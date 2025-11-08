
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Interpolater.H>
#include <AMReX_Interp_C.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <cmath>

#include <AmrCoreCNS.H>
#include <Kernels.H>
#include <Prob.H>
#include <Thermo.H>
#include <Tagging.H>
#include <bc_fill.H>

using namespace amrex;
AmrCoreCNS::AmrCoreCNS ()
{
    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    if (do_subcycle) {
        for (int lev = 1; lev <= max_level; ++lev) {
            nsubsteps[lev] = MaxRefRatio(lev-1);
        }
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    qcons_new.resize(nlevs_max);
    qcons_old.resize(nlevs_max);
    dq.resize(nlevs_max);
    qprims.resize(nlevs_max);
    phi.resize(nlevs_max);
    rhs.resize(nlevs_max);
    ebfactory.resize(nlevs_max);

    int bc_lo[AMREX_SPACEDIM];
    int bc_hi[AMREX_SPACEDIM];

    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        if (Geom(0).isPeriodic()[idim] == 1) {
            bc_lo[idim] = bc_hi[idim] = BCType::int_dir;  // periodic
        } else {
            bc_lo[idim] = bc_hi[idim] = BCType::foextrap;  // walls (Neumann)
        }
    }

    bcs.resize(1);     // Setup 1-component
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        // lo-side BCs
        if (bc_lo[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_lo[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_lo[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setLo(idim, bc_lo[idim]);
        }
        else {
            amrex::Abort("Invalid bc_lo");
        }

        // hi-side BCSs
        if (bc_hi[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_hi[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_hi[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setHi(idim, bc_hi[idim]);
        }
        else {
            amrex::Abort("Invalid bc_hi");
        }
    }

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] is never actually used in the reflux operation
    flux_reg.resize(nlevs_max+1);
}

AmrCoreCNS::~AmrCoreCNS ()
{
}

// initializes multilevel data
void
AmrCoreCNS::InitData ()
{
    if (restart_chkfile == "") {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        AverageDown(2); // average down from finest level to level 0

        if (chk_int > 0) {
            WriteCheckpointFile();
        }

    }
    else {
        // restart from a checkpoint
        ReadCheckpointFile();
    }

    if (plot_int > 0) {
        WritePlotFile();
    }
}

void AmrCoreCNS::MakeNewLevelFromScratch(int lev, Real time, const BoxArray& ba,
        const DistributionMapping& dm)
{
    const int ncomp_cons = NCONS; // [\rho, \rho u, \rho v, \rho w, E]^T
    const int ncomp_prims = NPRIM; // [\rho, u, v, w, P, T]
    const int nghost = 2;

    qcons_new[lev].define(ba, dm, ncomp_cons, nghost);
    qcons_old[lev].define(ba, dm, ncomp_cons, nghost);
    dq[lev].define(ba, dm, ncomp_cons, nghost);

    qprims[lev].define(ba, dm, ncomp_prims, nghost);
    phi[lev].define(ba, dm, 1, nghost);
    rhs[lev].define(ba, dm, 1, nghost);

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp_cons));
    }

    MultiFab& state = qprims[lev];
    MultiFab& cons  = qcons_new[lev];

    const auto problo = Geom(lev).ProbLoArray();
    const auto dx     = Geom(lev).CellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> fab = state[mfi].array();
        Array4<Real> cfab = cons[mfi].array();
        const Box& box = mfi.growntilebox(nghost);

        amrex::launch(box,
        [=] AMREX_GPU_DEVICE (Box const& tbx)
        {
            initdata(tbx, fab, problo, dx);
        });
    }

    // Call prims to cons to fill in qcons_new
    Prims2Cons();


    // Initialize phi as signed distance to a sphere centered at (0,0,0) with radius 0.5
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> pfab = phi[lev][mfi].array();
        const Box& box = mfi.tilebox();

        amrex::ParallelFor(box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real x = problo[0] + (static_cast<Real>(i) + Real(0.5)) * dx[0];
#if (AMREX_SPACEDIM >= 2)
            const Real y = problo[1] + (static_cast<Real>(j) + Real(0.5)) * dx[1];
#else
            const Real y = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
            const Real z = problo[2] + (static_cast<Real>(k) + Real(0.5)) * dx[2];
#else
            const Real z = Real(0.0);
#endif
            const Real r = std::sqrt(x*x + y*y + z*z);
            pfab(i,j,k,0) = r - Real(0.5);
        });
    }

    // Build EB geometry (sphere) and factory
    {
        // EB2::SphereIF sphere(Real(0.5), {AMREX_D_DECL(Real(0.0), Real(0.0), Real(0.0))}, false);
        // auto gshop = EB2::makeShop(sphere);
        // EB2::Build(gshop, Geom(lev), 0, 0, 4);
        const Vector<int> ng_ebs{1,1,1};
        ebfactory[lev] = amrex::makeEBFabFactory(Geom(lev), ba, dm, ng_ebs, EBSupport::full);
    }
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
AmrCoreCNS::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                                    const DistributionMapping& dm)
{
    const int ncomp_cons = qcons_new[lev-1].nComp(); // [\rho, \rho u, \rho v, \rho w, E]^T
    const int ncomp_prims = qprims[lev-1].nComp(); // [\rho, u, v, w, P, T]
    const int nghost = qcons_new[lev-1].nGrow();

    qcons_new[lev].define(ba, dm, ncomp_cons, nghost);
    qcons_old[lev].define(ba, dm, ncomp_cons, nghost);
    dq[lev].define(ba, dm, ncomp_cons, nghost);
    qprims[lev].define(ba, dm, ncomp_prims, nghost);
    phi[lev].define(ba, dm, 1, nghost);
    rhs[lev].define(ba, dm, 1, nghost);

    setVal(dq[lev], 0.0);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp_cons));
    }

    FillCoarsePatch(lev, 2, time, qprims[lev], 0, ncomp_prims); // only need to copy over qprims

    // Build EB geometry (sphere) and factory on this level as well
    {
        // EB2::SphereIF sphere(Real(0.5), {AMREX_D_DECL(Real(0.0), Real(0.0), Real(0.0))}, false);
        // auto gshop = EB2::makeShop(sphere);
        // EB2::Build(gshop, Geom(lev), 0, 0, 4);
        const Vector<int> ng_ebs{1,1,1};
        ebfactory[lev] = amrex::makeEBFabFactory(Geom(lev), ba, dm, ng_ebs, EBSupport::full);
    }
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
AmrCoreCNS::RemakeLevel (int lev, Real time, const BoxArray& ba,
                         const DistributionMapping& dm)
{
    const int ncomp_cons = qcons_new[lev].nComp();
    const int ncomp_prims = qprims[lev].nComp();
    const int nghost = qcons_new[lev].nGrow();

    MultiFab new_qcons(ba, dm, ncomp_cons, nghost);
    MultiFab new_qprims(ba, dm, ncomp_prims, nghost);

    FillPatch(lev, 1, time, new_qcons, 0, ncomp_cons);
    FillPatch(lev, 2, time, new_qprims, 0, ncomp_prims);

    std::swap(new_qcons, qcons_new[lev]);
    std::swap(new_qprims, qprims[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp_cons));
    }
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
AmrCoreCNS::ClearLevel (int lev)
{
    qcons_new[lev].clear();
    qcons_old[lev].clear();
    dq[lev].clear();
    qprims[lev].clear();
    phi[lev].clear();
    rhs[lev].clear();
    flux_reg[lev].reset(nullptr);
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
AmrCoreCNS::ErrorEst (int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/)
{
    static bool first = true;
    static Vector<Real> phierr;

    // only do this during the first call to ErrorEst
    if (first)
    {
        first = false;
        // read in an array of "phierr", which is the tagging threshold
        // in this example, we tag values of "phi" which are greater than phierr
        // for that particular level
        // in subroutine state_error, you could use more elaborate tagging, such
        // as more advanced logical expressions, or gradients, etc.
        ParmParse pp("adv");
        int n = pp.countval("phierr");
        if (n > 0) {
            pp.getarr("phierr", phierr, 0, n);
        }
    }

    if (lev >= phierr.size()) return;

//    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const MultiFab& state = qprims[lev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
    {

        for (MFIter mfi(state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx  = mfi.tilebox();
            const auto statefab = state.array(mfi);
            const auto tagfab  = tags.array(mfi);
            Real phierror = phierr[lev];

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                state_error(i, j, k, tagfab, statefab, phierror, tagval);
            });
        }
    }
}


// read in some parameters from inputs file
void
AmrCoreCNS::ReadParameters ()
{
    {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
    }

    {
        ParmParse pp("amr"); // Traditionally, these have prefix, amr.

        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        pp.query("chk_file", chk_file);
        pp.query("chk_int", chk_int);
        pp.query("restart",restart_chkfile);
    }

    {
        ParmParse pp("adv");

        pp.query("cfl", cfl);
        pp.query("do_reflux", do_reflux);
        pp.query("do_subcycle", do_subcycle);
    }
}

// set covered coarse cells to be the average of overlying fine cells
void
AmrCoreCNS::AverageDown (int opt)
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        if (opt == 1) {
            amrex::average_down(qcons_new[lev+1], qcons_new[lev],
                    geom[lev+1], geom[lev],
                    0, qcons_new[lev].nComp(), refRatio(lev));
        } else if (opt == 2) {
            amrex::average_down(qprims[lev+1], qprims[lev],
                    geom[lev+1], geom[lev],
                    0, qprims[lev].nComp(), refRatio(lev));
        }
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
AmrCoreCNS::AverageDownTo (int crse_lev, int opt)
{
    if (opt == 1) {
        amrex::average_down(qcons_new[crse_lev+1], qcons_new[crse_lev],
                geom[crse_lev+1], geom[crse_lev],
                0, qcons_new[crse_lev].nComp(), refRatio(crse_lev));
    } else {
        amrex::average_down(qprims[crse_lev+1], qprims[crse_lev],
                geom[crse_lev+1], geom[crse_lev],
                0, qprims[crse_lev].nComp(), refRatio(crse_lev));
    }
}

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
AmrCoreCNS::FillPatch (int lev, int opt, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, opt, time, smf, stime);

        if(Gpu::inLaunchRegion())
        {
            GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > physbc(geom[lev],bcs,gpu_bndry_func);
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                        geom[lev], physbc, 0);
        }
        else
        {
            CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev],bcs,bndry_func);
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                        geom[lev], physbc, 0);
        }
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(lev-1, opt, time, cmf, ctime);
        GetData(lev  , opt, time, fmf, ftime);

        Interpolater* mapper = &cell_cons_interp;

        if(Gpu::inLaunchRegion())
        {
            GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[lev],bcs,gpu_bndry_func);

            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                      mapper, bcs, 0);
        }
        else
        {
            CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
            PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);

            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                      mapper, bcs, 0);
        }
    }
}


// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
AmrCoreCNS::FillCoarsePatch (int lev, int opt, Real time, MultiFab& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, opt, time, cmf, ctime);
    Interpolater* mapper = &cell_cons_interp;

    if (cmf.size() != 1) {
        amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    if(Gpu::inLaunchRegion())
    {
        GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[lev],bcs,gpu_bndry_func);

        amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                     cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                     mapper, bcs, 0);
    }
    else
    {
        CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);

        amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                     cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                     mapper, bcs, 0);
    }
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void
AmrCoreCNS::GetData (int lev, int opt, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (opt == 1) {
        data.push_back(&qcons_new[lev]);
    } else if (opt == 2) {
        data.push_back(&qprims[lev]);
    }
    datatime.push_back(t_new[lev]);
}


// get plotfile name
std::string
AmrCoreCNS::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

// put together an array of multifabs for writing
Vector<const MultiFab*>
AmrCoreCNS::PlotFileMF () const
{
    Vector<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
        r.push_back(&qprims[i]);
    }
    return r;
}

// set plotfile variables names
Vector<std::string>
AmrCoreCNS::PlotFileVarNames () const
{
    Vector<std::string> names;
    names.push_back("density");
    names.push_back("x_velocity");
    names.push_back("y_velocity");
#if (AMREX_SPACEDIM == 3)
    names.push_back("z_velocity");
#endif
    names.push_back("pressure");
    names.push_back("temperature");
    names.push_back("Y01");
    names.push_back("Y02");
    return names;
}

// write plotfile to disk
void
AmrCoreCNS::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);

    // Build per-level MultiFabs that include qprims + phi + EB volfrac and bndry area magnitude
    Vector<MultiFab> out_mf;
    out_mf.reserve(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        const int ncomp_q = qprims[lev].nComp();
        const int ncomp_out = ncomp_q + 3; // phi, vfrac, barea_mag
        MultiFab tmp(grids[lev], dmap[lev], ncomp_out, 0);

        // copy qprims into [0, ncomp_q)
        MultiFab::Copy(tmp, qprims[lev], 0, 0, ncomp_q, 0);
        // copy phi
        MultiFab::Copy(tmp, phi[lev], 0, ncomp_q, 1, 0);

        // copy EB volume fraction if available
        if (true) {
        /* if (ebfactory[lev]) { */
            MultiFab::Copy(tmp, ebfactory[lev]->getVolFrac(), 0, ncomp_q+1, 1, 0);

            // boundary area magnitude: sum of face area fracs as a simple scalar diagnostic
            MultiFab barea_mag(grids[lev], dmap[lev], 1, 0);
            barea_mag.setVal(0.0);
//             auto area = ebfactory[lev]->getAreaFrac();
// #ifdef AMREX_USE_OMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//             for (MFIter mfi(barea_mag,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
//                 const Box& bx = mfi.tilebox();
//                 auto bav = barea_mag.array(mfi);
//                 AMREX_D_TERM(auto apx = area[0]->const_array(mfi);,
//                              auto apy = area[1]->const_array(mfi);,
//                              auto apz = area[2]->const_array(mfi);)
//                 amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
//                     Real s = Real(0.0);
//                     s += apx(i,j,k);
// #if (AMREX_SPACEDIM >= 2)
//                     s += apy(i,j,k);
// #endif
// #if (AMREX_SPACEDIM == 3)
//                     s += apz(i,j,k);
// #endif
//                     bav(i,j,k,0) = s;
//                 });
//             }
            MultiFab::Copy(tmp, barea_mag, 0, ncomp_q+2, 1, 0);
        } else {
            tmp.setVal(0.0, ncomp_q+1, 2, 0);
        }

        out_mf.emplace_back(std::move(tmp));
    }

    // Build pointers vector expected by WriteMultiLevelPlotfile
    Vector<const MultiFab*> mf_ptrs;
    mf_ptrs.reserve(out_mf.size());
    for (auto& mf : out_mf) {
        mf_ptrs.push_back(&mf);
    }

    // Variable names: existing qprims names + phi + vfrac + barea_mag
    auto varnames = PlotFileVarNames();
    varnames.push_back("phi");
    varnames.push_back("vfrac");
    varnames.push_back("barea_mag");

    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf_ptrs, varnames,
                                   Geom(), t_new[0], istep, refRatio());
}

// write checkpoint file to disk
void
AmrCoreCNS::WriteCheckpointFile () const
{
    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains the global data, namely,
    //                     the evolution of the state data, the grid
    //                     information, and the time
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories contain the MultiFab data at each level of refinement

    const std::string& checkpointname = amrex::Concatenate(chk_file,istep[0],5);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subdirPrefix_0 .. dirName/subdirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
    if (ParallelDescriptor::IOProcessor()) {

        std::string HeaderFileName(checkpointname + "/Header");
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                        std::ofstream::trunc |
                        std::ofstream::binary);
        if( ! HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);

        // write out title line
        HeaderFile << "Checkpoint file for CNS\n";

        // write out finest_level
        HeaderFile << finest_level << "\n";

        // write out array of istep
        for (int i = 0; i < istep.size(); ++i) {
            HeaderFile << istep[i] << " ";
        }
        HeaderFile << "\n";

        // write out array of dt
        for (int i = 0; i < dt.size(); ++i) {
            HeaderFile << dt[i] << " ";
        }
        HeaderFile << "\n";

        // write out array of t_new
        for (int i = 0; i < t_new.size(); ++i) {
            HeaderFile << t_new[i] << " ";
        }
        HeaderFile << "\n";

        // write the BoxArray at each level
        for (int lev = 0; lev <= finest_level; ++lev) {
            boxArray(lev).writeOn(HeaderFile);
            HeaderFile << '\n';
        }
    }

    // write the MultiFab data to, e.g., chk00010/Level_0/
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Write(qcons_new[lev],
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "qcons_new"));
    }
}

// read checkpoint file from disk
void
AmrCoreCNS::ReadCheckpointFile ()
{
    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    Gpu::streamSynchronize();

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            istep[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            t_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        Gpu::streamSynchronize();

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // set BoxArray grids and DistributionMapping dmap in AMReX_AmrCore.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFab and FluxRegister data
        int ncomp_cons = 5;
        int ncomp_prims = 6;
        int nghost = 2;

        qcons_new[lev].define(ba, dm, ncomp_cons, nghost);
        qcons_old[lev].define(ba, dm, ncomp_cons, nghost);
        dq[lev].define(ba, dm, ncomp_cons, nghost);
        qprims[lev].define(ba, dm, ncomp_prims, nghost);
        phi[lev].define(ba, dm, 1, nghost);
        rhs[lev].define(ba, dm, 1, nghost);

        if (lev > 0 && do_reflux) {
            flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp_cons));
        }
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(qcons_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "qcons_new"));
    }
}
