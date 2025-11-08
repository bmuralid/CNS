
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
#include <Thermo.H>
#include <Kernels.H>

#include <CNS_index_macros.H>
using namespace amrex;

/* time integration constants used `updatedq` */

void
AmrCoreCNS::AdvanceSingleStage (Real time, Real dt, int istage)

{
    const std::array<amrex::Real, 2> const1 = {1.0, 0.5};
    const std::array<amrex::Real, 2> const2 = {0.0, 1.0};
    Vector< Array <MultiFab, AMREX_SPACEDIM>> fluxes(finest_level + 1);
    pyro::pyro<double> const lpyro = thermo;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            BoxArray ba = grids[lev];
            ba.surroundingNodes(idim);
            fluxes[lev][idim] = MultiFab(ba, dmap[lev], NCONS, 0);
        }
    }

    /*----------------------------------------------------------------------*/
    for (int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab& mfprims =  qprims[lev];
        MultiFab& mfcons = qcons_new[lev];
        const auto dx = geom[lev].CellSizeArray();
        const amrex::Array<Real, 3> dx_local = {dx[0], dx[1], dx[2]};
        AMREX_D_TERM(Real dtdx = dt/dx[0];,
                     Real dtdy = dt/dx[1];,
                     Real dtdz = dt/dx[2]);
        // amrex::Print() << qprims[lev].boxArray() <<"\n";
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            FArrayBox tmpfab;
            for (MFIter mfi(qprims[lev], false); mfi.isValid(); ++mfi)
            {
                Array4<Real> cfab = mfcons[mfi].array();
                Array4<Real> pfab = mfprims[mfi].array();

                const Box& bx = mfi.growntilebox(1);
                FArrayBox derv(bx, NPRIM, The_Async_Arena());

                Box fbx = mfi.tilebox();

                Array4<Real> dfab = derv.array();

                AMREX_D_TERM(Array4<Real> fluxx = fluxes[lev][0].array(mfi);,
                             Array4<Real> fluxy = fluxes[lev][1].array(mfi);,
                             Array4<Real> fluxz = fluxes[lev][2].array(mfi));

                // Flux reconstruction in X-direction
                int idir = 0;
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                   derivative(i, j, k, idir, 1, pfab, dfab, dx_local[idir]);
                });

                // amrex::ParallelFor(fbx.grow(Direction::x, -1).surroundingNodes(Direction::x),
                //
                amrex::ParallelFor(amrex::surroundingNodes(fbx, Direction::x),                // amrex::ParallelFor(fbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {

                    flux_recon(i, j, k, fluxx, pfab, dfab, idir, lpyro, dt, dx[idir]);

                });

                // Flux reconstruction in Y-direction
                idir = 1;

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                   derivative(i, j, k, idir, 1, pfab, dfab, dx[idir]);
                });

                // amrex::Print()<<fluxy <<fbx.surroundingNodes(Direction::y) <<"\n";

                amrex::ParallelFor(amrex::surroundingNodes(fbx, Direction::y),                // amrex::ParallelFor(fbx,
                // amrex::ParallelFor(fbx.surroundingNodes(Direction::y),
                // amrex::ParallelFor(fbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_recon(i, j, k, fluxy, pfab, dfab, idir, lpyro, dt, dx[idir]);

                });

                // Flux reconstruction in Z-direction
                idir = 2;
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                   derivative(i, j, k, idir, 1, pfab, dfab, dx[idir]);
                });

                amrex::ParallelFor(amrex::surroundingNodes(fbx, Direction::z),                // amrex::ParallelFor(fbx,
                // amrex::ParallelFor(fbx.surroundingNodes(Direction::z),
                // amrex::ParallelFor(fbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_recon(i, j, k, fluxz, pfab, dfab, idir, lpyro, dt, dx[idir]);

                });
            }
        }
    }

    // =======================================================
    // Average down the fluxes before using them to update phi
    // =======================================================
    for (int lev = finest_level; lev > 0; lev--)
    {
       average_down_faces(amrex::GetArrOfConstPtrs(fluxes[lev  ]),
                          amrex::GetArrOfPtrs     (fluxes[lev-1]),
                          refRatio(lev-1), Geom(lev-1));
    }


    /*----------------------------------------------------------------------*/
    for (int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab& mfdq = dq[lev];
        MultiFab& mfold = qcons_old[lev];
        MultiFab& mfnew = qcons_new[lev];
        const auto dx = geom[lev].CellSizeArray();
        AMREX_D_TERM(Real dtdx = dt/dx[0];,
                     Real dtdy = dt/dx[1];,
                     Real dtdz = dt/dx[2]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            FArrayBox tmpfab;
            for (MFIter mfi(qprims[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Array4<Real> dqfab = mfdq[mfi].array();
                Array4<Real> qnew = mfnew[mfi].array();
                Array4<Real> qold = mfold[mfi].array();
                const Box& bx = mfi.tilebox();
                Box fbx = mfi.tilebox();
                AMREX_D_TERM(Array4<Real> fluxx = fluxes[lev][0].array(mfi);,
                             Array4<Real> fluxy = fluxes[lev][1].array(mfi);,
                             Array4<Real> fluxz = fluxes[lev][2].array(mfi));
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                   computedq(i, j, k, fluxx, fluxy, fluxz, dqfab, dx, dt);
                   updatedq(i, j, k, qnew, qold, dqfab, const1[istage], const2[istage]);

                });
            }
        }
    }
}



void
AmrCoreCNS::Advance (Real time, Real dt)
{

    for (int lev=0; lev <= finest_level; lev++){
        std::swap(qcons_old[lev], qcons_new[lev]);
    }

    for (int istage=0; istage<1; istage++){
        // Need to call BC function here for primitives
        AdvanceSingleStage(time, dt, istage);
        AverageDown(1);
        Cons2Prims(1);
    }

    /* to compute the fluxes with refluxing */
/*     for (int lev = 0; lev <= finest_level; lev++) { */
/*         MultiFab& mfprims =  qprims[lev]; */
/*         MultiFab& mfcons = qcons_new[lev]; */

/*         for (MFIter mfi(mfcons, TilingIfNotGPU()); mfi.isValid(); ++mfi) */
/*         { */
/*             Array4<Real> cfab = mfcons[mfi].array(); */
/*             Array4<Real> pfab = mfprims[mfi].array(); */



/*             const Box& bx = mfi.growntilebox(1); */

/*             FArrayBox derv(bx, 3, The_Async_Arena()); */

/*             Array4<Real> darr = derv.array(); */

/*             // NTBC: Check if this is needed */
/*             /1* const Box& gbx = amrex::grow(bx, 1); *1/ */

/*             ParallelFor(bx, */
/*             [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept */
/*             { */
/*                 cns_deriv(i, j, k, cfab, darr); */


/*             }); */


/*         } */
/*     } */
}
