
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

using namespace amrex;

void
AmrCoreCNS::Advance (Real time, Real dt)
{
    /* Compute conservatives to primitives */
    Cons2Prims();


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

