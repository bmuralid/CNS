
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

void
AmrCoreCNS::AdvanceSingleStage (Real time, Real dt)

{
    Vector< Array <MultiFab, AMREX_SPACEDIM>> fluxes(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        for (int idim = 0; idim <= AMREX_SPACEDIM; ++idim)
        {
            BoxArray ba = grids[lev];
            ba.surroundingNodes(idim);
            fluxes[lev][idim] = MultiFab(ba, dmap[lev], 1, 0);
        }
    }

    /*----------------------------------------------------------------------*/
    for (int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab& mfprims =  qprims[lev];
        MultiFab& mfcons = qcons_new[lev];
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
                int idir;
                idir = 0;
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                   derivative(i, j, k, idir, 1, pfab, dfab, dx[idir]);
                });

                amrex::ParallelFor(fbx.grow(Direction::x, -1).surroundingNodes(Direction::x),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_recon_x(i, j, k, fluxx, pfab, dfab, dt);

                });

                // Flux reconstruction in Y-direction
                idir = 1;
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                   derivative(i, j, k, idir, 1, pfab, dfab, dx[idir]);
                });

                amrex::ParallelFor(fbx.grow(Direction::y, -1).surroundingNodes(Direction::y),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_recon_y(i, j, k, fluxy, pfab, dfab, dt);

                });

                // Flux reconstruction in Z-direction
                idir = 2;
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                   derivative(i, j, k, idir, 1, pfab, dfab, dx[idir]);
                });

                amrex::ParallelFor(fbx.grow(Direction::z, -1).surroundingNodes(Direction::z),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_recon_z(i, j, k, fluxz, pfab, dfab, dt);

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
        MultiFab& mfcons = dq[lev];
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
                Array4<Real> dqfab = mfcons[mfi].array();
                const Box& bx = mfi.tilebox();
                Box fbx = mfi.tilebox();
                AMREX_D_TERM(Array4<Real> fluxx = fluxes[lev][0].array(mfi);,
                             Array4<Real> fluxy = fluxes[lev][1].array(mfi);,
                             Array4<Real> fluxz = fluxes[lev][2].array(mfi));
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                   computedq(i, j, k, fluxx, fluxy, fluxz, dqfab, dt);
                });

            }
        }
    }
}



void
AmrCoreCNS::Advance (Real time, Real dt)
{
    /* Compute conservatives to primitives */
    // Use Q_old for conervative to primitive
    Cons2Prims(1);

    AdvanceSingleStage(time, dt);

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

