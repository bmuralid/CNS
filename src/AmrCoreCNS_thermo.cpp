#include <cmath>

#include <AmrCoreCNS.H>
#include <Kernels.H>
#include <pyro.H>
#include <Thermo.H>

using namespace amrex;

void
AmrCoreCNS::Cons2Prims (int opt)
{
    /* Compute conservatives to primitives */
    pyro::pyro<double> const lpyro  = thermo;

    for (int lev = 0; lev <= finest_level; lev++) {
        MultiFab& mfprims =  qprims[lev];
        MultiFab& mfcons = opt == 0 ? qcons_old[lev] : qcons_new[lev];
        const int ng = mfcons.nGrow();
        // mfcons.FillBoundary(Geom(lev).periodicity());

        for (MFIter mfi(mfcons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<Real> cfab = mfcons[mfi].array();
            Array4<Real> pfab = mfprims[mfi].array();

            // const Box& bx = mfi.growntilebox(ng);
            const Box& bx = mfi.tilebox();

            ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cns_ctoprim(i, j, k, cfab, pfab, lpyro);

            });
        }
    }
}

void
AmrCoreCNS::Prims2Cons ()
{
    /* Compute primitives to primitives*/
    pyro::pyro<double> const lpyro  = thermo;

    for (int lev = 0; lev <= finest_level; lev++) {
        MultiFab& mfprims =  qprims[lev];
        MultiFab& mfcons = qcons_new[lev];

        const int ng = mfcons.nGrow();
        for (MFIter mfi(mfcons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<Real> cfab = mfcons[mfi].array();
            Array4<Real> pfab = mfprims[mfi].array();

            // const Box& bx = mfi.growntilebox(ng);
            const Box& bx = mfi.tilebox();

            ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cns_primtoc(i, j, k, pfab, cfab, lpyro);

            });
        }
    }
}
