#ifndef PROB_H_
#define PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
initdata(amrex::Box const& bx, amrex::Array4<amrex::Real> const& fab,
         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& problo,
         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    using namespace amrex;

    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real x = problo[0] + (i+0.5) * dx[0];
                Real y = problo[1] + (j+0.5) * dx[1];
#if (AMREX_SPACEDIM == 3)
                Real z = problo[2] + (k+0.5) * dx[2];
#else
                Real z = 0.0;
#endif

                // Initialize primitive variables: [rho, u, v, w, P, T]
                // Simple initialization - uniform density and pressure
                fab(i,j,k,0) = 1.225;  // density
                fab(i,j,k,1) = 0.0;  // u velocity
                fab(i,j,k,2) = 0.0;  // v velocity
#if (AMREX_SPACEDIM == 3)
                fab(i,j,k,3) = 0.0;  // w velocity
                fab(i,j,k,4) = 101325.0;  // pressure
                fab(i,j,k,5) = 300.0;  // temperature
                fab(i, j, k, 6) = 0.767;
                fab(i, j, k, 7) = 0.233; //
#else
                fab(i,j,k,3) = 101325.0;  // pressure
                fab(i,j,k,4) = 300.0;  // temperature
#endif
            }
        }
    }
}

#endif
