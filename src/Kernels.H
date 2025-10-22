#ifndef Kernels_H_
#define Kernels_H_

#include <AMReX_FArrayBox.H>
#include "CNS_index_macros.H"
#include "pyro.H"

/* #if (AMREX_SPACEDIM == 2) */
/* #include <compute_flux_2D_K.H> */
/* #else */
/* #include <compute_flux_3D_K.H> */
/* #endif */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
cns_computedt(int i, int j, int k,
        amrex::Array4<amrex::Real> const& prims,
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const&  dx,
        pyro::pyro<double> const& thermo) noexcept
{
    using amrex::Real;

    Real rho = prims(i, j, k, QRHO);
    Real ux = prims(i, j, k, QU);
    Real uy = prims(i, j, k, QV);
    Real uz = prims(i, j, k, QW);

    std::array<double, NSP> species_massfractions = {0.0};
    for (int n = 0; n < thermo.num_species; ++n) {
        species_massfractions[n] = prims(i, j, k, QY1 + n);
    }
    Real temperature = prims(i, j, k, QTEMP);

    Real c = sqrt(thermo.get_mixture_specific_heat_cp_mass(temperature, species_massfractions) *
                  thermo.get_specific_gas_constant(species_massfractions) /
                  thermo.get_mixture_specific_heat_cv_mass(temperature, species_massfractions) *
                  temperature);

    Real dt_min = Real(1.0e10);
#if (AMREX_SPACEDIM == 2)
    Real dt_x = dx[0] / (std::abs(ux) + c);
    Real dt_y = dx[1] / (std::abs(uy) + c);
    dt_min = amrex::min(dt_x, dt_y);
#else
    Real dt_x = dx[0] / (std::abs(ux) + c);
    Real dt_y = dx[1] / (std::abs(uy) + c);
    Real dt_z = dx[2] / (std::abs(uz) + c);
    dt_min = amrex::min(dt_x, dt_y);
    dt_min = amrex::min(dt_min, dt_z);
#endif
    return dt_min;
}
#endif
