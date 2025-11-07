#ifndef Kernels_H_
#define Kernels_H_

#include <AMReX_FArrayBox.H>
#include "CNS_index_macros.H"
#include "pyro.H"
#include <AMReX_Print.H>
#include "Thermo.H"
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


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
monocent(amrex::Real a, amrex::Real b, amrex::Real c) noexcept
{
    using amrex::Real;

    if ( a > Real(0.0) && b > Real(0.0) && c > Real(0.0) ) {
        return amrex::min( amrex::min(a, b), c );
    } else if ( a < Real(0.0) && b < Real(0.0) && c < Real(0.0) ) {
        return amrex::max( amrex::max(a, b), c );
    } else {
        return Real(0.0);
    }
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
derivative(int i, int j, int k,
        int idir,
        int isLimited,
        amrex::Array4<amrex::Real> const& prims,
        amrex::Array4<amrex::Real> const& deriv,
        amrex::Real dx) noexcept
{
    using amrex::Real;

    Real inv_2dx = Real(0.5) / dx;
    Real inv_dx = Real(1.0) / dx;
    Real delta_plus = Real(0.0);
    Real delta_central = Real(0.0);
    Real delta_minus = Real(0.0);

    const int i1 = idir== 0 ? 1 : 0;
    const int j1 = idir== 1 ? 1 : 0;
    const int k1 = idir== 2 ? 1 : 0;

    for (int icomp = 0; icomp < NPRIM; ++icomp){
        delta_minus = (prims(i, j, k, icomp) - prims(i-i1, j-j1, k-k1, icomp)) * inv_dx;
        delta_central = (prims(i+i1, j+j1, k+k1, icomp) - prims(i-i1, j-j1, k-k1, icomp)) * inv_2dx;
        delta_plus = (prims(i+i1, j+j1, k+k1, icomp) - prims(i, j, k, icomp)) * inv_dx;
        deriv(i, j, k, icomp) = isLimited == 1 ? monocent(delta_minus, delta_central, delta_plus) : delta_central;
    }
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
compute_riemann_flux(
        amrex::GpuArray<amrex::Real, NPRIM> const&  uleft,
        amrex::GpuArray<amrex::Real, NPRIM> const&  uright,
        const int idir,
        amrex::GpuArray<amrex::Real, NCONS> & flux,
        pyro::pyro<double> const& thermo) noexcept
{
    using amrex::Real;

    amrex::GpuArray<Real, NCONS> qleft, qright;
    // primtoc(uleft, qleft, thermo);
    // primtoc(uright, qright, thermo);

      /* call sos(ul, cl) */
      /* call sos(ur, cr) */

      /* if (idir == 1) then */
      /*    fl(1) = ul(irho)*ul(iu) */
      /*    fl(2) = ul(irho)*ul(iu)*ul(iu) + ul(ip) */
      /*    fl(3) = ul(irho)*ul(iu)*ul(iv) */
      /*    fl(4) = ul(iu)*(ql(iRhoE) + ul(ip)) */

      /*    fr(1) = ur(irho)*ur(iu) */
      /*    fr(2) = ur(irho)*ur(iu)*ur(iu) + ur(ip) */
      /*    fr(3) = ur(irho)*ur(iu)*ur(iv) */
      /*    fr(4) = ur(iu)*(qr(iRhoE) + ur(ip)) */
      /* else */
      /*    fl(1) = ul(irho)*ul(iv) */
      /*    fl(2) = ul(irho)*ul(iu)*ul(iv) */
      /*    fl(3) = ul(irho)*ul(iv)*ul(iv) + ul(ip) */
      /*    fl(4) = ul(iv)*(ql(iRhoE) + ul(ip)) */

      /*    fr(1) = ur(irho)*ur(iv) */
      /*    fr(2) = ur(irho)*ur(iu)*ur(iv) */
      /*    fr(3) = ur(irho)*ur(iv)*ur(iv) + ur(ip) */
      /*    fr(4) = ur(iv)*(qr(iRhoE) + ur(ip)) */
      /* end if */

      /* a1 = max(ul(irho + idir) + cl, ur(irho + idir) + cr, 0.0d0) */
      /* a2 = min(ul(irho + idir) - cl, ur(irho + idir) - cr, 0.0d0) */

      /* do l = 1, ncons */
      /*    f(l) = (a1*fl(l) - a2*fr(l) + a1*a2*(qr(l) - ql(l)))/(a1 - a2) */
      /* end do */

    /* Use the above Fortan code as reference to implement the Riemann solver */
    Real cl = Real(0.0); // Speed of sound on the left state
    Real cr = Real(0.0); // Speed of sound on the right state
    amrex::GpuArray<Real, NSP> species_massfractions = {0.0};
    for (int n = 0; n < thermo.num_species; ++n) {
        species_massfractions[n] = uleft[QTEMP + 1 + n];
    }
    Real temperature_left = uleft[QTEMP];
    // cl = sqrt(thermo.get_mixture_specific_heat_cp_mass(temperature_left, species_massfractions) *
    //           thermo.get_specific_gas_constant(species_massfractions) /
    //           thermo.get_mixture_specific_heat_cv_mass(temperature_left, species_massfractions) *
    //             temperature_left);
    // for (int n = 0; n < thermo.num_species; ++n) {
    //     species_massfractions[n] = uright[QTEMP + 1 + n];
    // }
    // Real temperature_right = uright[QTEMP];
    // cr = sqrt(thermo.get_mixture_specific_heat_cp_mass(temperature_right, species_massfractions) *
    //           thermo.get_specific_gas_constant(species_massfractions) /
    //           thermo.get_mixture_specific_heat_cv_mass(temperature_right, species_massfractions) *
    //             temperature_right);


    // Convert left and right primitive states to conservative states


    // Placeholder for Riemann solver implementation
    // This function should compute the flux based on the left and right states
    // For now, we set flux to zero
    for(int icomp = 0; icomp < NCONS; ++icomp){
        flux[icomp] = Real(0.0);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
flux_recon(int i, int j, int k,
        amrex::Array4<amrex::Real> const& fluxx,
        amrex::Array4<amrex::Real> const& prims,
        amrex::Array4<amrex::Real> const& deriv,
        const int idir,
        pyro::pyro<double> const& thermo,
        amrex::Real dt, amrex::Real dx) noexcept
{
    using amrex::Real;
    amrex::GpuArray<Real, NPRIM> uleft, uright;
    amrex::GpuArray<Real, NCONS> flux;

    const int i1 = idir== 0 ? 1 : 0;
    const int j1 = idir== 1 ? 1 : 0;
    const int k1 = idir== 2 ? 1 : 0;
    // Perform reconstruction of the left and right states using the derivative information
    for(int icomp = 0; icomp < NPRIM; ++icomp){
        uleft[icomp]  = prims(i-i1, j-j1, k-k1, icomp) + Real(0.5) * deriv(i-i1, j-j1, k-k1, icomp) * dx;
        uright[icomp] = prims(i, j, k, icomp) - Real(0.5) * deriv(i, j, k, icomp) * dx;
    }

    compute_riemann_flux(uleft, uright, idir, flux, thermo);

    // Store the computed fluxes
    for(int icomp = 0; icomp < NCONS; ++icomp){
        fluxx(i, j, k, icomp) = flux[icomp];
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
flux_recon_y(int i, int j, int k,
        amrex::Array4<amrex::Real> const& fluxy,
        amrex::Array4<amrex::Real> const& prims,
        amrex::Array4<amrex::Real> const& deriv,
        amrex::GpuArray<amrex::Real, NPRIM> const&  uleft,
        amrex::GpuArray<amrex::Real, NPRIM> const&  uright,
        amrex::GpuArray<amrex::Real, NCONS> const&  qleft,
        amrex::GpuArray<amrex::Real, NCONS> const&  qright,
        pyro::pyro<double> const& thermo,
        amrex::Real dt) noexcept
{
    using amrex::Real;

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
flux_recon_z(int i, int j, int k,
        amrex::Array4<amrex::Real> const& fluxz,
        amrex::Array4<amrex::Real> const& prims,
        amrex::Array4<amrex::Real> const& deriv,
        amrex::GpuArray<amrex::Real, NPRIM> const&  uleft,
        amrex::GpuArray<amrex::Real, NPRIM> const&  uright,
        amrex::GpuArray<amrex::Real, NCONS> const&  qleft,
        amrex::GpuArray<amrex::Real, NCONS> const&  qright,
        pyro::pyro<double> const& thermo,
        amrex::Real dt) noexcept
{
    using amrex::Real;

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
computedq(int i, int j, int k,
        amrex::Array4<amrex::Real> const& fluxx,
        amrex::Array4<amrex::Real> const& fluxy,
        amrex::Array4<amrex::Real> const& fluxz,
        amrex::Array4<amrex::Real> const& dq,
        amrex::Real dt) noexcept
{
    using amrex::Real;

}

#endif
