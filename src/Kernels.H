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
using namespace amrex;

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
    Real sumY = 0.0;
    for (int n = 0; n < thermo.num_species; ++n) {
        species_massfractions[n] = prims(i, j, k, QY1 + n);
        sumY += species_massfractions[n];
    }
    // Normalize mass fractions
    for (int n = 0; n < thermo.num_species; ++n) {
        species_massfractions[n] /= sumY;
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
compute_flux_3D(int idir,
        amrex::GpuArray<amrex::Real, NPRIM> const&  u,
        amrex::GpuArray<amrex::Real, NCONS> const&  q,
        amrex::GpuArray<amrex::Real, NCONS> & flux) noexcept
{
    using amrex::Real;
    flux[0] = q[URHO] * u[idir + 1];
    flux[1] = q[URHO] * u[1] * u[idir + 1];
    flux[2] = q[URHO] * u[2] * u[idir + 1];
    flux[3] = q[URHO] * u[3] * u[idir + 1];
    flux[4] = u[idir + 1] * (q[UEINT] + u[QPRES]);
    for (int n = 0; n < NSP; ++n) {
        flux[5 + n] = q[URHO] * u[idir + 1] * u[QTEMP + 1 + n];
    }
    flux[idir + 1] += u[QPRES];
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
    amrex::GpuArray<Real, NCONS> fluxleft, fluxright;
    primtoc(uleft, qleft, thermo);
    primtoc(uright, qright, thermo);

    compute_flux_3D(idir, uleft, qleft, fluxleft);
    compute_flux_3D(idir, uright, qright, fluxright);


    /* Use the above Fortan code as reference to implement the Riemann solver */
    Real cl = compute_sos(uleft, thermo);  // Speed of sound on the left state
    Real cr = compute_sos(uright, thermo); // Speed of sound on the right state
                                                  //
    Real a1 = amrex::max(uleft[idir + 1] + cl, uright[idir + 1] + cr, Real(0.0));
    Real a2 = amrex::min(uleft[idir + 1] - cl, uright[idir + 1] - cr, Real(0.0));
    for(int icomp = 0; icomp < NCONS; ++icomp){
        flux[icomp] = (a1 * fluxleft[icomp] - a2 * fluxright[icomp] +
                       a1 * a2 * (qright[icomp] - qleft[icomp])) / (a1 - a2);
    }

//     // This function should compute the flux based on the left and right states
//     for(int icomp = 0; icomp < NCONS; ++icomp){
//         flux[icomp] = Real(0.0);
//     }
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
        uleft[icomp]  = prims(i-i1, j-j1, k-k1, icomp); //+ Real(0.5) * deriv(i-i1, j-j1, k-k1, icomp) * dx;
        uright[icomp] = prims(i, j, k, icomp); //- Real(0.5) * deriv(i, j, k, icomp) * dx;
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
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const&  dx,
        amrex::Real dt) noexcept
{
    using amrex::Real;
    Real dxinv = (Real) 1.0/ dx[0];
    Real dyinv = (Real) 1.0/ dx[1];
    Real dzinv = (Real) 1.0/ dx[2];

    for(int icomp=0; icomp<NCONS; ++icomp) {
        dq(i, j, k, icomp) = -dt *( (fluxx(i+1, j, k, icomp) - fluxx(i, j, k, icomp))*dxinv +
                (fluxy(i, j+1, k, icomp) - fluxy(i, j, k, icomp)) * dyinv +
                (fluxz(i, j, k+1, icomp) - fluxz(i, j, k, icomp)) * dzinv);
        // dq(i, j, k, icomp) = (Real) 0.0;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
updatedq(int i, int j, int k,
        amrex::Array4<amrex::Real> const& qnew,
        amrex::Array4<amrex::Real const> const& qold,
        amrex::Array4<amrex::Real const> const& dq,
        amrex::Real const1,
        amrex::Real const2
        ) noexcept
{
    using amrex::Real;
    for(int icomp=0; icomp<NCONS; ++icomp) {
        qnew(i, j, k, icomp) = const1 * (const2 * qnew(i, j, k, icomp) + qold(i, j, k, icomp) + dq(i, j, k, icomp));
    }
}
#endif
