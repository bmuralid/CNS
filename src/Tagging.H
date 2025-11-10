#ifndef TAGGING_H
#define TAGGING_H

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include "CNS_index_macros.H"
#include <cmath>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
compute_normalized_curvature(int i, int j, int k,
                              amrex::Array4<amrex::Real const> const& state,
                              int comp,
                              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx) noexcept
{
    using amrex::Real;

    // Compute second derivatives (curvature) using central differences
    // Laplacian: ∇²f = f_xx + f_yy + f_zz

    Real f_xx = Real(0.0);
    Real f_yy = Real(0.0);
    Real f_zz = Real(0.0);

    // Second derivative in x-direction: f_xx ≈ (f(i+1) - 2*f(i) + f(i-1)) / dx²
    Real f_curr = state(i, j, k, comp);
    Real f_xp = state(i+1, j, k, comp);
    Real f_xm = state(i-1, j, k, comp);
    f_xx = (f_xp - Real(2.0) * f_curr + f_xm) / (dx[0] * dx[0]);

    // Second derivative in y-direction
    Real f_yp = state(i, j+1, k, comp);
    Real f_ym = state(i, j-1, k, comp);
    f_yy = (f_yp - Real(2.0) * f_curr + f_ym) / (dx[1] * dx[1]);

#if (AMREX_SPACEDIM == 3)
    // Second derivative in z-direction
    Real f_zp = state(i, j, k+1, comp);
    Real f_zm = state(i, j, k-1, comp);
    f_zz = (f_zp - Real(2.0) * f_curr + f_zm) / (dx[2] * dx[2]);
#endif

    // Laplacian (total curvature)
    Real laplacian = f_xx + f_yy;
#if (AMREX_SPACEDIM == 3)
    laplacian += f_zz;
#endif

    // Normalize by local field value to make it dimensionless
    // Use absolute value of field with small epsilon to avoid division by zero
    Real field_mag = std::abs(f_curr);
    Real epsilon = Real(1.0e-10);

    // Normalized curvature: |∇²f| * dx² / (|f| + epsilon)
    // This gives a measure of relative curvature that is scale-independent
    Real dx_min = amrex::min(dx[0], dx[1]);
#if (AMREX_SPACEDIM == 3)
    dx_min = amrex::min(dx_min, dx[2]);
#endif

    Real normalized_curvature = std::abs(laplacian) * dx_min * dx_min / (field_mag + epsilon);

    return normalized_curvature;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
state_error (int i, int j, int k,
             amrex::Array4<char> const& tag,
             amrex::Array4<amrex::Real const> const& state,
             amrex::Real phierr, char tagval,
             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx)
{
    using amrex::Real;

    // Compute normalized curvature for pressure and density
    /* Real curvature_pressure = compute_normalized_curvature(i, j, k, state, QPRES, dx); */
    /* Real curvature_density = compute_normalized_curvature(i, j, k, state, QRHO, dx); */

    // Tag cell if normalized curvature exceeds threshold for either field
    // phierr is now interpreted as the normalized curvature threshold
    /* if (curvature_pressure > phierr || curvature_density > phierr) */
    /* { */
    /*     tag(i, j, k) = tagval; */
    /* } */
    if (state(i, j, k, QRHO) > 1.8) {
        tag(i, j, k) = tagval;
    }
}

#endif
