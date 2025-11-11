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
    Real laplacian = f_xx * f_xx + f_yy * f_yy;
#if (AMREX_SPACEDIM == 3)
    laplacian += f_zz * f_zz;
#endif
    laplacian = std::sqrt(laplacian);

    // Normalize by local field value to make it dimensionless
    // Use absolute value of field with small epsilon to avoid division by zero
    Real field_mag = std::abs(f_xm + (Real)2.0 * f_curr + f_xp);
    field_mag += std::abs(f_ym + (Real) 2.0 * f_curr + f_yp);
#if (AMREX_SPACEDIM == 3)
    field_mag += std::abs(f_zm + (Real) 2.0 * f_curr + f_zp);
#endif
    Real epsilon = Real(1.0e-10);

    // Normalized curvature: |∇²f| * dx² / (|f| + epsilon)
    // This gives a measure of relative curvature that is scale-independent
    Real dx_min = amrex::min(dx[0], dx[1]);
#if (AMREX_SPACEDIM == 3)
    dx_min = amrex::min(dx_min, dx[2]);
#endif

    Real normalized_curvature = std::abs(laplacian) * dx_min/ (field_mag + epsilon);

    return normalized_curvature;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
compute_curvature_leslie(int i, int j, int k,
                         int idir,
                         amrex::Real eps, amrex::Real threshold,
                         amrex::Array4<amrex::Real const> const& state,
                         int comp,
                         amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx) noexcept
{
    using amrex::Real;

    // Second derivative in x-direction: f_xx ≈ (f(i+1) - 2*f(i) + f(i-1)) / dx²
    int i1 = idir == 0 ? 1 : 0;
    int j1 = idir == 1 ? 1 : 0;
    int k1 = idir == 2 ? 1 : 0;

    Real f_curr = state(i, j, k, comp);
    Real f_xp = state(i+i1, j+j1, k+k1, comp);
    Real f_xm = state(i-i1, j-j1, k-k1, comp);
    Real delta_xm = f_curr - f_xm;
    Real delta_xp = f_xp - f_curr;
    Real dif1 = std::abs(delta_xp - delta_xm);
    Real dif2 = std::abs(delta_xp) + std::abs(delta_xm);
    Real test1 = dif1 - eps * f_curr;
    Real test2 = dif1 - threshold * dif2;
    return (test1 > 0.0 && test2 > 0.0) ? 1.0 : 0.0;
}



AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
compute_sensor (int i, int j, int k,
                amrex::Array4<amrex::Real> const& tag,
                amrex::Array4<amrex::Real const> const& state,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx
               )
{
    using amrex::Real;

    // // Compute normalized curvature for pressure and density
    // Real curvature_pressure = compute_normalized_curvature(i, j, k, state, QPRES, dx);
    // Real curvature_density = compute_normalized_curvature(i, j, k, state, QRHO, dx);

    // // Tag cell if normalized curvature exceeds threshold for either field
    // const Real threshold = 0.1; // Example threshold value; adjust as needed
    // if (curvature_pressure >= threshold || curvature_density >= threshold)
    // {
    //     tag(i, j, k) = curvature_pressure; // Tag value
    // }
    Real epsP = 0.05;
    Real threshP = 0.5;
    Real s1 = compute_curvature_leslie(i, j, k, 0, epsP, threshP, state, QPRES, dx);
    Real s2 = compute_curvature_leslie(i, j, k, 1, epsP, threshP, state, QPRES, dx);
#if (AMREX_SPACEDIM == 3)
    Real s3 = compute_curvature_leslie(i, j, k, 2, epsP, threshP, state, QPRES, dx);
    tag(i, j, k) = amrex::max(s1, amrex::max(s2, s3));
#else
    tag(i, j, k) = amrex::max(s1, s2);
#endif

    Real epsRho = 0.05;
    Real threshRho = 0.125;
    s1 = compute_curvature_leslie(i, j, k, 0, epsRho, threshRho, state, QRHO, dx);
    s2 = compute_curvature_leslie(i, j, k, 1, epsRho, threshRho, state, QRHO, dx);
#if (AMREX_SPACEDIM == 3)
    s3 = compute_curvature_leslie(i, j, k, 2, epsRho, threshRho, state, QRHO, dx);
    tag(i, j, k) = amrex::max(tag(i, j, k), amrex::max(s1, amrex::max(s2, s3)));
#else
    tag(i, j, k) = amrex::max(tag(i, j, k), amrex::max(s1, s2));
#endif
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
state_error (int i, int j, int k,
             bool first,
             amrex::Array4<char> const& tag,
             amrex::Array4<amrex::Real const> const& state,
             amrex::Real phierr, char tagval, char clearval
            )
{
    using amrex::Real;

    if (state(i, j, k, 0) > 0.5) {
        tag(i, j, k) = tagval;
    } else {
        tag(i, j, k) = clearval;
    }
}

#endif
