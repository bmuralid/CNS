#ifndef THERMO_H
#define THERMO_H

#include <AMReX_Array4.H>
#include <array>
#include <AMReX_Print.H>
#include "pyro.H"
#include "CNS_index_macros.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
cns_ctoprim (int i, int j, int k,
        amrex::Array4<amrex::Real const> const& cons,
        amrex::Array4<amrex::Real> const& prims,
        pyro::pyro<double> const& thermo)
{
    using amrex::Real;

    // Real rho = cons(i, j, k, URHO);
    // Recompute denisty from species mass fractions for consistency
    Real rho = Real(0.0);
    for (int n = 0; n < thermo.num_species; ++n) {
        rho += cons(i, j, k, UEINT + 1 + n);
    }

    Real rhoinv = Real(0.0);
    if (rho > Real(0.0))
        rhoinv = Real(1.0)/rho;

    Real ux = cons(i, j, k, UMX) * rhoinv;
    Real uy = cons(i, j, k, UMY) * rhoinv;
    Real uz = cons(i, j, k, UMZ) * rhoinv;

    Real kineng = Real(0.5) * rho * (ux*ux + uy*uy + uz*uz);
    Real ei = (cons(i, j, k, UEINT) - kineng) * rhoinv;

    std::array<double, NSP> species_massfractions = {0.767, 0.233}; // default values
    Real sumY = Real(0.0);
    for (int n = 0; n < thermo.num_species; ++n) {
        species_massfractions[n] = cons(i, j, k, UEINT + 1 + n) * rhoinv;
        species_massfractions[n] = amrex::max(species_massfractions[n], Real(0.0));
        species_massfractions[n] = amrex::min(species_massfractions[n], Real(1.0));
        sumY += species_massfractions[n];
    }
    // Normalize mass fractions
    for (int n = 0; n < thermo.num_species; ++n) {
        species_massfractions[n] /= sumY;
    }
    Real temperature = thermo.get_temperature(ei, prims(i, j, k, QTEMP), species_massfractions);
    Real pressure = thermo.get_pressure(rho, temperature, species_massfractions);

    prims(i, j, k, QRHO) = rho;
    prims(i, j, k, QU)  = ux;
    prims(i, j, k, QV)  = uy;
    prims(i, j, k, QW)  = uz;
    prims(i, j, k, QPRES) = pressure;
    prims(i, j, k, QTEMP) = temperature;
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
cns_primtoc (int i, int j, int k,
        amrex::Array4<amrex::Real const> const& prims,
        amrex::Array4<amrex::Real> const& cons,
        pyro::pyro<double> const& thermo)
{
    using amrex::Real;

    Real rho = prims(i, j, k, QRHO);
    Real ux = prims(i, j, k, QU);
    Real uy = prims(i, j, k, QV);
    Real uz = prims(i, j, k, QW);
    Real pressure = prims(i, j, k, QPRES);
    Real temperature = prims(i, j, k, QTEMP);


    Real kineng = Real(0.5) * rho * (ux*ux + uy*uy + uz*uz);
    std::array<double, NSP> species_massfractions = {0.767, 0.233};
    Real sumY = 0.0;
    for (int n = 0; n < thermo.num_species; ++n) {
        species_massfractions[n] = prims(i, j, k, QTEMP + 1 + n);
        species_massfractions[n] = amrex::max(species_massfractions[n], Real(0.0));
        species_massfractions[n] = amrex::min(species_massfractions[n], Real(1.0));
        sumY += species_massfractions[n];
    }
    // Normalize mass fractions
    for (int n = 0; n < thermo.num_species; ++n) {
        species_massfractions[n] /= sumY;
    }
    // Recompute pressure from temperature to ensure consistency
    pressure = thermo.get_pressure(rho, temperature, species_massfractions);

    Real ei = thermo.get_mixture_internal_energy_mass(temperature, species_massfractions);
    Real eint = ei * rho;
    Real total_energy = eint + kineng;
    cons(i, j, k, URHO) = rho;
    cons(i, j, k, UMX)  = rho * ux;
    cons(i, j, k, UMY)  = rho * uy;
    cons(i, j, k, UMZ)  = rho * uz;
    cons(i, j, k, UEINT) = total_energy;
    for (int n = 0; n < thermo.num_species; ++n) {
        cons(i, j, k, UEINT + 1 + n) = rho * species_massfractions[n];
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
primtoc (amrex::GpuArray<amrex::Real, NPRIM> const& prims,
        amrex::GpuArray<amrex::Real, NCONS> & cons,
        pyro::pyro<double> const& thermo)
{
    using amrex::Real;

    Real rho = prims[QRHO];
    Real ux = prims[QU];
    Real uy = prims[QV];
    Real uz = prims[QW];
    Real pressure = prims[QPRES];
    Real temperature = prims[QTEMP];

    Real kineng = Real(0.5) * rho * (ux*ux + uy*uy + uz*uz);
    std::array<double, NSP> species_massfractions = {0.767, 0.233}; // default values
    Real sumY = (Real) 0.0;
    for (int n = 0; n < thermo.num_species; ++n) {
        species_massfractions[n] = prims[QTEMP + 1 + n];
        species_massfractions[n] = amrex::max(species_massfractions[n], Real(0.0));
        species_massfractions[n] = amrex::min(species_massfractions[n], Real(1.0));
        sumY += species_massfractions[n];
    }
    // Normalize mass fractions
    for (int n = 0; n < thermo.num_species; ++n) {
        species_massfractions[n] /= sumY;
    }
    // Recompute pressure from temperature to ensure consistency
    pressure = thermo.get_pressure(rho, temperature, species_massfractions);

    Real ei = thermo.get_mixture_internal_energy_mass(temperature, species_massfractions);
    Real eint = ei * rho;
    Real total_energy = eint + kineng;
    cons[URHO] = rho;
    cons[UMX]  = rho * ux;
    cons[UMY]  = rho * uy;
    cons[UMZ]  = rho * uz;
    cons[UEINT] = total_energy;
    for (int n = 0; n < thermo.num_species; ++n) {
        cons[UEINT + 1 + n] = rho * species_massfractions[n];
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
compute_sos(
        amrex::GpuArray<amrex::Real, NPRIM> const&  prims,
        pyro::pyro<double> const& thermo) noexcept
{
    using amrex::Real;
    Real rho = prims[QRHO];
    Real ux = prims[QU];
    Real uy = prims[QV];
    Real uz = prims[QW];
    Real temperature = prims[QTEMP];
    std::array<double, NSP> species_massfractions = {0.0};
    for (int n = 0; n < thermo.num_species; ++n) {
        species_massfractions[n] = prims[QTEMP + 1 + n];
    }
    Real c = sqrt(thermo.get_mixture_specific_heat_cp_mass(temperature, species_massfractions) *
                  thermo.get_specific_gas_constant(species_massfractions) /
                  thermo.get_mixture_specific_heat_cv_mass(temperature, species_massfractions) *
                  temperature);
    return c;
}


#endif
