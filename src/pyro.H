#pragma once

#if __cplusplus < 201103L
#error "Pyrometheus requires C++11 or later."
#endif

#include <array>
#include <cmath>
#include <string>

namespace pyro {

template <typename _DataTypeT = double, typename _ContainerT = _DataTypeT>
struct pyro
{
    constexpr static int num_species = 2;
    constexpr static int num_reactions = 1;
    constexpr static int num_falloff = 0;

    static const char* species_names[];
    static const char* element_names[];

    using DataTypeT  = _DataTypeT;
    using ContainerT = _ContainerT;
    using SpeciesT   = std::array<ContainerT, num_species>;
    using Species2T  = std::array<SpeciesT, num_species>;
    using ReactionsT = std::array<ContainerT, num_reactions>;

    static const DataTypeT molecular_weights[];
    static const DataTypeT inv_molecular_weights[];
    static const DataTypeT gas_constant;
    static const DataTypeT one_atm;

    static int nSpecies() { return num_species; };
    
    // Helper function to convert double literals to template type
    static constexpr ContainerT to_type(double value) {
        return static_cast<ContainerT>(value);
    }

    static std::string get_species_name(int index);
    static std::string get_element_name(int index);

    static int get_species_index(const std::string& name);
    static int get_element_index(const std::string& name);

    // Thermodynamic properties
    static ContainerT get_specific_gas_constant(SpeciesT const &mass_fractions);
    static ContainerT get_mix_molecular_weight(SpeciesT const &mass_fractions);
    static ContainerT get_mixture_molecular_weight(SpeciesT mass_fractions);
    
    static SpeciesT get_concentrations(ContainerT rho, SpeciesT const &mass_fractions);
    static SpeciesT get_mole_fractions(ContainerT mix_mol_weight, SpeciesT mass_fractions);
    
    static ContainerT get_mass_average_property(SpeciesT const &mass_fractions, SpeciesT const &spec_property);
    
    // Energy and thermodynamic state
    static ContainerT get_density(ContainerT p, ContainerT temperature, SpeciesT const &mass_fractions);
    static ContainerT get_pressure(ContainerT density, ContainerT temperature, SpeciesT const &mass_fractions);
    
    // Specific heats
    static ContainerT get_mixture_specific_heat_cp_mass(ContainerT temperature, SpeciesT const &mass_fractions);
    static ContainerT get_mixture_specific_heat_cv_mass(ContainerT temperature, SpeciesT const &mass_fractions);
    
    // Energy properties
    static ContainerT get_mixture_enthalpy_mass(ContainerT temperature, SpeciesT const &mass_fractions);
    static ContainerT get_mixture_internal_energy_mass(ContainerT temperature, SpeciesT const &mass_fractions);
    
    // Species-specific thermodynamic properties
    static SpeciesT get_species_specific_heats_r(ContainerT temperature);
    static SpeciesT get_species_enthalpies_rt(ContainerT temperature);
    static SpeciesT get_species_entropies_r(ContainerT temperature);
    static SpeciesT get_species_gibbs_rt(ContainerT temperature);
    
    // Chemical kinetics
    static ReactionsT get_equilibrium_constants(ContainerT temperature);
    static ReactionsT get_fwd_rate_coefficients(ContainerT temperature, SpeciesT const &concentrations);
    static ReactionsT get_net_rates_of_progress(ContainerT temperature, SpeciesT const &concentrations);
    static SpeciesT get_net_production_rates(ContainerT rho, ContainerT temperature, SpeciesT const &mass_fractions);
    
    // Temperature calculation
    static ContainerT get_temperature(ContainerT energy_or_enthalpy, ContainerT t_guess, 
                                     SpeciesT const &mass_fractions, bool const do_energy = true);
    
    // Transport properties
    static SpeciesT get_species_viscosities(ContainerT temperature);
    static SpeciesT get_species_thermal_conductivities(ContainerT temperature);
    static Species2T get_species_binary_mass_diffusivities(ContainerT temperature);
    
    // Mixture transport properties
    static ContainerT get_mixture_viscosity_mixavg(ContainerT temperature, SpeciesT mass_fractions);
    static ContainerT get_mixture_thermal_conductivity_mixavg(ContainerT temperature, SpeciesT mass_fractions);
    static SpeciesT get_species_mass_diffusivities_mixavg(ContainerT pressure, ContainerT temperature, SpeciesT mass_fractions);
};

} // namespace pyro
