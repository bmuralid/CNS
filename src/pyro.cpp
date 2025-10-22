#include "pyro.H"

namespace pyro {

template <typename _DataTypeT, typename _ContainerT>
std::string pyro<_DataTypeT, _ContainerT>::get_species_name(int index) {
    return species_names[index];
}

template <typename _DataTypeT, typename _ContainerT>
std::string pyro<_DataTypeT, _ContainerT>::get_element_name(int index) {
    return element_names[index];
}

template <typename _DataTypeT, typename _ContainerT>
int pyro<_DataTypeT, _ContainerT>::get_species_index(const std::string& name) {
    if (name == "O2") return 0;
    if (name == "N2") return 1;
    return -1;
}

template <typename _DataTypeT, typename _ContainerT>
int pyro<_DataTypeT, _ContainerT>::get_element_index(const std::string& name) {
    if (name == "O") return 0;
    if (name == "N") return 1;
    return -1;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ContainerT pyro<_DataTypeT, _ContainerT>::get_specific_gas_constant(SpeciesT const &mass_fractions)
{
    return gas_constant * (
                + inv_molecular_weights[0]*mass_fractions[0]
                + inv_molecular_weights[1]*mass_fractions[1]
            );
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ContainerT pyro<_DataTypeT, _ContainerT>::get_mix_molecular_weight(SpeciesT const &mass_fractions)
{
    return 1.0/(
        + inv_molecular_weights[0]*mass_fractions[0]
        + inv_molecular_weights[1]*mass_fractions[1]
    );
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::SpeciesT pyro<_DataTypeT, _ContainerT>::get_concentrations(
    ContainerT rho, SpeciesT const &mass_fractions)
{
    SpeciesT concentrations = {
            inv_molecular_weights[0]*mass_fractions[0]*rho,
            inv_molecular_weights[1]*mass_fractions[1]*rho,
    };
    return concentrations;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::SpeciesT pyro<_DataTypeT, _ContainerT>::get_mole_fractions(
    ContainerT mix_mol_weight, SpeciesT mass_fractions)
{
    return SpeciesT{
        inv_molecular_weights[0] * mass_fractions[0] * mix_mol_weight,
        inv_molecular_weights[1] * mass_fractions[1] * mix_mol_weight,
    };
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ContainerT pyro<_DataTypeT, _ContainerT>::get_mass_average_property(
    SpeciesT const &mass_fractions, SpeciesT const &spec_property)
{
    return (
            + mass_fractions[0]*
              spec_property[0]*
              inv_molecular_weights[0]
            + mass_fractions[1]*
              spec_property[1]*
              inv_molecular_weights[1]
    );
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ContainerT pyro<_DataTypeT, _ContainerT>::get_mixture_specific_heat_cv_mass(
    ContainerT temperature, SpeciesT const &mass_fractions)
{
    SpeciesT cp0_r = get_species_specific_heats_r(temperature);

    for (int i = 0; i < num_species; ++i)
        cp0_r[i] -= 1.0;

    const ContainerT cpmix = get_mass_average_property(mass_fractions, cp0_r);

    return gas_constant * cpmix;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ContainerT pyro<_DataTypeT, _ContainerT>::get_mixture_specific_heat_cp_mass(
    ContainerT temperature, SpeciesT const &mass_fractions)
{
    const SpeciesT cp0_r = get_species_specific_heats_r(temperature);
    const ContainerT cpmix = get_mass_average_property(mass_fractions, cp0_r);
    return gas_constant * cpmix;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ContainerT pyro<_DataTypeT, _ContainerT>::get_mixture_enthalpy_mass(
    ContainerT temperature, SpeciesT const &mass_fractions)
{
    SpeciesT h0_rt = get_species_enthalpies_rt(temperature);
    ContainerT hmix = get_mass_average_property(mass_fractions, h0_rt);
    return gas_constant * hmix * temperature;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ContainerT pyro<_DataTypeT, _ContainerT>::get_mixture_internal_energy_mass(
    ContainerT temperature, SpeciesT const &mass_fractions)
{
    SpeciesT e0_rt = get_species_enthalpies_rt(temperature);

    for (int i = 0; i < num_species; ++i)
        e0_rt[i] -= 1.0;

    const ContainerT emix = get_mass_average_property(mass_fractions, e0_rt);
    return gas_constant * emix * temperature;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ContainerT pyro<_DataTypeT, _ContainerT>::get_density(ContainerT p, ContainerT temperature,
        SpeciesT const &mass_fractions)
{
    ContainerT mmw = get_mix_molecular_weight(mass_fractions);
    ContainerT rt = gas_constant * temperature;
    return p * mmw / rt;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ContainerT pyro<_DataTypeT, _ContainerT>::get_pressure(
    ContainerT density, ContainerT temperature, SpeciesT const &mass_fractions)
{
    const double mmw = get_mix_molecular_weight(mass_fractions);
    const double rt = gas_constant * temperature;
    return density * rt / mmw;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ContainerT pyro<_DataTypeT, _ContainerT>::get_mixture_molecular_weight(SpeciesT mass_fractions) {
    return 1.0/(
        + inv_molecular_weights[0]*mass_fractions[0]
        + inv_molecular_weights[1]*mass_fractions[1]
    );
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::SpeciesT pyro<_DataTypeT, _ContainerT>::get_species_specific_heats_r(ContainerT temperature)
{
    SpeciesT cp0_r = {
        (temperature > to_type(1000.0) ? to_type(3.28253784) + to_type(0.00148308754) * temperature + to_type(2.09470555e-10) * static_cast<ContainerT>(pow(temperature, 3)) + to_type(-7.57966669e-07) * temperature * temperature + to_type(-2.16717794e-14) * static_cast<ContainerT>(pow(temperature, 4)) : to_type(3.78245636) + to_type(9.84730201e-06) * temperature * temperature + to_type(3.24372837e-12) * static_cast<ContainerT>(pow(temperature, 4)) + to_type(-9.68129509e-09) * static_cast<ContainerT>(pow(temperature, 3)) + to_type(-0.00299673416) * temperature),
        (temperature > to_type(1000.0) ? to_type(2.92664) + to_type(0.0014879768) * temperature + to_type(1.0097038e-10) * static_cast<ContainerT>(pow(temperature, 3)) + to_type(-6.753351e-15) * static_cast<ContainerT>(pow(temperature, 4)) + to_type(-5.68476e-07) * temperature * temperature : to_type(3.298677) + to_type(0.0014082404) * temperature + to_type(5.641515e-09) * static_cast<ContainerT>(pow(temperature, 3)) + to_type(-3.963222e-06) * temperature * temperature + to_type(-2.444854e-12) * static_cast<ContainerT>(pow(temperature, 4))),
        };
    return cp0_r;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::SpeciesT pyro<_DataTypeT, _ContainerT>::get_species_enthalpies_rt(ContainerT temperature)
{
    SpeciesT h0_rt = {
        (temperature > to_type(1000.0) ? to_type(3.28253784) + to_type(0.00074154377) * temperature + to_type(5.236763875e-11) * static_cast<ContainerT>(pow(temperature, 3)) + to_type(-4.33435588e-15) * static_cast<ContainerT>(pow(temperature, 4)) + to_type(-2.526555563333333e-07) * temperature * temperature + to_type(-1088.45772) / temperature : to_type(3.78245636) + to_type(6.48745674e-13) * static_cast<ContainerT>(pow(temperature, 4)) + to_type(3.282434003333333e-06) * temperature * temperature + to_type(-2.4203237725e-09) * static_cast<ContainerT>(pow(temperature, 3)) + to_type(-1063.94356) / temperature + to_type(-0.00149836708) * temperature),
        (temperature > to_type(1000.0) ? to_type(2.92664) + to_type(0.0007439884) * temperature + to_type(2.5242595e-11) * static_cast<ContainerT>(pow(temperature, 3)) + to_type(-922.7977) / temperature + to_type(-1.8949200000000001e-07) * temperature * temperature + to_type(-1.3506701999999999e-15) * static_cast<ContainerT>(pow(temperature, 4)) : to_type(3.298677) + to_type(0.0007041202) * temperature + to_type(1.41037875e-09) * static_cast<ContainerT>(pow(temperature, 3)) + to_type(-4.889707999999999e-13) * static_cast<ContainerT>(pow(temperature, 4)) + to_type(-1020.8999) / temperature + to_type(-1.3210739999999999e-06) * temperature * temperature),
        };
    return h0_rt;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::SpeciesT pyro<_DataTypeT, _ContainerT>::get_species_entropies_r(ContainerT temperature)
{
    SpeciesT s0_r = {
        (temperature > to_type(1000.0) ? to_type(5.45323129) + to_type(3.28253784) * static_cast<ContainerT>(log(temperature)) + to_type(0.00148308754) * temperature + to_type(6.982351833333333e-11) * static_cast<ContainerT>(pow(temperature, 3)) + to_type(-5.41794485e-15) * static_cast<ContainerT>(pow(temperature, 4)) + to_type(-3.789833345e-07) * temperature * temperature : to_type(3.78245636) * static_cast<ContainerT>(log(temperature)) + to_type(3.65767573) + to_type(8.109320925e-13) * static_cast<ContainerT>(pow(temperature, 4)) + to_type(4.923651005e-06) * temperature * temperature + to_type(-3.2270983633333334e-09) * static_cast<ContainerT>(pow(temperature, 3)) + to_type(-0.00299673416) * temperature),
        (temperature > to_type(1000.0) ? to_type(5.980528) + to_type(2.92664) * static_cast<ContainerT>(log(temperature)) + to_type(0.0014879768) * temperature + to_type(3.3656793333333334e-11) * static_cast<ContainerT>(pow(temperature, 3)) + to_type(-2.84238e-07) * temperature * temperature + to_type(-1.68833775e-15) * static_cast<ContainerT>(pow(temperature, 4)) : to_type(3.950372) + to_type(3.298677) * static_cast<ContainerT>(log(temperature)) + to_type(0.0014082404) * temperature + to_type(1.8805050000000002e-09) * static_cast<ContainerT>(pow(temperature, 3)) + to_type(-6.112135e-13) * static_cast<ContainerT>(pow(temperature, 4)) + to_type(-1.981611e-06) * temperature * temperature),
        };
    return s0_r;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::SpeciesT pyro<_DataTypeT, _ContainerT>::get_species_gibbs_rt(ContainerT temperature)
{
    SpeciesT h0_rt = get_species_enthalpies_rt(temperature);
    SpeciesT s0_r = get_species_entropies_r(temperature);
    SpeciesT g0_rt = {
    h0_rt[0] - s0_r[0],
    h0_rt[1] - s0_r[1],
    };
    return g0_rt;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ReactionsT pyro<_DataTypeT, _ContainerT>::get_equilibrium_constants(ContainerT temperature)
{
    ContainerT rt = gas_constant * temperature;
    ContainerT c0 = static_cast<ContainerT>(std::log(one_atm/rt));
    SpeciesT g0_rt = get_species_gibbs_rt(temperature);
    ReactionsT k_eq = {
    to_type(1.0) * g0_rt[1] + to_type(1.0) * g0_rt[0] + to_type(0) - (to_type(1.0) * g0_rt[1] + to_type(1.0) * g0_rt[0] + to_type(0)),
    };
    return k_eq;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ContainerT pyro<_DataTypeT, _ContainerT>::get_temperature(
    ContainerT energy_or_enthalpy,
    ContainerT t_guess,
    SpeciesT const &mass_fractions,
    bool const do_energy)
{
    ContainerT (*pv_fun)(ContainerT, SpeciesT const &);
    ContainerT (*he_fun)(ContainerT, SpeciesT const &);

    if (do_energy) {
        pv_fun = get_mixture_specific_heat_cv_mass;
        he_fun = get_mixture_internal_energy_mass;
    } else {
        pv_fun = get_mixture_specific_heat_cp_mass;
        he_fun = get_mixture_enthalpy_mass;
    }

    int num_iter = 500;
    ContainerT tol = to_type(1.0e-06);
    ContainerT iter_temp = t_guess;

    for(int iter = 0; iter < num_iter; ++iter){
        ContainerT iter_rhs   = energy_or_enthalpy - he_fun(iter_temp, mass_fractions);
        ContainerT iter_deriv = -pv_fun(iter_temp, mass_fractions);
        ContainerT iter_dt    = -iter_rhs/iter_deriv;
        iter_temp += iter_dt;
        if(std::fabs(iter_dt) < tol){ break; }
    }
    return iter_temp;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ReactionsT pyro<_DataTypeT, _ContainerT>::get_fwd_rate_coefficients(ContainerT temperature,
                                            SpeciesT const &concentrations)
{
    ReactionsT k_fwd = {
    static_cast<ContainerT>(exp(to_type(24.01910270295074) + to_type(0.0) * static_cast<ContainerT>(log(temperature)) - to_type(178.64293439206185) / temperature)),
    };

    return k_fwd;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ReactionsT pyro<_DataTypeT, _ContainerT>::get_net_rates_of_progress(
    ContainerT temperature, SpeciesT const &concentrations)
{
    ReactionsT k_fwd = get_fwd_rate_coefficients(temperature, concentrations);
    ReactionsT log_k_eq = get_equilibrium_constants(temperature);
    ReactionsT r_net = {
    k_fwd[0] * (concentrations[1] * concentrations[0] - static_cast<ContainerT>(exp(log_k_eq[0])) * concentrations[1] * concentrations[0]),
    };
    return r_net;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::SpeciesT pyro<_DataTypeT, _ContainerT>::get_net_production_rates(
    ContainerT rho, ContainerT temperature, SpeciesT const &mass_fractions)
{
    SpeciesT concentrations = get_concentrations(rho, mass_fractions);
    ReactionsT r_net = get_net_rates_of_progress(temperature, concentrations);
    SpeciesT omega = {
    (to_type(1.0) * r_net[0] + to_type(0) - (to_type(1.0) * r_net[0] + to_type(0))) * (to_type(1.0) + to_type(0) * r_net[0]),
    (to_type(1.0) * r_net[0] + to_type(0) - (to_type(1.0) * r_net[0] + to_type(0))) * (to_type(1.0) + to_type(0) * r_net[0]),
    };
    return omega;
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::SpeciesT pyro<_DataTypeT, _ContainerT>::get_species_viscosities(ContainerT temperature)
{
    ContainerT c1 = ContainerT(0.0036188245122444136);
    ContainerT c2 = ContainerT(-0.00618642807174783);
    ContainerT c3 = ContainerT(5.9160127096038374e-05);
    ContainerT c4 = ContainerT(-1.9049771784292983e-06);
    ContainerT c5 = ContainerT(-0.0006861983404476142);

    ContainerT c6 = ContainerT(0.0031249811202346835);
    ContainerT c7 = ContainerT(-0.0052325034588374736);
    ContainerT c8 = ContainerT(5.1908396950248406e-05);
    ContainerT c9 = ContainerT(-1.6850989392683053e-06);
    ContainerT c10 = ContainerT(-0.0005968857242770352);

    return SpeciesT{
        static_cast<ContainerT>(sqrt(temperature)) * (c1 * static_cast<ContainerT>(log(temperature)) + c2 + c3 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 3)) + c4 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 4)) + c5 * static_cast<ContainerT>(log(temperature)) * static_cast<ContainerT>(log(temperature))) * (c1 * static_cast<ContainerT>(log(temperature)) + c2 + c3 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 3)) + c4 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 4)) + c5 * static_cast<ContainerT>(log(temperature)) * static_cast<ContainerT>(log(temperature))),
        static_cast<ContainerT>(sqrt(temperature)) * (c6 * static_cast<ContainerT>(log(temperature)) + c7 + c8 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 3)) + c9 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 4)) + c10 * static_cast<ContainerT>(log(temperature)) * static_cast<ContainerT>(log(temperature))) * (c6 * static_cast<ContainerT>(log(temperature)) + c7 + c8 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 3)) + c9 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 4)) + c10 * static_cast<ContainerT>(log(temperature)) * static_cast<ContainerT>(log(temperature))),
    };
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::SpeciesT pyro<_DataTypeT, _ContainerT>::get_species_thermal_conductivities(ContainerT temperature)
{
    ContainerT c1 = ContainerT(0.10689552101630528);
    ContainerT c2 = ContainerT(0.014217756178712414);
    ContainerT c3 = ContainerT(5.09174438888175e-05);
    ContainerT c4 = ContainerT(-0.06376711343770658);
    ContainerT c5 = ContainerT(-0.0013908411531386269);

    ContainerT c6 = ContainerT(0.0026129214099176513);
    ContainerT c7 = ContainerT(0.0015932386446981589);
    ContainerT c8 = ContainerT(0.00016507154037197712);
    ContainerT c9 = ContainerT(-8.297315315373045e-06);
    ContainerT c10 = ContainerT(-0.000984277527739843);

    return SpeciesT{
        static_cast<ContainerT>(sqrt(temperature)) * (c1 + c2 * static_cast<ContainerT>(log(temperature)) * static_cast<ContainerT>(log(temperature)) + c3 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 4)) + c4 * static_cast<ContainerT>(log(temperature)) + c5 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 3))),
        static_cast<ContainerT>(sqrt(temperature)) * (c6 + c7 * static_cast<ContainerT>(log(temperature)) + c8 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 3)) + c9 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 4)) + c10 * static_cast<ContainerT>(log(temperature)) * static_cast<ContainerT>(log(temperature))),
    };
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::Species2T pyro<_DataTypeT, _ContainerT>::get_species_binary_mass_diffusivities(ContainerT temperature)
{
    ContainerT c1 = ContainerT(0.0016332321885344269);
    ContainerT c2 = ContainerT(-0.0031272899371257495);
    ContainerT c3 = ContainerT(2.3795154191077426e-05);
    ContainerT c4 = ContainerT(-7.135754459127752e-07);
    ContainerT c5 = ContainerT(-0.00029023244731892585);

    ContainerT c6 = ContainerT(0.0015963119455941247);
    ContainerT c7 = ContainerT(-0.00303005347627421);
    ContainerT c8 = ContainerT(2.3589003022566823e-05);
    ContainerT c9 = ContainerT(-7.128139183150021e-07);
    ContainerT c10 = ContainerT(-0.00028558271238060057);

    ContainerT c11 = ContainerT(0.001572491741915425);
    ContainerT c12 = ContainerT(-0.00295673852294946);
    ContainerT c13 = ContainerT(2.364377574630533e-05);
    ContainerT c14 = ContainerT(-7.213510385775298e-07);
    ContainerT c15 = ContainerT(-0.000283699905789111);

    return Species2T{
        SpeciesT{
            static_cast<ContainerT>(sqrt(temperature)) * temperature * (c1 * static_cast<ContainerT>(log(temperature)) + c2 + c3 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 3)) + c4 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 4)) + c5 * static_cast<ContainerT>(log(temperature)) * static_cast<ContainerT>(log(temperature))),
            static_cast<ContainerT>(sqrt(temperature)) * temperature * (c6 * static_cast<ContainerT>(log(temperature)) + c7 + c8 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 3)) + c9 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 4)) + c10 * static_cast<ContainerT>(log(temperature)) * static_cast<ContainerT>(log(temperature))),
        },
        SpeciesT{
            static_cast<ContainerT>(sqrt(temperature)) * temperature * (c6 * static_cast<ContainerT>(log(temperature)) + c7 + c8 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 3)) + c9 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 4)) + c10 * static_cast<ContainerT>(log(temperature)) * static_cast<ContainerT>(log(temperature))),
            static_cast<ContainerT>(sqrt(temperature)) * temperature * (c11 * static_cast<ContainerT>(log(temperature)) + c12 + c13 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 3)) + c14 * static_cast<ContainerT>(pow(static_cast<ContainerT>(log(temperature)), 4)) + c15 * static_cast<ContainerT>(log(temperature)) * static_cast<ContainerT>(log(temperature))),
        },
    };
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ContainerT pyro<_DataTypeT, _ContainerT>::get_mixture_viscosity_mixavg(
    ContainerT temperature, SpeciesT mass_fractions)
{
    ContainerT mmw = get_mixture_molecular_weight(mass_fractions);
    SpeciesT mole_fractions = get_mole_fractions(mmw, mass_fractions);
    SpeciesT viscosities = get_species_viscosities(temperature);
    SpeciesT mix_rule_f = {
        to_type(0) + (mole_fractions[1] * (static_cast<ContainerT>(sqrt(viscosities[0] / viscosities[1] * static_cast<ContainerT>(sqrt(to_type(0.8754922182636414))))) + to_type(1)) * (static_cast<ContainerT>(sqrt(viscosities[0] / viscosities[1] * static_cast<ContainerT>(sqrt(to_type(0.8754922182636414))))) + to_type(1))) / static_cast<ContainerT>(sqrt(to_type(17.137716855857786))) + (mole_fractions[0] * (static_cast<ContainerT>(sqrt(viscosities[0] / viscosities[0] * static_cast<ContainerT>(sqrt(to_type(1.0))))) + to_type(1)) * (static_cast<ContainerT>(sqrt(viscosities[0] / viscosities[0] * static_cast<ContainerT>(sqrt(to_type(1.0))))) + to_type(1))) / static_cast<ContainerT>(sqrt(to_type(16.0))),
        to_type(0) + (mole_fractions[1] * (static_cast<ContainerT>(sqrt(viscosities[1] / viscosities[1] * static_cast<ContainerT>(sqrt(to_type(1.0))))) + to_type(1)) * (static_cast<ContainerT>(sqrt(viscosities[1] / viscosities[1] * static_cast<ContainerT>(sqrt(to_type(1.0))))) + to_type(1))) / static_cast<ContainerT>(sqrt(to_type(16.0))) + (mole_fractions[0] * (static_cast<ContainerT>(sqrt(viscosities[1] / viscosities[0] * static_cast<ContainerT>(sqrt(to_type(1.1422146069822232))))) + to_type(1)) * (static_cast<ContainerT>(sqrt(viscosities[1] / viscosities[0] * static_cast<ContainerT>(sqrt(to_type(1.1422146069822232))))) + to_type(1))) / static_cast<ContainerT>(sqrt(to_type(15.00393774610913))),
    };

    return (to_type(0.0) +
        + mole_fractions[0]*viscosities[0]/mix_rule_f[0]
        + mole_fractions[1]*viscosities[1]/mix_rule_f[1]
    );
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::ContainerT pyro<_DataTypeT, _ContainerT>::get_mixture_thermal_conductivity_mixavg(
    ContainerT temperature, SpeciesT mass_fractions)
{
    ContainerT mmw = get_mixture_molecular_weight(mass_fractions);
    SpeciesT mole_fractions = get_mole_fractions(mmw, mass_fractions);
    SpeciesT conductivities = get_species_thermal_conductivities(temperature);

    ContainerT lhs = (to_type(0.0) +
        + mole_fractions[0]*conductivities[0]
        + mole_fractions[1]*conductivities[1]
    );

    ContainerT rhs = to_type(1.0) / (to_type(0.0) +
        + mole_fractions[0]/conductivities[0]
        + mole_fractions[1]/conductivities[1]
    );

    return to_type(0.5)*(lhs + rhs);
}

template <typename _DataTypeT, typename _ContainerT>
typename pyro<_DataTypeT, _ContainerT>::SpeciesT pyro<_DataTypeT, _ContainerT>::get_species_mass_diffusivities_mixavg(
    ContainerT pressure, ContainerT temperature, SpeciesT mass_fractions)
{
    ContainerT mmw = get_mixture_molecular_weight(mass_fractions);
    SpeciesT mole_fractions = get_mole_fractions(mmw, mass_fractions);
    Species2T bdiff_ij = get_species_binary_mass_diffusivities(temperature);

    SpeciesT x_sum = {
        mole_fractions[1] / bdiff_ij[1][0] + mole_fractions[0] / bdiff_ij[0][0] + 0,
        mole_fractions[1] / bdiff_ij[1][1] + mole_fractions[0] / bdiff_ij[0][1] + 0,
    };

    SpeciesT denom = {
        x_sum[0] - mole_fractions[0]/bdiff_ij[0][0],
        x_sum[1] - mole_fractions[1]/bdiff_ij[1][1],
    };

    return SpeciesT{
            denom[0] > 0.0 ?
                (mmw - mole_fractions[0] * molecular_weights[0])
                  / (pressure * mmw * denom[0])
              : (bdiff_ij[0][0] / pressure),
            denom[1] > 0.0 ?
                (mmw - mole_fractions[1] * molecular_weights[1])
                  / (pressure * mmw * denom[1])
              : (bdiff_ij[1][1] / pressure),
    };
}

// Static member definitions
template <typename _DataTypeT, typename _ContainerT>
const char* pyro<_DataTypeT, _ContainerT>::species_names[] = {
    "O2", "N2"
};

template <typename _DataTypeT, typename _ContainerT>
const char* pyro<_DataTypeT, _ContainerT>::element_names[] = {
    "O", "N"
};

template <typename _DataTypeT, typename _ContainerT>
const _DataTypeT pyro<_DataTypeT, _ContainerT>::molecular_weights[] = {
    31.998, 28.014
};

template <typename _DataTypeT, typename _ContainerT>
const _DataTypeT pyro<_DataTypeT, _ContainerT>::inv_molecular_weights[] = {
    0.03125195324707794, 0.03569643749553795
};

template <typename _DataTypeT, typename _ContainerT>
const _DataTypeT pyro<_DataTypeT, _ContainerT>::gas_constant = 8314.46261815324;

template <typename _DataTypeT, typename _ContainerT>
const _DataTypeT pyro<_DataTypeT, _ContainerT>::one_atm = 101325.0;

// Explicit instantiations for common types
template struct pyro<double, double>;
template struct pyro<float, float>;

} // namespace pyro
