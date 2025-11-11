
#pragma once

#if __cplusplus < 202002L
#error "Pyrometheus requires C++20 or later."
#endif

#include <array>
#include <cmath>
#include <string>
#include <AMReX_Gpu.H>
#if __cplusplus >= 202002L
#include <concepts>
#endif

namespace pyro {

template <typename T>
concept IntegralLike = requires(T a, T b) {
    { T(1.0)  } -> std::same_as<T>;    { a         } -> std::convertible_to<double>;
    { a + b   } -> std::same_as<T>;    { a - b     } -> std::same_as<T>;
    { a * b   } -> std::same_as<T>;    { a / b     } -> std::same_as<T>;
    { -a      } -> std::same_as<T>;    { +a        } -> std::same_as<T>;
    { a += b  } -> std::same_as<T&>;   { a -= b    } -> std::same_as<T&>;
    { a *= b  } -> std::same_as<T&>;   { a /= b    } -> std::same_as<T&>;
    { a == b  } -> std::same_as<bool>; { a != b    } -> std::same_as<bool>;
    { a < b   } -> std::same_as<bool>; { a > b     } -> std::same_as<bool>;
    { a <= b  } -> std::same_as<bool>; { a >= b    } -> std::same_as<bool>;
    { log(a)  } -> std::same_as<T>;    { exp(a)    } -> std::same_as<T>;
    { sqrt(a) } -> std::same_as<T>;    { pow(a, b) } -> std::same_as<T>;
};

template <IntegralLike _DataTypeT = double, typename _ContainerT = _DataTypeT>
struct pyro
{
    constexpr static int num_species = 2;
    constexpr static int num_reactions = 1;
    constexpr static int num_falloff = 0;

    constexpr static const char* species_names[] = {
        "O2", "N2"
    };

    constexpr static const char* element_names[] = {
        "O", "N"
    };

    using DataTypeT  = _DataTypeT;
    using ContainerT = _ContainerT;
    using SpeciesT   = std::array<ContainerT, num_species>;
    using Species2T  = std::array<SpeciesT, num_species>;
    using ReactionsT = std::array<ContainerT, num_reactions>;

    constexpr static DataTypeT molecular_weights[] =
        {31.998, 28.014};
    constexpr static DataTypeT inv_molecular_weights[] =
        {0.03125195324707794, 0.03569643749553795};
    constexpr static DataTypeT gas_constant = 8314.46261815324;
    constexpr static DataTypeT one_atm = 101325.0;

    static int nSpecies() { return num_species; };

    static std::string get_species_name(int index) { return species_names[index]; }
    static std::string get_element_name(int index) { return element_names[index]; }

    static int get_species_index(const std::string& name) {
            if (name == "O2") return 0;
            if (name == "N2") return 1;

        return -1;
    }

    static int get_element_index(const std::string& name) {
            if (name == "O") return 0;
            if (name == "N") return 1;

        return -1;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ContainerT get_specific_gas_constant(SpeciesT const &mass_fractions)
    {
        return gas_constant * (
                    + inv_molecular_weights[0]*mass_fractions[0]
                    + inv_molecular_weights[1]*mass_fractions[1]
                );
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ContainerT get_mix_molecular_weight(SpeciesT const &mass_fractions)
    {
        return 1.0/(
            + inv_molecular_weights[0]*mass_fractions[0]
            + inv_molecular_weights[1]*mass_fractions[1]
        );
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static SpeciesT get_concentrations(
        ContainerT rho, SpeciesT const &mass_fractions)
    {
        SpeciesT concentrations = {
                inv_molecular_weights[0]*mass_fractions[0]*rho,
                inv_molecular_weights[1]*mass_fractions[1]*rho,
        };
        return concentrations;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static SpeciesT get_mole_fractions(
        ContainerT mix_mol_weight, SpeciesT mass_fractions)
    {
        return SpeciesT{
            inv_molecular_weights[0] * mass_fractions[0] * mix_mol_weight,
            inv_molecular_weights[1] * mass_fractions[1] * mix_mol_weight,
        };
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ContainerT get_mass_average_property(
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

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ContainerT get_mixture_specific_heat_cv_mass(
        ContainerT temperature, SpeciesT const &mass_fractions)
    {
        SpeciesT cp0_r = get_species_specific_heats_r(temperature);

        for (int i = 0; i < num_species; ++i)
            cp0_r[i] -= 1.0;

        const ContainerT cpmix = get_mass_average_property(mass_fractions, cp0_r);

        return gas_constant * cpmix;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ContainerT get_mixture_specific_heat_cp_mass(
        ContainerT temperature, SpeciesT const &mass_fractions)
    {
        const SpeciesT cp0_r = get_species_specific_heats_r(temperature);
        const ContainerT cpmix = get_mass_average_property(mass_fractions, cp0_r);
        return gas_constant * cpmix;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ContainerT get_mixture_enthalpy_mass(
        ContainerT temperature, SpeciesT const &mass_fractions)
    {
        SpeciesT h0_rt = get_species_enthalpies_rt(temperature);
        ContainerT hmix = get_mass_average_property(mass_fractions, h0_rt);
        return gas_constant * hmix * temperature;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ContainerT get_mixture_internal_energy_mass(
        ContainerT temperature, SpeciesT const &mass_fractions)
    {
        SpeciesT e0_rt = get_species_enthalpies_rt(temperature);

        for (int i = 0; i < num_species; ++i)
            e0_rt[i] -= 1.0;

        const ContainerT emix = get_mass_average_property(mass_fractions, e0_rt);
        return gas_constant * emix * temperature;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ContainerT get_density(ContainerT p, ContainerT temperature,
            SpeciesT const &mass_fractions)
    {
        ContainerT mmw = get_mix_molecular_weight(mass_fractions);
        ContainerT rt = gas_constant * temperature;
        return p * mmw / rt;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ContainerT get_pressure(
        ContainerT density, ContainerT temperature, SpeciesT const &mass_fractions)
    {
        const double mmw = get_mix_molecular_weight(mass_fractions);
        const double rt = gas_constant * temperature;
        return density * rt / mmw;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ContainerT get_mixture_molecular_weight(SpeciesT mass_fractions) {
        return 1.0/(
            + inv_molecular_weights[0]*mass_fractions[0]
            + inv_molecular_weights[1]*mass_fractions[1]
        );
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static SpeciesT get_species_specific_heats_r(ContainerT temperature)
    {
        SpeciesT cp0_r = {
            (temperature > 1000.0 ? 3.28253784 + 0.00148308754 * temperature + (2.09470555e-10) * pow(temperature, 3) + (-7.57966669e-07) * temperature * temperature + (-2.16717794e-14) * pow(temperature, 4) : 3.78245636 + (9.84730201e-06) * temperature * temperature + (3.24372837e-12) * pow(temperature, 4) + (-9.68129509e-09) * pow(temperature, 3) + (-0.00299673416) * temperature),
            (temperature > 1000.0 ? 2.92664 + 0.0014879768 * temperature + (1.0097038e-10) * pow(temperature, 3) + (-6.753351e-15) * pow(temperature, 4) + (-5.68476e-07) * temperature * temperature : 3.298677 + 0.0014082404 * temperature + (5.641515e-09) * pow(temperature, 3) + (-3.963222e-06) * temperature * temperature + (-2.444854e-12) * pow(temperature, 4)),
            };
        return cp0_r;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static SpeciesT get_species_enthalpies_rt(ContainerT temperature)
    {
        SpeciesT h0_rt = {
            (temperature > 1000.0 ? 3.28253784 + 0.00074154377 * temperature + (5.236763875e-11) * pow(temperature, 3) + (-4.33435588e-15) * pow(temperature, 4) + (-2.526555563333333e-07) * temperature * temperature + (-1088.45772) / temperature : 3.78245636 + (6.48745674e-13) * pow(temperature, 4) + (3.282434003333333e-06) * temperature * temperature + (-2.4203237725e-09) * pow(temperature, 3) + (-1063.94356) / temperature + (-0.00149836708) * temperature),
            (temperature > 1000.0 ? 2.92664 + 0.0007439884 * temperature + (2.5242595e-11) * pow(temperature, 3) + (-922.7977) / temperature + (-1.8949200000000001e-07) * temperature * temperature + (-1.3506701999999999e-15) * pow(temperature, 4) : 3.298677 + 0.0007041202 * temperature + (1.41037875e-09) * pow(temperature, 3) + (-4.889707999999999e-13) * pow(temperature, 4) + (-1020.8999) / temperature + (-1.3210739999999999e-06) * temperature * temperature),
            };
        return h0_rt;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static SpeciesT get_species_entropies_r(ContainerT temperature)
    {
        SpeciesT s0_r = {
            (temperature > 1000.0 ? 5.45323129 + 3.28253784 * log(temperature) + 0.00148308754 * temperature + (6.982351833333333e-11) * pow(temperature, 3) + (-5.41794485e-15) * pow(temperature, 4) + (-3.789833345e-07) * temperature * temperature : 3.78245636 * log(temperature) + 3.65767573 + (8.109320925e-13) * pow(temperature, 4) + (4.923651005e-06) * temperature * temperature + (-3.2270983633333334e-09) * pow(temperature, 3) + (-0.00299673416) * temperature),
            (temperature > 1000.0 ? 5.980528 + 2.92664 * log(temperature) + 0.0014879768 * temperature + (3.3656793333333334e-11) * pow(temperature, 3) + (-2.84238e-07) * temperature * temperature + (-1.68833775e-15) * pow(temperature, 4) : 3.950372 + 3.298677 * log(temperature) + 0.0014082404 * temperature + (1.8805050000000002e-09) * pow(temperature, 3) + (-6.112135e-13) * pow(temperature, 4) + (-1.981611e-06) * temperature * temperature),
            };
        return s0_r;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static SpeciesT get_species_gibbs_rt(ContainerT temperature)
    {
        SpeciesT h0_rt = get_species_enthalpies_rt(temperature);
        SpeciesT s0_r = get_species_entropies_r(temperature);
        SpeciesT g0_rt = {
        h0_rt[0] - s0_r[0],
        h0_rt[1] - s0_r[1],
        };
        return g0_rt;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ReactionsT get_equilibrium_constants(ContainerT temperature)
    {
        ContainerT rt = gas_constant * temperature;
        ContainerT c0 = std::log(one_atm/rt);
        SpeciesT g0_rt = get_species_gibbs_rt(temperature);
        ReactionsT k_eq = {
        1.0 * g0_rt[1] + 1.0 * g0_rt[0] + 0 - (1.0 * g0_rt[1] + 1.0 * g0_rt[0] + 0),
        };
        return k_eq;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ContainerT get_temperature(
        ContainerT energy_or_enthalpy,
        ContainerT t_guess,
        SpeciesT const &mass_fractions,
        bool const do_energy = true)
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
        double tol = 1.0e-06;
        ContainerT iter_temp = t_guess;

        for(int iter = 0; iter < num_iter; ++iter){
            auto iter_rhs   = energy_or_enthalpy - he_fun(iter_temp, mass_fractions);
            auto iter_deriv = -pv_fun(iter_temp, mass_fractions);
            auto iter_dt    = -iter_rhs/iter_deriv;
            iter_temp += iter_dt;
            if(std::fabs(iter_dt) < tol){ break; }
        }
        return iter_temp;
    }


    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ReactionsT get_fwd_rate_coefficients(ContainerT temperature,
                                                SpeciesT const &concentrations)
    {
        ReactionsT k_fwd = {
        exp(24.01910270295074 + 0.0 * log(temperature) - 178.64293439206185 / temperature),
        };


        return k_fwd;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ReactionsT get_net_rates_of_progress(
        ContainerT temperature, SpeciesT const &concentrations)
    {
        ReactionsT k_fwd = get_fwd_rate_coefficients(temperature, concentrations);
        ReactionsT log_k_eq = get_equilibrium_constants(temperature);
        ReactionsT r_net = {
        k_fwd[0] * (concentrations[1] * concentrations[0] - exp(log_k_eq[0]) * concentrations[1] * concentrations[0]),
        };
        return r_net;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static SpeciesT get_net_production_rates(
        ContainerT rho, ContainerT temperature, SpeciesT const &mass_fractions)
    {
        SpeciesT concentrations = get_concentrations(rho, mass_fractions);
        ReactionsT r_net = get_net_rates_of_progress(temperature, concentrations);
        SpeciesT omega = {
        (1.0 * r_net[0] + 0 - (1.0 * r_net[0] + 0)) * (1.0 + 0 * r_net[0]),
        (1.0 * r_net[0] + 0 - (1.0 * r_net[0] + 0)) * (1.0 + 0 * r_net[0]),
        };
        return omega;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static SpeciesT get_species_viscosities(ContainerT temperature)
    {
        return SpeciesT{
            sqrt(temperature) * (0.0036188245122444136 * log(temperature) + -0.00618642807174783 + (5.9160127096038374e-05) * pow(log(temperature), 3) + (-1.9049771784292983e-06) * pow(log(temperature), 4) + (-0.0006861983404476142) * log(temperature) * log(temperature)) * (0.0036188245122444136 * log(temperature) + -0.00618642807174783 + (5.9160127096038374e-05) * pow(log(temperature), 3) + (-1.9049771784292983e-06) * pow(log(temperature), 4) + (-0.0006861983404476142) * log(temperature) * log(temperature)),
            sqrt(temperature) * (0.0031249811202346835 * log(temperature) + -0.0052325034588374736 + (5.1908396950248406e-05) * pow(log(temperature), 3) + (-1.6850989392683053e-06) * pow(log(temperature), 4) + (-0.0005968857242770352) * log(temperature) * log(temperature)) * (0.0031249811202346835 * log(temperature) + -0.0052325034588374736 + (5.1908396950248406e-05) * pow(log(temperature), 3) + (-1.6850989392683053e-06) * pow(log(temperature), 4) + (-0.0005968857242770352) * log(temperature) * log(temperature)),
        };
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static SpeciesT get_species_thermal_conductivities(ContainerT temperature)
    {
        return SpeciesT{
            sqrt(temperature) * (0.10689552101630528 + 0.014217756178712414 * log(temperature) * log(temperature) + (5.09174438888175e-05) * pow(log(temperature), 4) + (-0.06376711343770658) * log(temperature) + (-0.0013908411531386269) * pow(log(temperature), 3)),
            sqrt(temperature) * (0.0026129214099176513 + 0.0015932386446981589 * log(temperature) + 0.00016507154037197712 * pow(log(temperature), 3) + (-8.297315315373045e-06) * pow(log(temperature), 4) + (-0.000984277527739843) * log(temperature) * log(temperature)),
        };
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static Species2T get_species_binary_mass_diffusivities(ContainerT temperature)
    {
        return Species2T{
            SpeciesT{
                sqrt(temperature) * temperature * (0.0016332321885344269 * log(temperature) + -0.0031272899371257495 + (2.3795154191077426e-05) * pow(log(temperature), 3) + (-7.135754459127752e-07) * pow(log(temperature), 4) + (-0.00029023244731892585) * log(temperature) * log(temperature)),
                sqrt(temperature) * temperature * (0.0015963119455941247 * log(temperature) + -0.00303005347627421 + (2.3589003022566823e-05) * pow(log(temperature), 3) + (-7.128139183150021e-07) * pow(log(temperature), 4) + (-0.00028558271238060057) * log(temperature) * log(temperature)),
            },
            SpeciesT{
                sqrt(temperature) * temperature * (0.0015963119455941247 * log(temperature) + -0.00303005347627421 + (2.3589003022566823e-05) * pow(log(temperature), 3) + (-7.128139183150021e-07) * pow(log(temperature), 4) + (-0.00028558271238060057) * log(temperature) * log(temperature)),
                sqrt(temperature) * temperature * (0.001572491741915425 * log(temperature) + -0.00295673852294946 + (2.364377574630533e-05) * pow(log(temperature), 3) + (-7.213510385775298e-07) * pow(log(temperature), 4) + (-0.000283699905789111) * log(temperature) * log(temperature)),
            },
        };
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ContainerT get_mixture_viscosity_mixavg(
        ContainerT temperature, SpeciesT mass_fractions)
    {
        ContainerT mmw = get_mixture_molecular_weight(mass_fractions);
        SpeciesT mole_fractions = get_mole_fractions(mmw, mass_fractions);
        SpeciesT viscosities = get_species_viscosities(temperature);
        SpeciesT mix_rule_f = {
            0 + (mole_fractions[1] * (sqrt(viscosities[0] / viscosities[1] * sqrt(0.8754922182636414)) + 1) * (sqrt(viscosities[0] / viscosities[1] * sqrt(0.8754922182636414)) + 1)) / sqrt(17.137716855857786) + (mole_fractions[0] * (sqrt(viscosities[0] / viscosities[0] * sqrt(1.0)) + 1) * (sqrt(viscosities[0] / viscosities[0] * sqrt(1.0)) + 1)) / sqrt(16.0),
            0 + (mole_fractions[1] * (sqrt(viscosities[1] / viscosities[1] * sqrt(1.0)) + 1) * (sqrt(viscosities[1] / viscosities[1] * sqrt(1.0)) + 1)) / sqrt(16.0) + (mole_fractions[0] * (sqrt(viscosities[1] / viscosities[0] * sqrt(1.1422146069822232)) + 1) * (sqrt(viscosities[1] / viscosities[0] * sqrt(1.1422146069822232)) + 1)) / sqrt(15.00393774610913),
        };

        return (0.0 +
            + mole_fractions[0]*viscosities[0]/mix_rule_f[0]
            + mole_fractions[1]*viscosities[1]/mix_rule_f[1]
        );
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static ContainerT get_mixture_thermal_conductivity_mixavg(
        ContainerT temperature, SpeciesT mass_fractions)
    {
        ContainerT mmw = get_mixture_molecular_weight(mass_fractions);
        SpeciesT mole_fractions = get_mole_fractions(mmw, mass_fractions);
        SpeciesT conductivities = get_species_thermal_conductivities(temperature);

        ContainerT lhs = (0.0 +
            + mole_fractions[0]*conductivities[0]
            + mole_fractions[1]*conductivities[1]
        );

        ContainerT rhs = 1.0 / (0.0 +
            + mole_fractions[0]/conductivities[0]
            + mole_fractions[1]/conductivities[1]
        );

        return 0.5*(lhs + rhs);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static SpeciesT get_species_mass_diffusivities_mixavg(
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

};
}
