//
// Created by chlorinepentoxide on 30/1/24.
//

#include <functional>
#include <cmath>

#ifndef FASTSIM_LYNX_TYPES_H
#define FASTSIM_LYNX_TYPES_H

struct RocketParameters {
    double mass_static;
    double pressure_initial;
    double radius_chamber;
    double radius_nozzle;
    double height_chamber;
    double height_propellant;
    double density_propellant;
    double coefficient_drag;
};

struct EnvironmentParameters {
    const double gravity = 9.8;
    std::function<double (double, double)> density_air;
    const double pressure_atm = 101325;
    const double gas_gamma = 1.4;
    const double temperature = 298.15;
    const double gas_constant = 287;
};

struct CalculatedParameters {
    double mass_static;
    double pressure_initial;
    double radius_chamber;
    double radius_nozzle;
    double height_chamber;
    double height_propellant;
    double density_propellant;
    double coefficient_drag;

    double area_chamber;
    double area_nozzle;
    double volume_chamber;
    double volume_fluid;
    double transient_mass;
};

struct StateVector {
    double pressure{};
    double volume_fluid{};
    double temperature{};
    double exit_velocity{};
    double transient_mass{};
    double rocket_velocity{};
    double rocket_angle{};
    double rocket_height{};
    CalculatedParameters rocket{};
    EnvironmentParameters env;
};

CalculatedParameters calculate_rocket_parameters(RocketParameters& parameters, EnvironmentParameters& eniv) {
    return {
        parameters.mass_static,
        parameters.pressure_initial,
        parameters.radius_chamber,
        parameters.radius_nozzle,
        parameters.height_chamber,
        parameters.height_propellant,
        parameters.density_propellant,
        parameters.coefficient_drag,
        M_PI * pow(parameters.radius_chamber, 2),
        M_PI * pow(parameters.radius_nozzle, 2),
        M_PI * pow(parameters. radius_chamber, 2) * parameters.height_chamber,
        M_PI * pow(parameters. radius_chamber, 2) * parameters.height_propellant,
        M_PI * pow(parameters. radius_chamber, 2) * parameters.height_propellant * parameters.density_propellant
        + M_PI * pow(parameters. radius_chamber, 2) * (parameters.height_chamber - parameters.height_propellant) * eniv.density_air(parameters.pressure_initial, eniv.temperature),
    };
}

#endif //FASTSIM_LYNX_TYPES_H
