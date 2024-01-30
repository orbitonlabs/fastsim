//
// Created by chlorinepentoxide on 30/1/24.
//

#ifndef FASTSIM_LYNX_H
#define FASTSIM_LYNX_H

#include "lynx_types.h"

/* FUNCTION DECLARATIONS */

// 1 - GAS EXPANSION MODEL EQUATIONS
double dPdt(StateVector& state);
double dV1dt(StateVector& state);
double dVwdt(StateVector& state);
double dHdt(StateVector& state);

// 2 - FLUID EJECTION MODEL EQUATIONS
double exit_velocity(StateVector& state);
double fluid_height(StateVector& state);
double dmtdt(StateVector& state);
int create_temp_store(StateVector& state);
int clear_temp_store();

// 3 - ROCKET MODEL EQUATIONS
double thrust(StateVector& state);
double drag(StateVector& state);
double dvdt(StateVector& state);
double dhdt(StateVector& state);
double fluid_pseudoacceleration(StateVector& state);

/* FUNCTION IMPLEMENTATIONS */

/* 0 - AIR DENSITY */
const double R = 287.052874;

double air_density(double P, double T) {
    return P / (R * T);
}

/* 1 - GAS EXPANSION */

// IDEAL GAS ISOTHERMAL MODEL
double idealgas_isothermal_dPdt(StateVector& state) {
    return -1 * (state.pressure / (state.rocket.volume_chamber - state.volume_fluid)) * dV1dt(state);
}

// IDEAL GAS ADIABATIC MODEL
double idealgas_adiabatic_dPdt(StateVector& state) {
    return -1 * state.env.gas_gamma * (state.pressure / (state.rocket.volume_chamber - state.volume_fluid)) * dV1dt(state);
}

double dPdt(StateVector& state) {
#ifdef LYNX_IDEAL_GAS_ISOTHERMAL
    return idealgas_isothermal_dPdt(state);
#endif
#ifdef LYNX_IDEAL_GAS_ADIABATIC
    return idealgas_adiabatic_dPdt(state);
#endif
}

// GENERIC FUNCTIONS

double dV1dt(StateVector& state) {
    return -1 * dVwdt(state);
}

double dVwdt(StateVector& state) {
    return state.rocket.area_nozzle * exit_velocity(state);
}

double dHdt(StateVector& state) {
    return (1 / state.rocket.area_chamber) * dVwdt(state);
}

/* FLUID EJECTION */
double local_store_ve;

double exit_velocity(StateVector& state) {
    return local_store_ve;
}

double fluid_height(StateVector& state) {
    return state.volume_fluid / state.rocket.area_chamber;
}

double dmtdt(StateVector& state) {
    return state.rocket.density_propellant * state.rocket.area_nozzle * exit_velocity(state);
}

int clear_temp_store() {
    local_store_ve = 0;
    return 0;
}

int create_temp_store(StateVector& state) {
    local_store_ve = -1 * sqrt( ( (2 / state.rocket.density_propellant) * (state.pressure - state.env.pressure_atm) + 2 * (state.env.gravity +
            fluid_pseudoacceleration(state)) * fluid_height(state)) / (1 - pow(state.rocket.area_nozzle / state.rocket.area_chamber, 2)) );
    return 0;
}

/* 3 - ROCKET MODEL */
double thrust(StateVector& state) {
    return state.rocket.density_propellant * state.rocket.area_nozzle * pow(exit_velocity(state), 2);
}

double drag(StateVector& state) {
    return -0.5 * state.rocket.coefficient_drag * state.env.density_air(state.env.pressure_atm, state.env.temperature)
                * state.rocket.area_chamber * state.rocket_velocity * fabs(state.rocket_velocity);
}

double dvdt(StateVector& state) {
    return (thrust(state) + drag(state)) / (state.transient_mass + state.rocket.mass_static) - state.env.gravity;
}

double fluid_pseudoacceleration(StateVector& state) {
    return (thrust(state) + drag(state)) / (state.transient_mass);
}

double dhdt(StateVector& state) {
    return state.rocket_velocity;
}

/* SOLVER SYSTEMS */
void push_next_state(StateVector& state, double step) {
    create_temp_store(state);
    state.pressure += dPdt(state) * step;
    state.volume_fluid += dVwdt(state) * step;
    state.transient_mass += dmtdt(state) * step;
    state.rocket_height += dhdt(state) * step;
    state.rocket_velocity += dvdt(state) * step;
    clear_temp_store();
}

void simple_solver(const StateVector& initial, const std::function<int (StateVector, double)>& relay) {
    double stepping = 1e-6;
    StateVector state = initial;
    for(unsigned int i = 0;;i++) {
        double t = i * stepping;
        push_next_state(state, stepping);
        int r = relay(state, t);
        if(r != 0) break;
    }
}

#endif //FASTSIM_LYNX_H
