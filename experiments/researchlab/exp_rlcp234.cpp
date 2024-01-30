/*
 * ID       : EXP RESEARCH LAB CP 23 4
 * Author   : CHLORINE PENTOXIDE and STELLAR TERROR
 * DOBE     : 26 01 2024
 * STATUS   : OPEN
 * OBJECTIVE: TRANSIENT FLOW BERNOULLI WITH ADIABATIC AND PSEUDO-FORCE CORRECTION
 */

#include <cmath>
#include <cxxplot/cxxplot>
#include "../../lib/types.h"
#include "../../solvers/euler.h"

namespace plt = cxxplot;

const double R = 287.052874;

double rho_a(double P, double T) {
    return P / (R * T);
}

// STATE VECTOR
const double Pe = 10125;
const double P0 = 10 * Pe;
const double Temp = 273.15 + 28;
const double x0 = 0.55;
const double Lpc = 1;
const double area_exit = M_PI * pow(35.0 / 2000, 2);
const double area_chamber = M_PI * pow(75.0 / 2000, 2);
const double Vc = Lpc * area_chamber;
const double Vf = x0 * area_chamber;
const double rho_w = 1000;
const double g = 9.8;
const double adbgamma = 1.4;
const double m0 = 0.5;
const double ma = rho_a(P0, Temp) * (Vc - Vf);
const double mt = rho_w * Vf + ma;
const double CD = 0.5;
const double eta = 0.5 * CD * rho_a(Pe, Temp);

/*
 * State Vector Defintions
 */
const unsigned int PRESSURE_CHANNEL = 0;
const unsigned int FLUID_VOLUME_CHANNEL = 1;
const unsigned int FLUID_HEIGHT_CHANNEL = 2;
const unsigned int TRANSIENT_MASS_CHANNEL = 3;
const unsigned int EXIT_VELOCITY_CHANNEL = 4;
const unsigned int ROCKET_VELOCITY_CHANNEL = 5;
const unsigned int ROCKET_THRUST = 6;
int WATCH_CHANNEL = ROCKET_VELOCITY_CHANNEL;              // Change for plot view

/* Rocket Model Equations */

double dmtdt(vector_type& state) {
    return rho_w * area_exit * state[EXIT_VELOCITY_CHANNEL];
}

double thrust(vector_type& state) {
    return rho_w * area_exit * pow(state[EXIT_VELOCITY_CHANNEL], 2);
}

double drag(vector_type& state) {
    return -1 * eta * state[ROCKET_VELOCITY_CHANNEL] * fabs(state[ROCKET_VELOCITY_CHANNEL]);
}

double dvdt(vector_type& state) {
    return (thrust(state) + drag(state)) / (m0 + state[TRANSIENT_MASS_CHANNEL]) - g;
}

double fluid_pseudoacceleration(vector_type& state) {
    return (thrust(state) + drag(state)) / (state[TRANSIENT_MASS_CHANNEL] - ma) - g;
}

/* GAS EXPANSION EQUATIONS */

double dVwdt(vector_type& state) {
    return area_exit * state[EXIT_VELOCITY_CHANNEL];
}

double dV1dt(vector_type& state) {
    return -1 * dVwdt(state);
}

double dHdt(vector_type& state) {
    return (1 / area_chamber) * dVwdt(state);
}

double d2Hdt2_mu_cofactor(vector_type& state) {
    return (area_exit / area_chamber);
}

double dPdt(vector_type& state) {
    double P = state[PRESSURE_CHANNEL];
    double Va = Vc - state[FLUID_VOLUME_CHANNEL];
    return -1 * adbgamma * (P / Va) * dV1dt(state);
}

double d2Pdt2_lambda(vector_type& state) {
    return -1 * (adbgamma / (Vc - state[FLUID_VOLUME_CHANNEL])) * dV1dt(state) * dPdt(state)
           + (adbgamma / pow((Vc - state[FLUID_VOLUME_CHANNEL]), 2)) * state[PRESSURE_CHANNEL] * pow(dV1dt(state), 2);
}

double d2Pdt2_mu_cofactor(vector_type& state) {
    return area_exit * (adbgamma / (Vc - state[FLUID_VOLUME_CHANNEL])) * state[PRESSURE_CHANNEL];
}

/* STEADY FLOW EQUATIONS - GUESS MODEL */

double steadyflow_exitvelocity(vector_type& state) {
    double P = state[PRESSURE_CHANNEL];
    double H = state[FLUID_HEIGHT_CHANNEL];
    double pressure_head = (2 / rho_w) * (P - Pe);
    double potential_factor = 2 * (g + fluid_pseudoacceleration(state)) * H;
    double curvature_factor = (1 - pow(area_exit / area_chamber, 2));
    return -1 * sqrt((pressure_head + potential_factor) / (curvature_factor));
}

// del/delP (ve) - Calculated for Steady Flow
double steadyflow_alpha(vector_type& state) {
    double curvature_factor = (1 - pow(area_exit / area_chamber, 2));
    return 1 / (rho_w * curvature_factor * steadyflow_exitvelocity(state));
}

// del/delH (ve) - Calculated for Steady Flow
double steadyflow_beta(vector_type& state) {
    double H = state[FLUID_HEIGHT_CHANNEL];
    double curvature_factor = (1 - pow(area_exit / area_chamber, 2));
    return g / (rho_w * curvature_factor * steadyflow_exitvelocity(state));
}

// del/delap (ve) - Calculated for Steady Flow
double steadyflow_delta(vector_type& state) {

}

// 1st order term - Calculated at Steady Flow
double steadyflow_dvedt(vector_type& state) {
    return steadyflow_alpha(state) * dPdt(state) + steadyflow_beta(state) * dHdt(state); // + Sigma (Pseudo-force correction)
}

// 1st order guess for unsteady term - Calculated at Steady Flow; applied omega-approximation
double steadyflow_omega(vector_type& state) {
    return state[FLUID_HEIGHT_CHANNEL] * steadyflow_dvedt(state);
}

/* TRANSIENT FLOW EQUATIONS */

// Transient Flow Exit Velocity with Omega
double transientflow_exitvelocity(vector_type& state, double omega) {
    double P = state[PRESSURE_CHANNEL];
    double H = state[FLUID_HEIGHT_CHANNEL];
    double pressure_head = (2 / rho_w) * (P - Pe);
    double potential_factor = 2 * (g + fluid_pseudoacceleration(state)) * H;
    double curvature_factor = (1 - pow(area_exit / area_chamber, 2));
    double unsteady_term = 2 * omega;
    return -1 * sqrt((pressure_head + potential_factor + unsteady_term) / (curvature_factor));
}

// del/delP (ve) - Calculated for Transient Flow
double transientflow_alpha(vector_type& state, double omega) {
    double curvature_factor = (1 - pow(area_exit / area_chamber, 2));
    return 1 / (rho_w * curvature_factor * transientflow_exitvelocity(state, omega));
}

double transientflow_dalphadt_mu_cofactor(vector_type& state, double omega) {
    double curvature_factor = (1 - pow(area_exit / area_chamber, 2));
    double ve = transientflow_exitvelocity(state, omega);
    return -1 / (rho_w * curvature_factor * ve * ve);
}

double transientflow_dbetadt_mu_cofactor(vector_type& state, double omega) {
    double curvature_factor = (1 - pow(area_exit / area_chamber, 2));
    double ve = transientflow_exitvelocity(state, omega);
    return -(g) / (rho_w * curvature_factor * ve * ve);
}

double transientflow_dbetadt_lambda(vector_type& state, double omega) {
    return 0;
}

// del/delH (ve) - Calculated for Transient Flow
double transientflow_beta(vector_type& state, double omega) {
    double H = state[FLUID_HEIGHT_CHANNEL];
    double curvature_factor = (1 - pow(area_exit / area_chamber, 2));
    return (g) / (rho_w * curvature_factor * transientflow_exitvelocity(state, omega));
}

// 2nd order term - Calculated at Transient Flow
double transientflow_dvedt(vector_type& state, double omega) {
    double lambda1 = transientflow_alpha(state, omega) * dPdt(state) + transientflow_beta(state, omega) * dHdt(state);
    double curvature_factor = (1 - pow(area_exit / area_chamber, 2));
    double k = (-1 / (curvature_factor * transientflow_exitvelocity(state, omega)));
    double mu1 = k * (area_exit / area_chamber) * dHdt(state);
    double epsilon = state[FLUID_HEIGHT_CHANNEL] * (area_exit / area_chamber);
    double mu2 = k * epsilon * transientflow_dalphadt_mu_cofactor(state, omega) * dPdt(state);
    double lambda2 = k * epsilon * transientflow_alpha(state, omega) * d2Pdt2_lambda(state);
    double mu3 = k * epsilon * transientflow_alpha(state, omega) * d2Pdt2_mu_cofactor(state);
    double mu4 = k * epsilon * transientflow_dbetadt_mu_cofactor(state, omega) * dHdt(state);
    double lambda3 = k * epsilon * transientflow_dbetadt_lambda(state, omega) * dHdt(state);
    double mu5 = k * epsilon * transientflow_beta(state, omega) * d2Hdt2_mu_cofactor(state);
    double lambda = lambda1 + lambda2 + lambda3;
    double mu = mu1 + mu2 + mu3 + mu4 + mu5;
    return (lambda) / (1 - mu);
}

double transientflow_omega(vector_type state, double omega) {
    return state[FLUID_HEIGHT_CHANNEL] * transientflow_dvedt(state, omega);
}

// Convergence Solver for Exit Velocity
double convegence_tol = 1e-6;
int convergence_max = 25;
int converge_transient_exitvelocity(vector_type& state, double& converged_exitvelocity, double& converged_omega) {
    double ve = steadyflow_exitvelocity(state);
    state[EXIT_VELOCITY_CHANNEL] = ve;
    double omega = steadyflow_omega(state);
    double deviation = fabs(ve - transientflow_exitvelocity(state, omega));
    double prev_dev = deviation;
    for(int i = 1; i <= convergence_max; i++) {
        ve = transientflow_exitvelocity(state, omega);
        state[EXIT_VELOCITY_CHANNEL] = ve;
        omega = transientflow_omega(state, omega);
        deviation = fabs(ve - transientflow_exitvelocity(state, omega));

        converged_exitvelocity = ve;
        converged_omega = omega;
        if(deviation <= convegence_tol) return i;
        if(deviation - prev_dev > 0) return -2;
        prev_dev = deviation;
    }
    return -1;
}

int main(int argc, char* argv[]) {

    return cxxplot::exec(argc, argv, [&]() {
        std::vector<plt::point2d> data;
        std::vector<plt::point2d> data2;
        using namespace plt::named_parameters;
        auto w = plt::plot(
                data,
                window_title_
                        = "Orbiton Labs Experimental Model : CHLORINE PENTOXIDE - STELLAR TERROR MODELS",
                window_size_ = {800, 800},
                auto_redraw_ = true,
                line_color_ = plt::color::blue,
                show_legend_ = true,
                legend_alignment_ = plt::VerticalAlignment::Top | plt::HorizontalAlignment::Right);

        auto &f = w.figure(0);
        auto &g0 = f.graph(0);
        g0.name = "Transient Flow Model";

        auto &g1 = w.add_graph(data2,
                               line_color_ = plt::color::rgb(240, 50, 50));
        g1.name = "Steady Flow Model";

        double exitvel = 0;
        double omega = 0;

        vector_type state = {P0, Vf, x0, mt, 0, 0, 0 };
        int convergence = 0;
        double dt = 1e-6;
        for(double t = 0.0; t <= 10; t += dt) {
            if(state[FLUID_HEIGHT_CHANNEL] < 1e-3) break;
            convergence = converge_transient_exitvelocity(state, exitvel, omega);
            if(convergence < 0) break;
            state[ROCKET_THRUST] = thrust(state);
            g0.append_data(t, state[WATCH_CHANNEL]);
            vector_type diff_state = { 0, 0, 0, 0, 0, 0 };
            diff_state[PRESSURE_CHANNEL] = dPdt(state) * dt;
            diff_state[FLUID_VOLUME_CHANNEL] = dVwdt(state) * dt;
            diff_state[FLUID_HEIGHT_CHANNEL] = dHdt(state) * dt;
            diff_state[TRANSIENT_MASS_CHANNEL] = dmtdt(state) * dt;
            diff_state[EXIT_VELOCITY_CHANNEL] = 0;
            diff_state[ROCKET_VELOCITY_CHANNEL] = dvdt(state) * dt;
            state += diff_state;
        }
        if(convergence == -1) std::cout << "Convergence not attained within solver limits. Solver exited." << std::endl;
        else if(convergence == -2) std::cout << "Parameters are divergent. Solver exited." << std::endl;
        else std::cout << "Convergence was found. Solution obtained." << std::endl;

        state = {P0, Vf, x0, mt, 0, 0};
        for(double t = 0.0; t <= 10; t += dt) {
            if(state[FLUID_HEIGHT_CHANNEL] < 1e-3) break;
            state[EXIT_VELOCITY_CHANNEL] = steadyflow_exitvelocity(state);
            state[ROCKET_THRUST] = thrust(state);
            g1.append_data(t, state[WATCH_CHANNEL]);
            vector_type diff_state = { 0, 0, 0, 0, 0, 0 };
            diff_state[PRESSURE_CHANNEL] = dPdt(state) * dt;
            diff_state[FLUID_VOLUME_CHANNEL] = dVwdt(state) * dt;
            diff_state[FLUID_HEIGHT_CHANNEL] = dHdt(state) * dt;
            diff_state[TRANSIENT_MASS_CHANNEL] = dmtdt(state) * dt;
            diff_state[EXIT_VELOCITY_CHANNEL] = 0;
            diff_state[ROCKET_VELOCITY_CHANNEL] = dvdt(state) * dt;
            state += diff_state;
        }

        return 0;
    });
}