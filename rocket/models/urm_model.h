#ifndef FASTSIM_URM_MODEL_H
#define FASTSIM_URM_MODEL_H

#include <cmath>
#include <cxxplot/cxxplot>
#include "../../lib/types.h"
#include "../../solvers/euler.h"
#include "../../solvers/runge_kutta.h"

namespace plt = cxxplot;

namespace fastsim {

    struct RocketParameters {
        double dry_mass;
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
        const double air_density = 1.293;
        const double pressure_atm = 101325;
        const double gamma = 1.4;
        const double T = 298.15;
        const double R = 287;
    };

    namespace urm {

        // COMMUNICATION CHANNELS

        const unsigned int CH_MASS = 0; // Total Mass Channel
        const unsigned int CH_PRES = 1; // Pressure Channel
        const unsigned int CH_VLCT = 2; // Velocity Channel
        const unsigned int CH_HGHT = 3; // Height Channel

        // OPTIMIZED PARAMETERS
        struct OptimizedParameters {
            double Ae;
            double Ap;
            double RA;
            double ZP1;
            double ZP2;
            double CCK1;
            double ALPHA1;
            double ALPHA2;
            double ALPHA3;
            double ALPHA4;
            double P1F;
            double M1F;
            double CCK2;
            double M2F;
            double T2F;
            double P2F;
            double T3F;
            double GMF;
            double P3F;
            double CCK3;
            double BYO;
            double DRAG;
            double GR;
            double PATM;
            double INITM;
        };

        static double area(double r) {
            return M_PI * r * r;
        }

        static OptimizedParameters calculate_optimized_parameters(const RocketParameters& rocket, const EnvironmentParameters& environment) {
            double Ae = area(rocket.radius_nozzle);
            double Ap = area(rocket.radius_chamber);
            return {
                Ae,
                Ap,
                Ae / Ap,
                rocket.height_chamber,
                -1 * rocket.pressure_initial * (rocket.height_chamber - rocket.density_propellant),
                rocket.pressure_initial * (1 - rocket.height_propellant / rocket.height_chamber),
                2 / rocket.density_propellant,
                -2 * environment.pressure_atm / rocket.density_propellant + 2 * environment.gravity * rocket.height_chamber,
                -2 * environment.gravity * rocket.pressure_initial * (rocket.height_chamber - rocket.height_propellant),
                1 - pow(Ae / Ap, 2),
                rocket.pressure_initial * (rocket.height_chamber - rocket.height_propellant),
                rocket.density_propellant * Ae,
                environment.pressure_atm * pow(2.0 / (environment.gamma + 1), -1 * environment.gamma / (environment.gamma - 1)),
                Ae * sqrt((environment.gamma / (environment.R * environment.T)) * pow((2.0 / (environment.gamma + 1)), (environment.gamma + 1) / (environment.gamma - 1))),
                Ae * environment.gamma * pow((2.0 / (environment.gamma + 1)), (environment.gamma) / (environment.gamma + 1)),
                -1 * (Ae / (rocket.height_chamber * Ap)) * sqrt(environment.gamma * environment.R * environment.T * pow((2.0 / (environment.gamma + 1)), (environment.gamma + 1) / (environment.gamma - 1))),
                environment.pressure_atm * Ae * ((2 * environment.gamma) / (environment.gamma - 1)),
                (environment.gamma - 1) / environment.gamma,
                pow((environment.pressure_atm * Ae) / (rocket.height_chamber * Ap), 2) * ((environment.gamma * 2 * environment.R * environment.T)/(environment.gamma - 1)),
                environment.pressure_atm,
                environment.air_density * environment.gravity * Ap * rocket.height_chamber,
                0.5 * environment.air_density * rocket.coefficient_drag * Ap,
                environment.gravity,
                environment.pressure_atm,
                rocket.dry_mass + rocket.density_propellant * Ap * rocket.height_propellant + environment.air_density * (rocket.pressure_initial / environment.pressure_atm) * Ap * (rocket.height_chamber - rocket.height_propellant)
            };
        }

        static functional_vector generate_simfunc1(const OptimizedParameters& params) {
            auto dmdt = [&params] (vector_type state, double t) -> double {
                double P = state[CH_PRES];
                return -1 * params.M1F * sqrt((params.ALPHA1 * P + params.ALPHA2 + params.ALPHA3 / P) / (params.ALPHA4));
            };
            auto dPdt = [&params] (vector_type state, double t) -> double {
                double P = state[CH_PRES];
                return -1 * ((P * P) / params.P1F) * params.RA * sqrt((params.ALPHA1 * P + params.ALPHA2 + params.ALPHA3 / P) / (params.ALPHA4));
            };
            auto dvdt = [&params] (vector_type state, double t) -> double {
                double P = state[CH_PRES];
                double M = state[CH_MASS];
                double V = state[CH_VLCT];
                double thr = params.M1F * ((params.ALPHA1 * P + params.ALPHA2 + params.ALPHA3 / P) / (params.ALPHA4));
                return (thr + params.BYO - params.DRAG * fabs(V) * V) / M - params.GR;
            };
            auto dhdt = [&params] (vector_type state, double t) -> double {
                return state[CH_VLCT];
            };
            return {dmdt, dPdt, dvdt, dhdt};
        }

        static functional_vector generate_simfunc2(const OptimizedParameters& params) {
            auto dmdt = [&params] (vector_type state, double t) -> double {
                double P = state[CH_PRES];
                return -1 * params.M2F * P;
            };
            auto dPdt = [&params] (vector_type state, double t) -> double {
                double P = state[CH_PRES];
                return P * params.P2F;
            };
            auto dvdt = [&params] (vector_type state, double t) -> double {
                double P = state[CH_PRES];
                double M = state[CH_MASS];
                double V = state[CH_VLCT];
                double thr = params.T2F * P;
                return (thr + params.BYO - params.DRAG * fabs(V) * V) / M - params.GR;
            };
            auto dhdt = [&params] (vector_type state, double t) -> double {
                return state[CH_VLCT];
            };
            return {dmdt, dPdt, dvdt, dhdt};
        }

        static functional_vector generate_simfunc3(const OptimizedParameters& params) {
            auto dmdt = [&params] (vector_type state, double t) -> double {
                return 0;
            };
            auto dPdt = [&params] (vector_type state, double t) -> double {
                double P = state[CH_PRES];
                return -1 * sqrt(params.P3F * (1 - pow(params.PATM / P, params.GMF)));
            };
            auto dvdt = [&params] (vector_type state, double t) -> double {
                double P = state[CH_PRES];
                double M = state[CH_MASS];
                double V = state[CH_VLCT];
                double thr = params.T3F * (1 - pow(params.PATM / P, params.GMF));
                return (thr + params.BYO - params.DRAG * fabs(V) * V) / M - params.GR;
            };
            auto dhdt = [&params] (vector_type state, double t) -> double {
                return state[CH_VLCT];
            };
            return {dmdt, dPdt, dvdt, dhdt};
        }

        static functional_vector generate_simfunc4(const OptimizedParameters& params) {
            auto dmdt = [&params](vector_type state, double t) -> double {
                return 0;
            };
            auto dPdt = [&params](vector_type state, double t) -> double {
                return 0;
            };
            auto dvdt = [&params](vector_type state, double t) -> double {
                double M = state[CH_MASS];
                double V = state[CH_VLCT];
                return (params.BYO - params.DRAG * fabs(V) * V) / M - params.GR;
            };
            auto dhdt = [&params](vector_type state, double t) -> double {
                return state[CH_VLCT];
            };
            return {dmdt, dPdt, dvdt, dhdt};
        }

        static int simulate(int argc, char* argv[], RocketParameters& rocket, EnvironmentParameters& env, int watch_channel) {
            return cxxplot::exec(argc, argv, [ & ]( ) {
                std::vector<plt::point2d> data;
                using namespace plt::named_parameters;
                auto w = plt::plot(
                        data,
                        window_title_
                                = "Orbiton Labs URM Model",
                        window_size_ = {800, 800},
                        auto_redraw_ = true,
                        line_color_ = plt::color::red);

                auto &f = w.figure(0);
                auto &g0 = f.graph(0);

                OptimizedParameters opti = calculate_optimized_parameters(rocket, env);

                vector_type sync_state = { opti.INITM, rocket.pressure_initial, 0, 0 };
                double sync_time = 0;

                auto controller1 = [&] (vector_type& s, double t) -> int {
                    g0.append_data(t, s[watch_channel]);
                    if(s[CH_PRES] - opti.CCK1 < 1e-3) {
                        sync_state = s;
                        sync_time = t;
                        return -1;
                    }
                    return 0;
                };
                auto controller2 = [&] (vector_type& s, double t) -> int {
                    g0.append_data(t, s[watch_channel]);
                    if(s[CH_PRES] - opti.CCK2 < 1e-3) {
                        sync_state = s;
                        sync_time = t;
                        return -1;
                    }
                    return 0;
                };
                auto controller3 = [&] (vector_type& s, double t) -> int {
                    g0.append_data(t, s[watch_channel]);
                    if(s[CH_PRES] - opti.CCK3 < 1e-3) {
                        sync_state = s;
                        sync_time = t;
                        return -1;
                    }
                    return 0;
                };
                auto controller4 = [&] (vector_type& s, double t) -> int {
                    g0.append_data(t, s[watch_channel]);
                    if(s[CH_HGHT] < 1e-3) {
                        sync_state = s;
                        sync_time = t;
                        return -1;
                    }
                    return 0;
                };

                auto fv1 = generate_simfunc1(opti);
                auto fv2 = generate_simfunc2(opti);
                auto fv3 = generate_simfunc3(opti);
                auto fv4 = generate_simfunc4(opti);

                ode_qzdirk2(fv1, sync_state, sync_time, 1, 1e-6, controller1, 1);
                ode_qzdirk2(fv2, sync_state, sync_time, 1, 1e-6, controller2, 1);
                ode_be(fv3, sync_state, sync_time, 1, 1e-6, controller3, 1);
                ode_fe(fv4, sync_state, sync_time, 100, 1e-4, controller4, 5);
                return 0;
            });
        }

    }

}

#endif //FASTSIM_URM_MODEL_H
