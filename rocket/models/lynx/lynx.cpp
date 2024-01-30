//
// Created by chlorinepentoxide on 30/1/24.
//

#define LYNX_IDEAL_GAS_ADIABATIC

#include <cxxplot/cxxplot>
#include "lynx.h"

namespace plt = cxxplot;

int main(int argc, char* argv[]) {
    EnvironmentParameters env;
    env.density_air = air_density;
    RocketParameters parameters = {
        0.5,
        10 * env.pressure_atm,
        75.0 / 2000,
        35.0 / 2000,
        1,
        0.33,
        1000,
        0.5
    };

    CalculatedParameters cparams = calculate_rocket_parameters(parameters, env);
    StateVector init = {
            cparams.pressure_initial,
            cparams.volume_fluid,
            env.temperature,
            0,
            cparams.transient_mass,
            0,
            M_PI / 2,
            0,
            cparams,
            env
    };

    return cxxplot::exec(argc, argv, [&]() {
        std::vector<plt::point2d> data;
        std::vector<plt::point2d> data2;
        using namespace plt::named_parameters;
        auto w = plt::plot(
                data,
                window_title_
                        = "Orbiton Labs Lynx Project",
                window_size_ = {800, 800},
                auto_redraw_ = true,
                line_color_ = plt::color::blue,
                show_legend_ = true,
                legend_alignment_ = plt::VerticalAlignment::Top | plt::HorizontalAlignment::Right);

        auto &f = w.figure(0);
        auto &g0 = f.graph(0);
        g0.name = "Ideal Gas Isothermal Steady State";

        auto relay = [&] (StateVector state, double time) -> int {
            g0.append_data(time, state.rocket_velocity);
            if(state.volume_fluid < 1e-3) return 1;
            return 0;
        };

        simple_solver(init, relay);

        return 0;
    });
}