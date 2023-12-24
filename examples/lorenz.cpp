/**
 * Lorenz Attractor
 * Written by Chlorine Pentoxide,
 * Orbiton Labs (C) 2023-24
 *
 */

#include "../lib/types.h"
#include "../solvers/euler.h"

static const double sigma = 10;
static const double rho = 28;
static const double beta = 8.0 / 3;

namespace plt = cxxplot;

int main(int argc, char* argv[]) {
    auto lorenz_x = [] (vector_type v, double t) -> double {
        return sigma * (v[1] - v[0]);
    };
    auto lorenz_y = [] (vector_type v, double t) -> double {
        return v[0] * (rho - v[2]) - v[1];
    };
    auto lorenz_z = [] (vector_type v, double t) -> double {
        return v[0] * v[1] - beta * v[2];
    };
    return cxxplot::exec( argc, argv, [ & ]( ) {
        std::vector< plt::point2d > data;
        using namespace plt::named_parameters;
        auto w = plt::plot(
                data,
                window_title_
                        = "Orbiton Labs Demo : Lorenz Map (XZ Plane)",
                window_size_      = { 800, 800 },
                auto_redraw_      = true,
                line_color_       = plt::color::green);

        auto& f = w.figure( 0 );
        auto& g0 = f.graph( 0 );

        functional_vector functions { lorenz_x, lorenz_y, lorenz_z };
        vector_type initial_values { 0.5, 1, 20 };

        auto relay = [&g0] (vector_type s, double t) -> int {
            g0.append_data(s[0], s[2]);
            return 0;
        };

        ode_be(functions, initial_values, 0, 100, 1e-5, relay, 1e3);
        return 0;
    });
}