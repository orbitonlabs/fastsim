/**
 * Van der pol Oscillator Demo
 * Written by Chlorine Pentoxide,
 * Orbiton Labs (C) 2023-24
 *
 */

#include "../lib/types.h"
#include "../solvers/euler.h"

static const double mu = 1.7;
namespace plt = cxxplot;

int main(int argc, char* argv[]) {
    auto vanderpol_f1 = [] (vector_type v, double t) -> double {
        return v[1];
    };
    auto vanderpol_f2 = [] (vector_type v, double t) -> double {
        return mu * (1 - v[0] * v[0]) * v[1] - v[0];
    };
    return cxxplot::exec( argc, argv, [ & ]( ) {
        std::vector< plt::point2d > data;
        using namespace plt::named_parameters;
        auto w = plt::plot(
                data,
                window_title_
                        = "Orbiton Labs Demo : Vanderpol Poincare Map",
                window_size_      = { 800, 800 },
                auto_redraw_      = true,
                line_color_       = plt::color::green);

        auto& f = w.figure( 0 );
        auto& g0 = f.graph( 0 );
        g0.name  = "Forward Euler (x)";

        functional_vector functions { vanderpol_f1, vanderpol_f2 };
        vector_type initial_values { 1, 1 };

        auto relay = [&g0] (vector_type s, double t) -> int {
            g0.append_data(s[0], s[1]);
            return 0;
        };

        ode_fe(functions, initial_values, 0, 7.5, 1e-5, relay, 5e3);
        return 0;
    });
}