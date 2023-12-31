
#include "../../solvers/euler.h"

double f1h(vector_type state, double m) {
    return state[1];
}

double f2h(vector_type state, double m) {
    return 1.766 * (1 - state[0] * state[0]) * state[1] - state[0];
}

namespace plt = cxxplot;

int main1(int argc, char* argv[])
{
    return cxxplot::exec( argc, argv, [ & ]( ) {
        std::vector< plt::point2d > data;
        std::vector< plt::point2d > data2;

        using namespace plt::named_parameters;

        auto w = plt::plot(
                data,
                window_title_
                        = "Van der pol",
                window_size_      = { 650, 400 },
                auto_redraw_      = true,
                line_color_       = plt::color::rgb( 72, 171, 72 ),
                xlabel_           = "Matrix size",
                ylabel_           = "Execution time (ms)",
                show_legend_      = true,
                legend_alignment_ = plt::VerticalAlignment::Top | plt::HorizontalAlignment::Left );

        auto& f = w.figure( 0 );

        auto& g1 = w.add_graph( data2,
                                line_color_   = plt::color::rgb( 240, 50, 50 ) );

        g1.name = "Forward Euler (y)";

        auto& g0 = f.graph( 0 );
        g0.name  = "Forward Euler (x)";

        functional_vector functions{ f1h, f2h};
        vector_type init{ 1, 1 };

        auto func = [&g0, &g1] (vector_type s, double t) -> int {
            g0.append_data(t, s[0]);
            g1.append_data(t, s[1]);
            return 0;
        };

        auto t1 = std::chrono::high_resolution_clock::now();
        ode_fe(functions, init, 0, 50, 1e-6, func, 1e5);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto dt1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "w/Graphing (Freq: 5e5): " << dt1.count() << "ms" << std::endl;
        return 0;
    });
    return 0;
}
