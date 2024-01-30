/*
 * ID       : EXP RESEARCH LAB CP 23 3
 * Author   : CHLORINE PENTOXIDE and STELLAR TERROR
 * DOBE     : 23 01 2024
 * STATUS   : OPEN
 * OBJECTIVE: SLEEP
 */

#include <cmath>
#include <cxxplot/cxxplot>
#include "../../lib/types.h"
#include "../../solvers/euler.h"

namespace plt = cxxplot;

int main(int argc, char* argv[]) {

    // STATE VECTOR
    const double Pe = 10125;
    const double P0 = 10 * Pe;
    const double x0 = 0.33;
    const double Lpc = 1;
    const double a = M_PI * pow(35.0 / 2000, 2);
    const double A = M_PI * pow(75.0 / 2000, 2);
    const double Vc = Lpc * A;
    const double rho = 1000;
    const double g = 9.8;
    const double adbgamma = 1.4;
    const double m0 = 0.5;
    const double sound_speed = 300;

    auto dPdt_adiabatic = [&](vector_type v, double t) -> double {
        double P = v[0];
        double V = v[1];
        double x = Lpc - (P0 / P) * (Lpc - x0);
        double ve = sqrt(((2 * (P - Pe) / rho) + 2 * g * x) / (1 - pow(a / A, 2)));
        double dxdt = (a / A) * ve;
        return adbgamma * (P / V) * -1 * A * dxdt;
    };

    auto dVdt_adiabatic = [&](vector_type v, double t) -> double {
        double P = v[0];
        double x = Lpc - (P0 / P) * (Lpc - x0);
        double ve = sqrt(((2 * (P - Pe) / rho) + 2 * g * x) / (1 - pow(a / A, 2)));
        double dxdt = (a / A) * ve;
        return A * dxdt;
    };

    auto dPdt_isothermal = [&](vector_type v, double t) -> double {
        double P = v[0];
        double V = v[1];
        double x = Lpc - (P0 / P) * (Lpc - x0);
        double ve = sqrt(((2 * (P - Pe) / rho) + 2 * g * x) / (1 - pow(a / A, 2)));
        double dxdt = (a / A) * ve;
        return 1 * (P / V) * -1 * A * dxdt;
    };

    auto dVdt_isothermal = [&](vector_type v, double t) -> double {
        double P = v[0];
        double x = Lpc - (P0 / P) * (Lpc - x0);
        double ve = sqrt(((2 * (P - Pe) / rho) + 2 * g * x) / (1 - pow(a / A, 2)));
        double dxdt = (a / A) * ve;
        return A * dxdt;
    };

    auto dPdt_cpst_2 = [&](vector_type v, double t) -> double {
        double P = v[0];
        double V = v[1];
        double fluidmass = rho * (Vc - V);
        double x = (Vc - V) / A;
        double ve = sqrt(
                ((2 * (P - Pe) / rho) + 2 * g * x) / (1 - pow(a / A, 2) - (2 * rho * a * x) / (m0 + fluidmass)));
        double dxdt = (a / A) * ve;
        return 1 * (P / V) * -1 * A * dxdt;
    };

    auto dVdt_cpst_2 = [&](vector_type v, double t) -> double {
        double P = v[0];
        double V = v[1];
        double fluidmass = rho * (Vc - V);
        double x = (Vc - V) / A;
        double ve = sqrt(
                ((2 * (P - Pe) / rho) + 2 * g * x) / (1 - pow(a / A, 2) - (2 * rho * a * x) / (m0 + fluidmass)));
        double dxdt = (a / A) * ve;
        return A * dxdt;
    };

    auto dvdt_cpst_2 = [&](vector_type v, double t) -> double {
        double P = v[0];
        double V = v[1];
        double fluidmass = rho * (Vc - V);
        double x = (Vc - V) / A;
        double ve = sqrt(
                ((2 * (P - Pe) / rho) + 2 * g * x) / (1 - pow(a / A, 2) - (2 * rho * a * x) / (m0 + fluidmass)));
        double Ft = rho * A * ve * ve;
        return (Ft / fluidmass) - g;
    };

    auto cpst_3_converge_H = [&] (double Vw) -> double {
        double Vref = 0.0;
        double dz = 1e-4;
        for(double z = 0.0; z <= Lpc; z += dz) {
            Vref += A * dz;
            if(Vw - Vref < 1e-3) return z;
        }
        return 0;
    };

    auto dPdt_cpst_3 = [&](vector_type v, double t) -> double {
        double Pressure = v[0];
        double Volume_Air = v[1];
        double Exit_Velocity = v[3];
        return adbgamma * (Pressure / Volume_Air) * a * Exit_Velocity;
    };

    auto dVdt_cpst_3 = [&](vector_type v, double t) -> double {
        double Exit_Velocity = v[3];
        return -1 * a * Exit_Velocity;
    };

    auto dvdt_cpst_3 = [&](vector_type v, double t) -> double {
        double Pressure = v[0];
        double Volume_Air = v[1];
        double Exit_Velocity = v[3];
        double Force_Thrust = rho * a * Exit_Velocity * Exit_Velocity;
        double fluidmass = (Vc - Volume_Air) * rho;
        return Force_Thrust / (fluidmass + m0) - g;
    };

    auto dvedt_cpst_3 = [&](vector_type v, double t) -> double {
        double Pressure = v[0];
        double Volume_Air = v[1];
        double Exit_Velocity = v[3];
        double H = cpst_3_converge_H(Vc - Volume_Air);
        double fluidmass = (Vc - Volume_Air) * rho;

        double Alpha = 0.5 * (pow(a /  A, 2) - 1) + ((0 * rho * a * pow(Exit_Velocity, 2))/ fluidmass + g) * H;
        double Beta = (Pressure - Pe) / rho;
        double Omega = (a / A) * H;

        return -1 * (Alpha * pow(Exit_Velocity, 2) + Beta) / Omega;
    };

    return cxxplot::exec(argc, argv, [&]() {
        std::vector<plt::point2d> data;
        std::vector<plt::point2d> data2;
        std::vector<plt::point2d> data3;
        std::vector<plt::point2d> data4;
        using namespace plt::named_parameters;
        auto w = plt::plot(
                data,
                window_title_
                        = "Orbiton Labs Experimental Model : CHLORINE PENTOXIDE - STELLAR TERROR MODELS",
                window_size_ = {800, 800},
                auto_redraw_ = true,
                line_color_ = plt::color::green,
                show_legend_ = true,
                legend_alignment_ = plt::VerticalAlignment::Top | plt::HorizontalAlignment::Right);

        auto &f = w.figure(0);
        auto &g0 = f.graph(0);
        g0.name = "Stock Model";

        auto &g1 = w.add_graph(data2,
                               line_color_ = plt::color::rgb(240, 50, 50));
        g1.name = "CPST Mk. 1";

        auto &g2 = w.add_graph(data3,
                               line_color_ = plt::color::rgb(140, 25, 75));
        g2.name = "CPST Mk. 2";

        auto &g3 = w.add_graph(data4,
                               line_color_ = plt::color::rgb(45, 175, 10));
        g3.name = "CPST Mk. 3";

        int VIEW_CHANNEL = 2; // 0 - PRESSURE, 1 - VOLUME (AIR)

        functional_vector isothermal_fns{dPdt_isothermal, dVdt_isothermal};
        vector_type isoth_state = {P0, (Lpc - x0) * A};

        auto isoth_relay = [&](vector_type v, double t) -> int {
            double P = v[0];
            double x = Lpc - (P0 / P) * (Lpc - x0);
            g0.append_data(t, v[VIEW_CHANNEL]);
            if (x < 1e-5) return -1;
            return 0;
        };
        ode_fe(isothermal_fns, isoth_state, 0, 100, 1e-6, isoth_relay, 1);

        functional_vector adb_fns{dPdt_adiabatic, dVdt_adiabatic};
        vector_type adiabatic_state = {P0, (Lpc - x0) * A};

        auto adiabatic_relay = [&](vector_type v, double t) -> int {
            double P = v[0];
            double x = Lpc - (P0 / P) * (Lpc - x0);
            g1.append_data(t, v[VIEW_CHANNEL]);
            if (x < 1e-5) return -1;
            return 0;
        };
        ode_fe(adb_fns, adiabatic_state, 0, 100, 1e-6, adiabatic_relay, 1);

        functional_vector cpst_2_fns{dPdt_cpst_2, dVdt_cpst_2, dvdt_cpst_2};
        vector_type cpst_2_state = {P0, (Lpc - x0) * A, 0};

        auto cpst_2_relay = [&](vector_type v, double t) -> int {
            double P = v[0];
            double V = v[1];
            double fluidmass = rho * (Vc - V);
            double x = (Vc - V) / A;
            g2.append_data(t, v[VIEW_CHANNEL]);
            if (x < 1e-5) return -1;
            return 0;
        };
        ode_fe(cpst_2_fns, cpst_2_state, 0, 100, 1e-6, cpst_2_relay, 1);

        functional_vector cpst_3_fns{dPdt_cpst_3, dVdt_cpst_3, dvdt_cpst_3, dvedt_cpst_3};
        vector_type cpst_3_state = {P0, (Lpc - x0) * A, 0, 0};

        auto cpst_3_relay = [&](vector_type v, double t) -> int {
            double V = v[1];
            double x = cpst_3_converge_H(Vc - V);
            g3.append_data(t, v[VIEW_CHANNEL]);
            //std::cout << x << "\n";
            if (x < 1e-5) return -1;
            return 0;
        };
        ode_fe(cpst_3_fns, cpst_3_state, 0, 100, 1e-7, cpst_3_relay, 1);

        return 0;
    });
}