#include <iostream>
#include <valarray>
#include <fstream>
#include <chrono>

#include "lib/types.h"
#include "solvers/euler.h"
#include "solvers/runge_kutta.h"

double f1(vector_type state, double m) {
    return state[1];
}

double f2(vector_type state, double m) {
    return 1.766 * (1 - state[0] * state[0]) * state[1] - state[0];
}

int main_off() {
    std::cout << "OrbitonLabs FastSim 2024" << std::endl;
    std::cout << "Copyright(C) 2024, Simulation Department, Orbiton Labs" << std::endl;
    std::cout << "Licensed under MIT License" << std::endl;
    std::cout << std::endl;

    functional_vector functions{ f1, f2 };
    vector_type init{ 1, 1 };

    std::ofstream out("coupled_euler_solve.csv");
    out << "x,y,t" << "\n";
    auto t1 = std::chrono::high_resolution_clock::now();
    ode_fe(functions, init, 0, 7.5, 1e-4, out);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "FE: " << dt.count() << "ms" << std::endl;
    out.close();

    std::ofstream out2("coupled_rk4_solve.csv");
    out2 << "x,y,t" << "\n";
    auto t12 = std::chrono::high_resolution_clock::now();
    ode_rk4(functions, init, 0, 7.5, 1e-4, out2);
    auto t22 = std::chrono::high_resolution_clock::now();
    auto dt2 = std::chrono::duration_cast<std::chrono::milliseconds>(t22 - t12);
    std::cout << "RK4: " << dt2.count() << "ms" << std::endl;
    out2.close();

    _rk_odesolver_tolerance = 1e-6;
    std::ofstream out3("coupled_rkf45_solve.csv");
    out3 << "x,y,t" << "\n";
    auto t13 = std::chrono::high_resolution_clock::now();
    ode_rk45orig(functions, init, 0, 7.5, out3);
    auto t23 = std::chrono::high_resolution_clock::now();
    auto dt3 = std::chrono::duration_cast<std::chrono::milliseconds>(t23 - t13);
    std::cout << "RKF45: " << dt3.count() << "ms" << std::endl;
    out3.close();


    _rk_odesolver_tolerance = 1e-2;
    std::ofstream out5("coupled_rk45_solve.csv");
    out5 << "x,y,t" << "\n";
    auto t15 = std::chrono::high_resolution_clock::now();
    ode_rk45(functions, init, 0, 7.5, out5);
    auto t25 = std::chrono::high_resolution_clock::now();
    auto dt5 = std::chrono::duration_cast<std::chrono::milliseconds>(t25 - t15);
    std::cout << "RK45: " << dt5.count() << "ms" << std::endl;
    out5.close();

    _rk_odesolver_tolerance = 1e-6;
    std::ofstream out6("coupled_rkcp45_solve.csv");
    out6 << "x,y,t" << "\n";
    auto t16 = std::chrono::high_resolution_clock::now();
    ode_rk45cp(functions, init, 0, 7.5, out6);
    auto t26 = std::chrono::high_resolution_clock::now();
    auto dt6 = std::chrono::duration_cast<std::chrono::milliseconds>(t26 - t16);
    std::cout << "RKCP45: " << dt6.count() << "ms" << std::endl;
    out6.close();

    _rk_odesolver_tolerance = 1e-6;
    std::ofstream out7("coupled_rk12eh_solve.csv");
    out7 << "x,y,t" << "\n";
    auto t17 = std::chrono::high_resolution_clock::now();
    ode_rk12eh(functions, init, 0, 7.5, out7);
    auto t27 = std::chrono::high_resolution_clock::now();
    auto dt7 = std::chrono::duration_cast<std::chrono::milliseconds>(t27 - t17);
    std::cout << "RKHE12: " << dt7.count() << "ms" << std::endl;
    out6.close();

    _rk_odesolver_tolerance = 1e-4;
    std::ofstream out8("coupled_rk23bs_solve.csv");
    out8 << "x,y,t" << "\n";
    auto t18 = std::chrono::high_resolution_clock::now();
    ode_rk23bs(functions, init, 0, 7.5, out8);
    auto t28 = std::chrono::high_resolution_clock::now();
    auto dt8 = std::chrono::duration_cast<std::chrono::milliseconds>(t28 - t18);
    std::cout << "RKBS23: " << dt8.count() << "ms" << std::endl;
    out6.close();

    std::ofstream outA("coupled_backward_euler_solve.csv");
    outA << "x,y,t" << "\n";
    auto t1A = std::chrono::high_resolution_clock::now();
    ode_be(functions, init, 0, 7.5, 1e-2, outA);
    auto t2A = std::chrono::high_resolution_clock::now();
    auto dtA = std::chrono::duration_cast<std::chrono::milliseconds>(t2A - t1A);
    std::cout << "BE(+NR): " << dtA.count() << "ms" << std::endl;
    outA.close();

    std::ofstream outA2("coupled_cnrk2_solve.csv");
    outA2 << "x,y,t" << "\n";
    auto t1A2 = std::chrono::high_resolution_clock::now();
    ode_cnrk2(functions, init, 1e-2, 0, 7.5, outA2);
    auto t2A2 = std::chrono::high_resolution_clock::now();
    auto dtA2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2A2 - t1A2);
    std::cout << "CNRK2(+NR): " << dtA2.count() << "ms" << std::endl;
    outA2.close();

    std::ofstream outA3("coupled_qzdirk2_solve.csv");
    outA3 << "x,y,t" << "\n";
    auto t1A3 = std::chrono::high_resolution_clock::now();
    ode_qzdirk2(functions, init, 1e-2, 0, 7.5, outA3);
    auto t2A3 = std::chrono::high_resolution_clock::now();
    auto dtA3 = std::chrono::duration_cast<std::chrono::milliseconds>(t2A3 - t1A3);
    std::cout << "QZDIRK2(+NR): " << dtA3.count() << "ms" << std::endl;
    outA3.close();

    return 0;
}