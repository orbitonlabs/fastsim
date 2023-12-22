#pragma once

#include <functional>
#include <valarray>
#include <fstream>
#include "../lib/range.h"
#include "../lib/types.h"
#include "../solvers/newton_raphson.h"

using std::valarray;
using std::function;
using util::lang::indices;

double forwardEuler(double xn, double yn, double step, function<double(double, double)> func) {
    return step * func(xn, yn);
}

double backwardEuler(double xn, double yn, double ynnext, double step, function<double(double, double)> func) {
    return step * func(xn + step, ynnext) - ynnext + yn;
}

auto ode_fe(const valarray<vector_function_type>& differentials,
    const vector_type& initial_conditions,
    double domain_start, double domain_end, double stepping) {
    valarray<double> state = initial_conditions;
    for (double marker = domain_start; marker <= domain_end; marker += stepping) {
        valarray<double> diff(differentials.size());
        for (auto i : indices(diff)) diff[i] = differentials[i](state, marker) * stepping;
        state += diff;
    }
    return state;
}

auto ode_fe(const valarray<vector_function_type>& differentials,
    const vector_type& initial_conditions,
    double domain_start, double domain_end, double stepping, std::ofstream& output) {
    valarray<double> state = initial_conditions;
    for (double marker = domain_start; marker <= domain_end; marker += stepping) {
        valarray<double> diff(differentials.size());
        for (auto i : indices(diff)) diff[i] = differentials[i](state, marker) * stepping;
        state += diff;
        for (double x : state) output << x << ",";
        output << marker << "\n";
    }
}

double zero(vector_type v, double t) {
    return 0;
}

auto ode_be(const functional_vector& differentials,
    const vector_type& initial_conditions,
    double domain_start, double domain_end, double stepping) {
    vector_type state = initial_conditions;
    for (double marker = domain_start; marker <= domain_end; marker += stepping) {
        functional_vector ies(zero, differentials.size());
        for (int i = 0; i < differentials.size(); i++) {
            auto imp_eq = [stepping, i, differentials, state](vector_type vnext, double marker) -> double {
                return vnext[i] - stepping * differentials[i](vnext, marker + stepping) - state[i];
                };
            ies[i] = imp_eq;
        }
        state = newton_raphson(ies, state, marker, 1e-6, 5);
    }
}

auto ode_be(const functional_vector& differentials,
    const vector_type& initial_conditions,
    double domain_start, double domain_end, double stepping, std::ofstream& output) {
    vector_type state = initial_conditions;
    for (double marker = domain_start; marker <= domain_end; marker += stepping) {
        functional_vector ies(zero, differentials.size());
        for (int i = 0; i < differentials.size(); i++) {
            auto imp_eq = [stepping, i, differentials, state](vector_type vnext, double marker) -> double  {
                return vnext[i] - stepping * differentials[i](vnext, marker + stepping) - state[i];
                };
            ies[i] = imp_eq;
        }
        state = newton_raphson(ies, state, marker, 1e-6, 5);
        for (double x : state) output << x << ",";
        output << marker << "\n";
    }
}