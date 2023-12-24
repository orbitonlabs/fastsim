#pragma once

#include <functional>
#include <valarray>
#include <fstream>
#include <cxxplot/cxxplot>
#include "../lib/range.h"
#include "../lib/types.h"
#include "../solvers/newton_raphson.h"

using std::valarray;
using std::function;
using util::lang::indices;

static auto ode_fe(const valarray<vector_function_type>& differentials,
            const vector_type& initial_conditions,
            double domain_start, double domain_end, double stepping,
            const custom_parametric_function& custom_func, int freq) {
    valarray<double> state = initial_conditions;
    int count = 0;
    for (double marker = domain_start; marker <= domain_end; marker += stepping) {
        valarray<double> diff(differentials.size());
        for (auto i : indices(diff)) diff[i] = differentials[i](state, marker) * stepping;
        state += diff;
        if(count == freq) {
            auto s = custom_func(state, marker);
            if (s != 0) break;
            count = 0;
        } else {
            count++;
        }
    }
}

static auto ode_be(const functional_vector& differentials,
    const vector_type& initial_conditions,
    double domain_start, double domain_end, double stepping,
    const custom_parametric_function& custom_func, int freq) {
    vector_type state = initial_conditions;
    int count = 0;
    for (double marker = domain_start; marker <= domain_end; marker += stepping) {
        functional_vector ies(zero, differentials.size());
        for (size_t i = 0; i < differentials.size(); i++) {
            auto imp_eq = [stepping, i, differentials, state](vector_type vnext, double marker) -> double {
                return vnext[i] - stepping * differentials[i](vnext, marker + stepping) - state[i];
                };
            ies[i] = imp_eq;
        }
        state = newton_raphson(ies, state, marker, 1e-6, 5);
        if(count == freq) {
            auto s = custom_func(state, marker);
            if (s != 0) break;
            count = 0;
        } else {
            count++;
        }
    }
}