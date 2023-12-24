#pragma once


#include <functional>
#include <valarray>
#include <fstream>
#include "../lib/range.h"
#include "../lib/types.h"

using std::valarray;
using std::ostream;
using util::lang::indices;

static double _rk_odesolver_tolerance = 1e-2;
static double _rk_odesolver_default_step = 1e-6;

static auto ode_rk4(const valarray<vector_function_type>& differentials,
    const vector_type& initial_conditions,
    double domain_start, double domain_end, double stepping,
    const custom_parametric_function& pf, int freq) {
    vector_type state = initial_conditions;
    double mul = stepping / 6.0;
    int count = 0;
    auto inddiff = indices(differentials);
    for (double marker = domain_start; marker <= domain_end; marker += stepping) {
        vector_type diff1(differentials.size());
        vector_type diff2(differentials.size());
        vector_type diff3(differentials.size());
        vector_type diff4(differentials.size());
        for (auto i : inddiff) diff1[i] = differentials[i](state, marker);
        for (auto i : inddiff) diff2[i] = differentials[i](state + (stepping / 2) * diff1, marker + (stepping / 2));
        for (auto i : inddiff) diff3[i] = differentials[i](state + (stepping / 2) * diff2, marker + (stepping / 2));
        for (auto i : inddiff) diff4[i] = differentials[i](state + stepping * diff3, marker + stepping);
        state += mul * (diff1 + 2 * diff2 + 2 * diff3 + diff4);
        if(count == freq) {
            auto s = pf(state, marker);
            if (s != 0) break;
            count = 0;
        } else {
            count++;
        }
    }
    return state;
}

/* Generic RK solver with adaptive step control */

auto gensolveAB(const valarray<vector_function_type>& differentials,
    const vector_type& initial_conditions,
    double domain_start, double domain_end, double epsilon,
    const vector_type& a, const valarray<vector_type>& b, const vector_type& c, const vector_type& ch,
    int A, int B, double h_start,
    const custom_parametric_function& pf, int freq) {
    vector_type state = initial_conditions;
    auto inddiff = indices(differentials);
    auto ct = c - ch;
    double h = h_start;
    int count = 0;
    for (double marker = domain_start; marker <= domain_end;) {
        valarray<vector_type> k(B);
        for (int i = 0; i < B; i++) {
            vector_type kp = state;
            for (int j = 0; j < i; j++) kp += b[i][j] * k[j];
            vector_type kn(differentials.size());
            for (auto index : inddiff) kn[index] = h * differentials[index](kp, marker + h * a[i]);
            k[i] = kn;
        }
        vector_type new_state = state;
        vector_type errors(state.size());
        double te = 0;
        for (int i = 0; i < B; i++) {
            new_state += ch[i] * k[i];
            errors += ct[i] * k[i];
        }
        for (auto index : inddiff) te += errors[index] * errors[index];
        te = sqrt(te);
        double h_next = 0.9 * h * pow(epsilon / te, 1.0 / A);
        if (te > epsilon) {
            h = h_next;
            if(count == freq) {
                auto s = pf(state, marker);
                if (s != 0) break;
                count = 0;
            } else {
                count++;
            }
        }
        else {
            state = new_state;
            h = h_next;
            marker += h;
        }
    }
}

/* Original adaptive step control integrator based on Fehlberg */

const vector_type rk45orig_a = { 0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0 };
const valarray<vector_type> rk45orig_b = {
                                  {},
                                  {1.0 / 4.0 },
                                  {3.0 / 32, 9.0 / 32},
                                  {1932.0 / 2197, -7200.0 / 2197, 7296.0 / 2197},
                                  {439.0 / 216, -8.0, 3680.0 / 513, -845.0 / 4104},
                                  {-8.0 / 27, 2, -3544.0 / 2565, 1859.0 / 4104, -11.0 / 40}
};
const vector_type rk45orig_ch = { 16.0 / 135, 0, 6656.0 / 12825, 28561.0 / 56430, -9.0 / 50, 2.0 / 55 };
const vector_type rk45orig_c = { 25.0/216, 0, 1408.0/2565, 2197.0 / 4104, -1.0 / 5, 0 };

auto ode_rk45orig(const valarray<vector_function_type>& differentials,
    const vector_type& initial_conditions,
    double domain_start, double domain_end,
    const custom_parametric_function& pf, int freq) {
    gensolveAB(differentials, initial_conditions, domain_start, domain_end, _rk_odesolver_tolerance, rk45orig_a, rk45orig_b, rk45orig_c, rk45orig_ch, 5, 6, 1e-4, pf, freq);
}

/* Dormand-Prince RK45 adaptive step integrator */

const vector_type rk45dp_a = { 0, 0.2, 0.3, 0.8, 8.0 / 9, 1, 1 };
const valarray<vector_type> rk45dp_b = {
    {},
    {0.2},
    {3.0/40, 9.0/40},
    {44.0/45, -56.0/15, 32.0/9},
    {19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729},
    {9017.0/3168, -355.0/33, 46732.0/5247, 49.0/176, -5103.0/18656},
    {35.0/384, 0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84}
};
const vector_type rk45dp_c = { 35.0 / 384, 0, 500.0 / 113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84, 0 };
const vector_type rk45dp_ch = { 5179.0 / 57600, 0, 7571.0 / 16695, 393.0 / 640, -92097.0 / 339200, 187.0 / 2100, 1.0 / 40 };

auto ode_rk45(const valarray<vector_function_type>& differentials,
    const vector_type& initial_conditions,
    double domain_start, double domain_end,
    const custom_parametric_function& pf, int freq) {
    gensolveAB(differentials, initial_conditions, domain_start, domain_end, _rk_odesolver_tolerance, rk45dp_a, rk45dp_b, rk45dp_c, rk45dp_ch, 6, 7, 1e-4, pf, freq);
}

/* Cash-Karp RK45 adaptive step integrator */

const vector_type rk45cp_a = { 0, 0.2, 0.3, 0.6, 1, 7.0 / 8 };
const valarray<vector_type> rk45cp_b = {
    {},
    {0.2},
    {3.0/40, 9.0/40},
    {3.0/10, -9.0/10, 6.0/5.0},
    {-11.0/54, 5.0/2, -70.0/27, 35.0/27},
    {1631.0/55296, 175.0/512, 575.0/13824, 44275.0/110592, 253.0/4096}
};
const vector_type rk45cp_c = { 37.0 / 378, 0, 250.0 / 621, 125.0 / 594, 0, 512.0 / 1771 };
const vector_type rk45cp_ch = { 2825.0 / 27648, 0, 18575.0 / 48384, 13525.0 /55296, 277.0 / 14336, 1.0 / 4 };

auto ode_rk45cp(const valarray<vector_function_type>& differentials,
    const vector_type& initial_conditions,
    double domain_start, double domain_end,
    const custom_parametric_function& pf, int freq) {
    gensolveAB(differentials, initial_conditions, domain_start, domain_end, _rk_odesolver_tolerance, rk45cp_a, rk45cp_b, rk45cp_c, rk45cp_ch, 5, 6, 1e-4, pf, freq);
}

/* Euler-Heun adaptive integrator */

const vector_type rk12eh_a = { 0, 1 };
const valarray<vector_type> rk12eh_b = {
    {},
    {1}
};
const vector_type rk12eh_c = { 0.5, 0.5 };
const vector_type rk12eh_ch = { 1, 0 };

auto ode_rk12eh(const valarray<vector_function_type>& differentials,
    const vector_type& initial_conditions,
    double domain_start, double domain_end,
    const custom_parametric_function& pf, int freq) {
    gensolveAB(differentials, initial_conditions, domain_start, domain_end, _rk_odesolver_tolerance, rk12eh_a, rk12eh_b, rk12eh_c, rk12eh_ch, 1, 2, 1e-4, pf, freq);
}

/* Bogacki-Shampine adaptive integrator */

const vector_type rk23bs_a = { 0, 0.5, 0.75, 1 };
const valarray<vector_type> rk23bs_b = {
    {},
    {0.5},
    {0, 0.75},
    {2.0/9, 1.0/3, 4.0/9}
};
const vector_type rk23bs_ch = { 2.0 / 9, 1.0 / 3, 4.0 / 9 };
const vector_type rk23bs_c = { 7.0/24, 0.25, 1.0/3, 1.0/8 };

auto ode_rk23bs(const valarray<vector_function_type>& differentials,
    const vector_type& initial_conditions,
    double domain_start, double domain_end,
    const custom_parametric_function& pf, int freq) {
    gensolveAB(differentials, initial_conditions, domain_start, domain_end, _rk_odesolver_tolerance, rk23bs_a, rk23bs_b, rk23bs_c, rk23bs_ch, 2, 3, 1e-4, pf, freq);
}

/* Crank Nicholson */

auto ode_cnrk2(const valarray<vector_function_type>& differentials,
    const vector_type& initial_conditions, const double stepping,
    double domain_start, double domain_end,
    const custom_parametric_function& pf, int freq) {
    vector_type state = initial_conditions;
    int count = 0;
    for (double marker = domain_start; marker <= domain_end; marker += stepping) {
        auto k1 = stepping * eval(differentials, state, marker);
        functional_vector ies2(zero, differentials.size());
        for (int i = 0; i < differentials.size(); i++) {
            auto imp_eq = [stepping, i, differentials, state, k1](vector_type vnext, double marker) -> double {
                return vnext[i] - stepping * differentials[i](state + 0.5 * k1 + 0.5 * vnext, marker + stepping);
                };
            ies2[i] = imp_eq;
        }
        auto k2 = newton_raphson(ies2, state, marker, 1e-6, 5);
        state += 0.5 * k1 + 0.5 * k2;
        if(count == freq) {
            auto s = pf(state, marker);
            if (s != 0) break;
            count = 0;
        } else {
            count++;
        }
    }
}

/* Qin and Zhang's two-stage */

auto ode_qzdirk2(const valarray<vector_function_type>& differentials,
    const vector_type& initial_conditions, const double stepping,
    double domain_start, double domain_end,
    const custom_parametric_function& pf, int freq) {
    vector_type state = initial_conditions;
    int count = 0;
    for (double marker = domain_start; marker <= domain_end; marker += stepping) {
        functional_vector ies1(zero, differentials.size());
        for (int i = 0; i < differentials.size(); i++) {
            auto imp_eq = [stepping, i, differentials, state](vector_type vnext, double marker) -> double {
                return vnext[i] - stepping * differentials[i](state + 0.25 * vnext, marker + 0.25 * stepping);
                };
            ies1[i] = imp_eq;
        }
        auto k1 = newton_raphson(ies1, state, marker, 1e-6, 5);
        functional_vector ies2(zero, differentials.size());
        for (int i = 0; i < differentials.size(); i++) {
            auto imp_eq = [stepping, i, differentials, state, k1](vector_type vnext, double marker) -> double {
                return vnext[i] - stepping * differentials[i](state + 0.5 * k1 + 0.25 * vnext, marker + 0.75 * stepping);
                };
            ies2[i] = imp_eq;
        }
        auto k2 = newton_raphson(ies2, state, marker, 1e-6, 5);
        state += 0.5 * k1 + 0.5 * k2;
        if(count == freq) {
            auto s = pf(state, marker);
            if (s != 0) break;
            count = 0;
        } else {
            count++;
        }
    }
}