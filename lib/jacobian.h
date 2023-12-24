#pragma once

#include "types.h"

namespace linalg
{
	static double delta = 1e-6;

	static vector_type grad(const vector_function_type& f, const vector_type& x, double t) {
		vector_type g(0.0, x.size());
		for (size_t i = 0; i < x.size(); i++) {
			vector_type dx(0.0, x.size());
			dx[i] = delta;
			g[i] = (f(x + dx, t) - f(x, t)) / delta;
		}
		return g;
	}

	static matrix::Matrix compute_jacobian(const functional_vector& f, const vector_type& x, double t) {
		auto r = f.size(), c = x.size();
		matrix_type m({0.0}, r);
		for (size_t i = 0; i < r; i++) {
			m[i] = grad(f[i], x, t);
		}
		return {m, r, c};
	}
}