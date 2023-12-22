#pragma once

#include "../lib/types.h"
#include "../lib/jacobian.h"
#include <cmath>

static double newton_raphson(std::function<double(double)> f, double g, double tol, int depth) {
	double x = g;
	auto der = [f](double x) -> double {
		return (f(x + linalg::delta) - f(x)) / linalg::delta;
	};
	for (int i = 1; i <= depth; i++) {
		x -= f(x) / der(x);
		if (fabs(x) < tol) return x;
	}
	return x;
}

static vector_type newton_raphson(functional_vector f, vector_type g, double t, double tol, int depth) {
	vector_type root = g;
	for (int i = 1; i <= depth; i++) {
		auto jcb = linalg::compute_jacobian(f, root, t);
		vector_type fval = eval(f, root, t);
		matrix::Matrix jinv(1, 1);
		matrix::Matrix::geninv(jcb, jinv);
		auto dch = jinv.dot(fval);
		if (sqrt(vector::dot(dch, dch)) < tol) break;
		root = root - dch;
	}
	return root;
}