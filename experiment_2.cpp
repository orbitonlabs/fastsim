#include "lib/types.h"
#include "solvers/newton_raphson.h"

#include <iostream>

double fe11(vector_type v, double t) {
	return (v[0] - 1.1);
}

double fe12(vector_type v, double t) {
	return (v[1] - 1.99);
}

int e2main() {
	std::cout << "Experiment 2" << "\n";
	std::cout << "Copyright(C) 2023, P. Chlorine" << "\n";

	functional_vector fv = { fe11, fe12 };
	vector_type guess = { 1.7, 2 };
	vector_type root = newton_raphson(fv, guess, 2, 1e-6, 3);
	vector::print(root);
	return 0;
}