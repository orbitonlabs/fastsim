#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

auto dydt(double y, double t) {
	return 2 * t;
}

int e1main() {
	std::cout << "Experiment 1" << "\n";
	std::cout << "Copyright(C) 2023, P. Chlorine" << "\n";
	
	auto h = 1e-4;
	auto der_h = 1e-10;
	auto y0 = 1;

	double y = y0;
	for (double t = 0.0; t <= 1; t += h) {
		auto impliciteq = [h, t, y] (double y_next) -> double{
			return - y_next + dydt(y_next, t + h) * h + y;
		};
		/*std::ofstream output("output_" + std::to_string(t) + ".csv");
		output << "x,f(x)" << "\n";
		for (double ax = y - 2.0; ax <= y + 2.0; ax += 1e-4) {
			output << ax << "," << impliciteq(ax) << "\n";
		}
		output.close();*/
		double y_next = y;
		for (int i = 1; i <= 3; i++) {
			auto der = [der_h, impliciteq] (double y) -> double {
				return (impliciteq(y + der_h) - impliciteq(y)) / der_h;
			};
			y_next -= (impliciteq(y_next) / der(y_next));
		}
		std::cout << y << "\t" << (t * t + 1) << "\n";
		y = y_next;
	}
	return 0;
}