/*
 * Notes: MF doesn't even have a proper paper, missing expression for ve
 * DOI: 10.1088/1742-6596/2508/1/012053
 */

#include <cmath>
#include "../../lib/types.h"

const double propellant_density = 1e3;
const double drag_cd = 0.35;
const double gravity = 9.8;
const double air_density = 1.293;
const double pressure_atm = 101325;

struct RocketModel {
    double mb;
    double P0;
    double Rp;
    double Re;
    double height_chamber;
    double L0;
};

double area(double radius) {
    return M_PI * pow(radius, 2);
}

void run_simulation(RocketModel rm) {
    const double vc = area(rm.Rp) * rm.height_chamber;
    const double ma = air_density * area(rm.Rp) * (rm.height_chamber - rm.L0);
    double vw = area(rm.Rp) * rm.L0;

}

int main(int argc, char* argv[]) {

    RocketModel rm = { 0.5, 4 * pressure_atm, 75.0e-3 / 2, 35.0e-3 / 2, 1, 0.3 };

}