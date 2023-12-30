/*
 * DOI: TBF
 */

#include "urm_model.h"

int main(int argc, char* argv[]) {
    using namespace fastsim;
    EnvironmentParameters env;
    RocketParameters rocket = { 1, 8 * env.pressure_atm, 75.0 / 2000, 35.0 / 2000, 1.0, 0.33, 1e3, 0.3};
    return urm::simulate(argc, argv, rocket, env, urm::CH_HGHT);
}