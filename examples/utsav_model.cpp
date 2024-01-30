/*
 * DOI: TBF
 */

#include "../rocket/models/urm_model.h"

int main(int argc, char* argv[]) {
    using namespace fastsim;
    EnvironmentParameters env;
    RocketParameters rocket = { 1, 3 * env.pressure_atm, 120.0 / 2000, 36.0 / 2000, 0.75, 0.75 * 0.4, 1e3, 0.4};
    return urm::simulate(argc, argv, rocket, env, urm::CH_HGHT);
}