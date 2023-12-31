/*
 * ID       : EXP RESEARCH LAB CP 23 2
 * Author   : CHLORINE PENTOXIDE
 * DOBE     : 31 12 2023
 * STATUS   : ALIVE
 * OBJECTIVE: GET PRESSURE(STG 1) DATA ONLY
 */

#include <fstream>
#include "../../rocket/models/urm_model.h"

/*
 * Standard Rocket Model
 *  1kg 7atm 75/2mm 35/2mm 1m (1.0*pi)m 1e3kg/m3 0.3
 *
 * Experimental Modifications: DRAG OFF, BYO OFF
 * Run(0:40) PRES(6atm:4atm:0.1atm) PI(0.2:0.8)
 *
 */
int main() {
    using namespace fastsim;
    EnvironmentParameters env;
    int watch_channel = urm::CH_VLCT;
    for(int i = 60, count = 0; i >= 20; i--, count++) {
        std::ofstream out("exp_rlcp231_results/exp_rlcp231_run"+std::to_string(count)+".csv");
        out << "pressure,pi,max_velocity,max_velocity_time" << "\n";
        double pi = 0.0;
        experimental_coprocessor coprocessor = [] (vector_type&, double, int) -> int {
            return 0;
        };
        experimental_coprocessor stage_processor = [&watch_channel, &out, &pi, &i] (vector_type& state, double t, int stage) -> int {
            if(stage == 3) {
                out << (i / 10.0) << "," << pi << ",";
                out << state[watch_channel] << "," << t << "\n";
            }
            return 0;
        };
        for (pi = 0.2; pi <= 0.8; pi += 0.01) {
            std::cout << "Running PRES(" << i / 10.0 << "atm)" << " pi(" << pi << ") ... ";
            RocketParameters rocket = {
                    1,
                    i * env.pressure_atm / 10.0,
                    75.0 / 2000,
                    35.0 / 2000,
                    1.0,
                    1.0 * pi,
                    1e3,
                    0.3
            };
            urm::OptimizedParameters parameters = urm::calculate_optimized_parameters(rocket, env);
            parameters.DRAG = 0;
            parameters.BYO = 0;
            urm::simulate(parameters, 3, coprocessor, stage_processor);
            std::cout << "DONE" << "\n";
        }
        out.close();
    }
}