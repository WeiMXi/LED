#include "ledlib/Engine/ProbabilityEngine.h++"
#include "ledlib/IO/IO.h++" /* my input-output routines */

#include <cstdlib>
#include <globes/globes.h> /* GLoBES library */
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[]) {

    /* output files path*/
    const std::string NOvAFile1{"../data/prob/NOvAprob4D.dat"};
    const std::string NOvAFile2{"../data/prob/NOvAprobcir2.dat"};
    const std::string NOvAFile3{"../data/prob/NOvAprobcir5.dat"};
    const std::string NOvAFile4{"../data/prob/NOvAprobmuir0.01.dat"};
    const std::string NOvAFile5{"../data/prob/NOvAprobmuir0.02.dat"};
    const std::string T2KFile1{"../data/prob/T2Kprob4D.dat"};
    const std::string T2KFile2{"../data/prob/T2Kprobcir2.dat"};
    const std::string T2KFile3{"../data/prob/T2Kprobcir5.dat"};
    const std::string T2KFile4{"../data/prob/T2Kprobmuir0.01.dat"};
    const std::string T2KFile5{"../data/prob/T2Kprobmuir0.02.dat"};

    LED::IO::Output outputFiles;

    /* Initialize libglobes */
    glbInit(argv[0]);
    glbClearExperimentList();
    glbInitExperiment("test.glb", &glb_experiment_list[0],
                      &glb_num_of_exps); // this glb has 2 sampling points to make calculation faster, so this file can only be used to calculate probability.

    glbRegisterProbabilityEngine(13, /*Number of parameters*/
                                 &LED::CalProbability::my_probability_matrix,
                                 &LED::CalProbability::my_set_oscillation_parameters,
                                 &LED::CalProbability::my_get_oscillation_parameters,
                                 NULL);

    // /*parameters for T2K*/
    // const double theta12 = asin(sqrt(0.303));
    // const double theta13 = asin(sqrt(0.028)); // 2303.03222
    // const double theta23 = asin(sqrt(0.467)); // 2303.03222
    // const double deltacp = -2.22;             // 2303.03222
    // const double sdm = 7.53e-5;
    // const double ldm = 2.494e-3 + sdm; // 2303.03222

    /*parameters for NOvA*/
    double theta12 = asin(sqrt(0.307)); // 2108
    double theta13 = asin(sqrt(0.021)); // 2108
    double theta23 = asin(sqrt(0.57));  // 2108 result
    double deltacp = 0.82 * M_PI;       // 2108 result
    double sdm = 7.53e-5;               // 2108
    double ldm = 2.41e-3 + sdm;         // 2108 result

    glb_params central_values = glbAllocParams();
    glbDefineParams(central_values, theta12, theta13, theta23, deltacp, sdm, ldm);

    glbSetOscParams(central_values, 10, LED::CalProbability::GLB_R);
    glbSetOscParams(central_values, 4, LED::CalProbability::GLB_C1R);
    glbSetOscParams(central_values, -4, LED::CalProbability::GLB_C2R);
    glbSetOscParams(central_values, -4, LED::CalProbability::GLB_C3R);
    // glbSetOscParams(central_values, 0.01, LED::CalProbability::GLB_MU1R);
    // glbSetOscParams(central_values, 0.0275, LED::CalProbability::GLB_MU2R);
    // glbSetOscParams(central_values, 0.1603, LED::CalProbability::GLB_MU3R); // T2K
    // glbSetOscParams(central_values, 0.1598, LED::CalProbability::GLB_MU3R); // NOvA
    glbSetOscParams(central_values, 2, LED::CalProbability::GLB_MU1R);
    // std::cout << LED::CalProbability::CalLightestm2(10, 4, 0.5) << std::endl;
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(10, -4, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(10, -4, ldm), LED::CalProbability::GLB_MU3R);
    LED::CalProbability::SetModesCutoff(20);

    glbSetOscillationParameters(central_values);
    glbSetRates();

    const double baseline = 810; // 295 for T2K and 810 for NOvA

    outputFiles.InitOutput(NOvAFile1, "");
    for (int i = 1; i <= 2000; i++) {
        const double prob = glbFilteredConstantDensityProbability(0, 2, 1, 1, i * 0.002);
        outputFiles.AddToOutput2(i * 0.002, prob);
    }
    // glbSetOscParams(central_values, 42, LED::CalProbability::GLB_C1R);
    // glbSetOscParams(central_values, -42, LED::CalProbability::GLB_C2R);
    // glbSetOscParams(central_values, -42, LED::CalProbability::GLB_C3R);
    // glbSetOscillationParameters(central_values);
    // glbSetRates();
    // outputFiles.InitOutput(T2KFile2, "");
    // for (int i = 1; i <= 1000; i++) {
    //     const double prob = glbConstantDensityProbability(2, 2, 1, i * 0.004, baseline, 2.8);
    //     outputFiles.AddToOutput2(i * 0.004, prob);
    // }

    // glbSetOscParams(central_values, 45, LED::CalProbability::GLB_C1R);
    // glbSetOscParams(central_values, -45, LED::CalProbability::GLB_C2R);
    // glbSetOscParams(central_values, -45, LED::CalProbability::GLB_C3R);
    // glbSetOscillationParameters(central_values);
    // glbSetRates();
    // outputFiles.InitOutput(T2KFile3, "");
    // for (int i = 1; i <= 1000; i++) {
    //     const double prob = glbConstantDensityProbability(2, 2, 1, i * 0.004, baseline, 2.8);
    //     outputFiles.AddToOutput2(i * 0.004, prob);
    // }

    // glbSetOscParams(central_values, 40, LED::CalProbability::GLB_C1R);
    // glbSetOscParams(central_values, -40, LED::CalProbability::GLB_C2R);
    // glbSetOscParams(central_values, -40, LED::CalProbability::GLB_C3R);
    // glbSetOscParams(central_values, 0.01 + 0.01, LED::CalProbability::GLB_MU1R);
    // glbSetOscParams(central_values, 0.01 + 0.0275, LED::CalProbability::GLB_MU2R);
    // glbSetOscParams(central_values, 0.01 + 0.1603, LED::CalProbability::GLB_MU3R);
    // glbSetOscillationParameters(central_values);
    // glbSetRates();
    // outputFiles.InitOutput(T2KFile4, "");
    // for (int i = 1; i <= 1000; i++) {
    //     const double prob = glbConstantDensityProbability(2, 2, 1, i * 0.004, baseline, 2.8);
    //     outputFiles.AddToOutput2(i * 0.004, prob);
    // }

    // glbSetOscParams(central_values, 0.02 + 0.01, LED::CalProbability::GLB_MU1R);
    // glbSetOscParams(central_values, 0.02 + 0.0275, LED::CalProbability::GLB_MU2R);
    // glbSetOscParams(central_values, 0.02 + 0.1603, LED::CalProbability::GLB_MU3R);
    // glbSetOscillationParameters(central_values);
    // glbSetRates();
    // outputFiles.InitOutput(T2KFile5, "");
    // for (int i = 1; i <= 1000; i++) {
    //     const double prob = glbConstantDensityProbability(2, 2, 1, i * 0.004, baseline, 2.8);
    //     outputFiles.AddToOutput2(i * 0.004, prob);
    // }

    return EXIT_SUCCESS;
}
