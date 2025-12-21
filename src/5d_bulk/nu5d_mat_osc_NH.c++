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
    const std::string NOvAFile1{"../data/prob/NOvAprob4D-2-2.dat"};
    const std::string NOvAFile2{"../data/prob/NOvAprobcir4-2-2.dat"};
    const std::string NOvAFile3{"../data/prob/NOvAprobcir8-2-2.dat"};
    const std::string NOvAFile4{"../data/prob/NOvAprobmuir1-2-2.dat"};
    const std::string NOvAFile5{"../data/prob/NOvAprobmuir1.5-2-2.dat"};
    const std::string T2KFile1{"../data/prob/T2Kprob4D+2-1.dat"};
    const std::string T2KFile2{"../data/prob/T2Kprobcir4+2-1.dat"};
    const std::string T2KFile3{"../data/prob/T2Kprobcir8+2-1.dat"};
    const std::string T2KFile4{"../data/prob/T2Kprobmuir1+2-1.dat"};
    const std::string T2KFile5{"../data/prob/T2Kprobmuir1.5+2-1.dat"};

    LED::IO::Output outputFiles;

    /* Initialize libglobes */
    glbInit(argv[0]);
    glbClearExperimentList();
    glbInitExperiment("T2K_Prob.glb", &glb_experiment_list[0],
                      &glb_num_of_exps);
    glbInitExperiment("NOvA_Prob.glb", &glb_experiment_list[0],
                      &glb_num_of_exps);

    glbRegisterProbabilityEngine(13, /*Number of parameters*/
                                 &LED::CalProbability::my_probability_matrix,
                                 &LED::CalProbability::my_set_oscillation_parameters,
                                 &LED::CalProbability::my_get_oscillation_parameters,
                                 NULL);

    // /*parameters for T2K*/
    double theta12 = asin(sqrt(0.307)); // nu-fit 6.0
    double theta13 = asin(sqrt(0.02195));
    double theta23 = asin(sqrt(0.561));
    double deltacp = -1.97;
    double sdm = 7.49e-5;        // nu-fit 5.2
    double ldm = 2.494e-3 + sdm; // NO

    /*parameters for NOvA*/
    // double theta12 = asin(sqrt(0.307));   // 2108
    // double theta13 = asin(sqrt(0.02195)); // 2108
    // double theta23 = asin(sqrt(0.57));    // 2108 result
    // double deltacp = 0.82 * M_PI;         // 2108 result
    // double sdm = 7.49e-5;                 // 2108
    // double ldm = 2.41e-3 + sdm;           // 2108 result

    glb_params central_values = glbAllocParams();
    glbDefineParams(central_values, theta12, theta13, theta23, deltacp, sdm, ldm);

    glbSetOscParams(central_values, 1, LED::CalProbability::GLB_R);
    glbSetOscParams(central_values, 10, LED::CalProbability::GLB_C1R);
    glbSetOscParams(central_values, -10, LED::CalProbability::GLB_C2R);
    glbSetOscParams(central_values, -10, LED::CalProbability::GLB_C3R);
    glbSetOscParams(central_values, 0.1, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(1, -10, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(1, -10, ldm), LED::CalProbability::GLB_MU3R);
    LED::CalProbability::SetModesCutoff(20);

    glbSetOscillationParameters(central_values);
    glbSetRates();

    const double baseline = 810; // 295 for T2K and 810 for NOvA

    outputFiles.InitOutput(T2KFile1, "");
    for (int i = 1; i <= 2000; i++) {
        const double prob = glbFilteredConstantDensityProbability(0, 2, 1, 1, i * 0.002);
        outputFiles.AddToOutput2(i * 0.002, prob);
    }

    glbSetOscParams(central_values, 10, LED::CalProbability::GLB_R);
    glbSetOscParams(central_values, 4, LED::CalProbability::GLB_C1R);
    glbSetOscParams(central_values, -4, LED::CalProbability::GLB_C2R);
    glbSetOscParams(central_values, -4, LED::CalProbability::GLB_C3R);
    glbSetOscParams(central_values, 0.1, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(10, -4, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(10, -4, ldm), LED::CalProbability::GLB_MU3R);
    glbSetOscillationParameters(central_values);
    glbSetRates();
    outputFiles.InitOutput(T2KFile2, "");
    for (int i = 1; i <= 2000; i++) {
        const double prob = glbFilteredConstantDensityProbability(0, 2, 1, 1, i * 0.002);
        outputFiles.AddToOutput2(i * 0.002, prob);
    }

    glbSetOscParams(central_values, 10, LED::CalProbability::GLB_R);
    glbSetOscParams(central_values, 8, LED::CalProbability::GLB_C1R);
    glbSetOscParams(central_values, -8, LED::CalProbability::GLB_C2R);
    glbSetOscParams(central_values, -8, LED::CalProbability::GLB_C3R);
    glbSetOscParams(central_values, 0.1, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(10, -8, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(10, -8, ldm), LED::CalProbability::GLB_MU3R);
    glbSetOscillationParameters(central_values);
    glbSetRates();
    outputFiles.InitOutput(T2KFile3, "");
    for (int i = 1; i <= 2000; i++) {
        const double prob = glbFilteredConstantDensityProbability(0, 2, 1, 1, i * 0.002);
        outputFiles.AddToOutput2(i * 0.002, prob);
    }

    glbSetOscParams(central_values, 10, LED::CalProbability::GLB_R);
    glbSetOscParams(central_values, 10, LED::CalProbability::GLB_C1R);
    glbSetOscParams(central_values, -10, LED::CalProbability::GLB_C2R);
    glbSetOscParams(central_values, -10, LED::CalProbability::GLB_C3R);
    glbSetOscParams(central_values, 1, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(10, -10, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(10, -10, ldm), LED::CalProbability::GLB_MU3R);
    glbSetOscillationParameters(central_values);
    glbSetRates();
    outputFiles.InitOutput(T2KFile4, "");
    for (int i = 1; i <= 2000; i++) {
        const double prob = glbFilteredConstantDensityProbability(0, 2, 1, 1, i * 0.002);
        outputFiles.AddToOutput2(i * 0.002, prob);
    }

    glbSetOscParams(central_values, 10, LED::CalProbability::GLB_R);
    glbSetOscParams(central_values, 10, LED::CalProbability::GLB_C1R);
    glbSetOscParams(central_values, -10, LED::CalProbability::GLB_C2R);
    glbSetOscParams(central_values, -10, LED::CalProbability::GLB_C3R);
    glbSetOscParams(central_values, 1.5, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(10, -10, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(10, -10, ldm), LED::CalProbability::GLB_MU3R);
    glbSetOscillationParameters(central_values);
    glbSetRates();
    outputFiles.InitOutput(T2KFile5, "");
    for (int i = 1; i <= 2000; i++) {
        const double prob = glbFilteredConstantDensityProbability(0, 2, 1, 1, i * 0.002);
        outputFiles.AddToOutput2(i * 0.002, prob);
    }

    return EXIT_SUCCESS;
}
