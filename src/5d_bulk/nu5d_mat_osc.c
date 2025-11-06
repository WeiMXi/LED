#include "nu5d_mat_osc.h"

#include "myio.h" /* my input-output routines */

#include <globes/globes.h> /* GLoBES library */
#include <stdlib.h>
#include <string.h>

#define sq(x) ((x) * (x))

char NOvAFILE1[] = "../data/prob/NOvAprob4D.dat";
char NOvAFILE2[] = "../data/prob/NOvAprobcir2.dat";
char NOvAFILE3[] = "../data/prob/NOvAprobcir5.dat";
char NOvAFILE4[] = "../data/prob/NOvAprobmuir0.01.dat";
char NOvAFILE5[] = "../data/prob/NOvAprobmuir0.02.dat";
char T2KFILE1[] = "../data/prob/T2Kprob4D.dat";
char T2KFILE2[] = "../data/prob/T2Kprobcir2.dat";
char T2KFILE3[] = "../data/prob/T2Kprobcir5.dat";
char T2KFILE4[] = "../data/prob/T2Kprobmuir0.01.dat";
char T2KFILE5[] = "../data/prob/T2Kprobmuir0.02.dat";

int main(int argc, char* argv[]) {
    /* Initialize libglobes */
    glbInit(argv[0]);
    glbClearExperimentList();
    glbInitExperiment("test.glb", &glb_experiment_list[0],
                      &glb_num_of_exps); // this glb has 2 sampling points to make calculation faster, so this file can only be used to calculate probability.

    /* Initialize T2K*/

    glbRegisterProbabilityEngine(13, /*Number of parameters*/
                                 &my_probability_matrix,
                                 &my_set_oscillation_parameters,
                                 &my_get_oscillation_parameters,
                                 NULL);
    /* Define central values for the prior function (adopted from neutrino 2022 conference results) */
    double theta12 = 33.68 / 180;
    double theta13 = 8.52 / 180;
    double theta23 = 48.5 / 180;
    double deltacp = 0;
    double sdm = DMSQ21;
    double ldm = DMSQ31_NH;

    glb_params central_values = glbAllocParams();
    glbDefineParams(central_values, theta12, theta13, theta23, deltacp, sdm, ldm);

    glbSetOscParams(central_values, 10, GLB_R);
    glbSetOscParams(central_values, 40, GLB_C1R);
    glbSetOscParams(central_values, -40, GLB_C2R);
    glbSetOscParams(central_values, -40, GLB_C3R);
    glbSetOscParams(central_values, 0.01, GLB_MU1R);
    glbSetOscParams(central_values, 0.0275, GLB_MU2R);
    glbSetOscParams(central_values, 0.1603, GLB_MU3R);
    SetModesCutoff(50);

    glbSetOscillationParameters(central_values);
    glbSetRates();

    double baseline = 295; // 295 for T2K and 810 for NOvA

    InitOutput(T2KFILE1, "");
    for (int i = 1; i <= 1000; i++) {
        double prob = glbVacuumProbability(2, 2, 1, i * 0.004, baseline);
        AddToOutput2(i * 0.004, prob);
    }
    glbSetOscParams(central_values, 42, GLB_C1R);
    glbSetOscParams(central_values, -42, GLB_C2R);
    glbSetOscParams(central_values, -42, GLB_C3R);
    glbSetOscillationParameters(central_values);
    glbSetRates();
    InitOutput(T2KFILE2, "");
    for (int i = 1; i <= 1000; i++) {
        double prob = glbVacuumProbability(2, 2, 1, i * 0.004, baseline);
        AddToOutput2(i * 0.004, prob);
    }

    glbSetOscParams(central_values, 45, GLB_C1R);
    glbSetOscParams(central_values, -45, GLB_C2R);
    glbSetOscParams(central_values, -45, GLB_C3R);
    glbSetOscillationParameters(central_values);
    glbSetRates();
    InitOutput(T2KFILE3, "");
    for (int i = 1; i <= 1000; i++) {
        double prob = glbVacuumProbability(2, 2, 1, i * 0.004, baseline);
        AddToOutput2(i * 0.004, prob);
    }

    glbSetOscParams(central_values, 40, GLB_C1R);
    glbSetOscParams(central_values, -40, GLB_C2R);
    glbSetOscParams(central_values, -40, GLB_C3R);
    glbSetOscParams(central_values, 0.01 + 0.01, GLB_MU1R);
    glbSetOscParams(central_values, 0.01 + 0.0275, GLB_MU2R);
    glbSetOscParams(central_values, 0.01 + 0.1603, GLB_MU3R);
    glbSetOscillationParameters(central_values);
    glbSetRates();
    InitOutput(T2KFILE4, "");
    for (int i = 1; i <= 1000; i++) {
        double prob = glbVacuumProbability(2, 2, 1, i * 0.004, baseline);
        AddToOutput2(i * 0.004, prob);
    }

    glbSetOscParams(central_values, 0.02 + 0.01, GLB_MU1R);
    glbSetOscParams(central_values, 0.02 + 0.0275, GLB_MU2R);
    glbSetOscParams(central_values, 0.02 + 0.1603, GLB_MU3R);
    glbSetOscillationParameters(central_values);
    glbSetRates();
    InitOutput(T2KFILE5, "");
    for (int i = 1; i <= 1000; i++) {
        double prob = glbVacuumProbability(2, 2, 1, i * 0.004, baseline);
        AddToOutput2(i * 0.004, prob);
    }

    return EXIT_SUCCESS;
}
