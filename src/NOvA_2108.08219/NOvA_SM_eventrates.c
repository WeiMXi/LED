/* This program produces the data needed to plot the total events in the T2K far detector.
 * The output consists of five different spectra: CCQE nu_mu, CCQE nu_mu_bar, CCQE nu_e,
 * CCQE nu_e_bar and CC nu_e pi^+.
 *
 *
 * To do list:
 *   Step 1: Uncalibrated events with builtin smearing\\completed
 *   Step 2: Events after correction to normalization \\ignored
 *   Step 3: Events after smearing corrections\\no change
 *   Step 4: Events after data-driven corrections
 */

#include "NOvA_setup.h"
#include "myio.h" /* my input-output routines */

#include <glb_smear.h>
#include <glb_types.h>
#include <globes/globes.h> /* GLoBES library */
#include <gsl/gsl_sf_erf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char MYFILE[] = "../data/prob/NOvASMprob.dat";
char MYFILE1[] = "../data/NOvA/NOvA_numu_2021.dat";
char MYFILE1A[] = "../data/NOvA/NOvA_numu_data.dat";
char MYFILE2[] = "../data/NOvA/NOvA_anumu_2021.dat";
char MYFILE2A[] = "../data/NOvA/NOvA_anumu_data.dat";

char MYFILE3[] = "../data/NOvA/NOvA_nue_LCNN_2021.dat";
char MYFILE3A[] = "../data/NOvA/NOvA_nue_LCNN_data.dat";
char MYFILE4[] = "../data/NOvA/NOvA_nue_HCNN_2021.dat";
char MYFILE4A[] = "../data/NOvA/NOvA_nue_HCNN_data.dat";

char MYFILE5[] = "../data/NOvA/NOvA_anue_LCNN_2021.dat";
char MYFILE5A[] = "../data/NOvA/NOvA_anue_LCNN_data.dat";
char MYFILE6[] = "../data/NOvA/NOvA_anue_HCNN_2021.dat";
char MYFILE6A[] = "../data/NOvA/NOvA_anue_HCNN_data.dat";

int main(int argc, char* argv[]) {
    /* Initialize libglobes */
    glbInit(argv[0]);
    /* Initialize NOvA */
    InitializeNOvA(&glb_experiment_list[0], &glb_num_of_exps);

    /* Define standard oscillation parameters for NO in NOvA */
    double theta12 = asin(sqrt(0.307)); // NOvA
    double theta13 = asin(sqrt(0.021)); // NOvA
    double theta23 = asin(sqrt(0.57));  // NOvA
    double deltacp = 0.82 * M_PI;       // NOvA
    double sdm = 7.53e-5;               // NOvA
    double ldm = 2.41e-3 + sdm;         // NOvA

    /* Initialize the parameter vector */
    glb_params true_values = glbAllocParams();

    /* Set the parameter vector */
    glbDefineParams(true_values, theta12, theta13, theta23, deltacp, sdm, ldm);
    glbSetDensityParams(true_values, 1.0, GLB_ALL);

    /* Set central values */
    glbSetCentralValues(true_values);

    /* The oscillation probabilities are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();

    InitOutput(MYFILE, "");
    for (int i = 1; i <= 1000; i++) {
        double prob = glbConstantDensityProbability(2, 2, 1, i * 0.004, 810, 2.8);
        AddToOutput2(i * 0.004, prob);
    }

    /* Obtain lists for the energy bins in the considered samples */

    double *true_rates_N1, *true_rates_N2, *true_rates_N0, *true_rates_N4, *true_rates_N5;
    double coefficient_N1, coefficient_N2, coefficient_N3, coefficient_N4, coefficient_N5;
    double *numu_bin_centers_E, *numubar_bin_centers_E, *nue_bin_centers_E;
    double *numu_bin_widths, *numubar_bin_widths, *nue_bin_widths;

    int n_bins1 = glbGetNumberOfBins(0);
    int n_bins2 = glbGetNumberOfBins(1);

    numu_bin_centers_E = glbGetBinCentersListPtr(0);
    numubar_bin_centers_E = glbGetBinCentersListPtr(0);
    nue_bin_centers_E = glbGetBinCentersListPtr(1);

    numu_bin_widths = glbGetBinSizeListPtr(0);
    numubar_bin_widths = glbGetBinSizeListPtr(0);
    nue_bin_widths = glbGetBinSizeListPtr(1);

    /* Compute event rate spectra for the neutrino muon-like event */
    InitOutput(MYFILE1, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0);
    true_rates_N2 = glbGetBGRatePtr(0, 0);

    for (int i = 0; i <= n_bins1 - 1; i++) {
        AddToOutput(numu_bin_centers_E[i], (true_rates_N1[i] + true_rates_N2[i]) / numu_bin_widths[i] * 0.1, numu_bin_widths[i]);
    }
    InitOutput(MYFILE1A, "");
    for (int i = 0; i <= n_bins1 - 1; i++) AddToOutput2(numu_bin_centers_E[i], true_rates_N0[i]);

    /* Compute event rate spectra for the antineutrino muon-like event */
    InitOutput(MYFILE2, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 1);
    true_rates_N1 = glbGetSignalRatePtr(0, 1);
    true_rates_N2 = glbGetBGRatePtr(0, 1);

    for (int i = 0; i <= n_bins1 - 1; i++) {
        AddToOutput(numubar_bin_centers_E[i], (true_rates_N1[i] + true_rates_N2[i]) / numu_bin_widths[i] * 0.1, numubar_bin_widths[i]);
    }
    InitOutput(MYFILE2A, "");
    for (int i = 0; i <= n_bins1 - 1; i++) AddToOutput2(numubar_bin_centers_E[i], true_rates_N0[i]);

    /* Compute event rate spectra for the neutrino electron-like samples with low CNN */
    InitOutput(MYFILE3, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0);
    true_rates_N2 = glbGetBGRatePtr(1, 0);

    for (int i = 0; i <= n_bins2 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }
    InitOutput(MYFILE3A, "");
    for (int i = 0; i <= n_bins2 - 1; i++) AddToOutput2(nue_bin_centers_E[i], true_rates_N0[i]);

    /* Compute event rate spectra for the neutrino electron-like samples with high CNN */
    InitOutput(MYFILE4, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 1);
    true_rates_N1 = glbGetSignalRatePtr(1, 1);
    true_rates_N2 = glbGetBGRatePtr(1, 1);

    for (int i = 0; i <= n_bins2 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }
    InitOutput(MYFILE4A, "");
    for (int i = 0; i <= n_bins2 - 1; i++) AddToOutput2(nue_bin_centers_E[i], true_rates_N0[i]);

    /* Compute event rate spectra for the neutrino antielectron-like samples with low CNN */
    InitOutput(MYFILE5, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 2);
    true_rates_N1 = glbGetSignalRatePtr(1, 2);
    true_rates_N2 = glbGetBGRatePtr(1, 2);

    for (int i = 0; i <= n_bins2 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }
    InitOutput(MYFILE5A, "");
    for (int i = 0; i <= n_bins2 - 1; i++) AddToOutput2(nue_bin_centers_E[i], true_rates_N0[i]);

    /* Compute event rate spectra for the neutrino electron-like samples with high CNN */
    double NOvA_anue_HCNN_MC[6] = {
        1.48931,
        7.96550,
        6.62295,
        2.40363,
        1.06121,
        0.74725};
    InitOutput(MYFILE6, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 3);
    true_rates_N1 = glbGetSignalRatePtr(1, 3);
    true_rates_N2 = glbGetBGRatePtr(1, 3);

    for (int i = 0; i <= n_bins2 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }
    InitOutput(MYFILE6A, "");
    for (int i = 0; i <= n_bins2 - 1; i++) AddToOutput2(nue_bin_centers_E[i], true_rates_N0[i]);

    glbFreeParams(true_values);

    /* Clear experiment list */
    glbClearExperimentList();

    exit(0);
}
