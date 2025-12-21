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
#include "myio.h" /* my input-output routines */
#include "nu5d_mat_osc.h"

#include "T2K_setup.h"

#include <glb_smear.h>
#include <glb_types.h>
#include <globes/globes.h> /* GLoBES library */
#include <gsl/gsl_sf_erf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char MYFILE[] = "../data/prob/T2KSMprob-2-2.dat";
char MYFILE1[] = "../data/T2K/T2K_numu_FHC_2023.dat";
char MYFILE1A[] = "../data/T2K/T2K_numu_FHC_data.dat";
char MYFILE2[] = "../data/T2K/T2K_numu_RHC_2023.dat";
char MYFILE2A[] = "../data/T2K/T2K_numu_RHC_data.dat";
char MYFILE3[] = "../data/T2K/T2K_nue_FHC_2023.dat";
char MYFILE3A[] = "../data/T2K/T2K_nue_FHC_data.dat";
char MYFILE4[] = "../data/T2K/T2K_nue_RHC_2023.dat";
char MYFILE4A[] = "../data/T2K/T2K_nue_RHC_data.dat";
char MYFILE5[] = "../data/T2K/T2K_nueCCpi_2023.dat";
char MYFILE5A[] = "../data/T2K/T2K_nueCCpi_data.dat";

int main(int argc, char* argv[]) {
    /* Initialize libglobes */
    glbInit(argv[0]);
    glbDefineChiFunction(&ChiT2K, 4, "ChiT2K", NULL);
    /* Initialize T2K */
    InitializeT2K(&glb_experiment_list[0], &glb_num_of_exps);
    /* Define standard oscillation parameters for NO in T2K */
    // NH
    double theta12 = asin(sqrt(0.307));   // nu-fit 6.0
    double theta13 = asin(sqrt(0.02195)); // 2303
    double theta23 = asin(sqrt(0.561));   // 2303
    double deltacp = -2.22;               // 2303
    double sdm = 7.49e-5;                 // nu-fit 5.2
    double ldm = 2.495e-3 + sdm;          // 2303

    // // IH
    // double theta12 = asin(sqrt(0.307)); // nu-fit 6.0
    // double theta13 = asin(sqrt(0.02224));
    // double theta23 = asin(sqrt(0.563));
    // double deltacp = -1.44;
    // double sdm = 7.49e-5;    // nu-fit 6.0
    // double ldm = -2.4634e-3; // IO

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
    for (int i = 1; i <= 2000; i++) {
        double prob = glbFilteredConstantDensityProbability(0, 2, 2, -1, i * 0.002);
        AddToOutput2(i * 0.002, prob);
    }

    /* Obtain lists for the energy bins in the considered samples */
    int i;

    double *true_rates_N1, *true_rates_N2, *true_rates_N0, *true_rates_N4, *true_rates_N5;
    double coefficient_N1, coefficient_N2, coefficient_N3, coefficient_N4, coefficient_N5;
    double *numu_bin_centers_E, *numubar_bin_centers_E, *nue_bin_centers_E;
    double *numu_bin_widths, *numubar_bin_widths, *nue_bin_widths;

    int n_bins1 = glbGetNumberOfBins(0); /* Neutrino beam muon-like sample */
    int n_bins2 = glbGetNumberOfBins(1); /* Antineutrino beam muon-like sample */
    int n_bins3 = glbGetNumberOfBins(2); /* Neutrino and antineutrino beam electron-like samples */

    numu_bin_centers_E = glbGetBinCentersListPtr(0);
    numubar_bin_centers_E = glbGetBinCentersListPtr(1);
    nue_bin_centers_E = glbGetBinCentersListPtr(2);

    numu_bin_widths = glbGetBinSizeListPtr(0);
    numubar_bin_widths = glbGetBinSizeListPtr(1);
    nue_bin_widths = glbGetBinSizeListPtr(2);

    /* Compute event rate spectra for the neutrino muon-like sample in Super-Kamiokande */
    InitOutput(MYFILE1, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (i = 0; i <= n_bins1 - 1; i++) {
        AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }
    InitOutput(MYFILE1A, "");
    for (i = 0; i <= n_bins1 - 1; i++) AddToOutput2(numu_bin_centers_E[i], true_rates_N0[i]);
    printf("The events in 1Rmu is total of %f\n", glbTotalRuleRate(0, 0, GLB_ALL, GLB_W_EFF, GLB_WO_BG, GLB_W_COEFF, GLB_SIG));

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    InitOutput(MYFILE2, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (i = 0; i <= n_bins2 - 1; i++) {
        AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }
    InitOutput(MYFILE2A, "");
    for (i = 0; i <= n_bins2 - 1; i++) AddToOutput2(numubar_bin_centers_E[i], true_rates_N0[i]);
    printf("The events in 1Rmu_anti is total of %f\n", glbTotalRuleRate(1, 0, GLB_ALL, GLB_W_EFF, GLB_WO_BG, GLB_W_COEFF, GLB_SIG));

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    InitOutput(MYFILE3, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }
    InitOutput(MYFILE3A, "");
    for (i = 0; i <= n_bins3 - 1; i++) AddToOutput2(nue_bin_centers_E[i], true_rates_N0[i]);
    printf("The events in 1Re is total of %f\n", glbTotalRuleRate(2, 0, GLB_ALL, GLB_W_EFF, GLB_WO_BG, GLB_W_COEFF, GLB_SIG));

    InitOutput(MYFILE4, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }
    InitOutput(MYFILE4A, "");
    for (i = 0; i <= n_bins3 - 1; i++) AddToOutput2(nue_bin_centers_E[i], true_rates_N0[i]);
    printf("The events in 1Re_anti is total of %f\n", glbTotalRuleRate(2, 1, GLB_ALL, GLB_W_EFF, GLB_WO_BG, GLB_W_COEFF, GLB_SIG));

    InitOutput(MYFILE5, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 2);
    true_rates_N1 = glbGetSignalRatePtr(2, 2); // 1Re1de event
    true_rates_N2 = glbGetBGRatePtr(2, 2);
    for (i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }
    InitOutput(MYFILE5A, "");
    for (i = 0; i <= n_bins3 - 1; i++) AddToOutput2(nue_bin_centers_E[i], true_rates_N0[i]);
    printf("The events in 1Re1de_anti is total of %f\n", glbTotalRuleRate(2, 2, GLB_ALL, GLB_W_EFF, GLB_WO_BG, GLB_W_COEFF, GLB_SIG));

    glbFreeParams(true_values);

    /* Clear experiment list */
    glbClearExperimentList();

    exit(0);
}
