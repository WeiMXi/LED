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

char MYFILE1_4Dfinite[] = "../data/T2K/LED/T2K_numu_FHC_4Dfinite.dat";
char MYFILE2_4Dfinite[] = "../data/T2K/LED/T2K_numu_RHC_4Dfinite.dat";
char MYFILE3_4Dfinite[] = "../data/T2K/LED/T2K_nue_FHC_4Dfinite.dat";
char MYFILE4_4Dfinite[] = "../data/T2K/LED/T2K_nue_RHC_4Dfinite.dat";
char MYFILE5_4Dfinite[] = "../data/T2K/LED/T2K_nueCCpi_4Dfinite.dat";

char MYFILE1_probcir2[] = "../data/T2K/LED/T2K_numu_FHC_cir2.dat";
char MYFILE2_probcir2[] = "../data/T2K/LED/T2K_numu_RHC_cir2.dat";
char MYFILE3_probcir2[] = "../data/T2K/LED/T2K_nue_FHC_cir2.dat";
char MYFILE4_probcir2[] = "../data/T2K/LED/T2K_nue_RHC_cir2.dat";
char MYFILE5_probcir2[] = "../data/T2K/LED/T2K_nueCCpi_cir2.dat";

char MYFILE1_probcir5[] = "../data/T2K/LED/T2K_numu_FHC_cir5.dat";
char MYFILE2_probcir5[] = "../data/T2K/LED/T2K_numu_RHC_cir5.dat";
char MYFILE3_probcir5[] = "../data/T2K/LED/T2K_nue_FHC_cir5.dat";
char MYFILE4_probcir5[] = "../data/T2K/LED/T2K_nue_RHC_cir5.dat";
char MYFILE5_probcir5[] = "../data/T2K/LED/T2K_nueCCpi_cir5.dat";

char MYFILE1_muir1[] = "../data/T2K/LED/T2K_numu_FHC_muir0.01.dat";
char MYFILE2_muir1[] = "../data/T2K/LED/T2K_numu_RHC_muir0.01.dat";
char MYFILE3_muir1[] = "../data/T2K/LED/T2K_nue_FHC_ muir0.01.dat";
char MYFILE4_muir1[] = "../data/T2K/LED/T2K_nue_RHC_ muir0.01.dat";
char MYFILE5_muir1[] = "../data/T2K/LED/T2K_nueCCpi_ muir0.01.dat";

char MYFILE1_muir2[] = "../data/T2K/LED/T2K_numu_FHC_muir0.02.dat";
char MYFILE2_muir2[] = "../data/T2K/LED/T2K_numu_RHC_muir0.02.dat";
char MYFILE3_muir2[] = "../data/T2K/LED/T2K_nue_FHC_muir0.02.dat";
char MYFILE4_muir2[] = "../data/T2K/LED/T2K_nue_RHC_muir0.02.dat";
char MYFILE5_muir2[] = "../data/T2K/LED/T2K_nueCCpi_muir0.02.dat";

int main(int argc, char* argv[]) {
    /* Initialize libglobes */
    glbInit(argv[0]);
    glbDefineChiFunction(&ChiT2K, 4, "ChiT2K", NULL);
    /* Initialize T2K */
    InitializeT2K(&glb_experiment_list[0], &glb_num_of_exps);
    glbRegisterProbabilityEngine(13, /*Number of parameters*/
                                 &my_probability_matrix,
                                 &my_set_oscillation_parameters,
                                 &my_get_oscillation_parameters,
                                 NULL);

    /* Initialize the parameter vector */
    glb_params true_values = glbAllocParams();

    /* Define standard oscillation parameters for NO in T2K */
    double theta12 = asin(sqrt(0.303)); // nu-fit 5.2
    double theta13 = asin(sqrt(0.028)); // 2303
    double theta23 = asin(sqrt(0.467)); // 2303
    double deltacp = -2.22;             // 2303
    double sdm = 7.53e-5;               // nu-fit 5.2
    double ldm = 2.495e-3 + sdm;        // 2303

    glbSetOscParams(true_values, 10, GLB_R);
    glbSetOscParams(true_values, 40, GLB_C1R);
    glbSetOscParams(true_values, -40, GLB_C2R);
    glbSetOscParams(true_values, -40, GLB_C3R);
    glbSetOscParams(true_values, 0.01, GLB_MU1R);
    glbSetOscParams(true_values, 0.0275, GLB_MU2R);
    glbSetOscParams(true_values, 0.1603, GLB_MU3R); // T2K
    SetModesCutoff(50);

    /* Set the parameter vector */
    glbDefineParams(true_values, theta12, theta13, theta23, deltacp, sdm, ldm);
    glbSetDensityParams(true_values, 1.0, GLB_ALL);

    /* Set central values */
    glbSetCentralValues(true_values);

    /* The oscillation probabilities are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Obtain lists for the energy bins in the considered samples */
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
    InitOutput(MYFILE1_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    InitOutput(MYFILE2_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    InitOutput(MYFILE3_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    InitOutput(MYFILE4_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    InitOutput(MYFILE5_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 2);
    true_rates_N1 = glbGetSignalRatePtr(2, 2); // 1Re1de event
    true_rates_N2 = glbGetBGRatePtr(2, 2);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbSetOscParams(true_values, 42, GLB_C1R);
    glbSetOscParams(true_values, -42, GLB_C2R);
    glbSetOscParams(true_values, -42, GLB_C3R);
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like sample in Super-Kamiokande */
    InitOutput(MYFILE1_probcir2, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    InitOutput(MYFILE2_probcir2, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    InitOutput(MYFILE3_probcir2, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    InitOutput(MYFILE4_probcir2, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    InitOutput(MYFILE5_probcir2, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 2);
    true_rates_N1 = glbGetSignalRatePtr(2, 2); // 1Re1de event
    true_rates_N2 = glbGetBGRatePtr(2, 2);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbSetOscParams(true_values, 45, GLB_C1R);
    glbSetOscParams(true_values, -45, GLB_C2R);
    glbSetOscParams(true_values, -45, GLB_C3R);
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like sample in Super-Kamiokande */
    InitOutput(MYFILE1_probcir5, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    InitOutput(MYFILE2_probcir5, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    InitOutput(MYFILE3_probcir5, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    InitOutput(MYFILE4_probcir5, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    InitOutput(MYFILE5_probcir5, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 2);
    true_rates_N1 = glbGetSignalRatePtr(2, 2); // 1Re1de event
    true_rates_N2 = glbGetBGRatePtr(2, 2);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbSetOscParams(true_values, 40, GLB_C1R);
    glbSetOscParams(true_values, -40, GLB_C2R);
    glbSetOscParams(true_values, -40, GLB_C3R);
    glbSetOscParams(true_values, 0.01 + 0.01, GLB_MU1R);
    glbSetOscParams(true_values, 0.01 + 0.0275, GLB_MU2R);
    glbSetOscParams(true_values, 0.01 + 0.1603, GLB_MU3R);
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like sample in Super-Kamiokande */
    InitOutput(MYFILE1_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    InitOutput(MYFILE2_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    InitOutput(MYFILE3_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    InitOutput(MYFILE4_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    InitOutput(MYFILE5_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 2);
    true_rates_N1 = glbGetSignalRatePtr(2, 2); // 1Re1de event
    true_rates_N2 = glbGetBGRatePtr(2, 2);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbSetOscParams(true_values, 0.02 + 0.01, GLB_MU1R);
    glbSetOscParams(true_values, 0.02 + 0.0275, GLB_MU2R);
    glbSetOscParams(true_values, 0.02 + 0.1603, GLB_MU3R);
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like sample in Super-Kamiokande */
    InitOutput(MYFILE1_muir2, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    InitOutput(MYFILE2_muir2, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    InitOutput(MYFILE3_muir2, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    InitOutput(MYFILE4_muir2, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    InitOutput(MYFILE5_muir2, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 2);
    true_rates_N1 = glbGetSignalRatePtr(2, 2); // 1Re1de event
    true_rates_N2 = glbGetBGRatePtr(2, 2);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbFreeParams(true_values);

    /* Clear experiment list */
    glbClearExperimentList();

    exit(0);
}
