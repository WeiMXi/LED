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
extern "C" {
#include "T2K_setup.h"
}
#include "ledlib/Engine/ProbabilityEngine.h++"
#include "ledlib/IO/IO.h++" /* my input-output routines */

#include <glb_smear.h>
#include <glb_types.h>
#include <globes/globes.h> /* GLoBES library */
#include <gsl/gsl_sf_erf.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

std::string MYFILE1_4Dfinite = "../data/T2K/LED/T2K_numu_FHC_4Dfinite.dat";
std::string MYFILE2_4Dfinite = "../data/T2K/LED/T2K_numu_RHC_4Dfinite.dat";
std::string MYFILE3_4Dfinite = "../data/T2K/LED/T2K_nue_FHC_4Dfinite.dat";
std::string MYFILE4_4Dfinite = "../data/T2K/LED/T2K_nue_RHC_4Dfinite.dat";
std::string MYFILE5_4Dfinite = "../data/T2K/LED/T2K_nueCCpi_4Dfinite.dat";

std::string MYFILE1_cirp5 = "../data/T2K/LED/T2K_numu_FHC_cir+5.dat";
std::string MYFILE2_cirp5 = "../data/T2K/LED/T2K_numu_RHC_cir+5.dat";
std::string MYFILE3_cirp5 = "../data/T2K/LED/T2K_nue_FHC_cir+5.dat";
std::string MYFILE4_cirp5 = "../data/T2K/LED/T2K_nue_RHC_cir+5.dat";
std::string MYFILE5_cirp5 = "../data/T2K/LED/T2K_nueCCpi_cir+5.dat";

std::string MYFILE1_cirn5 = "../data/T2K/LED/T2K_numu_FHC_cir-5.dat";
std::string MYFILE2_cirn5 = "../data/T2K/LED/T2K_numu_RHC_cir-5.dat";
std::string MYFILE3_cirn5 = "../data/T2K/LED/T2K_nue_FHC_cir-5.dat";
std::string MYFILE4_cirn5 = "../data/T2K/LED/T2K_nue_RHC_cir-5.dat";
std::string MYFILE5_cirn5 = "../data/T2K/LED/T2K_nueCCpi_cir-5.dat";

std::string MYFILE1_muirp5 = "../data/T2K/LED/T2K_numu_FHC_muir+0.005.dat";
std::string MYFILE2_muirp5 = "../data/T2K/LED/T2K_numu_RHC_muir+0.005.dat";
std::string MYFILE3_muirp5 = "../data/T2K/LED/T2K_nue_FHC_muir+0.005.dat";
std::string MYFILE4_muirp5 = "../data/T2K/LED/T2K_nue_RHC_muir+0.005.dat";
std::string MYFILE5_muirp5 = "../data/T2K/LED/T2K_nueCCpi_muir+0.005.dat";

std::string MYFILE1_muirn5 = "../data/T2K/LED/T2K_numu_FHC_muir-0.005.dat";
std::string MYFILE2_muirn5 = "../data/T2K/LED/T2K_numu_RHC_muir-0.005.dat";
std::string MYFILE3_muirn5 = "../data/T2K/LED/T2K_nue_FHC_muir-0.005.dat";
std::string MYFILE4_muirn5 = "../data/T2K/LED/T2K_nue_RHC_muir-0.005.dat";
std::string MYFILE5_muirn5 = "../data/T2K/LED/T2K_nueCCpi_muir-0.005.dat";

LED::IO::Output outputFiles;

int main(int argc, char* argv[]) {
    /* Initialize libglobes */
    glbInit(argv[0]);
    glbDefineChiFunction(&ChiT2K, 4, "ChiT2K", NULL);
    /* Initialize T2K */
    InitializeT2K(&glb_experiment_list[0], &glb_num_of_exps);
    glbRegisterProbabilityEngine(13, /*Number of parameters*/
                                 &LED::CalProbability::my_probability_matrix,
                                 &LED::CalProbability::my_set_oscillation_parameters,
                                 &LED::CalProbability::my_get_oscillation_parameters,
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

    glbSetOscParams(true_values, 1, LED::CalProbability::GLB_R);
    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_C1R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C2R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C3R);
    glbSetOscParams(true_values, 0.1, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(true_values, LED::CalProbability::CalMuiR(1, -10, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(true_values, LED::CalProbability::CalMuiR(1, -10, ldm), LED::CalProbability::GLB_MU3R); // T2K
    LED::CalProbability::SetModesCutoff(50);

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
    outputFiles.InitOutput(MYFILE1_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE2_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE3_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE4_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE5_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 2);
    true_rates_N1 = glbGetSignalRatePtr(2, 2); // 1Re1de event
    true_rates_N2 = glbGetBGRatePtr(2, 2);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_R);
    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_C1R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C2R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C3R);
    glbSetOscParams(true_values, 0.1, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(true_values, LED::CalProbability::CalMuiR(10, -10, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(true_values, LED::CalProbability::CalMuiR(10, -10, ldm), LED::CalProbability::GLB_MU3R); // T2K
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE1_cirp5, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE2_cirp5, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE3_cirp5, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE4_cirp5, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE5_cirp5, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 2);
    true_rates_N1 = glbGetSignalRatePtr(2, 2); // 1Re1de event
    true_rates_N2 = glbGetBGRatePtr(2, 2);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_R);
    glbSetOscParams(true_values, 20, LED::CalProbability::GLB_C1R);
    glbSetOscParams(true_values, -20, LED::CalProbability::GLB_C2R);
    glbSetOscParams(true_values, -20, LED::CalProbability::GLB_C3R);
    glbSetOscParams(true_values, 0.1, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(true_values, LED::CalProbability::CalMuiR(10, -20, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(true_values, LED::CalProbability::CalMuiR(10, -20, ldm), LED::CalProbability::GLB_MU3R); // T2K
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE1_cirn5, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE2_cirn5, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE3_cirn5, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE4_cirn5, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE5_cirn5, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 2);
    true_rates_N1 = glbGetSignalRatePtr(2, 2); // 1Re1de event
    true_rates_N2 = glbGetBGRatePtr(2, 2);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_R);
    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_C1R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C2R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C3R);
    glbSetOscParams(true_values, 0.5, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(true_values, LED::CalProbability::CalMuiR(10, -10, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(true_values, LED::CalProbability::CalMuiR(10, -10, ldm), LED::CalProbability::GLB_MU3R); // T2K
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE1_muirp5, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE2_muirp5, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE3_muirp5, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE4_muirp5, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE5_muirp5, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 2);
    true_rates_N1 = glbGetSignalRatePtr(2, 2); // 1Re1de event
    true_rates_N2 = glbGetBGRatePtr(2, 2);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    // glbSetOscParams(true_values, -0.005 + 0.01, GLB_MU1R);
    // glbSetOscParams(true_values, -0.005 + 0.0275, GLB_MU2R);
    // glbSetOscParams(true_values, -0.005 + 0.1603, GLB_MU3R);
    // glbSetOscillationParameters(true_values);
    // glbSetRates();

    // /* Compute event rate spectra for the neutrino muon-like sample in Super-Kamiokande */
    // InitOutput(MYFILE1_muirn5, "");
    // true_rates_N0 = glbGetRuleRatePtr(0, 0);
    // true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    // true_rates_N2 = glbGetBGRatePtr(0, 0);
    // for (int i = 0; i <= n_bins1 - 1; i++) {
    //     AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    // }

    // /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    // InitOutput(MYFILE2_muirn5, "");
    // true_rates_N0 = glbGetRuleRatePtr(1, 0);
    // true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    // true_rates_N2 = glbGetBGRatePtr(1, 0);
    // for (int i = 0; i <= n_bins2 - 1; i++) {
    //     AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    // }

    // /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    // InitOutput(MYFILE3_muirn5, "");
    // true_rates_N0 = glbGetRuleRatePtr(2, 0);
    // true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    // true_rates_N2 = glbGetBGRatePtr(2, 0);
    // for (int i = 0; i <= n_bins3 - 1; i++) {
    //     AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    // }

    // InitOutput(MYFILE4_muirn5, "");
    // true_rates_N0 = glbGetRuleRatePtr(2, 1);
    // true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    // true_rates_N2 = glbGetBGRatePtr(2, 1);
    // for (int i = 0; i <= n_bins3 - 1; i++) {
    //     AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    // }

    // InitOutput(MYFILE5_muirn5, "");
    // true_rates_N0 = glbGetRuleRatePtr(2, 2);
    // true_rates_N1 = glbGetSignalRatePtr(2, 2); // 1Re1de event
    // true_rates_N2 = glbGetBGRatePtr(2, 2);
    // for (int i = 0; i <= n_bins3 - 1; i++) {
    //     AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    // }

    glbFreeParams(true_values);

    /* Clear experiment list */
    glbClearExperimentList();

    exit(0);
}
