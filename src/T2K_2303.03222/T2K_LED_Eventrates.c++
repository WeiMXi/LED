extern "C" {
#include "T2K_Setup.h"
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

const std::string MYFILE1_cir4 = "../data/T2K/LED/T2K_numu_FHC_cir4.dat";
const std::string MYFILE2_cir4 = "../data/T2K/LED/T2K_numu_RHC_cir4.dat";
const std::string MYFILE3_cir4 = "../data/T2K/LED/T2K_nue_FHC_cir4.dat";
const std::string MYFILE4_cir4 = "../data/T2K/LED/T2K_nue_RHC_cir4.dat";
const std::string MYFILE5_cir4 = "../data/T2K/LED/T2K_nueCCpi_cir4.dat";

const std::string MYFILE1_cir8 = "../data/T2K/LED/T2K_numu_FHC_cir8.dat";
const std::string MYFILE2_cir8 = "../data/T2K/LED/T2K_numu_RHC_cir8.dat";
const std::string MYFILE3_cir8 = "../data/T2K/LED/T2K_nue_FHC_cir8.dat";
const std::string MYFILE4_cir8 = "../data/T2K/LED/T2K_nue_RHC_cir8.dat";
const std::string MYFILE5_cir8 = "../data/T2K/LED/T2K_nueCCpi_cir8.dat";

const std::string MYFILE1_muir1 = "../data/T2K/LED/T2K_numu_FHC_muir1.dat";
const std::string MYFILE2_muir1 = "../data/T2K/LED/T2K_numu_RHC_muir1.dat";
const std::string MYFILE3_muir1 = "../data/T2K/LED/T2K_nue_FHC_muir1.dat";
const std::string MYFILE4_muir1 = "../data/T2K/LED/T2K_nue_RHC_muir1.dat";
const std::string MYFILE5_muir1 = "../data/T2K/LED/T2K_nueCCpi_muir1.dat";

const std::string MYFILE1_muir15 = "../data/T2K/LED/T2K_numu_FHC_muir1.5.dat";
const std::string MYFILE2_muir15 = "../data/T2K/LED/T2K_numu_RHC_muir1.5.dat";
const std::string MYFILE3_muir15 = "../data/T2K/LED/T2K_nue_FHC_muir1.5.dat";
const std::string MYFILE4_muir15 = "../data/T2K/LED/T2K_nue_RHC_muir1.5.dat";
const std::string MYFILE5_muir15 = "../data/T2K/LED/T2K_nueCCpi_muir1.5.dat";

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
    double theta12 = asin(sqrt(0.307));
    double theta13 = asin(sqrt(0.02195));
    double theta23 = asin(sqrt(0.561));
    double deltacp = -1.97;
    double sdm = 7.49e-5;
    double ldm = 2.495e-3 + sdm;

    glbSetOscParams(true_values, 1, LED::CalProbability::GLB_R);
    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_C1R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C2R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C3R);
    glbSetOscParams(true_values, 0.1, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(true_values, LED::CalProbability::CalculateMuiR(1, -10, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(true_values, LED::CalProbability::CalculateMuiR(1, -10, ldm), LED::CalProbability::GLB_MU3R); // T2K
    LED::CalProbability::SetModesCutoff(40);

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
    glbSetOscParams(true_values, 4, LED::CalProbability::GLB_C1R);
    glbSetOscParams(true_values, -4, LED::CalProbability::GLB_C2R);
    glbSetOscParams(true_values, -4, LED::CalProbability::GLB_C3R);
    glbSetOscParams(true_values, 0.1, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(true_values, LED::CalProbability::CalculateMuiR(10, -4, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(true_values, LED::CalProbability::CalculateMuiR(10, -4, ldm), LED::CalProbability::GLB_MU3R); // T2K
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE1_cir4, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE2_cir4, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE3_cir4, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE4_cir4, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE5_cir4, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 2);
    true_rates_N1 = glbGetSignalRatePtr(2, 2); // 1Re1de event
    true_rates_N2 = glbGetBGRatePtr(2, 2);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_R);
    glbSetOscParams(true_values, 8, LED::CalProbability::GLB_C1R);
    glbSetOscParams(true_values, -8, LED::CalProbability::GLB_C2R);
    glbSetOscParams(true_values, -8, LED::CalProbability::GLB_C3R);
    glbSetOscParams(true_values, 0.1, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(true_values, LED::CalProbability::CalculateMuiR(10, -8, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(true_values, LED::CalProbability::CalculateMuiR(10, -8, ldm), LED::CalProbability::GLB_MU3R); // T2K
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE1_cir8, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE2_cir8, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE3_cir8, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE4_cir8, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE5_cir8, "");
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
    glbSetOscParams(true_values, 1, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(true_values, LED::CalProbability::CalculateMuiR(10, -10, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(true_values, LED::CalProbability::CalculateMuiR(10, -10, ldm), LED::CalProbability::GLB_MU3R); // T2K
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE1_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE2_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE3_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE4_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE5_muir1, "");
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
    glbSetOscParams(true_values, 1.5, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(true_values, LED::CalProbability::CalculateMuiR(10, -10, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(true_values, LED::CalProbability::CalculateMuiR(10, -10, ldm), LED::CalProbability::GLB_MU3R); // T2K
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE1_muir15, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0); // Neutrino beam muon-like events rates
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numu_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE2_muir15, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0); // Neutrino beam antimuon-like events rates
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(numubar_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like sample in Super-Kamiokande */
    outputFiles.InitOutput(MYFILE3_muir15, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 0);
    true_rates_N1 = glbGetSignalRatePtr(2, 0); // 1Re FHC event
    true_rates_N2 = glbGetBGRatePtr(2, 0);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE4_muir15, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 1);
    true_rates_N1 = glbGetSignalRatePtr(2, 1); // 1Re RHC event
    true_rates_N2 = glbGetBGRatePtr(2, 1);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    outputFiles.InitOutput(MYFILE5_muir15, "");
    true_rates_N0 = glbGetRuleRatePtr(2, 2);
    true_rates_N1 = glbGetSignalRatePtr(2, 2); // 1Re1de event
    true_rates_N2 = glbGetBGRatePtr(2, 2);
    for (int i = 0; i <= n_bins3 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbFreeParams(true_values);

    /* Clear experiment list */
    glbClearExperimentList();

    exit(0);
}
