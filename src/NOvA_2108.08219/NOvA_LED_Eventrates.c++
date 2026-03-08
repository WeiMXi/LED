extern "C" {
#include "NOvA_Setup_New.h"
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

const std::string MYFILE1_4Dfinite = "../data/NOvA/LED/NOvA_numu_4Dfinite.dat";
const std::string MYFILE2_4Dfinite = "../data/NOvA/LED/NOvA_anumu_4Dfinite.dat";
const std::string MYFILE3_4Dfinite = "../data/NOvA/LED/NOvA_nue_LCNN_4Dfinite.dat";
const std::string MYFILE4_4Dfinite = "../data/NOvA/LED/NOvA_nue_HCNN_4Dfinite.dat";
const std::string MYFILE5_4Dfinite = "../data/NOvA/LED/NOvA_anue_LCNN_4Dfinite.dat";
const std::string MYFILE6_4Dfinite = "../data/NOvA/LED/NOvA_anue_HCNN_4Dfinite.dat";

const std::string MYFILE1_cir4 = "../data/NOvA/LED/NOvA_numu_cir4.dat";
const std::string MYFILE2_cir4 = "../data/NOvA/LED/NOvA_anumu_cir4.dat";
const std::string MYFILE3_cir4 = "../data/NOvA/LED/NOvA_nue_LCNN_cir4.dat";
const std::string MYFILE4_cir4 = "../data/NOvA/LED/NOvA_nue_HCNN_cir4.dat";
const std::string MYFILE5_cir4 = "../data/NOvA/LED/NOvA_anue_LCNN_cir4.dat";
const std::string MYFILE6_cir4 = "../data/NOvA/LED/NOvA_anue_HCNN_cir4.dat";

const std::string MYFILE1_cir8 = "../data/NOvA/LED/NOvA_numu_cir8.dat";
const std::string MYFILE2_cir8 = "../data/NOvA/LED/NOvA_anumu_cir8.dat";
const std::string MYFILE3_cir8 = "../data/NOvA/LED/NOvA_nue_LCNN_cir8.dat";
const std::string MYFILE4_cir8 = "../data/NOvA/LED/NOvA_nue_HCNN_cir8.dat";
const std::string MYFILE5_cir8 = "../data/NOvA/LED/NOvA_anue_LCNN_cir8.dat";
const std::string MYFILE6_cir8 = "../data/NOvA/LED/NOvA_anue_HCNN_cir8.dat";

const std::string MYFILE1_muir1 = "../data/NOvA/LED/NOvA_numu_muir1.dat";
const std::string MYFILE2_muir1 = "../data/NOvA/LED/NOvA_anumu_muir1.dat";
const std::string MYFILE3_muir1 = "../data/NOvA/LED/NOvA_nue_LCNN_muir1.dat";
const std::string MYFILE4_muir1 = "../data/NOvA/LED/NOvA_nue_HCNN_muir1.dat";
const std::string MYFILE5_muir1 = "../data/NOvA/LED/NOvA_anue_LCNN_muir1.dat";
const std::string MYFILE6_muir1 = "../data/NOvA/LED/NOvA_anue_HCNN_muir1.dat";

const std::string MYFILE1_muir15 = "../data/NOvA/LED/NOvA_numu_muir1.5.dat";
const std::string MYFILE2_muir15 = "../data/NOvA/LED/NOvA_anumu_muir1.5.dat";
const std::string MYFILE3_muir15 = "../data/NOvA/LED/NOvA_nue_LCNN_muir1.5.dat";
const std::string MYFILE4_muir15 = "../data/NOvA/LED/NOvA_nue_HCNN_muir1.5.dat";
const std::string MYFILE5_muir15 = "../data/NOvA/LED/NOvA_anue_LCNN_muir1.5.dat";
const std::string MYFILE6_muir15 = "../data/NOvA/LED/NOvA_anue_HCNN_muir1.5.dat";
LED::IO::Output outputFiles;

int main(int argc, char* argv[]) {
    /* Initialize libglobes */
    glbInit(argv[0]);
    /* Initialize NOvA */
    InitializeNOvA(&glb_experiment_list[0], &glb_num_of_exps);
    glbRegisterProbabilityEngine(11, /*Number of parameters*/
                                 &LED::CalProbability::my_probability_matrix,
                                 &LED::CalProbability::my_set_oscillation_parameters,
                                 &LED::CalProbability::my_get_oscillation_parameters,
                                 NULL);

    /* Define standard oscillation parameters for NO in NOvA */
    double theta12 = asin(sqrt(0.307));
    double theta13 = asin(sqrt(0.02195));
    double theta23 = asin(sqrt(0.57));
    double deltacp = 0.82 * M_PI;
    double sdm = 7.49e-5;
    double ldm = 2.41e-3 + sdm;

    /* Initialize the parameter vector */
    glb_params true_values = glbAllocParams();

    glbSetOscParams(true_values, 1, LED::CalProbability::GLB_R);
    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_C1R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C2R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C3R);
    double m02 = LED::CalProbability::CalLightestm2(1, 10, 0.1);
    glbSetOscParams(true_values, m02, LED::CalProbability::GLB_M0SQUARE);

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

    int n_bins1 = glbGetNumberOfBins(0);
    int n_bins2 = glbGetNumberOfBins(1);

    numu_bin_centers_E = glbGetBinCentersListPtr(0);
    numubar_bin_centers_E = glbGetBinCentersListPtr(0);
    nue_bin_centers_E = glbGetBinCentersListPtr(1);

    numu_bin_widths = glbGetBinSizeListPtr(0);
    numubar_bin_widths = glbGetBinSizeListPtr(0);
    nue_bin_widths = glbGetBinSizeListPtr(1);

    /* Compute event rate spectra for the neutrino muon-like event */
    outputFiles.InitOutput(MYFILE1_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0);
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numu_bin_centers_E[i], (true_rates_N1[i] + true_rates_N2[i]) / numu_bin_widths[i] * 0.1, numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like event */
    outputFiles.InitOutput(MYFILE2_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 1);
    true_rates_N1 = glbGetSignalRatePtr(0, 1);
    true_rates_N2 = glbGetBGRatePtr(0, 1);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numubar_bin_centers_E[i], (true_rates_N1[i] + true_rates_N2[i]) / numu_bin_widths[i] * 0.1, numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with low CNN */
    outputFiles.InitOutput(MYFILE3_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0);
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with high CNN */
    outputFiles.InitOutput(MYFILE4_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 1);
    true_rates_N1 = glbGetSignalRatePtr(1, 1);
    true_rates_N2 = glbGetBGRatePtr(1, 1);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino antielectron-like samples with low CNN */
    outputFiles.InitOutput(MYFILE5_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 2);
    true_rates_N1 = glbGetSignalRatePtr(1, 2);
    true_rates_N2 = glbGetBGRatePtr(1, 2);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with high CNN */
    outputFiles.InitOutput(MYFILE6_4Dfinite, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 3);
    true_rates_N1 = glbGetSignalRatePtr(1, 3);
    true_rates_N2 = glbGetBGRatePtr(1, 3);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }
    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_R);
    glbSetOscParams(true_values, 4, LED::CalProbability::GLB_C1R);
    glbSetOscParams(true_values, -4, LED::CalProbability::GLB_C2R);
    glbSetOscParams(true_values, -4, LED::CalProbability::GLB_C3R);
    m02 = LED::CalProbability::CalLightestm2(10, 4, 0.1);
    glbSetOscParams(true_values, m02, LED::CalProbability::GLB_M0SQUARE);
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like event */
    outputFiles.InitOutput(MYFILE1_cir4, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0);
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numu_bin_centers_E[i], (true_rates_N1[i] + true_rates_N2[i]) / numu_bin_widths[i] * 0.1, numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like event */
    outputFiles.InitOutput(MYFILE2_cir4, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 1);
    true_rates_N1 = glbGetSignalRatePtr(0, 1);
    true_rates_N2 = glbGetBGRatePtr(0, 1);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numubar_bin_centers_E[i], (true_rates_N1[i] + true_rates_N2[i]) / numu_bin_widths[i] * 0.1, numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with low CNN */
    outputFiles.InitOutput(MYFILE3_cir4, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0);
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with high CNN */
    outputFiles.InitOutput(MYFILE4_cir4, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 1);
    true_rates_N1 = glbGetSignalRatePtr(1, 1);
    true_rates_N2 = glbGetBGRatePtr(1, 1);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino antielectron-like samples with low CNN */
    outputFiles.InitOutput(MYFILE5_cir4, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 2);
    true_rates_N1 = glbGetSignalRatePtr(1, 2);
    true_rates_N2 = glbGetBGRatePtr(1, 2);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with high CNN */
    outputFiles.InitOutput(MYFILE6_cir4, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 3);
    true_rates_N1 = glbGetSignalRatePtr(1, 3);
    true_rates_N2 = glbGetBGRatePtr(1, 3);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_R);
    glbSetOscParams(true_values, 8, LED::CalProbability::GLB_C1R);
    glbSetOscParams(true_values, -8, LED::CalProbability::GLB_C2R);
    glbSetOscParams(true_values, -8, LED::CalProbability::GLB_C3R);
    m02 = LED::CalProbability::CalLightestm2(10, 8, 0.1);
    glbSetOscParams(true_values, m02, LED::CalProbability::GLB_M0SQUARE);
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like event */
    outputFiles.InitOutput(MYFILE1_cir8, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0);
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numu_bin_centers_E[i], (true_rates_N1[i] + true_rates_N2[i]) / numu_bin_widths[i] * 0.1, numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like event */
    outputFiles.InitOutput(MYFILE2_cir8, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 1);
    true_rates_N1 = glbGetSignalRatePtr(0, 1);
    true_rates_N2 = glbGetBGRatePtr(0, 1);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numubar_bin_centers_E[i], (true_rates_N1[i] + true_rates_N2[i]) / numu_bin_widths[i] * 0.1, numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with low CNN */
    outputFiles.InitOutput(MYFILE3_cir8, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0);
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with high CNN */
    outputFiles.InitOutput(MYFILE4_cir8, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 1);
    true_rates_N1 = glbGetSignalRatePtr(1, 1);
    true_rates_N2 = glbGetBGRatePtr(1, 1);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino antielectron-like samples with low CNN */
    outputFiles.InitOutput(MYFILE5_cir8, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 2);
    true_rates_N1 = glbGetSignalRatePtr(1, 2);
    true_rates_N2 = glbGetBGRatePtr(1, 2);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with high CNN */
    outputFiles.InitOutput(MYFILE6_cir8, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 3);
    true_rates_N1 = glbGetSignalRatePtr(1, 3);
    true_rates_N2 = glbGetBGRatePtr(1, 3);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_R);
    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_C1R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C2R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C3R);
    m02 = LED::CalProbability::CalLightestm2(10, 10, 1);
    glbSetOscParams(true_values, m02, LED::CalProbability::GLB_M0SQUARE);
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like event */
    outputFiles.InitOutput(MYFILE1_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0);
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numu_bin_centers_E[i], (true_rates_N1[i] + true_rates_N2[i]) / numu_bin_widths[i] * 0.1, numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like event */
    outputFiles.InitOutput(MYFILE2_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 1);
    true_rates_N1 = glbGetSignalRatePtr(0, 1);
    true_rates_N2 = glbGetBGRatePtr(0, 1);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numubar_bin_centers_E[i], (true_rates_N1[i] + true_rates_N2[i]) / numu_bin_widths[i] * 0.1, numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with low CNN */
    outputFiles.InitOutput(MYFILE3_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0);
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with high CNN */
    outputFiles.InitOutput(MYFILE4_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 1);
    true_rates_N1 = glbGetSignalRatePtr(1, 1);
    true_rates_N2 = glbGetBGRatePtr(1, 1);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino antielectron-like samples with low CNN */
    outputFiles.InitOutput(MYFILE5_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 2);
    true_rates_N1 = glbGetSignalRatePtr(1, 2);
    true_rates_N2 = glbGetBGRatePtr(1, 2);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with high CNN */
    outputFiles.InitOutput(MYFILE6_muir1, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 3);
    true_rates_N1 = glbGetSignalRatePtr(1, 3);
    true_rates_N2 = glbGetBGRatePtr(1, 3);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_R);
    glbSetOscParams(true_values, 10, LED::CalProbability::GLB_C1R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C2R);
    glbSetOscParams(true_values, -10, LED::CalProbability::GLB_C3R);
    m02 = LED::CalProbability::CalLightestm2(10, 10, 1.5);
    glbSetOscParams(true_values, m02, LED::CalProbability::GLB_M0SQUARE);
    glbSetOscillationParameters(true_values);
    glbSetRates();

    /* Compute event rate spectra for the neutrino muon-like event */
    outputFiles.InitOutput(MYFILE1_muir15, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 0);
    true_rates_N1 = glbGetSignalRatePtr(0, 0);
    true_rates_N2 = glbGetBGRatePtr(0, 0);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numu_bin_centers_E[i], (true_rates_N1[i] + true_rates_N2[i]) / numu_bin_widths[i] * 0.1, numu_bin_widths[i]);
    }

    /* Compute event rate spectra for the antineutrino muon-like event */
    outputFiles.InitOutput(MYFILE2_muir15, "");
    true_rates_N0 = glbGetRuleRatePtr(0, 1);
    true_rates_N1 = glbGetSignalRatePtr(0, 1);
    true_rates_N2 = glbGetBGRatePtr(0, 1);
    for (int i = 0; i <= n_bins1 - 1; i++) {
        outputFiles.AddToOutput(numubar_bin_centers_E[i], (true_rates_N1[i] + true_rates_N2[i]) / numu_bin_widths[i] * 0.1, numubar_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with low CNN */
    outputFiles.InitOutput(MYFILE3_muir15, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 0);
    true_rates_N1 = glbGetSignalRatePtr(1, 0);
    true_rates_N2 = glbGetBGRatePtr(1, 0);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with high CNN */
    outputFiles.InitOutput(MYFILE4_muir15, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 1);
    true_rates_N1 = glbGetSignalRatePtr(1, 1);
    true_rates_N2 = glbGetBGRatePtr(1, 1);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino antielectron-like samples with low CNN */
    outputFiles.InitOutput(MYFILE5_muir15, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 2);
    true_rates_N1 = glbGetSignalRatePtr(1, 2);
    true_rates_N2 = glbGetBGRatePtr(1, 2);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    /* Compute event rate spectra for the neutrino electron-like samples with high CNN */
    outputFiles.InitOutput(MYFILE6_muir15, "");
    true_rates_N0 = glbGetRuleRatePtr(1, 3);
    true_rates_N1 = glbGetSignalRatePtr(1, 3);
    true_rates_N2 = glbGetBGRatePtr(1, 3);
    for (int i = 0; i <= n_bins2 - 1; i++) {
        outputFiles.AddToOutput(nue_bin_centers_E[i], true_rates_N1[i] + true_rates_N2[i], nue_bin_widths[i]);
    }

    glbFreeParams(true_values);

    /* Clear experiment list */
    glbClearExperimentList();

    exit(0);
}
