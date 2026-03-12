#include <glb_smear.h>
#include <glb_types.h>
#include <globes/globes.h> // GLoBES library
#include <gsl/gsl_sf_erf.h>
#include <math.h>

// double ChiNOvAZero(int exp, int rule, int n_params, double* x, double* errors, void* user_data) {

// double chi2 = 0.0;

// return chi2;
// }

// double my_sigma(int j, double par1, double par2, const glb_smear* data) {
//     double energy = data->simbincenter[j];

// return par1 * sqrt(energy) + par2 * pow(energy, 2);
// }

// static double my_SmearMatrixElementA(int i, int j, double e_shift, double par1, double par2, const glb_smear* data) {
//     double p1, p2;
//     p1 = 0.5 * (1 + gsl_sf_erf(
//                         (glb_upper_bin_boundary(i, data) - glb_sbin_center(j, data) - e_shift) /
//                         (sqrt(2) * my_sigma(j, par1, par2, data))));
//     p2 = 0.5 * (1 + gsl_sf_erf(
//                         (glb_lower_bin_boundary(i, data) - glb_sbin_center(j, data) - e_shift) /
//                         (sqrt(2) * my_sigma(j, par1, par2, data))));

// return (p1 - p2);
// }

// void create_my_smear(int exp, int channel, double e_shift, double par1, double par2) // static double**
// {
//     int i, j;
//     double aux;

// struct glb_experiment* e = glb_experiment_list[exp];

// int rows = e->smear_data[channel]->numofbins + 1;
// int columns = e->smear_data[channel]->simbins + 1;

// for (int i = 0; i < e->smear_data[channel]->numofbins; i++) {
//     free(e->smear[channel][i]);
// }

// free(e->lowrange[channel]);
// free(e->uprange[channel]);

// e->lowrange[channel] = (int*)malloc(rows * sizeof(int));
// e->uprange[channel] = (int*)malloc(rows * sizeof(int));
// for (int i = 0; i < e->smear_data[channel]->numofbins; i++) {
//     e->smear[channel][i] = (double*)malloc(columns * sizeof(double));
// }

// for (i = 0; i < e->smear_data[channel]->numofbins; i++) {
//     e->lowrange[channel][i] = 0;
//     e->uprange[channel][i] = e->smear_data[channel]->simbins - 1;

// for (j = 0; j < e->smear_data[channel]->simbins; j++) {
//     aux = my_SmearMatrixElementA(i, j, e_shift, par1, par2, e->smear_data[channel]);
//     if (aux > 1.0e-6)
//         e->smear[channel][i][j] = aux;
//     else
//         e->smear[channel][i][j] = 0;
// }
// }

// e->lowrange[channel][i] = -1;
// e->uprange[channel][i] = -1;
// e->smear[channel][i] = NULL;
// }

inline double log_likelihood_NOVA(double true_rate, double fit_rate) {
    if (true_rate > 1.0e-3)
        return 2 * (fit_rate - true_rate * (1 + log(fit_rate / true_rate)));
    else
        return 0;
}

double ChiNOvAPeripheral(int exp, int rule, int n_params, double* x, double* errors, void* user_data) {

    double chi2 = 0.0;
    int ew_low, ew_high;

    double *signal_fit_rates, *true_rates;
    double* bg_fit_rates;
    double total_rates = 0.0;
    double BackTotal_rates = 0.0;
    double Data = 0.0;
    double cosmic = 0.0;

    signal_fit_rates = glbGetSignalFitRatePtr(exp, rule);
    bg_fit_rates = glbGetBGFitRatePtr(exp, rule);
    true_rates = glbGetRuleRatePtr(exp, rule);

    glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);

    if (rule == 4) {

        for (int i = ew_low; i <= ew_high; i++) {
            total_rates = total_rates + signal_fit_rates[i];     // 8.1
            BackTotal_rates = BackTotal_rates + bg_fit_rates[i]; // 4.7
            Data = 27.0;
            cosmic = 3.835616438;
        }
    }
    if (rule == 5) {

        for (int i = ew_low; i <= ew_high; i++) {
            total_rates = total_rates + signal_fit_rates[i];     // 2.819
            BackTotal_rates = BackTotal_rates + bg_fit_rates[i]; // 2.22
            Data = 3.0;
            cosmic = 0.909090909;
        }
    }

    chi2 = log_likelihood_NOVA(Data, (1 + x[0]) * total_rates + (1 + x[1]) * (BackTotal_rates) + cosmic);

    // Systematic part of chi^2
    for (int i = 0; i < n_params; i++)
        chi2 += (x[i] / errors[i]) * (x[i] / errors[i]);

    return chi2;
}

void InitializeNOvA(glb_exp* in, int* counter) {
    int Nexp, Nsmear;

    glbInitExperiment("../experiments/NOvA_2025.04361/NOvA_FD_numu_samples.glb", in, counter);

    // Nexp   = *counter - 1;
    // Nsmear = 0;

    // create_my_smear(Nexp, Nsmear, -0.1, 0.091, 0.00);

    // Nsmear = Nsmear + 1;
    // create_my_smear(Nexp, Nsmear, 0.00, 0.082, 0.00);

    // Nsmear = Nsmear + 1;
    // create_my_smear(Nexp, Nsmear, -0.03, 0.15, 0.00);

    // Nsmear = Nsmear + 1;
    // create_my_smear(Nexp, Nsmear, 0.00, 0.082, 0.00);

    glbDefineChiFunction(&ChiNOvAPeripheral, 2, "chiPeripheral", NULL);
    glbInitExperiment("../experiments/NOvA_2025.04361/NOvA_FD_nue_samples_Peripheral.glb", in, counter);
    glbInitExperiment("../experiments/NOvA_2025.04361/NOvA_FD_nue_lowEnergy.glb", in, counter);

    // Nexp   = *counter - 1;
    // Nsmear = 0;

    // create_my_smear(Nexp, Nsmear, -0.15, 0.107, 0.0);

    // Nsmear = Nsmear + 1;
    // create_my_smear(Nexp, Nsmear, -0.04, 0.02, 0.01);

    // Nsmear = Nsmear + 1;
    // create_my_smear(Nexp, Nsmear, -0.2, 0.107, 0.0);

    // Nsmear = Nsmear + 1;
    // create_my_smear(Nexp, Nsmear, -0.06, 0.017, 0.01);

    // Nsmear = Nsmear + 1;
    // create_my_smear(Nexp, Nsmear, -0.8, 0.1, 0.0);

    // Nsmear = Nsmear + 1;
    // create_my_smear(Nexp, Nsmear, -0.5, 0.01, 0.01);

    // Nsmear = Nsmear + 1;
    // create_my_smear(Nexp, Nsmear, -0.8, 0.1, 0.0);

    // Nsmear = Nsmear + 1;
    // create_my_smear(Nexp, Nsmear, -0.5, 0.01, 0.01);
}
