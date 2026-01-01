#include <glb_smear.h>
#include <glb_types.h>
#include <globes/globes.h> // GLoBES library
#include <gsl/gsl_sf_erf.h>
#include <math.h>

double log_likelihood_NOVA(double true_rate, double fit_rate) {
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
    double cosmic = 0;
    double Data = 0.0;

    signal_fit_rates = glbGetSignalFitRatePtr(exp, rule);
    bg_fit_rates = glbGetBGFitRatePtr(exp, rule);
    true_rates = glbGetRuleRatePtr(exp, rule);

    glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);

    if (rule == 4) {

        for (int i = ew_low; i <= ew_high; i++) {
            total_rates = total_rates + signal_fit_rates[i];     // 8.1
            BackTotal_rates = BackTotal_rates + bg_fit_rates[i]; // 4.7
            cosmic = 1.677777777;
            Data = 16;
        }
    }
    if (rule == 5) {

        for (int i = ew_low; i <= ew_high; i++) {
            total_rates = total_rates + signal_fit_rates[i];     // 2.819
            BackTotal_rates = BackTotal_rates + bg_fit_rates[i]; // 2.22
            Data = 2;
            cosmic = 0.993836672;
        }
    }

    chi2 = log_likelihood_NOVA(Data, (1 + x[0]) * total_rates + (1 + x[1]) * (BackTotal_rates) + cosmic);

    // Systematic part of chi^2
    for (int i = 0; i < n_params; i++) {
        chi2 += (x[i] / errors[i]) * (x[i] / errors[i]);
    }
    return chi2;
}

void InitializeNOvA(glb_exp* in, int* counter) {

    glbInitExperiment("../experiments/NOvA_2108.08219/NOvA_FD_numu_samples.glb", in, counter);

    glbDefineChiFunction(&ChiNOvAPeripheral, 2, "chiPeripheral", NULL);
    glbInitExperiment("../experiments/NOvA_2108.08219/NOvA_FD_nue_samples_Peripheral.glb", in, counter);
}
