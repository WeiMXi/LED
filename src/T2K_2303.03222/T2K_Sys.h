#ifndef T2K_SYS_H
#define T2K_SYS_H

#include <globes/globes.h> // GLoBES library
#include <math.h>
double chi_zero(int exp, int rule, int n_params, double* x, double* errors, void* user_data) {

    return 0.0;
}

/***************************************************************************
 * Function glb_likelihood                                                 *
 ***************************************************************************
 * Calculate the Poisson likelihood for a specific true and fitted event   *
 * rate.                                                                   *
 ***************************************************************************/
inline double log_fatorial(double in) {
    int i;
    double factorial = 0;
    for (i = 1; i <= in; i++) factorial += log(i);
    return factorial;
}

inline double glb_likelihood(double true_rate, double fit_rate) {
    double res;
    res = fit_rate - true_rate;
    if (true_rate > 0) {
        if (fit_rate <= 0.0)
            res = 1e100;
        else
            res += true_rate * log(true_rate / fit_rate);
    } else
        res = fabs(res);

    return 2.0 * res;
}
// There are some '-1' in data set, which means no count in this bin

/***************************************************************************
 * Function glb_prior                                                      *
 ***************************************************************************
 * Calculate prior term of the form ((x - x_center)/error)^2.              *
 ***************************************************************************/
inline double glb_prior(double x, double center, double sigma) {
    double tmp = (x - center) / sigma;
    return tmp * tmp;
}

/***************************************************************************
 *                     C H I ^ 2   F U N C T I O N S                       *
 ***************************************************************************/

/***************************************************************************
 * Function glbChiSpectrumTilt                                             *
 ***************************************************************************
 * chi^2 including the standard signal and background errors, as well as   *
 * spectral information                                                    *
 ***************************************************************************/
double ChiT2K(int exp, int rule, int n_params, double* x, double* errors, void* user_data) {
    double *signal_fit_rates, *true_rates;
    double* bg_fit_rates;
    double* bin_centers;
    double signal_norm, signal_tilt;
    double bg_norm_center, bg_tilt_center;
    double bg_norm, bg_tilt;
    int ew_low, ew_high;
    double emin, emax, ecenter;
    double fit_rate;
    double chi2 = 0.0;
    int i;

    double Devents1, Devents2;

    signal_fit_rates = glbGetSignalFitRatePtr(exp, rule);
    bg_fit_rates = glbGetBGFitRatePtr(exp, rule);
    true_rates = glbGetRuleRatePtr(exp, rule);
    bin_centers = glbGetBinCentersListPtr(exp);

    glbGetEminEmax(exp, &emin, &emax);
    ecenter = 0.5 * (emax + emin);
    glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
    glbGetBGCenters(exp, rule, &bg_norm_center, &bg_tilt_center);
    signal_norm = 1.0 + x[0];
    signal_tilt = x[1] / (emax - emin);
    bg_norm = bg_norm_center * (1.0 + x[2]);
    bg_tilt = x[3] / (emax - emin);
    for (i = ew_low; i <= ew_high; i++) {
        fit_rate = signal_norm * signal_fit_rates[i] + signal_tilt * (bin_centers[i] - ecenter) * signal_fit_rates[i] + bg_norm * bg_fit_rates[i] + bg_tilt * (bin_centers[i] - ecenter) * bg_fit_rates[i];
        chi2 += glb_likelihood(true_rates[i], fit_rate);
    }

    chi2 += glb_prior(x[0], 0.0, errors[0]) + glb_prior(x[1], 0.0, errors[1]) + glb_prior(bg_norm, bg_norm_center, errors[2]) + glb_prior(bg_tilt, bg_tilt_center, errors[3]);

    return chi2;
}
#endif
