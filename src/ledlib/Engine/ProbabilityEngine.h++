#pragma once

#ifdef HAVE_CUDA
#    include "CUDA/CUDAHermDiag.h++"
#endif
#include "globes/globes.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <iostream>
#include <numbers>
#include <vector>

namespace LED::CalProbability {

inline double sq(const double x) { return x * x; }

constexpr double TOLERANCE = 1.0e-10;
constexpr double DMSQ21 = 7.49e-5;      // [eV^2] NuFit 6.0
constexpr double DMSQ31_NH = 2.513e-3;  // [eV^2] NuFit 6.0
constexpr double DMSQ32_IH = -2.484e-3; // [eV^2] NuFit 6.0

constexpr double GLB_V_factor = 7.5e-14; // Conversion factor for matter potentials
constexpr double GLB_ne_mantle = 0.5;    // Effective electron numbers for calculation

enum ParIndexes {
    GLB_R = 6,
    GLB_C1R = 7,
    GLB_C2R = 8,
    GLB_C3R = 9,
    GLB_M0SQUARE = 10
};

/***************************************************************************
 *     U S E R - D E F I N E D   P R O B A B I L I T Y   E N G I N E       *
 ***************************************************************************/

/* Fundamental oscillation parameters */
static double th12, th13, th23; // Mixing angles
static double delta;            // Dirac CP phase
static double mq[3];            // Squared masses static
static double sigma_E;
static double deltacp;
static double sdm;
static double ldm;
static double R, c1R, c2R, c3R, mu1R, mu2R, mu3R; // 5d bulk
static double m0Sqare;                            // the sqaure of lightest nu mass
static int modes;                                 // modes cutoff in 5d bulk

double mum_to_eVinv(const double x) { return x * 5.06773; }
double eVinv_to_mum(const double x) { return x / 5.06773; }

inline double CalMuiR(const double r, const double ciR, const double mi2) {
    double r_eVinv = mum_to_eVinv(r);
    if (std::abs(ciR) > TOLERANCE) {
        return std::sqrt(mi2 * (r_eVinv) * (r_eVinv) / 2 / std::numbers::pi / std::abs(ciR));
    } else {
        return sqrt(mi2 * (r_eVinv) * (r_eVinv));
    }
}

inline double CalculateMuiR(const double r, const double ciR, const double mi2) {
    double y{}, muiR2{};
    double r_eVinv = mum_to_eVinv(r);
    if (ciR * ciR - r_eVinv * r_eVinv * mi2 >= 0) {
        y = std::sqrt(ciR * ciR - r_eVinv * r_eVinv * mi2);
        muiR2 = (-y * y + ciR * ciR) / (std::numbers::pi * y / std::tanh(std::numbers::pi * y) - std::numbers::pi * ciR);
    } else {
        y = std::sqrt(r_eVinv * r_eVinv * mi2 - ciR * ciR);
        muiR2 = (y * y + ciR * ciR) / (std::numbers::pi * y / std::tan(std::numbers::pi * y) - std::numbers::pi * ciR);
    }
    return std::sqrt(muiR2);
}

// KK masses equation

bool is_coth(double cR, double muR) {
    return sq(cR) + M_PI * cR * sq(muR) - sq(muR) >= 0;
}

double masseq_vac(double x, void* params) {
    double* p = (double*)params;
    double cR = p[0];
    double muR = p[1];
    if (x > sq(cR) + TOLERANCE)
        return M_PI * sq(muR) * sqrt(x - sq(cR)) / tan(M_PI * sqrt(x - sq(cR))) - x - M_PI * cR * sq(muR);
    else if (x < sq(cR) - TOLERANCE)
        return M_PI * sq(muR) * sqrt(sq(cR) - x) / tanh(M_PI * sqrt(sq(cR) - x)) - x - M_PI * cR * sq(muR);
    else
        return sq(muR) - M_PI * sq(muR) * cR - sq(cR);
}

double cot_y(double y, void* params) {
    double* p = (double*)params;
    double cR = p[0];
    return masseq_vac(sq(cR) + sq(y), params);
}

double solver(double (*eqn)(double, void*), double xl, double xr, void* p) {
    const gsl_root_fsolver_type* T = gsl_root_fsolver_brent; // Brent's method adopted
    gsl_root_fsolver* s = gsl_root_fsolver_alloc(T);
    gsl_function eq;
    eq.function = eqn;
    eq.params = p;

    int status;
    int iter = 0;
    int max_iter = 100;
    double root = -1.0;
    double epsilon = 1.0e-12;
    double eq_xl = eqn(xl, p);
    double eq_xr = eqn(xr, p);

    if (fabs(eq_xl) <= epsilon)
        root = xl;
    else if (fabs(eq_xr) <= epsilon)
        root = xr;
    else {
        gsl_root_fsolver_set(s, &eq, xl, xr);
        do {
            iter = iter + 1;
            status = gsl_root_fsolver_iterate(s);
            root = gsl_root_fsolver_root(s);
            xl = gsl_root_fsolver_x_lower(s);
            xr = gsl_root_fsolver_x_upper(s);
            status = gsl_root_test_interval(xl, xr, 0, epsilon);
            // if (status == GSL_SUCCESS)
            // printf("converged!\n");
            // printf("%5d [%2.10f, %2.10f], %2.10f, %2.10f\n", iter, xl, xr, root, xr-xl);
        } while (status == GSL_CONTINUE && iter < max_iter);
    }
    return root;
}

double find_bracket_bound(double (*eq)(double, void*), double xs, double delta, double bound, void* p) {
    double eqs = eq(xs, p);
    double x = xs + delta;
    if (delta > 0) {
        while (eqs * eq(x, p) > 0) {
            if (x > bound)
                break;
            else
                x += delta;
        }
    } else {
        while (eqs * eq(x, p) > 0) {
            if (x < bound)
                break;
            else
                x += delta;
        }
    }
    // printf("find_bracket_bound: bracket found\n");
    return x;
}

double solve_masseq_vac_coth(double* p) {
    double xl = 0.0;
    double cR = p[0];
    if (fabs(masseq_vac(xl, p)) < TOLERANCE)
        return xl;
    else {
        double delta = sq(cR) / 100;
        double xr = find_bracket_bound(masseq_vac, xl, delta, sq(cR) - TOLERANCE, p);
        return solver(masseq_vac, xl, xr, p);
    }
};

double solve_masseq_vac_cot(int n, double* p) {
    double cR = p[0];
    double yl = n + TOLERANCE;
    if (fabs(cot_y(yl, p)) < TOLERANCE)
        return sq(cR) + sq(yl);
    else {
        double delta = 1.0e-2;
        double yr = find_bracket_bound(cot_y, yl, delta, n + 1 - TOLERANCE, p);
        return sq(cR) + sq(solver(cot_y, yl, yr, p));
    }
};

double solve_masseq_vac(int n, double* p) {
    double cR = p[0];
    double muR = p[1];
    if (n == 0) {
        if (is_coth(cR, muR)) {
            return solve_masseq_vac_coth(p);
        } else {
            return solve_masseq_vac_cot(0, p);
        }
    } else
        return solve_masseq_vac_cot(n, p);
};

inline int my_set_oscillation_parameters(glb_params p, void* user_data) {
    th12 = glbGetOscParams(p, GLB_THETA_12);
    th13 = glbGetOscParams(p, GLB_THETA_13);
    th23 = glbGetOscParams(p, GLB_THETA_23);
    deltacp = glbGetOscParams(p, GLB_DELTA_CP);
    sdm = glbGetOscParams(p, GLB_DM_21);
    ldm = glbGetOscParams(p, GLB_DM_31);

    /* 5d bulk PARAMETERS*/
    R = glbGetOscParams(p, GLB_R);
    c1R = glbGetOscParams(p, GLB_C1R);
    c2R = glbGetOscParams(p, GLB_C2R);
    c3R = glbGetOscParams(p, GLB_C3R);
    m0Sqare = glbGetOscParams(p, GLB_M0SQUARE);

    return 0;
}

/***************************************************************************
 * Write oscillation parameters from internal data structures into p.      *
 ***************************************************************************/
inline int my_get_oscillation_parameters(glb_params p, void* user_data) {
    glbSetOscParams(p, th12, GLB_THETA_12);
    glbSetOscParams(p, th13, GLB_THETA_13);
    glbSetOscParams(p, th23, GLB_THETA_23);
    glbSetOscParams(p, deltacp, GLB_DELTA_CP);
    glbSetOscParams(p, sdm, GLB_DM_21);
    glbSetOscParams(p, ldm, GLB_DM_31);

    /* 5d bulk PARAMETERS*/
    glbSetOscParams(p, R, GLB_R);
    glbSetOscParams(p, c1R, GLB_C1R);
    glbSetOscParams(p, c2R, GLB_C2R);
    glbSetOscParams(p, c3R, GLB_C3R);
    glbSetOscParams(p, m0Sqare, GLB_M0SQUARE);

    return 0;
}

inline int GetModesCutoff() {
    return modes;
}

inline void SetModesCutoff(const int n) {
    modes = n;
}

inline double CalLightestm2(const double r, const double c1R, const double mu1R) {
    double r_eVinv = mum_to_eVinv(r);
    return std::exp(-2 * std::numbers::pi * std::abs(c1R)) * (2 * std::numbers::pi * std::abs(c1R)) * mu1R * mu1R / r_eVinv / r_eVinv;
}

inline double CalMu1R(const double r, const double c1R, const double lightestm2) {
    double abs_c1R = std::abs(c1R);
    double r_eVinv = mum_to_eVinv(r);
    double exp_term = std::exp(2 * std::numbers::pi * abs_c1R); // 1/e^{-...}=e^...
    double linear_term = 2 * std::numbers::pi * abs_c1R;
    return sqrt(exp_term * lightestm2 * r_eVinv * r_eVinv / linear_term);
}

inline void make_PMNS(gsl_matrix_complex* U, const int cp_sign) {
    if (cp_sign > 0) {
        gsl_matrix_complex_set(U, 0, 0, gsl_complex_rect(cos(th12) * cos(th13), 0));
        gsl_matrix_complex_set(U, 0, 1, gsl_complex_rect(sin(th12) * cos(th13), 0));
        gsl_matrix_complex_set(U, 0, 2, gsl_complex_polar(sin(th13), -deltacp));

        gsl_matrix_complex_set(U, 1, 0, gsl_complex_add(gsl_complex_rect(-sin(th12) * cos(th23), 0), gsl_complex_mul_real(gsl_complex_polar(sin(th13) * sin(th23) * cos(th12), deltacp), -1)));

        gsl_matrix_complex_set(U, 1, 1, gsl_complex_add(gsl_complex_rect(cos(th12) * cos(th23), 0), gsl_complex_mul_real(gsl_complex_polar(-sin(th13) * sin(th23) * sin(th12), deltacp), 1)));

        gsl_matrix_complex_set(U, 1, 2, gsl_complex_rect(sin(th23) * cos(th13), 0));

        gsl_matrix_complex_set(U, 2, 0, gsl_complex_add(gsl_complex_rect(sin(th12) * sin(th23), 0), gsl_complex_mul_real(gsl_complex_polar(-cos(th12) * cos(th23) * sin(th13), deltacp), 1)));

        gsl_matrix_complex_set(U, 2, 1, gsl_complex_add(gsl_complex_rect(-cos(th12) * sin(th23), 0), gsl_complex_mul_real(gsl_complex_polar(-sin(th12) * cos(th23) * sin(th13), deltacp), 1)));

        gsl_matrix_complex_set(U, 2, 2, gsl_complex_rect(cos(th23) * cos(th13), 0));

    } else { /* delta_CP -> -delta_CP */
        gsl_matrix_complex_set(U, 0, 0, gsl_complex_rect(cos(th12) * cos(th13), 0));
        gsl_matrix_complex_set(U, 0, 1, gsl_complex_rect(sin(th12) * cos(th13), 0));
        gsl_matrix_complex_set(U, 0, 2, gsl_complex_polar(sin(th13), deltacp));

        gsl_matrix_complex_set(U, 1, 0, gsl_complex_add(gsl_complex_rect(-sin(th12) * cos(th23), 0), gsl_complex_mul_real(gsl_complex_polar(-cos(th12) * sin(th23) * sin(th13), -deltacp), 1)));

        gsl_matrix_complex_set(U, 1, 1, gsl_complex_add(gsl_complex_rect(cos(th12) * cos(th23), 0), gsl_complex_mul_real(gsl_complex_polar(-sin(th12) * sin(th23) * sin(th13), -deltacp), 1)));

        gsl_matrix_complex_set(U, 1, 2, gsl_complex_rect(sin(th23) * cos(th13), 0));

        gsl_matrix_complex_set(U, 2, 0, gsl_complex_add(gsl_complex_rect(sin(th12) * sin(th23), 0), gsl_complex_mul_real(gsl_complex_polar(-cos(th12) * cos(th23) * sin(th13), -deltacp), 1)));

        gsl_matrix_complex_set(U, 2, 1, gsl_complex_add(gsl_complex_rect(-cos(th12) * sin(th23), 0), gsl_complex_mul_real(gsl_complex_polar(-sin(th12) * cos(th23) * sin(th13), -deltacp), 1)));

        gsl_matrix_complex_set(U, 2, 2, gsl_complex_rect(cos(th23) * cos(th13), 0));
    }
}

inline void make_RRH(gsl_matrix_complex* RRH, gsl_matrix_complex* U,
                     double E, double cp_sign, double density)
// E [GeV], masses [eV], R [micron], density [g/cm^3]
{
    double R_eVinv = mum_to_eVinv(R); // [eV^{-1}]
    double E_eV = E * 1e9;            // fix：should be 1e9 instead of 10e9
    double cR[3] = {c1R, c2R, c3R};
    if (ldm > 0) {
        mu1R = CalMu1R(R, c1R, m0Sqare);
        mu2R = CalculateMuiR(R, c2R, sdm);
        mu3R = CalculateMuiR(R, c3R, ldm);
    } else if (ldm < 0) {
        mu3R = CalMu1R(R, c3R, m0Sqare);
        mu2R = CalculateMuiR(R, c2R, m0Sqare + sdm - ldm);
        mu1R = CalculateMuiR(R, c1R, m0Sqare - ldm);
    }
    double muR[3] = {mu1R, mu2R, mu3R};

    double Vc = cp_sign * density * GLB_V_factor * GLB_ne_mantle;
    double Vn = cp_sign * density * GLB_V_factor * (1.0 - GLB_ne_mantle) / 2.0;

    std::vector<std::vector<double>> RY(3, std::vector<double>(modes));
    for (int i = 0; i < 3; i++) {
        if (std::fabs(cR[i]) < TOLERANCE) {
            RY[i][0] = muR[i];
        } else {
            RY[i][0] = muR[i] * std::sqrt(2 * std::numbers::pi * cR[i] / (std::exp(2 * std::numbers::pi * cR[i]) - 1));
        }

        for (int n = 1; n < modes; n++) {
            RY[i][n] = muR[i] * n * std::sqrt(2.0 / (n * n + cR[i] * cR[i]));
        }
    }

    double A[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; i++) {
        for (int n = 0; n < modes; n++) {
            A[i] += sq(RY[i][n]);
        }
    }

    for (int i = 0; i < 3; i++) {
        gsl_complex Ui = gsl_matrix_complex_get(U, 0, i);
        gsl_matrix_complex_set(RRH, i, i,
                               gsl_complex_rect(A[i] + 2 * E_eV * sq(R_eVinv) * (gsl_complex_abs2(Ui) * Vc - Vn), 0));

        for (int j = 0; j < 3; j++) {
            if (j != i) {
                gsl_complex Uj = gsl_matrix_complex_get(U, 0, j);
                gsl_matrix_complex_set(RRH, i, j,
                                       gsl_complex_mul_real(gsl_complex_mul(gsl_complex_conjugate(Ui), Uj), 2 * E_eV * sq(R_eVinv) * Vc));
            }
        }

        for (int n = 1; n < modes; n++) {
            double lsq = sq(cR[i]) + sq(n);
            gsl_complex B = gsl_complex_rect(RY[i][n] * std::sqrt(lsq), 0);

            gsl_matrix_complex_set(RRH, 3 * n + i, 3 * n + i, gsl_complex_rect(lsq, 0));
            gsl_matrix_complex_set(RRH, 3 * n + i, i, B);
            gsl_matrix_complex_set(RRH, i, 3 * n + i, B);
        }
    }
}

/**
 * @brief will be used when sigma > TOLERANCE, at the situation, low-pass filter will be applied. Multiple dispatch version of probability matrix calculation for LED. See more information in globes manual at 11.5.3.
 *
 * @param P probability matrix
 * @param cp_sign
 * @param E energy [GeV]
 * @param length baseline [km]
 * @param density [g/cm^3]
 * @param sigma low-pass filter sigma
 * @return int
 */
inline int prob_matrix_5dnu(double (*P)[3], const int cp_sign, const double E, const double length, const double density, const double sigma) {
    // E [GeV], length [km], density [g/cm^3]

    const double R_eVinv = mum_to_eVinv(R);       // [eV^{-1}]
    const double L_GeVinv = mum_to_eVinv(length); // [GeV^{-1}]

#ifndef HAVE_CUDA
    gsl_eigen_hermv_workspace* work = gsl_eigen_hermv_alloc(3 * modes);
#endif
    gsl_vector* e_val = gsl_vector_alloc(3 * modes);
    gsl_matrix_complex* U = gsl_matrix_complex_alloc(3, 3);
    gsl_matrix_complex* RRH = gsl_matrix_complex_alloc(3 * modes, 3 * modes);
    gsl_matrix_complex* e_vec = gsl_matrix_complex_alloc(3 * modes, 3 * modes);
    gsl_matrix_complex* exp_mat = gsl_matrix_complex_alloc(3 * modes, 3 * modes);
    gsl_matrix_complex* L = gsl_matrix_complex_alloc(3, 3 * modes);
    gsl_matrix_complex* UL = gsl_matrix_complex_alloc(3, 3 * modes);
    gsl_matrix_complex* ULexp = gsl_matrix_complex_alloc(3, 3 * modes);
    gsl_matrix_complex* A = gsl_matrix_complex_alloc(3, 3);

    gsl_vector_set_zero(e_val);
    gsl_matrix_complex_set_zero(U);
    gsl_matrix_complex_set_zero(RRH);
    gsl_matrix_complex_set_zero(e_vec);
    gsl_matrix_complex_set_zero(exp_mat);
    gsl_matrix_complex_set_zero(L);
    gsl_matrix_complex_set_zero(UL);
    gsl_matrix_complex_set_zero(ULexp);
    gsl_matrix_complex_set_zero(A);

    make_PMNS(U, cp_sign);
    make_RRH(RRH, U, E, cp_sign, density);

#ifdef HAVE_CUDA
    const int n = 3 * modes;
    cuda_eigen_hermv_gsl(RRH, e_val, e_vec, n);
#else
    gsl_eigen_hermv(RRH, e_val, e_vec, work); // RRH is destroyed here.
#endif

    // 1. convert eigenvectors to Flavor-Mass mix matrix UL
    // UL[alpha][k] means flavor eigenstate alpha mix with mass eigenstate k
    for (int i = 0; i < 3; i++) {
        for (int n = 0; n < 3 * modes; n++) {
            gsl_matrix_complex_set(L, i, n, gsl_matrix_complex_get(e_vec, i, n));
        }
    }
    // UL = U * L (PMNS * Eigenvectors)
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, U, L, GSL_COMPLEX_ZERO, UL);

    // 2. compute phase of every eigenstate
    std::vector<double> phases(3 * modes);
    for (int n = 0; n < 3 * modes; n++) {
        phases[n] = L_GeVinv * gsl_vector_get(e_val, n) / (2 * E * sq(R_eVinv));
    }

    // 3. compute probability, apply low-pass filter
    // formula: P_ab = Sum_{i,j} U_ai U*_bi U*_aj U_bj * exp(-i * dPhi) * Damping

    // compute constant term of filter coefficients: -0.5 * (sigma/E)^2
    const double filter_const = -0.5 * sq(sigma / E);

    for (int alpha = 0; alpha < 3; alpha++) {  // inital flavor
        for (int beta = 0; beta < 3; beta++) { // final flavor
            double prob_sum = 0.0;

            for (int i = 0; i < 3 * modes; i++) {
                // get mix matrix element U_alpha_i and U_beta_i
                const gsl_complex U_ai = gsl_matrix_complex_get(UL, alpha, i);
                const gsl_complex U_bi = gsl_matrix_complex_get(UL, beta, i);

                // diagonal elements need not filter
                // herm P = |Sum|^2 = Sum |term|^2 + 2 Re(Sum cross_terms)
                const double term_diag = gsl_complex_abs2(U_ai) * gsl_complex_abs2(U_bi);
                prob_sum += term_diag;

                // interference term (j > i)
                for (int j = i + 1; j < 3 * modes; j++) {
                    const gsl_complex U_aj = gsl_matrix_complex_get(UL, alpha, j);
                    const gsl_complex U_bj = gsl_matrix_complex_get(UL, beta, j);

                    // compute phase difference
                    // dPhi = phi_i - phi_j
                    double dPhi = phases[j] - phases[i];

                    // apply low-pass filter
                    // dPhi(\Phi_{ij} in manual) is already contains L/2E factor
                    // Damping = exp( -0.5 * (sigma/E * dPhi)^2 )
                    const double damping = std::exp(filter_const * sq(dPhi));

                    // if damping is too small (< 1e-9), skip this term to speed up
                    if (damping < 1e-9) continue;

                    // compute : J = U_ai * conj(U_bi) * conj(U_aj) * U_bj
                    // formula: P = |A|^2 = A * A*
                    // Term_ij = (U_ai U*_bi) * (U_aj U*_bj)* * e^{-i(phi_j - phi_i)}
                    //         = U_ai U*_bi U*_aj U_bj * (cos(dPhi) - i*sin(dPhi))

                    const gsl_complex J = gsl_complex_mul(
                        gsl_complex_mul(U_ai, gsl_complex_conjugate(U_bi)),
                        gsl_complex_mul(gsl_complex_conjugate(U_aj), U_bj));

                    const double cos_phi = std::cos(dPhi);
                    const double sin_phi = std::sin(dPhi);

                    // get real part: Re(J * (cos - i*sin)) = Re(J)*cos + Im(J)*sin
                    const double real_part = GSL_REAL(J) * cos_phi + GSL_IMAG(J) * sin_phi;

                    // double because j<i and j>i are conjugated
                    // times 2 filter factor
                    prob_sum += 2.0 * real_part * damping;
                }
            }
            P[alpha][beta] = prob_sum;
        }
    }

#ifndef HAVE_CUDA
    gsl_eigen_hermv_free(work);
#endif
    gsl_vector_free(e_val);
    gsl_matrix_complex_free(U);
    gsl_matrix_complex_free(RRH);
    gsl_matrix_complex_free(e_vec);
    gsl_matrix_complex_free(exp_mat);
    gsl_matrix_complex_free(L);
    gsl_matrix_complex_free(UL);
    gsl_matrix_complex_free(ULexp);
    gsl_matrix_complex_free(A);

    return EXIT_SUCCESS;
}

inline int prob_matrix_5dnu(double (*P)[3], const int cp_sign, const double E, const double length, const double density) {
    // E [GeV], length [km], density [g/cm^3]

    double R_eVinv = mum_to_eVinv(R);       // [eV^{-1}]
    double L_GeVinv = mum_to_eVinv(length); // [GeV^{-1}]

#ifndef HAVE_CUDA
    gsl_eigen_hermv_workspace* work = gsl_eigen_hermv_alloc(3 * modes);
#endif
    gsl_vector* e_val = gsl_vector_alloc(3 * modes);
    gsl_matrix_complex* U = gsl_matrix_complex_alloc(3, 3);
    gsl_matrix_complex* RRH = gsl_matrix_complex_alloc(3 * modes, 3 * modes);
    gsl_matrix_complex* e_vec = gsl_matrix_complex_alloc(3 * modes, 3 * modes);
    gsl_matrix_complex* exp_mat = gsl_matrix_complex_alloc(3 * modes, 3 * modes);
    gsl_matrix_complex* L = gsl_matrix_complex_alloc(3, 3 * modes);
    gsl_matrix_complex* UL = gsl_matrix_complex_alloc(3, 3 * modes);
    gsl_matrix_complex* ULexp = gsl_matrix_complex_alloc(3, 3 * modes);
    gsl_matrix_complex* A = gsl_matrix_complex_alloc(3, 3);

    gsl_vector_set_zero(e_val);
    gsl_matrix_complex_set_zero(U);
    gsl_matrix_complex_set_zero(RRH);
    gsl_matrix_complex_set_zero(e_vec);
    gsl_matrix_complex_set_zero(exp_mat);
    gsl_matrix_complex_set_zero(L);
    gsl_matrix_complex_set_zero(UL);
    gsl_matrix_complex_set_zero(ULexp);
    gsl_matrix_complex_set_zero(A);

    make_PMNS(U, cp_sign);
    make_RRH(RRH, U, E, cp_sign, density);

#ifdef HAVE_CUDA
    const int n = 3 * modes;
    cuda_eigen_hermv_gsl(RRH, e_val, e_vec, n);
#else
    gsl_eigen_hermv(RRH, e_val, e_vec, work); // RRH is destroyed here.
#endif

    for (int i = 0; i < 3; i++) {
        for (int n = 0; n < 3 * modes; n++) {
            gsl_matrix_complex_set(L, i, n, gsl_matrix_complex_get(e_vec, i, n));
        }
    }

    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, U, L, GSL_COMPLEX_ZERO, UL);

    for (int n = 0; n < 3 * modes; n++) {
        gsl_matrix_complex_set(exp_mat, n, n,
                               gsl_complex_polar(1.0, L_GeVinv * gsl_vector_get(e_val, n) / (2 * E * sq(R_eVinv))));
    }

    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, UL, exp_mat, GSL_COMPLEX_ZERO, ULexp);
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, ULexp, UL, GSL_COMPLEX_ZERO, A);

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            P[a][b] = gsl_complex_abs2(gsl_matrix_complex_get(A, a, b));
        }
    }

#ifndef HAVE_CUDA
    gsl_eigen_hermv_free(work);
#endif
    gsl_vector_free(e_val);
    gsl_matrix_complex_free(U);
    gsl_matrix_complex_free(RRH);
    gsl_matrix_complex_free(e_vec);
    gsl_matrix_complex_free(exp_mat);
    gsl_matrix_complex_free(L);
    gsl_matrix_complex_free(UL);
    gsl_matrix_complex_free(ULexp);
    gsl_matrix_complex_free(A);

    return EXIT_SUCCESS;
}

inline int my_probability_matrix(double P[3][3], const int cp_sign, const double E,
                                 int, const double* length, const double* density,
                                 double filter_sigma, void* user_data) {
    if (filter_sigma > TOLERANCE) {
        prob_matrix_5dnu(P, cp_sign, E, length[0], density[0], filter_sigma);
    } else {
        prob_matrix_5dnu(P, cp_sign, E, length[0], density[0]);
    }
    return EXIT_SUCCESS;
}

} // namespace LED::CalProbability