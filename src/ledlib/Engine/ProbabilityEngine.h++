#pragma once

#ifdef HAVE_CUDA
#    include "CUDA/CUDAHermDiag.h++"
#endif
#include "globes/globes.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <format>
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

// deprecated
constexpr int GLB_R = 6;
constexpr int GLB_C1R = 7;
constexpr int GLB_C2R = 8;
constexpr int GLB_C3R = 9;
constexpr int GLB_MU1R = 10;
constexpr int GLB_MU2R = 11;
constexpr int GLB_MU3R = 12;

enum ParIndexes {
    LedParR = 6,
    LedParC1R = 7,
    LedParC2R = 8,
    LedParC3R = 9,
    LedParMU1R = 10,
    LedParMU2R = 11,
    LedParMU3R = 12,
    LedParM0 = 10 // do not compatible with LedParMUXR
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
static double m0;                                 // the lightest nu mass
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

inline double compute_muiR(double r, double ciR, double mi_sq) {
    double r_eVinv = mum_to_eVinv(r);
    double abs_c = fabs(ciR);
    double pi = std::numbers::pi;
    double sinh_p = sinh(pi * abs_c);
    double coth_p = cosh(pi * abs_c) / sinh_p;
    double num_base = 2 * pi * abs_c * (abs_c * coth_p - ciR);
    double denom_base = 2 * abs_c;
    double denom_c = pi * coth_p;
    double denom_d = pi * pi * abs_c / (sinh_p * sinh_p);
    double y = r_eVinv * r_eVinv * mi_sq; // (R * m_i)^2
    double res = y * denom_base / (num_base + denom_d * y - denom_c * y);
    double muiR = sqrt(res);
    return muiR;
}

/**
 * @brief get the dirac mass form neutrino mass [eV] (IN C = 0!, with second-order correction)
 *
 * @param mui neutrino mass [eV]
 * @param r_mum compactification radius [mum]
 * @return double, dirac mass [eV]
 */
auto CalMuid(const double mui, const double r_mum) -> double {
    return mui / (1 - (std::pow(std::numbers::pi * mum_to_eVinv(r_mum) * mui, 2) / 6));
};

/**
 * @brief get the three muR form the lightest neutrino mass [eV] (IN C = 0!, with second-order correction)
 *
 * @param m0 the lightest neutrino mass [eV]
 * @param r_mum compactification radius [mum]
 * @return std::array<double, 3>, mu1R, mu2R, mu3R
 */
auto calThreeMuiR(const double m0, const double r_mum) -> std::array<double, 3> {
    const double r_eVinv{LED::CalProbability::mum_to_eVinv(r_mum)};
    return std::array<double, 3>{
        CalMuid(m0, r_mum) * r_eVinv,
        CalMuid(std::sqrt(sdm + m0 * m0), r_mum) * r_eVinv,
        CalMuid(std::sqrt(ldm + m0 * m0), r_mum) * r_eVinv};
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
    mu1R = glbGetOscParams(p, GLB_MU1R);
    mu2R = glbGetOscParams(p, GLB_MU2R);
    mu3R = glbGetOscParams(p, GLB_MU3R);

    return 0;
}

/**
 * @brief Set the lightest neutrino mass and compactification radius, the compute other parameters
 * WARRING: this function need to be modified if you want to set c != 0
 * @param p glb_params
 * @param user_data
 * @return int
 */
inline int SetRMParameters(glb_params p, void* user_data) {
    th12 = glbGetOscParams(p, GLB_THETA_12);
    th13 = glbGetOscParams(p, GLB_THETA_13);
    th23 = glbGetOscParams(p, GLB_THETA_23);
    deltacp = glbGetOscParams(p, GLB_DELTA_CP);
    sdm = glbGetOscParams(p, GLB_DM_21);
    ldm = glbGetOscParams(p, GLB_DM_31);

    /* 5d bulk PARAMETERS*/

    R = glbGetOscParams(p, LedParR);
    c1R = glbGetOscParams(p, LedParC1R);
    c2R = glbGetOscParams(p, LedParC2R);
    c3R = glbGetOscParams(p, LedParC3R);

    m0 = glbGetOscParams(p, LedParM0);

    if (std::abs(c2R) < TOLERANCE) {
        const auto [mu1RTemp, mu2RTemp, mu3RTemp] = calThreeMuiR(m0, R);
        mu1R = mu1RTemp;
        mu2R = mu2RTemp;
        mu3R = mu3RTemp;
    } else { // need modify to c != 0 case
        mu1R = CalMuid(m0, R) * mum_to_eVinv(R);
        mu2R = CalMuiR(R, c2R, sdm);
        mu3R = CalMuiR(R, c3R, ldm);
    }

    return EXIT_SUCCESS;
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
    glbSetOscParams(p, mu1R, GLB_MU1R);
    glbSetOscParams(p, mu2R, GLB_MU2R);
    glbSetOscParams(p, mu3R, GLB_MU3R);

    std::cout << std::format("m0: {:.6f} eV, mu1R: {:.6f}, mu2R: {:.6f}, mu3R: {:.6f}\n", m0, mu1R, mu2R, mu3R);

    return 0;
}

/**
 * @brief get the oscillation parameters, only for use m0 and R situation
 *
 * @param p glb_params
 * @param user_data
 * @return int
 */
inline int GetRMParameters(glb_params p, void* user_data) {
    glbSetOscParams(p, th12, GLB_THETA_12);
    glbSetOscParams(p, th13, GLB_THETA_13);
    glbSetOscParams(p, th23, GLB_THETA_23);
    glbSetOscParams(p, deltacp, GLB_DELTA_CP);
    glbSetOscParams(p, sdm, GLB_DM_21);
    glbSetOscParams(p, ldm, GLB_DM_31);

    /* 5d bulk PARAMETERS*/
    glbSetOscParams(p, R, LedParR);
    glbSetOscParams(p, c1R, LedParC1R);
    glbSetOscParams(p, c2R, LedParC2R);
    glbSetOscParams(p, c3R, LedParC3R);
    glbSetOscParams(p, m0, LedParM0);

    return EXIT_SUCCESS;
}

inline int GetModesCutoff() {
    return modes;
}

inline void SetModesCutoff(const int n) {
    modes = n;
}

double hi_zero(double cbar, double mubar) {
    const double pi = std::numbers::pi_v<double>;

    double abs_c = std::abs(cbar);
    if (abs_c < 1e-12) {
        return 0.0;
    }

    double coth = std::cosh(pi * abs_c) / std::sinh(pi * abs_c);
    double term = cbar - cbar * cbar * coth; // 注意 cbar² = (|c|)^2
    return pi * mubar * mubar * term;
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
 * @param E energy
 * @param length baseline
 * @param density
 * @param sigma low-pass filter sigma
 * @return int
 */
inline int prob_matrix_5dnu(double (*P)[3], const int cp_sign, const double E, const double length, const double density, const double sigma) {
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
    double filter_const = -0.5 * sq(sigma / E);

    for (int alpha = 0; alpha < 3; alpha++) {  // inital flavor
        for (int beta = 0; beta < 3; beta++) { // final flavor
            double prob_sum = 0.0;

            for (int i = 0; i < 3 * modes; i++) {
                // get mix matrix element U_alpha_i and U_beta_i
                gsl_complex U_ai = gsl_matrix_complex_get(UL, alpha, i);
                gsl_complex U_bi = gsl_matrix_complex_get(UL, beta, i);

                // diagonal elements need not filter
                // herm P = |Sum|^2 = Sum |term|^2 + 2 Re(Sum cross_terms)
                double term_diag = gsl_complex_abs2(U_ai) * gsl_complex_abs2(U_bi);
                prob_sum += term_diag;

                // interference term (j > i)
                for (int j = i + 1; j < 3 * modes; j++) {
                    gsl_complex U_aj = gsl_matrix_complex_get(UL, alpha, j);
                    gsl_complex U_bj = gsl_matrix_complex_get(UL, beta, j);

                    // compute phase difference
                    // dPhi = phi_i - phi_j
                    double dPhi = phases[i] - phases[j];

                    // apply low-pass filter
                    // dPhi(\Phi_{ij} in manual) is already contains L/2E factor
                    // Damping = exp( -0.5 * (sigma/E * dPhi)^2 )
                    double damping = std::exp(filter_const * sq(dPhi));

                    // if damping is too small (< 1e-9), skip this term to speed up
                    if (damping < 1e-9) continue;

                    // compute : J = U_ai * conj(U_bi) * conj(U_aj) * U_bj
                    // formula: P = |A|^2 = A * A*
                    // Term_ij = (U_ai U*_bi) * (U_aj U*_bj)* * e^{-i(phi_i - phi_j)}
                    //         = U_ai U*_bi U*_aj U_bj * (cos(dPhi) - i*sin(dPhi))

                    gsl_complex J = gsl_complex_mul(
                        gsl_complex_mul(U_ai, gsl_complex_conjugate(U_bi)),
                        gsl_complex_mul(gsl_complex_conjugate(U_aj), U_bj));

                    double cos_phi = std::cos(dPhi);
                    double sin_phi = std::sin(dPhi);

                    // get real part: Re(J * (cos - i*sin)) = Re(J)*cos + Im(J)*sin
                    double real_part = GSL_REAL(J) * cos_phi + GSL_IMAG(J) * sin_phi;

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