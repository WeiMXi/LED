#pragma once
#ifdef HAVE_CUDA
#    include "CUDA/CUDAHermDiag.h++"
#endif
#include "globes/globes.h"

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
#include <numbers>
#include <vector>

namespace LED::CalProbability {

inline double sq(const double x) { return x * x; }

constexpr double DELTA = 1.0e-10;
constexpr double DMSQ21 = 7.49e-5;      // [eV^2] NuFit 6.0
constexpr double DMSQ31_NH = 2.513e-3;  // [eV^2] NuFit 6.0
constexpr double DMSQ32_IH = -2.484e-3; // [eV^2] NuFit 6.0

constexpr double GLB_V_factor = 7.5e-14; // Conversion factor for matter potentials
constexpr double GLB_ne_mantle = 0.5;    // Effective electron numbers for calculation

constexpr int GLB_R = 6;
constexpr int GLB_C1R = 7;
constexpr int GLB_C2R = 8;
constexpr int GLB_C3R = 9;
constexpr int GLB_MU1R = 10;
constexpr int GLB_MU2R = 11;
constexpr int GLB_MU3R = 12;

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
static int modes;                                 // modes cutoff in 5d bulk

inline int my_set_oscillation_parameters(glb_params p, void* user_data) {
    th12 = glbGetOscParams(p, GLB_THETA_12);
    th13 = glbGetOscParams(p, GLB_THETA_13);
    th23 = glbGetOscParams(p, GLB_THETA_23);
    deltacp = glbGetOscParams(p, GLB_DELTA_CP);
    sdm = glbGetOscParams(p, GLB_DM_21) * 1.0e-18; /* Convert to GeV^2 */
    ldm = glbGetOscParams(p, GLB_DM_31) * 1.0e-18; /* Convert to GeV^2 */

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

/***************************************************************************
 * Write oscillation parameters from internal data structures into p.      *
 ***************************************************************************/
inline int my_get_oscillation_parameters(glb_params p, void* user_data) {
    glbSetOscParams(p, th12, GLB_THETA_12);
    glbSetOscParams(p, th13, GLB_THETA_13);
    glbSetOscParams(p, th23, GLB_THETA_23);
    glbSetOscParams(p, deltacp, GLB_DELTA_CP);
    glbSetOscParams(p, sdm * 1.0e18, GLB_DM_21); /* Convert to eV^2 */
    glbSetOscParams(p, ldm * 1.0e18, GLB_DM_31); /* Convert to eV^2 */

    /* 5d bulk PARAMETERS*/
    glbSetOscParams(p, R, GLB_R);
    glbSetOscParams(p, c1R, GLB_C1R);
    glbSetOscParams(p, c2R, GLB_C2R);
    glbSetOscParams(p, c3R, GLB_C3R);
    glbSetOscParams(p, mu1R, GLB_MU1R);
    glbSetOscParams(p, mu2R, GLB_MU2R);
    glbSetOscParams(p, mu3R, GLB_MU3R);
    return 0;
}

inline int GetModesCutoff() {
    return modes;
}

inline void SetModesCutoff(const int n) {
    modes = n;
}

inline double mum_to_eVinv(const double x) { return x * 5.06773; }
inline double eVinv_to_mum(const double x) { return x / 5.06773; }

inline double CalMuiR(const double r, const double ciR, const double m12, const double mi2) {
    double r_eVinv = mum_to_eVinv(r);
    return std::sqrt((mi2 - m12) * (r_eVinv) * (r_eVinv) / 2 / std::numbers::pi / std::abs(ciR));
}

inline double CalLightestm2(const double c1R, const double mu1R) {
    return std::exp(-2 * std::numbers::pi * std::abs(c1R)) * (2 * std::numbers::pi * std::abs(c1R)) / mu1R;
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
        if (std::fabs(cR[i]) < DELTA) {
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
    prob_matrix_5dnu(P, cp_sign, E, length[0], density[0]);
    return EXIT_SUCCESS;
}

} // namespace LED::CalProbability