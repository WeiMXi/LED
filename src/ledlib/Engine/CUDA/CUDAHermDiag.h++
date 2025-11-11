#pragma once

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector.h>
#include <cuComplex.h>

/*
 * Hermitian eigen-decomposition using cuSOLVER
 *
 * Inputs:
 *   RRH    : gsl_matrix_complex* (n x n) - input Hermitian matrix
 *   n      : dimension (should match RRH->size1)
 *
 * Outputs:
 *   e_val  : gsl_vector* (size n) - eigenvalues ascending
 *   e_vec  : gsl_matrix_complex* (n x n) - eigenvectors (columns)
 *
 * Returns 0 on success, non-zero on failure.
 */
int cuda_eigen_hermv_gsl(const gsl_matrix_complex* RRH,
                         gsl_vector* e_val,
                         gsl_matrix_complex* e_vec,
                         int n);
