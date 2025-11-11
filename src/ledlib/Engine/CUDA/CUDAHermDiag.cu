#include "CUDAHermDiag.h++"

#include <cstdio>
#include <cstdlib>
#include <cuComplex.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>

// Helper macros
#define CUDA_CHECK(x)                                                                               \
    do {                                                                                            \
        cudaError_t err = x;                                                                        \
        if (err != cudaSuccess) {                                                                   \
            fprintf(stderr, "CUDA Error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
            return -1;                                                                              \
        }                                                                                           \
    } while (0)
#define CUSOLVER_CHECK(x)                                                  \
    do {                                                                   \
        cusolverStatus_t s = x;                                            \
        if (s != CUSOLVER_STATUS_SUCCESS) {                                \
            fprintf(stderr, "CUSOLVER Error %s:%d\n", __FILE__, __LINE__); \
            return -2;                                                     \
        }                                                                  \
    } while (0)

int cuda_eigen_hermv_gsl(const gsl_matrix_complex* RRH,
                         gsl_vector* e_val,
                         gsl_matrix_complex* e_vec,
                         int n) {
    if (!RRH || !e_val || !e_vec) return -10;
    size_t matrix_size = size_t(n) * size_t(n);

    // Allocate host column-major array
    cuDoubleComplex* h_A = (cuDoubleComplex*)malloc(matrix_size * sizeof(cuDoubleComplex));
    if (!h_A) return -11;

    // Copy GSL row-major -> column-major
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            gsl_complex z = gsl_matrix_complex_get(RRH, i, j);
            h_A[j * n + i] = make_cuDoubleComplex(GSL_REAL(z), GSL_IMAG(z));
        }
    }

    // CUDA memory
    cuDoubleComplex* d_A;
    double* d_W;
    CUDA_CHECK(cudaMalloc((void**)&d_A, matrix_size * sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc((void**)&d_W, n * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_A, h_A, matrix_size * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));

    // cuSOLVER
    cusolverDnHandle_t handle;
    CUSOLVER_CHECK(cusolverDnCreate(&handle));

    int lwork = 0;
    CUSOLVER_CHECK(cusolverDnZheevd_bufferSize(
        handle,
        CUSOLVER_EIG_MODE_VECTOR,
        CUBLAS_FILL_MODE_UPPER,
        n,
        d_A,
        n,
        d_W,
        &lwork));

    cuDoubleComplex* d_work;
    CUDA_CHECK(cudaMalloc((void**)&d_work, lwork * sizeof(cuDoubleComplex)));
    int* devInfo;
    CUDA_CHECK(cudaMalloc((void**)&devInfo, sizeof(int)));

    CUSOLVER_CHECK(cusolverDnZheevd(
        handle,
        CUSOLVER_EIG_MODE_VECTOR,
        CUBLAS_FILL_MODE_UPPER,
        n,
        d_A,
        n,
        d_W,
        d_work,
        lwork,
        devInfo));

    // Copy results back
    double* W_host = (double*)malloc(n * sizeof(double));
    CUDA_CHECK(cudaMemcpy(W_host, d_W, n * sizeof(double), cudaMemcpyDeviceToHost));
    cuDoubleComplex* A_host = (cuDoubleComplex*)malloc(matrix_size * sizeof(cuDoubleComplex));
    CUDA_CHECK(cudaMemcpy(A_host, d_A, matrix_size * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));

    // Write back to GSL
    for (int i = 0; i < n; i++) {
        gsl_vector_set(e_val, i, W_host[i]);
        for (int j = 0; j < n; j++) {
            // gsl_matrix_complex_set(e_vec,i,j,A_host[j*n + i]);
            cuDoubleComplex z = A_host[j * n + i];
            gsl_complex gz = gsl_complex_rect(cuCreal(z), cuCimag(z));
            gsl_matrix_complex_set(e_vec, i, j, gz);
        }
    }

    // Cleanup
    free(h_A);
    free(W_host);
    free(A_host);
    CUDA_CHECK(cudaFree(d_A));
    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_work));
    CUDA_CHECK(cudaFree(devInfo));
    CUSOLVER_CHECK(cusolverDnDestroy(handle));

    return EXIT_SUCCESS;
}
