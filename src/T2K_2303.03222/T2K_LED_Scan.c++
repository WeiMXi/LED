/* This program produces the data needed to plot the total events in the T2K far detector.
 * The output consists of five different spectra: CCQE nu_mu, CCQE nu_mu_bar, CCQE nu_e,
 * CCQE nu_e_bar and CC nu_e pi^+.
 *
 *
 * To do list:
 *   Step 1: Uncalibrated events with builtin smearing\\completed
 *   Step 2: Events after correction to normalization \\ignored
 *   Step 3: Events after smearing corrections\\no change
 *   Step 4: Events after data-driven corrections
 */

#include "cstdlib"
#include "ledlib/Engine/ProbabilityEngine.h++"
#include "ledlib/IO/IO.h++" /* my input-output routines */
extern "C" {
#include "T2K_setup.h"
}
#include <glb_smear.h>
#include <glb_types.h>
#include <globes/globes.h> /* GLoBES library */
#include <gsl/gsl_sf_erf.h>
#include <iostream>
#include <math.h>
#include <mpi.h> // 添加MPI头文件
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const std::string MYFILEN1 = "../data/T2K/T2K_LED_Scan_NH.dat";
const std::string MYFILEI1 = "../data/T2K/T2K_LED_Scan_IH.dat";
LED::IO::Output outputFiles;

int main(int argc, char* argv[]) {
    /* Initialize libglobes */
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
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
    glb_params central_values = glbAllocParams();

    /* Define standard oscillation parameters for NO in T2K with Reactor Constraint */
    double theta12 = asin(sqrt(0.307)); // nu-fit 6.0
    double theta13 = asin(sqrt(0.02195));
    double theta23 = asin(sqrt(0.561));
    double deltacp = -1.97;
    double sdm = 7.49e-5;        // nu-fit 5.2
    double ldm = 2.494e-3 + sdm; // NO
    /* Set the parameter vector */
    glbSetOscParams(central_values, 1, LED::CalProbability::GLB_R);
    glbSetOscParams(central_values, 4, LED::CalProbability::GLB_C1R);
    glbSetOscParams(central_values, -4, LED::CalProbability::GLB_C2R);
    glbSetOscParams(central_values, -4, LED::CalProbability::GLB_C3R);
    glbSetOscParams(central_values, 0.01, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(10, 4, sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(10, 4, ldm), LED::CalProbability::GLB_MU3R); // T2K
    LED::CalProbability::SetModesCutoff(20);

    /*Obtained from T2K paper 2303.03222*/
    double theta12_error = 0.75 * M_PI / 180;
    double theta13_error = 1.91e-3;
    double theta23_error = 1.1 * M_PI / 180;
    double deltacp_error = 1.25;
    double sdm_error = 0.19e-5;
    double ldm_error = 0.058e-3;

    glb_params input_errors = glbAllocParams();
    glbDefineParams(input_errors, theta12_error, theta13_error, theta23_error, deltacp_error, sdm_error, ldm_error);
    glbSetOscParams(input_errors, 0, LED::CalProbability::GLB_R);
    glbSetOscParams(input_errors, 0, LED::CalProbability::GLB_C1R);
    glbSetOscParams(input_errors, 0, LED::CalProbability::GLB_C2R);
    glbSetOscParams(input_errors, 0, LED::CalProbability::GLB_C3R);
    glbSetOscParams(input_errors, 0, LED::CalProbability::GLB_MU1R);
    glbSetOscParams(input_errors, 0.01, LED::CalProbability::GLB_MU2R);
    glbSetOscParams(input_errors, 0.02, LED::CalProbability::GLB_MU3R);

    /* Initialize parameter and projection vector(s) */
    glb_params test_values = glbAllocParams();
    glb_params minimum = glbAllocParams();
    glbDefineParams(central_values, theta12, theta13, theta23, deltacp, sdm, ldm);
    glb_projection T2K_projection = glbAllocProjection();

    /* Set central values and Input_error*/
    glbSetDensityParams(central_values, 1.0, GLB_ALL);
    glbSetDensityParams(input_errors, 0.05, GLB_ALL);
    glbSetCentralValues(central_values);
    glbSetInputErrors(input_errors);

    /*Set up the Projection  */
    glbDefineProjection(T2K_projection, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED);
    glbSetProjectionFlag(T2K_projection, GLB_FIXED, LED::CalProbability::GLB_R);
    glbSetProjectionFlag(T2K_projection, GLB_FIXED, LED::CalProbability::GLB_C1R);
    glbSetProjectionFlag(T2K_projection, GLB_FIXED, LED::CalProbability::GLB_C2R);
    glbSetProjectionFlag(T2K_projection, GLB_FIXED, LED::CalProbability::GLB_C3R);
    glbSetProjectionFlag(T2K_projection, GLB_FIXED, LED::CalProbability::GLB_MU1R);
    glbSetProjectionFlag(T2K_projection, GLB_FIXED, LED::CalProbability::GLB_MU2R);
    glbSetProjectionFlag(T2K_projection, GLB_FIXED, LED::CalProbability::GLB_MU3R);
    glbSetDensityProjectionFlag(T2K_projection, GLB_FIXED, GLB_ALL);
    glbSetProjection(T2K_projection);

    /* The oscillation probabilities are computed */
    glbSetOscillationParameters(central_values);
    glbCopyParams(central_values, test_values);
    glbSetRates();

    //////////////////////////////////////////////
    /// 		Calculate data		      ///
    /////////////////////////////////////////////

    /* Initiate a parameter vector for the scan */
    glbCopyParams(central_values, test_values);

    double xmin = 4;
    double xmax = 10;
    int xsteps = 40;
    double ymin = 0.1;
    double ymax = 0.5;
    int ysteps = 20;
    int total_tasks = xsteps * ysteps;
    double dx = (xmax - xmin) / xsteps;
    double dy = (ymax - ymin) / ysteps;

    int quotient = total_tasks / size;
    int remainder = total_tasks % size;
    int start_task, num_tasks;
    if (rank < remainder) {
        start_task = rank * (quotient + 1);
        num_tasks = quotient + 1;
    } else {
        start_task = remainder * (quotient + 1) + (rank - remainder) * quotient;
        num_tasks = quotient;
    }
    if (num_tasks == 0) num_tasks = 1;
    /* result */
    double* local_res = (double*)malloc(num_tasks * sizeof(double));
    double local_x, local_y;
    double theAbsCR, themu1R, theLightestm2, theR, themu3R;
    double res;
    double local_chi_min = 1000000.0;
    double local_theta23_min = 0.0;
    double local_deltacp_min = 0.0;

    double start_time = MPI_Wtime();
    theR = 10;
    /* MPI */
    for (int t = 0; t < num_tasks; t++) {
        int task_idx = start_task + t;
        int x_idx = task_idx / ysteps;
        int y_idx = task_idx % ysteps;
        local_x = xmin + x_idx * dx;
        local_y = ymin + y_idx * dy;
        /* Set vector of test values */
        theAbsCR = local_x;
        themu1R = local_y;

        // double hc1 = theAbsCR * theAbsCR + (std::numbers::pi * theAbsCR - 1) * themu1R;
        // double h0 = LED::CalProbability::hi_zero(theAbsCR, themu1R);
        // std::cout << hc1 << " " << h0 << std::endl;

        glbSetOscParams(test_values, theR, LED::CalProbability::GLB_R);
        glbSetOscParams(test_values, theAbsCR, LED::CalProbability::GLB_C1R);
        glbSetOscParams(test_values, -theAbsCR, LED::CalProbability::GLB_C2R);
        glbSetOscParams(test_values, -theAbsCR, LED::CalProbability::GLB_C3R);

        glbSetOscParams(test_values, themu1R, LED::CalProbability::GLB_MU1R);
        glbSetOscParams(test_values, LED::CalProbability::CalculateMuiR(theR, -theAbsCR, sdm), LED::CalProbability::GLB_MU2R);
        glbSetOscParams(test_values, LED::CalProbability::CalculateMuiR(theR, -theAbsCR, ldm), LED::CalProbability::GLB_MU3R);
        // std::cout << sqrt(LED::CalProbability::CalLightestm2(theR, theAbsCR, themu1R)) * LED::CalProbability::mum_to_eVinv(theR) << std::endl;
        res = glbChiNP(test_values, minimum, GLB_ALL);
        printf("%f %f %f\n", theAbsCR, themu1R, res);
        local_res[t] = res;
        double local_elapsed = MPI_Wtime() - start_time;

        int local_completed = t + 1;
        int global_completed = local_completed * size;
        double tasks_per_sec = (double)global_completed / local_elapsed;
        int remaining_tasks = total_tasks - global_completed;
        double remaining_time = remaining_tasks / tasks_per_sec;

        if (rank == 0) {
            printf("Progress: %d/%d tasks completed. average time: %.2f /s, %.2f s left.\n",
                   global_completed, total_tasks, tasks_per_sec, remaining_time);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    double global_chi_min, global_theta23_min, global_deltacp_min;
    MPI_Allreduce(&local_chi_min, &global_chi_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    double* all_res = NULL;
    if (rank == 0) {
        all_res = (double*)malloc(total_tasks * sizeof(double));
    }

    int* all_num_tasks = NULL;
    if (rank == 0) {
        all_num_tasks = (int*)malloc(size * sizeof(int));
    }
    MPI_Gather(&num_tasks, 1, MPI_INT, all_num_tasks, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int *displs = NULL, *recvcounts = NULL;
    if (rank == 0) {
        displs = (int*)malloc(size * sizeof(int));
        recvcounts = (int*)malloc(size * sizeof(int));
        displs[0] = 0;
        recvcounts[0] = all_num_tasks[0];
        for (int i = 1; i < size; i++) {
            displs[i] = displs[i - 1] + all_num_tasks[i - 1];
            recvcounts[i] = all_num_tasks[i];
        }
    }
    MPI_Gatherv(local_res, num_tasks, MPI_DOUBLE, all_res, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        outputFiles.InitOutput(MYFILEN1, "");
        double chi_min = global_chi_min;
        double theta23_min = 0.0;
        double deltacp_min = 0.0;
        int z = 0;
        for (int idx = 0; idx < total_tasks; idx++) {
            int x_idx = idx / ysteps;
            int y_idx = idx % ysteps;
            double x_out = xmin + x_idx * dx;
            double y_out = ymin + y_idx * dy;
            res = all_res[idx];
            outputFiles.AddToOutput(x_out, y_out, res);
            if (res < chi_min) {
                chi_min = res;
                theta23_min = asin(sqrt(x_out));
                deltacp_min = y_out;
            }
            z++;
            if (z % 100 == 0) {
                printf("%d\n", z);
            }
        }
    }

    // Flip hierarchy
    glbSetOscParams(central_values, asin(sqrt(0.56)), GLB_THETA_23);
    glbSetOscParams(central_values, -1.44, GLB_DELTA_CP);
    glbSetOscParams(central_values, 1, LED::CalProbability::GLB_R);
    glbSetOscParams(central_values, -10, LED::CalProbability::GLB_C1R);
    glbSetOscParams(central_values, -10, LED::CalProbability::GLB_C2R);
    glbSetOscParams(central_values, 10, LED::CalProbability::GLB_C3R);

    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(1, -10, 2.463e-3), LED::CalProbability::GLB_MU1R);
    glbSetOscParams(central_values, LED::CalProbability::CalculateMuiR(1, -10, 2.463e-3 + sdm), LED::CalProbability::GLB_MU2R);
    glbSetOscParams(central_values, 0.01, LED::CalProbability::GLB_MU3R);
    glbSetOscillationParameters(central_values);
    glbSetRates();

    glbCopyParams(central_values, test_values);
    start_time = MPI_Wtime();
    local_chi_min = 1000000.0;
    local_theta23_min = 0.0;
    local_deltacp_min = 0.0;

    /* MPI */
    for (int t = 0; t < num_tasks; t++) {
        int task_idx = start_task + t;
        int x_idx = task_idx / ysteps;
        int y_idx = task_idx % ysteps;
        local_x = xmin + x_idx * dx;
        local_y = ymin + y_idx * dy;
        /* Set vector of test values */
        theAbsCR = local_x;
        themu3R = local_y;

        glbSetOscParams(test_values, theR, LED::CalProbability::GLB_R);
        glbSetOscParams(test_values, -theAbsCR, LED::CalProbability::GLB_C1R);
        glbSetOscParams(test_values, -theAbsCR, LED::CalProbability::GLB_C2R);
        glbSetOscParams(test_values, theAbsCR, LED::CalProbability::GLB_C3R);

        glbSetOscParams(test_values, LED::CalProbability::CalculateMuiR(theR, -theAbsCR, 2.463e-3), LED::CalProbability::GLB_MU1R);
        glbSetOscParams(test_values, LED::CalProbability::CalculateMuiR(theR, -theAbsCR, 2.463e-3 + sdm), LED::CalProbability::GLB_MU2R);
        glbSetOscParams(test_values, themu3R, LED::CalProbability::GLB_MU3R);
        // std::cout << sqrt(LED::CalProbability::CalLightestm2(10, theAbsCR, themu1R)) * LED::CalProbability::mum_to_eVinv(10) << std::endl;
        res = glbChiNP(test_values, minimum, GLB_ALL);

        printf("%f %f %f\n", theAbsCR, themu3R, res);
        local_res[t] = res;
        double local_elapsed = MPI_Wtime() - start_time;

        int local_completed = t + 1;
        int global_completed = local_completed * size;
        double tasks_per_sec = (double)global_completed / local_elapsed;
        int remaining_tasks = total_tasks - global_completed;
        double remaining_time = remaining_tasks / tasks_per_sec;

        if (rank == 0) {
            printf("Progress: %d/%d tasks completed. average time: %.2f /s, %.2f s left.\n",
                   global_completed, total_tasks, tasks_per_sec, remaining_time);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(&local_chi_min, &global_chi_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    all_res = NULL;
    if (rank == 0) {
        all_res = (double*)malloc(total_tasks * sizeof(double));
    }
    // 先gather num_tasks到rank0
    all_num_tasks = NULL;
    if (rank == 0) {
        all_num_tasks = (int*)malloc(size * sizeof(int));
    }
    MPI_Gather(&num_tasks, 1, MPI_INT, all_num_tasks, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 然后displs和recvcounts for irregular gather
    displs = NULL,
    recvcounts = NULL;
    if (rank == 0) {
        displs = (int*)malloc(size * sizeof(int));
        recvcounts = (int*)malloc(size * sizeof(int));
        displs[0] = 0;
        recvcounts[0] = all_num_tasks[0];
        for (int i = 1; i < size; i++) {
            displs[i] = displs[i - 1] + all_num_tasks[i - 1];
            recvcounts[i] = all_num_tasks[i];
        }
    }
    MPI_Gatherv(local_res, num_tasks, MPI_DOUBLE, all_res, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        outputFiles.InitOutput(MYFILEI1, "");
        double chi_min = global_chi_min;
        double theta23_min = 0.0;
        double deltacp_min = 0.0;
        int z = 0;
        for (int idx = 0; idx < total_tasks; idx++) {
            int x_idx = idx / ysteps;
            int y_idx = idx % ysteps;
            double x_out = xmin + x_idx * dx;
            double y_out = ymin + y_idx * dy;
            res = all_res[idx];
            outputFiles.AddToOutput(x_out, y_out, res);
        }
    }

    /* clean */
    free(local_res);
    if (rank == 0) {
        free(all_res);
        free(all_num_tasks);
        free(displs);
        free(recvcounts);
    }

    MPI_Finalize();
    // Destroy parameter vector(s)
    glbFreeParams(central_values);
    glbFreeParams(test_values);
    glbFreeParams(input_errors);
    glbFreeParams(minimum);
    glbFreeProjection(T2K_projection);

    exit(0);
}
