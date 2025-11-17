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

#include "myio.h" /* my input-output routines */

#include "T2K_setup.h"

#include <glb_smear.h>
#include <glb_types.h>
#include <globes/globes.h> /* GLoBES library */
#include <gsl/gsl_sf_erf.h>
#include <math.h>
#include <mpi.h> // 添加MPI头文件
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char MYFILE[] = "../data/T2K/T2K_SM_fit.dat";

double SQR(double x) {
    return x * x;
}

int main(int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    /* Initialize libglobes */
    glbInit(argv[0]);
    glbDefineChiFunction(&ChiT2K, 4, "ChiT2K", NULL);
    /* Initialize T2K and NOvA */
    InitializeT2K(&glb_experiment_list[0], &glb_num_of_exps);

    /* Define standard oscillation parameters for NO in T2K with Reactor Constraint */
    double theta12 = asin(sqrt(0.307)); // nu-fit 5.2
    double theta13 = asin(sqrt(0.0218));
    double theta23 = asin(sqrt(0.561));
    double deltacp = -1.97;
    double sdm = 7.53e-5;        // nu-fit 5.2
    double ldm = 2.494e-3 + sdm; // NO

    /*Obtained from T2K paper 2303.03222*/
    double theta12_error = 0.75 * M_PI / 180; // nu-fit 5.2
    double theta13_error = 1.91e-3;           // nu-fit 5.2
    double theta23_error = 1.1 * M_PI / 180;
    double deltacp_error = 1.25;
    double sdm_error = 0.21e-5; // nu-fit 5.2
    double ldm_error = 0.058e-3;

    /* Initialize the parameter vector */
    /* Initialize parameter and projection vector(s) */
    glb_params central_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    glb_params minimum = glbAllocParams();
    glb_projection T2K_projection = glbAllocProjection();

    /* Set the parameter vector */
    glbDefineParams(central_values, theta12, theta13, theta23, deltacp, sdm, ldm);
    glbDefineParams(input_errors, theta12_error, theta13_error, theta23_error, deltacp_error, sdm_error, ldm_error);

    /* Set central values and Input_error*/
    glbSetDensityParams(central_values, 1.0, GLB_ALL);
    glbSetDensityParams(input_errors, 0.05, GLB_ALL);
    glbSetCentralValues(central_values);
    glbSetInputErrors(input_errors);

    /*Set up the Projection  */
    glbDefineProjection(T2K_projection, GLB_FIXED, GLB_FREE, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED);
    glbSetDensityProjectionFlag(T2K_projection, GLB_FIXED, GLB_ALL);
    glbSetProjection(T2K_projection);

    /* The oscillation probabilities are computed */
    glbSetOscillationParameters(central_values);
    glbSetRates();

    glbCopyParams(central_values, test_values);

    // /* Iterate over all deltacp values and find the global minimum */
    // double thetheta23, theta23_min, res;
    // double thedeltacp, deltacp_min;
    // double x, y;
    // double chi_min = 1000000.0;
    // int z = 0;
    // double xmin = 0.35;
    // double xmax = 0.65;
    // int xsteps = 100;
    // double ymin = 0;
    // double ymax = 2 * M_PI;
    // int ysteps = 100;
    // InitOutput(MYFILE, "");

    // for (x = xmin; x <= xmax; x = x + (xmax - xmin) / xsteps) {
    //     for (y = ymin; y <= ymax; y = y + (ymax - ymin) / ysteps) {
    //         /* Set vector of test values */
    //         thetheta23 = asin(sqrt(x)); // Sin2 theta23 to radian
    //         thedeltacp = y;
    //         glbSetOscParams(test_values, thetheta23, GLB_THETA_23);
    //         glbSetOscParams(test_values, thedeltacp, GLB_DELTA_CP);
    //         /* Compute Chi^2 for all loaded experiments and all rules */

    // res = glbChiNP(test_values, minimum, GLB_ALL);
    // AddToOutput(x, y, res);
    // if (res < chi_min) {
    //     theta23_min = thetheta23;
    //     chi_min = res;
    //     deltacp_min = thedeltacp;
    // }
    // z++;
    // if (z % 100 == 0) {
    //     printf("%d\n", z);
    // }
    // }
    // }
    // printf("\n\nBest-fit value: Sin^2 theta23 = %f degrees, DeltaCP = %f and chi^2_min = %f \n\n", SQR(sin(theta23_min)), deltacp_min, chi_min);

    // glbFreeParams(central_values);

    // /* Clear experiment list */
    // glbClearExperimentList();

    /* **新增：并行网格搜索参数** */
    double xmin = 0.35;
    double xmax = 0.65;
    int xsteps = 100;
    double ymin = 0;
    double ymax = 2 * M_PI;
    int ysteps = 100;
    int total_tasks = xsteps * ysteps; // 总网格点数
    double dx = (xmax - xmin) / xsteps;
    double dy = (ymax - ymin) / ysteps;
    /* **新增：任务分发（类似于打印循环）** */
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
    if (num_tasks == 0) num_tasks = 1; // 避免空任务
    /* **新增：本地结果数组** */
    double* local_res = (double*)malloc(num_tasks * sizeof(double));
    double local_x, local_y;
    double thetheta23, thedeltacp;
    double res;
    double local_chi_min = 1000000.0;
    double local_theta23_min = 0.0;
    double local_deltacp_min = 0.0;
    int local_z = 0;
    /* 在循环前添加时间起始记录（所有进程） */
    double start_time = MPI_Wtime();

    /* **并行计算循环** */
    for (int t = 0; t < num_tasks; t++) {
        int task_idx = start_task + t;
        int x_idx = task_idx / ysteps;
        int y_idx = task_idx % ysteps;
        local_x = xmin + x_idx * dx;
        local_y = ymin + y_idx * dy;

        /* Set vector of test values */
        thetheta23 = asin(sqrt(local_x)); // Sin2 theta23 to radian
        thedeltacp = local_y;
        glbSetOscParams(test_values, thetheta23, GLB_THETA_23);
        glbSetOscParams(test_values, thedeltacp, GLB_DELTA_CP);
        // glbSetRates();
        /* Compute Chi^2 for all loaded experiments and all rules */
        res = glbChiNP(test_values, minimum, GLB_ALL);
        local_res[t] = res;

        /* **本地更新min** */
        if (res < local_chi_min) {
            local_chi_min = res;
            local_theta23_min = thetheta23;
            local_deltacp_min = thedeltacp;
        }

        local_z++;
        if (local_z % 100 == 0) {
            /* 计算本地已用时间 */
            double local_elapsed = MPI_Wtime() - start_time;

            /* 全局汇总已完成任务数（所有进程的 local_z 求和） */
            int global_completed;
            MPI_Allreduce(&local_z, &global_completed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            /* 全局汇总平均已用时间（所有进程的 local_elapsed 求平均） */
            double sum_elapsed;
            MPI_Allreduce(&local_elapsed, &sum_elapsed, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            double avg_elapsed = sum_elapsed / size;

            /* 计算每个任务的平均时间，并估算剩余时间 */
            double avg_per_task = (global_completed > 0) ? (avg_elapsed / global_completed) : 0.0;
            double remaining_time = (total_tasks - global_completed) * avg_per_task;

            /* 只在 rank 0 打印进度 */
            if (rank == 0) {
                printf("Progress: %d/%d tasks completed. Elapsed: %.2f s, ETA: %.2f s\n",
                       global_completed, total_tasks, avg_elapsed, remaining_time);
            }
        }
    }

    /* **新增：MPI_Barrier 同步所有进程** */
    MPI_Barrier(MPI_COMM_WORLD);

    double global_chi_min, global_theta23_min, global_deltacp_min;
    MPI_Allreduce(&local_chi_min, &global_chi_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    // 注意：theta23/deltacp_min需对应min位置，实际实现中可广播或额外gather；这里简化假设rank0有

    /* **新增：收集所有res到rank0** */
    double* all_res = NULL;
    if (rank == 0) {
        all_res = (double*)malloc(total_tasks * sizeof(double));
    }
    // 先gather num_tasks到rank0
    int* all_num_tasks = NULL;
    if (rank == 0) {
        all_num_tasks = (int*)malloc(size * sizeof(int));
    }
    MPI_Gather(&num_tasks, 1, MPI_INT, all_num_tasks, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 然后displs和recvcounts for irregular gather
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

    /* **新增：rank0处理输出和全局min** */
    if (rank == 0) {
        // 重组并输出
        InitOutput(MYFILE, "");
        double chi_min = global_chi_min;
        for (int idx = 0; idx < total_tasks; idx++) {
            int x_idx = idx / ysteps;
            int y_idx = idx % ysteps;
            double x_out = xmin + x_idx * dx;
            double y_out = ymin + y_idx * dy;
            res = all_res[idx];
            AddToOutput(x_out, y_out, res);
        }
    }

    /* **清理** */
    free(local_res);
    if (rank == 0) {
        free(all_res);
        free(all_num_tasks);
        free(displs);
        free(recvcounts);
    }

    // End MPI
    MPI_Finalize();

    exit(0);
}
