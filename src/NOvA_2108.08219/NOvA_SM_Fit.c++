/* ################################################################################## *
 * 	This program fits sin2th23 and deltaCP to NOvA far detector data alone	 *
 * ################################################################################## *
 */

// #define _GNU_SOURCE
#include <math.h>
#include <mpi.h> // 添加MPI头文件
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define atoa(x) #x

extern "C" {
#include "NOvA_Setup_New.h"
}
#include "ledlib/Engine/ProbabilityEngine.h++"
#include "ledlib/IO/IO.h++" /* my input-output routines */

#include <globes/globes.h> // GLoBES library

char MYFILEN1[] = "../data/NOvA/NOvA_fit_data_normal.dat";
char MYFILEI1[] = "../data/NOvA/NOvA_fit_data_inverted.dat";
LED::IO::Output outputFiles;

static int n_params = 6;

double square(double x) {
    return x * x;
}

double my_prior(const glb_params in, void* user_data) {
    glb_projection p = glbAllocProjection();
    glbGetProjection(p);
    double aux;
    double pv = 0.0;
    double fitvalue, centralvalue, inputerror;

    double th13Prior = asin(sqrt(0.022));

    if (glbGetProjectionFlag(p, GLB_THETA_13) == GLB_FREE) {
        fitvalue = glbGetOscParams(in, GLB_THETA_13);
        fitvalue = square(sin(2 * fitvalue));
        centralvalue = square(sin(2 * th13Prior));
        inputerror = 0.0019;
        pv += square((centralvalue - fitvalue) / inputerror);
    }
    if (glbGetProjectionFlag(p, GLB_DELTA_CP) == GLB_FREE) {
        fitvalue = glbGetOscParams(in, GLB_DELTA_CP);
        centralvalue = 0.82 * M_PI;
        inputerror = 0.27 * M_PI;
        pv += square((centralvalue - fitvalue) / inputerror);
    }

    glbFreeProjection(p);
    return pv;
}

void percent_bar(double* check_bar, double DeltaP, int current, double print_percent) {
    *check_bar = *check_bar + DeltaP;
    int aux_percent;
    if (*check_bar > print_percent) {
        printf("#");

        printf("](%i %%)", current);
        if (current < 10)
            printf("\b\b\b\b\b\b");
        else if (current <= 100)
            printf("\b\b\b\b\b\b\b");
        else
            printf("\b\b\b\b\b\b\b\b");
        fflush(stdout);
        *check_bar = 0.0;
    }
}

void set_osc_params_zero(glb_params in_params) {
    for (int i = 0; i < n_params; i++)
        glbSetOscParams(in_params, 0.0, i);
}

///////////////////////////////////////////////
///		Define projection           ///
//////////////////////////////////////////////

/* Set all proj flags to GLB_FIXED */
void set_proj_all_fixed(glb_projection in_proje) {
    for (int i = 0; i < n_params; i++)
        glbSetProjectionFlag(in_proje, GLB_FIXED, i);

    glbSetDensityProjectionFlag(in_proje, GLB_FIXED, GLB_ALL);
}

/* Set up a projection on thet23-deltaCP plane */
void set_proj_SM_TH23_DCP_2D(glb_projection in_proje) {
    // set_proj_all_fixed(in_proje);
    /* THETA_12, THETA_13, THETA_23, DELTA_CP,  DM_21,     DM_31 */
    glbDefineProjection(in_proje, GLB_FIXED, GLB_FREE, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FREE);
    glbSetDensityProjectionFlag(in_proje, GLB_FIXED, GLB_ALL);
}

////////////////////////////////////////////////////////////
///		       MAIN PROGRAM 		           //
///////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    /* Initialize libglobes */
    glbInit(argv[0]);

    /* Initialize experiments */
    InitializeNOvA(&glb_experiment_list[0], &glb_num_of_exps);

    /* Select prior function (defined above) */
    glbRegisterPriorFunction(my_prior, NULL, NULL, NULL);

    /* Define standard oscillation parameters */
    double theta12 = asin(sqrt(0.307));
    double theta13 = asin(sqrt(0.02195));
    double theta23 = asin(sqrt(0.547));
    double deltacp = 0.87 * M_PI;
    double sdm = 7.49e-5;
    double ldm = 2.441e-3 + sdm;

    /* Initialize parameter vectors */
    glb_params central_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    glb_params minimum = glbAllocParams();

    /* Set central values */
    set_osc_params_zero(central_values);
    glbDefineParams(central_values, theta12, theta13, theta23, deltacp, sdm, ldm);
    glbSetDensityParams(central_values, 1.0, GLB_ALL);

    /* Set prior values (no priors, only 5% uncertainty for matter density) */
    set_osc_params_zero(input_errors);
    glbSetDensityParams(input_errors, 0.05, GLB_ALL);

    /* Insert central values and input errors into GLoBES */
    glbSetCentralValues(central_values);
    glbSetInputErrors(input_errors);

    /* The simulated data is computed to acquire correct background */
    glbSetOscillationParameters(central_values);
    glbSetRates();

    /* Set projection as previously defined */
    glb_projection my_projection = glbAllocProjection();
    set_proj_SM_TH23_DCP_2D(my_projection); /* SM fit onto TH23-DCP plane */

    glbSetProjection(my_projection);

    //////////////////////////////////////////////
    /// 		Calculate data		      ///
    /////////////////////////////////////////////

    /* Initiate a parameter vector for the scan */
    glbCopyParams(central_values, test_values);

    double xmin = 0.35;
    double xmax = 0.65;
    int xsteps = 101;
    double ymin = 0;
    double ymax = 2 * M_PI;
    int ysteps = 101;
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
    double thetheta23, thedeltacp;
    double res;
    double local_chi_min = 1000000.0;
    double local_theta23_min = 0.0;
    double local_deltacp_min = 0.0;

    double start_time = MPI_Wtime();

    /* MPI */
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

        double local_elapsed = MPI_Wtime() - start_time;

        int local_completed = t + 1;
        int global_completed = local_completed * size;
        double tasks_per_sec = (double)local_completed / local_elapsed;
        int remaining_tasks = total_tasks - global_completed;
        double remaining_time = remaining_tasks / tasks_per_sec / size;

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
    glbSetOscParams(central_values, asin(sqrt(0.02224)), GLB_THETA_13);
    glbSetOscParams(central_values, asin(sqrt(0.541)), GLB_THETA_23);
    glbSetOscParams(central_values, 1.53 * M_PI, GLB_DELTA_CP);
    glbSetOscParams(central_values, -2.481e-3 + sdm, GLB_DM_31); // DM31 = DM32 + DM21
    glbSetOscillationParameters(central_values);
    glbSetRates();
    set_proj_SM_TH23_DCP_2D(my_projection);
    glbSetProjection(my_projection);

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

        double local_elapsed = MPI_Wtime() - start_time;

        int local_completed = t + 1;
        int global_completed = local_completed * size;
        double tasks_per_sec = (double)local_completed / local_elapsed;
        int remaining_tasks = total_tasks - global_completed;
        double remaining_time = remaining_tasks / tasks_per_sec / size;

        if (rank == 0) {
            printf("Progress: %d/%d tasks completed. average time: %.2f /s, %.2f s left.\n",
                   global_completed, total_tasks, tasks_per_sec, remaining_time);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // double global_chi_min, global_theta23_min, global_deltacp_min;
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
    displs = NULL;
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
    glbFreeProjection(my_projection);

    exit(0);
}
