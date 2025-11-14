/* ################################################################################## *
 * 	This program fits sin2th23 and deltaCP to NOvA far detector data alone	 *
 * ################################################################################## *
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define atoa(x) #x

#include "NOvA_setup.h"

#include <globes/globes.h> // GLoBES library

static int global_counter = 0;
static int nexp = 1, nexpmu = 1;

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
        centralvalue = square(sin(2 * th13Prior)); // 0.0857;
        inputerror = 0.0046;
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

/* Set up an all-free projection (except solar parameters) */
void set_proj_SM(glb_projection in_proje) {
    set_proj_all_fixed(in_proje);
    /* THETA_12, THETA_13, THETA_23, DELTA_CP,  DM_21,     DM_31 */
    glbDefineProjection(in_proje, GLB_FIXED, GLB_FIXED, GLB_FREE, GLB_FREE, GLB_FIXED, GLB_FREE);
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
    /* Initialize libglobes */
    glbInit(argv[0]);

    /* Initialize experiments */
    InitializeNOvA(&glb_experiment_list[0], &glb_num_of_exps);

    /* Select prior function (defined above) */
    glbRegisterPriorFunction(my_prior, NULL, NULL, NULL);

    /* Define standard oscillation parameters */
    double theta12 = asin(sqrt(0.307));
    double theta13 = asin(sqrt(0.0210));
    double theta23 = asin(sqrt(0.57));
    double deltacp = 0.82 * M_PI;
    double sdm = 7.53e-5;
    double ldm = 2.41e-3 + sdm;

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
    glbSetOscParams(input_errors, 0.27, GLB_DELTA_CP);
    glbSetOscParams(input_errors, 1 / (sqrt(1 - 0.307)) * 1 / (2 * sqrt(0.307)) * 0.0011, GLB_THETA_13);
    glbSetOscParams(input_errors, 0.07e-3, GLB_DM_31);
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

    /* Iterate over all theta23 and deltacp values */
    double dcp, sin2th23, res;
    int Nsteps = 101;
    double bar_aux = 0.0;
    double percentage = 100.0 * 1.0 / Nsteps;
    int aux_percent = 0;

    FILE* fp;

    char str[100];

    double dcp_ini = 0.0;
    double dcp_fin = 2.0;
    double sin2th23_ini = 0.28;
    double sin2th23_fin = 0.72;

    double dcpmin, t23min, d31min, chimin;

    // Start the run
    printf("\n Commencing fit to NOvA data...\n");
    strcpy(str, "../data/NOvA/NOvA_fit_data_normal.dat");
    fp = fopen(str, "w+");

    fflush(stdout);

    for (int icp = 0; icp < Nsteps; icp++) {
        aux_percent = (int)100.0 * icp / (Nsteps - 1);
        percent_bar(&bar_aux, percentage, aux_percent, 1.0);

        dcp = dcp_ini + (dcp_fin - dcp_ini) * icp / (Nsteps - 1);
        glbSetOscParams(test_values, dcp * M_PI, GLB_DELTA_CP);

        for (int ith = 0; ith < Nsteps; ith++) {
            sin2th23 = sin2th23_ini + (sin2th23_fin - sin2th23_ini) * ith / (Nsteps - 1);
            glbSetOscParams(test_values, asin(sqrt(sin2th23)), GLB_THETA_23);
            res = glbChiNP(test_values, minimum, GLB_ALL);

            fprintf(fp, "%.4f\t%.4f\t%.4f", dcp, sin2th23, res);

            // this is to remove the space in the last line.
            if (ith < Nsteps - 1 || icp < Nsteps - 1)
                fprintf(fp, "\n");
        }
    }

    percent_bar(&bar_aux, 100, aux_percent, 1.0);
    printf("\n \t Fitting completed for Normal Ordering.\n");
    fclose(fp);

    set_proj_SM(my_projection);
    glbSetProjection(my_projection);
    chimin = glbChiNP(central_values, minimum, GLB_ALL);
    dcpmin = glbGetOscParams(minimum, GLB_DELTA_CP) / M_PI;
    t23min = square(sin(glbGetOscParams(minimum, GLB_THETA_23)));
    d31min = glbGetOscParams(minimum, GLB_DM_31);
    printf("\t Best-fit point: deltacp = %.4f, sin^2th23 = %.4f, dm32^2 = %.4g \t minimum chi2: %.4f\n", dcpmin, t23min, d31min + sdm, chimin);

    // Flip hierarchy
    glbSetOscParams(central_values, asin(sqrt(0.56)), GLB_THETA_23);
    glbSetOscParams(central_values, 1.52 * M_PI, GLB_DELTA_CP);
    glbSetOscParams(central_values, -2.45e-3 + sdm, GLB_DM_31); // DM31 = DM32 + DM21
    glbSetOscillationParameters(central_values);
    glbSetRates();
    set_proj_SM_TH23_DCP_2D(my_projection);
    glbSetProjection(my_projection);

    glbCopyParams(central_values, test_values);

    strcpy(str, "../data/NOvA/NOvA_fit_data_inverted.dat");
    fp = fopen(str, "w+");

    for (int icp = 0; icp < Nsteps; icp++) {
        aux_percent = (int)100.0 * icp / (Nsteps - 1);
        percent_bar(&bar_aux, percentage, aux_percent, 1.0);

        dcp = dcp_ini + (dcp_fin - dcp_ini) * icp / (Nsteps - 1);
        glbSetOscParams(test_values, dcp * M_PI, GLB_DELTA_CP);

        for (int ith = 0; ith < Nsteps; ith++) {
            sin2th23 = sin2th23_ini + (sin2th23_fin - sin2th23_ini) * ith / (Nsteps - 1);
            glbSetOscParams(test_values, asin(sqrt(sin2th23)), GLB_THETA_23);
            res = glbChiNP(test_values, minimum, GLB_ALL);

            fprintf(fp, "%.4f\t%.4f\t%.4f", dcp, sin2th23, res);

            // this is to remove the space in the last line.
            if (ith < Nsteps - 1 || icp < Nsteps - 1)
                fprintf(fp, "\n");
        }
    }

    percent_bar(&bar_aux, 100, aux_percent, 1.0);
    printf("\n \t Fitting completed for Inverted Ordering.\n");
    fclose(fp);

    set_proj_SM(my_projection);
    glbSetProjection(my_projection);
    chimin = glbChiNP(central_values, minimum, GLB_ALL);
    dcpmin = glbGetOscParams(minimum, GLB_DELTA_CP) / M_PI;
    t23min = square(sin(glbGetOscParams(minimum, GLB_THETA_23)));
    d31min = glbGetOscParams(minimum, GLB_DM_31);
    printf("\t Best-fit point: deltacp = %.4f, sin^2th23 = %.4f, dm32^2 = %.4g \t minimum chi2: %.4f\n", dcpmin, t23min, d31min + sdm, chimin);

    // Destroy parameter vector(s)
    glbFreeParams(central_values);
    glbFreeParams(test_values);
    glbFreeParams(input_errors);
    glbFreeParams(minimum);
    glbFreeProjection(my_projection);

    exit(0);
}
