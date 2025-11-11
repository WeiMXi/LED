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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char MYFILE[] = "../data/T2K/T2K_fit.dat";

double SQR(double x) {
    return x * x;
}

int main(int argc, char* argv[]) {
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

    /* Iterate over all deltacp values and find the global minimum */
    double thetheta23, theta23_min, res;
    double thedeltacp, deltacp_min;
    double x, y;
    double chi_min = 1000000.0;
    int z = 0;
    double xmin = 0.35;
    double xmax = 0.65;
    int xsteps = 100;
    double ymin = 0;
    double ymax = 2 * M_PI;
    int ysteps = 100;
    InitOutput(MYFILE, "");

    for (x = xmin; x <= xmax; x = x + (xmax - xmin) / xsteps) {
        for (y = ymin; y <= ymax; y = y + (ymax - ymin) / ysteps) {
            /* Set vector of test values */
            thetheta23 = asin(sqrt(x)); // Sin2 theta23 to radian
            thedeltacp = y;
            glbSetOscParams(test_values, thetheta23, GLB_THETA_23);
            glbSetOscParams(test_values, thedeltacp, GLB_DELTA_CP);
            /* Compute Chi^2 for all loaded experiments and all rules */

            res = glbChiNP(test_values, minimum, GLB_ALL);
            AddToOutput(x, y, res);
            if (res < chi_min) {
                theta23_min = thetheta23;
                chi_min = res;
                deltacp_min = thedeltacp;
            }
            z++;
            if (z % 100 == 0) {
                printf("%d\n", z);
            }
        }
    }
    printf("\n\nBest-fit value: Sin^2 theta23 = %f degrees, DeltaCP = %f and chi^2_min = %f \n\n", SQR(sin(theta23_min)), deltacp_min, chi_min);

    glbFreeParams(central_values);

    /* Clear experiment list */
    glbClearExperimentList();

    exit(0);
}
