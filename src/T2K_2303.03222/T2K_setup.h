#ifndef T2Ksetup_h
#define T2Ksetup_h

#include "sys_T2K.h"

#include <glb_smear.h>
#include <glb_types.h>
#include <globes/globes.h> // GLoBES library
#include <gsl/gsl_sf_erf.h>
#include <math.h>
#include <stdlib.h>

double T2K_sigma(int j, double par1, double par2, const glb_smear* data) {
    double energy = data->simbincenter[j];

    return par1 * sqrt(energy) + par2 * pow(energy, 2);
}

static double T2K_SmearMatrixElementA(int i, int j, double e_shift, double par1, double par2, const glb_smear* data) {
    // double e_shift = -0.2; /////// This is the energy shift!!
    double p1, p2;
    p1 = 0.5 * (1 + gsl_sf_erf(
                        (glb_upper_bin_boundary(i, data) - glb_sbin_center(j, data) - e_shift) /
                        (sqrt(2) * T2K_sigma(j, par1, par2, data))));
    p2 = 0.5 * (1 + gsl_sf_erf(
                        (glb_lower_bin_boundary(i, data) - glb_sbin_center(j, data) - e_shift) /
                        (sqrt(2) * T2K_sigma(j, par1, par2, data))));

    return (p1 - p2);
}

void create_T2K_smear(int exp, int channel, double e_shift, double par1, double par2) // static double**
{

    int i, j;
    double aux;

    struct glb_experiment* e = glb_experiment_list[exp];

    int rows = e->smear_data[channel]->numofbins + 1;
    int columns = e->smear_data[channel]->simbins + 1;

    for (int i = 0; i < e->smear_data[channel]->numofbins; i++) {
        free(e->smear[channel][i]);
    }
    free(e->lowrange[channel]);
    free(e->uprange[channel]);

    e->lowrange[channel] = (int*)malloc(rows * sizeof(int));
    e->uprange[channel] = (int*)malloc(rows * sizeof(int));
    for (int i = 0; i < e->smear_data[channel]->numofbins; i++) {

        e->smear[channel][i] = (double*)malloc(columns * sizeof(double));
    }

    for (i = 0; i < e->smear_data[channel]->numofbins; i++) {
        e->lowrange[channel][i] = 0;
        e->uprange[channel][i] = e->smear_data[channel]->simbins - 1;

        for (j = 0; j < e->smear_data[channel]->simbins; j++) {
            aux = T2K_SmearMatrixElementA(i, j, e_shift, par1, par2, e->smear_data[channel]);
            if (aux > 1.0e-6)
                e->smear[channel][i][j] = aux;
            else
                e->smear[channel][i][j] = 0;
        }
    }

    e->lowrange[channel][i] = -1;
    e->uprange[channel][i] = -1;
    e->smear[channel][i] = NULL;
}

void InitializeT2K(glb_exp* in, int* counter) {
    double rule_eshift, rule_par1, rule_par2;
    // Initialize experiment list
    glbInitExperiment("../experiments/T2K_2303.03222/T2K_numu_2023.glb", in, counter);

    // This will create a smearing metrix of globes type "A", on exp, channel,  with energy shifted e_shift and parameters par1 and par2
    // Check T2K_SmearMatrixElementA, e_shift variable to see the energy shift.
    // the variables are: create_T2K_smear(exp, channel, e_shift, par1, par2)
    int Nexp = *counter - 1;
    int Nsmear = 0;

    create_T2K_smear(Nexp, Nsmear, 0.035, 0.16, 0.0); // 0.035, 0.15

    Nsmear = Nsmear + 1;
    create_T2K_smear(Nexp, Nsmear, -0.05, 0.11, 0.07);

    Nsmear = Nsmear + 1;
    create_T2K_smear(Nexp, Nsmear, -0.3, 0.0001, 0.0001);

    glbInitExperiment("../experiments/T2K_2303.03222/T2K_anumu_2023.glb", in, counter);

    Nexp = *counter - 1;
    Nsmear = 0;
    rule_eshift = 0;
    rule_par1 = 0;
    rule_par2 = 0;
    create_T2K_smear(Nexp, Nsmear, 0.03 + rule_eshift, 0.1 + rule_par1, 0.05 + rule_par2); // 0.035, 0.15

    Nsmear = Nsmear + 1;
    create_T2K_smear(Nexp, Nsmear, -0.05 + rule_eshift, 0.12 + rule_par1, 0.11 + rule_par2);

    Nsmear = Nsmear + 1;
    create_T2K_smear(Nexp, Nsmear, -0.3 + rule_eshift, 0.0001 + rule_par1, 0.0001 + rule_par2);

    glbInitExperiment("../experiments/T2K_2303.03222/T2K_nue_2023.glb", in, counter);

    Nexp = *counter - 1;
    Nsmear = 0;

    create_T2K_smear(Nexp, Nsmear, -0.05, 0.13, 0.1); // nueCCQE

    Nsmear = Nsmear + 1;
    create_T2K_smear(Nexp, Nsmear, -0.25, 0.12, 0.11); // nue_NC

    Nsmear = Nsmear + 1;

    create_T2K_smear(Nexp, Nsmear, -0.025, 0.14, 0.0); // anueCCQE

    Nsmear = Nsmear + 1;

    create_T2K_smear(Nexp, Nsmear, -0.05, 0.13, 0.1); // anueCCQEWrong

    Nsmear = Nsmear + 1;
    create_T2K_smear(Nexp, Nsmear, -0.35, 0.12, 0.11); // anue_NC

    Nsmear = Nsmear + 1;
    create_T2K_smear(Nexp, Nsmear, -0.01, 0.04, 0.1); // nueCCQE // nueCCQE -0.01

    Nsmear = Nsmear + 1;
    create_T2K_smear(Nexp, Nsmear, -0.25, 0.12, 0.11); // nue_NC
}
#endif
