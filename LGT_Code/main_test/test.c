
/*******************************************************************************
*
* File test.c

* Copyright (C) 2021 Alessandro Sciarra
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Program under construction for testing lattice gauge theory code
*
* Syntax: test -i <input-filename>
*
* For usage instructions see the file README.main
*
*******************************************************************************/

#define MAIN_C

#include "ranlxd.h"
#include "modules.h"
#include "test_utils.h"
#include <assert.h>
#include <stdio.h>

static double generateSU2InHeatbath_Creutz(su2mat *X, double alpha)
{
    double ud1[2], ud2[3], x[4], norm;
    int count = 0;
    double y_min = exp(-alpha);
    double y_max = exp(+alpha);
    double y;

    do
    {
        ranlxd(ud1, 2);
        y = y_min + (y_max - y_min) * ud1[0];

        x[0] = log(y) / alpha;
        count++;

    } while (ud1[1] > sqrt(1.0 - pow(x[0], 2.0)));

    do
    {
        ranlxd(ud2, 3);
        x[1] = 2.0 * ud2[0] - 1.0;
        x[2] = 2.0 * ud2[1] - 1.0;
        x[3] = 2.0 * ud2[2] - 1.0;

        norm = pow(x[1], 2.0) + pow(x[2], 2.0) + pow(x[3], 2.0);

    } while (norm < 0.0000000001 || norm > 1.0);

    norm = sqrt((1.0 - pow(x[0], 2.0)) / norm);

    x[1] *= norm;
    x[2] *= norm;
    x[3] *= norm;

    mk_dble_array_su2(x, *X);
    return count;
}

static double generateSU2InHeatbath_Improved(su2mat *X, double alpha)
{
    int count = 0;
    double ud[4], delta, phi, theta, length, x[4], angles[2];

    double X1, X2, C, A;

    do
    {
        ranlxd(ud, 4);

        X1 = -log(1.0 - ud[0]) / alpha;
        X2 = -log(1.0 - ud[1]) / alpha;
        C = pow(cos(2.0 * M_PI * ud[2]), 2.0);
        A = X1 * C;
        delta = X2 + A;

        // delta = -1.0 / alpha * (log(1. - ud[0]) + pow(cos(2. * PI * (1. - ud[1])), 2) * log(1. - ud[2]));
        count++;
    } while (pow(ud[3], 2.0) > 1.0 - (0.5 * delta));

    x[0] = 1.0 - delta;

    ranlxd(angles, 2);

    phi = 2.0 * PI * angles[0];
    theta = acos(1.0 - (2.0 * angles[1]));

    length = sqrt(1.0 - pow(x[0], 2.0));
    x[1] = length * sin(theta) * cos(phi);
    x[2] = length * sin(theta) * sin(phi);
    x[3] = length * cos(theta);

    mk_dble_array_su2(x, *X);

#if DEBUG == 1
    double det_1, det_2;
    su2mat X_dag, res;
    su2_dag(X_dag, X);
    su2_mat_mul(res, X_dag, X);
    su2_det(det_2, res);
    su2_det(det_1, X);
    printf("%3.20f   %3.20f\n", det_1, det_2);
#endif
    return count;
}

int main(int argc, char *argv[])
{
    double alpha = 7.0;
    double temp;
    su2mat X;
    int myInt;
    unsigned long int N = 100000000;
    unsigned long int M = 200;
    unsigned long int i;
    unsigned long int histo1[M];
    unsigned long int histo2[M];

    rlxd_init(2, 1);

    /*
    for(i = 0; i < N; i++)
        myInt = generateSU2InHeatbath_Improved(&X, alpha);
    */

    for (i = 0; i < M; i++)
    {
        histo1[i] = 0;
        histo2[i] = 0;
    }

    for (i = 0; i < N; i++)
    {
        myInt = generateSU2InHeatbath_Improved(&X, alpha);
        myInt = (X.c0 + 1.0) / 2.0 * M;
        histo1[myInt]++;
        myInt = generateSU2InHeatbath_Creutz(&X, alpha);
        myInt = (X.c0 + 1.0) / 2.0 * M;
        histo2[myInt]++;
    }

    FILE *fout;
    fout = fopen("improved.dat", "w");
    for (i = 0; i < M; i++)
    {
        fprintf(fout, "%+1.10f\t%1.10f\n", 2.0 * (double)i / M - 1.0 + 1.0 / M, (double)histo1[i] * (1.0 * M / N));
    }

    fclose(fout);
    fout = fopen("oldschool.dat", "w");
    for (i = 0; i < M; i++)
    {
        fprintf(fout, "%+1.10f\t%1.10f\n", 2.0 * (double)i / M - 1.0 + 1.0 / M, (double)histo2[i] * (1.0 * M / N));
    }
    fclose(fout);

    return 0;
}
