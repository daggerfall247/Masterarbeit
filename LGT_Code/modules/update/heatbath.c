
/*******************************************************************************
 *
 * File heatbath.c
 *
 * Copyright (C) 2022 Marco Stilger
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Includes the routines to perform heatbath updates for pure gauge theory in SU(2).
 *
 * Externally accessible functions:
 *
 * void staples(int n,int dir,sun_mat *stap)
 *      Computes the staples for the link starting at point n in direction
 *      dir. The staples matrix is returned via the pointer stap.
 *
 * double localHeatbathUpdate(int n,int dir,int m)
 *        Performs m local Heatbath update steps for the link U_dir(n).
 *        On return it hands back the fraction of successful updates.
 *        For now only works for SUN=2.
 *
 *******************************************************************************/

#define HEATHBATH_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ranlxd.h"
#include "headers.h"
#include "modules.h"

/*
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
*/

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
        C = pow(cos(2.0 * PI * ud[2]), 2.0);
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

#if (SUN == 2)
double localHeatbathUpdate(uint64_t n, int dir, int m)
{
    int count = 0;
    double sqrt_det, alpha;
    su2mat A, X, U;

    staples(n, dir, &A);
    su2_det_sqrt(sqrt_det, A);

    if (sqrt_det <= __DBL_EPSILON__)
    {
        su2RandomMatrix(&U);
        *pu[n][dir] = U;
        return 0.0;
    }

    su2_dble_div(A, sqrt_det);

    alpha = sqrt_det * runParams.beta;

    count += generateSU2InHeatbath_Improved(&X, alpha);
    su2_mat_mul_dag(U, X, A); /*U = X*V^dag*/
    *pu[n][dir] = U;

    return 1.0 / count;
}
#elif (SUN == 3)
double localHeatbathUpdate(uint64_t n, int dir, int m)
{
    su3mat heatbath_su3, old_link, staple_sum, Temp_1, Temp_2;
    su2mat heatbath_su2, generated, new_link_su2;

    staples(n, dir, &staple_sum);
    old_link = *pu[n][dir];

    for (int i = 1; i <= 3; i++)
    {
        su3_mat_mul(heatbath_su3, old_link, staple_sum);
        switch (i)
        {
        case 1:
            extr_sub1(heatbath_su2, heatbath_su3);
            break;
        case 2:
            extr_sub2(heatbath_su2, heatbath_su3);
            break;
        case 3:
            extr_sub3(heatbath_su2, heatbath_su3);
            break;
        }
        su2_det_sqrt(sqrt_det, heatbath_su2);

        if (sqrt_det <= __DBL_EPSILON__)
        {
            /*su3RandomMatrix(pu[n][dir]);*/
            return 0.0;
        }
        su2_dble_div(heatbath_su2, sqrt_det);

        alpha = sqrt_det * runParams.beta * 2.0 / 3.0;

        count += generateSU2InHeatbath(&generated, alpha);
        su2_mat_mul_dag(new_link_su2, generated, heatbath_su2);

        switch (i)
        {
        case 1:
            create_su3_sub1(Temp_1, new_link_su2);
            break;
        case 2:
            create_su3_sub2(Temp_1, new_link_su2);
            break;
        case 3:
            create_su3_sub3(Temp_1, new_link_su2);
            break;
        }

        su3_mat_mul(Temp_2, Temp_1, heatbath_su3);
        heatbath_su3 = Temp_2;
        su3_mat_mul(Temp_2, Temp_1, old_link);
        old_link = Temp_2;
    }

    *pu[n][dir] = old_link;

    return 3.0 / count;
}
#endif
