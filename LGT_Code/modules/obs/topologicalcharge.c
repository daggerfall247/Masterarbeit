//
//  topologicalcharge.c
//  glueballs
//
//  Created by Carolin Riehl on 14.05.20.
//  Copyright © 2020 Carolin Riehl. All rights reserved.
//

#define TOPCHARGE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gauge.h"
#include "headers.h"
#include "modules.h"

static void computeLine(sun_mat *u, uint64_t *n, int dir, int l, int flag)
{
    sun_mat temp1, temp2;
    int i;

    sun_unit(temp2);
    temp1 = *u;
    if (flag == 0)
    {
        for (i = 0; i < l; i++) // pos dir
        {
            temp2 = temp1;
            sun_mul(temp1, temp2, *pu[*n][dir]);
            *n = neib[*n][dir];
        }
    }
    else if (flag == 1)
    {
        for (i = 0; i < l; i++) // neg dir
        {
            *n = neib[*n][dir + DIM];
            temp2 = temp1;
            sun_mul_dag(temp1, temp2, *pu[*n][dir]);
        }
    }
    *u = temp1;
}

static void computeCloverAverages(sun_mat *u, uint64_t n, int dir1, int dir2, int e1, int e2)
{
    uint64_t n1;
    int i, l1, l2;
    sun_mat temp1, sum;

    sun_zero(sum);

    n1 = n;

    l1 = e1;
    l2 = e2;
    for (i = 0; i < 2; i++)
    {
        if (e1 != e2 && i == 1)
        {
            l1 = e2;
            l2 = e1;
        }
        else if (e1 == e2 && i == 1)
        {
            break;
        }
        //+dir1,+dir2
        sun_unit(temp1);
        computeLine(&temp1, &n1, dir1, l1, 0);
        computeLine(&temp1, &n1, dir2, l2, 0);
        computeLine(&temp1, &n1, dir1, l1, 1);
        computeLine(&temp1, &n1, dir2, l2, 1);
        sun_self_add(sum, temp1);

        //+dir1,-dir2
        sun_unit(temp1);
        computeLine(&temp1, &n1, dir2, l2, 1);
        computeLine(&temp1, &n1, dir1, l1, 0);
        computeLine(&temp1, &n1, dir2, l2, 0);
        computeLine(&temp1, &n1, dir1, l1, 1);
        sun_self_add(sum, temp1);

        //-dir1,+dir2
        sun_unit(temp1);
        computeLine(&temp1, &n1, dir2, l2, 0);
        computeLine(&temp1, &n1, dir1, l1, 1);
        computeLine(&temp1, &n1, dir2, l2, 1);
        computeLine(&temp1, &n1, dir1, l1, 0);
        sun_self_add(sum, temp1);

        //-dir1,-dir2
        sun_unit(temp1);
        computeLine(&temp1, &n1, dir1, l1, 1);
        computeLine(&temp1, &n1, dir2, l2, 1);
        computeLine(&temp1, &n1, dir1, l1, 0);
        computeLine(&temp1, &n1, dir2, l2, 0);
        sun_self_add(sum, temp1);
    }
    if (e1 != e2)
    {
        sun_dble_div(sum, 8.);
    }
    else
    {
        sun_dble_div(sum, 4.);
    }
    *u = sum;
}

static void hermitianTraceless(sun_mat *U)
{
    sun_mat X, X_dag;

    X = *U;
    sun_dag(X_dag, X);
    sun_self_sub(X, X_dag);
    sun_dble_div(X, 2.);
    *U = X;
}

static double computeCloverProducts(uint64_t x, int d1, int d2, int d3, int d4)
{
    double tr;
    sun_mat U1, U2, U3, F_munu, F_rhosigma;
    double k1 = 1.5, k2 = -0.15, k3 = 1. / 90.;

    
    
    computeCloverAverages(&U1, x, d1, d2, 1, 1);
    computeCloverAverages(&U2, x, d1, d2, 2, 2);
    computeCloverAverages(&U3, x, d1, d2, 3, 3);

    sun_dble_mul(U1, k1);
    sun_dble_mul(U2, k2);
    sun_dble_mul(U3, k3);

    sun_add(F_munu, U1, U2);
    sun_self_add(F_munu, U3);
    
    hermitianTraceless(&F_munu);

    
    computeCloverAverages(&U1, x, d3, d4, 1, 1);
    computeCloverAverages(&U2, x, d3, d4, 2, 2);
    computeCloverAverages(&U3, x, d3, d4, 3, 3);

    sun_dble_mul(U1, k1);
    sun_dble_mul(U2, k2);
    sun_dble_mul(U3, k3);

    sun_add(F_rhosigma, U1, U2);
    sun_self_add(F_rhosigma, U3);

    hermitianTraceless(&F_rhosigma);

    sun_mul(U1, F_munu, F_rhosigma);
    sun_trace(tr, U1);

    return tr;
}

double topologicalCharge(void)
{
    uint64_t n;
    double tstart1, tend1;
    double qtop;
    checkpoint("meas_topologicalcharge");

    tstart1 = getTime();
    qtop = 0.0;
    logging("calculating topological charge ..\n");
#if (MASTER_FIELD == 0)
    for (n = 0; n < VOL; n++)
    {
        qtop += 8. * computeCloverProducts(n, 0, 1, 2, 3) + 8. * computeCloverProducts(n, 0, 2, 3, 1) + 8. * computeCloverProducts(n, 0, 3, 1, 2);
    }
#elif (MASTER_FIELD == 1)
    uint64_t k;
    char cnfg_file_temp[FULL_PATH_SIZE + 12];
    sprintf(cnfg_file_temp, "%s_temp", CNFG_FILE);
    for (k = 0; k < VOL_SL; k++)
    {
        updateArrayOfIndexMF(sl[k]);
        readConfig(cnfg_file_temp);
        for (n = 0; n < VOL; n++)
        {
            qtop += 8. * computeCloverProducts(n, 0, 1, 2, 3) + 8. * computeCloverProducts(n, 0, 2, 3, 1) + 8. * computeCloverProducts(n, 0, 3, 1, 2);
        }
    }

#endif
    qtop /= (32. * PI * PI);
    tend1 = getTime();
    logging("measurement done (took %.3e sec)\n", tend1 - tstart1);

    return qtop;
}

static void U_plpl(sun_mat *u, int ii, int jj, int kk)
{
    int n;
    sun_mat *un;

    un = malloc(3 * sizeof(sun_mat));

    // Compute Plaquette U_mu nu.
    un[0] = *pu[ii][jj];
    n = neib[ii][jj];
    un[1] = *pu[n][kk];
    sun_mul(un[2], un[0], un[1]);
    n = neib[ii][kk];
    sun_dag(un[0], *pu[n][jj]);
    sun_mul(un[1], un[2], un[0]);
    sun_dag(un[2], *pu[ii][kk]);
    sun_mul(un[0], un[1], un[2]);

    u[0] = un[0];
    free(un);
}

static void U_mipl(sun_mat *u, int ii, int jj, int kk)
{
    int n;
    sun_mat *un;

    un = malloc(3 * sizeof(sun_mat));

    // Compute Plaquette U_mu nu.
    n = neib[ii][jj + DIM];
    sun_dag(un[0], *pu[n][jj]);
    un[1] = *pu[n][kk];
    sun_mul(un[2], un[0], un[1]);

    n = neib[n][kk];

    sun_mul(un[1], un[2], *pu[n][jj]);
    sun_dag(un[2], *pu[ii][kk]);
    sun_mul(un[0], un[1], un[2]);

    u[0] = un[0];
    free(un);
}

static void U_mimi(sun_mat *u, int ii, int jj, int kk)
{
    int n;
    sun_mat *un;

    un = malloc(3 * sizeof(sun_mat));

    // Compute Plaquette U_ -mu -nu.
    n = neib[ii][jj + DIM];
    sun_dag(un[0], *pu[n][jj]);
    n = neib[n][kk + DIM];
    sun_mul_dag(un[2], un[0], *pu[n][kk]);
    sun_mul(un[1], un[2], *pu[n][jj]);
    n = neib[ii][kk + DIM];
    sun_mul(un[0], un[1], *pu[n][kk]);

    u[0] = un[0];
    free(un);
}

static void U_plmi(sun_mat *u, int ii, int jj, int kk)
{
    int n;
    sun_mat *un;

    un = malloc(3 * sizeof(sun_mat));

    // Compute Plaquette U_ +mu -nu.
    un[0] = *pu[ii][jj];
    n = neib[ii][jj];
    n = neib[n][kk + DIM];
    sun_mul_dag(un[2], un[0], *pu[n][kk]);
    n = neib[ii][kk + DIM];
    sun_dag(un[0], *pu[n][jj]);
    sun_mul(un[1], un[2], un[0]);
    sun_mul(un[0], un[1], *pu[n][kk]);

    u[0] = un[0];
    free(un);
}

void clover(sun_mat *clov, int x, int d1, int d2)
{
    sun_mat *un, *p;

    // checkpoint("compute_clover");

    un = malloc(3 * sizeof(sun_mat));
    p = malloc(2 * sizeof(sun_mat));
    sun_zero(p[0]);
    sun_zero(p[1]);

    U_plpl(un, x, d1, d2);
    sun_add(p[1], un[0], p[0]);

    U_mipl(un, x, d2, d1);
    sun_add(p[0], un[0], p[1]);
    U_mimi(un, x, d1, d2);
    sun_add(p[1], un[0], p[0]);
    U_plmi(un, x, d2, d1);
    sun_add(p[0], un[0], p[1]);

    un[0] = p[0];

    // Im(C_{mu nu}).
    sun_dag(un[1], un[0]);
    sun_sub(p[0], un[0], un[1]);

    sun_dble_div(p[0], 8.);
    *clov = p[0];
}

double compute_clover_products(int x, int d1, int d2, int d3, int d4)
{

    double tr;
    sun_mat *un, *unn, *p;

    // checkpoint("compute_plaquette");

    un = malloc(3 * sizeof(sun_mat));
    unn = malloc(3 * sizeof(sun_mat));
    p = malloc(2 * sizeof(sun_mat));
    sun_zero(p[0]);
    sun_zero(p[1]);

    U_plpl(un, x, d1, d2);

    sun_add(p[1], un[0], p[0]);

    U_mipl(un, x, d2, d1);
    sun_add(p[0], un[0], p[1]);
    U_mimi(un, x, d1, d2);
    sun_add(p[1], un[0], p[0]);
    U_plmi(un, x, d2, d1);
    sun_add(p[0], un[0], p[1]);

    un[0] = p[0];
    sun_dble_div(un[0], 4.);

    sun_zero(p[0]);
    sun_zero(p[1]);

    U_plpl(unn, x, d3, d4);

    sun_add(p[1], unn[0], p[0]);
    U_mipl(unn, x, d4, d3);
    sun_add(p[0], unn[0], p[1]);
    U_mimi(unn, x, d3, d4);
    sun_add(p[1], unn[0], p[0]);
    U_plmi(unn, x, d4, d3);
    sun_add(p[0], unn[0], p[1]);
    unn[0] = p[0];
    sun_dble_div(unn[0], 4.);

    // compute trace of product of imaginary part of plaquettes C_{mu nu} C_{rho sigma}.

    // compute Im (C) * Im (C).
    sun_zero(p[0]);
    sun_zero(p[1]);
    sun_zero(un[1]);
    sun_zero(unn[1]);
    sun_zero(un[2]);

    // Im(C_{mu nu}).
    sun_dag(un[1], un[0]);
    sun_sub(p[0], un[0], un[1]);

    // Im(C_{rho sigma}).
    sun_dag(unn[1], unn[0]);
    sun_sub(p[1], unn[0], unn[1]);

    sun_dble_div(p[0], 2.);
    sun_dble_div(p[1], 2.);
    // ImC*ImC.
    sun_mul(un[2], p[0], p[1]);
    sun_trace(tr, un[2]);

    free(un);
    free(unn);
    free(p);

    return tr;
}

void meas_topologicalcharge(double *q)
{
    int n = 0;

    checkpoint("meas_topologicalcharge");

    double q2 = 0;

    //********************* Reduced formula (faster version)  *******************************************

    for (n = 0; n < VOL; n++)
    {
        q2 += -compute_clover_products(n, 0, 1, 3, 2) - compute_clover_products(n, 1, 2, 3, 0) - compute_clover_products(n, 2, 0, 3, 1);
    }

    // printf("Time for computation of topological charge = %e\n", tend1-tstart1);

    *q = q2 / (4. * PI * PI);
}
