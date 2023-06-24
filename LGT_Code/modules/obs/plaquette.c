
/*******************************************************************************
 *
 * File plaquette.c
 *
 * Copyright (C) 2013 Bastian Brandt
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Includes the routines to calculate the plaquette and the staples.
 *
 * double plaquette(void)
 *      Returns the average value of the plaquette on the gauge field
 *      configuration.
 *
 * double gaugeAction(void)
 *      Returns the Wilson gauge action measured on the gauge field
 *      configuration.
 *******************************************************************************/

#define PLAQUETTE_C

#include <stdlib.h>
#include <math.h>
#include "gauge.h"
#include "headers.h"
#include "modules.h"

#if (MASTER_FIELD == 1)
double plaquette(void)
{
    uint64_t ii, jj, kk, n;
    double plaq, tr;
    sun_mat *un;

    checkpoint("plaquette");

    plaq = 0.0;
    un = malloc(3 * sizeof(sun_mat));
    for (ii = 0; ii < VOL; ii++)
    {
        for (jj = 0; jj < (DIM - 1); jj++)
        {
            for (kk = jj + 1; kk < DIM; kk++)
            {
                un[0] = *pu[i[ii]][jj];
                n = neib[i[ii]][jj];
                un[1] = *pu[n][kk];
                sun_mul(un[2], un[0], un[1]);
                n = neib[i[ii]][kk];
                sun_dag(un[0], *pu[n][jj]);
                sun_mul(un[1], un[2], un[0]);
                sun_dag(un[2], *pu[i[ii]][kk]);
                sun_mul(un[0], un[1], un[2]);
                sun_trace(tr, un[0]);
                plaq += tr;
            }
        }
    }
    plaq /= (double)(SUN * NPLAQ * VOL);
    free(un);

    return plaq;
}
#elif (MASTER_FIELD == 0)
double plaquette(void)
{
    uint64_t ii, jj, kk, n;
    double plaq, tr;
    sun_mat *un;

    checkpoint("plaquette");

    plaq = 0.;
    un = malloc(3 * sizeof(sun_mat));
    for (ii = 0; ii < VOL; ii++)
    {
        for (jj = 0; jj < (DIM - 1); jj++)
        {
            for (kk = jj + 1; kk < DIM; kk++)
            {
                un[0] = *pu[ii][jj];
                n = neib[ii][jj];
                un[1] = *pu[n][kk];
                sun_mul(un[2], un[0], un[1]);
                n = neib[ii][kk];
                sun_dag(un[0], *pu[n][jj]);
                sun_mul(un[1], un[2], un[0]);
                sun_dag(un[2], *pu[ii][kk]);
                sun_mul(un[0], un[1], un[2]);
                sun_trace(tr, un[0]);
                plaq += tr;
            }
        }
    }
    plaq /= (double)(SUN * NPLAQ * VOL);
    free(un);

    return plaq;
}
#endif

double gaugeAction(void)
{
    double plaq;

    plaq = plaquette();
    return runParams.beta * NPLAQ * VOL * (1. - plaq);
}
