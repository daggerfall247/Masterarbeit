
/*******************************************************************************
 *
 * File wilson.c
 *
 * Copyright (C) 2019 Francesca Cuteri
 * Copyright (C) 2016 Bastian Brandt
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Includes the routines to calculate Wilson loops.
 *
 *******************************************************************************/

#define WILSON_C

#include "modules.h"
#include <stdio.h>

static double compute_wloop(uint64_t n, int d1, int d2, int l1, int l2)
{
    uint64_t ii, nn;
    double tr;
    sun_mat u[2];

    sun_unit(u[0]);
    nn = n;
    for (ii = 0; ii < l1; ii++) // positive d1 direction
    {
        u[1] = u[0];
        sun_mul(u[0], u[1], *pu[nn][d1]); // iu1
        nn = neib[nn][d1];
    }
    for (ii = 0; ii < l2; ii++) // positive d2 direction
    {
        u[1] = u[0];
        sun_mul(u[0], u[1], *pu[nn][d2]); // iu2
        nn = neib[nn][d2];
    }
    for (ii = 0; ii < l1; ii++) // negative d1 direction
    {
        u[1] = u[0];
        nn = neib[nn][d1 + DIM];
        sun_mul_dag(u[0], u[1], *pu[nn][d1]); // iu1
    }
    for (ii = 0; ii < l2; ii++) // negative d2 direction
    {
        u[1] = u[0];
        nn = neib[nn][d2 + DIM];
        sun_mul_dag(u[0], u[1], *pu[nn][d2]); // iu2
    }

    sun_trace(tr, u[0]);
    return tr;
}

#if (MASTER_FIELD == 1)
void measureWilsonLoop(double *w)
{
    uint64_t it, ir, nr, nt, nn, inn;
    uint64_t in, dir, ii, index;
    uint64_t nloop;
    char cnfg_file_temp[FULL_PATH_SIZE + 12];

    uint64_t volumeFactor_small[DIM], volumeFactor[DIM], volumeFactor_mf[DIM];
    uint64_t numSmallLattices, coord[DIM];
    uint64_t latticeExtent_small[DIM], latticeExtent_mf[DIM];
    uint64_t VOL_small;

    nt = ((measParams.tf - measParams.ts) / measParams.dt) + 1;
    nr = ((measParams.rf - measParams.rs) / measParams.dr) + 1;
    nloop = nt * nr;
    for (nn = 0; nn < nloop; nn++)
        w[nn] = 0.;

#if (DIM > 3)
    error(measParams.rf >= latParams.linearExtent[3], "measureWilsonLoop (mf)", "Sublattice too small in z direction to measure loops");
#endif
#if (DIM > 2)
    error(measParams.rf >= latParams.linearExtent[2], "measureWilsonLoop (mf)", "Sublattice too small in y direction to measure loops");
#endif
    error(measParams.rf >= latParams.linearExtent[1], "measureWilsonLoop (mf)", "Sublattice too small in x direction to measure loops");
    error(measParams.tf >= latParams.linearExtent[0], "measureWilsonLoop (mf)", "Sublattice too small in t direction to measure loops");

#if (DIM == 2)
    latticeExtent_mf[0] = LENGT_MF;
    latticeExtent_mf[1] = LENGS1_MF;
#if (FRACT == 1)
    latticeExtent_small[0] = LENGT;
#else
    latticeExtent_small[0] = LENGT - measParams.tf;
#endif
#if (FRACS1 == 1)
    latticeExtent_small[1] = LENGS1;
#else
    latticeExtent_small[1] = LENGS1 - measParams.rf;
#endif
    volumeFactor[0] = LENGS1;
    volumeFactor[1] = 1;
    volumeFactor_mf[0] = LENGS1_MF;
    volumeFactor_mf[1] = 1;
#elif (DIM == 3)
    latticeExtent_mf[0] = LENGT_MF;
    latticeExtent_mf[1] = LENGS1_MF;
    latticeExtent_mf[2] = LENGS2_MF;
#if (FRACT == 1)
    latticeExtent_small[0] = LENGT;
#else
    latticeExtent_small[0] = LENGT - measParams.tf;
#endif
#if (FRACS1 == 1)
    latticeExtent_small[1] = LENGS1;
#else
    latticeExtent_small[1] = LENGS1 - measParams.rf;
#endif
#if (FRACS2 == 1)
    latticeExtent_small[2] = LENGS2;
#else
    latticeExtent_small[2] = LENGS2 - measParams.rf;
#endif
    volumeFactor[0] = LENGS1 * LENGS2;
    volumeFactor[1] = LENGS2;
    volumeFactor[2] = 1;
    volumeFactor_mf[0] = LENGS1_MF * LENGS2_MF;
    volumeFactor_mf[1] = LENGS2_MF;
    volumeFactor_mf[2] = 1;
#elif (DIM == 4)
    latticeExtent_mf[0] = LENGT_MF;
    latticeExtent_mf[1] = LENGS1_MF;
    latticeExtent_mf[2] = LENGS2_MF;
    latticeExtent_mf[3] = LENGS3_MF;
#if (FRACT == 1)
    latticeExtent_small[0] = LENGT;
#else
    latticeExtent_small[0] = LENGT - measParams.tf;
#endif
#if (FRACS1 == 1)
    latticeExtent_small[1] = LENGS1;
#else
    latticeExtent_small[1] = LENGS1 - measParams.rf;
#endif
#if (FRACS2 == 1)
    latticeExtent_small[2] = LENGS2;
#else
    latticeExtent_small[2] = LENGS2 - measParams.rf;
#endif
#if (FRACS3 == 1)
    latticeExtent_small[3] = LENGS3;
#else
    latticeExtent_small[3] = LENGS3 - measParams.rf;
#endif
    volumeFactor[0] = LENGS1 * LENGS2 * LENGS3;
    volumeFactor[1] = LENGS2 * LENGS3;
    volumeFactor[2] = LENGS3;
    volumeFactor[3] = 1;
    volumeFactor_mf[0] = LENGS1_MF * LENGS2_MF * LENGS3_MF;
    volumeFactor_mf[1] = LENGS2_MF * LENGS3_MF;
    volumeFactor_mf[2] = LENGS3_MF;
    volumeFactor_mf[3] = 1;
#endif

    numSmallLattices = 1;
    for (dir = 0; dir < DIM; dir++)
    {
        while (latticeExtent_mf[dir] % latticeExtent_small[dir] > 0)
        {
            latticeExtent_small[dir] -= 1;
        }
        // logging("LatticeExtent_small[%d] = %d\n", dir, latticeExtent_small[dir]);
        numSmallLattices *= (latticeExtent_mf[dir] / latticeExtent_small[dir]);

        coord[dir] = 0;
    }
    // logging("numSmallLattices = %d\n\n", numSmallLattices);

#if (DIM == 2)
    VOL_small = latticeExtent_small[0] * latticeExtent_small[1];
    volumeFactor_small[0] = latticeExtent_small[1];
    volumeFactor_small[1] = 1;
#elif (DIM == 3)
    VOL_small = latticeExtent_small[0] * latticeExtent_small[1] * latticeExtent_small[2];
    volumeFactor_small[0] = latticeExtent_small[1] * latticeExtent_small[2];
    volumeFactor_small[1] = latticeExtent_small[2];
    volumeFactor_small[2] = 1;
#elif (DIM == 4)
    VOL_small = latticeExtent_small[0] * latticeExtent_small[1] * latticeExtent_small[2] * latticeExtent_small[3];
    volumeFactor_small[0] = latticeExtent_small[1] * latticeExtent_small[2] * latticeExtent_small[3];
    volumeFactor_small[1] = latticeExtent_small[2] * latticeExtent_small[3];
    volumeFactor_small[2] = latticeExtent_small[3];
    volumeFactor_small[3] = 1;
#endif

    uint64_t pos_small[numSmallLattices];
    index = 0;
    for (in = 0; in < numSmallLattices; in++)
    {
        pos_small[in] = index;
        // logging("pos_small[%d] = %d\n", in, pos_small[in]);
        for (dir = DIM - 1; dir >= 0; dir--)
        {
            if (coord[dir] + latticeExtent_small[dir] < latticeExtent_mf[dir])
            {
                coord[dir] += latticeExtent_small[dir];
                break;
            }
            else
            {
                coord[dir] = 0;
            }
        }
        index = 0;
        for (dir = 0; dir < DIM; dir++)
        {
            index += volumeFactor_mf[dir] * coord[dir];
        }
    }

    uint64_t i_small[VOL_small];
    for (ii = 0; ii < VOL_small; ii++) // TODO: maybe test, but seems ok
    {
        index = ii;
        for (dir = 0; dir < DIM; dir++) // calc small coordinates
        {
            coord[dir] = index / volumeFactor_small[dir];
            index = index % volumeFactor_small[dir];
        }
        index = 0;
        for (dir = 0; dir < DIM; dir++) // calc normal index
        {
            index += coord[dir] * volumeFactor[dir];
        }
        i_small[ii] = index;

        // logging("i_small[%d] = %d\n", ii, i_small[ii]);
    }

    sprintf(cnfg_file_temp, "%s_temp", CNFG_FILE);
    for (in = 0; in < numSmallLattices; in++) // TODO: choose in smart!!!!!
    {
        index = pos_small[in];
        updateArrayOfIndexMF(index);
        readConfig(cnfg_file_temp);

        for (ii = 0; ii < VOL_small; ii++)
        {
            nn = 0;
            for (ir = measParams.rs; ir <= measParams.rf; ir += measParams.dr, nn++)
            {
                inn = nn * nt;
                for (it = measParams.ts; it <= measParams.tf; it += measParams.dt, inn++)
                {
                    for (dir = 1; dir < DIM; dir++)
                    {
                        w[inn] += compute_wloop(i[i_small[ii]], 0, dir, it, ir);
                    }
                }
            }
        }
    }

    for (nn = 0; nn < nloop; nn++)
    {
        w[nn] /= (double)((DIM - 1) * VOL_MF * SUN);
    }
}
#elif (MASTER_FIELD == 0)
void measureWilsonLoop(double *w)
{
    uint64_t nr, nt, nn, inn;
    uint64_t in, dir;
    uint64_t lbord, nloop;

    if(measParams.mwil_mode <= 0)
    {
        nt = ((measParams.tf - measParams.ts) / measParams.dt) + 1;
        nr = ((measParams.rf - measParams.rs) / measParams.dr) + 1;
        nloop = nt * nr;
    }
    else
    {
        nt = measParams.mwil_mode;
        nr = measParams.mwil_mode;
        nloop = measParams.mwil_mode;
    }
    
    for (nn = 0; nn < nloop; nn++)
        w[nn] = 0.;

    lbord = latParams.linearExtent[0] - 1;
    error(lbord < measParams.tf, "meas_wilson", "Lattice too small to measure loops!");

    //for(inn = 0; inn < nloop; inn++)
    //{
    //    printf("r,t = %d, %d\n", measParams.rExtents[inn], measParams.tExtents[inn] );
    //}
    
    for (in = 0; in < VOL; in++)
    {
        for (dir = 1; dir < DIM; dir++)
        {
            for(inn = 0; inn < nloop; inn++)
            {
                w[inn] += compute_wloop(in, 0, dir, measParams.tExtents[inn], measParams.rExtents[inn]);
            }
            /*
            for (ir = 0; ir < nr; ir++, nn++)
            {
                inn = nn * nt;
                for (it = 0; it < nt; it++, inn++)
                {
                    w[inn] += compute_wloop(in, 0, dir, measParams.tExtents[ir*nt + it], measParams.rExtents[ir*nt + it]);
                }
            }
            */
        }
    }
    

    for (nn = 0; nn < nloop; nn++)
    {
        w[nn] /= (double)((DIM - 1) * VOL * SUN);
    }
}
#endif
