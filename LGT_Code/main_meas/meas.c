
/*******************************************************************************
 *
 * File qcd.c
 *
 * Copyright (C) 2019 Francesca Cuteri
 * Copyright (C) 2016 Bastian Brandt
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Program for
 *
 * Syntax: meas -i [input-filename]
 *
 * For usage instructions see the file README.main
 *
 *******************************************************************************/

#define MAIN_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ranlxd.h"
#include "modules.h"

// static int init = 0;

void setGaugeFieldToCold(void)
{
    int n, dir;
    for (n = 0; n < VOL; n++)
        for (dir = 0; dir < DIM; dir++)
            sun_unit(*pu[n][dir]);
}

int main(int argc, char *argv[])
{
    int n, ir, it, nr, nt, nloop;
    int seed;
    int sconf, fconf, econf;
    char out_dir[NAME_SIZE];
    char cnfg_file[FULL_PATH_SIZE + 12];
    double *wil = NULL;
    double qtop;
    double t1, t2;

    readInputFileForMeas(&seed, &sconf, &fconf, &econf, out_dir, argc, argv);
    rlxd_init(2, seed);

    setupOutputFilesForMeas(runParams.idForOutputFilesName, out_dir);
    initArrayOfNeighbours();
    initGaugeField(1);

    if (measParams.mwil_mode <= 0)
    {
        nt = ((measParams.tf - measParams.ts) / measParams.dt) + 1;
        nr = ((measParams.rf - measParams.rs) / measParams.dr) + 1;
        nloop = nt * nr;
    }
    else
    {
        nloop = measParams.mwil_mode;
    }
    if (runParams.mwil)
    {
        wil = malloc(nloop * sizeof(double));
    }
    if (runParams.mcorrs)
        allocateFermionFieldsFor2ptFunctions(measParams.stype);

    /*
    if (init == 0)
    {
       allocateGaugeField(iu1);
       allocateGaugeField(iu2);
       init = 1;
    }
    */

    for (n = sconf; n <= fconf; n += econf)
    {

        sprintf(cnfg_file, "%s_n%d", CNFG_FILE, n);
        // sprintf(cnfg_file, "../../SU2_Marc/lattice_SU2/configs/16x16x16x16/conf.%4d", n);

        logging("\nMeasuring on configuration:\n %s\n", cnfg_file);
        t1 = getTime();

        // readConfig_MarcStyle(cnfg_file);

        if (runParams.mtopcharge)
        {
            readConfig(cnfg_file);
            smearing_APE_all(measParams.ism_tc, measParams.smpar_tc, pu);
            qtop = topologicalCharge();
            logging("QTOP %.8e\n", qtop);
        }

        if (runParams.mwil)
        {
            // copyGaugeField(pu, iu1);
            // copyGaugeField(pu, iu2);
            readConfig(cnfg_file);
            smearing_APE_temporal(measParams.ismt, measParams.smpart, pu); // iu1
            smearing_APE_spatial(measParams.isms, measParams.smpars, pu);  // iu2

            measureWilsonLoop(wil);
            int inn;
            for (inn = 0; inn < nloop; inn++)
            {
                logging("WIL %d %d\t%.8e\n", measParams.rExtents[inn], measParams.tExtents[inn], wil[inn]);
            }
            /*
            for (ir = 0; ir < nr; ir++, nn++)
            {
                inn = nn * nt;
                for (it = 0; it < nt; it++, inn++)
                {
                    logging("WIL %d %d\t%.8e\n", measParams.rExtents[nt * ir + it], measParams.tExtents[nt * ir + it], wil[inn]);
                    // logging("WIL %d %d\t%.8e\n", measParams.rs + ir * measParams.dr, measParams.ts + it * measParams.dt, wil[inn]);
                }
            }
            */
        }
        if (runParams.mcorrs)
        {
            measure2ptFunctions(measParams.source, measParams.stype, measParams.n2pt, measParams.g1, measParams.g2);
        }
        t2 = getTime();
        logging("\nMeasurement complete (took %.3e sec)\n\n", t2 - t1);
    }

    /*
    if (init == 1)
    {
       deallocateGaugeField(iu1);
       deallocateGaugeField(iu2);
       init = 0;
    }
    */

    if (runParams.mwil)
    {
        free(wil);
        free(measParams.rExtents);
        free(measParams.tExtents);
    }
    if (runParams.mcorrs)
    {
        deallocateFermionFieldsFor2ptFunctions();
        free(measParams.g1);
        free(measParams.g2);
    }

    releaseGaugeField();

    return 0;
}
