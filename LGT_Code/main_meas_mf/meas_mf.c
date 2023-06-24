
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

#if (MASTER_FIELD == 1)
int main(int argc, char *argv[])
{
    int n, ir, it, nr, nt, il;
    uint64_t k;
    int seed;
    int sconf, fconf, econf;
    char out_dir[NAME_SIZE];
    char cnfg_file[FULL_PATH_SIZE + 12], cnfg_file_temp[FULL_PATH_SIZE + 12], cnfg_file_temp2[FULL_PATH_SIZE + 12];
    double *wil = NULL;
    double qtop;
    double t1, t2;

    readInputFileForMeas(&seed, &sconf, &fconf, &econf, out_dir, argc, argv);
    rlxd_init(1, seed);

    setupOutputFilesForMeas(runParams.idForOutputFilesName, out_dir);
    initGaugeField(1);
    initGlobalArrays();

    nt = ((measParams.tf - measParams.ts) / measParams.dt) + 1;
    nr = ((measParams.rf - measParams.rs) / measParams.dr) + 1;
    if (runParams.mwil)
        wil = malloc(nr * nt * sizeof(double));
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

    sprintf(cnfg_file_temp, "%s_temp", CNFG_FILE);
    sprintf(cnfg_file_temp2, "%s_temp2", CNFG_FILE);
    prepareConfig(cnfg_file_temp);
    prepareConfig(cnfg_file_temp2);

    for (n = sconf; n <= fconf; n += econf)
    {
        sprintf(cnfg_file, "%s_n%d", CNFG_FILE, n);
        // sprintf(cnfg_file, "../../SU2_Marc/lattice_SU2/configs/16x16x16x16/conf.%4d", n);

        logging("\nMeasuring on configuration:\n %s\n", cnfg_file);
        t1 = getTime();

        if (runParams.mtopcharge)
        {
            for (k = 0; k < VOL_SL; k++)
            {
                updateArrayOfIndexMF(sl[k]);
                // readConfig_MarcStyle(cnfg_file);
                readConfig(cnfg_file);
                writeConfig(cnfg_file_temp);
            }
            for (it = 0; it < measParams.ism_tc; it++)
            {
                for (k = 0; k < VOL_SL; k++)
                {
                    updateArrayOfIndexMF(sl[k]);
                    readConfig(cnfg_file_temp);
                    smearing_APE_all(1, measParams.smpar_tc, pu); // iu1
                    writeConfig(cnfg_file_temp2);
                }
                for (k = 0; k < VOL_SL; k++)
                {
                    updateArrayOfIndexMF(sl[k]);
                    readConfig(cnfg_file_temp2);
                    writeConfig(cnfg_file_temp);
                }
            }

            measureTopologicalCharge(&qtop);
            logging("QTOP %.8e\n", qtop);
        }

        if (runParams.mwil)
        {
            for (k = 0; k < VOL_SL; k++)
            {
                updateArrayOfIndexMF(sl[k]);
                // readConfig_MarcStyle(cnfg_file);
                readConfig(cnfg_file);
                writeConfig(cnfg_file_temp);
            }

            logging("APE-Smearing on configuration:\n %s\n", cnfg_file);
            t1 = getTime();

            for (it = 0; it < measParams.ismt; it++)
            {
                for (k = 0; k < VOL_SL; k++)
                {
                    updateArrayOfIndexMF(sl[k]);
                    readConfig(cnfg_file_temp);
                    smearing_APE_temporal(1, measParams.smpart, pu); // iu1
                    writeConfig(cnfg_file_temp2);
                }
                for (k = 0; k < VOL_SL; k++)
                {
                    updateArrayOfIndexMF(sl[k]);
                    readConfig(cnfg_file_temp2);
                    writeConfig(cnfg_file_temp);
                }
            }
            for (ir = 0; ir < measParams.isms; ir++)
            {
                for (k = 0; k < VOL_SL; k++)
                {
                    updateArrayOfIndexMF(sl[k]);
                    readConfig(cnfg_file_temp);
                    smearing_APE_spatial(1, measParams.smpars, pu); // iu2
                    writeConfig(cnfg_file_temp2);
                }
                for (k = 0; k < VOL_SL; k++)
                {
                    updateArrayOfIndexMF(sl[k]);
                    readConfig(cnfg_file_temp2);
                    writeConfig(cnfg_file_temp);
                }
            }

            t2 = getTime();
            logging("\nAPE-Smearing complete (took %.3e sec)\n\n", t2 - t1);

            measureWilsonLoop(wil);
            int nn = 0, inn;
            for (ir = 0; ir < nr; ir++, nn++)
            {
                inn = nn * nt;
                for (it = 0; it < nt; it++, inn++)
                {
                    logging("WIL %d %d\t%.8e\n", measParams.rs + ir * measParams.dr, measParams.ts + it * measParams.dt, wil[inn]);
                }
            }
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
        free(wil);
    if (runParams.mcorrs)
    {
        deallocateFermionFieldsFor2ptFunctions();
        free(measParams.g1);
        free(measParams.g2);
    }

    releaseGaugeField();

    return 0;
}
#elif (MASTER_FIELD == 0)
int main(void)
{
    printf("Abort: use non-masterfield version of this main in _mf folder");
}
#endif