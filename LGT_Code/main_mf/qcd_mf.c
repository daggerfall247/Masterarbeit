/*******************************************************************************
 *
 * File qcd.c
 *
 * Copyright (C) 2022 Marco Stilger
 * Copyright (C) 2020 Alessandro Sciarra
 * Copyright (C) 2019 Francesca Cuteri
 * Copyright (C) 2016 Bastian Brandt
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Program under construction for lattice gauge theory simulations
 *
 * Syntax: qcd -i [input-filename]
 *
 * For usage instructions see the file ../README.txt
 *
 *******************************************************************************/

#define MAIN_C

#include "ranlxd.h"
#include "modules.h"
#include <stdio.h>
#include <stdbool.h>

int main(int argc, char *argv[])
{
    int n, nw;
    uint64_t k;
    bool write;
    int seed, rconf, tconf, start;
    char out_dir[NAME_SIZE];
    char cnfg_file[FULL_PATH_SIZE + 12], cnfg_file_temp[FULL_PATH_SIZE + 12];
    double plaq;

    readInputFile(&seed, out_dir, &start, &rconf, &tconf, cnfg_file, argc, argv);
    rlxd_init(2, seed);
    initProgram(start);

    printf("PI after init Program: PI = %1.15f\n\n", PI);

    setupOutputFiles(runParams.idForOutputFilesName, out_dir);
    printStartupInfo(seed, rconf, cnfg_file);

    initGaugeField(0);
    initGlobalArrays();

   
    /*
    if (rconf)
        readConfig(cnfg_file);
    */

    
    checkForErrors(1, 0);

    double t1, t2;

    sprintf(cnfg_file_temp, "%s_temp", CNFG_FILE);
    prepareConfig(cnfg_file_temp); // auskommentieren, wenn man mittendrin anfängt

    for (n = 0, nw = 0; n < runParams.numConfs; n++) // nw anpassen wenn man mittendrin anfängt
    {
        t1 = getTime();
        plaq = 0.0;
        write = (((n - runParams.numThermConfs + 1) % runParams.writeConfsFreq) == 0) && (n > runParams.numThermConfs) && (runParams.writeConfsFreq > 0);

        if (write)
        {
            nw++;
            sprintf(cnfg_file, "%s_n%d", CNFG_FILE, nw);
            writeHeaderToConfig(cnfg_file);
        }

        for (k = 0; k < VOL_SL; k++)
        {
            updateArrayOfIndexMF(sl[k]);
            readConfig(cnfg_file_temp);
            gaugefieldUpdate(n + 1, 1, SIM_TYPE, SWEEP_TYPE); // adjust number sweeps
            writeConfig(cnfg_file_temp);
            if (write)
            {
                writeConfig(cnfg_file);
            }
        }

        plaq = 0.0;
        for (k = 0; k < VOL_SL; k++)
        {
            updateArrayOfIndexMF(sl[k]);
            readConfig(cnfg_file_temp);
            plaq += plaquette();
        }
        plaq /= VOL_SL;
        
        logging("plaq\t%.6e", plaq);
        writePlaquetteToConfig(cnfg_file_temp, plaq);
        if (write)
        {
            writePlaquetteToConfig(cnfg_file, plaq);
        }

        checkForErrors(1, 0);
        t2 = getTime();
        printf("\nSweep %d done (took %.3e sec)!\n", n, t2 - t1);
    }
    

    releaseGaugeField();

    return 0;
}
