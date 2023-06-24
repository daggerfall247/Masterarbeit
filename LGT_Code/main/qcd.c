/*******************************************************************************
 *
 * File qcd.c
 *
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

int main(int argc, char *argv[])
{
    int n, nw;
    int seed, rconf, tconf, start;
    char out_dir[NAME_SIZE];
    char cnfg_file[FULL_PATH_SIZE + 12];
    //double qtop;
    double plaq, t1, t2;

    readInputFile(&seed, out_dir, &start, &rconf, &tconf, cnfg_file, argc, argv);
    rlxd_init(2, seed);
    initProgram(1);

    printf("PI after init Program: PI = %1.15f\n\n", PI);

    setupOutputFiles(runParams.idForOutputFilesName, out_dir);
    printStartupInfo(seed, rconf, cnfg_file);
    initArrayOfNeighbours();
    initGaugeField(start);

    if (rconf)
        readConfig(cnfg_file);

    if (tconf)
        readConfigForThermalisation(cnfg_file);

    checkForErrors(1, 0);

    for (n = 0, nw = 0; n < runParams.numConfs; n++)
    {
        t1 = getTime();
        gaugefieldUpdate(n + 1, 1, SIM_TYPE, SWEEP_TYPE);
        t2 = getTime();
        plaq = plaquette();
        logging("plaq\t%.8f\n", plaq);
        //qtop = topologicalCharge();
        //logging("qtop\t%.8f\n", qtop);

        if ((((n - runParams.numThermConfs + 1) % runParams.writeConfsFreq) == 0) && (n > runParams.numThermConfs) && (runParams.writeConfsFreq > 0))
        {
            nw++;
            sprintf(cnfg_file, "%s_n%d", CNFG_FILE, nw);
            writeConfig(cnfg_file);
        }

       
        printf("\nSweep %d done (took %.3e sec)!\n", n, t2 - t1);
        checkForErrors(1, 0);
    }

    releaseGaugeField();

    return 0;
}
