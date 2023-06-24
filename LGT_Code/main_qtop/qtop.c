//
//  meas_topcharge.c
//  glueballs
//
//  Created by Carolin Riehl on 15.05.20.
//  Copyright © 2020 Carolin Riehl. All rights reserved.
//

#define MAIN_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ranlxd.h"
#include "headers.h"
#include "modules.h"

int main(int argc, char *argv[])
{
    int n;
    int N_APE;
    int conf_start, conf_finish, conf_step;
    int N_APE_start, N_APE_finish, N_APE_step;

    char out_dir[NAME_SIZE];
    char cnfg_file[FULL_PATH_SIZE + 12];

    double t1, t2;
    double qtop;

    double alpha = 0.1;
    // sprintf(out_dir, "/scratch/mesonqcd/stilger/MA/configs/48x48x48x48_SU2_b3.09/temp");
    sprintf(out_dir, "../Configs/16x16x16x16_SU2_b2.300");

    if (argc == 2)
    {
        strcpy(out_dir, argv[1]);
        printf("%s\n\n", out_dir);
    }

    conf_start = 1;
    conf_finish = 1;
    conf_step = 1;

    N_APE_start = 0;
    N_APE_finish = 20;
    N_APE_step = 1;

    runParams.beta = 2.3;

    // ###############################
    runParams.decorSteps = 1;
    runParams.numConfs = 1;
    runParams.numThermConfs = 0;
    runParams.writeConfsFreq = 1;
    // ###############################

    rlxd_init(2, 1337);

    initProgram(1);
    setupOutputFilesForMeas(1, out_dir);
    initArrayOfNeighbours();
    initGaugeField(0);

    printf("PI after init Program: PI = %1.15f\n\n", PI);

    for (n = conf_start; n <= conf_finish; n += conf_step)
    {

        t1 = getTime();
        printf("Configuration : %04d\n", n);
        sprintf(cnfg_file, "%s_n%d", CNFG_FILE, n);
        printf("%s\n", cnfg_file);
        logging("measuring on Configuration %s", cnfg_file);

        readConfig(cnfg_file);
        smearing_APE_all(N_APE_start, alpha, pu);
        for (N_APE = N_APE_start; N_APE <= N_APE_finish; N_APE += N_APE_step)
        {
            // Measurement.
            qtop = topologicalCharge();
            // meas_topologicalcharge(&qtop);

            // Printing.
            logging("QTOP %d %d %e %+.8e\n", n, N_APE, alpha, qtop);
            // printf("QTOP %d %d %e %+.8e\n", n, N_APE, alpha, qtop);

            // Smearing.
            if(N_APE != N_APE_finish)
                smearing_APE_all(N_APE_step, alpha, pu);
        }
        t2 = getTime();
        logging("measurement took %.3e sec\n\n", t2 - t1);
    }

    releaseGaugeField();

    return 0;
}
