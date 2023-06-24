/*******************************************************************************
 *
 * Copyright (C) 2020 Alessandro Sciarra
 * Copyright (C) 2019 Francesca Cuteri
 * Copyright (C) 2016 Sebastian Schmalzbauer
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Integrated autocorrelation time estimate for data read from a (one column) file
 *
 *  compile like:
 *  cc autocorr.c -o autocorr -lm
 *
 *  run:
 *  ./autocorr <data.txt>
 *
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double average(double *array, int start, int end) {
    double average = 0.0;
    for(int i = start; i < end; i++)
        average += array[i];
    return average / (end - start);
}


double autocorrelationFunction(double *array, int N, int k) {
    // See section 13.7 and equations (13.48) (13.49) of
    //   http://www.physics.ntua.gr/konstant/ComputationalPhysics/C++/Book/ComputationalPhysicsKNA2ndEd_nocover.pdf
    double c_k = 0.0;
    double avg_0 = average(array, 0, N - k);
    double avg_k = average(array, k, N);
    for(int i = 0; i < N - k; i++)
        c_k += (array[i] - avg_0) * (array[i + k] - avg_k);
    return c_k / (N - k);
}


double integratedAutocorrelationTime(double *array, int N) {
    double gamma_0 = autocorrelationFunction(array, N, 0);
    double t = 0.5 * gamma_0;
    double gamma_k;
    for(int k = 1; k < N; k++) {
        gamma_k = autocorrelationFunction(array, N, k);
        if(gamma_k < 0.0)
            break;
        t += gamma_k;
    }
    return t / gamma_0;
}


double naiveError(double *array, int N) {
    return sqrt(autocorrelationFunction(array, N, 0)/(N - 1));
}


void analyseDataAndPrintResult(double *array, int N) {
    double tauInt = integratedAutocorrelationTime(array, N);
    int nTauIntToCompareWithStatistics = (int)ceil(100 * tauInt);
    if(nTauIntToCompareWithStatistics > N) {
        printf("Not enough measurements to resolve the estimated tau_int!\n");
        exit(1);
    }
    else {
        double trueError = sqrt(2.0 * tauInt) * naiveError(array, N);
        printf("%f  +-  %f    with tauInt = %f\n", average(array, 0, N), trueError, 2.0 * tauInt);
    }
}


int checkThatInputFileExistsAndGetNumberOfLines(char *filename) {
    char line[BUFSIZ];
    int nlines = 0;

    FILE *check = fopen(filename, "r");
    if (!check) {
        printf("Error: could not open %s!\n", filename);
        exit(1);
    }

    while(fgets(line, BUFSIZ, check) != NULL) {
        nlines++;
    }
    return nlines;
    fclose(check);
}


void readDataFromFile(char *filename, double *data) {
    int i = 0;
    int dummy = 0;
    FILE *read = fopen(filename, "r");
    while (fscanf(read, "%lf\n", &data[i]) != EOF)
        i++;
    fclose(read);
}


int main(int argc,char *argv[]) {
    if(argc != 2) {
        printf("Specify parameters: ./%s <data.txt>\n", argv[0]);
        exit(1);
    }

    int N = checkThatInputFileExistsAndGetNumberOfLines(argv[1]);
    printf("File %s has %i lines\n", argv[1], N);

    double *data = (double*)malloc(N * sizeof(double));
    if (data == NULL) {
        printf("ERROR : malloc failed ! Aborting ...\n");
        exit(1);
    }

    readDataFromFile(argv[1], data);
    analyseDataAndPrintResult(data, N);

    free(data);
    return 0;
}
