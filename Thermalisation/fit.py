import math
import numpy as np

fileCold = "24x24x24x24_b3.09_WL_7_7_cold.dat"
fileHot = "24x24x24x24_b3.09_WL_7_7_hot.dat"

file_out = "fit_10_24x24x24x24_b3.09_WL_7_7.dat"
intervall_length = 10

def jackknife(values):
    mean = np.mean(values)
    N = len(values)
    sum = 0.0
    for ii in range(N): #for all subsets
        subset = [values[x] for x in range(N) if x is not ii]
        subset_mean = np.mean(subset)
        sum += (subset_mean - mean)**2

    error = math.sqrt((N-1)/N * sum)
    return mean, error


if __name__ == "__main__":
    valuesCold = []
    valuesHot = []

    for line in open(fileCold):
        valuesCold.append(float(line))

    for line in open(fileHot):
        valuesHot.append(float(line))
    
    outfile = open(file_out, "w")
    for ii in range(0,50000-intervall_length+1):
        valuesForFitCold = []
        valuesForFitHot = []
        for jj in range(ii, ii+intervall_length):
            valuesForFitCold.append(valuesCold[jj])
            valuesForFitHot.append(valuesHot[jj])

        meanCold, errorCold = jackknife(valuesForFitCold)
        meanHot, errorHot = jackknife(valuesForFitHot)

        print(f'{ii:6d} {meanCold:+.8e} {meanHot:+.8e} {errorCold:+.8e} {errorHot:+.8e}', file=outfile)

    outfile.close()

    
