import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t

fileCold = "24x24x24x24_b3.09_WL_1_1_cold.dat"
fileHot = "24x24x24x24_b3.09_WL_1_1_hot.dat"

#file_out = "variance_24x24x24x24_b3.09_WL_9_9.dat"
intervall_length = 1000
diff_intervall_length = 10
xrange = 50000
t_int = 2.0

def variance(values):
    mean = 0.0
    #mean = np.mean(values)
    N = len(values)

    summe = 0.0
    for ii in range(N):
        summe += (values[ii] - mean)*(values[ii] - mean)

    
    return summe/ (N-1)


def jackknife(values):
    mean = np.mean(values)
    N = len(values)

    summe = 0.0
    for ii in range(N):
        summe += (values[ii] - mean)*(values[ii] - mean)

    error = math.sqrt(summe/ (N*(N-1)))
    
    return mean, error

def jackknife2(values):
    mean = np.mean(values)
    N = len(values)

    summe = 0.0
    for ii in range(N):
        summe += (values[ii] - mean)*(values[ii] - mean)

    error = math.sqrt(t_int * summe/ (N*(N-1)))
    
    return mean, error


if __name__ == "__main__":
    valuesCold = []
    valuesHot = []

    for line in open(fileCold):
        valuesCold.append(float(line))

    for line in open(fileHot):
        valuesHot.append(float(line))

    
    diff = []
    variances = []
    for ii in range(0,50000-intervall_length+1,1):
        valuesForFitCold = []
        valuesForFitHot = []
        
        for jj in range(ii, ii+intervall_length,1):
            valuesForFitCold.append(valuesCold[jj])
            valuesForFitHot.append(valuesHot[jj])

        meanHot, errorHot = jackknife(valuesForFitHot)
        meanCold, errorCold = jackknife(valuesForFitCold)
        diff.append((meanCold-meanHot)/math.sqrt((errorCold**2+errorHot**2)))

    print("\ncheckpoint0\n")

    fig1, axs = plt.subplots(2)
    axs[0].set_xlim([0,xrange])
    axs[0].set_ylim([-6,6])
    axs[0].plot(range(0,50000-intervall_length+1,1), diff)
    axs[0].axhline(y = 0.0, color = 'black', linestyle = '-')
    axs[0].fill_between(range(0,xrange-intervall_length+1,1), -1, 1, color='orange', alpha=0.6)
    axs[0].fill_between(range(0,xrange-intervall_length+1,1), -2, 2, color='orange', alpha=0.4)
    axs[0].fill_between(range(0,xrange-intervall_length+1,1), -3, 3, color='orange', alpha=0.2)

    print("\ncheckpoint1\n")

    for ii in range(0,len(diff)-diff_intervall_length+1,1):
        diffForVariance = []
        for jj in range(ii, ii+diff_intervall_length,1):
            diffForVariance.append(diff[jj])

        variances.append(variance(diffForVariance))

    print("\ncheckpoint2\n")

    axs[1].set_xlim([0,xrange])
    axs[1].set_ylim([0,3])
    axs[1].plot(range(0,len(diff)-diff_intervall_length+1,1), variances)
    axs[1].axhline(y = 1.0, color = 'black', linestyle = '-')
    axs[1].axhline(y = math.sqrt(2.0), color = 'black', linestyle = '-')
    fig1.show()

    plt.show()

    
