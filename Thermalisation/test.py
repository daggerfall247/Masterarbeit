import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t
from scipy.stats import ttest_ind

fileCold = "24x24x24x24_b3.09_WL_5_5_cold.dat"
fileHot = "24x24x24x24_b3.09_WL_5_5_hot.dat"

#file_out = "variance_24x24x24x24_b3.09_WL_9_9.dat"
intervall_length = 10
diff_intervall_length = 10
xrange = 2000

def variance(values):
    #mean = 0.0
    mean = np.mean(values)
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


if __name__ == "__main__":
    valuesCold = []
    valuesHot = []

    for line in open(fileCold):
        valuesCold.append(float(line))

    for line in open(fileHot):
        valuesHot.append(float(line))

    
    diff = [(valuesCold[ii] - valuesHot[ii]) for ii in range(len(valuesCold))]

    print("\ncheckpoint0\n")

    fig1, axs = plt.subplots(2)
    axs[0].set_xlim([0,xrange])
    axs[0].set_ylim([-0.0025,0.02])
    axs[0].plot(range(len(valuesCold)), diff)
    axs[0].axhline(y = 0.0, color = 'black', linestyle = '-')

    axs[1].set_xlim([0,xrange])
    axs[1].set_ylim([0.025,0.075])
    axs[1].plot(range(len(valuesCold)), valuesCold)
    axs[1].plot(range(len(valuesHot)), valuesHot)


    print("\ncheckpoint1\n")

    fig1.show()

    plt.show()

    
