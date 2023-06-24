import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t
import scipy.stats as sc

intervall_lengths = [6, 11, 21, 51, 101, 201, 501, 1001, 2001]

fileCold = ["24x24x24x24_b3.09_WL_1_1_cold.dat", "24x24x24x24_b3.09_WL_3_3_cold.dat", "24x24x24x24_b3.09_WL_5_5_cold.dat", "24x24x24x24_b3.09_WL_7_7_cold.dat", "24x24x24x24_b3.09_WL_9_9_cold.dat", "24x24x24x24_b3.09_WL_10_16_cold.dat"]
fileHot = ["24x24x24x24_b3.09_WL_1_1_hot.dat", "24x24x24x24_b3.09_WL_3_3_hot.dat", "24x24x24x24_b3.09_WL_5_5_hot.dat", "24x24x24x24_b3.09_WL_7_7_hot.dat", "24x24x24x24_b3.09_WL_9_9_hot.dat", "24x24x24x24_b3.09_WL_10_16_hot.dat"]
loops = ["WL(1,1)", "WL(3,3)", "WL(5,5)", "WL(7,7)", "WL(9,9)", "WL(10,16)"]
t_int_hot = [1.656759, 4.718209, 7.002166, 4.738989, 1.579215, 1.000000]
t_int_cold = [1.765300, 5.257504, 9.178870, 7.109655, 2.366152, 1.033835]

alpha = 0.01

xrange = 50000

numLoops = len(fileCold)
numIntervals = len(intervall_lengths)

def r1(values):
    N = len(values)
    mean = np.mean(values)
    nom = 0.0
    denom = 0.0

    for t in range(N-1):
        nom += (values[t] - mean)*(values[t+1]-mean)

    for t in range(N):
        denom += (values[t] - mean) * (values[t] - mean)

    result = nom / denom
    return result

def t_value_1(values1, values2, t_int1, t_int2):
    mean1 = np.mean(values1)
    mean2 = np.mean(values2)
    s1 = variance(values1)
    s2 = variance(values2)
    #s1 = variance(values1)
    #s2 = variance(values2)
    N1 = len(values1)/t_int1
    N2 = len(values2)/t_int2

    #s1 = 0.0
    #for ii in range(len(values1)):
    #    s1 += (values1[ii] - mean1)**2

    #s2 = 0.0
    #for ii in range(len(values2)):
    #    s2 += (values2[ii] - mean2)**2

    #print(N1)
    #print(N2)
    error = math.sqrt((s1+s2)/2) * math.sqrt(1.0/N1 + 1.0/N2)

    result = (mean1-mean2)/error
    return result


def variance(values):
    #mean = 0.0
    mean = np.mean(values)
    N = len(values)

    summe = 0.0
    for ii in range(len(values)):
        summe += (values[ii] - mean)**2

    error = summe / (N-1)

    return error


def jackknife(values):
    mean = np.mean(values)
    N = len(values)

    summe = 0.0
    for ii in range(N):
        summe += (values[ii] - mean)*(values[ii] - mean)

    error = math.sqrt(summe / (N*(N-1)))

    return mean, error


def jackknife2(values):
    mean = np.mean(values)
    N = len(values)

    summe = 0.0
    for ii in range(N):
        summe += (values[ii] - mean)*(values[ii] - mean)

    error = math.sqrt(summe / (N*(N-1)))

    return mean, error


if __name__ == "__main__":

    figs = []
    for ll in range(numIntervals):
        fig, ax1 = plt.subplots(nrows=numLoops, ncols=1, sharex=True, sharey=True, figsize=(20,20))
        fig.set_dpi(400)
        fig.supxlabel(r"$\tau$ (MC time)")
        fig.supylabel(r"$\frac{\left<WL\right>_{cold} - \left<WL\right>_{hot}}{\sigma_p} \quad \sigma_p = \sqrt{\frac{\sigma^2_{cold} + \sigma^2_{hot}}{2}} \cdot \sqrt{\frac{1}{N_{1,e}}+\frac{1}{N_{2,e}}}$")
        plt.setp(ax1, xticks=[x for x in range(0,50001,2500)], yticks=[y for y in range(-8,9,2)], xlim=[0,xrange], ylim=[-7,7])
        fig.subplots_adjust(left=0.06, right=0.97, bottom=0.04, top=0.92, wspace=0.001, hspace=0.2)
        fig.suptitle(f"two sample t-test. H_0: <WL>_cold = <WL>_hot. H_0: <WL>_cold != <WL>_hot.\n Each t value is calculated on an interval of size {intervall_lengths[ll]} starting at the position of the t value. Shaded area is confidence interval with alpha = {alpha}")
        
        for kk in range(numLoops):

            valuesCold = []
            valuesHot = []
            valuesHot.clear()
            valuesCold.clear()

            for line in open(fileCold[kk]):
                valuesCold.append(float(line))

            for line in open(fileHot[kk]):
                valuesHot.append(float(line))


            ts = []
            for ii in range(0, 50000-intervall_lengths[ll]+1, 1):
                valuesForFitCold = []
                valuesForFitHot = []
                valuesForFitCold.clear()
                valuesForFitHot.clear()

                for jj in range(ii, ii+intervall_lengths[ll], 1):
                    valuesForFitCold.append(valuesCold[jj])
                    valuesForFitHot.append(valuesHot[jj])

                ts.append(t_value_1(valuesForFitCold, valuesForFitHot, t_int_cold[kk], t_int_cold[kk]))
        
            quantil = t.ppf(1-alpha/2, 2*intervall_lengths[ll]-2)

            #print(t.ppf(1-0.001/2, 20))

            ax1[kk].set_title(f"{loops[kk]}", y=1.0, pad=-14)
            #ax1[kk].set_xlabel(r"$\tau$ MC time")
            #ax1[kk].set_ylabel("")
            ax1[kk].plot(range(0, 50000-intervall_lengths[ll]+1, 1), ts, color = 'blue', label='t value')
            ax1[kk].fill_between(range(0,50000-intervall_lengths[ll]+1,1), -quantil, quantil, color='orange', alpha=0.3, label=f"shaded area y = [{-quantil, quantil}]")
            ax1[kk].axhline(y = 0, color = 'black', linestyle = '-')

        plt.savefig(f"WL_{intervall_lengths[ll]}_alt1.pdf", format='pdf')
        
