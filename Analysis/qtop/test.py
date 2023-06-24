import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

file1 = "24x24x24x24_SU2_b3.09_qtop_30.meas"
file2 = "48x48x48x48_SU2_b3.09_qtop_30.meas"
file3 = "96x96x96x96_SU2_b3.09_qtop_30.meas"

def variance(list):
    mean = np.mean(list)
    N = len(list)
    sum = 0
    for ii in range(N):
        sum += (mean - list[ii])**2

    sum = sum / (N)
    return sum

if __name__ == "__main__":

    V_24 = []
    V_48 = []
    V_96 = []
    
    
    for line in open(file1):
        if line.startswith("QTOP"):
            splitted_line = line.split()
            V_24.append(float(splitted_line[4]))

    
    for line in open(file2):
        if line.startswith("QTOP"):
            splitted_line = line.split()
            V_48.append(float(splitted_line[4]))
                
    
    for line in open(file3):
        if line.startswith("QTOP"):
            splitted_line = line.split()
            V_96.append(float(splitted_line[4]))

    var24 = variance(V_24)
    var48 = variance(V_48)
    var96 = variance(V_96)

    print(var24)
    print(var48)
    print(var96)