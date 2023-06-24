import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

file1 = "V_24.dat"
file2 = "V_48.dat"
file3 = "V_96.dat"

loop_r_min = 1
loop_r_max = 12
loop_dr = 1

num_r_loops = loop_r_max-loop_r_min+1

if __name__ == "__main__":

    V_24 = [[] for rr in range(num_r_loops)]
    V_48 = [[] for rr in range(num_r_loops)]
    V_96 = [[] for rr in range(num_r_loops)]
    
    for rr in range(num_r_loops):
        for line in open(file1):
            if line.startswith("%d\t" % (rr*loop_dr + loop_r_min)):
                splitted_line = line.split()
                V_24[rr].append(float(splitted_line[2]))
                V_24[rr].append(float(splitted_line[3]))

    for rr in range(num_r_loops):
        for line in open(file2):
            if line.startswith("%d\t" % (rr*loop_dr + loop_r_min)):
                splitted_line = line.split()
                V_48[rr].append(float(splitted_line[2]))
                V_48[rr].append(float(splitted_line[3]))
                
    for rr in range(num_r_loops):
        for line in open(file3):
            if line.startswith("%d\t" % (rr*loop_dr + loop_r_min)):
                splitted_line = line.split()
                V_96[rr].append(float(splitted_line[2]))
                V_96[rr].append(float(splitted_line[3]))

    ratio = []
    ratio_err = []

    for rr in range(num_r_loops):
        ratio.append((V_24[rr][0]-V_48[rr][0])/(V_48[rr][0]-V_96[rr][0]))
        ratio_err.append(math.sqrt((V_24[rr][0]-V_48[rr][0])**2/(V_48[rr][0]-V_96[rr][0])**4 * V_96[rr][1]**2 + (V_96[rr][0]-V_24[rr][0])**2/(V_48[rr][0]-V_96[rr][0])**4 * V_48[rr][1]**2 + 1/(V_48[rr][0]-V_96[rr][0])**2 * V_24[rr][1]**2))

    for rr in range(num_r_loops):
        print("%d\t%f\t%f" % (rr*loop_dr + loop_r_min, ratio[rr], ratio_err[rr]))