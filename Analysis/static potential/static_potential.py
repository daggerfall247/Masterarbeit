import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

file = "24x24x24x24_SU2_b3.09_NAPE_400.meas"
#file = "48x48x48x48_SU2_b3.09_NAPE_400.meas"
#file = "96x96x96x96_SU2_b3.09_NAPE_400.meas"

Vol = 24**4
#Vol = 48**4
#Vol = 96**4

num_confs = 1000
num_bins = 20

loop_r_min = 1
loop_r_max = 12
loop_dr = 1

loop_t_min = 1
loop_t_max = 20
loop_dt = 1

num_t_loops = loop_t_max-loop_t_min+1
num_r_loops = loop_r_max-loop_r_min+1

def fitting_constant(f, s):
    if(len(f) != len(s)):
        exit("ERROR")

    N = len(f)

    sigma = 0.0
    for ii in range(N):
        sigma += 1.0/(s[ii]*s[ii])

    a = 0.0
    for ii in range(N):
        a += f[ii]/(sigma*s[ii]**2)

    return a

def jackknife1(values):
    mean = np.mean(values)
    N = len(values)
    sum = 0.0
    for ii in range(N): #for all subsets
        subset = [values[x] for x in range(N) if x is not ii]
        subset_mean = np.mean(subset)
        sum += (subset_mean - mean)**2

    error = math.sqrt((N-1)/N * sum)
    return mean, error

def binData(array):
    N = num_confs
    L = N//num_bins
    print(N)
    print(L)
    return np.array([np.mean(array[N%L:])] + [np.mean(np.delete(array[N%L:], np.s_[L*k:L*(k+1)], 0)) for k in range(num_bins)])

if __name__ == "__main__":

    WL = [[[] for x in range(num_t_loops)] for y in range(num_r_loops)]
    
    for rr in range(num_r_loops):
        for tt in range(num_t_loops):
            for line in open(file):
                if line.startswith("WIL %d %d\t" % (rr*loop_dr + loop_r_min, tt*loop_dt + loop_t_min)):
                    splitted_line = line.split()
                    WL[rr][tt].append(float(splitted_line[3]))
            
            if (len(WL[rr][tt]) != num_confs):
                exit()

            #print(len(WL[rr][tt]))


    WL_jack_averages = [[[] for x in range(num_t_loops)] for y in range(num_r_loops)] # Wilson loop averages of reduced jackknife samples
    WL_averages = [[np.mean(WL[rr][tt]) for tt in range(num_t_loops)] for rr in range(num_r_loops)] # Wilson averages of complete sample
    for rr in range(num_r_loops):
        for tt in range(num_t_loops):
            for kk in range(num_bins):
                WL_jack_averages[rr][tt].append(np.mean(np.delete(WL[rr][tt], np.s_[num_confs//num_bins*kk:num_confs//num_bins*(kk+1)], 0)))
            
            #print(len(WL_jack_averages[rr][tt]))


    V_eff_jack = [[[] for tt in range(num_t_loops-1)] for rr in range(num_r_loops)] # effective potential of jackknife samples
    V_eff = [[math.log(WL_averages[rr][tt]/WL_averages[rr][tt+1]) for tt in range(num_t_loops-1)] for rr in range(num_r_loops)] # effective potential of complete sample
    for rr in range(num_r_loops):
        for tt in range(num_t_loops-1):
            for ii in range(num_bins):
                V_eff_jack[rr][tt].append(math.log(WL_jack_averages[rr][tt][ii]/WL_jack_averages[rr][tt+1][ii]))
            
            #print(len(V_eff_jack[rr][tt]))

    V_eff_jack_averages = [[np.mean(V_eff_jack[rr][tt]) for tt in range(num_t_loops-1)] for rr in range(num_r_loops)] # averages of effective potential over all jackknife samples
    sigma_V_eff = [[[] for tt in range(num_t_loops-1)] for rr in range(num_r_loops)] # error of effective potentials
    for rr in range(num_r_loops):
        for tt in range(num_t_loops-1):
            sum = 0
            for ii in range(num_bins):
                sum += (V_eff_jack[rr][tt][ii] - V_eff[rr][tt])**2
            
            sigma_V_eff[rr][tt] = math.sqrt((num_bins-1)/num_bins * sum)


    for rr in range(num_r_loops):
        for tt in range(num_t_loops-1):
            V_eff[rr][tt] = V_eff[rr][tt] - (num_bins - 1)*(V_eff_jack_averages[rr][tt] - V_eff[rr][tt]) # correct for bias
            print('%d\t%d\t%.6e' % (rr*loop_dr+loop_r_min, tt*loop_dt+loop_t_min, V_eff[rr][tt]))
        print('\n')


    fit_jack = [[] for rr in range(num_r_loops)] # fitted constant for jackknife samples
    V = []
    t_fit_start = [7,7,7,7,7,7,7,7,7,8,8,8]
    t_fit_end = [17,17,17,17,17,17,17,17,17,17,17,17]
    for rr in range(num_r_loops):
        X_for_fit = range(t_fit_start[rr], t_fit_end[rr])
        for kk in range(num_bins):
            Y_for_fit = [V_eff_jack[rr][ii][kk] for ii in range(t_fit_start[rr], t_fit_end[rr])]
            Y_err_for_fit = [sigma_V_eff[rr][ii] for ii in range(t_fit_start[rr], t_fit_end[rr])]
            a = fitting_constant(Y_for_fit, Y_err_for_fit)
            fit_jack[rr].append(a)

        Y_for_fit = [V_eff[rr][ii] for ii in range(t_fit_start[rr], t_fit_end[rr])]
        Y_err_for_fit = [sigma_V_eff[rr][ii] for ii in range(t_fit_start[rr], t_fit_end[rr])]
        a = fitting_constant(Y_for_fit, Y_err_for_fit)
        V.append(a)
        #print(len(tertiary[rr]))

    sigma_V = []
    for rr in range(num_r_loops):
        sum = 0
        for ii in range(num_bins):
            sum += (fit_jack[rr][ii] - V[rr])**2
            
        sigma_V.append(math.sqrt((num_bins-1)/num_bins * sum))

    fit_jack_averages = [np.mean(fit_jack[rr]) for rr in range(num_r_loops)]
    for rr in range(num_r_loops):
        V[rr] = V[rr] - (num_bins-1)*(fit_jack_averages[rr] - V[rr])
    

    for rr in range(num_r_loops):
        print('%d\t%.6e\t%.6e\t%.6e' % (rr*loop_dr+loop_r_min, 1.0/Vol, V[rr], sigma_V[rr]))




    fig2, ax2 = plt.subplots(num_r_loops, figsize=(10,60))
    fit_results = []
    fit_results_err = []
    X = range(loop_t_min, loop_t_max, loop_dt)
    for rr in range(num_r_loops):
        Y = [V_eff[rr][ii] for ii in range(num_t_loops-1)]
        Y_err = [sigma_V_eff[rr][ii] for ii in range(num_t_loops-1)]

        ax2[rr].errorbar(X, Y, yerr=Y_err)
        ax2[rr].errorbar(X, [V[rr] for ii in range(len(X))])
    
    fig2.savefig("V_eff_split.pdf", format='pdf')

    fm_in_GeV = 5.068
    a = 0.012899
    fig3, ax3 = plt.subplots()
    X = [range(loop_r_min, (loop_r_max+1), loop_dr)[rr]*a for rr in range(num_r_loops)]
    Y = [V[rr]/(a*fm_in_GeV) for rr in range(num_r_loops)]
    Y_err = [sigma_V[rr]/(a*fm_in_GeV) for rr in range(num_r_loops)]
    ax3.errorbar(X, Y, marker='x', yerr=Y_err, label=r'Test $\beta=3.09$')
    

    a_carolin = 0.026
    #offset = 0
    #offset = 2.325e6
    offset = 2.7225
    X_carolin = [range(2,12,1)[ii]*a_carolin for ii in range(10)]
    Y_raw = [0.353848, 0.394722, 0.417829, 0.433611, 0.445685, 0.455557, 0.463982, 0.471474, 0.478364, 0.484684]
    Y_carolin = [(Y_raw[ii]/(a_carolin*fm_in_GeV) + offset) for ii in range(10)]
    ax3.plot(X_carolin, Y_carolin, '-', label=r'Referenz $\beta=2.85$')

    ax3.set_xlabel(r"R [fm]")
    ax3.set_ylabel(r"V(R) [GeV]")
    ax3.legend(loc='lower right')

    fig3.savefig("V.pdf", format='pdf')