# -*- coding: utf-8 -*-
"""
Thermal Project Code:
    Using data from 'Cu3AuData.xslx', which contains the temperature (K) and
    background-subtracted specific molar heat capactity (R, gas constant) of
    Cu3Au, extract the excess entropy associated with the data.
"""
import math
import numpy as np
import scipy.signal as ss
import matplotlib.pyplot as plt


def plotter(ylabel, xlabel, filename, x, y, yerr, xerr, peak):
    """
    Plotting Function
    Creates a figure and plots (x,y), includes axis labels
    Parameters
    ----------
    ylabel : y_axis label, str.
    xlabel : x_axis label, str.
    filename : Name of file which saves the figure as, str.
    x : list of data points for the x coordinates, list.
    y : list of data points for the y coordinates, list.
    yerr : error in the y coordinates, list.
    xerr : error in the x coordinates, list.
    peak: Peak of graph, int.
    """
    fig, ax = plt.subplots()
    plt.plot(x, y)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.errorbar(x, y, yerr=yerr, xerr=xerr)
    ax.annotate(f"T_c = {x[peak]} K", (x[peak], y[peak]), (x[peak]-75, y[peak-5]), arrowprops=dict(shrink=0.05))
    plt.show()
    fig.savefig(filename)


def integration(c, t1, t2):
    """
    Integration function


    Parameters
    ----------
    c : specific molar heat per T at i
        float
    t1 : Temperature at i
        float
    t2 : Temperature at i + 1
        float

    Returns
    -------
    s : Molar Entropy change element at i
        float

    """
    s = c*(t2 - t1)
    return s


def uncertainty_T(t):
    """
    Using the uncertainty equation from "article" Appendix A 
    for the temperature measurement:
        sigma(T) = T*(0.01 + T*(9*10**-7))

    Parameters
    ----------
    t : Temperature measurement, in Kelvin
        float

    Returns
    -------
    sig : Uncertainty in Temperature measurement, in Kelvin
        float

    """
    sig = t*(0.01+t*9e-7)
    return sig


def uncertainty_CbyT(t, c, s):
    """
    Uncertainty calculation for specific molar heat capacity per temperature,
    in R/K
    Using sigma(Cp) = Cp*(-0.198 + T*(9.4*10**-4) - (t**2)*(9.2*10**-7))
    and sigma(f)/f = sqrt[ (sigma(x)/x)**2 + (sigma(y)/y)**2] for f = x*y

    Parameters
    ----------
    t : Temperature measurement, in Kelvin
        float
    c : Specific molar heat capacity measurement, in R
        float
    s : Specific molar heat per temeprature value, in R/K
        float

    Returns
    -------
    sig : Uncertainty in specific molar heat per temperature value
        float

    """
    # uncertainty in Cp, sig = Cp*(-0.198+t*9.4e-4 - (t**2)*9.2e-7)
    # f = a*b, (sig_f)/f = sqrt[ (sig_a/a)**2 + (sig_b/b)**2 ]
    sig_c = c*(-0.198 + t*(9.4*10**-4) - (t**2)*(9.2*10**-7))
    sig_t = uncertainty_T(t)
    sig = s*math.sqrt((sig_c/c)**2 + (sig_t/t)**2)
    return sig


def uncertainty_s(t, c, s1, s2):
    """
    Uncertainty in Entropy calculation.
    This is for a single element of the entropy calculation.

    Parameters
    ----------
    t : Temperature measurement, in Kelvin.
    c : Specific molar heat measurement, in R.
    s1 : Specific molar heat per kelvin value, in R/K.
    s2 : Entropy value, in R.

    Returns
    -------
    sig : Uncertainty in entropy calculation for one instance, in R.

    """
    sig_c = uncertainty_CbyT(t, c, s1)
    sig_t = uncertainty_T(t)
    sig = s2*math.sqrt((sig_c/c)**2 + (sig_t/t)**2)
    return sig


def main():
    # Plot Cp against T and Cp/T against T
    infile = "Thermal_data.txt"
    data = open(infile, "r")
    t_data = []  # temperature data
    t_err = []  # temperature uncertainty
    c_data = []  # specific molar heat data
    y_data = []  # specific molar heat per temp data
    y_err = []  # specific molar heat per temperature uncertainty
    for line in data.readlines():
        if not line.startswith("#"):
            tokens = line.split()
            x = float(tokens[0])
            y = float(tokens[1])
            z = y/x
            t_data.append(x)
            t_err.append(uncertainty_T(x))
            c_data.append(y)
            y_data.append(z)
            y_err.append(uncertainty_CbyT(x, y, z))
    s = 0  # Molar entropy change
    s_err = []  # error in entropy calculation
    start = 1137
    end = 2145
    y_new = []
    x_new = []
    for i in range(start, end):
        s += integration(y_data[i], t_data[i], t_data[i+1])
        s_err.append(uncertainty_s(t_data[i], c_data[i], y_data[i], s))
        y_new.append(y_data[i])
        x_new.append(t_data[i])

    # Uncertainty in s = sqrt(sum(ds**2))
    ds = []
    for i in range(len(s_err)):
        ds.append(s_err[i]**2)
    s_err = sum(ds)
    s_err = math.sqrt(s_err)
    print(f"Molar Entropy change, S = {s:.02} + {s_err:.01} R")
    print("Integration bounds:")
    print(f"start = {t_data[start]} K, end = {t_data[end]} K")

    # Find Critical Temperature
    peak,_ = ss.find_peaks(y_data, height=0.025)
    peak = np.array(peak)
    print(f"Critical Temperature: T_c = {t_data[peak[0]]} K.")
    
    # Figure 1 shows Cp against T
    fig1 = plotter("Background subtracted specific molar heat capacity (R)",
                   "Temperature (K)", "Thermal_fig_1", t_data, c_data, 0, 0, peak)

    # Figure 2 shows Cp/T against T
    fig2 = plotter("Specific Molar heat capacity per temperature (R/K)",
                   "Temperature (K)", "Thermal_fig_2", t_data, y_data, 0, 0, peak)

    # Figure 3 shows Cp/T against T within the integration bounds (start, end)
    fig3 = plotter("Specific Molar heat capacity per temperature (R/K)",
                   "Temperature (K)", "Thermal_fig_3", x_new, y_new, 0, 0, peak-start)


main()
