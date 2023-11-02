# -*- coding: utf-8 -*-
"""
Thermal Project Code:
    Using data from 'Cu3AuData.xslx', which contains the temperature (K) and
    background-subtracted specific molar heat capactity (R, gas constant) of
    Cu3Au, extract the excess entropy associated with the data.
"""
import math
import matplotlib.pyplot as plt


def plot(i, ylabel, xlabel, filename, x, y, yerr, xerr):
    """
    Plotting Function
    Creates a figure and plots (x,y), includes axis labels
    Parameters
    ----------
    i : Figure Number, int.
    ylabel : y_axis label, str.
    xlabel : x_axis label, str.
    filename : Name of file which saves the figure as, str.
    x : list of data points for the x coordinates, list.
    y : list of data points for the y coordinates, list.
    yerr : error in the y coordinates, list.
    xerr : error in the x coordinates, list.
    """
    fig = plt.figure(i)
    plt.plot(x, y)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    fig.savefig(filename)
    plt.errorbar(x, y, yerr=yerr, xerr=xerr)
    plt.show()


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
    s = 0  # Molar entropy change
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
    s_err = []  # error in entropy calculation
    for i in range(len(t_data)-1):
        s += integration(y_data[i], t_data[i], t_data[i+1])
        s_err.append(uncertainty_s(t_data[i], c_data[i], y_data[i], s))

    # Figure 1 shows Cp against T
    fig1 = plot(1, "Background subtracted specific molar heat capacity (R)",
                "Temperature (K)", "Thermal_fig_1", t_data, c_data, 0, t_err)

    # Figure 2 shows Cp/T against T
    fig2 = plot(2, "Specific Molar heat capacity per temperature (R/K)",
                "Temperature (K)", "Thermal_fig_2", t_data, y_data, y_err, t_err)

    # Uncertainty in s = sqrt(sum(ds**2))
    ds = []
    for i in range(len(s_err)):
        ds.append(s_err[i]**2)
    s_err = sum(ds)
    s_err = math.sqrt(s_err)
    print(f"Molar Entropy change, S = {s:.02} ({round(10*s_err)}) R")
    print("Value given in article, S = 0.45 (4) R")
    diff = abs(0.45-s)/s_err
    print(f"Calculated value is different by {diff:.03} error values")


main()