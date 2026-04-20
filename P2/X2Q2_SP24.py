# region imports
from math import sin
import math
import numpy as np
import scipy as sp
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
# endregion

# region function definitions
def odeSystem(t, X, *args):
    """
    this is the odeSystem callback I'm using for solve_ivp().  The *args contains a callback for the input voltage.
    :param X: the current values of the state variables
    :param t: the current time
    :param args: fn(t), L, R, C
    :return: list of derivatives of state variables
    """
    # ChatGPT helped write this function.

    # unpack X into convenient variables
    fn, L, R, C = args

    # assign friendly variable names
    i1 = X[0]
    i2 = X[1]

    # calculate the current input voltage
    vt = fn(t)

    # From the circuit:
    # iL = i1 + i2
    # vR = R*(i1 - i2) = vC
    # L*d(i1+i2)/dt = v(t) - vR
    # i2 = C*dvC/dt = C*d(R*(i1-i2))/dt
    #
    # This gives:
    # i1dot + i2dot = (vt - R*(i1 - i2))/L
    # i1dot - i2dot = i2/(R*C)

    a = (vt - R * (i1 - i2)) / L
    b = i2 / (R * C)

    i1dot = 0.5 * (a + b)
    i2dot = 0.5 * (a - b)

    return [i1dot, i2dot]


def simulate(L=20, R=10, C=0.05, A=20, f=20, p=0, t=10, pts=500):
    """
    For simulating transient behavior of circuit.
    :param L: Inductance (H)
    :param R: Resistance (ohm)
    :param C: Capacitance (F)
    :param A: Amplitude (V)
    :param f: frequency (Hz)
    :param p: phase (deg)
    :param t: simulation time (s)
    :param pts: number of time points
    :return: solve_ivp solution object
    """
    # ChatGPT helped write this function.
    w = f * 2 * math.pi  # frequency in radians/sec
    phi = p * math.pi / 180.0  # phase in radians
    vin = lambda t: A * sin(w * t + phi)
    myargs = (vin, L, R, C)
    x0 = [0, 0]  # initial values state variables
    tList = np.linspace(0, t, int(pts))  # time vector
    I = solve_ivp(odeSystem, t_span=[0, t], y0=x0, t_eval=tList, args=myargs)
    return I


def doPlot(*args, ax=None):
    """
    Re-written on 4/21/2022 to adapt to plotting on GUI if ax is not None
    :param args: contains ((R, list of time values, and results of solve_ivp))
    :param ax:
    :return:
    """
    if ax is None:
        ax = plt.subplot()
        QTPlotting = False
    else:
        QTPlotting = True

    R, tList, I = args[0]
    ax.clear()
    ax.plot(tList, I.y[0], linestyle='solid', color='k', label=r'$i_1(t)$')
    ax.plot(tList, I.y[1], linestyle='dashed', color='k', label=r'$i_2(t)$')
    ax.set_xlim(0, max(tList))
    minI = min(min(I.y[0]), min(I.y[1]))
    maxI = max(max(I.y[0]), max(I.y[1]))
    rangeI = abs(maxI - minI)
    ax.set_ylim(minI - 0.01 * rangeI, maxI + 0.01 * rangeI)
    ax.tick_params(axis='both', which='both', direction='in', top=True, labelsize=12)
    ax.tick_params(axis='both', grid_linewidth=1, grid_linestyle='solid', grid_alpha=0.5)
    ax.tick_params(axis='both', which='minor')
    ax.grid(True)
    ax.set_xlabel('t (s)', fontsize=12)
    ax.set_ylabel(r'$i_1, i_2$ (A)', fontsize=12)

    ax1 = ax.twinx()
    yvals = R * (I.y[1] - I.y[0])
    yrange = abs(max(yvals) - min(yvals))
    ax1.plot(tList, yvals, linestyle='dotted', color='k', label=r'$v_c(t)$')
    ax1.set_ylim(min(yvals) - yrange * 0.01, max(yvals) + yrange * 0.01)
    ax1.tick_params(axis='y', which='both', direction='in', top=True, right=True, labelsize=12)
    ax.legend(fontsize=12)
    ax1.legend(loc='lower right', fontsize=12)
    ax1.set_ylabel(r'$V_c(t)$ (V)', fontsize=12)

    if not QTPlotting:
        plt.show()


def main():
    """
    For solving problem 2 on exam.
    :return:
    """
    I = simulate(L=20, R=10, C=0.05, A=20, f=20, p=0, t=10, pts=500)
    doPlot((10, I.t, I))
# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion