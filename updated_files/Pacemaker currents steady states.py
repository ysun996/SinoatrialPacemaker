import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def m_inf(Vm, max_a, max_b):
    a_m = max_a*((Vm+47)/(1-np.exp(-(Vm+47)/10)))
    b_m = max_b*np.exp(-(Vm + 72)/17.86)
    minf = a_m/(a_m+b_m)
    return minf

def h_inf(Vm, max_a, max_b):
    a_h = max_a*np.exp(-(Vm+71)/5.43)
    b_h = max_b/(1+np.exp(-(Vm+10)/12.2))
    hinf = a_h/(a_h + b_h)
    return hinf

def d_inf(Vm, max_a, max_b):
    a_d = max_a*((Vm+34)/(1-np.exp(-(Vm+34)/10)))
    b_d = max_b * np.exp(-(Vm+34)/6.67)
    dinf = a_d/(a_d + b_d)
    return dinf

def f_inf(Vm, max_a, max_b):
    a_f = max_a*np.exp(-(Vm+47)/20)
    b_f = max_b/(1 + np.exp(-(Vm+13)/11.49))
    finf = a_f/(a_f+b_f)
    return finf

def x1_inf(Vm, max_a, max_b):
    a_x1 = (max_a*np.exp((Vm+30)/12.11))/(1+np.exp((Vm+30)/50))
    b_x1 = max_b * (np.exp(-(Vm-20)/16.67))/(1 + np.exp(-(Vm - 20)/25))
    x1inf = a_x1/(a_x1+b_x1)
    return x1inf

def y_inf(Vm, max_a, max_b):
    a_y = max_a * np.exp(-(Vm + 17)/14.93)
    b_y = max_b * (Vm + 17)/(1 - np.exp(-(17 + Vm)/14.93)) ###Get your shit together why are your values messed up??
    yinf = a_y/(a_y+b_y)
    return yinf

def y_t(Vm, max_a, max_b):
    a_y = max_a * np.exp(-(Vm + 17)/14.93)
    b_y = max_b * (Vm + 17)/(1 - np.exp(-(17 + Vm)/14.93)) ###Get your shit together why are your values messed up??
    yinf = 1/(a_y+b_y)
    return yinf

#Max rate constants
max_a = [1, 0.0085, 0.005, 0.002468, 0.0025, 0.0000167]

max_b = [40, 2.5, 0.05, 0.05, 0.0065, 0.00033]

#Initial condition
initcon = np.arange(-100, 50, 0.5)

plotvar = [m_inf(initcon, max_a[0], max_b[0]),
           h_inf(initcon, max_a[1], max_b[1]),
           d_inf(initcon, max_a[2], max_b[2]),
           f_inf(initcon, max_a[3], max_b[3]),
           x1_inf(initcon, max_a[4], max_b[4]),
           y_inf(initcon, max_a[5], max_b[5])]


speccon = np.arange(-100, 0, 0.5)
yt = y_t(speccon, max_a[5], max_b[5])



count = 1


for i in plotvar:
    plt.plot(initcon, i, label = str(count))
    plt.xlim(-100, 50)
    plt.ylim(0, 1)
    count +=1

plt.legend(['m','h','d','f','x1','y'])