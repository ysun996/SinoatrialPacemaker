import numpy as np
from scipy.integrate import odeint
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt

#comment

#####Sodium current gating constants
def am(Vm):
    max_a_m = 1

    a_m = max_a_m * ((Vm + 47) / (1 - np.exp(-(Vm + 47) / 10)))
    return a_m

def bm(Vm):
    max_b_m = 40
    b_m = max_b_m * np.exp(-(Vm + 72) / 17.86)


    return b_m

def ah(Vm):
    max_a_h = 0.0085
    a_h = max_a_h * np.exp(-(Vm + 71) / 5.43)


    return a_h

def bh(Vm):
    max_b_h = 2.5
    b_h = max_b_h / (1 + np.exp(-(Vm + 10) / 12.2))

    return b_h

#####Slow inward gating constants

def ad(Vm):
    max_a_d = 0.005
    a_d = max_a_d * ((Vm + 34) / (1 - np.exp(-(Vm + 34) / 10)))
    return a_d

def bd(Vm):
    max_b_d = 0.05
    b_d = max_b_d * (np.exp(-(Vm + 34) / 6.67))

    return b_d

def not_d(Vm):
    ### d' variable needed for the slow inward current
    not_d = 1 / (1 + np.exp(-(Vm + 15) / 6.67))

    return not_d

def af(Vm):
    max_a_f = 0.002468
    a_f = max_a_f * np.exp(-(Vm + 47) / 20)

    return a_f

def bf(Vm):
    max_b_f = 0.05
    b_f = max_b_f / (1 + np.exp(-(Vm + 13) / 11.49))


    return b_f

##### Constants for delayed rectifier (I_x1)

def ax1(Vm):
    max_a_x1 = 0.0025
    a_x1 = (max_a_x1 * np.exp((Vm + 30) / 12.11)) / (1 + np.exp((Vm + 30) / 50))


    return a_x1

def bx1(Vm):
    max_b_x1 = 0.0065
    b_x1 = max_b_x1 * (np.exp(-(Vm - 20) / 16.67)) / (1 + np.exp(-(Vm - 20) / 25))

    return b_x1

def l_bar_x1(Vm):
    l_bar = 2.25 * ((np.exp((Vm + 95) / 25) - 1) / (np.exp((Vm + 45) / 25)))

    return l_bar

#####Constants for pacemaker current (I_y)

def ay(Vm):
    max_a_y = 0.00005/3.0811907104922107
    a_y = max_a_y * np.exp(-(Vm + 14.25) / 14.93)

    return a_y

def by(Vm):
    max_b_y = 0.001/3.0811907104922107
    b_y = max_b_y * (Vm + 14.25) / (1 - np.exp(-(Vm + 14.25) / 14.93))

    return b_y


#####Time-independent background potassium current (I_k1)

def i_k1(Vm):
    IPot = 1.3 * (np.exp((Vm + 110) / 25) - 1) / (np.exp((Vm + 60) / 12.5) + np.exp((Vm + 60) / 25)) + \
           0.26 * (Vm + 30) / (1 - np.exp(-(Vm + 30) / 25))
    return IPot

#Current functions
def Isi(gSi, gnSi, d, f, v, ESi):
    si = (gSi * d * f + gnSi * not_d(v)) * (v - ESi)
    return si

def Ina(gNa, m, h, gnNa, v, ENa):
    na = (gNa * (m**3) * h + gnNa) * (v - ENa)
    return na

def Ix1(x1, v):
    xi = x1 * l_bar_x1(v)
    return xi

def If(y, gfNa, ENa, EK, v):
    f = ((y**2) * gfNa) * (v - ENa) + (v-EK) * (y**2) * (-120/EK)
    return f

def Iall(gSi, gnSi, gNa, gnNa, gfNa, ESi, ENa, EK, v, d, f, m, h, x1, y, Ibias, Istim):
    all = Isi(gSi, gnSi, d, f, v, ESi) + Ina(gNa, m, h, gnNa, v, ENa) + Ix1(x1, v) + i_k1(v) + If(y, gfNa, ENa, EK, v) + \
         Ibias + Istim
    return all

####Base model: condition = 0
####Return currents: condition = 1
####Return voltage clamp: condition = 2

def PacemakerODE(state, t, parameters, condition):
    gNa = parameters[0]
    gnNa = parameters[1]
    gSi = parameters[2]
    gnSi = parameters[3]
    gfNa = parameters[4]
    ENa = parameters[5]
    ESi = parameters[6]
    EK = parameters[7]
    Cm = parameters[8]
    Ibias = parameters[9]
    Istim = parameters[10]

    v = state[0]
    d = state[1]
    f = state[2]
    m = state[3]
    h = state[4]
    x1 = state[5]
    y = state[6]

    if Ibias < 0:
        Is = 0.78 * Isi(gSi, gnSi, d, f, v, ESi)
        Ix = 0.8 * Ix1(x1, v)
    else:
        Is = Isi(gSi, gnSi, d, f, v, ESi)
        Ix = Ix1(x1, v)

    In = Ina(gNa, m, h, gnNa, v, ENa)
    Ik = i_k1(v)
    Iy = If(y, gfNa, ENa, EK, v)

    Ii = Is + In + Ix + Ik + Iy + Ibias + Istim

    if condition <= 1:
        dv = -Ii/Cm

    else:
        dv = 0

    dd = (ad(v) * (1-d)) - (bd(v) * d)
    df = (af(v) * (1-f)) - (bf(v) * f)

    dm = (am(v) * (1-m)) - (bm(v) * m)
    dh = (ah(v) * (1-h)) - (bh(v) * h)

    dx1 = (ax1(v) * (1-x1)) - (bx1(v) * x1)

    dy = (ay(v) * (1-y)) - (by(v) * y)

    rstate = [dv, dd, df, dm, dh, dx1, dy]
    retcurrent = [Is, In, Ix, Ik, Iy, Ii]

    if condition != 1:
        return rstate
    else:
        return retcurrent
#Values
##Parameters
parametersdep = [4.4, 0.066, 0.5175, 0.161, 1.2, 40, 70, -93, 6, -2, 0]
parametershyp = [4.4, 0.066, 0.5175, 0.161, 1.2, 40, 70, -93, 6, 0.7, 0]
parametersbase = [4.4, 0.066, 0.5175, 0.161, 1.2, 40, 70, -93, 6, 0, 0]

##Time
tmax = 10000
dt = 0.2
t = np.arange(0, tmax, dt)
