import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

###Sodium gating

def am(v):
    amax = 1
    a_m = amax * ((v + 47)/(1 - np.exp(-(v + 47)/10)))

    return a_m

def bm(v):
    bmax = 40
    b_m = bmax * np.exp(-(v + 72)/17.86)

    return b_m

def ah(v):
    amax = 0.0085
    a_h = amax * np.exp(-(v + 71)/5.43)

    return a_h

def bh(v):
    bmax = 2.5
    b_h = bmax/(1 + np.exp(-(v + 10)/12))

    return b_h

###I_si

def ad(v):
    amax = 0.005
    a_d = amax * ((v + 34)/(1-np.exp((v + 34)/10)))

    return a_d

def bd(v):
    bmax = 0.05
    b_d = bmax * np.exp(-(v + 34)/6.67)

    return b_d

def not_d(v):
    notd = 1/(1+ np.exp(-(v+15)/6.67))

    return notd

def af(v):
    amax = 0.002468
    a_f = amax * np.exp(-(v + 47))

    return a_f

def bf(v):
    bmax = 0.05
    b_f = bmax/(1 + np.exp(-(v + 13)/11.49))

    return b_f

###I_xi

def ax1(v):
    amax = 0.0025
    a_x1 = (amax * np.exp(v + 30)/12.11)/(1 + np.exp((v + 30)/50))

    return a_x1

def bx1(v):
    bmax = 0.0065
    b_x1 = bmax * (np.exp(-(v - 20)/16.67))/(1 + np.exp(-(v - 20)/25))

    return b_x1

def lx1max(v):
    lx1 = 2.25 * ((np.exp((v+95)/25) - 1)/np.exp((v + 45)/25))

    return lx1

###I_y

def ay(v):
    amax = 0.0000167
    a_y = amax * np.exp((v + 32)/14.93)

    return a_y

def by(v):
    bmax = 0.00033
    b_y = bmax * (v + 32) / (1 - np.exp(-(v + 32) / 5))

    return b_y

def I_k1(Vm):
    # ik1pt1 = 1.3 * ((np.exp((v + 110)/25)-1))/(np.exp((v + 60)/12.5) + np.exp((v + 60)/25))
    # ik1pt2 = 0.26 * (v + 30)/(1 - np.exp(-(v + 30)/25))
    # ik = ik1pt1 + ik1pt2
    ik = 1.3*(np.exp((Vm+110)/25)-1)/(np.exp((Vm+60)/12.5)+np.exp((Vm+60)/25))+\
         0.26*(Vm+30)/(1-np.exp(-(Vm+30)/25))

    return ik


########Parameters

gNa = 4.4
gnNa = 0.066
gSi = 0.5175
gnSi = 0.161
gfNa = 1.2
Cm = 6
ENa = 40
ESi = 70
EK = -93
I = 0
# I = np.arange(20, 42, 2)
parameters = [gNa, gnNa, gSi, gnSi, gfNa, Cm, ENa, ESi, EK, I]

########Pacemaker ODE

def pacemakerODE(state, t, parameters):

    gNa = parameters[0]
    gnNa = parameters[1]
    gSi = parameters[2]
    gnSi = parameters[3]
    gfNa = parameters[4]
    Cm = parameters[5]
    ENa = parameters[6]
    ESi = parameters[7]
    EK = parameters[8]
    I = parameters[9]

    v = state[0]
    m = state[1]
    h = state[2]
    d = state[3]
    f = state[4]
    x1 = state[5]
    y = state[6]

    dno = not_d(v)
    gna = gNa * (m**3) * h + gnNa
    gsi = gSi * d * f + gnSi * dno
    gfna = (y **2) * gfNa
    gfk = (-120/EK) * (y **2)

    Ina = gna * (v - ENa)
    Isi = gsi * (v - ESi)
    Ix1 = x1 * lx1max(v)
    Ik1 = I_k1(v)
    If = gfna * (v - ENa) + gfk * (v - EK)

    Ii = I + Ina + Isi + Ix1 + Ik1 + If

    dv = -Ii/Cm
    dm = am(v) * (1 - m) - bm(v) * m
    dh = ah(v) * (1 - h) - bh(v) * h

    dd = ad(v) * (1 - d) - bd(v) * d
    df = af(v) * (1 - f) - bf(v) * f

    dx1 = ax1(v) * (1 - x1) - bx1(v) * x1

    dy = ay(v) * (1 - y) - by(v) * y

    rstate = [dv, dm, dh, dd, df, dx1, dy]

    return rstate

tmax = 2000
dt = 0.2
t = np.arange(0, tmax, dt)

biascurrent = np.arange(-70, -40, 1)

openstate = np.arange(0, 1.1, 0.1)
state0 = [-40, 0, 0, 0, 0, 0, 0]
# for i in biascurrent:
#     state0 = [i, 0, 1, 0, 1, 0, 1]
#     state = odeint(pacemakerODE, state0, t, args = (parameters,))
#     plt.plot(t, state[:,0], label = str(i))

# for i in I:
#     parameters = [gNa, gnNa, gSi, gnSi, gfNa, Cm, ENa, ESi, EK, i]
#     state0 = [-45, 0, 1, 0, 1, 0, 1]
#     state = odeint(pacemakerODE, state0, t, args=(parameters,))
#     plt.plot(t, state[:, 0], label=str(i))
# state0 = [-10, 0, 0, 0, 0, 0, 0]
state = odeint(pacemakerODE, state0, t, args=(parameters,))
plt.plot(t, state[:, 0])
plt.legend()
