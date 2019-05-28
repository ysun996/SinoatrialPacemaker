import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.fftpack import fft
from scipy.signal import argrelextrema
from common import *


#Parameters
gNa = 4.4
gnNa = 0.066
gSi = 0.5175
gnSi = 0.161
gfNa = 1.2
ENa = 40
ESi = 70
EK = -93
Cm = 6
I = 0
parameters = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, I]

#Time
tmax = 4000
dt = 0.2
t = np.arange(400, tmax, dt)

def PacemakerODE(state, t, parameters):
    gNa = parameters[0]
    gnNa = parameters[1]
    gSi = parameters[2]
    gnSi = parameters[3]
    gfNa = parameters[4]
    ENa = parameters[5]
    ESi = parameters[6]
    EK = parameters[7]
    Cm = parameters[8]
    I = parameters[9]

    v = state[0]
    d = state[1]
    f = state[2]
    m = state[3]
    h = state[4]
    x1 = state[5]
    y = state[6]

    Isi = (gSi * d * f + gnSi * not_d(v)) * (v - ESi)
    Ina = (gNa * (m**3) * h + gnNa) * (v - ENa)
    Ix1 = x1 * l_bar_x1(v)
    Ik1 = i_k1(v)
    If = ((y**2) * gfNa) * (v - ENa) + (v-EK) * (y**2) * (-120/EK)

    Ii = Isi + Ina + Ix1 + Ik1 + If + I

    dv = -Ii/Cm

    dd = (ad(v) * (1-d)) - (bd(v) * d)
    df = (af(v) * (1-f)) - (bf(v) * f)

    dm = (am(v) * (1-m)) - (bm(v) * m)
    dh = (ah(v) * (1-h)) - (bh(v) * h)

    dx1 = (ax1(v) * (1-x1)) - (bx1(v) * x1)
    #print(dx1)

    dy = (ay(v) * (1-y)) - (by(v) * y)

    rstate = [dv, dd, df, dm, dh, dx1, dy]

    return rstate

state0 = [-66, 0, 0, 0, 0, 0, 0]

state = odeint(PacemakerODE, state0, t, args=(parameters,))
dv = []
for points in state:
    pacemaker = PacemakerODE(points, t, parameters)
    dv.append(pacemaker[0])
# state1 = state[190]
# plt.plot(state[:,0], state[:,0]/t)
plt.plot(t, state[:,0])
# plt.figure()
#plt.plot(state[:,0], dv)

#######Fast Fourier transform
volts = np.zeros(len(state[:,0]))
for ind in range(len(volts)):
    volts[ind] = state[ind, 0]
CL = argrelextrema(volts, np.greater) ####Cycle length array, numbers refer to indices
timpoints = CL[0] * 0.2








