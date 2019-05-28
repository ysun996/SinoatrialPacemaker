import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.fftpack import fft
from scipy.signal import argrelextrema
from common import *
from commonIV import *

parametersd = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, -2, Istim]

def memcurrent(state, parameters):
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

    Isi = 0.78 * (gSi * d * f + gnSi * not_d(v)) * (v - ESi)
    Ina = (gNa * (m ** 3) * h + gnNa) * (v - ENa)
    Ix1 = 0.8 * x1 * l_bar_x1(v)
    Ik1 = i_k1(v)
    If = ((y ** 2) * gfNa) * (v - ENa) + (v - EK) * (y ** 2) * (-120 / EK)

    Ii = Isi + Ina + Ix1 + Ik1 + If + Ibias + Istim
    memcurrent = Ii

    return memcurrent

currentarrayd = []

for v in vclamp:
    state0 = [v, 0, 0, 0, 0, 0, 0]
    state = odeint(PacemakerODE, state0, t, args=(parametersd,))
    currentarrayd.append(memcurrent(state[-1], parametersd))

plt.plot(vclamp, currentarrayd)
plt.xlabel('Voltage (mV)')
plt.ylabel('Current (uA/cm^2)')
plt.title('Depressed model')