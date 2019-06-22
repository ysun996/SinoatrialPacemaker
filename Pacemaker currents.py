import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.fftpack import fft
from scipy.signal import argrelextrema
from common import *

condition = 0

state0 = [-66, 0, 0, 0, 0, 0, 0]

state = odeint(PacemakerODE, state0, t, args=(parametersbase,))
dv = []
for points in state:
    pacemaker = PacemakerODE(points, t, parametersbase, condition)
    dv.append(pacemaker[0])
plt.plot(t, state[:,0])
plt.figure()
plt.plot(state[:,0], dv)









