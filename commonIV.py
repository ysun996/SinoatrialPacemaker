from common import *

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

    if Ibias < 0:
        Ii = 0.78 * Isi(gSi, gnSi, d, f, v, ESi) + Ina(gNa, m, h, gnNa, v, ENa) + 0.8 * Ix1(x1, v) + i_k1(v) + \
             If(y, gfNa, ENa, EK, v) + Ibias + Istim
    else:
        Ii = Isi(gSi, gnSi, d, f, v, ESi) + Ina(gNa, m, h, gnNa, v, ENa) + Ix1(x1, v) + i_k1(v) + \
             If(y, gfNa, ENa, EK, v) + Ibias + Istim

    return Ii

#Values:
##voltage clamp
vclamp = np.arange(-60, -10, 0.1)
condition = 2