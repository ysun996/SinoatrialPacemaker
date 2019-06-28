from scipy.integrate import odeint

from dynamics.ode import *


def default_parameters(state='base'):
    gNa = 4.4
    gnNa = 0.066
    gSi = 0.5175
    gnSi = 0.161
    gfNa = 1.2
    ENa = 40
    ESi = 70
    EK = -93
    Cm = 6
    i_stim = 0

    if state == 'base':
        i_bias = 0
    elif state == 'dep':
        i_bias = -2
    elif state == 'hyp':
        i_bias = 0.7
    else:
        raise(ValueError)

    parameters = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, i_bias, i_stim]

    return parameters


def solve(parameters, condition=0, tmax=1500, dt=0.4, clamp_ix=None,
          state0=[-1.26312652e+01, 6.73685294e-01, 5.28447665e-01, 9.60815694e-01, 4.83891944e-07, 3.70080101e-02,
                  1.24208141e-01]):
    t = np.arange(0, tmax, dt)
    state = odeint(lambda s, t, p: PacemakerODE(s, t, p, condition=condition, clamp_ix=clamp_ix), state0, t,
                   args=(parameters,))

    return t, state
