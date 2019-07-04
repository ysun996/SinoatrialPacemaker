from scipy.integrate import odeint
from scipy.signal import argrelmax

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
        raise (ValueError)

    parameters = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, i_bias, i_stim]

    return parameters


def solve(parameters, condition=0, t0=0, tmax=1500, dt=0.4, clamp_ix=None,
          state0=[-1.26312652e+01, 6.73685294e-01, 5.28447665e-01, 9.60815694e-01, 4.83891944e-07, 3.70080101e-02,
                  1.24208141e-01]):
    t = np.arange(t0, tmax, dt)
    state = odeint(lambda s, t, p: PacemakerODE(s, t, p, condition=condition, clamp_ix=clamp_ix), state0, t,
                   args=(parameters,))

    return t, state


def perturb(stim_time, parameters, amp, width, ic=[-1.26312652e+01, 6.73685294e-01, 5.28447665e-01, 9.60815694e-01, 4.83891944e-07, 3.70080101e-02,
          1.24208141e-01], tmax=1500):
    all_t = []
    all_y = []
    times = [stim_time, width, tmax - stim_time - 50]
    for ix, t_length in enumerate(times):
        if ix == 1:
            parameters[-1] = amp
        else:
            parameters[-1] = 0
        [partial_t, partial_y] = solve(parameters, dt=0.2, tmax=t_length, state0=ic)
        all_t.append(partial_t[:-1])
        all_y.append(partial_y[:-1, 0])
        ic = partial_y[-1, :]

    all_t = np.concatenate((all_t[0], all_t[1] + all_t[0][-1], all_t[2] + all_t[0][-1] + all_t[1][-1]))
    all_y = np.concatenate((all_y[0], all_y[1], all_y[2]))

    return all_t, all_y


def prc(parameters, stim_resolution, amplitude, width):
    t, state = solve(parameters, condition=0, tmax=2000)
    pks = argrelmax(state[:, 0])[0]

    base_cycle = t[pks[1]] - t[pks[0]]
    stim_time_list = np.linspace(t[pks[0]], t[pks[1]], stim_resolution)
    period = []

    for stim_time in stim_time_list:
        t, y = perturb(stim_time, parameters, amplitude, width)

        pks = argrelmax(y, order=10)[0]
        pks = pks[y[pks] > 0]
        period.append(t[pks[1]] - t[pks[0]])

    percent_change = 100 * (period - base_cycle) / base_cycle
    return stim_time_list, percent_change
