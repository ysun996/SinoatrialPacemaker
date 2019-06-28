from dynamics.currents import *

"""
####Base model: condition = 0
####Return currents: condition = 1
####Return voltage clamp: condition = 2


parametersdep = [4.4, 0.066, 0.5175, 0.161, 1.2, 40, 70, -93, 6, -2, 0]
parametershyp = [4.4, 0.066, 0.5175, 0.161, 1.2, 40, 70, -93, 6, 0.7, 0]
parametersbase = [4.4, 0.066, 0.5175, 0.161, 1.2, 40, 70, -93, 6, 0, 0]

tmax = 10000
dt = 0.2
t = np.arange(0, tmax, dt)
"""


def PacemakerODE(state, t, parameters, condition=0, clamp_ix=None):
    g_na = parameters[0]
    g_n_na = parameters[1]
    g_si = parameters[2]
    g_n_si = parameters[3]
    g_f_na = parameters[4]
    e_na = parameters[5]
    e_si = parameters[6]
    e_k = parameters[7]
    cm = parameters[8]
    i_bias = parameters[9]
    stimulus = parameters[10]

    v, d, f, m, h, x1, y = state
    i_s, i_n, i_x, i_k, i_y, i_bias, stimulus = current(state, t, parameters)

    i_total = i_s + i_n + i_x + i_k + i_y + i_bias + stimulus

    dv = -i_total / cm

    dd = (ad(v) * (1 - d)) - (bd(v) * d)
    df = (af(v) * (1 - f)) - (bf(v) * f)
    dm = (am(v) * (1 - m)) - (bm(v) * m)
    dh = (ah(v) * (1 - h)) - (bh(v) * h)
    dx1 = (ax1(v) * (1 - x1)) - (bx1(v) * x1)
    dy = (ay(v) * (1 - y)) - (by(v) * y)

    rstate = [dv, dd, df, dm, dh, dx1, dy]
    if clamp_ix is not None:
        rstate[clamp_ix] = 0

    retcurrent = [i_s, i_n, i_x, i_k, i_y, i_total]

    if condition != 1:
        return rstate
    else:
        return retcurrent


def current(state, t, parameters):
    g_na = parameters[0]
    g_n_na = parameters[1]
    g_si = parameters[2]
    g_n_si = parameters[3]
    g_f_na = parameters[4]
    e_na = parameters[5]
    e_si = parameters[6]
    e_k = parameters[7]
    cm = parameters[8]
    i_bias = parameters[9]
    stimulus = parameters[10]

    v, d, f, m, h, x1, y = state

    if i_bias < 0:
        i_s = 0.78 * i_si(g_si, g_n_si, d, f, v, e_si)
        i_x = 0.8 * i_x1(x1, v)
    else:
        i_s = i_si(g_si, g_n_si, d, f, v, e_si)
        i_x = i_x1(x1, v)

    i_n = i_na(g_na, m, h, g_n_na, v, e_na)
    i_ff = i_f(y, g_f_na, e_na, e_k, v)
    i_k = i_k1(v)

    return i_s, i_n, i_x, i_k, i_ff, i_bias, stimulus


def total_current(state, t, parameters):
    return np.sum(current(state, t, parameters))
