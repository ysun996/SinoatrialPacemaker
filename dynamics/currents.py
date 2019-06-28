from dynamics.gating import *


def i_si(g_si, g_n_si, d, f, v, e_si):
    si = (g_si * d * f + g_n_si * not_d(v)) * (v - e_si)
    return si


def i_na(g_na, m, h, g_n_na, v, e_na):
    na = (g_na * (m ** 3) * h + g_n_na) * (v - e_na)
    return na


def i_x1(x1, v):
    xi = x1 * l_bar_x1(v)
    return xi


def i_f(y, g_f_na, e_na, e_k, v):
    f = ((y ** 2) * g_f_na) * (v - e_na) + (v - e_k) * (y ** 2) * (-120 / e_k)
    return f


def membrane_current(g_si, g_n_si, g_na, g_n_na, g_f_na, e_si, e_na, e_k, v, d, f, m, h, x1, y, i_bias, stimulus):
    return i_si(g_si, g_n_si, d, f, v, e_si) + i_na(g_na, m, h, g_n_na, v, e_na) + i_x1(x1, v) + i_k1(v) + \
           i_f(y, g_f_na, e_na, e_k, v) + i_bias + stimulus
