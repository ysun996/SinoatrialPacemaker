import numpy as np

"""
Sodium current gating constants
"""


def am(v, max_a_m=1):
    return max_a_m * ((v + 47) / (1 - np.exp(-(v + 47) / 10)))


def bm(v, max_b_m=40):
    return max_b_m * np.exp(-(v + 72) / 17.86)


def ah(v, max_a_h=0.0085):
    return max_a_h * np.exp(-(v + 71) / 5.43)


def bh(v, max_b_h=2.5):
    return max_b_h / (1 + np.exp(-(v + 10) / 12.2))


"""
Slow inward gating constants
"""


def ad(v, max_a_d=0.005):
    return max_a_d * ((v + 34) / (1 - np.exp(-(v + 34) / 10)))


def bd(v, max_b_d=0.05):
    return max_b_d * (np.exp(-(v + 34) / 6.67))


def not_d(v):  # TODO: wat?
    nd = 1 / (1 + np.exp(-(v + 15) / 6.67))  # d' variable needed for the slow inward current
    return nd


def af(v, max_a_f=0.002468):
    return max_a_f * np.exp(-(v + 47) / 20)


def bf(v, max_b_f=0.05):
    return max_b_f / (1 + np.exp(-(v + 13) / 11.49))


"""
Constants for delayed rectifier (I_x1)
"""


def ax1(v, max_a_x1=0.0025):
    return (max_a_x1 * np.exp((v + 30) / 12.11)) / (1 + np.exp((v + 30) / 50))


def bx1(v, max_b_x1=0.0065):
    return max_b_x1 * (np.exp(-(v - 20) / 16.67)) / (1 + np.exp(-(v - 20) / 25))


def l_bar_x1(v):
    return 2.25 * ((np.exp((v + 95) / 25) - 1) / (np.exp((v + 45) / 25)))


"""
Time-independent background potassium current (I_k1)
"""


def i_k1(v):
    return 1.3 * (np.exp((v + 110) / 25) - 1) / (np.exp((v + 60) / 12.5) + np.exp((v + 60) / 25)) + 0.26 * (v + 30) / (
            1 - np.exp(-(v + 30) / 25))


"""
Constants for pacemaker current (I_y)
"""


# def ay(Vm):
#     max_a_y = 0.00005/3.0811907104922107
#     a_y = max_a_y * np.exp(-(Vm + 14.25) / 14.93)
#
#     return a_y
#
# def by(Vm):
#     max_b_y = 0.001/3.0811907104922107
#     b_y = max_b_y * (Vm + 14.25) / (1 - np.exp(-(Vm + 14.25) / 14.93))
#
#     return b_y


def ay(v, w=0.03375000000000002, max_a_y=0.00005 / 3.0811907104922107):
    # Polynomial Regression

    if isinstance(v, float):
        if v < 0:
            a = [1.4557319134e-12, 4.0945641782e-10, 4.6549818992e-08, 2.4903140216e-06, 6.1460577425e-05,
                 4.7453248494e-04, \
                 2.5019715465e-03]
            a_y = a[0] * v ** 6 + a[1] * v ** 5 + a[2] * v ** 4 + a[3] * v ** 3 + a[4] * v ** 2 + a[5] * v + a[6]
            if a_y <= 0:
                a_y = 0.00001
        else:
            a_y = max_a_y * np.exp(-(v + 14.25) / 14.93)

    else:
        a = [1.4557319134e-12, 4.0945641782e-10, 4.6549818992e-08, 2.4903140216e-06, 6.1460577425e-05,
             4.7453248494e-04, 2.5019715465e-03]
        a_y = np.zeros(v.shape)
        a_y[v < 0] = a[0] * v[v < 0] ** 6 + a[1] * v[v < 0] ** 5 + a[2] * v[v < 0] ** 4 + a[3] * v[v < 0] ** 3 + a[4] * \
                     v[v < 0] ** 2 + a[5] * v[v < 0] + a[6]
        a_y[a_y < 0] = 0.00001
        a_y[v >= 0] = max_a_y * np.exp(-(v[v >= 0] + 14.25) / 14.93)

    a_y = a_y * w
    # Modified version of the original equation
    # max_a_y = 0.00005
    # a_y = (max_a_y * np.exp(-(Vm + 14.25) / 14.93))/3
    return a_y


def by(v, w=0.03375000000000002):
    # Polynomial Regression

    if isinstance(v, float):

        a = [3.5607174324e-13, 3.9587887660e-11, -6.9345321240e-09, -8.8541673551e-07, 4.5605591007e-05,
             9.4190808268e-03, 3.3771510156e-01]
        b_y = a[0] * v ** 6 + a[1] * v ** 5 + a[2] * v ** 4 + a[3] * v ** 3 + a[4] * v ** 2 + a[5] * v + a[6]
        if b_y <= 0:
            b_y = 0.00001
    else:
        a = [3.5607174324e-13, 3.9587887660e-11, -6.9345321240e-09, -8.8541673551e-07, 4.5605591007e-05,
             9.4190808268e-03, 3.3771510156e-01]
        b_y = a[0] * v ** 6 + a[1] * v ** 5 + a[2] * v ** 4 + a[3] * v ** 3 + a[4] * v ** 2 + a[5] * v + a[6]
        for i in range(len(v)):
            if b_y[i] <= 0:
                b_y[i] = 0.00001
    b_y = b_y * w

    # Modified version of the original equation
    # max_b_y = 0.001
    # b_y = (max_b_y * (Vm + 14.25) / (1 - np.exp(-(Vm + 14.25) / 5)))/3

    return b_y


def get_a_b(current):
    a = {'m': am, 'h': ah, 'f': af, 'y': ay, 'x1': ax1, 'd': ad}
    b = {'m': bm, 'h': bh, 'f': bf, 'y': by, 'x1': bx1, 'd': bd}

    return a[current], b[current]


def x_inf(v, current):
    a, b = get_a_b(current)

    return a(v) / (a(v) + b(v))


def tau_x(v, current):
    a, b = get_a_b(current)

    return 1 / (a(v) + b(v))
