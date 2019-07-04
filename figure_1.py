import matplotlib.pyplot as plt
import numpy as np

from dynamics.gating import x_inf, tau_x, ay, by

plt.figure()

"""
x_inf value
"""

plt.subplot(211)

v = np.linspace(-100, 50, 100)
current_names = ['m', 'h', 'f', 'y', 'x1', 'd']
for current in current_names:
    plt.plot(v, x_inf(v, current))

plt.legend(current_names)

"""
I_f gating
"""

plt.subplot(212)

v = np.linspace(-100, 0, 100)
plt.plot(v, x_inf(v, 'y'))
plt.plot(v, 0.001 * tau_x(v, 'y'))  # change to seconds
plt.plot(v, 25 * ay(v))  # extra scale
plt.plot(v, 25 * by(v))  # extra scale
