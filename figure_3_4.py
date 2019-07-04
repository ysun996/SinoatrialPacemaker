import matplotlib.pyplot as plt

from solver.solver import *

parameters = default_parameters()

for pulse_amplitude in [[-0.6, -1.3, -1.9, -2.8], [0.5, 2.0, 3.5]]:
    plt.figure()
    for amp in pulse_amplitude:
        period = []
        stim_time, percent_change = prc(parameters, 30, amp, 50)
        break_point = len(percent_change)
        try:
            break_point = np.where(np.abs(np.diff(percent_change)) > 20)[0][0] + 1
            plt.plot(stim_time[break_point:], percent_change[break_point:])
        except IndexError:
            pass
        plt.plot(stim_time[:break_point], percent_change[:break_point])
