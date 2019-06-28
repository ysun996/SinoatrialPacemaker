
from solver.solver import *
import matplotlib.pyplot as plt

parameters = default_parameters()
[t, state] = solve(parameters)
plt.plot(state[:, 0], np.gradient(state[:, 0],  0.4))
plt.ylabel('V')
plt.xlabel('t')
#############################################################