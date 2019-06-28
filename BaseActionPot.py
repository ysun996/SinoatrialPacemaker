
from solver.solver import *
import matplotlib.pyplot as plt

parameters = default_parameters()
[t, state] = solve(parameters)
plt.plot(t, state[:, 0])
plt.ylabel('V')
plt.xlabel('t')
plt.xlim(250, 1250)
#############################################################