from solver.solver import *
import matplotlib.pyplot as plt
############################################################
# Currents

parameters = default_parameters()
t, state = solve(parameters)
currents = np.array(list(map(lambda x: current(x, t, parameters), state)))

Ilabels = ['Isi', 'Ina', 'Ix1', 'Ik1', 'If', 'Ii']
for w in range(0, 5):
    plt.plot(t, currents[:, w], label=Ilabels[w])

plt.legend()
plt.ylabel('Current')
plt.xlabel('t')
plt.xlim(250, 1250)

############################################################
