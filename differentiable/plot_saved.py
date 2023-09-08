#!/usr/bin/env python3

from recorder import SimulationReader
import matplotlib.pyplot as plt
import numpy as np

nDoF = 6
reader = SimulationReader(nDoF)

history = reader.get_all_history()

history = np.array(history)
plt.plot(history[3], label="x2")
plt.plot(history[4], label="y2")
plt.plot(history[5], label="z2")
plt.plot(history[9], label="vx2")
plt.plot(history[10], label="vy2")
plt.plot(history[11], label="vz2")
plt.legend()
plt.show()

plt.plot(history.T, label="xd")
plt.legend()
plt.show()
