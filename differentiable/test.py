#!/usr/bin/env python3

from symulathon import Simulation
import numpy as np
import matplotlib.pyplot as plt

stiffness = 8
bend_stiffness = 0.1
tilt_angle = 0

sim = Simulation(stiffness, bend_stiffness, tilt_angle, False)
sim.fill_containers()

x = sim.getPosition()
dx0dp = sim.getInitialPositionJacobian().todense().T[3].A1
dx0dp = np.array(dx0dp)

# print(x)
print(dx0dp)
print(dx0dp.shape)

dp = 0.01
sim = Simulation(stiffness, bend_stiffness, tilt_angle + dp, False)
sim.fill_containers()

x2 = sim.getPosition()

dxdp = (x2-x)/dp
print(dxdp.shape)

plt.plot(dxdp, 'x')
plt.plot(-dx0dp, '.')
plt.legend()
plt.show()
