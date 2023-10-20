#!/usr/bin/env python3

from symulathon import Simulation
import numpy as np

sim = Simulation()
ndof = sim.getDoF()

initial_velocities = np.zeros(ndof)

sim.reset_simulation(initial_velocities)
