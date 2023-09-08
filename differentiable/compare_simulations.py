#!/usr/bin/env python3
from symulathon import Simulation
from recorder import SimulationReader
from sys import argv
from meansquareloss import MeanSquareLoss
import numpy as np

if (len(argv) != 3):
    print("It is necessary to pass 2 files!")
    exit()

nDoF = Simulation(1, 1).getDoF()

file1 = argv[1]
file2 = argv[2]
reader1 = SimulationReader(nDoF, file1)
reader2 = SimulationReader(nDoF, file2)
N_FRAMES = 400
loss = MeanSquareLoss()
g = 0
for i in range(N_FRAMES):
    x1, v1 = reader1.get_next_state()
    x2, v2 = reader2.get_next_state()
    loss.update_containers(x1, v1, x2, v2)
    g += loss.evaluate()

print(f"Loss g = {g}")
