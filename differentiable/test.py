#!/usr/bin/env python3

from solve_system import solve_system
from symulathon import Simulation
from simulation_functions import simulate, nFlex, nBend
import numpy as np

k_value = 40
k_bend_value = 0.1

bp = simulate(k_value, k_bend_value, 100)
print(bp.get_dgdp())
print(f'Single parameter g =\t {bp.get_g()}')

tension_springs = [k_value]*nFlex
bending_springs = [k_bend_value]*nBend
bp = simulate(tension_springs, bending_springs, 100)
base_g = bp.get_g()
print(f'Multiple parameters g =\t {base_g}')
dgdp = bp.get_dgdp()

index = 33
dk = 0.01
tension_springs[index] += dk
bp = simulate(tension_springs, bending_springs, 100)
new_g = bp.get_g()
dgdp2 = bp.get_dgdp()

print(dgdp[index])
print(dgdp2[index])
print((base_g - new_g) / dk)
