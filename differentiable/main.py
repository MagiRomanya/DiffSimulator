#!/usr/bin/env python3

from symulathon import Simulation
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from simulation_functions import simulate

plt.rcParams['text.usetex'] = True
plt.rcParams["savefig.bbox"] = 'tight'
plt.rcParams["savefig.dpi"] = 96*2
plt.rcParams["savefig.format"] = 'pdf'
plt.rcParams["figure.figsize"] = (20,3)

if __name__ == "__main__":
    sim = Simulation(1, 1)
    nDoF = sim.getDoF()
    mass = sim.getMassMatrix()
    h = sim.getTimeStep()
    DIFF_FRAMES = 100

    k_values = np.linspace(0.01, 40, 3)
    g_values = []
    dgdp_values = []
    for k in tqdm(k_values):
        bp = simulate(k, 0.1, DIFF_FRAMES)
        g_values.append(bp.get_g())
        dgdp_values.append(bp.get_dgdp()[0])

    # Calculate finite differences
    dgdp_finite = []
    for i in range(len(g_values)-1):
        value = (g_values[i+1]-g_values[i]) / (k_values[i+1] - k_values[i])
        dgdp_finite.append(value)
    dgdp_finite.append(dgdp_finite[-1])

    plt.plot(k_values, g_values, "-", label="Loss Function")
    plt.plot(k_values, dgdp_finite, ".", label="Finite $\\frac{\mathrm{d}g}{\mathrm{d}p}$")
    plt.plot(k_values, dgdp_values, "x", label="Back propagation $\\frac{\mathrm{d}g}{\mathrm{d}p}$")

    plt.legend()
    plt.xlabel("$k_s$")
    plt.grid()
    plt.show()
    print("Finished succesfully")
