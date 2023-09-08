#!/usr/bin/env python3
"""Computes finite derivative"""
from symulathon import Simulation, count_springs
from simulation_functions import simulate
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    # Initialize important constant variables
    DIFF_FRAMES = 10
    sim = Simulation(1, 1)
    nFlex, nBend = count_springs()
    print(f"nFlex = {nFlex}, nBend = {nBend}, total_diff = {nFlex + nBend}")
    nDoF = sim.getDoF()
    mass = sim.getMassMatrix()
    h = sim.getTimeStep()

    # Define finite der
    k_list = np.ones(nFlex) * 1
    k_bend_list = np.ones(nBend) * 0.1
    parameters = np.concatenate((k_list, k_bend_list))
    dgdp_finite = np.zeros(len(parameters))
    bp = simulate(k_list, k_bend_list, DIFF_FRAMES)

    # Calculate gradient and loss at point
    g0 = bp.get_g()
    dgdp = bp.get_dgdp()

    DeltaK = 0.01
    for i in tqdm(range(len(parameters))):
        dk_vec = np.zeros(len(parameters))
        dk_vec[i] = DeltaK
        dparameters = parameters + dk_vec
        dk_list, dk_bend_list = (dparameters[:nFlex], dparameters[nFlex:])
        bp = simulate(dk_list, dk_bend_list, DIFF_FRAMES)
        gi = bp.get_g()
        dgdp_finite[i] = (gi - g0) / DeltaK

    error = (dgdp - dgdp_finite)**2
    print("Backpropagation:", dgdp)
    print("Finite diff:", dgdp_finite)
    print("Error term:", error)

    plt.plot(dgdp, marker='x', label="backpropagation")
    plt.plot(dgdp_finite, marker='+', label="finite differences")
    plt.plot(error, label="Error")
    plt.legend()
    plt.show()
