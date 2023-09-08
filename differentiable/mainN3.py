#!/usr/bin/env python3
"""Computes plot of 2 interpolated parameters"""

from symulathon import Simulation
from simulation_functions import simulate
import numpy as np
from interpolate_parameters import generate_parameters
from tqdm import tqdm
import matplotlib.pyplot as plt


if __name__ == "__main__":
    # Initialize important constant variables
    sim = Simulation(1, 1)
    nDoF = sim.getDoF()
    mass = sim.getMassMatrix()
    h = sim.getTimeStep()
    DIFF_FRAMES = 100

    # Define the region studied
    n_points = 11
    # k_values = np.linspace(3.90, 4.1, n_points)
    # k_bend_values = np.linspace(0, 0.5, n_points)
    # print(k_values)
    # print(k_bend_values)
    # X, Y = np.meshgrid(k_values, k_bend_values)
    k_valuesX = np.linspace(3.95, 4.05, n_points)
    k_valuesY = np.linspace(3.95, 4.05, n_points)
    kbend4 = 0.1 * np.ones(4)
    X, Y = np.meshgrid(k_valuesX, k_valuesY)

    g_values = X.tolist()

    # Computation of the surface and gradients
    # for i in tqdm(range(len(k_values))):
    #     for j in range(len(k_bend_values)):
    #         k4 = np.ones(4) * k_values[i]
    #         kbend4 = np.ones(4) * k_bend_values[j]
    #         # k4 = np.array([k_valuesX[i], k_valuesX[i],
    #         #                k_valuesY[j], k_valuesY[j]])
    #         k, kb = generate_parameters(k4, kbend4)
    #         bp = simulate(k, kb, DIFF_FRAMES)
    #         g_values[j][i] = bp.get_g()

    for i in tqdm(range(len(k_valuesX))):
        for j in range(len(k_valuesY)):
            k4 = np.array([k_valuesX[i], k_valuesX[j],
                           k_valuesY[i], k_valuesY[j]])
            # k4 = np.array([k_valuesX[i], k_valuesX[i],
            #                k_valuesY[j], k_valuesY[j]])
            k, kb = generate_parameters(k4, kbend4)
            bp = simulate(k, kb, DIFF_FRAMES)
            g_values[j][i] = bp.get_g()

    # Plotting
    # 3D plot
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    g_values = np.array(g_values)

    surf = ax.plot_surface(X, Y, g_values)
    ax.set(
        xlabel="k values",
        ylabel="k bend values",
    )

    plt.show()
