#!/usr/bin/env python3
"""Plot loss landscape for 2 differentiable parameters."""

from symulathon import Simulation
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from tqdm import tqdm
from simulation_functions import simulate
from concurrent import futures

if __name__ == "__main__":
    # Initialize important constant variables
    sim = Simulation(1, 1)
    nDoF = sim.getDoF()
    mass = sim.getMassMatrix()
    h = sim.getTimeStep()
    DIFF_FRAMES = 100

    # Define the region studied
    n_points = 21
    k_values = np.linspace(3.5, 4, n_points)
    # k_values = np.linspace(6, 9, n_points)
    k_bend_values = np.linspace(2, 2.5, n_points-1)
    # k_values = np.linspace(3.95, 4.05, n_points)
    # k_bend_values = np.linspace(0, 0.2, n_points)
    # k_values = np.linspace(3.45, 3.55, n_points)
    # k_bend_values = np.linspace(0.35, 0.4, n_points)
    X, Y = np.meshgrid(k_values, k_bend_values)

    g_values = X.tolist()
    dgdk_values = X.tolist()
    dgdk_bend_values = X.tolist()
    ones = np.zeros(X.shape)

    # Computation of the surface and gradients
    # for i in tqdm(range(len(k_values))):
    #     for j in range(len(k_bend_values)):
    #         bp = simulate(k_values[i], k_bend_values[j], DIFF_FRAMES)
    #         dgdp = bp.get_dgdp()
    #         dgdk_values[j][i] = dgdp[0]
    #         dgdk_bend_values[j][i] = dgdp[1]
    #         g_values[j][i] = bp.get_g()

    pool = futures.ProcessPoolExecutor()
    for i in tqdm(range(len(k_values))):
        k_values_i = [k_values[i]] * len(k_bend_values)
        DIFF_FRAMES_i = [DIFF_FRAMES] * len(k_bend_values)
        bps = pool.map(simulate, k_values_i, k_bend_values, DIFF_FRAMES_i)
        for j, bp in enumerate(bps):
            dgdp = bp.get_dgdp()
            dgdk_values[j][i] = dgdp[0]
            dgdk_bend_values[j][i] = dgdp[1]
            g_values[j][i] = bp.get_g()
    # Finite differences
    dk = k_values[1]-k_values[0]
    dk_bend = k_bend_values[1]-k_bend_values[0]
    dgdk_bend_values_finite, dgdk_values_finite = np.gradient(g_values, dk_bend, dk)

    # Plotting
    # 3D plot
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    g_values = np.array(g_values)
    dgdk_values = np.array(dgdk_values)
    dgdk_bend_values = np.array(dgdk_bend_values)
    diff = np.sqrt(
        (dgdk_values - dgdk_values_finite)**2 +
        (dgdk_bend_values - dgdk_bend_values_finite)**2)

    magnitudes = np.sqrt(dgdk_bend_values**2 + dgdk_values**2)

    diff_rel = diff / magnitudes

    colors = cm.Greys(diff)
    # surf = ax.plot_surface(X, Y, g_values, facecolors=colors)
    surf = ax.plot_surface(X, Y, g_values)

    gradient = ax.quiver(X, Y, g_values, dgdk_values, dgdk_bend_values, ones,
                         arrow_length_ratio=0.001,
                         length=0.5,
                         normalize=True)

    gradient_finite = ax.quiver(X, Y, g_values,
                                dgdk_values_finite,
                                dgdk_bend_values_finite,
                                ones,
                                arrow_length_ratio=0.001,
                                color="red",
                                length=0.5,
                                normalize=True)
    ax.set(
        xlabel="k values",
        ylabel="k bend values",
        zlabel="loss",
    )
    # Add a color bar which maps values to colors.
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

    plt.title("Real magnitudes")
    plt.contourf(X, Y, g_values, levels=1000)
    plt.colorbar()
    plt.quiver(X, Y, dgdk_values, dgdk_bend_values,
               label="Backpropagation", color="blue", angles='xy')
    plt.quiver(X, Y, dgdk_values_finite, dgdk_bend_values_finite,
               label="Finite diff", color="red", angles='xy')
    plt.legend()
    plt.show()

    plt.title("Real magnitudes")
    plt.contourf(X, Y, g_values, levels=1000)
    plt.colorbar()
    plt.legend()
    plt.show()

    plt.title("Normalized magnitudes")
    plt.contourf(X, Y, g_values, levels=200)
    plt.colorbar()
    plt.quiver(X, Y, dgdk_values/magnitudes, dgdk_bend_values/magnitudes,
               label="Backpropagation", color="blue", angles='xy')
    plt.quiver(X, Y, dgdk_values_finite/magnitudes, dgdk_bend_values_finite/magnitudes,
               label="Finite diff", color="red", angles='xy')
    plt.legend()
    plt.show()

    plt.title("Relative differences")
    plt.contourf(X, Y, diff_rel, levels=200)
    plt.colorbar()
    plt.quiver(X, Y, dgdk_values/magnitudes, dgdk_bend_values/magnitudes,
               label="Backpropagation", color="blue", angles='xy')
    plt.quiver(X, Y, dgdk_values_finite/magnitudes, dgdk_bend_values_finite/magnitudes,
               label="Finite diff", color="red", angles='xy')
    plt.legend()
    plt.show()

    # save the data to disk
    np.save("data/x_data.npy", X)
    np.save("data/y_data.npy", Y)
    np.save("data/g_data.npy", g_values)
    np.savez("data/dgdk_data.npy",
             dgdk_values,
             dgdk_bend_values,
             dgdk_values_finite,
             dgdk_bend_values_finite)
