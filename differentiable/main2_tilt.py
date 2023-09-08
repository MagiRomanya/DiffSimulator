#!/usr/bin/env python3

""" Plot using 2 simulation parameters one of them being tilt."""
from symulathon import Simulation
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from tqdm import tqdm
from simulation_functions import simulate_tilt, simulate
from concurrent import futures

if __name__ == "__main__":
    # Initialize important constant variables
    sim = Simulation(1, 1)
    nDoF = sim.getDoF()
    mass = sim.getMassMatrix()
    h = sim.getTimeStep()
    DIFF_FRAMES = 100

    # Define the region studied
    n_points = 41
    k_values = np.linspace(6, 10, n_points)
    tilt_angle_values = np.linspace(20, -20, n_points-1)

    X, Y = np.meshgrid(k_values, tilt_angle_values)

    g_values = X.tolist()
    dgdk_values = X.tolist()
    dgdtilt_values = X.tolist()
    ones = np.zeros(X.shape)

    pool = futures.ProcessPoolExecutor()
    for i in tqdm(range(len(k_values))):
        k_values_i = [k_values[i]] * len(tilt_angle_values)
        DIFF_FRAMES_i = [DIFF_FRAMES] * len(tilt_angle_values)
        bps = pool.map(simulate_tilt, k_values_i, tilt_angle_values, DIFF_FRAMES_i)
        for j, bp in enumerate(bps):
            dgdp = bp.get_dgdp()
            dgdk_values[j][i] = dgdp[0]     # ( k value )
            dgdtilt_values[j][i] = dgdp[2]  # ( tilt angle value )
            g_values[j][i] = bp.get_g()

    # Finite differences
    dk = k_values[1]-k_values[0]
    dtilt = tilt_angle_values[1] - tilt_angle_values[0]
    dgdtilt_values_finite, dgdk_values_finite = np.gradient(g_values, dtilt, dk)

    # Plotting
    # 3D plot
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    g_values = np.array(g_values)
    dgdk_values = np.array(dgdk_values)
    dgdtilt_values = np.array(dgdtilt_values)
    diff = np.sqrt(
        (dgdk_values - dgdk_values_finite)**2 +
        (dgdtilt_values - dgdtilt_values_finite)**2)

    magnitudes = np.sqrt(dgdtilt_values**2 + dgdk_values**2)

    diff_rel = diff / magnitudes

    colors = cm.Greys(diff)
    # surf = ax.plot_surface(X, Y, g_values, facecolors=colors)
    surf = ax.plot_surface(X, Y, g_values)

    gradient = ax.quiver(X, Y, g_values, dgdk_values, dgdtilt_values, ones,
                         arrow_length_ratio=0.001,
                         length=0.5,
                         normalize=True)

    gradient_finite = ax.quiver(X, Y, g_values,
                                dgdk_values_finite,
                                dgdtilt_values_finite,
                                ones,
                                arrow_length_ratio=0.001,
                                color="red",
                                length=0.5,
                                normalize=True)
    ax.set(
        xlabel="k values",
        ylabel="tilt angle values",
        zlabel="loss",
    )
    # Add a color bar which maps values to colors.
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

    plt.title("Real magnitudes")
    plt.contourf(X, Y, g_values, levels=40)
    plt.colorbar()
    plt.quiver(X, Y, dgdk_values, dgdtilt_values,
               label="Backpropagation", color="blue", angles='xy')
    plt.quiver(X, Y, dgdk_values_finite, dgdtilt_values_finite,
               label="Finite diff", color="red", angles='xy')
    plt.legend()
    plt.show()

    plt.title("Normalized magnitudes")
    plt.contourf(X, Y, g_values, levels=200)
    plt.colorbar()
    plt.quiver(X, Y, dgdk_values/magnitudes, dgdtilt_values/magnitudes,
               label="Backpropagation", color="blue", angles='xy')
    plt.quiver(X, Y, dgdk_values_finite/magnitudes, dgdtilt_values_finite/magnitudes,
               label="Finite diff", color="red", angles='xy')
    plt.legend()
    plt.show()

    plt.title("Relative differences")
    plt.contourf(X, Y, diff_rel, levels=200)
    plt.colorbar()
    plt.quiver(X, Y, dgdk_values/magnitudes, dgdtilt_values/magnitudes,
               label="Backpropagation", color="blue", angles='xy')
    plt.quiver(X, Y, dgdk_values_finite/magnitudes, dgdtilt_values_finite/magnitudes,
               label="Finite diff", color="red", angles='xy')
    plt.legend()
    plt.show()

    # save the data to disk
    np.save("data/x_data.npy", X)
    np.save("data/y_data.npy", Y)
    np.save("data/g_data.npy", g_values)
    np.savez("data/dgdtilt_data.npy",
             dgdk_values,
             dgdtilt_values,
             dgdk_values_finite,
             dgdtilt_values_finite)
