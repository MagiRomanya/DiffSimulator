#!/usr/bin/env python3
"""Minimize function 2 differentiable parameters."""

from symulathon import Simulation
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from colorama import Fore, Style
from simulation_functions import simulate, get_user_plots
from animation_plots import AnimatedPlot
from scipy.optimize import minimize, Bounds


def simulation_wrapper(x):
    bp = simulate(x[0], x[1], DIFF_FRAMES)
    g = bp.get_g()
    dgdp = bp.get_dgdp()
    # print(f"g = {g}, params = {x}")
    return g, dgdp


if __name__ == "__main__":
    # Initialize important constant variables
    sim = Simulation(1, 1)
    nDoF = sim.getDoF()
    mass = sim.getMassMatrix()
    h = sim.getTimeStep()
    DIFF_FRAMES = 100

    #load g values
    X = np.load("data/x_data.npy")
    Y = np.load("data/y_data.npy")
    g = np.load("data/g_data.npy")
    maxvalue= 2e6
    plt.contourf(X, Y, g, 1000, vmax=maxvalue)
    plt.colorbar()
    scatter = plt.scatter([1],[1])
    arrow = plt.arrow(1, 1, 0.1, 0.1,color="white")

    x0 = np.array([2, 2])
    # x0 = np.array([3.5, 0.35])
    MAX_ITERATIONS = 1000
    ALPHA = 1e-7
    previous_g = 0
    bool_plots = get_user_plots()
    res = minimize(simulation_wrapper,
                   x0,
                   bounds=Bounds(0, 1000),
                   jac=True,
                   options={"disp": True})
    print(res)
    if bool_plots:
        plot = AnimatedPlot()
    # for i in range(MAX_ITERATIONS):
    #     bp = simulate(x0[0], x0[1], DIFF_FRAMES)
    #     g = bp.get_g()
    #     if bool_plots:
    #         plot.get_data(g)
    #     diff = g - previous_g
    #     previous_g = g
    #     diff_str = ""
    #     if (diff < 0):
    #         diff_str = Fore.GREEN + str(diff) + Style.RESET_ALL
    #     else:
    #         diff_str = Fore.RED + str(diff) + Style.RESET_ALL

    #     print(f"Iteration #{i}\t loss = {g:10.4e},",
    #           "Difference = "+diff_str)
    #     dgdp = bp.get_dgdp()
    #     scatter.set_offsets(x0)
    #     dx, dy = -ALPHA*dgdp*100
    #     arrow.set_data(x=x0[0], y=x0[1], dx=dx, dy=dy)
    #     plt.pause(0.01)
    #     x0 = x0 - ALPHA * dgdp
    #     print(f"New parameters {x0}")
