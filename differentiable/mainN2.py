#!/usr/bin/env python3
"""Computes gradient descent for a interpolated parameter cloth."""

from symulathon import Simulation, count_springs
from animation_plots import AnimatedPlot
import numpy as np
from interpolate_parameters import generate_parameters, get_parameter_indices
from simulation_functions import simulate, get_user_plots
from scipy.optimize import minimize, Bounds


def simulation_wrapper(x):
    k4 = np.array(x)
    k_bend4 = np.ones(4)*0.1
    k_values, k_bend_values = generate_parameters(k4, k_bend4)
    bp = simulate(k_values, k_bend_values, DIFF_FRAMES)
    g = bp.get_g()
    # print(f"Loss = {g:10.4e}")
    dgdp = bp.get_dgdp()
    dgdp = dgdp[indx][:4]
    return g, dgdp


if __name__ == "__main__":
    # Initialize important constant variables
    sim = Simulation(1, 1)
    nFlex, nBend = count_springs()
    print(f"nFlex = {nFlex}, nBend = {nBend}, total_diff = {nFlex + nBend}")
    nDoF = sim.getDoF()
    mass = sim.getMassMatrix()
    h = sim.getTimeStep()

    # Define inital state
    indx = get_parameter_indices()
    # Minimization process
    MAX_ITER = 1000
    ALPHA = 0.1
    last_loss = 0

    DIFF_FRAMES = 100
    x0 = np.ones(4) * 1
    # x0 = np.array([3.989, 3.967, 4.079,  4.160])
    res = minimize(simulation_wrapper, x0, jac=True,
                   bounds=Bounds(lb=0, ub=200),
                   tol=1e-10,
                   options={
                       "disp": True,
                       "maxfun": 100,
                       "maxls": 30,
                       "maxiter": 100},
                   )
    print(res)
