#!/usr/bin/env python3

from solve_system import solve_system
from recorder import SimulationReader
from backpropagation import Backpropagation
from symulathon import Simulation
from interpolate_parameters import generate_parameters
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from tqdm import tqdm


def newton_iteration(sim: Simulation, x0, v0, xi, vi):
    A = sim.getEquationMatrix()
    dfdx = sim.getForcePositionJacobian()
    f = sim.getForce()
    # M v0 + h fi + h df/dx (x0 – xi) - (M – h² df/dx) vi
    b = mass * v0 + h * f + h * dfdx @ (x0 - xi) - (mass - h**2 * dfdx) @ vi
    delta_v = solve_system(A, b)
    v1 = vi + delta_v
    x1 = x0 + h * v1
    return x1, v1


def simulate(k, k_bend):
    reader = SimulationReader(nDoF)
    backpropagation = Backpropagation(mass, h)
    sim = Simulation(k, k_bend)
    sim.fill_containers()
    for i in range(DIFF_FRAMES+1):
        ##################################
        # Record step for backpropagation
        ##################################
        x = sim.getPosition()
        v = sim.getVelocity()
        x_t, v_t = reader.get_next_state()
        A = sim.getEquationMatrix()
        dfdp = sim.getParameterJacobian()
        dfdx = sim.getForcePositionJacobian()
        backpropagation.step(x, v, x_t, v_t, A, dfdp, dfdx)

        ##################################
        # Newton Iterations
        ##################################
        iterations = 3
        xi = x
        vi = v
        for it in range(iterations):
            xi, vi = newton_iteration(sim, x, v, xi, vi)
            sim.set_state(xi, vi)
            sim.fill_containers()

    return backpropagation


if __name__ == "__main__":
    # Initialize important constant variables
    sim = Simulation(1, 1)
    nDoF = sim.getDoF()
    mass = sim.getMassMatrix()
    h = sim.getTimeStep()
    DIFF_FRAMES = 10

    # Define the region studied
    n_points = 20
    # k_values = np.linspace(0.01, 10, n_points)
    # k_bend_values = np.linspace(-10, 500, n_points-1) / 100
    k_values = np.linspace(0.5, 5, n_points)
    k_bend_values = np.linspace(1, 500, n_points-1) / 100
    X, Y = np.meshgrid(k_values, k_bend_values)

    g_values = X.tolist()

    # Computation of the surface and gradients
    for i in tqdm(range(len(k_values))):
        for j in range(len(k_bend_values)):
            k4 = k_values[i] * np.ones(4)
            kbend4 = k_bend_values[j] * np.ones(4)
            k, kb = generate_parameters(k4, kbend4)
            bp = simulate(k, kb)
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
