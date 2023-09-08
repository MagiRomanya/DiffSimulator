#!/usr/bin/env python3

from solve_system import solve_system
from recorder import SimulationReader
from backpropagation import Backpropagation
from symulathon import Simulation
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from joblib import Parallel, delayed


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


def simulate(k):
    reader = SimulationReader(nDoF)
    backpropagation = Backpropagation(mass, h)
    sim = Simulation(k, k/100)
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
        iterations = 2
        xi = x
        vi = v
        for it in range(iterations):
            xi, vi = newton_iteration(sim, x, v, xi, vi)
            sim.set_state(xi, vi)
            sim.fill_containers()
            # delta_v = vi - v
            # f = symulathon.get_force()
            # print(f"{it+1}th iteration convergence check = {sum(mass @ delta_v - h * f)}")

    dgdp = backpropagation.get_dgdp()
    g = backpropagation.get_g()
    return (g, dgdp)


if __name__ == "__main__":
    sim = Simulation(1, 1)
    nDoF = sim.getDoF()
    mass = sim.getMassMatrix()
    h = sim.getTimeStep()
    DIFF_FRAMES = 30

    k_values = np.linspace(0.01, 10, 200)

    results = Parallel(n_jobs=10)(delayed(simulate)(k) for k in tqdm(k_values))
    g_values = [results[i][0] for i in range(len(results))]
    dgdp_values = [results[i][1][0] for i in range(len(results))]

    # Calculate finite differences
    dgdp_finite = []
    for i in range(len(g_values)-1):
        value = (g_values[i+1]-g_values[i]) / (k_values[i+1] - k_values[i])
        dgdp_finite.append(value)
    dgdp_finite.append(dgdp_finite[-1])

    plt.plot(k_values, g_values, "-", label="Loss Function")
    plt.plot(k_values, dgdp_finite, ".", label="Finite dgdp")
    plt.plot(k_values, dgdp_values, "x", label="Backpropagation dgdp")

    plt.legend()
    plt.xlabel("k value")
    plt.grid()
    plt.show()
    print("Finished succesfully")
