#!/usr/bin/env python3

from solve_system import solve_system
from recorder import SimulationRecorder
# from interpolate_parameters import generate_parameters
from symulathon import Simulation
import numpy as np


def newton_iteration(sim: Simulation, x0, v0, xi, vi):
    """Update xi and vi to a new Newton-Rhapson iteration."""
    A = sim.getEquationMatrix()
    dfdx = sim.getForcePositionJacobian()
    f = sim.getForce()
    # M v0 + h fi + h df/dx (x0 – xi) - (M – h² df/dx) vi
    b = mass * v0 + h * f + h * dfdx @ (x0 - xi) - (mass - h**2 * dfdx) @ vi
    delta_v = solve_system(A, b)
    v1 = vi + delta_v
    x1 = x0 + h * v1
    return x1, v1


if __name__ == "__main__":
    rec = SimulationRecorder()
    K_VALUE = 3
    RECORD_FRAMES = 500
    sim = Simulation()

    n_tension, n_bending = sim.getSpringNumbers()
    param = np.genfromtxt("paramters_gioconda_xd.csv", delimiter=",")
    print(len(param))
    k = 0
    k_bend = 0
    vel = param

    sim = Simulation(0,0, True)
    sim.reset_simulation(k, k_bend, vel)
    h = sim.getTimeStep()
    frames_count = 0
    while (frames_count <= RECORD_FRAMES):
        frames_count += 1

        sim.fill_containers()
        x = sim.getPosition()
        v = sim.getVelocity()
        mass = sim.getMassMatrix()
        rec.record_timestep(x, v)

        A = sim.getEquationMatrix()
        b = sim.getEquationVector()

        # delta_v = solve_system(A, b)
        # v_new = v + delta_v
        # x_new = x + h * v_new
        # sim.set_state(x_new, v_new)

        iterations = 10
        xi = x
        vi = v
        for it in range(iterations):
            xi, vi = newton_iteration(sim, x, v, xi, vi)
            sim.set_state(xi, vi)
            sim.fill_containers()

        sim.render_state()

    print("Done")
