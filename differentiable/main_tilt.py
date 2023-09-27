#!/usr/bin/env python3
"""
STUDY: Computing the loss derivative wrt tilt parameter.
Single parameter derivative scan.
"""
from symulathon import Simulation
from simulation_functions import newton_iteration, get_user_weather_graphics
from recorder import SimulationReader
from backpropagation import Backpropagation
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

k = 20
k_bend = 0.1
sim = Simulation(k, k_bend, 0, get_user_weather_graphics())

def simulate(tilt_angle: float, diff_frames: int):
    nDoF = sim.getDoF()
    mass = sim.getMassMatrix()
    h = sim.getTimeStep()
    reader = SimulationReader(nDoF)
    sim.reset_simulation(k, k_bend, tilt_angle)
    dx0dp = sim.getInitialPositionJacobian()
    dv0dp = sim.getInitialVelocityJacobian()
    backpropagation = Backpropagation(mass, dx0dp, dv0dp, h)
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

        sim.render_state()
    return backpropagation

if __name__ == "__main__":
    nDoF = sim.getDoF()
    mass = sim.getMassMatrix()
    h = sim.getTimeStep()
    DIFF_FRAMES = 30

    angle_values = np.linspace(0, 360, 300)
    g_values = []
    dgdp_values = []
    for angle in tqdm(angle_values):
        bp = simulate(np.radians(angle), DIFF_FRAMES)
        g_values.append(bp.get_g())
        dgdp_values.append(bp.get_dgdp()[3])

    # Calculate finite differences
    dgdp_finite = []
    for i in range(len(g_values)-1):
        value = (g_values[i+1]-g_values[i]) / (angle_values[i+1] - angle_values[i])
        dgdp_finite.append(value)
    dgdp_finite.append(dgdp_finite[-1])

    plt.plot(angle_values, g_values, "-", label="Loss Function")
    plt.plot(angle_values, dgdp_finite, ".", label="Finite dgdp")
    plt.plot(angle_values, dgdp_values, "x", label="Backpropagation dgdp")

    plt.legend()
    plt.xlabel("Tilt angle value")
    plt.grid()
    plt.show()
    print("Finished succesfully")
