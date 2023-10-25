#!/usr/bin/env python3

from symulathon import Simulation
import numpy as np
from tqdm.notebook import tqdm
from scipy.optimize import minimize, Bounds
from simulation_functions import newton_iteration
from recorder import SimulationReader
from backpropagation import Backpropagation
k = 100
k_bend = 0.1
graphics=True
sim = Simulation(k, k_bend, 0, graphics)
nDoF = sim.getDoF()
DIFF_FRAMES = 100

def simulate(initial_velocities: np.array, diff_frames: int):
    nDoF = sim.getDoF()
    mass = sim.getMassMatrix()
    h = sim.getTimeStep()
    reader = SimulationReader(nDoF,"record_shape.csv")
    sim.reset_simulation(initial_velocities)
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

def simulate2(k: float, k_bend: float, initial_velocities: np.array, diff_frames: int):
    nDoF = sim.getDoF()
    mass = sim.getMassMatrix()
    h = sim.getTimeStep()
    reader = SimulationReader(nDoF,"record_shape.csv")
    sim.reset_simulation(k, k_bend, initial_velocities)
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

initial_velocities = np.zeros(nDoF)
initial_tension = np.array([100])
initial_bending = np.array([0.1])
parameters = np.concatenate((initial_tension, initial_bending, initial_velocities))
print(parameters)

def simulation_wrapper(initial_velocities):
    bp = simulate2(0, 0, initial_velocities, DIFF_FRAMES)
    g = bp.get_g()
    dgdp = bp.get_dgdp()[3:]
    return g, dgdp

def simulation_wrapper2(paramters):
    k = paramters[0]
    k_bend = paramters[1]
    ivel = paramters[2:]
    bp = simulate2(k, k_bend, ivel, DIFF_FRAMES)
    g = bp.get_g()
    dgdp = bp.get_dgdp()
    dgdp = np.delete(dgdp, 2) # remove garbage paramter
    return g, dgdp

print(len(parameters))
res = minimize(simulation_wrapper, initial_velocities, jac=True,
               bounds=Bounds(-100, 100),
               options={"disp": True,
                        "maxiter": 500,
                       })
print(res)
print(res.x)
