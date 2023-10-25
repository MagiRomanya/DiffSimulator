#!/usr/bin/env python3

from symulathon import Simulation
import numpy as np
from recorder import SimulationRecorder

sim = Simulation(0,0,True)
rec = SimulationRecorder("record_shape.csv")
dof = sim.getDoF()
N, M = sim.getGridDimensions()

x = sim.getPosition()

delta = x[3] - x[0]
start_y = x[1] + 10
side_length = 5
n_nodes = 20
sphere_radius = side_length / (2.*np.pi)


def sphere_point(angle):
    return sphere_radius * np.cos(angle), start_y + sphere_radius * np.sin(angle)

# We are working in the x-y plane, z axis is depth
def make_cilinder():
    dangle = 2.*np.pi / 20
    angle = np.pi / 2. + dangle / 2.
    for i in range(0, len(x), 3):
        print(f"x={x[i]}, y={x[i+1]}, z={x[i+2]}")

        xx, yy = sphere_point(angle)
        x[i] = xx
        x[i+1] = yy
        angle += dangle

def rotate_vertices(x: np.array, angle: float):
    y_offset = 10
    c = np.cos(angle)
    s = np.sin(angle)
    Rx = np.array([
        [1, 0, 0],
        [0, c, -s],
        [0, s, c]
    ])
    Ry = np.array([
        [c, 0, s],
        [0, 1, 0],
        [-s, 0, c]
    ])
    Rz = np.array([
        [c, -s, 0],
        [s, c, 0],
        [0,0,1]
    ])
    for i in range(0, len(x), 3):
        pos = np.array([x[i], x[i+1], x[i+2]])
        new_pos = Rx @ pos
        x[i] = new_pos[0]
        x[i+1] = new_pos[1] +y_offset
        x[i+2] = new_pos[2]


make_cilinder()
# rotate_vertices(x, np.pi/2)
sim.set_state(x, np.zeros(dof))

for i in range(100):
    rec.record_timestep(np.zeros(dof), np.zeros(dof))

rec.record_timestep(x, np.zeros(dof))
while not sim.window_should_close():
    sim.render_state()
