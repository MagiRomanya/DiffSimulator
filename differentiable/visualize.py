#!/usr/bin/env python3

from symulathon import Simulation
from recorder import SimulationReader
import numpy as np

sim = Simulation(0,0,True)
dof = sim.getDoF()

reader = SimulationReader(dof)
trajectory = np.array(reader.get_all_history()).T
n_states = len(trajectory)

print(f"Number of DoF : {dof}")
print(f"Number of states : {n_states}")

current_state_index = 0
KEY_COMMA = 44
KEY_PERIOD = 46
KEY_P = 80
KEY_R = 82
GAME_PAUSED = True

def set_state(index: int):
    global current_state_index
    current_state_index = index
    x, v = trajectory[index][:dof], trajectory[index][dof:]
    sim.set_state(x, v)

def next_state():
    global current_state_index
    current_state_index += 1
    current_state_index = min(n_states-1, current_state_index)
    set_state(current_state_index)

def unext_state():
    global current_state_index
    current_state_index -= 1
    current_state_index = max(0, current_state_index)
    set_state(current_state_index)

set_state(current_state_index)

while not sim.window_should_close():
    sim.render_state()
    if GAME_PAUSED:
        if (sim.is_key_pressed(KEY_PERIOD)):
            next_state()
        elif (sim.is_key_pressed(KEY_COMMA)):
            unext_state()
    else:
        next_state()

    if sim.is_key_pressed(KEY_P):
        GAME_PAUSED = not GAME_PAUSED

    if sim.is_key_pressed(KEY_R):
        set_state(0)

    print(f"Frame = {current_state_index}")
