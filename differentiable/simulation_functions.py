from solve_system import solve_system
from recorder import SimulationReader
from backpropagation import Backpropagation
from symulathon import Simulation
import sys
import getopt


sim = Simulation()
nFlex, nBend = sim.getSpringNumbers()
nDoF = sim.getDoF()
mass = sim.getMassMatrix()
h = sim.getTimeStep()


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


def get_user_weather_graphics():
    """Return weather the user wants a graphical simulation."""
    enable_graphics = False
    argv = sys.argv[1:]
    opts, args = getopt.getopt(argv, "gp")
    for opt, arg in opts:
        if opt in ['-g']:
            enable_graphics = True
    return enable_graphics


def get_user_plots():
    """Return weather the user wants a plot of the gradient descent convergence."""
    enable_plots = False
    argv = sys.argv[1:]
    opts, args = getopt.getopt(argv, "gp")
    for opt, arg in opts:
        if opt in ['-p']:
            enable_plots = True
    return enable_plots


def simulate(k_list, k_bend_list, DIFF_FRAMES: int):
    """Make a simulation with the specified parameters and return the backpropagation information."""
    reader = SimulationReader(nDoF)
    sim.reset_simulation(k_list, k_bend_list)
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
        iterations = 1
        xi = x
        vi = v
        for it in range(iterations):
            xi, vi = newton_iteration(sim, x, v, xi, vi)
            sim.set_state(xi, vi)
            sim.fill_containers()

        sim.render_state()
    return backpropagation


def simulate_tilt(k_value, tilt_angle, DIFF_FRAMES: int):
    """Make a simulation with the specified parameters and return the backpropagation information."""
    reader = SimulationReader(nDoF)
    k_bend_value = 0.1
    sim = Simulation(k_value, k_bend_value, tilt_angle, get_user_weather_graphics())
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
