#!/usr/bin/env python3

"""
CAUTION:
 Modifying the tolerance and the maximum iterations of the solver could
 potentially destablize the simulation and the simulation derivatives.

A small tolerance could make the solver not work by imposing too much accuracy to the solver.
A high tolerance could make the simulation not accurate enough and this would be reflected in both
forward and backward passes.

I have found that a value of tol~1e-5 is good enough for the simulations we do here.
Actually when introducing contact it is best to have lower tolerances 1e-6, 5000 iterations.

The maximum iterations of the solver is also important but not as much as the tolerance.
THEORY: If the system is not defined positive, imposing a strict tolerance will put strain in
the solver making it diverge as it is not designed for this kind of solves.
By using somewhat big tolerances, the solver can converge into a good enough solution, and we
do not force it to go in the divergent regime.
"""

import scipy


def check_cg_convergence(convergence: int):
    """
    Check weather the conjugate gradient method has converged or not.

    Helper function which outputs a warning to the console when the conjugate
    gradient method has had some kind of problem.
    """
    if convergence > 0:
        print(f"Warning: conjugate gradient did not converge\
        ({convergence} iterations)")
    elif convergence < 0:
        print("Warning: conjugate gradient illegal input")


def solve_system(eq_mat, eq_vec, threshold=1e-6, maxiter=500):
    """
    Solves the sparse linear system defined by one eq_mat & eq_vec.

    The function uses scipy's conjugate gradient method to approximate
    the results. The user can adjust the threshold & maximum iterations.
    """
    result, convergence = scipy.sparse.linalg.cg(eq_mat,
                                                 eq_vec,
                                                 tol=threshold,
                                                 maxiter=maxiter)
    check_cg_convergence(convergence)

    # result = scipy.sparse.linalg.spsolve(eq_mat, eq_vec)
    return result
