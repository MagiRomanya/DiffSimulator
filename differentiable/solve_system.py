#!/usr/bin/env python3

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


def solve_system(eq_mat, eq_vec, threshold=1e-6, maxiter=200):
    """
    Solves the sparse linear system defined by one eq_mat & eq_vec.

    The function uses scipy's conjugate gradient method to approximate
    the results. The user can adjust the threshold & maximum iterations.
    """
    result, convergence = scipy.sparse.linalg.cg(eq_mat,
                                                 eq_vec,
                                                 tol=threshold,
                                                 maxiter=maxiter)
    # result = scipy.sparse.linalg.spsolve(eq_mat, eq_vec)
    check_cg_convergence(convergence)
    return result
