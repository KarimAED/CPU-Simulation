# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 09:09:05 2020

@author: Karim

This checks for all of the manual linear algebra functions working properly.
"""
import numpy as np
import time

import defaults as d
import solver
import linalg

def jac_check_consistency():
    """
    A test for the jacobian matrix equation solver.

    Returns
    -------
    passed : bool
        If the solver is self-consistent.

    """
    print("\n\nTest: Jacobian Self Consistency")
    
    # random diagonal matrix is known to be solvable
    vec = np.ones(30)
    diag = np.random.randint(1, 30, 30)*np.random.choice([-1, 1], 30)
    mat = np.diag(diag)
    
    # check for differences between parameter and prediction
    x = linalg.jac_solve(mat, vec)
    vec_guess = mat.dot(x)
    e = np.abs(vec-vec_guess)
    passed = np.all(e < 1e-5)
    
    print("Passed {}".format(passed))
    return passed

def lu_check_consistency():
    """
    A test for the LU inversion solver.

    Returns
    -------
    passed : bool
        If the LU-inverter is self-consistent.

    """
    print("\n\nTest: LU inversion Self Consistency")
    
    # random diagonal matrix is known to be invertible
    diag = np.random.randint(1, 30, 30)*np.random.choice([-1, 1])
    mat = np.diag(diag)
    inv = linalg.inv(linalg.lu_solve, mat)
    
    
    # check for differences in error and value
    id_guess = inv.dot(mat)
    e = np.abs(np.identity(30)-id_guess)
    passed = np.all(e < 1e-7)
    print("Passed {}".format(passed))
    return passed

def linalg_compare_numpy():
    """
    This is a function comparing the self-written linear algebra results with
    the numpy solver results. This is done on the system without heat sink and
    with 'NATURAL' convection.
    
    It shows that the error is sufficiently small to
    ensure that the self-written algorithm works. It furthermore highlights,
    why the numpy solver is used instead for the rest of the implementation, 
    as it's runtime is much lower.

    Returns
    -------
    bool
        If the absolute error in peak temperature between the two algorithms
        is larger or smaller than 0.01K.

    """
    print("\n\nTest: Comparing manual linalg solver to numpy.\n\n")
    
    # set up system
    no_sink_system = d.no_hs_system()

    print("NUMPY:")
    
    # get prediction with numpy linear algebra
    sol = solver.NewtonSolver(no_sink_system)
    start = time.time()
    sol.solve(1e-3)
    end = time.time()
    m_temp_np = np.max(sol.temp)
    np_time = end-start


    print("\n\nMANUAL:")
    
    # get prediction with manual linear algebra
    sol2 = solver.NewtonSolver(no_sink_system)
    start = time.time()
    sol2.solve(1e-3, kp_linalg=True, verbose=True)
    end = time.time()
    m_temp_kp = np.max(sol2.temp)
    man_time = end-start

    err = np.abs(m_temp_np-m_temp_kp)
    passed = err < 1e-2

    print("\n\nNumpy max temp: {} - time: {} / Manual max temp: {} - time: {} \
          / error: {}\n\nPassed: {}".format(
        round(m_temp_np, 3), np_time, round(m_temp_kp, 3), man_time,
        round(np.abs(m_temp_np-m_temp_kp), 3), passed))

    return passed