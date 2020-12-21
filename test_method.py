# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 10:03:58 2020

@author: Karim
"""
import numpy as np

import solver
import defaults as d


def compare_methods():
    """
    A test comparing the different solvers.

    Returns
    -------
    passed : bool
        If the solvers are consistent with each other.

    """
    print("\n\nTest: Method cross-consistency.")
    
    # set up system (no heat sink)
    s = d.no_hs_system()
    
    # set up solvers with forced convection
    ex_solv = solver.ExactSolver(s, v=10)
    rel_solv = solver.RelaxSolver(s, "F", v=10)
    newt_solv = solver.NewtonSolver(s, "F", v=10)
    
    print("\nExact Solver:")
    ex_pred = ex_solv.solve()
    print("\nRelaxation Solver:")
    # Takes longer to converge, needs lower threshold
    rel_pred = rel_solv.solve(1e-5, verbose=True)
    print("\nNewton Solver:")
    newt_pred = newt_solv.solve(1e-3)
    
    
    # check for differences in predictions
    preds = np.array([ex_pred, rel_pred, newt_pred])
    preds -= ex_pred
    print(preds)
    
    passed = np.all(preds < 1e-2)
    print("Passed: {} - Difference: <{}".format(passed, 1e-2))
    return passed