# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 10:52:12 2020

@author: Karim
"""

import defaults as d
import solver


def check_bcs():
    """
    This function compares the forward and central difference scheme.
    It clearly shows that the central difference scheme returns lower
    temperatures.
    
    This is expected behaviour, as explained in the report (link in README.md)

    Returns
    -------
    passed : TYPE
        DESCRIPTION.

    """
    print("\n\nTest: Do cds BCs return lower temperatures than fds BCs")
    hs_sys = d.hs_system(20, 20, 2)
    sol =solver.NewtonSolver(hs_sys)
    
    print("\nCentral Difference Scheme:")
    cds_t = sol.solve(1e-3, cds=True)
    
    print("\nForward Difference Scheme:")
    fds_t = sol.solve(1e-3, cds=False)
    
    passed = cds_t < fds_t
    
    print("\n\n\nPassed: {} - Difference: {}".format(passed,
                                                     round(cds_t-fds_t, 2)))
    return passed