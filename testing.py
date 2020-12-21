# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 00:50:30 2020

@author: Karim
"""

import numpy as np

import test_linalg as lat
import test_boundaries as bct
import test_method as mt
import test_resolution as rt
import test_system as st

def test_all():
    """
    Function to perform all of the tests in the project

    Returns
    -------
    passed : bool
        If all tests were passed.

    """
    
    statement = """\n\nThis runs 7 tests:
    1. Jacobian method consistency
    2. LU inversion consistency
    3. Internal linear algebra and numpy linear algebra consistency.
    4. CDS vs FDS boundary conditions
    5. Exact, Relaxation and Newton method consistency
    6. Resolution consistency
    7. System setup conditions
"""
    
    print(statement)
    input("press enter to start tests.")
    l1 = lat.jac_check_consistency()
    l2 = lat.lu_check_consistency()
    l3 = lat.linalg_compare_numpy()
    b1 = bct.check_bcs()
    m1 = mt.compare_methods()
    r1 = rt.check_res()
    s1 = st.check_mat_hs_sys()
    
    results = np.array([l1, l2, l3, b1, m1, r1, s1])
    
    conf = np.sum(results)
    passed = np.all(results)
    
    print("\n\nPassed: {} --- Passed all: {}".format(conf, passed))
    return passed