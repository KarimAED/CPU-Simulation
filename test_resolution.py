# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 10:12:30 2020

@author: Karim
"""
import numpy as np
import matplotlib.pyplot as plt

import defaults as d
import solver

def check_res(plot=False):
    """
    Test to check if the resolution above 2 pxl/mm makes a significant
    difference.

    Returns
    -------
    bool
        If the difference between all of the different resolutions is less than
        1°C.

    """
    print("\n\nTest: Resolution invariance.")
    m_temp = [] # set up array to store temperatures
    
    # resolution range
    resolutions = np.arange(2, 11)
    for res in resolutions:
        print("\nResolution: {} pxl/mm".format(res))
        # set up system
        no_sink_system = d.no_hs_system(res=int(res))
        
        # solve system
        sol = solver.NewtonSolver(no_sink_system)
        m_temp.append(sol.solve(1e-3))
        
    # get differences
    m_temp = np.array(m_temp)
    m_temp -= m_temp[0]
    
    
    # plot differences
    if plot:
        plt.xlabel("Resolution in pxl/mm")
        plt.ylabel("Difference(°C) to Resolution = 2 pxl/mm")
        plt.grid()
        plt.plot(resolutions, m_temp)
        plt.show()
    
    # accurate to within 1°C
    passed = np.all(m_temp < 1)
    print("Passed: {}".format(passed))
    
    return passed