# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 18:58:46 2020

@author: Karim
"""
import solver as s
import system as sys


def no_hs_system(show=False, res=2):
    """
    This sets up a default system without a heat sink.
    
    Resolution can be specified, as well as optional plotting of the system.

    Parameters
    ----------
    show : bool, optional
        If the system should be plotted after setup. The default is False.
    res : int, optional
        The resolution at which to generate the system. The default is 2.
        The minimum is 2.

    Returns
    -------
    no_hs_sys : sys.System
        A sys.System dataclass instance without a sys.HeatSink object.

    """
    
    # set up materials
    air = sys.Material("Air")
    processor = sys.Material("Processor", power=0.5, conductivity=0.15)
    case = sys.Material("Case", conductivity=0.23)
    
    materials = [air, processor, case]
    
    # set up the system
    no_hs_sys = sys.System("System (without Heat Sink)",materials,
                           resolution=res)
    
    # plot layout if required
    if show:
        no_hs_system.show()
    
    return no_hs_sys



def hs_system(fin_len, fin_ct, fin_dist, show=False):
    """
    This sets up a default system with a heat sink.
    
    Resolution is fixed at 2, but plotting can be specified.

    Parameters
    ----------
    fin_len : int
        Length of heat sink fins in mm.
    fin_ct : int
        Number of heat sink fins.
    fin_dist : int
        Distance between heat sink fins in mm.
    show : bool, optional
        If the system layout should be plotted. The default is False.

    Returns
    -------
    hs_sys : sys.System
        A sys.System instance containing a sys.HeatSink object.

    """
    # fixed resolution
    res = 2
    
    # set up materials
    air = sys.Material("Air")
    processor = sys.Material("Processor", power=0.5, conductivity=0.15)
    case = sys.Material("Case", conductivity=0.23)
    alum = sys.Material("Heat Sink", conductivity=0.25)
    
    materials = [air, processor, case, alum]
    
    # set up heat sink
    hs = sys.HeatSink(fin_len, fin_dist, fin_ct)
    
    # set up system
    hs_sys = sys.System("System (with Heat Sink)", materials,
                        resolution=res, heat_sink=hs)
    
    # show layout if required
    if show:
        hs_sys.show()
    
    return hs_sys
    

def no_hs_modelling(err=1e-3, kp_linalg=False):
    """
    Perform default modelling of a system without a heat sink and natural
    convection.
    
    Uses solver.NewtonSolver.
    
    Automatically plots a heat map of the system.

    Parameters
    ----------
    err : float, optional
        The error threshold at which the newton solver terminates.
        The default is 1e-3.
    kp_linalg: bool, optional
        If internal linear algebra should be used instead of numpy.
        The default is False.

    Returns
    -------
    float
        The average temperature across the processor.

    """
    system = no_hs_system()
    solver = s.NewtonSolver(system)
    return solver.solve(err, verbose=True, show=True, kp_linalg=kp_linalg)


def hs_nat_modelling(fin_len, fin_ct, fin_dist, err=1e-3, 
                     verbose=True, show=False, kp_linalg=False):
    """
    Performs a default modelling of a system with a heat sink and natural
    convection.
    
    Uses solver.NewtonSolver.

    Parameters
    ----------
    fin_len : int
        Length of heat sink fins in mm.
    fin_ct : int
        Number of heat sink fins.
    fin_dist : int
        Distance between heat sink fins in mm.
    err : float, optional
        The error threshold at which the newton solver terminates.
        The default is 1e-3.
    verbose : bool, optional
        If verbose output should be printed. The default is True.
    show : bool, optional
        If the system heat map should be shown after solving.
        The default is False.
    kp_linalg: bool, optional
        If internal linear algebra should be used instead of numpy.
        The default is False.

    Returns
    -------
    float
        The average temperature across the processor.

    """
    
    system = hs_system(fin_len, fin_ct, fin_dist)
    solver = s.NewtonSolver(system)
    
    # use verbose for longer runs
    return solver.solve(err, verbose=verbose, show=show, kp_linalg=kp_linalg) 


def hs_for_modelling(fin_len, fin_ct, fin_dist, err=1e-3, vel=20, 
                     verbose=True, show=False, kp_linalg=False):
    """
    Performs a default modelling of a system with a heat sink and forced
    convection.
    
    Uses solver.NewtonSolver.

    Parameters
    ----------
    fin_len : int
        Length of heat sink fins in mm.
    fin_ct : int
        Number of heat sink fins.
    fin_dist : int
        Distance between heat sink fins in mm.
    err : float, optional
        The error threshold at which the newton solver terminates.
        The default is 1e-3.
    vel: float, optional
        The wind speed for the forced convection in m/s.
        The default is 20.
    verbose : bool, optional
        If verbose output should be printed. The default is True.
    show : bool, optional
        If the system heat map should be shown after solving.
        The default is False.
    kp_linalg: bool, optional
        If internal linear algebra should be used instead of numpy.
        The default is False.

    Returns
    -------
    float
        The average temperature across the processor.

    """
    
    system = hs_system(fin_len, fin_ct, fin_dist)
    solver = s.NewtonSolver(system, "F", v=vel)
    
    # use verbose for longer runs
    return solver.solve(err, verbose=verbose, show=show, kp_linalg=kp_linalg)
