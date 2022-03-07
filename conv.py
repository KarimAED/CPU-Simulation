# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 18:47:25 2020

@author: Karim
"""


def h_0(T):
    """
    Natural Convection in W/mm^2K

    Parameters
    ----------
    T : float
        The temperature difference between the surface and ambient temperature.
        In Kelvin.

    Returns
    -------
    float
        The amount of natural convection for the given T.

    """
    if not isinstance(T, (float, int)):
        raise TypeError("T must be a float or an int.")
    return 1.31*(T**(1/3))*1e-6


def h_f(v):
    """
    Forced Convection in W/mm^2K

    Parameters
    ----------
    v : float
        Wind speed in meters per second.

    Returns
    -------
    float
        The amount of forced convection at this wind speed.

    """
    if not (isinstance(v, float) or isinstance(v, int)):
        raise TypeError("v must be a float or an int.")
    return (11.4+5.7*v)*1e-6


def h(T, v):
    """
    The convection function in W/mm^2K

    Parameters
    ----------
    T : float
        Surface temperature difference to ambient temperature, in Kelvin.
    v : float
        Wind speed in m/s. Set to 0 to evaluate with natural convection.

    Returns
    -------
    float
        convection parameter.

    """
    if not isinstance(v, (float, int)):
        raise TypeError("v must be a float or an int.")
    if v == 0:
        return h_0(T)
    return h_f(v)
