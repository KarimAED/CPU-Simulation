# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 19:39:21 2020

@author: Karim
"""

import defaults as d
import optimization as opt
import testing

in_str = \
"""
This is the solution to the ICL Computational Physics course project 4.
Written by K. Alaa El-Din.
"""

loop_str = \
"""
Please enter a function you want to use:
    q, quit, exit - Quits the program
    t, test - Runs all tests (may take a long time)
    grid, opt - Accesses the heat sink optimiser (may also take a long time)
    sys, solve, s - Accesses the regular system solving
    p, pre, preset - Calls systems with all of the presets shown in the report
    h, help - Prints an extended help message


"""

invalid_str = \
"""
That input does not match any of the available options.
Please type "h" or "help" for more information.
"""

help_str = \
"""
The following features are available through this interface:
    
    
    1. Running all Tests, including linear algebra, resolution and methods,
    2. Optimising a heat sink design,
    and
    3. Setting up a single simulation,
    4. Running all of the simulations highlighted in the report.
    
    To do 1., simply enter "t" or "test" and all the tests will be run
    and output their result. A final line will show if all tests were success-
    ful or if they failed.
    
    To do 2., enter "grid", "gs", "grid_search", "opt" or "optimise" and you
    will be able to choose all of the different parameters to perform a grid
    search over the heat sink parameters. The 5 lowest processor temperature
    parameters will be returned.
    
    To do 3., enter "sys", "solve", or "s". You will then be able to specify
    the exact parameters for the system you wish to run or choose from a pre-
    defined systems which are plotted in the report.
    
    To do 4., enter "p", "pre", or "preset". You will then be able to choose
    to run all of the simulations from the report or to exclude the huge slow
    system.
    
    Systems and grid searches beyond 50 fins of 50mm length are highly discour-
    aged. They take extremely long even with optimal parameters.
"""

# used to understand user input
yes = ("y", "yes", "1")
no = ("n", "no", "0")


def cl():
    """
    Clears the console
    """
    print("\n"*100)


def manage_inp(label, defaults):
    """
    Shows user label as input description, defaults as default values.
    
    Then gets user input, splits it, checks for correct types and finally
    has the option of returning defaults.

    Parameters
    ----------
    label : str
        The label shown to the user as what they should input.
    defaults : list[]
        The default values for the properties.

    Returns
    -------
    None.

    """
    print("\n\nEnter "+label+".")
    inp = input("Enter value(s). For more than one value, \
enter seperated by a comma as in '2,4,5' \
(defaults to "+str(defaults)+"): ")

    inp_form = inp.split(",")
    if inp == "":
            out = defaults
            print("Defaulted to "+str(defaults))
    else:
        out = []
        for i, default in enumerate(defaults):
            tp = type(default)
            try:
                out.append(tp(inp_form[i]))
            except:
                raise TypeError("The user input had the wrong type.")
        print("User chose "+str(out))
    
    return out


def run_system():
    """
    Gets the user input, processes it and runs a single system.

    Returns
    -------
    None.

    """
    cl()
    
    print("Select system Parameters.")
    
    inp_mode = input("Should the system use internal linear algebra \
(much slower with heat sink) (y/n)? ")
    
    if (inp_mode not in yes) and (inp_mode not in no):
        return
    
    if inp_mode in yes:
        kp_linalg = True
    elif inp_mode in no:
        kp_linalg = False
    
    inp = input("Should the system use a heat sink (y/n)? ")
    
    # if the user input is unknown
    if (inp not in yes) and (inp not in no):
        return
    
    if inp in no:
        fin_params = []
        mode = "N"

    elif inp in yes:
        label = "Heat sink params: fin length(mm), count, and distance(mm)"
        fin_params = manage_inp(label, [20, 20, 2])
        
        lab_mode = "Convection mode: 'Natural'('N') or 'Forced'('F')"
        mode, = manage_inp(lab_mode, ["N"])
    
    label = "Convergence limit/error"
    err, = manage_inp(label, [1e-3])
    
    if not fin_params:
        d.no_hs_modelling(err, kp_linalg)
    else:
        if mode.upper() in ("NATURAL", "NAT", "N"):
            d.hs_nat_modelling(*fin_params, err, True, True, kp_linalg)
    
        elif mode.upper() in ("FORCED", "FOR", "F"):
            lab_vel = "Convection velocity (m/s)"
            v, = manage_inp(lab_vel, [20])
            
            d.hs_for_modelling(*fin_params, err, v, True, True, kp_linalg)
        
        else:
            raise ValueError("Invalid convection mode identifier.")


def optimise():
    """
    Gets the user input, processes it and generates a grid search.

    Raises
    ------
    ValueError
        If the convection mode identifier is invalid.

    Returns
    -------
    None.

    """
    cl()
    
    print("Choose parameters for grid search.")
    lab_mode = "Convection mode: 'Natural'('N') or 'Forced'('F')"
    mode, = manage_inp(lab_mode, ["N"])

    label = "Fin Length params: start, stop and step (in mm)"
    len_params = manage_inp(label, [10, 50, 10])
    
    label = "Fin count params: start, stop and step (in mm)"
    ct_params = manage_inp(label, [20, 40, 2])

    label = "Fin dist params: start, stop and step (in mm)"
    dist_params = manage_inp(label, [1, 4, 1])
    
    label = "Success threshold temperature"
    temp, = manage_inp(label, [80])

    label = "Newton Method Convergence Error"
    err, = manage_inp(label, [1e-3])
    
    if mode.upper() in ("NATURAL", "NAT", "N"):
        fun = d.hs_nat_modelling
        opt.grid_search_hs(fun, len_params, ct_params, dist_params, temp,
                        [err], plot_best=False, plot_surface=True)

    elif mode.upper() in ("FORCED", "FOR", "F"):
        fun = d.hs_for_modelling
        lab_vel = "Convection velocity (m/s)"
        v, = manage_inp(lab_vel, [20])
        opt.grid_search_hs(fun, len_params, ct_params, dist_params, temp,
                        [err, v], False, True)
    else:
        raise ValueError("Invalid convection mode identifier.")


def run_presets():
    """
    A function to run the simulation for all of the systems plotted in the
    report.

    Returns
    -------
    None.

    """
    cl()
    
    pres_statement = """Run the following presets:
    1. No heat sink and natural convection
    2. (optional, very slow) heat sink with 50 100mm long fins 3mm apart.
    3. Perform the grid search for the forced convection
    4. Run the smallest successful system with forced convection (28 fins,
      1mm apart and 30mm long)
"""
    
    print(pres_statement)
    
    inp = input("Run with extremely slow options (y/n)? ")
    
    if inp in yes:
        slow = True
    elif inp in no:
        slow = False
    else:
        return
    
    print("Running no heat sink system...")
    d.no_hs_modelling(1e-9)
    
    if slow:
        print("Running huge heat sink system with natural convection...")
        d.hs_nat_modelling(100, 50, 3, show=True)
    
    print("Running grid search for forced convection...")
    opt.grid_search_hs(d.hs_for_modelling,
                       [10, 50, 10],
                       [10, 40, 2],
                       [1, 4, 1],
                       args=[1e-3],
                       plot_surface=True)
    
    print("Running smallest successful system for forced convection...")
    d.hs_for_modelling(30, 28, 1, show=True)


def main():
    """
    A main loop for user interactions.

    Returns
    -------
    None.

    """
    i = 0
    ex_con = ("q", "quit", "exit")
    t_con = ("t", "test")
    gs_con = ("grid", "gs", "grid_search", "opt", "optimise")
    sys_con = ("sys", "solve", "s")
    help_con = ("h", "help")
    pre_con = ("p", "preset", "pre")
    print(in_str)
    while i==0:
        cl()
        print(loop_str)
        inp = input("Enter what you wish to do: ").lower()
        cl()
        if inp in ex_con:
            i = 1
        elif inp in t_con:
            print("Testing everything (This may take a while):")
            testing.test_all()
        elif inp in gs_con:
            optimise()
        elif inp in sys_con:
            run_system()
        elif inp in pre_con:
            run_presets()
        elif inp in help_con:
            print(help_str)
        else:
            print(invalid_str)
        input("Press enter to continue.")


if __name__ == "__main__":
    main()
