# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 18:56:05 2020

@author: Karim
"""

from abc import ABC, abstractmethod
import time as t

import numpy as np
import matplotlib.pyplot as plt

import linalg
from system import System
from conv import h

def neighbours(i_j):
    """
    
    A function to find the coordinate neighbours of a 2 coordinate Tupel.

    Parameters
    ----------
    i_j : int[2]
        Tupel containing a given elements coordinates.

    Returns
    -------
    list[int[2]]
        List of tupels containing all neighbouring element coordinates.

    """
    i, j = i_j
    return [(i-1, j), # bottom 
            (i, j-1), # left
            (i+1, j), # top
            (i, j+1)] # right




class Solver(ABC):
    """
    Abstract base class implementing the capability all other solvers have in
    common. Distinguishes between Natural and Forced convection on a class
    level.
    
    Generates the weighting stencil matrix, maps between 2d and 1d and 
    finds boundary locations.
    
    Initialises temperature across a given system and knows how to plot itself.
    
    Temperature map will first be in relative temperatures but after solving
    be returned in °C.
    """

    # This Class Attribute is used to determine the convection conditions a 
    # Solver should use.
    modes = ("NATURAL", "FORCED", "NAT", "FOR", "N", "F")

    @classmethod
    def check_mode(cls, mode):
        """
        Class Method checking if mode is valid.
        
        Parameters
        ----------
        mode: str
            string identifying the convection mode for which we want to solve.
        
        Returns
        ----------
        None
        
        Raises
        ----------
        TypeError
            If mode is not 'natural' or 'forced' identifier.
        
        """
        if mode.upper() not in cls.modes:
            raise ValueError("Mode must be 'NATURAL' or 'FORCED' identifier.")



    def __init__(self, sys, mode="NATURAL", velocity=0):
        """
        The solver constructor.
        
        Parameters
        ----------
        sys : system.System
            The system for which we want to solve the steady-state temperature.
        mode : str, optional
            Which kind of convection applies. The default is "NATURAL".
        velocity : float, optional
            The wind speed for forced convection. The default is 0.

        Raises
        ------
        TypeError
            If sys is not a System.
        ValueError
            If mode and velocity specifications don't match.

        Returns
        -------
        None.

        """
        
        # check if the mode used is valid, i.e. natural or forced
        self.check_mode(mode)

        # check type of system
        if not isinstance(sys, System):
            raise TypeError("sys is not a System.")
        
        # Match velocity and mode conditions (Velocity only works with forced
        # convection)
        if velocity != 0 and mode not in ("FORCED", "FOR", "F"):
            raise ValueError("v!=0 allowed only in mode 'FORCED'.")
        elif mode in ("FORCED", "FOR", "F") and velocity <= 0:
            raise ValueError("v must be positive in mode 'FORCED'.")
        
        # store the parameters
        self.sys = sys
        self.mode = mode
        self.velocity = velocity
        
        # select the region which is the processor
        self.core = np.where(sys.layout == 1, True, False)
        
        # create the empty array used as a 2d to 1d map
        self.map = []
        
        # calculate the step size from the resolution
        self.h = 1/sys.resolution
        
        # generate the 2d to 1d map
        self.gen_map()
        
        # set up the initial temperature map (in T-T_air)
        self.init_temp()
        
        # generate the matrix and RHS vector whose solution is the temperature
        # map of the system (or variation thereof).
        self.gen_lin_alg()
        
        # used in the relaxation solver to store the solution matrix inverse
        self.inv = None

    def gen_map(self):
        """
        This generates the 2d to 1d mapping for non-air points in the 2d
        system.

        Returns
        -------
        list[int[2]]
            1d list of 2d-coordinate tupels used as a mapping for the matrix
            eqn.

        """
        # for all points on 2d grid, add coordinate tupel to 1d mapping array
        for i in range(self.sys.layout.shape[0]):
            for j in range(self.sys.layout.shape[1]):
                if self.sys.layout[i, j] != 0:
                    self.map.append((i, j))
        return self.map

    def init_temp(self):
        """
        Initialises the heat map for the system as all zeros (processor 
        temperature matching ambient temperature).

        Returns
        -------
        float[][]
            2d heat-map containing all zeros.

        """
        # initialise at zero temperature relative to the air
        self.temp = np.zeros(self.sys.layout.shape)
        return self.temp


    def gen_lin_alg(self):
        """
        This generates the linear algebra components required for all methods.
        This includes the matrix of weights for the curvature stencil as well
        as the steady-state heat equation power conditions.
        It also finds all of the boundary points of the system and stores them
        for easy access.
        
        The matrix includes the thermal conductivity values and averages over
        boundaries. This allows it to deal with non-homogenous k and interfaces
        with varying k.

        Returns
        -------
        None.

        """
        
        #=====================================================================
        #The matrix is designed to include the thermal conductivities.
        #This allows it to average thermal conductivities and appropriately
        #scale the corresponding contributions to the matrix.
        #=====================================================================
        
        # set up empty arrays to fill up
        self.mat = []
        self.outer_nb_infos = []
        self.b_vector = []

        # for each point in the system
        for i in range(len(self.map)):

            # get the position and material and store to use when looking for
            # outer edges
            pos = self.map[i]
            matter = self.sys.materials[self.sys.layout[pos]]
            outer_nb_info = [pos, matter]
            
            # calculate the power value for the RHS vector and append it
            b_value = -matter.power*self.h*self.h
            self.b_vector.append(b_value)
            
            # initialise matrix row at zeros, as it will be a sparse matrix
            mat_row = np.zeros(len(self.map))
            
            # consider all of the neighbours to the point
            for nb in neighbours(pos):
                
                # if the neighbour is not air
                if nb in self.map:
                    
                    # get it's material
                    nb_material = self.sys.materials[self.sys.layout[nb]]
                    
                    # calculate mean conductivity and store it as neighbour
                    # conductivity
                    k = (nb_material.conductivity+matter.conductivity)/2
                    mat_row[self.map.index(nb)] = k
                
                # if the neighbour is air
                else:
                    # assume ficticiousextension with the same material as
                    # here and store that.
                    k = matter.conductivity
                    outer_nb_info.append(nb)
                
                # each of the neighbours contributes once to the matrix value
                # on the diagonal, either with the k at pos or the avg. val
                mat_row[i] -= k
            
            # fill up arrays
            self.mat.append(mat_row)
            self.outer_nb_infos.append(outer_nb_info)
        
        # return proper ndarrays
        self.mat = np.stack(self.mat)
        self.b_vector = np.array(self.b_vector)

    def gen_boundary_temp(self, ij, material, deriv=False,
                          cds_outer_pos=None):
        """
        Generates the ficticious temperature for points on the air side of the
        boundary at position ij with material "material". Can also generate
        derivatives used for the Newton-Raphson method.
        
        Uses the central or forward difference scheme, dependent on the param-
        eter cds_outer_pos.
        
        DEPRECATED: central difference scheme.

        Parameters
        ----------
        ij : int[2]
            The x and y positions whose neighbours temperature is to be eval-
            uated.
        material : system.Material
            The material at location ij.
        deriv : bool, optional
            If the derivative with respect to the temperatures should be
            evaluated instead (For NewtonSolver). The default is False.
        cds_outer_pos: int[2], optional
            If the cds_outer_pos is given, the central difference scheme will
            be used, otherwise it uses the forward difference scheme.
            The default is None.

        Returns
        -------
        float or float[]
            The fictitious temperature at the neighbouring point of ij.
            For deriv=True and cds_outer_pos!=none returns two derivatives.

        """
        # if the forward difference scheme is used
        if not cds_outer_pos:
            # get the convection contribution across the boundary
            grad = h(self.temp[ij], self.velocity)*self.h/material.conductivity
            
            # The ficticious temperature according to forward difference scheme
            if not deriv:
                return material.conductivity*self.temp[ij]*(1-grad)
            
            # The derivatives are used in the Newton-Raphson method
            elif deriv and not self.velocity:
                # The derivative for natural convection
                return material.conductivity*(1-4*grad/3)
            
            elif deriv and self.velocity:
                # The derivative for forced convection
                return material.conductivity*(1-grad)
        
        # for the central difference scheme
        else:
            # find the inner neighbour of the point
            d_i = ij[0]-cds_outer_pos[0]
            d_j = ij[1]-cds_outer_pos[1]
            ij_in = (ij[0]+d_i, ij[1]+d_j)
            
            # the convection contribution is now scaled            
            grad = 2*h(self.temp[ij], self.velocity)*self.h \
                    / material.conductivity
            
            # return ficticious temperature with central difference scheme
            if not deriv:
                return material.conductivity* \
                    (self.temp[ij_in]-self.temp[ij]*grad)
            
            # for central difference scheme we have two partial derivatives
            # for the point itself and its inner neighbour
            elif deriv and not self.velocity:
                # for natural convection, returns [point_deriv, nb_deriv]
                return -material.conductivity*4*grad/3, material.conductivity
            
            
            elif deriv and self.velocity:
                # for natural convection, returns [point_deriv, nb_deriv]
                return -material.conductivity*grad, material.conductivity

    def show(self, err):
        """
        A function allowing the solver to show the generated heatmap.
        
        Parameters
        ----------
        
        err: float
            The precision to which the current heat map has been evaluated.

        Returns
        -------
        None.

        """
        # Set up heat map figure
        fig = plt.figure()
        plt.title(self.sys.name+" Heat Map")
        plt.xlabel("x in mm")
        plt.ylabel("y in mm")
        
        #adjust axes
        x_y = np.array(self.sys.layout.shape)/self.sys.resolution
        ext = [0, x_y[1], 0, x_y[0]]
        
        # show plot and add temperature reference
        im = plt.imshow(self.temp, cmap="hot", extent=ext)
        cb = plt.colorbar(im)
        cb.ax.set_ylabel("°C", rotation=0, va="bottom")
        
        # add label giving processor temperature and relative error
        label = r"$T_{p}$ = "+str(round(np.mean(self.temp[self.core]), 2)) \
                + r"$\pm$"+str(err)+"°C"
        plt.text(2, 2, 
                 label,
                 bbox=dict(facecolor='white', alpha=0.8))
        fig.show()

    @abstractmethod
    def solve(self):
        """
        An abstract method overwritten by all children solvers.

        Returns
        -------
        None.

        """
        pass




class ExactSolver(Solver):
    """
    DEPRECATED
    This solver uses the fact that for FD-scheme BCs and forced convection, 
    we can formulate an exact linear algebra equation which solves the system.
    
    This is equivalent to using the Newton solver, but less versatile,
    hence deprecated.
    """

    def __init__(self, sys, mode="FORCED", v=0):
        """
        Constructor for ExactSolver, calls Solver constructor, after checking
        for the only allowed mode: 'FORCED'. Central difference scheme also
        not available for this solver.

        Parameters
        ----------
        sys : system.System
            The system for which we want to solve the steady-state temperature.
        mode : str, optional
            Which kind of convection applies. The default is "NATURAL".
        velocity : float, optional
            The wind speed for forced convection. The default is 0.

        Raises
        ------
        TypeError
            If sys is not a System.
        ValueError
            If mode is not 'FORCED'.

        Returns
        -------
        None.

        """
        if mode.upper() not in ("FORCED", "FOR", "F"):
            raise ValueError("Exact Solver only works with forced convection.")
        super().__init__(sys, mode, v)

    def solve(self, show=False, kp_linalg=False, err=1e-7):
        """
        This is the method solving the system temperature for the 'FORCED'
        convection case. It relies on the FD-Scheme for Boundary conditions.

        Parameters
        ----------
        show : bool, optional
            If the system is to be plotted after solving. The default is False.
        err: float, optional
            Manual override for the threshold to which accuracy the solver con-
            verges.
        kp_linalg: bool, optional
            If karim's python linear algebra should be used.
            The default is False.

        Returns
        -------
        float.
            The mean steady-state temperature across the processor.

        """
        start = t.time()
        mat = np.copy(self.mat)
        
        # loop over all the points (the RHS vector)
        for i in range(len(self.map)):
            pos = self.map[i]
            matter = self.sys.materials[self.sys.layout[pos]]
            # For each neighbour that is air, add a forward difference scheme
            # forced convection contribution
            for _ in self.outer_nb_infos[i][2:]:
                mat[i, i] -= h(0, self.velocity)*self.h
                mat[i, i] += matter.conductivity
        
        # personally implemented linear algebra (slower)
        if kp_linalg:
            # This needs much higher accuracy as its a non-iterative method,
            # cannot tolerate initial errors
            T = linalg.jac_solve(mat, self.b_vector, err)
        
        # numpy linear algebra (orders of magnitude faster)
        else:
            T = np.linalg.solve(mat, self.b_vector)
        
        # update the solver temperatures based on the solution
        for i in range(len(self.map)):
            self.temp[self.map[i]] = T[i]
        
        end = t.time()
        dt = end-start
        
        # turn relative temperatures into °C
        self.temp += 20 # return to °C
        
        # plot the heat map if required
        if show:
            self.show(err)
        print("\n\nRuntime: {}\
              \nFinal Processor Temp: {} / Final Error: {}".format(
              dt, np.mean(self.temp[self.core]), err))
        
        return np.mean(self.temp[self.core])




class RelaxSolver(Solver):
    """
    DEPRECATED
    A child of the Solver class implementing solving the heat equation by
    iterative relaxation. It works for both 'NATURAL' and 'FORCED' convection, but is
    much slower than the NewtonSolver and is hence deprecated.
    """

    def __init__(self, sys, mode="NATURAL", v=0):
        """
        This calls the Solver constructor.

        Parameters
        ----------
        sys : system.System
            The system for which we want to solve the steady-state temperature.
        mode : str, optional
            Which kind of convection applies. The default is "NATURAL".
        velocity : float, optional
            The wind speed for forced convection. The default is 0.

        Raises
        ------
        TypeError
            If sys is not a System.
        ValueError
            If mode and velocity specifications don't match.

        Returns
        -------
        None.

        """
        super().__init__(sys, mode, v)

    def update_temp(self, cds=False):
        """
        The method updating the temperature according to the relaxation
        algorithm once.
        
        cds option discouraged (see report).
        
        Parameters
        ---------
        cds: bool, optional
            if the central difference scheme should be used instead of the
            forward difference scheme. default is False.

        Returns
        -------
        None.

        """
        
        # get a copy of the constant parts (power) of the RHS vector
        b = np.copy(self.b_vector)
        
        # add contributions for all points with air neighbours
        for i in range(len(self.outer_nb_infos)):
            
            # for each air neighbour
            for nb in self.outer_nb_infos[i][2:]:
                # fake boundary temperature according to fds
                if not cds:
                    boundary_t = \
                        self.gen_boundary_temp(self.outer_nb_infos[i][0], 
                                               self.outer_nb_infos[i][1])
                # fake boundary temperature according to cds
                else:
                    boundary_t = \
                        self.gen_boundary_temp(self.outer_nb_infos[i][0], 
                                               self.outer_nb_infos[i][1],
                                               cds_outer_pos=nb)
                # add (subtract) the fake temperature contributions.
                b[i] -= boundary_t
        
        # use the previously determined matrix inverse to calculate the temp
        T = self.inv.dot(b)
        
        # update the solver 2d temperature array
        for i in range(len(self.map)):
            self.temp[self.map[i]] = T[i]


    def solve(self, err, init_temp=0, verbose=False, show=False, 
              kp_linalg=False, cds=False):
        """
        The iterative method for solving the heat equation across the system
        with the relaxation method.
        Convergence criterion is the updating step across the area of the 
        system with material 1 (processor). 200 additional iterations ensure
        robustness against extrema due to overshooting.
        
        Initial temperatures can improve convergence if the guess is good.

        Parameters
        ----------
        err : float
            The mean temperature difference across processor between two
            iterations which is considered sufficient for convergence.
        init_temp : float, optional
            The temperature to which the non-air components are initialised.
            The default is 0.
        verbose : bool, optional
            If True, print status every 100 iterations. The default is False.
        show : bool, optional
            If True, plot system temp after solving. The default is False.
        kp_linalg: bool, optional
            If karim's python linear algebra should be used.
            The default is False.
        cds: bool, optional, DEPRECATED(Not Recommended, see report)
            If the central difference scheme should be used instead of the
            forward difference scheme. The default is False.

        Returns
        -------
        float.
            The mean steady-state temperature across the processor.

        """
        
        #=====================================================================
        #For this method it is fastest to first invert the matrix, as the
        #same matrix is always used over and over again.
        #=====================================================================
        
        if kp_linalg:
            # manually implemented solver, uses LU-decomposition and inversion
            self.inv = linalg.inv(linalg.jac_solve, self.mat)
        else:
            # exact numpy inversion
            self.inv = np.linalg.inv(self.mat)
            
        # keep a copy of temperature to use as convergence reference
        t_0 = self.temp
        
        # initial temperatures improve convergence in some cases
        for i in self.map:
            self.temp[i]=init_temp
            
        # initialise one higher to start iterating
        t_1 = t_0+1
        i = 0
        j = 0 # use to prevent convergence to an extremum
        
        # get the initial error
        e = np.mean(np.abs(t_1[self.core]-t_0[self.core]))

        start_time = t.time()

        # even if updating error low, keep going for 200 iters to avoid
        # convergence
        while j < 200:
            
            # keep copies and update temperature
            t_0 = np.copy(t_1)
            self.update_temp(cds)
            t_1 = self.temp

            i += 1
            
            #get error
            e = np.mean(np.abs(t_1[self.core]-t_0[self.core]))
            
            # user output every 100 iterations as this iterates A LOT
            if verbose and i%100==0:
                print("Iteration: {} \t Updating/Error: {} \t Max Temp: {}"
                      .format(i, round(e, 3), round(np.max(t_0)+20, 3)))

            # start counting up j once error threshold is reached
            if e < err:
                j += 1
            
            # if we do diverge again, reset j to 0
            else:
                j = 0
        
        # calculations done in relative temperature, results returned in °C
        self.temp += 20 # Return to °C
        dt = t.time()-start_time
        
        # show heat map if required
        if show:
            self.show(err)
        print("\n\nRuntime: {} / Iterations: {} / Time/Iteration: {} \
              \nFinal Processor Temp: {} / Final Updating Error: {}".format(
              dt, i, dt/i, np.mean(self.temp[self.core]), e))
        
        return np.mean(self.temp[self.core])
        


class NewtonSolver(Solver):
    """
    Child class of Solver implementing the n-dimensional newton raphson method
    to solve for the system temperature.
    Equivalent to the ExactSolver for 'FORCED' convection and faster than the
    relaxation method for either convection mode.
    Works with FD-difference scheme only.
    """
    
    def __init__(self, sys, mode="NATURAL", v=0):
        """
        Internally calls the Solver constructor.

        Parameters
        ----------
        sys : system.System
            The system for which we want to solve the steady-state temperature.
        mode : str, optional
            Which kind of convection applies. The default is "NATURAL".
        velocity : float, optional
            The wind speed for forced convection. The default is 0.

        Raises
        ------
        TypeError
            If sys is not a System.
        ValueError
            If mode and velocity specifications don't match.

        Returns
        -------
        None.

        """
        super().__init__(sys, mode, v)
    
    
    def update_temp(self, kp_linalg=False, cds=False):
        """
        The method updating the temperature according to the newton-raphson-
        method once.
        
        Parameters
        ----------
        kp_linalg: bool, optional
            If karim's python linear algebra should be used.
            The default is False.
        cds: bool, optional, DEPRECATED(Not Recommended, see report)
            If the central difference scheme should be used instead of the
            forward difference scheme. The default is False.

        Returns
        -------
        None.

        """
        
        #=====================================================================
        #This method first creates a vector equation that equals zero,
        #then considers it's partial derivatives and uses that to gain better
        #estimates on the solution.
        #
        #Here, A will refer to the derivative matrix and f to the matrix
        #equation that should equal zero.
        #=====================================================================
        
        temp_vec = []
        
        
        # copy the matrix and the b vector initally
        A = np.copy(self.mat)
        b = np.copy(self.b_vector)
        
        
        # we modify the curvature stencil matrix and the RHS vector with power
        # contributions to get f and A
        
        # loop over all non-air points
        for i in range(len(self.map)):
            
            # find the position and material 
            pos = self.map[i]
            matter = self.sys.materials[self.sys.layout[pos]]
            
            # add the temperatures to the current temperature vector estimate
            temp_vec.append(self.temp[pos])
            
            # for each outer neighbour
            for nb in self.outer_nb_infos[i][2:]:
                
                # for forward difference scheme
                if not cds:
                    # add derivative contributions to diagonal of deriv. matrix
                    A[i, i] += self.gen_boundary_temp(pos, matter, deriv=True)
                    
                    # get the ficticious temperature at the boundary point
                    boundary_t = \
                        self.gen_boundary_temp(self.outer_nb_infos[i][0], 
                                               self.outer_nb_infos[i][1])
                
                # for cds, add derivative contributions to diagonal and inner
                # (opposite to the air point) neighbours
                else:
                    
                    # find inner neighbour position and index
                    diff = (pos[0]-nb[0], pos[1]-nb[1])
                    in_nb = (pos[0]+diff[0], pos[1]+diff[1])
                    ind = self.map.index(in_nb)
                    
                    # get the derivatives
                    derivs = self.gen_boundary_temp(pos, matter, deriv=True,
                                                    cds_outer_pos=nb)
                    
                    # update the diagonal and inner neighbour pos in the deriv
                    # matrix
                    A[i, i] += derivs[0]
                    A[i, ind] += derivs[1]
                    
                    # generate the boundary temperature (ficticious)
                    boundary_t = \
                        self.gen_boundary_temp(self.outer_nb_infos[i][0], 
                                               self.outer_nb_infos[i][1],
                                               cds_outer_pos=nb)
                
                # add ficticious temperature contributions to RHS vector copy
                b[i] -= boundary_t
                
        # get the current temperature estimates as an array
        temp_vec = np.array(temp_vec)
        
        # create f, must equal 0 to find soln.
        f = self.mat.dot(temp_vec)-b
        
        # Solve for difference between current and better temperature estimate
        if kp_linalg:
            # manual, less accurate jacobian solver
            delta_T = linalg.jac_solve(A, -f)
        else:
            # numpy, exact and fast solver
            delta_T = np.linalg.solve(A, -f)
        
        # add updates to the solver heat map
        for i in range(len(self.map)):
            self.temp[self.map[i]] += delta_T[i]
        
    
    def solve(self, err, verbose=False, show=False,
              kp_linalg=False, cds=False):
        """
        The iterative method for solving the heat equation across the system
        with the newton-raphson method.
        Convergence criterion is the updating step across the area of the 
        system with material 1 (processor). 

        Parameters
        ----------
        err : float
            Updating average across the material 1 (processor) considered 
            sufficient to terminate the algorithm. Values around 1e-3
            recommended.
        verbose : bool, optional
            If output should be generated for every iteration.
            The default is True.
        show : bool, optional
            If the temperature should be plotted after solving.
            The default is False.
        kp_linalg: bool, optional
            If karim's python linear algebra should be used.
            The default is False.
        cds: bool, optional, DEPRECATED(Not Recommended, see report)
            If the central difference scheme should be used instead of the
            forward difference scheme. The default is False.

        Returns
        -------
        float.
            Mean steady-state temperature across the processor.

        """
        
        # formatted output
        if verbose:
            print("Solving the System:\n\n", self.sys)
            print("\n\nSetting up...")
        
        # keep temperature copies for convergence reference
        t_0 = self.temp
        t_1 = t_0+1
        i = 0
        
        # get initial error (1 by means of initialisation)
        e = np.mean(np.abs(t_1[self.core]-t_0[self.core]))

        start_time = t.time()
        if verbose:
            print("Done! Now starting to iterate...")
        
        # no overshooting guard needed here, this always converges well
        while e > err:
            
            # update temperatures and keep references
            t_0 = np.copy(t_1)
            self.update_temp(kp_linalg, cds)
            t_1 = self.temp
            
            # calculate updating error
            e = np.mean(np.abs(t_1[self.core]-t_0[self.core]))
            
            i += 1
            
            if verbose:
                print("Iteration: {} \t Updating/Error: {} \t Max Temp: {}"
                      .format(i, round(e, 3), round(np.max(t_0)+20, 3)))
        
        # return temperatures in °C and not relative values
        self.temp += 20 # return to °C
        
        dt = t.time() - start_time
        
        
        # plot heat map if required
        if show:
            self.show(err)
        
        
        print("\n\nRuntime: {} / Iterations: {} / Time/Iteration: {} \
              \nFinal Processor Temp: {} / Final Updating Error: {}".format(
              dt, i, dt/i, np.mean(self.temp[self.core]), e))
        
        return np.mean(self.temp[self.core])