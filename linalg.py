# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 22:28:02 2020

@author: Karim
"""

import numpy as np


def jac_solve(mat, vec, e=1e-5):
    """
    Implements the jacobian method to slove matrix equations.

    Parameters
    ----------
    mat : np.ndarray
        The matrix in the equation.
    vec : np.array
        The vector in the equation.
    e : float, optional
        The convergence error after which the equation is considered solved.
        The default is 1e-5.

    Returns
    -------
    x_1 : np.array
        The vector solving the equation.

    """
    
    m = np.copy(mat)

    diag_inv = []
    
    # invert the diagonal and remove it from the matrix
    for i in range(len(vec)):
        temp = np.zeros(m.shape[1])
        temp[i] = 1/m[i, i]
        diag_inv.append(temp)
        m[i, i] = 0
    
    diag_inv = np.stack(diag_inv)
    
    # To ensure no divide by zero initially
    x = np.ones(len(vec))
    x_1 = np.zeros(len(vec))+1e-10
    
    # get initial error estimate
    err = np.abs(np.sum(x_1)-np.sum(x))/np.abs(np.sum(x))
    while err > e:
        
        # find new guess and keep copy from previously
        x = np.copy(x_1)
        x_1 = -diag_inv.dot(m.dot(x))+diag_inv.dot(vec)
        err = np.abs(np.sum(x_1)-np.sum(x))/np.abs(np.sum(x))
    return x_1


def lu(mat):
    """
    Performs the LU-decomposition of a matrix.

    Parameters
    ----------
    mat : np.ndarray
        The matrix to decompose.

    Raises
    ------
    ValueError
        Determinant of mat is 0.
    Exception
        Self validation failed.

    Returns
    -------
    l : np.ndarray
        Lower diagonal matrix.
    u : np.ndarray
        Upper diagonal matrix.

    """
    l = np.zeros(mat.shape)
    u = np.zeros(mat.shape)
    
    # Doolittle choice
    for i in range(mat.shape[0]):
        l[i, i] = 1
    
    # Crouts algorithm
    for j in range(mat.shape[1]):
        for i in range(j+1):
            s = 0
            for k in range(i):
                s += l[i, k]*u[k, j]
            u[i, j] = mat[i, j]-s
        
        for i in range(j+1, mat.shape[0]):
            s = 0
            for k in range(j):
                s += l[i, k]*u[k, j]
            l[i, j] = mat[i, j]-s
            if u[j, j] == 0:
                raise ValueError("Det(mat) = 0, cannot divide by 0.")
            l[i, j] = l[i, j]/u[j, j]
    
    # check if results make sense
    if not np.all(np.abs(l.dot(u)-mat) < 1e-10):
        raise Exception("LU cross-validation failed.")
    
    return l, u


def f_sub(mat, b):
    """
    Implements forward substitution of a lower diagonal matrix.

    Parameters
    ----------
    mat : np.ndarray
        lower diagonal matrix.
    b : np.array
        vector to solve for.

    Raises
    ------
    Exception
        Self validation failed.

    Returns
    -------
    x : np.array
        The result of the forward substitution.

    """
    x = np.array([])
    
    # row-wise forward substitution
    for i, b_v in enumerate(b):
        x_new = b_v - mat[i, :i].dot(x)
        x = np.append(x, [x_new])
    # check if the results match
    if not np.all(np.abs(mat.dot(x)-b) < 1e-10):
        raise Exception("Forward substitution failed.")
    return x


def b_sub(mat, b):
    """
    Implements backwards substitution of an upper diagonal matrix.

    Parameters
    ----------
    mat : np.ndarray
        upper diagonal matrix.
    b : np.array
        vector to solve for.

    Raises
    ------
    ValueError
        Matrix determinant is 0
    Exception
        Self validation failed.

    Returns
    -------
    x : np.array
        The result of the backward substitution.

    """
    x = np.array([])
    
    # row-wise backwards substitution
    for i in range(len(b)):
        j = len(b)-i-1
        x_new = b[j] - mat[j, j+1:].dot(x)
        if mat[j, j] == 0:
            raise ValueError("Cannot perform back-substitution with matrix \
                             with zero determinant.")
        x_new /= mat[j, j]
        x = np.append([x_new], x)
    
    # check if the results make sense
    if not np.all(np.abs(mat.dot(x)-b) < 1e-10):
        raise Exception("Backward substitution failed.")
    return x


def lu_solve(mat, b):
    """
    combines LU-decomposition, forward and backward substitution to
    solve a matrix equation.

    Parameters
    ----------
    mat : np.ndarray
        Matrix to decompose.
    b : np.array
        Vector to solve for.

    Raises
    ------
    Exception
        Self validation failed.

    Returns
    -------
    x : np.array
        Vector that solves the equation.

    """
    # set up L and U matrix
    l, u = lu(mat)
    
    # first forward substitute with L and the RHS
    y = f_sub(l, b)
    # Then use that in backwards substitution with U
    x = b_sub(u, y)
    
    # check for consistency
    if not np.all(np.abs(mat.dot(x)-b) < 1e-10):
        raise Exception("LU solving failed.")
    return x


def inv(sol, mat):
    """
    Generalised inversion method. Takes a solver function as an argument and
    constructs the inverse with it.

    Parameters
    ----------
    sol : function
        A method that solves a matrix equation.
    mat : np.ndarray
        The matrix to invert.
    
    Raises
    ------
    Exception
        self validation failed.

    Returns
    -------
    np.ndarray
        The inverse of mat.

    """
    # set up identity matrix
    i = np.identity(mat.shape[0])
    x = []
    
    # solve for x column by column
    for col in i:
        x.append(sol(mat, col))
    
    # make ndarray
    inv = np.stack(x)

    # check for consistency
    if np.any(mat.dot(inv)-i > 1e-10):
        raise Exception("Matrix inversion failed")
    
    return np.stack(x)
