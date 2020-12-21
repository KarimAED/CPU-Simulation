# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 17:55:42 2020

@author: Karim
"""
import numpy as np
import matplotlib.pyplot as plt


# very slow but gives sufficient results for heat sink analysis
def grid_search_hs(fun, l_params, ct_params, dist_params, t_thresh=80,
                   args=[], plot_best=True, plot_surface=False):
    """
    The function to perform a grid search over the heat sink parameters.
    It takes user input for the different ranges to check or alternatively
    defaults to predefined values.
    
    Also calculates heat sink area and adds it to the plot

    Parameters
    ----------
    fun : function
        The function wrapping the solver and returning the processor temp.
    l_params: int[]
        start, end and step for the fin lengths to loop over (in mm).
    ct_params: int[]
        start, end and step for the fin number to loop over.
    dist_params: int[]
        start, end and step for the fin distances to loop over (in mm).
    t_thresh: float, optional
        Temperature threshold below which the values in the grid plot will be
        marked green. The default is 80°C.
    args : list, optional
        Additional args such as error or show to pass to the function.
        The default is [].
    plot_best : bool, optional
        If the lowest temp result should be plotted. The default is True.
    plot_surface : bool, optional
        If a heatmap/surface plot of all the grid search options should be
        generated. The default is False.

    Returns
    -------
    min_vals : float[][]
        A list of the different temperatures, their parameters and area.
        Sorted ascending from lowest to highest temperature in the range.

    """
    
    
    # unpack values
    l_min, l_max, l_step = l_params
    ct_min, ct_max, ct_step = ct_params
    d_min, d_max, d_step = dist_params
    
    
    # Create ranges over which to perform the grid search
    lengths = np.arange(l_min, l_max, l_step)
    counts = np.arange(ct_min, ct_max, ct_step)
    dists = np.arange(d_min, d_max, d_step)
    
    # set up array to store all the temperatures in
    temps = []
    
    # get progress indicators
    i = 1
    m = lengths.size*counts.size*dists.size
    
    #A nested loop for each of the parameters
    for dist in dists:
        #2d array slices to fill up
        temp_2d = []
        for count in counts:
            #1d array slices to fill up
            temp_1d = []
            for l in lengths:
                print("\n\n\n\nGridsearch run {} of {}".format(i, m))
                
                #Adds the temperature returned to the array slice
                temp_1d.append(fun(int(l), int(count), int(dist), *args))
                i += 1
            temp_2d.append(temp_1d)
        temps.append(np.array(temp_2d))
    
    # to turn into ndarray
    temps = np.stack(temps)

    # Make the grid search heat maps
    if plot_surface:
        
        # One heat map for each fin dist value
        for i in range(temps.shape[0]):
            
            fig = plt.figure()
            
            # Set proper limits on the heat map and plot it
            ext = [lengths[0]-l_step/2, lengths[-1]+l_step/2,
                   counts[0]-ct_step/2, counts[-1]+ct_step/2]
            im = plt.imshow(temps[i], cmap="hot", extent=ext,
                               aspect="equal", interpolation="none")
            
            # Label axes and plot fin dist
            plt.title("Fin Separation: {}mm".format(dists[i]))
            plt.ylabel("Fin Number")
            plt.xlabel("Fin Length in mm")
            
            
            # Add area notation to each of the points
            for l in range(lengths.size):
                for c in range(counts.size):
                    
                    # get the area
                    area = lengths[l]*counts[c]*(1+dists[i])/1000
                    
                    #Conditional Formatting for points under temperature
                    #threshold
                    if temps[i, c, l] < t_thresh:
                        col = "g"
                    else:
                        col = "w"
                    
                    # Add the text
                    plt.text(lengths[l], counts[c], round(area, 2),
                       ha="center", va="center", color=col)
            
            # Add a reference for temperature levels and plot
            cb = plt.colorbar(im)
            cb.ax.set_ylabel("°C", rotation=0, va="bottom")
            fig.show()
    
    # for sorting the results
    min_vals = []
    
    # for all values in temps
    for i in range(temps.size):
        
        # Find new minimum
        m = np.amin(temps).copy()
        pos = np.where(temps==m)
        
        # Store minimum information
        vals = [int(lengths[pos[2]][0]),
                int(counts[pos[1]][0]),
                int(dists[pos[0]][0]),
                m]
        
        # Add area to the minimum info
        area = vals[0]*vals[1]*(1+vals[2])
        vals.append(area)
        
        # Set minimum value to inf to find next minimum
        temps[pos] = np.inf
        min_vals.append(vals)
    
    # Find the absolute minimum and plot it
    minimum = min_vals[0]
    if plot_best:
        fig = plt.figure()
        fun(*minimum[:-2], *args, show=True)

    
    return min_vals

    