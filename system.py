# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 18:52:11 2020

@author: Karim
"""

from dataclasses import dataclass, field
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpat



@dataclass(frozen=True) # frozen as it's properties should not change
class Material:
    """
    The dataclass that stores all of the information for a given material.
    It takes a name, power and thermal conductivity to fully define it's
    properties.
    """
    name: str
    power: float = 0.0
    conductivity: float = 0.0

    def invalid(self):
        """
        A function to enforce typing on all attributes of this dataclass.

        Returns
        -------
        ret : str[]
            A list of all attributes for which the typing failed. None is
            considered an acceptable value

        """
        ret = []
        
        # Get all of the fields names and definitions
        for f_name, f_def in self.__dataclass_fields__.items():
            
            # Get the type of the assigned value
            true_type = type(getattr(self, f_name))
            
            # If the type is wrong and the value is not None, mark for error
            if true_type != f_def.type and getattr(self, f_name) != None:
                ret.append(f_name)
        
        return ret
    
    def __post_init__(self):
        """
        A Dataclass function that gets called after the default constructor.
        Here used to ensure proper typing.

        Raises
        ------
        TypeError
            If the typing through invalid() method fails for any attribute.

        Returns
        -------
        None.

        """
        # get attributes with invalid type
        inv = self.invalid()
        
        # if inv not empty, raise error with the relevant field names
        if inv:
            raise TypeError("Attributes "+" ".join(inv)+" have invalid types.")
    
    def __str__(self):
        """
        Returns
        -------
        str
            A String representation of the material dataclass.

        """
        return "Material {} with thermal conductivity {} (W/mK) and power \
            generation {} (W/mm^3)" \
            .format(self.name,self.power, self.conductivity)
        
  

@dataclass(frozen=True) # frozen as it's values should not change
class HeatSink:
    """
    The dataclass encoding all of the information about the heat sink. It also
    creates further data such as width and height of the actual heat sink.
    
    The width of each fin is fixed at 1mm.
    """
    fin_len: int = 20
    fin_dist: int = 2
    fin_count: int = 10
    base_height: int = 4 # base height defaults to 4 and is usually unchanged
    width: int = None
    height: int = None
    
    def invalid(self):
        """
        A function to enforce typing on all attributes of this dataclass.

        Returns
        -------
        ret : str[]
            A list of all attributes for which the typing failed. None is
            considered an acceptable value

        """
        ret = []
        
        # Get all of the fields names and definitions
        for f_name, f_def in self.__dataclass_fields__.items():
            
            # Get the type of the assigned value
            true_type = type(getattr(self, f_name))
            
            # If the type is wrong and the value is not None, mark for error
            if true_type != f_def.type and getattr(self, f_name) != None:
                ret.append(f_name)
        
        return ret
    
    def __post_init__(self):
        """
        A Dataclass function that gets called after the default constructor.
        Here used to ensure proper typing and to calculate width and height
        of the HeatSink before freezing the data class.

        Raises
        ------
        TypeError
            If the typing through invalid() method fails for any attribute.

        Returns
        -------
        None.

        """
        #Evaluate width and height
        wd = self.fin_count*(self.fin_dist + 1)
        ht = self.fin_len + self.base_height
        
        #Still need to use this to change frozen class attributes
        super().__setattr__("width", wd)
        super().__setattr__("height", ht)
        
        #Find invalid types and raise error
        inv = self.invalid()
        if inv:
            raise TypeError("Attributes "+" ".join(inv)+" have invalid types.")
    
    def __str__(self):
        """
        Returns
        -------
        str
            A String representation of the material dataclass.

        """
        return "Heat sink (width: {}mm, height: {}mm) with {} fins of length \
            \n{}mm and distance {}mm".format(self.width, self.height, 
            self.fin_count, self.fin_len, self.fin_dist)



@dataclass()
class System:
    """
    A dataclass encoding the information about the whole system.
    Takes a name, a list of materials to construct it from, a resolution
    and optionally a heat sink to generate the layout. The shape of the ceramic
    holder and chip can also be varied if required.
    
    This assumes the following Material instances in materials:
            materials[0] - Air
            materials[1] - Processor
            materials[2] - Ceramic Case
            materials[3] - Heat Sink
    """
    name: str
    materials: list
    resolution: int = 2 # in pxl/mm
    heat_sink: HeatSink = None
    base_params: list = field(default_factory=list)
    layout: np.ndarray = None # gets generated automatically
    
    def invalid(self):
        """
        A function to enforce typing on all attributes of this dataclass.

        Returns
        -------
        ret : str[]
            A list of all attributes for which the typing failed. None is
            considered an acceptable value

        """
        ret = []
        
        # Get all of the fields names and definitions
        for f_name, f_def in self.__dataclass_fields__.items():
            
            # Get the type of the assigned value
            true_type = type(getattr(self, f_name))
            
            # If the type is wrong and the value is not None, mark for error
            if true_type != f_def.type and getattr(self, f_name) != None:
                ret.append(f_name)
        
        return ret
    
    def __post_init__(self):
        """
        The post init function is called after the constructor.
        It is used to check for wrong arguments and to generate the layout.

        Raises
        ------
        ValueError
            If the resolution is in the wrong range or the number of materials
            does not match the requirements (see help(System) for more info).
        TypeError
            If any of the typing through invalid() fails.

        Returns
        -------
        None.

        """
        
        # Resolution minimum at 2pxl/mm, to ensure both forward and central
        # difference schemes work properly
        if self.resolution < 2:
            raise ValueError("Resolution must be at lest 2 pxl/mm.")
        
        # 3 materials without heat sink and 4 with heat sink required
        if (len(self.materials)!=3 and not self.heat_sink) or \
        (len(self.materials)!=4 and self.heat_sink):
            raise ValueError("Wrong number of materials defined.")
        
        
        # Sets default parameters for shape of chip and ceramic case
        if not self.base_params:
            self.base_params = [14, # processor width
                                1,  # processor height
                                20, # case width
                                2]  # case height
        
        # Ensures that the parameters for the base are properly given
        if len(self.base_params)!=4:
            raise ValueError("Wrong number of base parameters given.")
        
        
        # generates the 2d numpy array representing the spatial layout of the
        # system
        self.generate_layout()
        
        # checks for invalid types
        inv = self.invalid()
        if inv:
            raise TypeError("Attributes "+" ".join(inv)+" have invalid types.")
    
    def generate_layout(self):
        """
        This function uses the heat sink and materials supplied to generate
        a spatial map of the system. It is called right after the constructor.
        
        There is a 1mm edge of air inserted around the system.
        
        This assumes the following Material instances:
            materials[0] - Air
            materials[1] - Processor
            materials[2] - Ceramic Case
            materials[3] - Heat Sink

        Returns
        -------
        None.

        """
        #====================================================================
        #This builds the layout in mm and then scales to match resolution
        #====================================================================
        
        # calculate the total height of the base
        base_height = self.base_params[1] + self.base_params[3]
        
        # to get the total height and width of the system + air gap
        if self.heat_sink:
            # with heat sink find maximum of either width of case or heat sink
            width = max(self.heat_sink.width+1, self.base_params[2]) + 2
            # height adds with the heat sink height
            height = base_height + self.heat_sink.height + 2
        
        # without a heat sink just add the 2mm for the air gap
        else:
            width = self.base_params[2] + 2
            height = base_height + 2
        
        # find the total shape and center of the system
        shape = (width, height)
        center = (shape[0]//2, shape[1]//2)
        
        # initialise empty layout
        layout = np.zeros(shape)

        
        # to keep track of the different regions along the height of the system
        x = 1
        
        # set up the processor region
        for i in range(x, self.base_params[1]+x):
            # centers it properly
            for j in range(center[0]-self.base_params[0]//2, \
                           center[0]+self.base_params[0]//2):
                # processor at index 1 in materials
                layout[j, i] = 1
        
        # end of processor
        x += self.base_params[1]
        
        # set up case region
        for i in range(x, self.base_params[3]+x):
            for j in range(center[0]-self.base_params[2]//2, \
                           center[0]+self.base_params[2]//2):
                # case at index 2 in materials
                layout[j, i] = 2
        
        # end of case
        x += self.base_params[3]
        
        
        # without heat sink, we are already done
        if self.heat_sink:
            
            # set up base of heat sink
            for i in range(x, self.heat_sink.base_height+x):
                for j in range(center[0]-self.heat_sink.width//2, \
                               center[0]+self.heat_sink.width//2):
                    # heat sink at index 3 in materials
                    layout[j, i] = 3
            
            # end of heat sink base
            x += self.heat_sink.base_height
            
            
            # generate fins of heat sink
            for i in range(x, self.heat_sink.fin_len+x):
                for j in range(center[0]-self.heat_sink.width//2, \
                               center[0]+self.heat_sink.width//2):
                    # find left-most point where fin starts
                    start = center[0]-self.heat_sink.width//2
                    
                    # the spatial period at which fins repeat
                    fin_period = self.heat_sink.fin_dist+1
                    # only set the parts where we are in the 1st mm of the
                    # period to be a part of the heat sink (the fin)
                    if (j-start)%(fin_period) < 1:
                        layout[j, i] = 3
        
        #====================================================================
        #Now convert from mm to pixels
        #====================================================================
        
        # create larger layout shape
        pxl_shape = (shape[0]*self.resolution, shape[1]*self.resolution)
        self.layout = np.zeros(pxl_shape)
        
        # loop over the larger layout
        for i in range(pxl_shape[0]):
            for j in range(pxl_shape[1]):
                
                # fill with corresponding elements in mm array
                self.layout[i, j] = \
                    layout[i//self.resolution, j//self.resolution]
        
        # typecast into int to use for indexing
        # then transpose to match desired array structure
        self.layout = self.layout.astype("int32").T
    
    
    def show(self):
        """
        A function to show the layout of the whole system as a heat map,
        by calling matplotlib.pyplot.imshow.

        Returns
        -------
        None.

        """
        
        # set up figure
        fig = plt.figure()
        
        # set up annotations
        plt.title(self.name)
        plt.xlabel("x in mm")
        plt.ylabel("y in mm")
        
        # get indices for materials for legend
        vals = np.arange(0, len(self.materials))
        
        # create x and y limits to scale the heat map
        x_y = np.array(self.layout.shape)/self.resolution
        ext = [0, x_y[1], 0, x_y[0]]
        
        # create the heat map
        img = plt.imshow(self.layout, extent=ext, origin="upper")
        
        # match the colors to the materials
        col = [img.cmap(img.norm(val)) for val in vals]
        
        # create patches to place in legend and show legend
        pat = [mpat.Patch(color=col[i], label=self.materials[i].name) \
               for i in vals]
        plt.legend(handles=pat, bbox_to_anchor=(1., 1),
                   loc=2, borderaxespad=0.)
        fig.show()
    
    def __str__(self):
        """
        Returns
        -------
        str
            A String representation of the material dataclass.

        """
        string = ""
        for i in self.materials:
            string += "\n"+repr(i)
        return "{} with resolution {}. Materials are:\n\n {}\n\nHeat Sink:\n{} \
                ".format(self.name, self.resolution,
                        string, str(self.heat_sink))