#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 13:04:31 2021

@author: Mohammad Asif Zaman

"""

import numpy as np
import pylab as py

from parameters import *
from geometry_def_DEP import *
from force_DEP import *


# Repeated definitions from force_DEP. Remove these later
# Parameters:
# =============================================================================

# Electrode array geometry parameters (units in um)
elec_width = 15
elec_spacing = 50
# =============================================================================


# Set particle colors
def particle_color(p):
    cl_p = cl_dred
 
    return cl_p







def draw_static_geo(ax_xy, ax_yz, ax_xz):
    Np, ro, tfinal, x_lim, y_lim, z_lim, xi_lim, yi_lim, zi_lim = time_pos_ax_limits()
    py.gcf().sca(ax_xy)
    substrate_xy = py.Rectangle((x_lim[0], y_lim[0]),(x_lim[1]-x_lim[0]), (y_lim[1]-y_lim[0]),fc=cl_lgrey)
    py.gca().add_patch(substrate_xy)
    for kk in range(-2,3):
        rectangle = py.Rectangle((x_lim[0]/2, -elec_width/2+kk*elec_spacing),x_lim[1],elec_width,fc=cl_royal_blue)
        py.gca().add_patch(rectangle)
        # ax.add_patch(rectangle)
    
    py.gcf().sca(ax_yz)
    substrate_yz = py.Rectangle((y_lim[0], z_lim[0]),(x_lim[1]-x_lim[0]), abs(z_lim[0]),fc=cl_dgrey, ec='k')
    py.gca().add_patch(substrate_yz)
     
    py.gcf().sca(ax_xz)
    substrate_xz = py.Rectangle((x_lim[0], z_lim[0]),(x_lim[1]-x_lim[0]), abs(z_lim[0]),fc=cl_dgrey, ec='k')
    py.gca().add_patch(substrate_xz)
   
    
   
    # Draw fluidic/reflecting geometries
    py.gcf().sca(ax_xy)
    w_type = geo_element_types()
    
    # Loop over all types of geometry segments
    for wt in range(w_type.shape[0]):
       
        fname = 'geo_points' + str(w_type[wt]) + '()'  # name of the corresponding geometry function
       
        # Evaluate geo_points*()function by string name
        p_w,z_w = eval(fname)
            
        Nwalls = p_w.shape[0]
        # Draw all sements of a particular type of geometry
        for m in range(Nwalls):
            
             f_wall_xy = py.Polygon(p_w[m],closed = True, fc = cl_army_green, ec = 'k')        
             py.gca().add_patch(f_wall_xy)
    
    return 0

    
    

# Draw electrodes for the animation
def draw_geo(tm, ax_xy, ax_yz, ax_xz):
    # tm is a dummy argument for this case (to make it consistent with other force functions)
    
    # The flag_source_state variable is used to draw/erase the source geometry only once
    # This is necessary to speed up the animation.
    
    if 'str1' not in globals():
        global str1
        global str2
        str1 = ''
        str2 = ''
    
    
    r_active, str1 = active_electrode(tm)   # call the active_electrode function to populate the str1 global variable
    
    if 'flag_source_state' not in globals():
        global flag_source_state     # Make this variable global so that the assigned value remains saved globally as t changes
        flag_source_state = 0        # initialize with OFF state
        print('Defining global flag for source geometry \n')
        
    
    
    
    if flag_source_state == 0:
        print('Drawing static geometry ...  ')
        
        draw_static_geo(ax_xy, ax_yz, ax_xz)
                
        print('Done \n')
        flag_source_state = 1

    
    return str1, str2



