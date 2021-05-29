#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 13:04:31 2021

@author: Mohammad Asif Zaman

"""

import numpy as np
import pylab as py

from parameters import *
from geometry_def_OT import *

from force_OT import *


# Set particle colors
def particle_color(p):
    cl_p = cl_dred
 
    return cl_p







def draw_static_geo(ax_xy, ax_yz, ax_xz):
     Np, ro, tfinal, x_lim, y_lim, z_lim, xi_lim, yi_lim, zi_lim = time_pos_ax_limits()
     py.gcf().sca(ax_yz)
     substrate_yz = py.Rectangle((y_lim[0], z_lim[0]),(y_lim[1] - y_lim[0]), abs(z_lim[0]),fc='#d4d4d4', ec='k')
     py.gca().add_patch(substrate_yz)
    
     py.gcf().sca(ax_xz)
     substrate_xz = py.Rectangle((x_lim[0], z_lim[0]),(x_lim[1] - x_lim[0]), abs(z_lim[0]),fc='#d4d4d4', ec='k')
     py.gca().add_patch(substrate_xz)
    
     py.gcf().sca(ax_xy)
     substrate_xy = py.Rectangle((x_lim[0], y_lim[0]), (x_lim[1] - x_lim[0]), (y_lim[1] - y_lim[0]),fc='#f9f9f9')
     py.gca().add_patch(substrate_xy)
     
     
     
     # Draw fluidic/reflecting geometries
    
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
 
    

def draw_dynamic_geo(ax_xy, ax_yz, ax_xz):
    n_order, w_well, A_well, rbase = geo_force_params()
    
    patch_spot_xy = py.Circle((0, 0), 0.5*w_well*1e6, fc='#ff8c00',alpha = 0.8)
    # patch_spot_yz = plt.patches.Arc((0, 0), 0.5*w_well*1e6, 0.5*w_well*1e6,0, 0, 180, fc='#ff8c00',alpha = 0.8)
    
    
    py.gcf().sca(ax_xy)
    py.gca().add_patch(patch_spot_xy)
    
    # py.sca(ax_yz)
    # py.gca().add_patch(patch_spot_yz)
    
    
    return 0


    
    
    
def draw_geo(tm, ax_xy, ax_yz, ax_xz):
    # March 7, 2021
    
    # The flag_source_state variable is used to draw/erase the source geometry only once
    # This is necessary to speed up the animation.
    global flag_source_state_1  # Make this variable global so that the assigned value remains saved globally as t changes
    global flag_source_state_2
    global str1
    global str2
 
    if 'flag_source_state_1' not in globals():
        global flag_source_state     # Make this variable global so that the assigned value remains saved globally as t changes
        flag_source_state_1 = 0        # initialize with OFF state
        print('Defining global flag for dynamic geometry \n')
    
    if 'flag_source_state_2' not in globals():
        global flag_source_state       # Make this variable global so that the assigned value remains saved globally as t changes
        flag_source_state_2 = 0        # initialize with OFF state
        print('Defining global flag for static geometry \n')
      
    
    # Draw static geometry (only once)
    if flag_source_state_2 < 1:
        draw_static_geo(ax_xy, ax_yz, ax_xz)
        flag_source_state_2 = 1 
        
    
    # Draw source
    if (tm > 1) & (tm < 8) & (flag_source_state_1 < 1):
        draw_dynamic_geo(ax_xy, ax_yz, ax_xz)
        str1 = 'Optical beam ON'
        str2 = ''
        # text_string2.set_text(str2)  
        # ax_xy.text(0.05, 0.8, 'Optical spot ON',color = '#FF0000',transform=ax_xy.transAxes)
        flag_source_state_1 = 1
        print('Drawing source\n')
        
        
        
        
    # Erase source (draw a white circle)   
    if (tm > 8) & (flag_source_state_1 == 1): 
        n_order, w_well, A_well, rbase = geo_force_params()
        patch_spot = py.Circle((0, 0), 0.51*w_well*1e6, fc='#f9f9f9',alpha = 1)
        py.gca().add_patch(patch_spot)
        str1 = 'Optical beam OFF'
        str2 = ''
        # ax_xy.text(0.05, 0.8, 'Optical spot ON',color = '#f9f9f9', transform=ax_xy.transAxes)
        # ax_xy.text(0.05, 0.8, 'Optical spot OFF',color = '#FF0000',transform=ax_xy.transAxes)
        print('Erasing source\n')
        flag_source_state_1 = 0
    
    if 'str1' not in globals():
        str1 = ''
        str2 = ''
    
    return str1, str2
   
