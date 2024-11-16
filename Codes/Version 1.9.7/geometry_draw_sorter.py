#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 13:04:31 2021

@author: Mohammad Asif Zaman

"""

import numpy as np
import pylab as py

from parameters import *
from geometry_def_sorter import *

Le = 350
Wg1 = 46
Wg2 = 70
Ws = 20
We = 16

fluid_wall_y = 120
fluid_wall_yW = 30
fluid_wall_zW = 100


# Edit on Aug 11, 2022. The new python updata creats issues with the sca function. So, copy-pasting
# the old sca function here. Ref: https://github.com/matplotlib/matplotlib/issues/19380. Use if error occrus.

# def sca(ax):
#     """
#     Set the current Axes instance to *ax*.

#     The current Figure is updated to the parent of *ax*.
#     """
#     managers = _pylab_helpers.Gcf.get_all_fig_managers()
#     for m in managers:
#         if ax in m.canvas.figure.axes:
#             _pylab_helpers.Gcf.set_active(m)
#             m.canvas.figure.sca(ax)
#             return
#     raise ValueError("Axes instance argument was not found in a figure")


# Set particle colors
def particle_color(p):
    if (p % 2) == 0:
         cl_p = cl_ored   # red color for even indexed particles
    else:
         cl_p = cl_pblue  # blue color for odd indexed particles
    
    return cl_p





def draw_static_geo(ax_xy, ax_yz, ax_xz):
     Np, ro, tfinal, x_lim, y_lim, z_lim, xi_lim, yi_lim, zi_lim = time_pos_ax_limits()

     
     
     # Define electrodes
     poly_cnt = np.array([ [-Le/2, Wg1/2], [Le/2, Wg2/2], [Le/2, -Wg2/2], [-Le/2, -Wg1/2] ])
     poly_top = np.array([ [-Le/2, Wg1/2+Ws+We], [Le/2, Wg2/2+Ws+We], [Le/2, Wg2/2+Ws], [-Le/2, Wg1/2+Ws] ])
     poly_bottom = np.array([ [-Le/2, -Wg1/2-Ws-We], [Le/2, -Wg2/2-Ws-We], [Le/2, -Wg2/2-Ws], [-Le/2, -Wg1/2-Ws] ])
     
     center_electrode = py.Polygon(poly_cnt,closed = True, fc = cl_gold, ec = 'k')        
     top_electrode = py.Polygon(poly_top,closed = True, fc = cl_gold, ec = 'k')        
     bottom_electrode = py.Polygon(poly_bottom,closed = True, fc = cl_gold, ec = 'k')        
     
     # Draw some fluid walls (not automated yet)
     # fluid_wall_1_yz = py.Rectangle((-fluid_wall_y, 0),  -fluid_wall_yW, fluid_wall_zW, fc=cl_army_green, ec = 'k')
     # fluid_wall_2_yz = py.Rectangle((fluid_wall_y, 0),  fluid_wall_yW, fluid_wall_zW, fc=cl_army_green, ec = 'k')
     
     
     py.gcf().sca(ax_yz)
     substrate_yz = py.Rectangle((y_lim[0], z_lim[0]),(y_lim[1] - y_lim[0]), abs(z_lim[0]),fc=cl_dgrey, ec='k')
     py.gca().add_patch(substrate_yz)
     # py.gca().add_patch(fluid_wall_1_yz)
     # py.gca().add_patch(fluid_wall_2_yz)
     
    
     py.gcf().sca(ax_xz)
     substrate_xz = py.Rectangle((x_lim[0], z_lim[0]),(x_lim[1] - x_lim[0]), abs(z_lim[0]),fc=cl_dgrey, ec='k')
     py.gca().add_patch(substrate_xz)
    
     py.gcf().sca(ax_xy)
     substrate_xy = py.Rectangle((x_lim[0], y_lim[0]), (x_lim[1] - x_lim[0]), (y_lim[1] - y_lim[0]),fc=cl_lgrey)
     py.gca().add_patch(substrate_xy)
     py.gca().add_patch(center_electrode)
     py.gca().add_patch(top_electrode)
     py.gca().add_patch(bottom_electrode)
     
     
     
     
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
    patch_spot_xy = py.Circle((0, 0), 0.5*w_well*1e6, fc=cl_yellow,alpha = 0.8)
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
        print('Defining global flag for source geometry \n')
    
    if 'flag_source_state_2' not in globals():
        global flag_source_state       # Make this variable global so that the assigned value remains saved globally as t changes
        flag_source_state_2 = 0        # initialize with OFF state
        print('Defining global flag for static geometry \n')
      
    
    # Draw static geometry (only once)
    if flag_source_state_2 < 1:
        draw_static_geo(ax_xy, ax_yz, ax_xz)
        flag_source_state_2 = 1 
        print('Static geometry drawn \n')
        
    

        
        
        

    
    if 'str1' not in globals():
        str1 = '' #'Live cell'
        str2 = '' #'Dead cell'
    
    return str1, str2
   
 
