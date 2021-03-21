#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 6, 2021

@author: asif
"""

import numpy as np
import pylab as py
import matplotlib as plt
from parameters import *

# Parameters:
# =============================================================================

# Particle parameters (number and raidus array)
Np = 4                             # Number of particles
ro =  np.zeros((Np,1)) + 2e-6
ro[0] = 1.5e-6
ro[1] = 2.5e-6

# Time parameter
tfinal = 12


# Axes parameter
xrange_limit = 30e-6    # Max and min of x axis range for plotting animation
zlow_limit = -10e-6
zhigh_limit = 30e-6


# Force parameters
r_active = 0
n_order = 1             # Order of the Gaussian potential = 2n
w_well = 10e-6          # 1/e *max width of the potential well
A_well = 4000*k_B*T   # well depth
# =============================================================================






def draw_geo(tm, ax_xy, ax_yz, ax_xz):
    # March 7, 2021
    
    # The flag_source_state variable is used to draw/erase the source geometry only once
    # This is necessary to speed up the animation.
    global flag_source_state_1  # Make this variable global so that the assigned value remains saved globally as t changes
    global flag_source_state_2
        
    if 'flag_source_state_1' not in globals():
        global flag_source_state     # Make this variable global so that the assigned value remains saved globally as t changes
        flag_source_state_1 = 0        # initialize with OFF state
        print('Defining global flag for source geometry \n')
    
    if 'flag_source_state_2' not in globals():
        global flag_source_state       # Make this variable global so that the assigned value remains saved globally as t changes
        flag_source_state_2 = 0        # initialize with OFF state
        print('Defining global flag for source geometry \n')
      
    
    # Draw static geometry (only once)
    if flag_source_state_2 < 1:
        py.sca(ax_yz)
        substrate_yz = py.Rectangle((-xrange_limit*1e6, zlow_limit*1e6),2*xrange_limit*1e6, abs(zlow_limit)*1e6,fc='#d4d4d4', ec='k')
        py.gca().add_patch(substrate_yz)
        
        py.sca(ax_xz)
        substrate_xz = py.Rectangle((-xrange_limit*1e6, zlow_limit*1e6),2*xrange_limit*1e6, abs(zlow_limit)*1e6,fc='#d4d4d4', ec='k')
        py.gca().add_patch(substrate_xz)
        
        py.sca(ax_xy)
        substrate_xy = py.Rectangle((-xrange_limit*1e6, -xrange_limit*1e6),2*xrange_limit*1e6,2*xrange_limit*1e6,fc='#f9f9f9')
        py.gca().add_patch(substrate_xy)
        flag_source_state_2 = 1
        
        
    
    # Draw source
    if (tm > 1) & (tm < 8) & (flag_source_state_1 < 1):
        patch_spot_xy = py.Circle((0, 0), 0.5*w_well*1e6, fc='#ff8c00',alpha = 0.8)
        # patch_spot_yz = plt.patches.Arc((0, 0), 0.5*w_well*1e6, 0.5*w_well*1e6,0, 0, 180, fc='#ff8c00',alpha = 0.8)
        
        
        py.sca(ax_xy)
        py.gca().add_patch(patch_spot_xy)
        
        # py.sca(ax_yz)
        # py.gca().add_patch(patch_spot_yz)
        
        flag_source_state_1 = 1
        print('Drawing source\n')
        
    # Erase source (draw a white circle)   
    if (tm > 8) & (flag_source_state_1 == 1): 
        patch_spot = py.Circle((0, 0), 0.51*w_well*1e6, fc='#f9f9f9',alpha = 1)
        py.gca().add_patch(patch_spot)
        print('Erasing source\n')
        flag_source_state_1 = 0
   
    

# def draw_yz(tm):
#     substrate_yz = py.Rectangle((-xrange_limit*1e6, zlow_limit*1e6),2*xrange_limit*1e6, abs(zlow_limit)*1e6,fc='#d4d4d4', ec='k')
#     py.gca().add_patch(substrate_yz)
    

# def draw_xz(tm):
#     substrate_xz = py.Rectangle((-xrange_limit*1e6, zlow_limit*1e6),2*xrange_limit*1e6, abs(zlow_limit)*1e6,fc='#d4d4d4', ec='k')
#     py.gca().add_patch(substrate_xz)


# This is function that is called from the main program
# Simplified spring force model
def force_profile(r_in, t):
    
    Np = r_in[0,:].size
    
    fm = np.zeros((3,Np))
    
    

    r_norm = np.linalg.norm(r_in, axis = 0) + 1e-30
    
    g = A_well*np.exp(-(r_norm/w_well)**(2*n_order))
        
    if (t > 1) & (t<8):
        fm[0,:] = -2*n_order*r_in[0,:]/(r_norm**2) * (r_norm/w_well)**(2*n_order) * g
        fm[1,:] = -2*n_order*r_in[1,:]/(r_norm**2) * (r_norm/w_well)**(2*n_order) * g
        fm[2,:] = -2*n_order*r_in[2,:]/(r_norm**2) * (r_norm/w_well)**(2*n_order) * g
    
    
    # fm[:,2] = 0
    # fm[:,3] = 0
    # fm[:,4] = 0
    # fm[:,5] = 0
    # fm[:,6] = 0
     
    
    return fm


def force_plot():
    Np = 1
    rin = np.zeros((3,Np)) 

    r_in = np.tile(np.linspace(-xrange_limit,xrange_limit,200),(3,1))

    F = force_profile(r_in,2)
    
    py.figure()
    py.plot(r_in[0,:]*1e6,F[0,:]*1e12, label = '$F_x$')
    # py.plot(r_in[1,:]*1e6,F[1,:]*1e12,'.', label = '$F_y$')
    # py.plot(r_in[2,:]*1e6,F[2,:]*1e12,'x', label = '$F_z$')
    py.xlabel('$x$ ($\mu$m)')
    py.ylabel('Force (pN)')
    py.legend()
    
# force_plot()

# draw_source(9)

