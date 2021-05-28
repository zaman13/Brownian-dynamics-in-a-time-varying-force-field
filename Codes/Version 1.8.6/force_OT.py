#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 6, 2021

@author: Mohammad Asif Zaman


- April 10, 2021
        - Added optical spot ON-OFF text
- April 11, 2021
        - Changed axis limit units to microns

- May 27, 2021
        - Moved the geometry drawing to a different module (geometry_def_X and geometry_draw_X)
        
"""

import numpy as np
import pylab as py
import matplotlib as plt


from parameters import *
from geometry_def_OT import *

# Module global Parameters:
# =============================================================================
# Force parameters
r_active = 0
n_order = 1             # Order of the Gaussian potential = 2n
w_well = 10e-6          # 1/e *max width of the potential well
A_well = 500*k_B*T      # well depth
rbase = 1e-6



# Setting fluid flow velocity
vx_flow = 0
vy_flow = 0
vz_flow = 0





# =============================================================================


    


# This is function that is called from the main program
# Simplified Gaussian potential-well model
def force_profile(r_in, t):
    
    # Np = r_in[0,:].size
    Np, ro, tfinal, x_lim, y_lim, z_lim, xi_lim, yi_lim, zi_lim = time_pos_ax_limits()
    
    
    
    fm = np.zeros((3,Np))
    
    

    r_norm = np.linalg.norm(r_in, axis = 0) + 1e-30
    
    g = A_well*np.exp(-(r_norm/w_well)**(2*n_order))/rbase**3
    
    
        
    if (t > 1) & (t<8):
        fm[0,:] = -2*n_order*r_in[0,:]/(r_norm**2) * (r_norm/w_well)**(2*n_order) * g *(ro)**3
        fm[1,:] = -2*n_order*r_in[1,:]/(r_norm**2) * (r_norm/w_well)**(2*n_order) * g *(ro)**3
        fm[2,:] = -2*n_order*r_in[2,:]/(r_norm**2) * (r_norm/w_well)**(2*n_order) * g *(ro)**3
        # fm[0,:] = np.gradient(g,r_in[0,:]) 
        # fm[1,:] = np.gradient(g,r_in[1,:]) 
        # fm[2,:] = np.gradient(g,r_in[2,:]) 
    

    
    # fm[:,2] = 0
    # fm[:,3] = 0
    # fm[:,4] = 0
    # fm[:,5] = 0
    # fm[:,6] = 0
     
    
    return fm



# Plotting function, can be called to plot the force profile in 1D. 
def force_plot():
    Np, ro, tfinal, x_lim, y_lim, z_lim, xi_lim, yi_lim, zi_lim = time_pos_ax_limits()
    
    
    
    r_in = np.linspace(x_lim[0]*.7e-6,x_lim[1]*.7e-6,200)
    F = np.zeros(200)
    
    for m in range(r_in.size):
        temp = r_in[m]
        temp2 = np.tile(temp, (3,Np))
        temp3 = force_profile(temp2, 2)
        F[m] = temp3[0,1]
    
    py.figure()
    py.plot(r_in*1e6,F*1e12/(ro[1]*1e6)**3)
    # py.plot(r_in[1,:]*1e6,F[1,:]*1e12,'.', label = '$F_y$')
    # py.plot(r_in[2,:]*1e6,F[2,:]*1e12,'x', label = '$F_z$')
    py.xlabel('Distance, $r$ ($\mu$m)')
    py.ylabel('$F_r/r_o^3$ ,  (pN/$\mu m^3$)')
    # py.legend()

# force_plot()



