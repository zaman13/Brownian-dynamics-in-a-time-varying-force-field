#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 17:44:46 2021

@author: Mohammad Asif Zaman

- April 10, 2021
        - Added active electrode text
- April 11, 2021
        - Changed axis limit units to microns
        
"""

import numpy as np
import pylab as py


# Parameters:
# =============================================================================

def time_pos_ax_limits():
    # Particle parameters (number and raidus array)
    Np = 3                             # Number of particles
    # ro =  np.zeros((Np,1)) + 10e-6
    ro =  np.zeros(Np) + 10e-6
    # ro[0] = 5e-6
    # ro[1] = 8e-6
    
    # Time parameters
    tfinal = 38
    
    # Axes parameters (in microns)
    
    x_lim = [-250, 250]
    y_lim = [-250, 250]
    z_lim = [-20,  150]
    
    # Limit of initial particle positions
    xi_lim = [-80e-6, 80e-6]
    yi_lim = [-80e-6, 80e-6]
    zi_lim = [max(ro)*1.5, 80e-6]
    
    return Np, ro, tfinal, x_lim, y_lim, z_lim, xi_lim, yi_lim, zi_lim

# Electrode array geometry parameters (units in um)
elec_width = 15
elec_spacing = 50
# =============================================================================



def draw_static_geo(ax_xy, ax_yz, ax_xz):
    Np, ro, tfinal, x_lim, y_lim, z_lim, xi_lim, yi_lim, zi_lim = time_pos_ax_limits()
    py.gcf().sca(ax_xy)
    substrate_xy = py.Rectangle((x_lim[0], y_lim[0]),(x_lim[1]-x_lim[0]), (y_lim[1]-y_lim[0]),fc='#f9f9f9')
    py.gca().add_patch(substrate_xy)
    for kk in range(-2,3):
        rectangle = py.Rectangle((x_lim[0]/2, -elec_width/2+kk*elec_spacing),x_lim[1],elec_width,fc='#002366')
        py.gca().add_patch(rectangle)
        # ax.add_patch(rectangle)
    
    py.gcf().sca(ax_yz)
    substrate_yz = py.Rectangle((y_lim[0], z_lim[0]),(x_lim[1]-x_lim[0]), abs(z_lim[0]),fc='#d4d4d4', ec='k')
    py.gca().add_patch(substrate_yz)
     
    py.gcf().sca(ax_xz)
    substrate_xz = py.Rectangle((x_lim[0], z_lim[0]),(x_lim[1]-x_lim[0]), abs(z_lim[0]),fc='#d4d4d4', ec='k')
    py.gca().add_patch(substrate_xz)
    
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

        draw_static_geo(ax_xy, ax_yz, ax_xz)
                
        print('Drawing source\n')
        flag_source_state = 1

    
    return str1, str2




# # Read force data from data file
# Mdata = np.genfromtxt('xsweep_ro=10u,zo=25u.csv',delimiter=',',skip_header=5)
# xdata = Mdata[:,0]*1e-6
# Fxdata = Mdata[:,3]*1e-12
# Fydata = 6*Mdata[:,2]*1e-12
# Fzdata = Mdata[:,4]*1e-12
# # Note that the axis has been changed. The axis from the csv data file is different.




# # Interpolation function when using force value from data file    
# def force_interp(ri,ri_active,Np):
#     # Interpolate function for the imported force data(csv file)
#     fi = np.zeros((3,Np))
#     fi[0,:] = np.interp((ri[1,:]-ri_active),xdata,Fxdata)  
#     fi[1,:] = np.interp((ri[1,:]-ri_active),xdata,Fydata)
#     fi[2,:] = np.interp((ri[1,:]-ri_active),xdata,Fzdata)
#     return fi







# Simplified spring force model
def force_model(rm,r_active,Np):
    k =  .3e-6    # spring constant
    r_mag = 0  # magnitude of the random force component
    fm = np.zeros((3,Np))
    
    fm[0,:] = 0
    fm[1,:] = -k*(rm[1,:]-r_active)
    fm[2,:] = -k*(rm[2,:])     
    return fm

# Random force component
def force_random(rr,r_active,Np,t):
    # This force component comes from random origins with random magnitude
    k =  .3e-6    # spring constant
    r_mag = 1/50  # magnitude of the random force component
    fr = np.zeros((3,Np))
    if t> 0.1 and 200*t%1 == 0:
        fr[0,:] = -5*r_mag*k*(rr[0,:]-np.random.normal(0,1,Np)*x_lim[1])*np.random.normal(0,1,Np)
        fr[1,:] = -r_mag*k*(rr[1,:]-np.random.normal(0,1,Np)*x_lim[1])*np.random.normal(0,1,Np)     # random y force 1/20th the magnitude of the x component
        fr[2,:] = 0     # setting the z component of the force to be 0
    else:
        fr[0,:] = 0
        fr[1,:] = 0
        fr[2,:] = 0
    
    
    return fr




# Force profile centered around the origin    
def force_origin(r,r_active,t,Np):
     # return force_interp(r,r_active,Np) + force_random(r,r_active,Np,t)
    return force_model(r,r_active,Np) + force_random(r,r_active,Np,t)




def active_electrode(t):
    
    x_e5 = elec_spacing*2e-6
    x_e4 = elec_spacing*1e-6
    x_e3 = elec_spacing*0
    x_e2 = -elec_spacing*1e-6
    x_e1 = -elec_spacing*2e-6
    
    
    r_active = x_e5
    ts = np.linspace(0,38,20) 
    strn = 'Active electrode = 5'
    # Assigning active electrode location based on time. 
    # Note the a if statement can be overridden by the subsequent if statement    

        
    if t < ts[0]:
        r_active = x_e5
        strn = 'Active electrode = 5'
    if t >= ts[1]:
        r_active = x_e4
        strn = 'Active electrode = 4'
    if t >= ts[2]:
        r_active = x_e3
        strn = 'Active electrode = 3'
    if t >= ts[3]:
        r_active = x_e2
        strn = 'Active electrode = 2'
    if t >= ts[4]:
        r_active = x_e1
        strn = 'Active electrode = 1'
    if t >= ts[5]:
        r_active = x_e2
        strn = 'Active electrode = 2'
    if t >= ts[6]:
        r_active = x_e3
        strn = 'Active electrode = 3'
    if t >= ts[7]:
        r_active = x_e4
        strn = 'Active electrode = 4'
    if t >= ts[8]:
        r_active = x_e5
        strn = 'Active electrode = 5'
    if t >= ts[9]:
        r_active = x_e4
        strn = 'Active electrode = 4'
    if t >= ts[10]:
        r_active = x_e3
        strn = 'Active electrode = 3'
    if t >= ts[11]:
        r_active = x_e2
        strn = 'Active electrode = 2'
    if t >= ts[12]:
        r_active = x_e1
        strn = 'Active electrode = 1'
    if t >= ts[13]:
        r_active = x_e2
        strn = 'Active electrode = 2'
    if t >= ts[14]:
        r_active = x_e3
        strn = 'Active electrode = 3'
    if t >= ts[15]:
        r_active = x_e4
        strn = 'Active electrode = 4'
    if t >= ts[16]:
        r_active = x_e5
        strn = 'Active electrode = 5'
    if t >= ts[17]:
        r_active = x_e4
        strn = 'Active electrode = 4'
    if t >= ts[18]:
        r_active = x_e3
        strn = 'Active electrode = 3'
    if t >= ts[19]:
        r_active = x_e2
        strn = 'Active electrode = 2'
        
    return r_active, strn


# This is function that is called from the main program
# Time dependent force field. Implementation of the switching sequence
def force_profile(r,t):
    
    # define switching instances
    r_active, str1 = active_electrode(t)
    Np = r[0,:].size
    
    return force_origin(r,r_active,t,Np)