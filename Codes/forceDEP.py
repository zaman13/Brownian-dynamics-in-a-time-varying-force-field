#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 17:44:46 2021

@author: asif
"""

import numpy as np

xrange_limit = 250e-6





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
    fm[2,:] = 0     
    return fm

# Random force component
def force_random(rr,r_active,Np,t):
    # This force component comes from random origins with random magnitude
    k =  .3e-6    # spring constant
    r_mag = 1/50  # magnitude of the random force component
    fr = np.zeros((3,Np))
    if t> 0.1 and 200*t%1 == 0:
        fr[0,:] = -5*r_mag*k*(rr[0,:]-np.random.normal(0,1,Np)*xrange_limit)*np.random.normal(0,1,Np)
        fr[1,:] = -r_mag*k*(rr[1,:]-np.random.normal(0,1,Np)*xrange_limit)*np.random.normal(0,1,Np)     # random y force 1/20th the magnitude of the x component
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



# Time dependent force field. Implementation of the switching sequence
def force_movingDEP(r,t, elec_spacing,Np):
    
    # define switching instances
   
    ts = np.linspace(0,38,20)  
  
    # define active electrode position at switching
    x_e5 = elec_spacing*2e-6
    x_e4 = elec_spacing*1e-6
    x_e3 = elec_spacing*0
    x_e2 = -elec_spacing*1e-6
    x_e1 = -elec_spacing*2e-6
   
    r_active = x_e5 

    # Assigning active electrode location based on time. 
    # Note the a if statement can be overridden by the subsequent if statement    
    if t < ts[0]:
        r_active = x_e5
    if t >= ts[1]:
        r_active = x_e4
    if t >= ts[2]:
        r_active = x_e3
    if t >= ts[3]:
        r_active = x_e2
    if t >= ts[4]:
        r_active = x_e1
    if t >= ts[5]:
        r_active = x_e2
    if t >= ts[6]:
        r_active = x_e3
    if t >= ts[7]:
        r_active = x_e4
    if t >= ts[8]:
        r_active = x_e5
    if t >= ts[9]:
        r_active = x_e4
    if t >= ts[10]:
        r_active = x_e3
    if t >= ts[11]:
        r_active = x_e2
    if t >= ts[12]:
        r_active = x_e1
    if t >= ts[13]:
        r_active = x_e2
    if t >= ts[14]:
        r_active = x_e3
    if t >= ts[15]:
        r_active = x_e4
    if t >= ts[16]:
        r_active = x_e5
    if t >= ts[17]:
        r_active = x_e4
    if t >= ts[18]:
        r_active = x_e3
    if t >= ts[19]:
        r_active = x_e2

    
    
    
    return force_origin(r,r_active,t,Np)