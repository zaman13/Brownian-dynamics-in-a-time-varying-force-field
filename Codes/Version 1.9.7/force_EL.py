#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 17:44:46 2021

@author: Mohammad Asif Zaman

- April 10, 2021
        - Added active electrode text
- April 11, 2021
        - Changed axis limit units to microns
- May 28, 2021
        - Functionalized fluid velocity      
"""

import numpy as np
import pylab as py
from scipy import interpolate

from parameters import *


# Setting fluid flow velocity
def fluid_vel(r_in, t):
    
    Np = r_in[0,:].size
    # xin = r_in[0,:]
    # yin = r_in[1,:]
    # zin = r_in[2,:]
    
    v_fluid = np.zeros((3,Np))
    v_fluid[1,:] = 15e-6
    
    return v_fluid




# Parameters:
# =============================================================================

# Electrode array geometry parameters (units in um)
def geo_force_params():
    elec_width = 15
    elec_spacing = 50
    circ_rad = 15
    circ_offset = 30
    
    return elec_width, elec_spacing, circ_rad, circ_offset
# =============================================================================






# Interpolation function when using force value from data file    
def force_interp(ri,ri_active,Np):
    # Interpolate function for the imported force data(csv file)
    
    # Read force data from data file
    
    xdata = np.genfromtxt('Data_EL/x1.csv',delimiter=',')*1e-6
    ydata = np.genfromtxt('Data_EL/y1.csv',delimiter=',')*1e-6
    
    Fxdata = np.genfromtxt('Data_EL/Fsx.csv',delimiter=',')*1e-12
    Fydata = np.genfromtxt('Data_EL/Fsy.csv',delimiter=',')*1e-12
    Fzdata = np.genfromtxt('Data_EL/Fsz.csv',delimiter=',')*1e-12
    # Note that the axis has been changed. The axis from the csv data file is different.
    
    
    fx_fun = interpolate.interp2d(xdata, ydata, Fxdata, kind='cubic')
    fy_fun = interpolate.interp2d(xdata, ydata, Fydata, kind='cubic')
    fz_fun = interpolate.interp2d(xdata, ydata, Fzdata, kind='cubic')
    
    fi = np.zeros((3,Np))
    
    
    
    if Np > 1:
        fi[0,:] = fx_fun(ri[0,:]-ri_active,ri[1,:]).diagonal()
        fi[1,:] = fy_fun(ri[0,:]-ri_active,ri[1,:]).diagonal()
        fi[2,:] = fz_fun(ri[0,:]-ri_active,ri[1,:]).diagonal()
    else:
        fi[0,:] = fx_fun(ri[0,:]-ri_active,ri[1,:]).squeeze()
        fi[1,:] = fy_fun(ri[0,:]-ri_active,ri[1,:]).squeeze()
        fi[2,:] = fz_fun(ri[0,:]-ri_active,ri[1,:]).squeeze()
        
    
    fi[0,:] = 1.3*fi[0,:]
    fi[1,:] = fi[1,:]
    fi[2,:] = 0.1*fi[2,:]
    # Fx_interpolator = interp.CloughTocher2DInterpolator(np.array([xdata,ydata]).T, Fxdata)
    # Fy_interpolator = interp.CloughTocher2DInterpolator(np.array([xdata,ydata]).T, Fydata)
    # Fz_interpolator = interp.CloughTocher2DInterpolator(np.array([xdata,ydata]).T, Fzdata)
    
    # fi[0,:] =  Fx_interpolator((ri[0,:]-ri_active),ydata,Fxdata)  
    # fi[1,:] =  Fx_interpolator((ri[1,:]-ri_active),ydata,Fydata)  
    # fi[2,:] =  Fx_interpolator((ri[1,:]-ri_active),ydata,Fzdata)  
    return fi





# Simplified spring force model
def force_model(rm,r_active,Np):
    

    

    k =  0 #.3e-6    # spring constant
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
    return force_interp(r,r_active,Np) 
    # return force_model(r,r_active,Np) + force_random(r,r_active,Np,t)




def active_electrode(t):
    
    elec_width, elec_spacing, circ_rad, circ_offset = geo_force_params()
    x_e6 = elec_spacing*3e-6
    x_e5 = elec_spacing*2e-6
    x_e4 = elec_spacing*1e-6
    x_e3 = elec_spacing*0
    x_e2 = -elec_spacing*1e-6
    x_e1 = -elec_spacing*2e-6
    x_e0 = -elec_spacing*3e-6
    
    r_active = x_e6
    ts = np.linspace(0,25,20) 
    strn = 'Active electrode = 5'
    # Assigning active electrode location based on time. 
    # Note the a if statement can be overridden by the subsequent if statement    

    
    t = t-ts[12]*int(t/ts[12])              # Switching will repeat with period = ts[12]
        
        
    if t < ts[0]:
        r_active = x_e6
        strn = 'Active electrode = 6'
    if t >= ts[1]:
        r_active = x_e5
        strn = 'Active electrode = 5'
    if t >= ts[2]:
        r_active = x_e4
        strn = 'Active electrode = 4'
    if t >= ts[3]:
        r_active = x_e3
        strn = 'Active electrode = 3'
    if t >= ts[4]:
        r_active = x_e2
        strn = 'Active electrode = 2'
    if t >= ts[5]:
        r_active = x_e1
        strn = 'Active electrode = 1'
    if t >= ts[6]:
        r_active = x_e0
        strn = 'Active electrode = 0'
    if t >= ts[7]:
        r_active = x_e1
        strn = 'Active electrode = 1'
    if t >= ts[8]:
        r_active = x_e2
        strn = 'Active electrode = 2'
    if t >= ts[9]:
        r_active = x_e3
        strn = 'Active electrode = 3'
    if t >= ts[10]:
        r_active = x_e4
        strn = 'Active electrode = 4'
    if t >= ts[11]:
        r_active = x_e5
        strn = 'Active electrode = 5'
    # if t >= ts[12]:
    #     r_active = x_e6
    #     strn = 'Active electrode = 6'
    # if t >= ts[13]:
    #     r_active = x_e5
    #     strn = 'Active electrode = 5'
    # if t >= ts[14]:
    #     r_active = x_e4
    #     strn = 'Active electrode = 4'
    # if t >= ts[15]:
    #     r_active = x_e3
    #     strn = 'Active electrode = 3'
    # if t >= ts[16]:
    #     r_active = x_e2
    #     strn = 'Active electrode = 2'
    # if t >= ts[17]:
    #     r_active = x_e1
    #     strn = 'Active electrode = 1'
    # if t >= ts[18]:
    #     r_active = x_e0
    #     strn = 'Active electrode = 0'
    # if t >= ts[19]:
    #     r_active = x_e1
    #     strn = 'Active electrode = 1'
        
    return r_active, strn


# This is function that is called from the main program
# Time dependent force field. Implementation of the switching sequence
def force_profile(r,t):
    
    # define switching instances
    r_active, str1 = active_electrode(t)
    Np = r[0,:].size
    
    return force_origin(r,r_active,t,Np)


def force_plotter(xlow,xhigh,ylow,yhigh):
    pNx = 20
    pNy = 20
    
    xi = np.linspace(xlow,xhigh,pNx)
    yi = np.linspace(ylow,yhigh,pNy)
    fx = np.zeros((pNy,pNx))
    fy = np.zeros((pNy,pNx))
    fz = np.zeros((pNy,pNx))
    
  
    
    ri = np.zeros((3,1))
  
    
    
    tx = np.zeros((pNx,pNy))
    
    for m in range(pNx):
        for n in range(pNy):
            ri[0,:] = xi[m]
            ri[1,:] = yi[n]
            
            tx[n, m] = force_interp(ri,60e-6,1)[0]
            # tx[n, m] = fi[0,:]n(xi[m],yi[n])
            
    
    py.contourf(xi,yi,tx)
    py.colorbar()
    
    
    return 0
            