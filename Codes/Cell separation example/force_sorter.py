#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 6, 2021

@author: Mohammad Asif Zaman


- May 28, 2021
        - Functionalized fluid velocity        
"""

import numpy as np
import pylab as py
import matplotlib as plt
import time

from scipy import interpolate

from parameters import *
from geometry_def_sorter import *

# Module global Parameters:
# =============================================================================


# Setting fluid flow velocity
# vx_flow = 120e-6
# vy_flow = 0e-6
# vz_flow = 0e-6



def vel_lm(xin, yin, zin):   # linear model of fluid velocity
    
    Np = xin.size
    
    v_out = np.zeros((3,Np))
    
    v_out[0,:] = 150e-6
    v_out[1,:] = 0
    v_out[2,:] = 0
    
    
    v_out[1,np.where( (xin > 170e-6) & (yin >= 30e-6 ) )]  = 30e-6
    v_out[1,np.where( (xin > 170e-6) & (yin <= -30e-6 ) )] = -30e-6
      
    return v_out



# def fluid_vel(r_in, t):
    
#     Np = r_in[0,:].size
#     xin = r_in[0,:]
#     yin = r_in[1,:]
#     zin = r_in[2,:]
    
#     v_fluid = np.zeros((3,Np))
    
#     v_fluid[0,:] = 120e-6
#     v_fluid[1,:] = 0
#     v_fluid[2,:] = 0
    
    
#     v_fluid[1,np.where( (xin > 170e-6) & (yin >= 30e-6 ) )]  = 30e-6
#     v_fluid[1,np.where( (xin > 170e-6) & (yin <= -30e-6 ) )] = -30e-6
    
    
    
    
#     return v_fluid

def fluid_vel(r_in, t):
    
    Np = r_in[0,:].size
    xi = r_in[0,:]
    yi = r_in[1,:]
    zi = r_in[2,:]
    d = 10e-6
    
    # temporary variables
    v1 = np.zeros((3,Np))
    v2 = np.zeros((3,Np))
    
    # Moving average smoothing of fluid velocity 
    
    # number of points per axis over which to average the velocity predicted by the linear velocity model
    N_avg_points = 15  # must be an odd number
    for m in range(int((N_avg_points-1)/2 + 1)):
        v1 = v1 + vel_lm(xi+d*m,yi,zi) + vel_lm(xi-d*m,yi,zi) if m > 0 else v1 + vel_lm(xi,yi,zi)
        v2 = v2 + vel_lm(xi,yi+d*m,zi) + vel_lm(xi,yi-d*m,zi) if m > 0 else v2 + vel_lm(xi,yi,zi)


    v_fluid = (v1 + v2)/(2*N_avg_points)
    
     
    return v_fluid
  
    
  
    
 # Read force data from data file
Mdata = np.genfromtxt('Fy_XY_grid2.csv',delimiter=',',skip_header=9)
xdata = Mdata[:,0]*1e-6
ydata = Mdata[:,1]*1e-6
points = np.array( (xdata, ydata)).T
Fydata = Mdata[:,2]*1e-12







# This is function that is called from the main program
# Simplified spring force model
def force_profile(r_in, t):
    
    # Np = r_in[0,:].size
    Np, ro, tfinal, x_lim, y_lim, z_lim, xi_lim, yi_lim, zi_lim = time_pos_ax_limits()
    
    temp = [.7,-0.25]
    od_ev = np.array(int(Np/2)*[temp]).flatten()
    
    
    xin = r_in[0,:]
    yin = r_in[1,:]
    
    fy = interpolate.griddata(points,Fydata,(xin,abs(yin)),method='linear',fill_value = 0)*np.sign(yin)
    fz = -od_ev*.3e-12
    fm = np.zeros((3,Np))
    fm[1,:] = fy*od_ev
    fm[2,:] = fz
    
    return fm



# force_plot()


# Np = 1
# # # xin = [1, 4, 2, 3]
# # # xin = np.array(xin)
# # # v_temp = np.zeros(Np)
# # # v_temp[np.where(xin > 2)]  = 7 

# r = np.random.rand(3,Np,2)
# rin = r[:,:,0]
# t = 0

# rin[0] = 176e-6
# rin[1] = 50e-6

# vf = fluid_vel(rin,t)
# # print(rin)
# print(vf)

# # vf = fluid_vel([[170e-6,40e-6,2e-6]],t)
# # print(vf)





# # interpolation speed test
# start_time = time.time()
# print('\n\n===========================================\n')

# Np = 24
# r = np.random.rand(3,Np,2)
# rin = r[:,:,0]
# t = 0

# tt= force_profile(rin,t)


# print("Execution time = %1.2f seconds \n" % (time.time() - start_time))
# print('\n===========================================\n')

