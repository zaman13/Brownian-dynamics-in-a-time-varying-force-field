#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 11:31:12 2021

@author: asif
"""

import numpy as np

xrange_limit = 100e-6    # Max and min of x axis range for plotting animation
r_active = 0


def draw_source():
    return 0

# This is function that is called from the main program
# Simplified spring force model
def force_profile(r_in, t):
    
    Np = r_in[0,:].size
    k =  .3e-6    # spring constant
    
    fm = np.zeros((3,Np))
    
    
    
    if (t > 1) :
        fm[0,:] = -k*(r_in[0,:]-r_active)
        fm[1,:] = -k*(r_in[1,:]-r_active)
        fm[2,:] =  -k*(r_in[2,:]-r_active)     

    
    
    
    # fm[:,2] = 0
    # fm[:,3] = 0
    # fm[:,4] = 0
    # fm[:,5] = 0
    # fm[:,6] = 0
     
    
    return fm