#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 11:31:12 2021

@author: asif
"""

import numpy as np

# Simplified spring force model
def force_trap(rm,r_active, t,Np):
    k =  .3e-6    # spring constant
    r_mag = 0  # magnitude of the random force component
    fm = np.zeros((3,Np))
    
    
    
    if (t > 1) :
        fm[0,:] = -k*(rm[0,:]-r_active)
        fm[1,:] = -k*(rm[1,:]-r_active)
        fm[2,:] =  -k*(rm[2,:]-r_active)     
    else:
        fm[0,:] = 0
        fm[1,:] = 0
        fm[2,:] = 0
    
    
    
    
    # fm[:,2] = 0
    # fm[:,3] = 0
    # fm[:,4] = 0
    # fm[:,5] = 0
    # fm[:,6] = 0
     
    
    return fm