#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 16:14:05 2021

@author: asif
"""


import numpy as np


def wall_collision_detect(r_in,ro):
    
    
    ind = []
    ind = np.where(r_in[2,:] <= ro) 
    
    
    return np.squeeze(ind,axis = 0)




def wall_collision_adjust(r_in,v_in,ro,damping_factor):
    
    v_out = np.copy(v_in)
    
    ro = 1.05*ro
    
    ind = wall_collision_detect(r_in,ro)
    print('No of particle-substrate collisions detected = %i \n' % ind.size)  
    
    for m in range(ind.size):
        v_out[2,ind[m]] = np.abs(v_in[2,ind[m]])
        
        
        
    return v_out
    
    
# Np = 40


# ro = .1


# r = np.random.rand(3,Np,2)
# # v = np.zeros((3,Np,2))
# v = np.random.rand(3,Np,2)

# p1 = r[:,:,0]
# v1 = v[:,:,0]

# i1 = wall_collision_detect(p1,ro)
# vo = substrate_collision_adjust(p1,v1,.1,.2)
# print(vo)