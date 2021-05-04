#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 16:14:05 2021

@author: Mohammad Asif Zaman

Input arguments:
    - r_in is a (3,Np) array
    - ro is a (Np,) array. For (Np,1) array there are some issues. 
      So, a squeeze command maybe required when calling the function with (Np,1) vector.
      For version 1.7.2 and up, this should not be any issue as ro and mo vectors are redefined as (Np,)
      vectors instead of (Np,1) vectors

- March 19, 2021
        - Adjustments made so that the code works when ro is an array



"""


import numpy as np


def wall_collision_detect(r_in,ro):
    
    
    ind = []
    ind = np.where(r_in[2,:] - ro <= 0) 
        
    return np.squeeze(ind,axis = 0)
  



def wall_collision_adjust(r_in,v_in,ro,damping_factor,*argv):
    
    
    # read variable arguments and store them in a list -> numpy array 
    n_sidewalls = []
    count = 0
    for arg in argv:
        n_sidewalls.append(arg)
    
    n_sidewalls = np.array(n_sidewalls)
    # print(n_sidewalls)
    # print(n_sidewalls[0,:])
    
    v_out = np.copy(v_in)
    
    ro = 1.05*ro
    
    ind = wall_collision_detect(r_in,ro)
    print('No of particle-substrate collisions detected = %i \n' % ind.size)  
    
    for m in range(ind.size):
        v_out[2,ind[m]] = np.abs(v_in[2,ind[m]])
        
        
        
    return v_out
    
    
# Np = 40
# # ro = 0.1

# ro = np.zeros((Np)) + 0.1


# r = np.random.rand(3,Np,2)
# # v = np.zeros((3,Np,2))
# v = np.random.rand(3,Np,2)

# p1 = r[:,:,0]
# v1 = v[:,:,0]

# i1 = wall_collision_detect(p1,ro)
# # print(i1)
# vo = wall_collision_adjust(p1,v1,.1,.2)
# print(vo)