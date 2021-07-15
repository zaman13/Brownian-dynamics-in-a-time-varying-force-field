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
      
    - n_wall_set is a (3,Nwall) array of vectors (one 3 vector per wall)
      Each reflecting wall is defined by a plane with equation ax + by + cz = d
      The n_wall_set vectors denote the coefficients a,b,c of the planes
    
    - d_wall_set is a (1,Nwall) array
      It denotes the right hand term of the plane equation  

- March 19, 2021
        - Adjustments made so that the code works when ro is an array
- May 4, 2021
        - A vector formula for calculating euclidean distance between the walls and particles

        

"""


import numpy as np



def wall_collision_detect(r_in,ro,p_wall,d_wall,lim_wall):
    
    
    ind = []
    euc_dist = np.abs(np.dot(p_wall, r_in) - d_wall)/np.linalg.norm(p_wall)
    # euc_dist is a (Np,) array. Row m represents particle m's distance from the wall
    
    # Check Euclidean distance
    ind1 = np.where( euc_dist - ro <= 0 )  
    # Check if the particle is within the bounds of the walls
    xmn = min(lim_wall[0],lim_wall[1])
    xmx = max(lim_wall[0],lim_wall[1])
    ymn = min(lim_wall[2],lim_wall[3])
    ymx = max(lim_wall[2],lim_wall[3])
    zmn = min(lim_wall[4],lim_wall[5])
    zmx = max(lim_wall[4],lim_wall[5])
    
    
    
    
    ind2 = np.where((xmn-ro <= r_in[0,:]) & (r_in[0,:] <= xmx+ro) & (ymn-ro <= r_in[1,:]) & (r_in[1,:] <= ymx+ro) & (zmn-ro <= r_in[2,:]) & (r_in[2,:] <= zmx+ro)) 
    
    # Index of particles satisfying both conditions
    ind = np.intersect1d(ind1,ind2)
    
    
    
    # return np.squeeze(ind,axis = 0)
    return ind
  



def wall_collision_adjust(r_in,v_in,ro,damping_factor,n_wall_set,d_wall_set, lim_wall_set):
    
    
    # read variable arguments and store them in a list -> numpy array 
    # print(n_wall_set[0]) 
    # print(n_walls)
    # print(n_walls[0,:])
    
    v_out = np.copy(v_in)
    
    ro = 1.05*ro
    
    
    Nwall = n_wall_set.shape[0]    # Detect number of walls
    
    # Loop over all walls
    for wall in range(Nwall):
        
        p_wall = n_wall_set[wall]           # coefficient vector for the wall
        d_wall = d_wall_set[wall]*1e-6      # offset vector for the wall
        lim_wall = lim_wall_set[wall]*1e-6  # 6 x 1 array containing limits of the walls [xmin, xmax, ymin, ymax, zmin, zmax] format.
        
        n_wall = p_wall/np.linalg.norm(p_wall)
        
        ind = wall_collision_detect(r_in,ro,p_wall,d_wall,lim_wall)   
        print('No of particle-substrate collisions detected = %i \n' % ind.size) 
        
        
        
        for m in range(ind.size):
            # v_out[2,ind[m]] = np.abs(v_in[2,ind[m]])
            # v_out[2,ind[m]] = -v_in[2,ind[m]]
            v_out[:,ind[m]] = v_in[:,ind[m]] - 2*np.dot(v_in[:,ind[m]], n_wall)*n_wall
            if wall == 0:    # wall == 0 indicates the z=0 wall. This is treated specially to avoid any error.
                v_out[2,ind[m]] = np.abs(v_out[2,ind[m]])   # force v_z to be positive when it hits the bottom surface. This force-avoids unphysical cases
            
        
    return v_out



   
    
# Np = 40
# # ro = 0.1
# Np = 3
# ro = np.zeros((Np)) + 0.1
# fluid_wall_y = 120


# r = np.random.rand(3,Np,2)
# # v = np.zeros((3,Np,2))
# v = np.random.rand(3,Np,2)

# p1 = r[:,:,0]
# # v1 = v[:,:,0]
# n_wall_set = [[0,0,1], [0,1,0], [0,1,0] ]                          # plane coefficient vector array
# d_wall_set = [0, fluid_wall_y*1e-6, -fluid_wall_y*1e-6]            # plane right-hand term array


# # Limits on the extent of the walls. Format: [xmin, xmax, ymin, ymax, zmin, zmax] per wall
# lim_wall_set = [[-np.inf, np.inf, -np.inf, np.inf, 0, np.inf],     
#                 [-np.inf, np.inf, -np.inf, np.inf, 0, np.inf],
#                 [-np.inf, np.inf, -np.inf, np.inf, 0, np.inf],                
#                 ]


# n_wall_set = np.array(n_wall_set)                                  # convert to numpy array
# d_wall_set = np.array(d_wall_set)

# lim_wall_set = np.array(lim_wall_set)
# i1 = wall_collision_detect(p1,ro,n_wall_set[0],d_wall_set[0],lim_wall_set[0])
# # print(i1)
# vo = wall_collision_adjust(p1,v1,.1,.2)
# print(vo)