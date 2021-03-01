#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 10:06:24 2021

@author: asif


Comments added on Feb. 28, 2021

Inputs:
    - r_in is a (3,Np) matrix of vectors containing the 3D cartesian coordinates of Np particles
    - ro is the radius of the particles (assuming all particles have the same size for now)
    - v_in is a (3,Np) matrix of vectors indication the velocity vectors of Np particles

Function euclidean_dist(r_in):
    - Calculates the Euclidean distance between every pair of particles. 
    - The index of the output matrix correspond to the particle-pairs whose distances are to be calculated

Function collision_detect(r_in,ro):
    - Finds the index of the particles that are closer than the sum of their diameters
    - The output are two array of indices: ind_fix1, indfix2
    - A given element of ind_fix1[] collides with the corresponding element of ind_fix2[]
    

"""

import numpy as np


def euclidean_dist(r_in):
    
    Np = r_in.shape[1]
    
    
    # Decomposing the input matrix in to Cartesian components 
    xp1 = np.tile(r_in[0,:],(Np,1))
    yp1 = np.tile(r_in[1,:],(Np,1))
    zp1 = np.tile(r_in[2,:],(Np,1))   # Ignoring the z-component to limit motion in xy plane only
    
    # Subtracting the matrix with its transpose
    deq = np.square(xp1-np.transpose(xp1)) + np.square(yp1-np.transpose(yp1)) + np.square(zp1-np.transpose(zp1))    
  
    
    return np.sqrt(deq)


def collision_detect(r_in,ro):
    deq = euclidean_dist(r_in)
    ind_mat = np.where(deq < 2.0*ro)
    ind1 = ind_mat[0]
    ind2 = ind_mat[1]
    
    # Defining empty numpy arrays
    ind_fix1 = np.array([])
    ind_fix2 = np.array([])
    
    
    # Eliminating repeated entries
    for m in range(ind1.size):
        if ind1[m] != ind2[m]:
            ind_fix1 = np.append(ind_fix1, ind1[m])
            ind_fix2 = np.append(ind_fix2, ind2[m])
    
    

    # print('No of collision detected = %i \n' % ind_fix1.size)        
    return ind_fix1, ind_fix2


def velocity_adjust(r_in,v_in,ro):
    ro = 1.05*ro
    v_out = np.copy(v_in)
    
    ind1, ind2 = collision_detect(r_in,ro)

    # # Reset velocities for the particles that are colliding
    # for m in range(ind1.size):
    #     i1 = int(ind1[m])
    #     v_out[:,i1] = 0
    

    for m in range(ind1.size):
        i1 = int(ind1[m])
        i2 = int(ind2[m])
        
        # The following equation holds when both particles have the same mass
        r_norm = np.linalg.norm(r_in[:,i1] - r_in[:,i2])
        n_vect = r_in[:,i1] - r_in[:,i2]
        
        if r_norm < 2*ro:
              
            print('Particle distance less than diamter! \n')
            print('No of collision detected = %i \n' % ind1.size)     
            r_norm = 2*ro
            n_vect = 2*ro*n_vect/r_norm
            
        v_out[:,i1] = v_in[:,i1] - n_vect * np.dot(v_in[:,i1] - v_in[:,i2], n_vect) / r_norm**2
        # v_out[:,i1] = v_out[:,i1] - n_vect * np.dot(v_out[:,i1] - v_out[:,i2], n_vect) / r_norm**2
        # v_out[:,i1] = v_out[:,i1] + v_in[:,i1] - n_vect * np.dot(v_in[:,i1] - v_in[:,i2], n_vect) / r_norm**2
        
        
    return v_out
        
        
    

  


