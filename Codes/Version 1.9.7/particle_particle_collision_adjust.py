#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 10:06:24 2021

@author: Mohammad Asif Zaman


Comments added on Feb. 28, 2021

Inputs:
    - r_in is a (3,Np) matrix of vectors containing the 3D cartesian coordinates of Np particles
    - ro is a (Np,) array containing the radius of Np particles (no longer assuming all particles have the same radius)
    - mass_o is a (Np,) array containing the mass of Np particles (no longer assuming all particles have the same mass)
    - v_in is a (3,Np) matrix of vectors indication the velocity vectors of Np particles

Function euclidean_dist(r_in):
    - Calculates the Euclidean distance between every pair of particles. 
    - The index of the output matrix correspond to the particle-pairs whose distances are to be calculated

Function collision_detect(r_in,ro):
    - Finds the index of the particles that are closer than the sum of their diameters
    - The output are two array of indices: ind_fix1, indfix2
    - A given element of ind_fix1[] collides with the corresponding element of ind_fix2[]


- March 20, 2021
    - Modifications made to work when ro and mo are arrays
    - Particles with different mass and radius can now be simulated

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
  
    # We return only the upper triangle matrix as that contains distances of all paired points. The lower
    # triangle matrix contains the same information. Selecting one avoids repeatation. 
    return np.triu(np.sqrt(deq))


def collision_detect(r_in,ro):
    deq = euclidean_dist(r_in)
    
    # Make the sum of radius matrix
    Np = r_in.shape[1]
    ro_sum_M = np.tile(ro, (Np,1)) 
    ro_sum_M = np.triu(ro_sum_M + np.transpose(ro_sum_M))
    
    ind_mat = np.where((deq - ro_sum_M <= 0) & (deq != 0))  
    # Some elements are avoided during the search
    # deq = 0 represents self-distance. 
    # deq = 0 also Or elements of the lower traingle matrix which was truncated in euclidean_dist(.)
    
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


def particle_collision_adjust(r_in,v_in,ro, mass_o, damping_factor):
    ro = 1.05*ro
    v_out = np.copy(v_in)
    
    ind1, ind2 = collision_detect(r_in,ro)

    
    # Randomize collision order
    # m_ord = np.arange(ind1.size)

    for m in range(ind1.size):
    # for m in m_ord:
        i1 = int(ind1[m])
        i2 = int(ind2[m])
        
        # The following equation holds when both particles have the same mass
        r_norm = np.linalg.norm(r_in[:,i1] - r_in[:,i2])
        n_vect = r_in[:,i1] - r_in[:,i2]
        
        if r_norm < ro[i1] + ro[i2]:
              
            print('Particle distance less than diamter! \n')
            print('No of particle-particle collisions detected = %i \n' % ind1.size)     
            
            # r_norm = 2*ro
            # n_vect = 2*ro*n_vect/r_norm
        
        Mass_total = mass_o[i1] + mass_o[i2]
        
        v_out[:,i1] = v_in[:,i1] - (2*mass_o[i2]/Mass_total) * n_vect * np.dot(v_in[:,i1] - v_in[:,i2], n_vect) / r_norm**2
        v_out[:,i2] = v_in[:,i2] - (2*mass_o[i1]/Mass_total) * n_vect * np.dot(v_in[:,i2] - v_in[:,i1], n_vect) / r_norm**2
        
        v_in = np.copy(v_out)   # Make sure that the updated velocity is used for the next collision

        # v_out[:,i1] = v_out[:,i1] - n_vect * np.dot(v_out[:,i1] - v_out[:,i2], n_vect) / r_norm**2
        # v_out[:,i1] = v_out[:,i1] + v_in[:,i1] - n_vect * np.dot(v_in[:,i1] - v_in[:,i2], n_vect) / r_norm**2
        
        
    return v_out*(1-damping_factor)
        
        

# Np = 10


# # ro = .1
# ro = np.zeros(Np) + 0.1
# ro[2] = 3


# r = np.random.rand(3,Np,2)
# # v = np.zeros((3,Np,2))
# v = np.random.rand(3,Np,2)

# p1 = r[:,:,0]
# v1 = v[:,:,0]


# # p1 = np.array([[0, 0],[0,1],[0, 0]])
# # v1 = np.array([[.4, 0],[1, -1],[0, 0]])


# ind1,ind2 = collision_detect(p1,ro)

# vo = particle_collision_adjust(p1,v1,ro,0)
# deq = euclidean_dist(p1)
# v[:,:,0] = vo

# print(v1)
# print(vo)

  


