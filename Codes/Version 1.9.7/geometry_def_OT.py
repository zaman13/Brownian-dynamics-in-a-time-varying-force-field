#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  8 13:18:12 2021

@author: Mohammad Asif Zaman

Date: May 15, 2021

This code defines geoemtrical elements by a set of points. It also calculates the 
normal vector for each surface of that geometry. 

Definition parameters:
    - In geo_pointsX() functions, the point array set p_wall[] contains the edge point of each geometry element
    - All numbers in the geo_pointsX() functions are in microns (um)

Outputs:
    - n_vector_set: Normal vectors for all the walls n = (a,b,c), n_vector_set = [(a1,b1,c1), (a2,b2,c2)....]. Size = (Nwalls,3)
    - d_set       : Set of d coefficients defining the wall planes ax + by + cz = d. Size (Nwalls,)
    - p_lim_set   : x,y,z limist for each wall. p_lim_1 = [xlow_1, xhigh_1, ylow_1, yhigh_1, zlow_1, zhigh_1], p_lim_set = [p_lim_1,p_lim_2....]. Size (Nwalls,6)   
    
    - Here, Nwalls = number of walls

- May 17, 2021: 
    - Added different geometry segment types (rectangular, triangular etc.)
    - Implemented call function by string 
    - Functionalized different parts



"""

import numpy as np



# This function states which elements are present in the definition of the geometry
def geo_element_types():
    # Types of geometrical elements used
    # 4 = rectangular geometry segments (geo_points4() function must exist)
    # 3 = triangular geometry segments  (geo_points3() function must exist)
    # .... and so on
    
    return np.array([])   # return empty array (no walls defined). Only the default z = 0 reflecting wall will be considered.
    

   


# # Triangular geometry segments
# def geo_points3():
        
#     # Points defining the geometry. Format 
#     #  [ [ [x1,y1], [x2,y2], [x3,y3], [x4,y4]  ]     # first geometry segment 
#     #    [ [X1,Y1], [X2,Y2], [X3,Y3], [X4,Y4]  ]     # second geometry segment 
#     #    [                                     ]     # .... (so on)
#     #  ]
#     p_wall   = [        [  [2,0], [8,0], [8,8]         ],
#                         [  [-20,-20], [-10,-20], [-10,-26]      ],
                        
#                       ]
    
    
#     # z range of each of the points
#     z_wall   = [ [0,140],
#                   [0,140],
                 
#                 ]
    
#     p_wall   = np.array(p_wall)

#     return p_wall, z_wall
    





# Calculates normal vector,d coeff and limits of x,y,z for a specific geometry set (example: geo_points4() )
def normal_v(pset,zset):
    n_vector = [[]]       # Set of normal vectors for each reflecting surface. Each element: n = (a,b,c)
    d_coeff = []          # Corresponding d coefficient for ax + by + cz = d plane equation
    p_lims =[[]]          # x,y,z extent/limit of the wall. Format: [xlow, xhigh, ylow, yhigh, zlow, zhigh]   
    Nw = pset.shape[0]    # Number of walls
    
    for q in range(Nw):   # goes through every wall geometry
        
        ps = pset[q]      # Select one geometry segment. ps contains all points defining that geometry segment. 
        zs = zset[q]
        M = ps.shape[0]   # Number of points in that geometry segment.For rectangle, it's 4 points.
         
        i2 = 0
        
        for m in range(M):  # goes through all the points in qth wall geometry
            
            # indices for pair of points
            i1 = m        
            i2 = m + 1 if i2 < M-1 else 0   # for the last point, select i2 = 0 as the index of the second point
            
            # making 3d coordiantes from 2d using z = 0 value for p1, p2 and z = 1 value for p3
            # 3 points are defined
            p1 = np.append(ps[i1],zs[0])
            p2 = np.append(ps[i2],zs[1])
            p3 = np.append(ps[i2],zs[0])    # a third arbitrary point at a different z value
            
            # min-max limit of the three points
            # format: xmin, xmax, ymin, ymax, zmin, zmax
            temp = [min(p1[0],p2[0],p3[0]),max(p1[0],p2[0],p3[0]),
                    min(p1[1],p2[1],p3[1]),max(p1[1],p2[1],p3[1]),
                    min(p1[2],p2[2],p3[2]),max(p1[2],p2[2],p3[2])]
            
            
            # Note: We are assuming that the geoemtry segments are uniform in the z direction (they extrude up from the xy plane in the z direction)
            # For generalized case, the z points should be defined within the definition of the points. However, for most LOC microfluidic devices, 
            # there are no variations along the z direction and this assumption holds. 
            
            nv = np.cross(p3-p2,p2-p1)  # normal vector to the plane containing the three points
            d = np.dot(nv,p1)           # corresponding d coefficient
            
            d = d/np.linalg.norm(nv)
            nv = nv/np.linalg.norm(nv)
            
           
            
            n_vector = np.append(n_vector,[nv], axis = 1) if m+q == 0 else np.append(n_vector,[nv], axis = 0)  
            p_lims = np.append(p_lims,[temp], axis = 1) if m+q == 0 else np.append(p_lims,[temp], axis = 0)  
            # Note for n_vector if statement: When appending the first element within the empty
            # n_vector, the axis is 1. For other cases, axis = 0. This ensures that the data structure is 
            # correct and the array of arrays don't blend into a 1D array. The first element occurs at 
            # q = 0 and m = 0 (thus, when m+q = 0). 
            
            d_coeff = np.append(d_coeff,d)
            # print(n_vector)
    
    return n_vector,d_coeff, p_lims
    # return ps


# Calculate normal vectors, d coeff and (x,y,z) limits for all geometries (combines results for diff type of geometries)



# Main function that will return normal vectors for all the geometry segments, d coefficients and (x,y,z) limits
def geo_define():
    
    el_type = geo_element_types()       # Array containing types of geometries used
    n_vector_set = np.array([[0,0,1]])            # List for storing normal vectors with z = 0 plane inputted as the 0th wall (substrate)
    d_set = np.array([0])                         # List for storing d coefficient with z = 0 plane inputted as the 0th wall 
    p_lim_set =np.array([[-np.inf, np.inf, -np.inf, np.inf, -np.inf, np.inf]])   # List for storing the limits/extent of each wall with an infinite z = 0 substrate already inputted
    
    # # Use these empty list initialization if z = 0 substrate is not present
    # n_vector_set = [[]]            # List for storing normal vectors 
    # d_set = []                         # List for storing d coefficient 
    # p_lim_set =[[]]   # List for storing the limits/extent of each wall 
    
    
    
    for m in range(el_type.shape[0]):             # looping over all geometry element types
        
        fname = 'geo_points' + str(el_type[m]) + '()'  # name of the corresponding geometry function
        
        # Evaluate geo_points*()function by string name
        p_w1,z_w1 = eval(fname)
        
        
        nv,d,pp = normal_v(p_w1, z_w1)
        
        
        n_vector_set = np.append(n_vector_set,nv, axis = 0)  
        d_set = np.append(d_set,d)  
        p_lim_set = np.append(p_lim_set,pp, axis = 0)  
        
        
        # # Use the following if you wish to initialize with an empty list instead of z = 0 default plane
        # n_vector_set = nv if m == 0 else np.append(n_vector_set,nv, axis = 0)  
        # d_set = d if m == 0 else np.append(d_set,d)  
        # p_lim_set = pp if m == 0 else np.append(p_lim_set,pp, axis = 0)  
        
    return n_vector_set, d_set, p_lim_set
    


# Particle parameters (number and raidus array)
def time_pos_ax_limits():
    Np = 4                             # Number of particles
    # ro =  np.zeros((Np,1)) + 2e-6
    ro =  np.zeros(Np) + 2e-6
    ro[0] = 1.5e-6
    ro[1] = 2.5e-6
    
    # Time parameter
    tfinal = 30  #12
    
    
    # Axes parameter (in microns)
    
    x_lim = [-30, 30]
    y_lim = [-30, 30]
    # z_lim = [-5,  25]
    z_lim = [-5,  55]
    
    # Limit of initial particle positions (in meters)
    xi_lim = [-10e-6, 10e-6]
    yi_lim = [-10e-6, 10e-6]
    zi_lim = [max(ro)*1.5, 20e-6]
    
    return Np, ro, tfinal, x_lim, y_lim, z_lim, xi_lim, yi_lim, zi_lim


# n_vector_set, d_set, p_lim_set = geo_define()
    
# print(n_vector_set)
# print(d_set)
# print(p_lim_set)





