#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 16:47:35 2021

@author: Mohammad Asif Zaman

- May 6, 2021
    - Added color library that can be accessed in other modules
- June 2, 2021
    - Defined N_plot_max (the maximum number of particle data to be plotted)

"""


# Parameters (All SI units)
# =============================================================================
k_B = 1.38e-23             # Boltzmann constant
T = 300                    # Temperature (K)
eta = 8.9e-4               # Dynamic viscosity of water (PaS) 


frame_rate = 30                 # Animation frame rate in fps
Nt =frame_rate*40 # 300 #1501   # Number of time steps


rho = 1055                      # density of polystyrene beads in kg/m3
damping_factor = 0              # Fraction of the velocity that is lost during each collision

N_plot_max = 5                  # Maximum number of particle data to plot

# Np = number of particles, defined inside the force_() functions
# ro = (Np,1) raidus vector, defined inside the force_() functions
# tfinal = final simulation time, defined inside the force_() functions


# Colors
# =============================================================================
cl_gold = '#d4af37'            # Golded yellow
cl_army_green = '#5a6b34'      # Dark green
cl_lgrey = '#f9f9f9'           # light grey
cl_dgrey = '#d4d4d4'           # Dark grey
cl_yellow = '#ff8c00'          # yellow
cl_royal_blue = '#002366'      # royal blue
cl_dred = '#C21808'            # dark red
cl_ored = '#fd3b01'            # orange red
cl_bblue = '#add8e6'           # baby blue
cl_pblue = '#78a2cc'           # pastel blue