#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 16:47:35 2021

@author: Mohammad Asif Zaman
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


# Np = number of particles, defined inside the force_() functions
# ro = (Np,1) raidus vector, defined inside the force_() functions
# tfinal = final simulation time, defined inside the force_() functions