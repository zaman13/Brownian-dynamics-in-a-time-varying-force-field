# -*- coding: utf-8 -*-
"""


@author: Mohammad Asif Zaman
Date: Dec. 2020

Modeling Brownian motion of a colloidal particle



Coordinate system notes: 
    - r is the 3D position vector (Cartesian). 
    - r[0] = x component 
    - r[1] = y component
    - r[2] = z component


Updates
- Dec. 21, 2020
        - Added force_random() function
- Dec. 22, 2020
        - Using a full force-field data file (spanning both positive and negative x)
        - Changed the axis. Electrodes are now horizontal. Motion is vertical. Note that the force data file 
          recorded assuming a differet coordinate (vertical electrodes and xsweep). So, the x and the y axis 
          has been reversed.
- Feb 18, 2021
        - Added xlabel and ylabel on the animation 
        - Added more comments
        - Uploaded as v1.0
- Feb 21. 2021
        - Working on multiparticle system. Working for Np = 2 case
        - Animation for 2 particles added
        - Generalization required
- Feb 22, 2021
        - Array index order changed to [space dimension, Np, Nt]
        - Animation working for arbitrary Np values
        - The force_interp() function needs to be rewritten. force_model() and force_rand() has been modified.
        - Uploaded as v1.4
- Feb 23, 2021
        - Moved the force related functions to a separate file (forceDEP.py)
        - Interpolation function for the experimental force data has been modified accordingly
        - Executation/elapsed time added
        - Status messages added
        - uploaded as v1.5
- Feb 24-28, 2021
        - Added particle-particle elastic collision dynaimics
        - Works fine for two simultaneous collisions. For larger number of simultaneous collisions, some modifications are needed.                - 

- March 2, 2021
        - Decoupled the simulation time step to animation frame rate. They can be independently defined now.
- March 3, 2021
        - Moved the draw_electrode() function outside the init() function
- March 3, 2021
        - Removed excess parameters from the force functions
        - Made the different force functions, i.e. forceDEP.py and force_spring_trap.py more consistent with each other
          Only the import call and final time needs to be adjusted for when switching the force functions. The main output
          from both imports are function_profile(r,t).
- March 7, 2021
        - Added zorder for animating the beads (beads append). This ensures that the beads are always in the foreground (compared to source geometry)        
        - Added draw_source(t) function inside the animate() function. Also, made the draw_source() as a funciton of time for dynamic
          manipulation of the source geometry (e.g. turning ON/OFF optical excitation for optical trapping demo.)
- March 8-10, 2021
        - Added x,y,z position vs time plots
        - Added animation for yz, and zx planes along with xy plane
- March 11, 2021
        - Added wall collision mechanics
        - Added substrate graphics
        - Renamed some of the modules
- March 14, 2021
        - Streamlined the plots
        - Fixed random distribution of initial position
- March 17, 2021
        - Steamlined the draw_geo() function within the force files
        - Moved the definition of ro and tfinal from the main file to the force files
- March 19-20, 2021
        - Generalized the code for working with non-homogeneous particles (particles with different radius and different mass)     
        - Released as version 1.7.0
        
"""


from __future__ import print_function    

import time
import math
import numpy as np
import pylab as py

import matplotlib as plt
from matplotlib import animation, rc




from parameters import *
from particle_particle_collision_adjust import *
from particle_wall_collision_adjust import *

from forceDEP import *
# from force_spring_trap import *
# from force_raised_gaussian import *


# Derived paremeters:
# =============================================================================
gamma0 = 6*np.pi*eta*ro
gamma = np.zeros((3,Np,Nt)) 

gamma[:,range(Np),:]= gamma0

D0 = k_B*T/gamma
mo = rho*(4/3)*np.pi*ro**3      # mass of a polystyrene bead in kg
# =============================================================================


start_time = time.time()
print('\n\n===========================================\n')



   

# Animation function. Reads out the positon coordinates sequentially

def init():
    # Initialization function for the animate function
  
    for np in range(Np):
        ax1.add_patch(beads_xy[np]) 
        ax2.add_patch(beads_yz[np]) 
        ax3.add_patch(beads_xz[np]) 
  
    
    ax1.set_xlabel('$x$ ($\mu$m)')
    ax1.set_ylabel('$y$ ($\mu$m)')
    ax2.set_xlabel('$y$ ($\mu$m)')
    ax2.set_ylabel('$z$ ($\mu$m)')
    ax3.set_xlabel('$x$ ($\mu$m)')
    ax3.set_ylabel('$z$ ($\mu$m)')
    
    
    return beads_xy

def animate(i):
    # Animation function
    # This function is called frame_rate*tfinal times
    # Simulation has Nt steps
    # So, animation should be called every Nt/(frame_rate*tfinal) steps.
    
    fct_adj = int(Nt/(frame_rate*tfinal))
    
    for np in range(Np):
        
        x1 = r_um[0,np,i*fct_adj]
        y1 = r_um[1,np,i*fct_adj]
        z1 = r_um[2,np,i*fct_adj]
        beads_xy[np].center = (x1,y1)
        beads_yz[np].center = (y1,z1)
        beads_xz[np].center = (x1,z1)
        
        # print(i)
    #py.text(0,100,'time, t =',fontsize=12) 
    time_string.set_text(time_template % (i*fct_adj*delta_t))  # Adjust time display by fct_adj factor
    py.sca(ax1)
    draw_geo(i*fct_adj*delta_t, ax1, ax2, ax3)
     
    return beads_xy


def my_plot(fig_no, xp,yp,xlbl,ylbl,lgnd):
    py.figure(fig_no)
    py.plot(xp, yp, label = lgnd)
    py.xlabel(xlbl)
    py.ylabel(ylbl)
    py.legend()

        
    
    







print('Number of particles = %i \n' % Np)
print('Number of time steps = %i \n' % Nt)
print('Final time = %1.2f seconds \n' % tfinal)





# Variables (initialization)
# =============================================================================
t = np.linspace(0,tfinal,Nt)         # time variable for the simulation
r = np.zeros((3,Np,Nt))              # position vector

v = np.zeros((3,Np,Nt))              # velocity vector
v_temp = np.zeros((3,Np))
force_temp = np.zeros((3,Np))
v_drift = np.zeros((3,Np))           # Drift velocity (constant, time independent, same for all particles)
# v_drift[0,:] = -2e-6                 # Drift adjust
w = np.random.normal(0,1,(3,Np,Nt))  # Random white noise term
D = D0            # Diffustion tensor
# =============================================================================


# Dependent parameters 
# =============================================================================
delta_t = t[2]-t[1]     # time step
# fct1 = np.sqrt(2*D*delta_t)
fct11 = np.sqrt(2*D/delta_t)

# fct = np.multiply(fct1,w)
v_br = np.multiply(fct11,w)

fct2 = mo/(gamma*delta_t)
# =============================================================================



# Time evolution of the Langevin equation
# =============================================================================

# r[:,:,0] = np.transpose(np.random.uniform(low = [-xrange_limit,-xrange_limit,max(ro)*1.5], high = [xrange_limit, xrange_limit, zhigh_limit/2], size =(Np,3)))  # Random white noise term

# r[2,:,0] = 0
# r[0,0,0] = 0e-6
# r[1,0,0] = 0e-6
# r[0,1,0] = -60e-6
# r[1,1,0] = 20e-6

r[0,0,0] = 0e-6
r[1,0,0] = -100e-6
r[0,1,0] = 5e-6
r[1,1,0] = 40e-6
r[0,2,0] = 20e-6
r[1,2,0] = 50e-6
r[2,:,0] = 100e-6

print('Solving Langevin equation ... \n')

for m in range(Nt-1):
       
    # force_temp = force_trap(r[:,:,m],0,t[m], Np)
    force_temp = force_profile(r[:,:,m],t[m])
    
    v_temp = 1/(1+fct2[:,:,m])*(fct2[:,:,m]*v[:,:,m] + force_temp/gamma[:,:,m] + v_br[:,:,m])
    # v_temp = v_temp*(-1.0)
    v_temp = particle_collision_adjust(r[:,:,m], v_temp, ro.squeeze(), mo.squeeze(), damping_factor) + v_drift
    v[:,:,m+1] = wall_collision_adjust(r[:,:,m], v_temp, ro.squeeze(), damping_factor)
    # v[:,:,m+1] = v_temp
    r[:,:, m+1] = r[:,:,m] + v[:,:,m+1]*delta_t             # Euler-Cromer method

    



# =============================================================================    
print('Done. \n')
    



r_um = r*1e6   #  um dimensions for plotting
v_um = v*1e6   # velocity in um/s


# Plotting time vs position data and velocity data
py.figure()

for np in range(Np):
    my_plot(1, t,r_um[0, np,:], '$t$ (s)', '$x$ ($\mu$m)', 'Particle %s' % np)         # x-position, partilce np, all time
    my_plot(2, t,r_um[1, np,:], '$t$ (s)', '$y$ ($\mu$m)', 'Particle %s' % np)         # y-position, partilce np, all time
    my_plot(3, t,r_um[2, np,:], '$t$ (s)', '$z$ ($\mu$m)', 'Particle %s' % np)         # z-position, partilce np, all time
    my_plot(4, t,v_um[1, np,:], '$t$ (s)', '$v_y$ ($\mu$m/s)', 'Particle %s' % np)     # vy-velocity, partilce np, all time


# Plotting position scatter data
py.figure()
py.plot(r_um[0,0,:],r_um[1,0,:],'r.')   # x,y position, particle 0, all time

py.xlabel('$x$ ($\mu$m)')
py.ylabel('$y$ ($\mu$m)')


#py.figure(3)



# Animation call
# =============================================================================

print('Processing animation ... \n')

fig = py.figure()
# Define a 2 x 2 grid
gs = py.GridSpec(2,2) # 2 rows, 2 columns


# Define 3 subplots. First one spanning both rows of column 1, the rest two taking the ramining 2 subplots
ax1 = fig.add_subplot(gs[:,0],aspect = 1)       # xy plane (column left, both rows)
ax2 = fig.add_subplot(gs[0,1],aspect = 1)       # yz plane (column right, top row)
ax3 = fig.add_subplot(gs[1,1],aspect = 1)       # xz plane (column right, bottoom row)



# Set figure dpi and size
fig.set_dpi(100)
fig.set_size_inches(14, 6)

# Lower and upper limit of the axes of the three subplots

ax1.set_xlim(-1e6*xrange_limit, 1e6*xrange_limit)
ax1.set_ylim(-1e6*xrange_limit, 1e6*xrange_limit)
ax2.set_xlim(-1e6*xrange_limit, 1e6*xrange_limit)
ax2.set_ylim(1e6*zlow_limit, 1e6*zhigh_limit)
ax3.set_xlim(-1e6*xrange_limit, 1e6*xrange_limit)
ax3.set_ylim(1e6*zlow_limit, 1e6*zhigh_limit)


# Initialize bead patches
beads_xy = []
beads_yz = []
beads_xz = []


for p in range(Np):
     beads_xy.append(plt.patches.Circle((r_um[0,p,0], r_um[1,p,0]), ro[p]*1e6, fc='#C21808',ec = 'k',zorder = 10))
     beads_yz.append(plt.patches.Circle((r_um[1,p,0], r_um[2,p,0]), ro[p]*1e6, fc='#C21808',ec = 'k',zorder = 10))
     beads_xz.append(plt.patches.Circle((r_um[0,p,0], r_um[2,p,0]), ro[p]*1e6, fc='#C21808',ec = 'k',zorder = 10))
    

time_template = 'Time = %.1f s'
time_string = ax1.text(0.05, 0.9, '', transform=ax1.transAxes)





anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=int(frame_rate*tfinal), interval=100, blit=True)



anim.save('br_v1.7.mp4', fps=frame_rate, extra_args=['-vcodec', 'libx264'])




print('Done.\n')
# =============================================================================



print("Execution time = %1.2f seconds \n" % (time.time() - start_time))
print('\n===========================================\n')





