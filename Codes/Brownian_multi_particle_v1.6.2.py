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
- March7, 2021
        - Added zorder for animating the beads (beads append). This ensures that the beads are always in the foreground (compared to source geometry)        
        - Added draw_source(t) function inside the animate() function. Also, made the draw_source() as a funciton of time for dynamic
          manipulation of the source geometry (e.g. turning ON/OFF optical excitation for optical trapping demo.)
        
        
"""


from __future__ import print_function    

import time
import math
import numpy as np
import pylab as py

import matplotlib as plt

from matplotlib import animation, rc

from collision_adjust import *


from forceDEP import *
# from force_spring_trap import *
# from force_raised_gaussian import *



start_time = time.time()
print('\n\n===========================================\n')



   

# Animation function. Reads out the positon coordinates sequentially

def init():
    # Initialization function for the animate function
  
    for np in range(Np):
        ax.add_patch(beads[np])    
  
    
    ax.set_xlabel('$x$ ($\mu$m)')
    ax.set_ylabel('$y$ ($\mu$m)')
    
    return beads

def animate(i):
    # Animation function
    # This function is called frame_rate*tfinal times
    # Simulation has Nt steps
    # So, animation should be called every Nt/(frame_rate*tfinal) steps.
    
    fct_adj = int(Nt/(frame_rate*tfinal))
    
    for np in range(Np):
        
        x1 = r_um[0,np,i*fct_adj]
        y1 = r_um[1,np,i*fct_adj]
        beads[np].center = (x1,y1)
        # print(i)
    #py.text(0,100,'time, t =',fontsize=12) 
    time_string.set_text(time_template % (i*fct_adj*delta_t))  # Adjust time display by fct_adj factor
    draw_source(i*fct_adj*delta_t)
     
    return beads




        
    
    

# Parameters (All SI units)
# =============================================================================
k_B = 1.38e-23
T = 300
eta = 8.9e-4
ro = 10e-6;
gamma = 6*np.pi*eta*ro
D0 = k_B*T/gamma


frame_rate = 30   # Animation frame rate in fps

tfinal =  38  #38 # 6 # 38
Nt = frame_rate*40 # 300 #1501   # Number of time steps
Np = 3       # Number of particles
# xplt_limit = 250e-6;


rho = 1055                      # density of polystyrene beads in kg/m3
mo = rho*(4/3)*np.pi*ro**3      # mass of a polystyrene bead in kg
damping_factor = 0              # Fraction of the velocity that is lost during each collision


# =============================================================================





print('Number of particles = %i \n' % Np)
print('Number of time steps = %i \n' % Nt)
print('Final time = %1.2f seconds \n' % tfinal)





# Variables (initialization)
# =============================================================================
t = np.linspace(0,tfinal,Nt)      # time variable for the simulation
r = np.zeros((3,Np,Nt))              # position vector

v = np.zeros((3,Np,Nt))              # velocity vector
v_temp = np.zeros((3,Np))
force_temp = np.zeros((3,Np))
v_drift = np.zeros(3)                # Drift velocity (constant, time independent, same for all particles)
#v_drift[0] = -2e-6                   # Drift adjust
w = np.random.normal(0,1,(3,Np,Nt))  # Random white noise term
D = np.ones((3,Np,Nt))*D0            # Diffustion tensor
# =============================================================================


# Dependent parameters 
# =============================================================================
delta_t = t[2]-t[1]     # time step
fct1 = np.sqrt(2*D*delta_t)
fct11 = np.sqrt(2*D/delta_t)

fct = np.multiply(fct1,w)
v_br = np.multiply(fct11,w)

fct2 = mo/(gamma*delta_t)
# =============================================================================



# Time evolution of the Langevin equation
# =============================================================================

# r[:,:,0] = np.random.normal(0,xrange_limit/2,(3,Np))  # Random white noise term
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

print('Solving Langevin equation ... \n')

for m in range(Nt-1):
    # v_temp[:,:,m] = fct[:,:,m]/delta_t + force_movingDEP(r[:,:,m],t[m], elec_spacing,Np)/gamma 
    # v_temp[:,:,m] = fct[:,:,m]/delta_t + force_trap(r[:,:,m],0,t[m], Np)/gamma 
    
    # tmp = (fct[:,:,m]/delta_t + force_trap(r[:,:,m],0,t[m], Np)/gamma - velocity_adjust(r[:,:,m],fct[:,:,m]/delta_t + force_trap(r[:,:,m],0,t[m], Np)/gamma, ro)).max()
    # print(tmp)
    
    # force_temp = force_trap(r[:,:,m],0,t[m], Np)
    force_temp = force_profile(r[:,:,m],t[m])
    
    v_temp = 1/(1+fct2)*(fct2*v[:,:,m] + force_temp/gamma + v_br[:,:,m])
    # v_temp = v_temp*(-1.0)
    v[:,:,m+1] = velocity_adjust(r[:,:,m], v_temp, ro, damping_factor)
    r[:,:, m+1] = r[:,:,m] + v[:,:,m+1]*delta_t             # Euler-Cromer method

    # r[:,:, m+1] = r[:,:,m] + (v[:,:,m] + v[:,:,m+1])*delta_t*0.5 
    # r[2,:,:] = 0    # Fix z position to 0 to limit motion in 2D only
    
    # v[:,:,m] = velocity_adjust(r[:,:,m], fct[:,:,m]/delta_t + force_trap(r[:,:,m],0,t[m], Np)/gamma, ro)
    #v[:,:,m] = velocity_adjust(r[:,:,m],fct[:,:,m]/delta_t + force_movingDEP(r[:,:,m],t[m], elec_spacing,Np)/gamma, ro)
    
    # r[:,:, m+1] = r[:,:,m] + v[:,:,m]*delta_t 

    
 



# =============================================================================    
print('Done. \n')
    


r_um = r*1e6   #  um dimensions for plotting
v_um = v*1e6   # velocity in um/s
#print(fct)
#print(x)


# Plotting time vs position data
py.figure(1)
#py.plot(t,x_um[:,0])
for np in range(Np):
    py.plot(t,r_um[1, np,:], label = 'Particle %s' % np)     # y-position, partilce np, all time

# py.plot(t,x_um[:,2])
py.xlabel('$t$ (s)')
py.ylabel('$y$ ($\mu$m)')
py.legend()


py.figure()
#py.plot(t,x_um[:,0])
for np in range(Np):
    py.plot(t,r_um[2, np,:], label = 'Particle %s' % np)     # y-position, partilce np, all time

# py.plot(t,x_um[:,2])
py.xlabel('$t$ (s)')
py.ylabel('$z$ ($\mu$m)')
py.legend()


# Plotting time vs position data
py.figure()
#py.plot(t,x_um[:,0])
for np in range(Np):
    py.plot(t,v_um[1, np,:], label = 'Particle %s' % np)     # y-position, partilce np, all time

# py.plot(t,x_um[:,2])
py.xlabel('$t$ (s)')
py.ylabel('$v_y$ ($\mu$m/s)')
py.legend()



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
fig.set_dpi(100)
fig.set_size_inches(7, 7)

ax = py.axes(xlim=(-1e6*xrange_limit, 1e6*xrange_limit), ylim=(-1e6*xrange_limit, 1e6*xrange_limit))
beads = []

for np in range(Np):
     beads.append(plt.patches.Circle((r_um[0,np,0], r_um[1,np,0]), ro*1e6, fc='#C21808',ec = 'k',zorder = 10))
    

time_template = 'Time = %.1f s'
time_string = ax.text(0.05, 0.9, '', transform=ax.transAxes)


# draw_source()



anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=int(frame_rate*tfinal), interval=100, blit=True)



anim.save('br_v1.6.mp4', fps=frame_rate, extra_args=['-vcodec', 'libx264'])


# anim = animation.FuncAnimation(fig, animate, init_func=init,
                                # frames=Nt, interval=100, blit=True)
# anim.save('br_v1.4.mp4', fps=Nt/tfinal, extra_args=['-vcodec', 'libx264'])

print('Done.\n')
# =============================================================================



print("Execution time = %1.2f seconds \n" % (time.time() - start_time))
print('\n===========================================\n')





