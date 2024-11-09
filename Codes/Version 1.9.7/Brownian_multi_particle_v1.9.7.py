# -*- coding: utf-8 -*-
"""


@author: Mohammad Asif Zaman
Date: Dec. 2020

Modeling Brownian motion of a colloidal particle



Coordinate system notes: 
    - r[3,Np,Nt] is the 3D position vector array for Np particles at Nt time instances (Cartesian). 
    - r[0,:,:] = x components 
    - r[1,:,:] = y components
    - r[2,:,:] = z components


External function (defined in other modules) locations:
    - time_pos_ax_limits()                          -- geometry_def_X.py
    - geo_define()                                  -- geometry_def_X.py
    - particle_color(a1)                            -- geometry_draw_X.py 
    - force_profile(a1,a2)                          -- force_X.py module
    - fluid_vel(a1,a2)                              -- force_X.py module
    - particle_collision_adjust(a1,a2,a3,a4,a5)     -- particle_particle_collision_adjust.py
    - wall_collision_adjust(a1,a2,a3,a4,a5,a6,a7)   -- particle_wall_collision_adjust.py




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
- March 24-25, 2021
        - Updated the Diffusion tensor + hydrodynamic interactions
        - Redefined ro as a (Np,) vector instead of (Np,1) vector
        - Updated plot font size
        - Released as version 1.7.2
- March 29, 2021
        - Bug fix: matplotlib version issues: Replaced py.sca() to py.gcf().sca()
- March 31, 2021
        - Added feature to plot specific frames of the animation separately (useful for saving images to file)
- April 10, 2021
        - Added parameters defining initial particle position range (xi_lim, yi_lim and zi_lim) inside the force functions
        - Extra text feature on the animation (text_string1 and text_string2)
- April 11, 2021
        - Improved axes limit parameter definitions
- May 3, 2021
        - Added fluid-velocity term
- May 6, 2021
        - Added particle_color() function         
- May 3-8, 2021      
        - Generalized wall collision dynamics
- May 9-27, 2021
        - Generalized node based geometry definition 
        - Integration of geometry drawing with geometry definition
        - Integration of collision dynamics with geometry definition
        - Reducing the number of particle data plots to reduce clutter
- May 28, 2021
        - Functionalized fluid velocity        
- Aug 11, 2022
        - Nt definition bug. Nt = frame_rate*40 was incorrect. It should also be multiplied by tfinal
        - sca() function may cause issues. It might be good be replace it with old version 
          explicitly defined within geometry_draw module
"""


from __future__ import print_function    

import time
import math
import numpy as np  
import pylab as py

import matplotlib as plt
from matplotlib import animation, rc

py.close('all')


from parameters import *

from particle_particle_collision_adjust import *
from particle_wall_collision_adjust import *


from force_OT import *
from geometry_def_OT import *
from geometry_draw_OT import *


# from force_sorter import *
# from geometry_def_sorter import *
# from geometry_draw_sorter import *


# from force_DEP import *
# from geometry_def_DEP import *
# from geometry_draw_DEP import *



# from force_EL import *
# from geometry_def_EL import *
# from geometry_draw_EL import *

# from force_k2 import *
# from geometry_def_k2 import *
# from geometry_draw_k2 import *







# Control parameters
# =============================================================================
save_video = 'y'
save_data = 'y'
data_thin_fct = 1     # at what step is the position data saved in the csv file


# Time, position and axis paremeters:
# =============================================================================
Np, ro, tfinal, x_lim, y_lim, z_lim, xi_lim, yi_lim, zi_lim = time_pos_ax_limits()   # Defined in geometry_def_X.py module
Nt =frame_rate*tfinal*10      # 40 # 300 #1501   # Number of time steps
# Geometry paremeters:
# =============================================================================
n_vector_set, d_set, p_lim_set = geo_define() # Defined in geometry_def_X.py module



# Derived paremeters:
# =============================================================================
gamma0 = 6*np.pi*eta*ro         # friction coefficient (Np,) array
D0 = k_B*T/gamma0               # Free space diffusion coefficien. (Np,) array
mo = rho*(4/3)*np.pi*ro**3      # mass of a polystyrene bead in kg. (Np,) array
scale_fct = 1

# =============================================================================


start_time = time.time()
print('\n\n===========================================\n')


py.rcParams['font.size'] = 14
   



# =============================================================================
# Animation function. Reads out the positon coordinates sequentially

def init():
    # Initialization function for the animate function
    
    # Each axis represents a subplot within the animation figure
    for np in range(Np):
        ax_xy.add_patch(beads_xy[np]) 
        ax_yz.add_patch(beads_yz[np]) 
        ax_xz.add_patch(beads_xz[np]) 
  
    
    # Set axis labels
    ax_xy.set_xlabel('$x$ ($\mu$m)')
    ax_xy.set_ylabel('$y$ ($\mu$m)')
    ax_yz.set_xlabel('$y$ ($\mu$m)')
    ax_yz.set_ylabel('$z$ ($\mu$m)')
    ax_xz.set_xlabel('$x$ ($\mu$m)')
    ax_xz.set_ylabel('$z$ ($\mu$m)')
    
    
    return beads_xy
# =============================================================================




# =============================================================================
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
        
      
    time_string.set_text(time_template % (i*fct_adj*delta_t))  # Adjust time display by fct_adj factor
    
    # py.sca(ax_xy)
    str1, str2 = draw_geo(i*fct_adj*delta_t, ax_xy, ax_yz, ax_xz)
    text_string1.set_text(str1)  
    # text_string2.set_text(str2)  
    
    
    return beads_xy
# =============================================================================


def anim_setup():
    fig = py.figure()
    # Define a 2 x 2 grid
    gs = py.GridSpec(2,2) # 2 rows, 2 columns
    
    
    
    # Define 3 subplots. First one spanning both rows of column 1, the rest two taking the ramining 2 subplots
    # Original
    ax_xy = fig.add_subplot(gs[:,0],aspect = 1)       # xy plane (column left, both rows)
    ax_yz = fig.add_subplot(gs[0,1],aspect = 1)       # yz plane (column right, top row)
    ax_xz = fig.add_subplot(gs[1,1],aspect = 1)       # xz plane (column right, bottoom row)
   
    # Set figure dpi and size
    fig.set_dpi(150)
    fig.set_size_inches(14, 9)
    
    # Nov. 1, 2024: Vertical stack (short format)
    # gs = py.GridSpec(3,1) # 3 rows, 1 columns
    # ax_xy = fig.add_subplot(gs[0,0],aspect = 1)       # xy plane (column 1, row 1)
    # ax_yz = fig.add_subplot(gs[1,0],aspect = 1)       # yz plane (column 1, row 2)
    # ax_xz = fig.add_subplot(gs[2,0],aspect = 1)       # xz plane (column 1, row 3)
    # fig.set_dpi(150)
    # fig.set_size_inches(9,16)

    
    
    # Lower and upper limit of the axes of the three subplots
    
    ax_xy.set_xlim(x_lim[0], x_lim[1])
    ax_xy.set_ylim(y_lim[0], y_lim[1])
    ax_yz.set_xlim(y_lim[0], y_lim[1])
    ax_yz.set_ylim(z_lim[0], z_lim[1])
    ax_xz.set_xlim(x_lim[0], x_lim[1])
    ax_xz.set_ylim(z_lim[0], z_lim[1])
    
    
    # Initialize bead patches
    beads_xy = []
    beads_yz = []
    beads_xz = []
    
    # Nov. 6, 2024
    # To find front/back order when drawing/plotting beads, it is necessary to find the lowest/highest position of the bead
    # along different axis directions. For example, in xy plane, the particle with the highest z value should appear on top of 
    # all other particles. 
    ind_sort_x = np.argsort(r_um[0,:,0])  # sorting index from lowest x position to highest x position 
    ind_sort_y = np.argsort(r_um[1,:,0])  # sorting index from lowest y position to highest y position
    ind_sort_z = np.argsort(r_um[2,:,0])  # sorting index from lowest z position to highest z position
    
    for p in range(Np):
         cl_p = particle_color(p)   # Defined in geometry_draw_X.py module
         
         # Nov. 6, 2024: Note the zorder modification. Take a default value of 10 and add the sorting index of the particle (along
         # (x, y or z direction). This should put all the particles in order. The order being, the particle with position value
         # largest in the axis direction perpendicular to the plane is plotted on the front most. 
         beads_xy.append(plt.patches.Circle((r_um[0,p,0], r_um[1,p,0]), ro[p]*1e6, fc=cl_p,ec = 'k',zorder = 10 + ind_sort_z[p]))
         beads_yz.append(plt.patches.Circle((r_um[1,p,0], r_um[2,p,0]), ro[p]*1e6, fc=cl_p,ec = 'k',zorder = 10 + ind_sort_x[p]))
         beads_xz.append(plt.patches.Circle((r_um[0,p,0], r_um[2,p,0]), ro[p]*1e6, fc=cl_p,ec = 'k',zorder = 10 + Np -  ind_sort_y[p])) # reverse the order for xz axis
        
    
    time_template = 'Time = %.1f s'
    time_string  = ax_xy.text(0.05, 0.90, '', transform=ax_xy.transAxes)
    text_string1 = ax_xy.text(0.05, 0.85, '', transform=ax_xy.transAxes)
    text_string2 = ax_xy.text(0.05, 0.80, '', transform=ax_xy.transAxes)
    
    return fig,ax_xy, ax_yz, ax_xz, beads_xy, beads_yz, beads_xz, time_template, time_string, text_string1, text_string2
    


# =============================================================================
def my_plot(fig_no, xp,yp,xlbl,ylbl,lgnd):
    py.figure(fig_no)
       
    py.plot(xp, yp, label = lgnd, )
    py.xlabel(xlbl)
    py.ylabel(ylbl)
    py.legend(prop={'size': 10})
# =============================================================================

        
    
    

# =============================================================================
def Diffusion_tensor(r_in, ro_in):
    # r_in is a (3,Np) vector
    # ro_in is (Np,) vector
    
    D_tensor = np.zeros((3,Np))
    
    z_in = r_in[2,:]    # should be a (Np,) vector
    ro_z = ro_in/z_in   # (Np,) vector
    
    H_parallel = 1 - (9/16)*(ro_z) + (1/8)*(ro_z)**3 - (45/256)*(ro_z)**4 - (1/16)*(ro_z)**5
    H_perp = (6*z_in**2 + 2*ro_in*z_in)/(6*z_in**2 + 9*ro_in*z_in + 2*ro_in**2)
    
    # (k_B*T/(6*np.pi*eta*ro_in))
    D_tensor[0,:] = D0 * H_parallel
    D_tensor[1,:] = D0 * H_parallel
    D_tensor[2,:] = D0 * H_perp
    
    
    return D_tensor
# =============================================================================
    




print('Number of particles = %i \n' % Np)
print('Number of time steps = %i \n' % Nt)
print('Final time = %1.2f seconds \n' % tfinal)





# Variables (initialization)
# =============================================================================
t = np.linspace(0,tfinal,Nt)         # time variable for the simulation, Nt defined in parameters.py
r = np.zeros((3,Np,Nt))              # position vector

v = np.zeros((3,Np,Nt))              # velocity vector
a = np.zeros((3,Np,Nt))              # acceleration vector
v_temp = np.zeros((3,Np))

force_ext = np.zeros((3,Np))
v_fluid = np.zeros((3,Np))           # Fluid-flow velocity (same for all particles)
           
w = np.random.normal(0,1,(3,Np,Nt))  # Random white noise term
# D = D0            # Diffustion tensor
# =============================================================================


# Dependent parameters 
# =============================================================================
delta_t = t[2]-t[1]          # time step
delta_t_us = delta_t*1e6     # time step in micro seconds

# =============================================================================


print('Time step = %1.8f us \n' % delta_t_us)


  
# Time evolution of the Langevin equation
# =============================================================================

r[:,:,0] = np.transpose(np.random.uniform(low = [xi_lim[0], yi_lim[0], zi_lim[0]], high = [xi_lim[1], yi_lim[1], zi_lim[1]], size =(Np,3)))  

# Forcing initial positions to specific y values for the particle sorter program
# r[0,0,0] = 500e-9
# r[1,0,0] = -20e-9
# r[2,0,0] = 200e-9


# v[0,:,0] = np.zeros(Np) - 100e-6 


print('Solving Langevin equation ... \n')

for m in range(Nt-1):
    D = Diffusion_tensor(r[:,:,0],ro.squeeze())        # (3,Np) vector array. Function defined in this file.
    
    # All of the following factors are (3,Np) vector array
    Lambda = (mo*D)/(k_B*T*delta_t)  
    # Lambda = 0*D + 1e2        
              
    fct_1 = np.sqrt(2*D/delta_t) * scale_fct
    fct_2 = D/(k_B*T) 
    
    # force_temp = force_trap(r[:,:,m],0,t[m], Np)
    force_ext = force_profile(r[:,:,m],t[m])         # Function defined in force_X.py module
    v_fluid   = fluid_vel(r[:,:,m],t[m])             # Function defined in force_X.py module
    # v_fluid[1,:] = -10e-6 
    
    v_temp = (Lambda/(1+Lambda)) * v[:,:,m] + (1/(1+Lambda)) * ( v_fluid + fct_1 * w[:,:,m]  +  fct_2 * force_ext )
    
    # Take into accoutn particle-particle collisions. Aug 15, 2022: This is only needed when there is more than one particle.
    if Np > 1:
        v_temp = particle_collision_adjust(r[:,:,m], v_temp, ro, mo, damping_factor)    # Function defiend in particle_particle_collsion_adjust.py module
    
    # Take into accoutn particle-wall collisions
    v[:,:,m+1] = wall_collision_adjust(r[:,:,m], v_temp, ro, damping_factor, n_vector_set, d_set, p_lim_set) # Function defiend in particle_wall_collsion_adjust.py module
    
    a[:,:,m+1] = (v[:,:,m+1] - v[:,:,m])/delta_t            # Acceleration
    
    # Update position
    # v[:,:,m+1] = v_temp
    r[:,:, m+1] = r[:,:,m] + v[:,:,m+1]*delta_t             # Euler-Cromer method

    
# =============================================================================    

print('Done. \n')
max_diff_x = max(r[0,0,1:-1] - r[0,0,0:-2])*1e6
max_diff_y = max(r[1,0,1:-1] - r[1,0,0:-2])*1e6
max_diff_z = max(r[2,0,1:-1] - r[2,0,0:-2])*1e6

print('Max x displacement = %1.2f micron \n' % max_diff_x)
print('Max y displacement = %1.2f micron \n' % max_diff_y)
print('Max z displacement = %1.2f micron \n' % max_diff_z)



# Data plotting
# =============================================================================

r_um = r*1e6   #  um dimensions for plotting
v_um = v*1e6   # velocity in um/s
a_um = a*1e6   # acceleration in um/s^2

# Plotting time vs position data and velocity data
py.figure()

# May 27, 2021 and June 2, 2021
# Plot N_plot_max particle data or less (with large number of particles, plotting all of them becomes crowded)
# N_plot_max defined in module parameters.py
for p in range(min(Np,N_plot_max)):
    my_plot(1, t,r_um[0, p,:], '$t$ (s)', '$x$ ($\mu$m)', 'Particle %s' % p)         # x-position, partilce np, all time
    my_plot(2, t,r_um[1, p,:], '$t$ (s)', '$y$ ($\mu$m)', 'Particle %s' % p)         # y-position, partilce np, all time
    my_plot(3, t,r_um[2, p,:], '$t$ (s)', '$z$ ($\mu$m)', 'Particle %s' % p)         # z-position, partilce np, all time
    my_plot(4, t,v_um[0, p,:], '$t$ (s)', '$v_x$ ($\mu$m/s)', 'Particle %s' % p)     # vx-velocity, partilce np, all time
    my_plot(5, t,v_um[1, p,:], '$t$ (s)', '$v_y$ ($\mu$m/s)', 'Particle %s' % p)     # vy-velocity, partilce np, all time
    my_plot(6, t,a_um[1, p,:], '$t$ (s)', '$a_y$ ($\mu$m/s)', 'Particle %s' % p)     # ay-velocity, partilce np, all time
    # my_plot(6, t,a[1, p,:]*mo[p]*1e12, '$t$ (s)', '$m_o a_y$ (pN)', 'Particle %s' % p)     # ay-velocity, partilce np, all time
    
# Plotting position scatter data
py.figure()
py.plot(r_um[0,0,:],r_um[1,0,:],'r.')   # x,y position, particle 0, all time

py.xlabel('$x$ ($\mu$m)')
py.ylabel('$y$ ($\mu$m)')
# =============================================================================



# Animation call
# =============================================================================

if save_video == 'y':
    
    print('Processing animation ... ### \n')
    
    fig_anim, ax_xy, ax_yz, ax_xz, beads_xy, beads_yz, beads_xz, time_template, time_string, text_string1, text_string2 = anim_setup()   # Function defined in this file
    
    anim = animation.FuncAnimation(fig_anim, animate, init_func=init,
                                frames=int(frame_rate*tfinal), interval=100, blit=True)


    print('Saving video file ... \n')
    anim.save('br_v1.9.mp4', fps=frame_rate, extra_args=['-vcodec', 'libx264'])
    print('Video file saved. \n')
    
    print('Animation processing done. ### \n')





# =============================================================================




# # Plot a specific frame of the animation
# # =============================================================================
# # t_snap_set = np.array([0.2, 0.8, 1.4, 2, 3, 4, 4.5])
# t_snap_set = np.array([0])


# for m in range(t_snap_set.shape[0]):
#     t_snap = t_snap_set[m]
#     print('Plotting time instance close to = %1.1f ... ' % (t_snap))
#     fig_ex, ax_xy, ax_yz, ax_xz, beads_xy, beads_yz, beads_xz, time_template, time_string,text_string1, text_string2 = anim_setup()
#     draw_static_geo(ax_xy, ax_yz, ax_xz)
#     init()
#     animate(int(t_snap*frame_rate))
#     print('Done.\n')



# # =============================================================================
if save_data == 'y':
    print('Saving data file... \n')
    np.savetxt('p0_position.csv', np.transpose(r_um[:,0,0:-1:data_thin_fct]),delimiter=",")
    np.savetxt('p0_velocity.csv', np.transpose(v_um[:,0,0:-1:data_thin_fct]),delimiter=",")
    np.savetxt('t.csv', t[0:-1:data_thin_fct],delimiter=",")
    print('Data file saved. \n')


print("Execution time = %1.2f seconds \n" % (time.time() - start_time))
print('\n===========================================\n')





