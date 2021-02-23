# -*- coding: utf-8 -*-
"""


@author: Mohammad Asif Zaman
Date: Dec. 2020

Modeling Brownian motion of a colloidal particle


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
-Feb 21. 2021
        - Working on multiparticle system. Working for Np = 2 case
        - Animation for 2 particles added
        - Generalization required
-Feb 22, 2021
        - Array index order changed to [space dimension, Np, Nt]
        - Animation working for arbitrary Np values
        - The force_interp() function needs to be rewritten. force_model() and force_rand() has been modified.
        - Uploaded as v1.4
        

"""


from __future__ import print_function    

import time
import math
import numpy as np
import pylab as py

import matplotlib as plt

from matplotlib import animation, rc


# Coordinate system notes: 
#    x is the 3D position vector. 
#   r[0] = x component 
#   r[1] = y component
#   r[2] = z component


# Read force data from data file
Mdata = np.genfromtxt('xsweep_ro=10u,zo=25u.csv',delimiter=',',skip_header=5)
xdata = Mdata[:,0]*1e-6
Fxdata = Mdata[:,3]*1e-12
Fydata = 6*Mdata[:,2]*1e-12
Fzdata = Mdata[:,4]*1e-12
# Note that the axis has been changed. The axis from the csv data file is different.


# Animation function. Reads out the positon coordinates sequentially

def init():
    # Initialization function for the animate function
    # Draw electrodes
    for kk in range(-2,3):
        rectangle = py.Rectangle((-xplt_limit*1e6/2, -elec_width/2+kk*elec_spacing),xplt_limit*1e6,elec_width,fc='b')
        py.gca().add_patch(rectangle)
  
    for np in range(Np):
        ax.add_patch(beads[np])    
    
    # patch1.center = (0, 0)
    # ax.add_patch(patch1)
    
    # patch2.center = (0, 0)
    # ax.add_patch(patch2)
    
    ax.set_xlabel('$x$ ($\mu$m)')
    ax.set_ylabel('$y$ ($\mu$m)')
    
    return beads

def animate(i):
    # Animation function
    # x1, y1 = patch1.center
    
    

    for np in range(Np):
        
        x1 = r_um[0,np,i]
        y1 = r_um[1,np,i]
        beads[np].center = (x1,y1)
    
    
        
    #py.text(0,100,'time, t =',fontsize=12) 
    time_string.set_text(time_template % (i*delta_t))
    
     
    return beads


# def animation_main(i):
#     for np in range(0,Np):
#         animate(i,np)
    
#     return []
    

# Simplified spring force model
def force_model(rm,xactive):
    k =  .3e-6    # spring constant
    r_mag = 0  # magnitude of the random force component
    fm = np.zeros((3,Np))
    
    fm[0,:] = 0
    fm[1,:] = -k*(rm[1,:]-xactive)
    fm[2,:] = 0     
    return fm

# Random force component
def force_random(rr,xactive,t):
    # This force component comes from random origins with random magnitude
    k =  .3e-6    # spring constant
    r_mag = 1/50  # magnitude of the random force component
    fr = np.zeros((3,Np))
    if t> 0.1 and 200*t%1 == 0:
        fr[0,:] = -5*r_mag*k*(rr[0,:]-np.random.normal(0,1,Np)*xplt_limit)*np.random.normal(0,1,Np)
        fr[1,:] = -r_mag*k*(rr[1,:]-np.random.normal(0,1,Np)*xplt_limit)*np.random.normal(0,1,Np)     # random y force 1/20th the magnitude of the x component
        fr[2,:] = 0     # setting the z component of the force to be 0
    else:
        fr[0,:] = 0
        fr[1,:] = 0
        fr[2,:] = 0
    
    
    return fr


# Interpolation function when using force value from data file    
# def force_interp(x,xactive):
#     # Interpolate function for the imported force data(csv file)
#     f = np.zeros(3)
#     f[0] = np.interp((r[0,1]-xactive),xdata,Fxdata)  
#     f[1] = np.interp((r[0,1]-xactive),xdata,Fydata)
#     f[2] = np.interp((r[0,1]-xactive),xdata,Fzdata)
#     return f



# Force profile centered around the origin    
def force_origin(x,xactive,t):
    # return force_interp(x,xactive) + force_random(x,xactive,t)
    return force_model(x,xactive) + force_random(x,xactive,t)



# Time dependent force field. Implementation of the switching sequence
def force_switch(x,t):
    
    # define switching instances
   
    ts = np.linspace(0,38,20)  
  
    # define active electrode position at switching
    x_e5 = elec_spacing*2e-6
    x_e4 = elec_spacing*1e-6
    x_e3 = elec_spacing*0
    x_e2 = -elec_spacing*1e-6
    x_e1 = -elec_spacing*2e-6
   
    xactive = x_e5 

    # Assigning active electrode location based on time. 
    # Note the a if statement can be overridden by the subsequent if statement    
    if t < ts[0]:
        xactive = x_e5
    if t >= ts[1]:
        xactive = x_e4
    if t >= ts[2]:
        xactive = x_e3
    if t >= ts[3]:
        xactive = x_e2
    if t >= ts[4]:
        xactive = x_e1
    if t >= ts[5]:
        xactive = x_e2
    if t >= ts[6]:
        xactive = x_e3
    if t >= ts[7]:
        xactive = x_e4
    if t >= ts[8]:
        xactive = x_e5
    if t >= ts[9]:
        xactive = x_e4
    if t >= ts[10]:
        xactive = x_e3
    if t >= ts[11]:
        xactive = x_e2
    if t >= ts[12]:
        xactive = x_e1
    if t >= ts[13]:
        xactive = x_e2
    if t >= ts[14]:
        xactive = x_e3
    if t >= ts[15]:
        xactive = x_e4
    if t >= ts[16]:
        xactive = x_e5
    if t >= ts[17]:
        xactive = x_e4
    if t >= ts[18]:
        xactive = x_e3
    if t >= ts[19]:
        xactive = x_e2

    
    
    
    return force_origin(x,xactive,t)
        
    
    

# Parameters (All SI units)
# =============================================================================
k_B = 1.38e-23
T = 300
eta = 8.9e-4
ro = 10e-6;
gamma = 6*np.pi*eta*ro
D0 = k_B*T/gamma
tfinal = 38
Nt =  1501   # Number of time steps
Np = 3       # Number of particles
#x_trap = 70e-6  # center of the trap for the test function
xplt_limit = 250e-6;


# =============================================================================



# Electrode array geometry (units in um)
# =============================================================================
elec_width = 15
elec_spacing = 50
# =============================================================================



# Variables (initialization)
# =============================================================================
t = np.linspace(0,tfinal,Nt)      # time variable for the simulation
r = np.zeros((3,Np,Nt))              # position vector
v = np.zeros((3,Np,Nt))              # velocity vector
v[0,:,:] = v[0,:,:] #- 2e-6            # Drift adjust
w = np.random.normal(0,1,(3,Np,Nt))  # Random white noise term
D = np.ones((3,Np,Nt))*D0            # Diffustion tensor
# =============================================================================


# Dependent parameters 
# =============================================================================
delta_t = t[2]-t[1]     # time step
fct1 = np.sqrt(2*D*delta_t)
fct = np.multiply(fct1,w)
# =============================================================================



# Time evolution of the Langevin equation
# =============================================================================
r[0,0,0] = 100e-6
r[1,0,0] = -100e-6
r[0,1,0] = -60e-6
r[1,1,0] = 40e-6



for m in range(Nt-1):
    r[:,:, m+1] = r[:,:,m] + fct[:,:,m] + delta_t*force_switch(r[:,:,m],t[m])/gamma + v[:,:,m]*delta_t
# =============================================================================    
    


r_um = r*1e6   #  um dimensions for plotting
#print(fct)
#print(x)


# Plotting time vs position data
py.figure(1)
#py.plot(t,x_um[:,0])
for np in range(Np):
    py.plot(t,r_um[1, np,:], label = 'particle %s' % np)     # y-position, partilce np, all time

# py.plot(t,x_um[:,2])
py.xlabel('$t$ (s)')
py.ylabel('$y$ ($\mu$m)')
py.legend()

# Plotting position scatter data
py.figure(2)
py.plot(r_um[0,0,:],r_um[1,0,:],'r.')   # x,y position, particle 0, all time

py.xlabel('$x$ ($\mu$m)')
py.ylabel('$y$ ($\mu$m)')


#py.figure(3)



# Animation call
# =============================================================================


fig = py.figure()
fig.set_dpi(100)
fig.set_size_inches(7, 7)

ax = py.axes(xlim=(-1e6*xplt_limit, 1e6*xplt_limit), ylim=(-1e6*xplt_limit, 1e6*xplt_limit))
# patch1 = py.Circle((0, 0), ro*1e6, fc='r')
# patch2 = py.Circle((0, 0), ro*1e6, fc='r')
beads = []

for np in range(Np):
     beads.append(plt.patches.Circle((r_um[0,np,0], r_um[1,np,0]), ro*1e6, fc='r'))
    

time_template = 'Time = %.1f s'
time_string = ax.text(0.05, 0.9, '', transform=ax.transAxes)

anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=Nt, interval=100, blit=True)
anim.save('br_v1.4.mp4', fps=Nt/tfinal, extra_args=['-vcodec', 'libx264'])

# =============================================================================




