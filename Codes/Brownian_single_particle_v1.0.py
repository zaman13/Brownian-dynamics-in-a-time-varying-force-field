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

"""


from __future__ import print_function    

import time
import math
import numpy as np
import pylab as py

from matplotlib import animation, rc


# Coordinate system notes: 
#    x is the 3D position vector. 
#   x[0] = x component 
#   x[1] = y component
#   x[2] = z component


# Read force data from data file
# Mdata = np.genfromtxt('xsweep_ro=10u,zo=25u.csv',delimiter=',',skip_header=5)
# xdata = Mdata[:,0]*1e-6
# Fxdata = Mdata[:,3]*1e-12
# Fydata = 6*Mdata[:,2]*1e-12
# Fzdata = Mdata[:,4]*1e-12
# Note that the axis has been changed. The axis from the csv data file is different.


# Animation function. Reads out the positon coordinates sequentially

def init():
    # Initialization function for the animate function
    # Draw electrodes
    for kk in range(-2,3):
        rectangle = py.Rectangle((-xplt_limit*1e6/2, -elec_width/2+kk*elec_spacing),xplt_limit*1e6,elec_width,fc='b')
        py.gca().add_patch(rectangle)
  
    
    patch.center = (0, 0)
    ax.add_patch(patch)
    ax.set_xlabel('$x$ ($\mu$m)')
    ax.set_ylabel('$y$ ($\mu$m)')
    
    return patch,

def animate(i):
    # Animation function
    x, y = patch.center    
    x = x_um[i,0]
    y = x_um[i,1]
    patch.center = (x, y)
    #py.text(0,100,'time, t =',fontsize=12) 
    time_string.set_text(time_template % (i*delta_t))
    return patch,




# Simplified spring force model
def force_model(x,xactive):
    k =  .3e-6    # spring constant
    r_mag = 0  # magnitude of the random force component
    f = np.zeros(3)
    f[1] = -k*(x[1]-xactive)
    f[0] = 0
    f[2] = 0     
    return f

# Random force component
def force_random(x,xactive,t):
    # This force component comes from random origins with random magnitude
    k =  .3e-6    # spring constant
    r_mag = 1/30  # magnitude of the random force component
    f = np.zeros(3)
    if t> 0.1 and 200*t%1 == 0:
        f[0] = -5*r_mag*k*(x[0]-np.random.normal(0,1)*xplt_limit)*np.random.normal(0,1)
        f[1] = -r_mag*k*(x[1]-np.random.normal(0,1)*xplt_limit)*np.random.normal(0,1)     # random y force 1/20th the magnitude of the x component
        f[2] = 0     # setting the z component of the force to be 0
    else:
        f[0] = 0
        f[1] = 0
        f[2] = 0
    
    
    return f


# Interpolation function when using force value from data file    
def force_interp(x,xactive):
    # Interpolate function for the imported force data(csv file)
    f = np.zeros(3)
    f[0] = np.interp((x[1]-xactive),xdata,Fxdata)  
    f[1] = np.interp((x[1]-xactive),xdata,Fydata)
    f[2] = np.interp((x[1]-xactive),xdata,Fzdata)
    return f



# Force profile centered around the origin    
def force_origin(x,xactive,t):
#    return force_interp(x,xactive) + force_random(x,xactive,t)
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
tfinal = 26 #38
Nt =  26*20 #5001   # Number of time steps
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
x = np.zeros((Nt,3))              # position vector
v = np.zeros((Nt,3))              # velocity vector
#v[:,0] = v[:,0] - 2e-6          # Drift adjust
w = np.random.normal(0,1,(Nt,3))  # Random white noise term
D = np.ones((Nt,3))*D0            # Diffustion tensor
# =============================================================================


# Dependent parameters 
# =============================================================================
delta_t = t[2]-t[1]     # time step
fct1 = np.sqrt(2*D*delta_t)
fct = np.multiply(fct1,w)
# =============================================================================



# Time evolution of the Langevin equation
# =============================================================================
x[0,1] = 100e-6

for m in range(Nt-1):
    x[m+1,:] = x[m,:] + fct[m,:] + delta_t*force_switch(x[m,:],t[m])/gamma + v[m,:]*delta_t
# =============================================================================    
    


x_um = x*1e6   #  um dimensions for plotting
#print(fct)
#print(x)


# Plotting time vs position data
py.figure(1)
#py.plot(t,x_um[:,0])
py.plot(t,x_um[:,1])
# py.plot(t,x_um[:,2])
py.xlabel('$t$ (s)')
py.ylabel('$y$ ($\mu$m)')


# Plotting position scatter data
py.figure(2)
py.plot(x_um[:,0],x_um[:,1],'r.')
py.xlabel('$x$ ($\mu$m)')
py.ylabel('$y$ ($\mu$m)')


#py.figure(3)



# Animation call
# =============================================================================


fig = py.figure()
fig.set_dpi(100)
fig.set_size_inches(7, 7)

ax = py.axes(xlim=(-1e6*xplt_limit, 1e6*xplt_limit), ylim=(-1e6*xplt_limit, 1e6*xplt_limit))
patch = py.Circle((0, 0), 10, fc='r')


time_template = 'Time = %.1f s'
time_string = ax.text(0.05, 0.9, '', transform=ax.transAxes)

anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=Nt, interval=100, blit=True)
anim.save('br_m1.mp4', fps=Nt/tfinal, extra_args=['-vcodec', 'libx264'])

# =============================================================================




