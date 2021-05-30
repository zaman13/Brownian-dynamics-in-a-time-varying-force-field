# Brownian-dynamics-in-a-time-varying-force-field

<p float="left">
<a href = "https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/tree/main/Codes"> <img src="https://img.shields.io/badge/Language-Python-blue" alt="alt text"> </a>
<a href = "https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/LICENSE"> <img src="https://img.shields.io/badge/license-MIT-green" alt="alt text"></a>
<a href = "https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/tree/main/Codes"> <img src="https://img.shields.io/badge/version-1.8-red" alt="alt text"> </a>
</p>

<img align = "right" src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Brownina_moving_force.gif" alt="alt text" width="360">

A python code to calculate the Brownian motion of a colloidal particle in a time varying force field. The code was developed for the research work of the <a href = "https://hesselink-lab.stanford.edu"> Hesselink research group </a>  at Stanford University. The particle trajectory computed by the code can be useful for designing lab-on-a-chip devices. Please cite this repository when using this code.

The current example code shows the transport of a micro-particle along a moving micro-electrode array producing moving dielectrophoretic force. As the micro-electrodes are excited in sequence, the micro-particle follows the position of the active electrode.

The Brownian motion code solves the Langevin equation in discrete time. The code is general. Although the example considers a dielectrophoretic force-field, the code would work for any type of external force-field. The code would also work for nano-particles instead of micro-particles. The Brownian vibrations are more prominent for nano-particles.  

The python script was tested with spyder IDE. 

## Package requirements
  - <b>NumPy</b>: <a href = "https://docs.scipy.org/doc/numpy/reference/"> <img src="https://img.shields.io/badge/Pkg-NumPy-FF4500" alt="alt text"> </a>
  - <b>Pylab, Matplotlib</b>: (sparse matrices, sparse linear algebra) <a href = "https://scipy.github.io/old-wiki/pages/PyLab"> <img src="https://img.shields.io/badge/Pkg-Pylab-FF7F50" alt="alt text"> </a> <a href = "https://www.tutorialspoint.com/matplotlib/matplotlib_pylab_module.htm"> <img src="https://img.shields.io/badge/Pkg-Matplotlib-FF7F50" alt="alt text"> </a>
  - <b>FuncAnimation()</b>: for animation <a href = "https://matplotlib.org/stable/api/_as_gen/matplotlib.animation.FuncAnimation.html"> <img src="https://img.shields.io/badge/Function-FuncAnimation-00783f" alt="alt text"> </a>
- <b>Ffmpeg</b>: for saving animation to video file <a href = "https://anaconda.org/conda-forge/ffmpeg"> <img src="https://img.shields.io/badge/Pkg-Ffmpeg-00783f" alt="alt text"> </a> 


## Incorporated Physics
- Brownian motion
- Hindered diffusion (hydrodynamic interactions)
- Elastic collisions (particle-particle collisions, particle-wall collisions)

## Features
- Can simulate multi-particle systems
- Each particle can have its own distinct size/mass
- Can save the animation of the motion of the particles as a video file
- Can create arbitrary polygon based geometry of walls (the collision mechanics are automatically formulated for any defined wall) 

## Theory
The Brownian motion of a colloidal particle in a low Reynolds number environment can be modeled by the Langevin equation: 

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{\dot{r}}(t) = \frac{\overset{\leftrightarrow}{\mathbf{D}}(\mathbf{r})}{k_B T} \mathbf{F}_\text{ext}(\mathbf{r},t) %2B \sqrt{2} \overset{\leftrightarrow}{\mathbf{D}}_{1/2}(\mathbf{r})  \mathbf{W}(t)">

Here,  <img src="https://render.githubusercontent.com/render/math?math=\mathbf{r} = (x_o,y_o,z_o)"> is the position of the center of the particle,  <img src="https://render.githubusercontent.com/render/math?math=\mathbf{F}_\text{ext}"> is the external trapping/manipulation force acting on the particle,  <img src="https://render.githubusercontent.com/render/math?math=k_B"> is the Boltzmann constant,  <img src="https://render.githubusercontent.com/render/math?math=T"> is the temperature,  <img src="https://render.githubusercontent.com/render/math?math=\overset{\leftrightarrow}{\mathbf{D}}"> is the diffusion tensor, and  <img src="https://render.githubusercontent.com/render/math?math=\mathbf{W}(t)"> is a vector white noise term. The tensor  <img src="https://render.githubusercontent.com/render/math?math=\overset{\leftrightarrow}{\mathbf{D}}_{1/2}"> is defined as the element-wise square root of  <img src="https://render.githubusercontent.com/render/math?math=\overset{\leftrightarrow}{\mathbf{D}}">. Each Cartesian component of  <img src="https://render.githubusercontent.com/render/math?math=\mathbf{W}(t)"> is a random process unit zero mean and unit variance. 

The Langevin equation can be converted into the following discrete form:
 <img src="https://render.githubusercontent.com/render/math?math=\mathbf{r}_{i%2B1} = \mathbf{r}_i %2B \frac{\overset{\leftrightarrow}{\mathbf{D}}(\mathbf{r}_i)}{k_B T} \mathbf{F}_\text{ext}(\mathbf{r}_i,t_i) %2B \sqrt{2 \Delta t} \overset{\leftrightarrow}{\mathbf{D}}_{1/2}(\mathbf{r}_i)  \mathbf{W}_i">

The equaiton is numerically solved using the Euler-Maruyama method.








## Sample output
### Sample output for moving DEP force profile:
The blue horizonal lines represent electrodes which are excited in sequence. Only one electrode is active at a time (the rest are grounded). A solid surface is assumed to be located along z = 0 plane. The red spheres are polystyrene micro-particles, each having a radius of 10 micron. 
<p float="left">
<img src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Figs_moving_DEP/br_v1.6.gif" alt="alt text" width="720">
</p>

The position and velocity profiles are:
<p float="left">
<img src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Figs_moving_DEP/Fig_DEP_x.png" alt="alt text" width="240">
<img src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Figs_moving_DEP/Fig_DEP_y.png" alt="alt text" width="240">
</p>

<p float="left">
<img src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Figs_moving_DEP/Fig_DEP_z.png" alt="alt text" width="240">
<img src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Figs_moving_DEP/Fig_DEP_vy.png" alt="alt text" width="240">
</p>


### Sample output for optical trap (Gaussian potential well):
A optical trap that is active during the time interval 1s < t < 8s is simulated. The optical spot is assumed to create a Gaussian potential well. The resulting force profile is given by the gradient of the potential well. A solid surface is assumed to be located along z = 0 plane (representing the substrate). Particles hitting the solid surface experience elastic collision. 

<p float="left">
<img src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Figs_Gaussian/br_v1.6.gif" alt="alt text" width="720">
</p>

It can be noted that particles outside capturing range of the trap (i.e. beyond the range where the trapping force is significant) do not get trapped.

The position and velocity profiles are:
<p float="left">
<img src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Figs_Gaussian/x.png" alt="alt text" width="240">
<img src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Figs_Gaussian/y.png" alt="alt text" width="240">
</p>

<p float="left">
<img src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Figs_Gaussian/z.png" alt="alt text" width="240">
<img src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Figs_Gaussian/vy.png" alt="alt text" width="240">
</p>




## References

