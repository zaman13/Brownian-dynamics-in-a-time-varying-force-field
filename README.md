# Brownian-dynamics-in-a-time-varying-force-field

<p float="left">
<a href = "https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/tree/main/Codes"> <img src="https://img.shields.io/badge/Language-Python-blue" alt="alt text"> </a>
<a href = "https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/LICENSE"> <img src="https://img.shields.io/badge/license-MIT-green" alt="alt text"></a>
<a href = "https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/tree/main/Codes"> <img src="https://img.shields.io/badge/version-1.9-red" alt="alt text"> </a>
</p>

<img align = "right" src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Figs_sorter/cell_sorter.gif" alt="alt text" width="480">

A python code to calculate the Brownian motion of a colloidal particle in a time varying force field. The code was developed for the research work of the <a href = "https://hesselink-lab.stanford.edu"> Hesselink research group </a>  at Stanford University. The particle trajectory computed by the code can be useful for designing lab-on-a-chip devices. Please cite this repository and the two papers listed in the reference section when using this code.

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

<img src="https://render.githubusercontent.com/render/math?math=m\frac{\partial\mathbf{v}(\mathbf{r},t)}{\partial t} = \frac{k_B T}{\overset{\leftrightarrow}{\mathbf{D}}(\mathbf{r})} \Bigg[\mathbf{v}_f(\mathbf{r},t)-\mathbf{v}(\mathbf{r},t)  %2B \sqrt{2} \overset{\leftrightarrow}{\mathbf{D}}_{\frac{1}{2}}(\mathbf{r})  \mathbf{W}(t) \Bigg] \\ %2B  \mathbf{F}_\text{ext}(\mathbf{r},t) \,.">

Here,  <img src="https://render.githubusercontent.com/render/math?math=\mathbf{r} = (x_o,y_o,z_o)"> is the position of the center of the particle, <img src="https://render.githubusercontent.com/render/math?math=\mathbf{v}"> is the particle velocity, <img src="https://render.githubusercontent.com/render/math?math=\mathbf{v}_f"> is the fluid velocity, <img src="https://render.githubusercontent.com/render/math?math=\mathbf{F}_\text{ext}"> is the external trapping/manipulation force acting on the particle,  <img src="https://render.githubusercontent.com/render/math?math=k_B"> is the Boltzmann constant,  <img src="https://render.githubusercontent.com/render/math?math=T"> is the temperature,  <img src="https://render.githubusercontent.com/render/math?math=\overset{\leftrightarrow}{\mathbf{D}}"> is the diffusion tensor, and  <img src="https://render.githubusercontent.com/render/math?math=\mathbf{W}(t)"> is a vector white noise term. The tensor  <img src="https://render.githubusercontent.com/render/math?math=\overset{\leftrightarrow}{\mathbf{D}}_{1/2}"> is defined as the element-wise square root of  <img src="https://render.githubusercontent.com/render/math?math=\overset{\leftrightarrow}{\mathbf{D}}">. Each Cartesian component of  <img src="https://render.githubusercontent.com/render/math?math=\mathbf{W}(t)"> is a random process unit zero mean and unit variance. 

The Langevin equation can be converted into the following discrete form:

 <img src="https://render.githubusercontent.com/render/math?math=m\frac{\mathbf{v}_{i%2B1} - \mathbf{v}_i}{\Delta t} = \frac{k_B T}{\overset{\leftrightarrow}{\mathbf{D}}(\mathbf{r}_i)} \left[\mathbf{v}_{f,i}-\mathbf{v}_{i%2B1} %2B \sqrt{\frac{2}{\Delta t}} \overset{\leftrightarrow}{\mathbf{D}}_{\frac{1}{2}}(\mathbf{r}_i)  \mathbf{w}_i \right]  %2B  \mathbf{F}_{\text{ext},i}(\mathbf{r}_i) \,.">

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

A large optical spot and a strong gradient force model is used here for illustration purposes. It is possible to simulate optical tweezers with arbitrary spot size and traping potential depth.

### Sample output for dielectrophoretic cell sorter/separator:
<p>
A microfluidic device that can sort/separate live and dead yeast cells is simulated. The device uses electrodes to apply dielectrophoretic forces. Since the material properties and hence the Clausius-Mossotti factor of live cells (read circles) differ from that of dead cells (blue circles), they experience different dielectrophoretic forces. Using this principle, the cell sorting/separation operation is accomplished.  
</p>

<p float="left">
<img align = "left" src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Figs_sorter/cell_sorter.gif" alt="alt text" width="720">
</p>



<p float="left">
<img src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Figs_sorter/y.svg" alt="alt text" width="480">
</p>

<p>
The cell separation mechanism can be analyzed by observing the y trajectory of different cells.
</p>

## Acknowledgement
This work is partially supported by the National Institute of Health (NIH) Grant R01GM138716. 


## References
1. Zaman, Mohammad Asif, et al. "Modeling Brownian Microparticle Trajectories in Lab-on-a-Chip Devices with Time Varying Dielectrophoretic or Optical Forces." Micromachines 12.10 (2021): 1265.
<a href = "https://doi.org/10.1063/5.0049126"> https://doi.org/10.1063/5.0049126 </a>
2. Zaman, Mohammad Asif, et al. "Microparticle transport along a planar electrode array using moving dielectrophoresis." Journal of Applied Physics 130.3 (2021): 034902.
<a href = "https://doi.org/10.3390/mi12101265"> https://doi.org/10.3390/mi12101265 </a>
