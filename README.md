# Brownian-dynamics-in-a-time-varying-force-field

<p float="left">
<a href = "https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/tree/main/Codes"> <img src="https://img.shields.io/badge/Language-Python-blue" alt="alt text"> </a>
<a href = "https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/LICENSE"> <img src="https://img.shields.io/badge/license-MIT-green" alt="alt text"></a>
<a href = "https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/tree/main/Codes"> <img src="https://img.shields.io/badge/version-1.0-red" alt="alt text"> </a>
</p>

<img align = "right" src="https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/blob/main/Brownina_moving_force.gif" alt="alt text" width="360">

A python code to calculate the Brownian motion of a colloidal particle in a time varying force field. The code was developed for the research work of the <a href = "https://hesselink-lab.stanford.edu"> Hesselink research group </a>  at Stanford University. The particle trajectory computed by the code can be useful for designing lab-on-a-chip devices. Please cite this repository when using this code.

The current example code shows the transport of a micro-particle along a moving micro-electrode array producing moving dielectrophoretic force. As the micro-electrodes are excited in sequence, the micro-particle follows the position of the active electrode.

The Brownian motion code solves the Langevin equation in discrete time. The code is general. Although the example considers a dielectrophoretic force-field, the code would work for any type of external force-field. The code would also work for nano-particles instead of micro-particles. The Brownian vibrations are more prominent for nano-particles.  

The python script was tested with spyder IDE. 

## Theory
The Brownian motion of a colloidal particle in a low Reynolds number environment can be modeled by the Langevin equation: 

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{\dot{r}}(t) = \frac{\overset{\leftrightarrow}{\mathbf{D}}(\mathbf{r})}{k_B T} \mathbf{F}_\text{ext}(\mathbf{r},t) %2B \sqrt{2} \overset{\leftrightarrow}{\mathbf{D}}_{1/2}(\mathbf{r})  \mathbf{W}(t)">

Here,  <img src="https://render.githubusercontent.com/render/math?math=\mathbf{r} = (x_o,y_o,z_o)"> is the position of the center of the particle,  <img src="https://render.githubusercontent.com/render/math?math=\mathbf{F}_\text{ext}"> is the external trapping/manipulation force acting on the particle,  <img src="https://render.githubusercontent.com/render/math?math=k_B"> is the Boltzmann constant,  <img src="https://render.githubusercontent.com/render/math?math=T"> is the temperature,  <img src="https://render.githubusercontent.com/render/math?math=\overset{\leftrightarrow}{\mathbf{D}}"> is the diffusion tensor, and  <img src="https://render.githubusercontent.com/render/math?math=\mathbf{W}(t)"> is a vector white noise term. The tensor  <img src="https://render.githubusercontent.com/render/math?math=\overset{\leftrightarrow}{\mathbf{D}}_{1/2}"> is defined as the element-wise square root of  <img src="https://render.githubusercontent.com/render/math?math=\overset{\leftrightarrow}{\mathbf{D}}">. Each Cartesian component of  <img src="https://render.githubusercontent.com/render/math?math=\mathbf{W}(t)"> is a random process unit zero mean and unit variance. 

The Langevin equation can be converted into the following discrete form:
 <img src="https://render.githubusercontent.com/render/math?math=\mathbf{r}_{i%2B1} = \mathbf{r}_i %2B \frac{\overset{\leftrightarrow}{\mathbf{D}}(\mathbf{r}_i)}{k_B T} \mathbf{F}_\text{ext}(\mathbf{r}_i,t_i) %2B \sqrt{2 \Delta t} \overset{\leftrightarrow}{\mathbf{D}}_{1/2}(\mathbf{r}_i)  \mathbf{W}_i">

The equaiton is numerically solved using the Euler-Maruyama method.


## Requirements






## Sample output


## References

