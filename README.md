# 2D Fluid Simulatons in Python

## Lattice-Boltzmann-Fluid-Simulation
Fluid flow simulation around a cyclinder using Lattice Boltzmann method. A script designed to simulate fluid flow around a cylinder utilizing the Lattice Boltzmann method, an efficient and straightforward approach for fluid dynamics simulation. Rather than directly solving the Navier-Stokes equations, this method models the behavior of microscopic particles on a lattice through streaming and collision processes.

The fluid is represented by a set of particles that follow specific directions and update their positions at each time step. When the particles collide, they adjust their velocities based on certain rules, which helps simulate the behavior of the fluid. This method is particularly useful for simulating complex fluid behavior, such as flows with intricate boundaries or multiple interacting phases, and it is computationally efficient, making it suitable for problems where traditional methods might be too slow or difficult to apply.

The Velocity Gradient:

![velo_grad](https://github.com/user-attachments/assets/08b555a1-5dd4-4c61-b510-246fba2877b9)

The Curl Gradient:

![curl_grad](https://github.com/user-attachments/assets/0b0f5a18-fdc0-43ef-9b31-bcadb710dc28)

The cylinder can be replaced with any other arbitrary shaped obstacle as long as it is defined in grid format.

Developed with the help of:

https://medium.com/swlh/create-your-own-lattice-boltzmann-simulation-with-python-8759e8b53b1c

https://www.youtube.com/watch?v=JFWqCQHg-Hs


## Pressure Driven Internal Pipe Flow
Developing Hagen-Poiseuille prabola by solving incompressible Navier-Stokes equations for an internal pipe flow. The pipe is designed to be periodic to simulate an infititely long pipe and to demonstrate stability condition.

The flow is forming boundary layers close to solid wall boundary conditions as a result of the viscous effects of the Navier-Stokes equations. In this video, we examine a pipe flow where the domain's top end and bottom edge are wall boundary conditions. The initially uniform velocity profile is forced to take on a parabolic shape by this symmetry. One of these uncommon situations is the pipe flow caused by a pressure gradient, for which the law of Hagen Poiseuille, an analytic solution to the Navier-Stokes equations, is actually available. 

Velocity Gradient and Boundary Layer Formation:

![periodic_pipe](https://github.com/user-attachments/assets/911b8db7-1c7a-420a-b169-22683728a24b)

Developed with the help of:

https://www.youtube.com/watch?v=PtfAPqIDaMI

