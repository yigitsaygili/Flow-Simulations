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

![periodic_pipe](https://github.com/user-attachments/assets/d836f719-209c-4a04-81ee-80877f40d8b6)

Developed with the help of:

https://www.youtube.com/watch?v=PtfAPqIDaMI


## Airfoil Generation with Joukowsky Transformation
Joukowsky airfoils are a kind of airfoils around which the transform is employed to solve for the two-dimensional potential flow in aerodynamics. Applying the Joukowsky transform to a circle creates a Joukowsky airfoil in the complex plane. The shape of the resulting airfoil can be altered by changing the coordinates of the circle's center.

By applying the Joukowski mapping to a circle in the Argand diagram, one may easily represent the cross section of an airfoil. Except for the locations where the complex derivative is 0, the map is conformal. A number of airfoil forms can be created by dragging the center of the circle, but it must pass through one of these points and either enclose or pass through the other.

Pressure Disrtibution:

![juk_air](https://github.com/user-attachments/assets/835d7728-ebc8-42be-9fa1-1b321ce75073)

Developed with the help of:

VBBV (2025). Matlab program for Joukowski Airfoil (https://www.mathworks.com/matlabcentral/fileexchange/65430-matlab-program-for-joukowski-airfoil), MATLAB Central File Exchange. Retrieved February 9, 2025.

