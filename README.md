# Lattice-Boltzmann-Fluid-Simulation
A 2D Fluid Simulation Using Lattice Boltzmann Method

A script designed to simulate fluid flow around a cylinder utilizing the Lattice Boltzmann method, an efficient and straightforward approach for fluid dynamics simulation. Rather than directly solving the Navier-Stokes equations, this method models the behavior of microscopic particles on a lattice through streaming and collision processes.

The fluid is represented by a set of particles that follow specific directions and update their positions at each time step. When the particles collide, they adjust their velocities based on certain rules, which helps simulate the behavior of the fluid. This method is particularly useful for simulating complex fluid behavior, such as flows with intricate boundaries or multiple interacting phases, and it is computationally efficient, making it suitable for problems where traditional methods might be too slow or difficult to apply.

The Velocity Gradient:

![velo_grad](https://github.com/user-attachments/assets/08b555a1-5dd4-4c61-b510-246fba2877b9)

The Curl Gradient:

![velo_grad](https://github.com/user-attachments/assets/dcce8200-1854-48b5-b20a-020ddd4671a5)


Developed with the help of:
https://medium.com/swlh/create-your-own-lattice-boltzmann-simulation-with-python-8759e8b53b1c
https://www.youtube.com/watch?v=JFWqCQHg-Hs
