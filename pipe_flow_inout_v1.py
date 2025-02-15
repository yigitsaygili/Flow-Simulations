# Transient internal pipe flow simulation
# Includes basic SIMPLE algorithm

# LIBRARIES
import numpy as np
import matplotlib.pyplot as plt

# CONSTANTS
N_POINTS_Y = 15                                                     # number of points in y direction
AR = 10                                                             # aspect ratio of the pipe
MU = 0.01                                                           # kinematic viscosity
TIME_STEP = 0.001                                                   # time step length
N_TIME_STEPS = 5000                                                 # total time steps
PLOT_EVERY = 50                                                     # plot frames per time step
N_POISSON = 50                                                      # number of pressure poisson iterations

# MAIN FUNCTION FOR THE SIMULATION
def main():
    cell_lenght = 1.0 / (N_POINTS_Y-1)                              # cell lenght
    n_points_x = (N_POINTS_Y-1) * AR + 1                            # number of points in x direction
    x_range = np.linspace(0.0, 1.0*AR, n_points_x)                  # x coordinates
    y_range = np.linspace(0.0, 1.0, N_POINTS_Y)                     # x coordinates
    coordinates_x, coordinates_y = np.meshgrid(x_range, y_range)    # 2D mesh domain

    # INITIAL CONDITIONS
    velocity_x_prev = np.ones((N_POINTS_Y+1, n_points_x))           # initial velocity in x direction
    velocity_x_prev[0, :] = -velocity_x_prev[1, :]                  # upper wall boundary condition
    velocity_x_prev[-1, :] = -velocity_x_prev[-2, :]                # lower wall boundary condition

    velocity_y_prev = np.ones((N_POINTS_Y, n_points_x+1))           # initial velocity in y direction
    pressure_prev = np.zeros((N_POINTS_Y+1, n_points_x+1))          # initial uniform zero pressure

    # PRE-ALLOCATING ARRAYS
    velocity_x_tent = np.zeros_like(velocity_x_prev)                # pre-allocated tentative x velocity
    velocity_x_next = np.zeros_like(velocity_x_prev)                # pre-allocated next x velocity
    velocity_y_tent = np.zeros_like(velocity_y_prev)                # pre-allocated tentative y velocity
    velocity_y_next = np.zeros_like(velocity_y_prev)                # pre-allocated next y velocity

    plt.figure(figsize=(1.5*AR, 4))

    # MAIN TIME LOOP
    for iter in range(N_TIME_STEPS):
        # UPDATING INTERIOR X VELOCITY WITH MOMENTUM EQUATION
        diffusion_x = MU * (                                        # definition of diffusion equation
            (
                velocity_x_prev[1:-1, 2: ]                          # forward stencil point x direction
                +
                velocity_x_prev[2: , 1:-1]                          # backward stencil point x direction
                +
                velocity_x_prev[1:-1, :-2]                          # forward stencil point y direction
                +
                velocity_x_prev[ :-2, 1:-1]                         # backward stencil point y direction
                - 4 *
                velocity_x_prev[1:-1, 1:-1]                         # interior velocity field
            ) / (
                cell_lenght**2                                      # cell lenght squared
            )
        )

        convection_x = (                                            # definition of convection equation
            (
                velocity_x_prev[1:-1, 2: ]**2                       # forward stencil point x direction
                -
                velocity_x_prev[1:-1, :-2]**2                       # backward stencil point x direction
            ) / (
                2 * cell_lenght
            )
            +
            (
                velocity_y_prev[1: ,1:-2]                           # top left stencil point y velocity
                +
                velocity_y_prev[1: ,2:-1]                           # top right stencil point y velocity
                +
                velocity_y_prev[ :-1, 1:-2]                         # bottom left stencil point y velocity
                +
                velocity_y_prev[ :-1, 2:-1]                         # bottom right stencil point y velocity
            ) / 4
            *
            (
                velocity_x_prev[2: , 1:-1]                          # forward stencil point y direction
                -
                velocity_x_prev[ :-2, 1:-1]                         # backward stencil point y direction
            ) / (
                2 * cell_lenght                                     # cell lenght squared
            )
        )

        pressure_gradient_x = (                                     # definition of pressure gradient
            (
                pressure_prev[1:-1, 2:-1]                           # forward stencil interior pressure in x
                -
                pressure_prev[1:-1, 1:-2]                           # backward stencil interior pressure in x
            ) / (
                cell_lenght                                         # cell lenght
            )
        )
    
        velocity_x_tent[1:-1, 1:-1] = (                             # interior tentative velocity in x direction
            velocity_x_prev[1:-1, 1:-1]                             # previous tentative velocity in x direction
            +
            TIME_STEP                                               # per time step
            *
            (
                -pressure_gradient_x                                # negative pressure gradient
                +
                diffusion_x                                         # effect of diffusion
                -
                convection_x                                        # effect of convection
            )
        )

if __name__ == "__main__":
    main()