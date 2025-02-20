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

        velocity_x_tent[1:-1, 0] = 1.0                              # left edge boundary condition
        velocity_x_tent[1:-1, -1] = velocity_x_tent[1:-1, -2]       # right edge boundary condition
        velocity_x_tent[0, :] = -velocity_x_tent[1, :]              # bottom edge boundary condition
        velocity_x_tent[-1, :] = -velocity_x_tent[-2, :]            # top edge boundary condition

        # UPDATING INTERIOR X VELOCITY WITH MOMENTUM EQUATION
        diffusion_y = MU * (                                        # definition of diffusion equation
            (
                velocity_y_prev[1:-1, 2: ]                          # forward stencil point x direction
                +
                velocity_y_prev[2: , 1:-1]                          # backward stencil point x direction
                +
                velocity_y_prev[1:-1, :-2]                          # forward stencil point y direction
                +
                velocity_y_prev[ :-2, 1:-1]                         # backward stencil point y direction
                - 4 *
                velocity_y_prev[1:-1, 1:-1]                         # interior velocity field
            ) / (
                cell_lenght**2                                      # cell lenght squared
            )
        )

        convection_y = (                                            # definition of convection equation
            (
                velocity_x_prev[2:-1, 1: ]                          # prefactor forward based on x velocity
                +
                velocity_x_prev[2:-1, :-1]                          # prefactor backward based on x velocity
                +
                velocity_x_prev[1:-2, 1: ]                          # prefactor forward based on x velocity
                +
                velocity_x_prev[1:-2, :-1]                          # prefactor backward based on x velocity
            ) / 4
            *
            (
                velocity_y_prev[1:-1, 2: ]                          # forward stencil point x direction
                -
                velocity_y_prev[1:-1, :-2]                          # backward stencil point x direction
            ) / (
                2 * cell_lenght                                     # cell lenght
            )
            +
            (
                velocity_y_prev[2: ,1:-1]**2                        # top left stencil point y velocity
                -
                velocity_y_prev[ :-2 ,1:-1]**2                      # top right stencil point y velocity
            ) / (
                2* cell_lenght                                      # cell lenght
            )
        )

        pressure_gradient_y = (                                     # definition of pressure gradient
            (
                pressure_prev[2:-1, 1:-1]                           # forward stencil interior pressure in x
                -
                pressure_prev[1:-2, 1:-1]                           # backward stencil interior pressure in x
            ) / (
                cell_lenght                                         # cell lenght
            )
        )

        velocity_y_tent[1:-1, 1:-1] = (                             # interior tentative velocity in x direction
            velocity_y_prev[1:-1, 1:-1]                             # previous tentative velocity in x direction
            +
            TIME_STEP                                               # per time step
            *
            (
                -pressure_gradient_y                                # negative pressure gradient
                +
                diffusion_y                                         # effect of diffusion
                -
                convection_y                                        # effect of convection
            )
        )

        velocity_y_tent[1:-1, 0] = -velocity_y_tent[1:-1, -1]       # left edge boundary condition
        velocity_y_tent[1:-1, -1] = velocity_y_tent[1:-1, -2]       # right edge boundary condition
        velocity_y_tent[0, :] = 0.0                                 # bottom edge boundary condition
        velocity_y_tent[-1, :] = 0.0                                # top edge boundary condition

        # COMPUTING DIVERGENCE FOR PRESSURE POISSON PROBLEM
        divergence = (                                                      # definition of divergence
            (
                velocity_x_tent[1:-1, 1:]                                   # tentative velocity interior forward in x direction
                -
                velocity_x_tent[1:-1, :-1]                                  # tentative velocity interior backward in x direction
            ) / (
                cell_lenght                                                 # cell lenght
            )
            +
            (
                velocity_y_tent[1: ,1:-1]                                   # tentative velocity interior forward in y direction
                -
                velocity_y_tent[ :-1, 1:-1]                                 # tentative velocity interior backward in x direction
            ) / (
                cell_lenght                                                 # cell lenght
            )
        )

        density = 1                                                         # fluid density
        pressure_poisson_rhs = divergence * density  / TIME_STEP            # pressure poisson equation right hand side

        # SOLVING PRESSURE CORRECTION POISSON EQUATION
        pressure_corr_prev = np.zeros_like(pressure_prev)                   # previous pressure correction

        for _ in range(N_POISSON):                                          # iteration for specific poisson iteration count
            pressure_corr_next = np.zeros_like(pressure_corr_prev)          # next pressure correction
            pressure_corr_next[1:-1, 1:-1] = 1/4 * (
                pressure_corr_prev[1:-1, 2: ]                               # interior pressure forward in x direction
                +
                pressure_corr_prev[2: , 1:-1]                               # interior pressure forward in y direction
                +
                pressure_corr_prev[1:-1, :-2]                               # interior pressure backward in x direction
                +
                pressure_corr_prev[ :-2, 1:-1]                              # interior pressure backward in y direction
                -
                cell_lenght**2                                              # cell lenght squared
                *
                pressure_poisson_rhs                                        # pressure poisson right hand side
            )

            pressure_corr_next[1:-1, 0] = pressure_corr_next[1:-1, 1]       # homogenous pressure left boundary condition
            pressure_corr_next[1:-1, -1] = -pressure_corr_next[1:-1, -2]    # homogenous pressure right boundary condition
            pressure_corr_next[0, :] = pressure_corr_next[1, :]             # homogenous pressure top boundary condition
            pressure_corr_next[-1, :] = pressure_corr_next[-2, :]           # homogenous pressure bottom boundary condition

            pressure_corr_prev = pressure_corr_next                         # advance in smoothing to enforce incompressibility

        pressure_next = pressure_prev + pressure_corr_next                  # updatingthe pressure
        
        # UPDATE VELOCITY
        pressure_corr_grad_x = (                                            # pressure correction gradient in x
            (
                pressure_corr_next[1:-1, 2:-1]                              # interior pressure forward in x direction
                -
                pressure_corr_next[1:-1, 1:-2]                              # interior pressure backward in x direction
            ) / (
                cell_lenght                                                 # cell lenght
            )
        )

        velocity_x_next[1:-1, 1:-1] = (                                     # velocity ıpdate in x direction
            velocity_x_tent[1:-1, 1:-1]                                     # tentative velocity in x
            -
            TIME_STEP                                                       # time step lenght
            *
            pressure_corr_grad_x                                            # pressure correction gradient in x direction
        )

        pressure_corr_grad_y = (                                            # pressure correction gradient in y
            (
                pressure_corr_next[2:-1, 1:-1]                              # interior pressure forward in y direction
                -
                pressure_corr_next[1:-2, 1:-1]                              # interior pressure backward in y direction
            ) / (
                cell_lenght                                                 # cell lenght
            )
        )

        velocity_y_next[1:-1, 1:-1] = (                                     # velocity ıpdate in y direction
            velocity_y_tent[1:-1, 1:-1]                                     # tentative velocity in y
            -
            TIME_STEP                                                       # time step lenght
            *
            pressure_corr_grad_y                                            # pressure correction gradient in y direction
        )
        
        velocity_x_next[1:-1, 0] = 1.0                              # left edge velocity boundary condition
        inflow_mass_rate_next = np.sum(velocity_x_next[1:-1, 0])    # inlet total velocity
        outflow_mass_rate_next = np.sum(velocity_x_next[1:-1, -2])  # outlet total velocity
        mrr = inflow_mass_rate_next / outflow_mass_rate_next        # mass rate ratio
        velocity_x_next[1:-1, -1] = velocity_x_next[1:-1, -2] * mrr # right edge velocity boundary condition
        velocity_x_next[0, :] = -velocity_x_next[1, :]              # bottom edge velocity boundary condition
        velocity_x_next[-1, :] = -velocity_x_next[-2, :]            # top edge velocity boundary condition

        velocity_y_next[1:-1, 0] = -velocity_y_next[1:-1, -1]       # left edge velocity boundary condition
        velocity_y_next[1:-1, -1] = velocity_y_next[1:-1, -2]       # right edge velocity boundary condition
        velocity_y_next[0, :] = 0.0                                 # bottom edge velocity boundary condition
        velocity_y_next[-1, :] = 0.0                                # top edge velocity boundary condition

        # REPEATING THE ITERATIONS
        velocity_x_prev = velocity_x_next                           # advance in time for x velocity
        velocity_y_prev = velocity_y_next                           # advance in time for y velocity
        pressure_prev = pressure_next                               # advance in time for pressure

        # VISUALIZATION
        if iter % PLOT_EVERY == 0:
            velocity_x_vertex_centered = (                          # vetices to be plotted in x
                (
                    velocity_x_next[1: , :]                         # left values
                    +
                    velocity_x_next[ :-1, :]                        # right values
                ) / 2                                               # calculating the mean value
            )

            velocity_y_vertex_centered = (                          # vetices to be plotted in y
                (
                    velocity_y_next[:, 1:]                          # top values
                    +
                    velocity_y_next[:, :-1]                         # bottom values
                ) / 2                                               # calculating the mean value
            )

            plt.contourf(                                           # plotting thr velocity field
                coordinates_x,
                coordinates_y,
                velocity_x_vertex_centered,
                levels=10,
                cmap='coolwarm'
            )
            plt.colorbar()

            plt.quiver(                                             # plotting the x velocity components
                coordinates_x[:, ::6],
                coordinates_y[:, ::6],
                velocity_x_vertex_centered[:, ::6],
                velocity_y_vertex_centered[:, ::6],
                alpha=0.4,
            )

            plt.plot(                                               # velocity parabola
                5*cell_lenght + velocity_x_vertex_centered[: ,5],
                coordinates_y[:, 5],
                color='white'
            )

            plt.plot(                                               # velocity parabola
                40*cell_lenght + velocity_x_vertex_centered[: ,40],
                coordinates_y[:, 40],
                color='white'
            )

            plt.plot(                                               # velocity parabola
                80*cell_lenght + velocity_x_vertex_centered[: ,80],
                coordinates_y[:, 80],
                color='white'
            )

            plt.title(f'Iteration {iter}/{N_TIME_STEPS}')
            plt.xlabel("Position along the Pipe")
            plt.ylabel("Cross section of the Pipe")

            plt.draw()
            plt.pause(0.05)
            plt.clf()

if __name__ == "__main__":
    main()
