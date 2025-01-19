# Pressure driven internal pipe flow simulation
# Developing Hagen-Poiseuille prabola by solving incompressible Navier-Stokes equations

# LIBRARIES
import numpy as np
import matplotlib.pyplot as plt

# MAIN FUNCTION FOR THE SIMULATION
def main():
    # CONSTANTS
    N = 11                                      # number of points
    mu = 0.01                                   # kinemiatic viscosity
    dt = 0.2                                    # time step size
    Nt = 100                                    # total time steps

    Pressure_Gradient = np.array([-1.0, 0.0])   # pressure gradient

    element_length = 1.0/(N-1)                  # cell length
    x_range = np.linspace(0.0, 1.0, N)          # x distances
    y_range = np.linspace(0.0, 1.0, N)          # y distances
    cx, cy = np.meshgrid(x_range, y_range)      # 2D mesh coordinates

    # DISCRETIZED SPATIAL DERIVATIVES
    def central_difference(field):              # function for periodic central difference approximation
        diff = (
            (
            np.roll(field, shift=1, axis=1)     # step forward in x direction
            -                                   # subtraction
            np.roll(field, shift=-1, axis=1)    # step backward in x direction
            ) / (                               # division
                2 * element_length              # central difference approximation
            )
        )
        return diff
    
    def laplace(field):                         # function for periodic laplace operator
        diff = (
            (
            np.roll(field, shift=1, axis=1)     # step forward in x direction
            +                                   # addition
            np.roll(field, shift=1, axis=0)     # step forward in y direction
            +                                   # addition
            np.roll(field, shift=-1, axis=1)    # step backward in x direction
            +                                   # addition
            np.roll(field, shift=-1, axis=0)    # step backward in y direction
            -                                   # subtraction
            4 * field                           # field itself
            ) / (                               # division
                element_length**2               # five point stencil
            )
        )
        return diff
    
    # INITIAL CONDITIONS
    Vx_prev = np.ones((N, N))                   # initial velocity as 1 in domain
    Vx_prev[0, :] = 0.0                         # boundary condition at upper wall
    Vx_prev[-1, :] = 0.0                        # boundary condition at lower wall
    
    # TIME LOOP
    for iter in range(Nt):
        convection_x = (                        # calculating convection term
            Vx_prev                             # previous velocity
            *                                   # multiplication
            central_difference(Vx_prev)         # difference term
            )
        
        diffusion_x = (                         # calculating diffusion term
            mu                                  # kinematic viscosity
            *                                   # multiplication
            laplace(Vx_prev)                    # laplace term
            )

        Vx_next = (                             # calculating updated velocity
            Vx_prev                             # previous velocity
            +                                   # addition
            dt                                  # time step
            *                                   # multiplication
            (
                -                               # negative
                Pressure_Gradient[0]            # pressure gradient
                +                               # addition
                diffusion_x                     # diffusion
                -                               # subtraction
                convection_x                    # convection
            )
        )

        Vx_next[0, :] = 0.0                     # boundary condition at upper wall
        Vx_next[-1, :] = 0.0                    # boundary condition at lower wall

        Vx_prev = Vx_next                       # advancing in time
    
        # VISUALIZATION
        plt.contourf(cx, cy, Vx_next, levels=50, cmap='coolwarm')
        plt.colorbar()
        plt.quiver(cx, cy, Vx_next, np.zeros_like(Vx_next))
        plt.xlabel("Position along the Pipe")
        plt.ylabel("Cross section of the Pipe")

        plt.twiny()
        plt.plot(Vx_next[:, 1], cy[:, 1], color="white")
        plt.xlabel("Flow Velocity")

        plt.draw()
        plt.pause(0.05)
        plt.clf()   

if __name__ == "__main__":
    main()