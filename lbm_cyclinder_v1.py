# Lattice-Boltzmann Fluid Simulation

# LIBRARIES
import numpy as np
from matplotlib import pyplot 

# FUNCTION FOR DISTANCE CALCULATION
def distance(x1, y1, x2, y2):
    return np.sqrt((x2-x1)**2 + (y2-y1)**2)

# MAIN FUNCTION FOR THE SIMULATION
def main():
    # CONSTANTS
    Nx = 400                                                                # number of lattices in x direction
    Ny = 100                                                                # number of lattices in y direction
    tau = 0.53                                                              # kinematic viscoty
    Nt = 3000                                                               # number of iterations

    # LATTICE SPEEDS AND WEIGHTS
    NL = 9                                                                  # number of discrete velocities
    cxs = np.array([0, 0, 1, 1, 1, 0,-1,-1,-1])                             # discrete velocity locations in x
    cys = np.array([0, 1, 1, 0,-1,-1,-1, 0, 1])                             # discrete velocity locations in y
    weights = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])   # discrete velocity weights

    # INITIAL CONDITIONS
    F = np.ones((Ny, Nx, NL)) + 0.1*np.random.randn(Ny, Nx, NL)             # mesoscopic velocity
    F[:, :, 3] = 2.3                                                        # constant right hand side velocity

    # BOUNDARIES
    cyclinder = np.full((Ny, Nx), False)                                    # definition of non boundary cells

    for y in range(0, Ny):                                                  # loop all cells in y direction
        for x in range(0, Nx):                                              # loop all cells in x direction
            if (distance(Nx//4, Ny//2, x, y)<13):                           # measure distance
                cyclinder[y][x] = True                                      # definition of boundary cells

    # MAIN LOOP
    for it in range(Nt):                                                    # time interval loop
        # print(it)                                                           # print iterations for debugging

        # ABSORBING BOUNDARIES
        F[:, -1, [6, 7, 8]] = F[:, -2, [6, 7, 8]]                           # right end velocity absorbtion
        F[:, 0, [2, 3, 4]] = F[:, 1, [2, 3, 4]]                             # right end velocity absorbtion

        # STREAMING
        for i, cx, cy in zip(range(NL), cxs, cys):                          # streaming loop into neighbour cells
            F[:, :, i] = np.roll(F[:, :, i], cx, axis=1)                    # rolling velocity in x direction
            F[:, :, i] = np.roll(F[:, :, i], cy, axis=0)                    # rolling velocity in y direction

        # BOUNDARY
        bndryF = F[cyclinder, :]                                            # checking the collision
        bndryF = bndryF[:, [0, 5, 6, 7, 8, 1, 2, 3, 4]]                     # inversing the velocity

        # FLUID VARIABLES
        rho = np.sum(F, 2)                                                  # calculating density
        ux = np.sum(F*cxs, 2) / rho                                         # calculating x velocity
        uy = np.sum(F*cys, 2) / rho                                         # calculating y velocity

        F[cyclinder, :] = bndryF                                            # setting velocities inside the boundary
        ux[cyclinder] = 0                                                   # omiting x velocities inside the boundary
        uy[cyclinder] = 0                                                   # omiting x velocities inside the boundary

        # COLLISION
        Feq = np.zeros(F.shape)                                             # defining force equillibrium
        for i, cx, cy, w in zip(range(NL), cxs, cys, weights):              # iterating each cell
            Feq[:, :, i] = rho * w * (                                      # calculating force equillibrium
                1 + 3*(cx*ux+cy*uy) + (9*(cx*ux+cy*uy)**2)/2 - 3*(ux**2+uy**2)/2
            )
        F = F + -(1/tau) * (F-Feq)                                          # calculating force in next step
        
        # PLOTTING
        if (it%50 == 0):                                                    # plot in periodic steps
            # CURL GRADIENT
            dfydx = ux[2: ,1:-1] - ux[0:-2, 1:-1]                           # calculating difference in x direction
            dfxdy = uy[1:-1, 2:] - uy[1:-1, 0:-2]                           # calculating difference in y direction
            curl = dfydx - dfxdy                                            # calculating curl
            pyplot.imshow(curl, cmap="bwr")                                 # plotting curl
            pyplot.title(f'Iteration {it}/{Nt}')                            # iteration counter
            pyplot.pause(0.01)                                              # pausing the visualization
            pyplot.cla()                                                    # clear the image

            # VELOCITY GRADIENT
            # pyplot.imshow(np.sqrt(ux**2 + uy**2))                           # plotting velocity
            # pyplot.title(f'Iteration {it}/{Nt}')                            # iteration counter
            # pyplot.pause(0.01)                                              # pausing the visualization
            # pyplot.cla()                                                    # clear the image

if __name__ == "__main__":
    main()
