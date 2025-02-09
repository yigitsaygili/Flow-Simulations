import numpy as np
import matplotlib.pyplot as plt

def jukowsky_airfoil(s_x, s_y, AoA):
    # FLOW PROPERTIES
    rho = 1.225
    v_inf = 20
    v = v_inf / v_inf
    theta = AoA * np.pi / 180
    
    # CIRCLE DEFINITION
    s = s_x + 1j * s_y
    r = 0.5
    lambda_ = r - s
    
    # CIRCULATION
    beta = theta
    k = 2 * r * v * np.sin(beta)
    Gamma = k / (2 * np.pi)
    
    # COMPLEX ASYMPTOTIC SPEED 
    w = v * np.exp(1j * theta)
    
    # GENERATING MESH
    toll = +5e-2
    x = np.arange(-2, 2.1, 0.1)
    y = x
    X, Y = np.meshgrid(x, y)
    z = X + 1j * Y
    
    # CALCULATIONS
    # Inside-circle points are Excluded!
    for a in range(len(x)):
        for b in range(len(y)):
            if abs(z[a, b] - s) <= r - toll:
                z[a, b] = np.nan
    
    # AERODYNAMIC POTENTIAL
    z_minus_s = z - s
    z_minus_s = np.where(np.abs(z_minus_s) < 1e-10, np.nan, z_minus_s)
    
    f = w * z + (v * np.exp(-1j * theta) * r**2) / z_minus_s + 1j * k * np.log(z_minus_s)
    
    # JOUKOWSKI TRANSFORMATION
    z = np.where(np.abs(z) < 1e-10, np.nan, z)
    J = z + lambda_**2 / z
    
    # GRAPHIC - Circle and Joukowski Airfoil
    angle = np.arange(0, 2 * np.pi, 0.1)
    z_circle = r * (np.cos(angle) + 1j * np.sin(angle)) + s
    z_airfoil = z_circle + lambda_**2 / z_circle
    
    # KUTTA JOUKOWSKI THEOREM
    L = v_inf * rho * Gamma
    L_str = str(L)

    x_coord = np.concatenate((np.array([max(np.real(z_airfoil)) + np.real(z_airfoil)]).T, np.array([[max(np.real(z_airfoil)) + np.real(z_airfoil[0])]])))
    y_coord = np.concatenate((np.imag(z_airfoil).reshape(-1, 1), [[np.imag(z_airfoil[0])]]))
    
    # VISUALIZATION
    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    plt.contourf(np.real(z), np.imag(z), np.imag(f), np.arange(-2, 2.2, 0.2), cmap='RdYlBu')  # Filled contours
    plt.fill(np.real(z_circle), np.imag(z_circle), 'k')
    plt.axis('equal')
    plt.axis([-2, 2, -2, 2])
    plt.title(f'Flow Around a Circle.   Lift:  {float(L_str):.3f}  [N/m]')

    plt.subplot(1, 2, 2)
    plt.contourf(np.real(J), np.imag(J), np.imag(f), np.arange(-2, 2.2, 0.2), cmap='RdYlBu')
    plt.fill(np.real(z_airfoil), np.imag(z_airfoil), 'k')
    plt.axis('equal')
    plt.axis([-2, 2, -2, 2])
    plt.title(f'Flow Around the Corresponding Airfoil.   Lift:  {float(L_str):.3f}  [N/m]')
    
    plt.tight_layout()
    plt.show()

    return x_coord, y_coord

s_x = 0.05
s_y = 0.05
AoA = 10

x_coord, y_coord = jukowsky_airfoil(s_x, s_y, AoA)
