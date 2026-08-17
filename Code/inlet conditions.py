import numpy as np
from scipy.optimize import fsolve

# Inputs
mdot = 11.67          # kg/s
Pt = 344740.0         # Pa
Tt = 709.0            # K

rh = 0.32576          # m
rt = 0.36576          # m

gamma = 1.333
R = 287.0

# Annulus area
A = np.pi * (rt**2 - rh**2)

def residual(Ca):
    # Static temperature
    T = Tt - Ca**2/(2 * gamma * R / (gamma - 1))

    # Static pressure
    P = Pt * (T/Tt)**(gamma/(gamma-1))

    # Density
    rho = P/(R*T)

    # Mass flow error
    return rho * A * Ca - mdot

# Initial guess
Ca0 = 120.0

Ca_solution = fsolve(residual, Ca0)[0]

# Compute final properties
T = Tt - Ca_solution**2/(2 * gamma * R / (gamma - 1))
P = Pt * (T/Tt)**(gamma/(gamma-1))
rho = P/(R*T)

print(f"Required axial velocity = {Ca_solution:.3f} m/s")
print(f"Static temperature = {T:.2f} K")
print(f"Static pressure = {P:.2f} Pa")
print(f"Density = {rho:.4f} kg/m³")