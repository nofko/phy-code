import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 50.0  # Spatial domain
T = 1.0   # Total simulation time
Nx = 500  # Number of spatial points
Nt = 500  # Number of time steps
dt = T / Nt
dx = L / Nx

# Initialize arrays
x = np.linspace(0, L, Nx)
t = np.linspace(0, T, Nt)
u = np.zeros((Nt, Nx))

# Initial condition: Single soliton
A = 0.2
x0 = L / 2.0
u[0, :] = A * np.exp(-(x - x0)**2 / (2.0 * A**2))

# Time-stepping using split-step method
for n in range(Nt - 1):
    # First half-step (linear part)
    u_half = u[n, :] - 0.5 * (np.roll(u[n, :], -1) - np.roll(u[n, :], 1)) * (dt / dx)

    u_half = np.fft.fft(u_half)

    # Second half-step (nonlinear part)
    k = 2.0 * np.pi / L * np.arange(-Nx / 2, Nx / 2)
    u_half *= np.exp(-1j * k**3 * dt)

    # Inverse Fourier transform
    u[n + 1, :] = np.fft.ifft(u_half).real

# Plot the results
plt.figure(figsize=(10, 5))
plt.contourf(x, t, u, cmap='viridis')
plt.colorbar()
plt.title('Korteweg-de Vries (KdV) Equation Simulation')
plt.xlabel('Spatial Domain')
plt.ylabel('Time Domain')
plt.show()
