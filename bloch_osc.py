# BLOCH OSCILLATIONS
# Autor: Aaron Trinidad
# GitHub user: aaron-trinidad

import numpy as np
import matplotlib.pyplot as plt
from ssfourier import ssf

# Hamiltonian for this problem: H = -\frac{\hbar^2}{2m}\frac{d^2}{dx^2}+V_0 cos(2*pi x / a)+Fx

# Definition of Parameters and Initial Conditions
hbar = 1
m = 1
dx = 0.004
dt = 0.001
x = np.arange(-10, 10 + dx, dx)  # Real space
t = np.arange(0, 0.75 + dt, dt)
Nx = len(x)
Nt = len(t)
k = 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(Nx, d=dx))  # reciprocal space
dk = np.abs(k[1] - k[0])

# Initial conditions of the wave
kx = 100
sigma = 0.9
x0 = 0.0
C = 1.0 / (sigma * np.sqrt(np.pi))
psi0 = (
    np.sqrt(C) * np.exp(-((x - x0) ** 2) / (2.0 * sigma**2)) * np.exp(1j * kx * x)
)  # Initial wave packet


# Periodic Potential
def V(V0f, xf):
    return np.cos(2 * np.pi * xf) * V0f


V0 = 0.1
v = V(V0, x)
# Uncomment the force you want
# F = 500000  # F strong
F = 40000  # F weak
# F = 200  # No oscillations
Fx = F * x

Potencial = v + Fx  # Full potential = Periodic potential plus External Force

# Time evolution of the initial condition
Psi_xt = ssf(Nx, Nt, dt, k, psi0, Potencial, hbar, m)
prob_Psi_xt = np.abs(Psi_xt) ** 2
Psi_kt = np.fft.fft(Psi_xt, axis=0)  # Solution in reciprocal space
Psi_kt = np.fft.fftshift(Psi_kt, axes=0)
prob_Psi_kt = np.abs(Psi_kt) ** 2

expect_k = np.sum(k[:, np.newaxis] * prob_Psi_kt, axis=0) * dk

# 2D Plots

time_indices = [20, 50, 80]

# Creating subplots
fig, axs = plt.subplots(1, 3, figsize=(10, 4))

# Subplots for Psi(x,t) --real space--
for i, idx in enumerate(time_indices):
    axs[i].plot(x, prob_Psi_xt[:, idx], color="b")
    axs[i].set_title(f"$\\Psi(x, t={t[idx]:.2f})$")
    axs[i].set_xlabel("Position (x)")
    axs[i].set_ylabel(r"$|\Psi(x,t)|^2$")
    axs[i].set_xlim(x[0], x[-1])
    axs[i].set_ylim(0, np.max(prob_Psi_xt) * 1.1)

plt.tight_layout()
# plt.savefig("images/evol_xspace.png", dpi=300)

# Subplots for Psi(k,t) --reciprocal-space--
fig, axs_k = plt.subplots(1, 3, figsize=(10, 4))

for i, idx in enumerate(time_indices):
    axs_k[i].plot(k, prob_Psi_kt[:, idx], color="r")
    axs_k[i].set_title(f"$\\Psi(k, t={t[idx]:.2f})$")
    axs_k[i].set_xlabel("Wave number (k)")
    axs_k[i].set_ylabel(r"$|\Psi(k,t)|^2$")
    axs_k[i].set_xlim(k[0], k[-1])
    axs_k[i].set_ylim(0, np.max(prob_Psi_kt) * 1.1)

plt.tight_layout()
# plt.savefig("images/evol_kspace.png", dpi=300)

# Density Probability
# Heat map
extend = (np.min(t), np.max(t), np.min(x), np.max(x))
plt.figure()
plt.imshow(prob_Psi_xt, aspect="auto", cmap="viridis", origin="lower", extent=extend)
plt.colorbar(label="$|\\Psi(x, t)|^2$")
plt.title("Bloch Oscillations")
plt.xlabel("t")
plt.ylabel("x")
# plt.savefig("images/densityheat.png", dpi=300)

# Expected value of the moment as a function of time
plt.figure()
plt.plot(t, expect_k, "--")
plt.title(r"$\langle k \rangle$")
plt.xlabel("t")
# plt.savefig("images/expect_k.png", dpi=300)

plt.show()
