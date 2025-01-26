import numpy as np


def step_ssf(dt, k, psi0, V, hbar, m):
    """
    Function that takes a step in SSF
    :param dt: temporal step
    :param k: vector in reciprocal space
    :param psi0: initial wave function for propagation
    :param V: potential in x space
    :param hbar: hbar constant
    :param m: mass constant
    :return: Temporal evolution of wave function in one step
    """
    psi = psi0.astype(complex)
    psi *= np.exp(-0.5 * (1j / hbar) * V * dt)
    psi_k = np.fft.fftshift(np.fft.fft(psi))
    psi_k *= np.exp(-0.5 * (1j * k**2 * dt) / m)
    psi = np.fft.ifft(np.fft.ifftshift(psi_k))
    psi = psi * np.exp(-0.5 * (1j / hbar) * V * dt)
    return psi


def ssf(Nx, Nt, dt, k, psi0, V, hbar, m):
    """
    Function to calculate the evolution of a wave function over a time range
    :param Nx: number of points in x space
    :param Nt: number of ponts in time
    :param dt: temporal step (resolution in time)
    :param k: vector in reciprocal space
    :param psi0: initial wave function for propagation
    :param V: potential in x space
    :param hbar: hbar constant
    :param m: mass constant
    :return: Matrix in which each column is the wave function at each time step
    """
    psi = psi0
    Psi_t = np.zeros((Nx, Nt), dtype=complex)
    Psi_t[:, 0] = psi
    for i in range(1, Nt):
        Psi_t[:, i] = step_ssf(dt, k, Psi_t[:, i - 1], V, hbar, m)
    return Psi_t
