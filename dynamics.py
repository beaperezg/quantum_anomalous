import qutip as qt
import numpy as np
import scipy as sp
from typing import Optional, Tuple
from hamiltonians import htotal
from tqdm.auto import tqdm


def _evolution_operator(time: np.ndarray, hop: np.ndarray, g: float,
                        Omega: float, Ncells: int, nmax: int) -> qt.Qobj:
    dt = time[1] - time[0]
    h = htotal(hop, g, Omega, Ncells, nmax)
    U_dt = (-1j * h * dt).expm()

    return U_dt


def calculate_dynamics(psi_0_elec: qt.Qobj, psi_0_phot: qt.Qobj, time: np.ndarray, hop: np.ndarray, g: float,
                       Omega: float, Ncells: int, nmax: int) -> Tuple[np.ndarray, np.ndarray]:
    U_dt = _evolution_operator(time, hop, g, Omega, Ncells, nmax)
    psi_initial = qt.tensor(psi_0_phot, psi_0_elec)

    prob_pho = [abs(np.array(psi_0_phot).flatten()) ** 2]
    prob_elec = [abs(np.array(psi_0_elec).flatten()) ** 2]
    psi_temp = psi_initial

    pbar = tqdm(total=len(time) - 1, desc='Computing dynamics')
    for i in range(1, len(time)):
        psi_temp = U_dt * psi_temp
        prob_pho.append(np.diag(psi_temp.ptrace(0)).real)
        prob_elec.append(np.diag(psi_temp.ptrace(1)).real)
        pbar.update()
    pbar.close()

    return np.array(prob_pho), np.array(prob_elec)


def discrete_fourier_transform(curve: np.ndarray, time: np.ndarray, height: Optional[int] = 450) -> Tuple[
    np.ndarray, np.ndarray, np.ndarray]:
    signal = np.real(curve)
    fourier = np.abs(np.fft.rfft(signal))
    n = signal.size
    d = time[1] - time[0]
    freq = np.fft.rfftfreq(n, d)
    per = 1 / freq
    peaks, _ = sp.signal.find_peaks(fourier, height=height)

    return per, fourier, peaks
