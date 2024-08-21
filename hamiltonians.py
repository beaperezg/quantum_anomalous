import numpy as np
import qutip as qt

"""
################################## ARTIFICIAL COUPLING TO THE CAVITY ##########################################
############################################### REAL SPACE ####################################################
"""


def _ssh(hop: np.ndarray, Ncells: int) -> qt.Qobj:
    diag1 = [np.mod(i + 1, 2) * hop[0] + np.mod(i, 2) * hop[1] for i in range(2 * Ncells - 1)]
    h = np.diag(diag1, 1) + np.diag(diag1, -1)
    return qt.Qobj(h)


def _sshCavity(hop: np.ndarray, Ncells: int, nmax: int) -> qt.Qobj:
    return qt.tensor(qt.identity(nmax + 1), _ssh(hop, Ncells))


def _vint(Ncells: int) -> qt.Qobj:
    diag1 = [np.mod(i + 1, 2) * 1 + np.mod(i, 2) * (-1) for i in range(2 * Ncells - 1)]
    h = np.diag(diag1, 1) + np.diag(diag1, -1)
    return qt.Qobj(h)


def _hinteraction(g: float, Ncells: int, nmax: int) -> qt.Qobj:
    return qt.tensor(g * (qt.create(nmax + 1) + qt.destroy(nmax + 1)), _vint(Ncells))


def _cavity(Omega: float, Ncells: int, nmax: int) -> qt.Qobj:
    return qt.tensor(Omega * qt.create(nmax + 1) * qt.destroy(nmax + 1), qt.identity(2 * Ncells))


def htotal(hop: np.ndarray, g: float, Omega: float, Ncells: int, nmax: int) -> qt.Qobj:
    return _sshCavity(hop, Ncells, nmax) + _cavity(Omega, Ncells, nmax) + _hinteraction(g, Ncells, nmax)


"""
################################## ARTIFICIAL COUPLING TO THE CAVITY ##########################################
################################################ K SPACE ######################################################
"""


def sshk(hop: np.ndarray, k: float) -> qt.Qobj:
    hsshk = np.zeros([2, 2], dtype='complex')
    hsshk[0, 1] = hop[0] + hop[1] * np.exp(1j * k)
    hsshk[1, 0] = hop[0] + hop[1] * np.exp(-1j * k)
    return qt.Qobj(hsshk)


def _sshkCavity(hop: np.ndarray, k: float, nmax: int) -> qt.Qobj:
    return qt.tensor(qt.identity(nmax + 1), sshk(hop, k))


def _vintk(k: float) -> qt.Qobj:
    matrixK = np.zeros([2, 2], dtype='complex')
    matrixK[0, 1] = 1 - np.exp(1j * k)
    matrixK[1, 0] = 1 - np.exp(-1j * k)
    return qt.Qobj(matrixK)


def _hinteractionk(g: float, k: float, nmax: int) -> qt.Qobj:
    return qt.tensor(g * (qt.create(nmax + 1) + qt.destroy(nmax + 1)), _vintk(k))


def _cavityk(Omega: float, nmax: int) -> qt.Qobj:
    return qt.tensor(Omega * qt.create(nmax + 1) * qt.destroy(nmax + 1), qt.identity(2))


def htotalK(hop: np.ndarray, k: float, g: float, Omega: float, nmax: int) -> qt.Qobj:
    return _sshkCavity(hop, k, nmax) + _cavityk(Omega, nmax) + _hinteractionk(g, k, nmax)
