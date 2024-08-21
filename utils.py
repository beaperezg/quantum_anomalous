import numpy as np
from typing import Optional


def probability(Ncells: int, nmax: int, slist: np.ndarray, v: np.ndarray) -> np.ndarray:
    sumat = np.zeros((np.size(slist), 2 * Ncells * (nmax + 1)))

    for q in range(np.size(slist)):
        s = slist[q]
        for k in range(2 * Ncells * (nmax + 1)):
            suma = 0
            for l in range(nmax):
                suma += np.abs(v[l * 2 * Ncells + s, k] ** 2)
            sumat[q, k] = suma

    return np.sum(sumat, axis=0)


def nphotons_operator(Ncells: int, nmax: int) -> np.ndarray:
    nph = np.zeros([2 * Ncells * (nmax + 1), 2 * Ncells * (nmax + 1)])

    for i in range(nmax + 1):
        nph[2 * Ncells * i: 2 * Ncells * (i + 1), 2 * Ncells * i: 2 * Ncells * (i + 1)] = i * np.identity(2 * Ncells)

    return nph


def expval(st: np.ndarray, op: np.ndarray) -> float:
    # expected value of an operator 'op' in a given state 'st'
    return np.dot(np.conjugate(st), np.dot(op, st))


def calculate_phase(avec: np.ndarray, nphases: Optional[int] = None, window: Optional[float] = 1e-10) -> np.ndarray:
    """
    This function calculates the Zak phase, with the result being in the interval [-pi, pi)
    """
    nk = np.shape(avec)[0]
    nh = np.shape(avec)[1]

    if nphases is None:
        nphases = nh
    products = []

    for j in range(nphases):
        prod = 1
        for k in range(nk - 1):
            prod *= np.dot(np.conjugate(avec[k, :, j]), avec[k + 1, :, j])
        products.append(prod)

    phases = - np.imag(np.log(products))
    cond1 = np.pi * (1 - window) < phases
    cond2 = phases < np.pi * (1 + window)
    phases[cond1 & cond2] = -np.pi

    return phases
