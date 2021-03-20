from scipy import integrate, special
from scipy.interpolate import interp1d
import numpy as np
from msynchro.units import unit


class lookup:
    """
    DEPRECATED: wrapper for F(X) lookup table

    :meta private:
    """

    def __init__(self, load=True, fname="Fxlookup.npy"):

        self.x_array = np.logspace(-10, 10, 10000)
        if load:
            self.ff = np.load(fname)
        else:

            self.ff = np.zeros_like(self.x_array)
            for i, x in enumerate(self.x_array):
                self.ff[i] = F(x)

            np.save(fname, self.ff)

        self.interp_func = interp1d(
            self.x_array, self.ff, kind="quadratic", fill_value="extrapolate"
        )


def fx_approximation(x):
    """
    DEPRECATED: An approximate form of F(x) as given by Aharonian et al.,

    :meta private:
    """
    f = 2.15 * (x ** (1.0 / 3.0))
    f *= (1 + 3.06 * x) ** (1.0 / 6.0)
    num = 1.0 + (0.884 * x ** (2.0 / 3.0)) + (0.471 * x ** (4.0 / 3.0))
    denom = 1.0 + (1.64 * x ** (2.0 / 3.0)) + (0.974 * x ** (4.0 / 3.0))
    f *= num / denom
    f *= np.exp(-x)
    return f


def F(x):
    """
    DEPRECATED: This is F(x) defined in equation 6.31c in Rybicki and Lightman.

    F(x) = x * integral(K_5/3(x)dx) from x to infinity.

    where K_5/3 is the modified Bessel function of order 5/3.

    Parameters:
        x       float
                nu/nu_c, frequency normalised to critical frequency

    :meta private:
    """

    # at large x, return eq 6.34b from Rybicki and Lightman (could perhaps use this earlier)
    if x > 1e5:
        answer = np.sqrt(np.pi / 2.0) * np.exp(-x) * np.sqrt(x)
    else:
        answer = (
            x
            * integrate.quad(lambda i: special.kv(5.0 / 3, i), x, np.inf, limit=200)[0]
        )
    return answer


def nu_crit(gamma, B):
    '''
    Calculate the critical frequency of a synchrotron electron

    Parameters:
        gamma       array-like, float
                    Lorentz factors of electrons

        B           float
                    magnetic field strength

    Returns:
        nu_c        float 
                    critical frequency of a synchrotron electron in Hz
    '''
    nu_c = gamma * gamma * nu_cyclotron(B)
    return (nu_c)

def nu_cyclotron(B):
    '''
    Calculate the cyclotron frequency or gyrofrequency

    Parameters:
        B           float
                    magnetic field strength

    Returns:
        nu          float 
                    cyclotron frequency in Hz

    '''
    nu = unit.e * B / 2.0 / np.pi / unit.melec / unit.c
    return (nu)

def psynch(gamma, nu, B):
    """
    equation 13 from Chiaberge & Ghisellini. This is the single
    particle synchrotron emissivity j_nu
    averaged over an isotropic distribution of pitch angles.

    Parameters:
        gamma       array-like
                    Lorentz factors of electrons

        nu          float
                    frequency in Hz

        B           float
                    magnetic field in Gauss
    """
    nu_B = unit.e * B / 2.0 / np.pi / unit.melec / unit.c
    t = nu / (3.0 * gamma * gamma * nu_B)

    x = 3.0 * np.sqrt(3.0) / np.pi * unit.thomson * unit.c * B * B / 8.0 / np.pi
    x *= t * t / nu_B

    # get the modified Bessel functions
    K13 = special.kv(1.0 / 3.0, t)
    K43 = special.kv(4.0 / 3.0, t)
    K43sq = K43 * K43
    K13sq = K13 * K13
    kterm = (K13 * K43) - (0.6 * t * (K43sq - K13sq))

    return x * kterm


def Ptot(nus, energies, ne, Bfield):
    """
    Get synchrotron spectrum for a given set of frequencies 
    from a differential spectrum of electrons ne=dN/dE.


    Parameters:
        nus         array-like
                    frequencies in Hz

        energies    array-like
                    energies of electrons - units must match ne array

        ne          array-like
                    differential energy spectrum - units must match energies array

        B           float
                    magnetic field in Gauss

    Returns:
        Ptot        array-like
                    synchrotron spectrum with same shape as 
                    nus input array
    """
    # array to store spectrum
    pnu = np.zeros_like(nus)

    # convert differential spectrum to CGS units
    dn_by_dE_cgs = ne / unit.ev

    # convert to dn/dgamma
    dn_by_dgamma = dn_by_dE_cgs * unit.melec_csq
    gamma = energies * unit.ev / unit.melec_csq

    for inu, nu in enumerate(nus):
        # get power for each gamma bin
        power = psynch(gamma, nu, Bfield)

        # integrate over distribution
        pnu[inu] = -np.trapz(gamma, dn_by_dgamma * power)

    return pnu
