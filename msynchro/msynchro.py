from scipy import integrate,special
from scipy.interpolate import interp1d
import numpy as np 
from msynchro.units import unit

class lookup:
    '''
    wrapper for F(X) lookup table
    DEPRECATED
    '''
    def __init__(self, load = True, fname = "Fxlookup.npy"):
        
        self.x_array = np.logspace(-10,10,10000)
        if load:
            self.ff = np.load(fname)
        else:

            self.ff = np.zeros_like(self.x_array)
            for i, x in enumerate(self.x_array):
                self.ff[i] = F(x)

            np.save(fname, self.ff)

        self.interp_func = interp1d(self.x_array, self.ff, kind="quadratic", fill_value="extrapolate")


def fx_approximation(x):
    ''' 
    An approximate form of F(x) as given by Aharonian et al.,
    '''
    f = 2.15 * (x**(1./3.))
    f *= (1 + 3.06 * x) ** (1./6.)
    num = 1.0 + (0.884  * x**(2./3.)) + (0.471  * x**(4./3.))
    denom = 1.0 + (1.64  * x**(2./3.)) + (0.974  * x**(4./3.))
    f *= num / denom 
    f *= np.exp(-x)
    return (f)


def F(x):
    """ 
    This is F(x) defined in equation 6.31c in Rybicki and Lightman.
        
    F(x) = x * integral(K_5/3(x)dx) from x to infinity.

    where K_5/3 is the modified Bessel function of order 5/3.

    Parameters:
        x       float
                nu/nu_c, frequency normalised to critical frequency

    DEPRECATED
    """

    # at large x, return eq 6.34b from Rybicki and Lightman (could perhaps use this earlier)
    if x > 1e5: 
        answer = np.sqrt(np.pi / 2.0) * np.exp(-x) * np.sqrt(x)
    else:
        answer = x * integrate.quad(lambda i: special.kv(5./3, i), x, np.inf, limit=200)[0]
    return (answer)


# def nu_crit(B, gamma):
#     '''
#     Get the critical frequency of an electron with
#     Lorentz factor gamma in a magnetic field B, in microGauss
#     '''
#     nu_c = 3.0 * gammas[i] * gammas[i] * unit.e * Bfield / unit.melec / unit.c / 2.0 * sinalpha / 2.0 / np.pi

def psynch_slower (energies, nu, Bfield, lookup = None):
    alphas = np.linspace(np.pi/100.0,np.pi,100)
    sinalpha = np.sin(alphas)

    gammas = energies * unit.ev / unit.melec_csq
    # nu_c = 3.0 * gammas * gammas * E * Bfield / MELEC / C / 2.0 * sinalpha
    term1 = np.sqrt(3)
    # term2 = E * E * E * Bfield * sinalpha / MELEC / C / C 
    integral = np.zeros_like(gammas)

    for i in range(len(gammas)):
        nu_c = 3.0 * gammas[i] * gammas[i] * unit.e * Bfield / unit.melec / unit.c / 2.0 * sinalpha / 2.0 / np.pi
        term2 = unit.e * unit.e * unit.e * Bfield * sinalpha / unit.melec_csq 
        #def array_map(x):
        x = nu/nu_c
        if lookup == None:
            term3 = np.array(list(map(F, x)))
        else:
            term3 = lookup.interp_func(x)
        #term3 = x * fx_approximation(x)

        # //integrand = term2 * term3
        integral[i] = np.trapz(alphas, 0.5 * sinalpha * term2 * term3)
    #term3 = np.array([F(nu/nuc) for nuc in nu_c])
    #print (integral)
    return integral * term1


def psynch (gamma, nu, B):
    '''
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
    '''
    nu_B = unit.e * B / 2.0 / np.pi / unit.melec / unit.c
    t = nu / (3.0  * gamma * gamma * nu_B)

    x = 3.0 * np.sqrt(3.0) / np.pi * unit.thomson * unit.c * B * B / 8.0 / np.pi 
    x *= t * t / nu_B

    # get the modified Bessel functions 
    K13 = special.kv(1./3., t)
    K43 = special.kv(4./3., t)
    K43sq = K43 * K43 
    K13sq = K13 * K13 
    kterm = (K13 * K43) - (0.6 * t * (K43sq - K13sq))

    return (x * kterm)




def Ptot(nus, energies, ne, Bfield, lookup=None):
    '''
    get synchrotron spectrum for a given 
    '''
    # array to store spectrum 
    pnu = np.zeros_like(nus)

    # convert differential spectrum to CGS units 
    dn_by_dE_cgs = ne / unit.ev

    # convert to dn/dgamma 
    dn_by_dgamma = dn_by_dE_cgs * unit.melec_csq 
    gamma = energies * unit.ev / unit.melec_csq 

    for inu, nu in enumerate(nus):
        # get power for each gamma bin 
        power = psynch (gamma, nu, Bfield)

        # integrate over distribution 
        #power = psynch(energies, nu, Bfield, lookup = lookup)
        #power = jsynch(energies, ncr, nus, Bfield, np.max(energies), len(energies))
        #print (gamma, dn_by_dgamma, power)
        pnu[inu] = -np.trapz(gamma, dn_by_dgamma * power) 


    return pnu