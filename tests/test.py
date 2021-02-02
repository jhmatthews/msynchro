import numpy as np 
import matplotlib.pyplot as plt 
import constants as const 
import msynchro
from msynchro.units import unit

def set_mpl_defaults():
    ## FIGURE
    plt.rcParams["text.usetex"] = "True"
    ## FONT
    plt.rcParams['font.serif']=['cm']
    plt.rcParams['font.family']='serif' 
    plt.rcParams['text.latex.preamble']=r'\usepackage{amsmath}'
    plt.rcParams['font.size']=18
    plt.rcParams['xtick.labelsize']=15
    plt.rcParams['ytick.labelsize']=15
    plt.rcParams['legend.fontsize']=14
    plt.rcParams['axes.titlesize']=16
    plt.rcParams['axes.labelsize']=16
    plt.rcParams['axes.linewidth']=2
    plt.rcParams["lines.linewidth"] = 2.2
    ## TICKS
    plt.rcParams['xtick.top']='True'
    plt.rcParams['xtick.bottom']='True'
    plt.rcParams['xtick.minor.visible']='True'
    plt.rcParams['xtick.direction']='out'
    plt.rcParams['ytick.left']='True'
    plt.rcParams['ytick.right']='True'
    plt.rcParams['ytick.minor.visible']='True'
    plt.rcParams['ytick.direction']='out'
    plt.rcParams['xtick.major.width']=1.5
    plt.rcParams['xtick.minor.width']=1
    plt.rcParams['xtick.major.size']=4
    plt.rcParams['xtick.minor.size']=3
    plt.rcParams['ytick.major.width']=1.5
    plt.rcParams['ytick.minor.width']=1
    plt.rcParams['ytick.major.size']=4
    plt.rcParams['ytick.minor.size']=3

def cooling_rate(gammas, B, B_CMB=3.24e-6):
    '''
    Get the electron cooling rate from synchrotron and inverse Compton.
    returns CGS units. 

    Parameters:
        gammas  array-like
                Lorentz factors of electrons 

        B       float 
                magnetic field in Gauss 

        B_CMB   float 
                equivalent magnetic field strength of the CMB radiation energy density

    Returns:
        Cooling rate dE/dt in erg/s
    '''
    Utot = (B**2 + B_CMB**2) / 8.0 / np.pi

    x = 4.0 / 3.0 * unit.thomson * unit.c * Utot * gammas * gammas

    return (x)

def run_powerlaw_test():

    B = 6e-6
    MCSQ = unit.melec * unit.c * unit.c
    rest_mass_ev = unit.melec * unit.c * unit.c / unit.ev

        # define energy ranges 
    Emin = np.log10(10.0 * rest_mass_ev)
    Emax = np.log10(1e9 * rest_mass_ev)
    energy_edges = np.logspace(Emin,Emax,10001)
    E1s = energy_edges[:-1]
    E2s = energy_edges[1:]

    # energies centres and bins
    energies = 0.5 * (E1s + E2s)
    Ebins = E2s - E1s


    # times 
    delta_t = 0.01 * unit.myr 
    tmax = 5.0 * unit.myr 
    time = 0.0

    betas = np.arange(1,3.5,0.5)

    plt.figure(figsize=(7,5))

    for ibeta, BETA in enumerate(betas):
        ne0 = energies ** -BETA 
        ne = ne0


        gamma_edges = energy_edges / rest_mass_ev
        energy_loss_rate = cooling_rate(gamma_edges, B, 0.0) / unit.ev
        loss_rate_cen = cooling_rate(energies / rest_mass_ev, B, 0.0) / unit.ev

        time = 0.0

        while time < tmax:

            minimum_ne = 1e-100
            select = (ne > minimum_ne)
            ne[~select] = 0.0
            dndt = ne[select] / Ebins[select] * loss_rate_cen[select]
            delta_t = 0.4 * np.min(ne[select] / dndt)
            ne = msynchro.evolve.particle_evolve(energy_edges, energy_loss_rate, 0.0, 0.0, ne, delta_t)

            time += delta_t

        plt.loglog(energies, ne/ne0, ls="--", alpha=0.9, lw=2, c="C" + str(ibeta))
        C_E = loss_rate_cen / energies / energies

        ne_analytic = (1 - energies * C_E * time) ** (BETA - 2)
        #ne = ne_analytic / ne0
        ne_analytic[(energies > 1.0 / C_E / time)] = 0.0
        plt.loglog(energies, ne_analytic, ls="-", alpha=0.5, lw=3, c="C" + str(ibeta), label="$p={}$".format(BETA))

    plt.legend()
    plt.title("Analytic (solid) v TDMA (dashed)")
    plt.ylabel("$dn/dE(t) / dn/dE(t_0)$", fontsize=14)
    plt.xlabel("$E$ (eV)", fontsize=14)
    plt.ylim(1e-3,1e3)
    plt.xlim(1e9,1e12)
    plt.savefig("relative_test.png", dpi=300)


def run_delta_test():
    '''
    Run a delta function injection test.
    '''

    B = 6e-6
    MCSQ = unit.melec * unit.c * unit.c
    rest_mass_ev = unit.melec * unit.c * unit.c / unit.ev
    BETA = 2.1

    # define energy ranges 
    Emin = np.log10(10.0 * rest_mass_ev)
    Emax = np.log10(1e9 * rest_mass_ev)
    energy_edges = np.logspace(Emin,Emax,10001)
    E1s = energy_edges[:-1]
    E2s = energy_edges[1:]

    # energy centres and bin sizes
    energies = 0.5 * (E1s + E2s)
    Ebins = E2s - E1s

    # times 
    delta_t = 0.01 * unit.myr 
    tmax = 11.0 * unit.myr 
    times_write = np.arange(0,12,1) * unit.myr
    time = 0.0
    iwrite = 0

    # initialise the distribution with a delta function at 1e12 eV 
    ne = np.zeros_like(energies)
    E_inject = 1e12
    iarg = np.argmin (np.fabs(energies - E_inject))
    ne[iarg] = 1

    # get lorentz factors and energy loss rates 
    gamma_edges = energy_edges / rest_mass_ev
    energy_loss_rate = cooling_rate(gamma_edges, B, 0.0) / unit.ev
    loss_rate_cen = cooling_rate(energies / rest_mass_ev, B, 0.0) / unit.ev

    plt.figure(figsize=(7,5))

    # evolve the population for 11 Myr, plotting every 1 Myr 
    while time < tmax:

        if time >= times_write[iwrite]:
            plt.loglog(energies / E_inject, ne, drawstyle="steps")
            iwrite+=1

        minimum_ne = 1e-50
        select = (ne > minimum_ne)
        ne[~select] = 0.0
        dndt = ne[select] / Ebins[select] * loss_rate_cen[select]

        # two possible choices for delta_t 
        delta_t = 0.4 * np.min(ne[select] / dndt)
        #delta_t = 0.4 * np.min(Ebins[select] / loss_rate_cen[select])

        # evolve particle distribution 
        ne = msynchro.evolve.particle_evolve(energy_edges, energy_loss_rate, 0.0, 0.0, ne, delta_t)

        time += delta_t


    E = E_inject
    iwrite = 0
    Enew = np.zeros_like(times_write)
    dt = 0.0001 * unit.myr
    for t in np.arange(0, tmax, dt):
        E = E - (cooling_rate(E/rest_mass_ev, B, 0.0) / unit.ev * dt)
        if t >= times_write[iwrite]:
            Enew[iwrite] = E
            iwrite += 1

    g0 = E_inject * unit.ev / unit.melec_csq
    D0 = cooling_rate(g0, B, B_CMB = 0.0) / unit.melec_csq / g0**2

    ##print (D0)
    #D0 = 4.0 * unit.thomson * B * B / (24.0 * np.pi * unit.melec * unit.c)
    #print (D0)

    g0 = E_inject * unit.ev / unit.melec_csq

    #print (D0 * g0 * times_write)
    for i, t in enumerate(times_write[:-1]):
        g = g0 / (1 + (D0 * g0 * t))
        E = g * unit.melec_csq / unit.ev
        plt.vlines([E / 1e12], 1e-3, 1, ls= "--", color="C"+str(i))

    # adjust plot 
    plt.ylim(1e-3,1)
    plt.xlim(1e-2,2)
    plt.xlabel("$E / E_0$")
    plt.ylabel("$dN/dE$")
    plt.savefig("delta_test.png")

if __name__ == "__main__":
    set_mpl_defaults()
    run_delta_test()
    run_powerlaw_test()
