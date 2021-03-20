import numpy as np 
class units:
    '''
    class containing some units. Should probably use astropy units 
    but I find them a bit annoying. Example to get the rest mass energy
    of an electron in eV would be

        from msynchro.units import unit 
        E = unit.melec_csq / unit.ev

    Includes the following quantities in CGS units:

    kpc         1 kiloparsec 
    pc          1 parsec
    c           speed of light in vacuum
    yr          1 year
    myr         1 Megayear
    kyr         1 kyr
    radian      degrees in 1 radian
    msol        solar mass
    mprot       proton mass
    melec       electron mass
    melec_csq   electron rest mass energy
    mprot_csq   proton rest mass energy
    e           fundamental charge 
    ev          electron volt
    kb          boltzmann constant
    h           plank constant 
    hbar        reduced plank constant (h/2pi)
    g           gravitational constant G
    hbar_c      hbar * c
    alpha       fine structure constant
    thomson     Thomson cross section
    '''
    def __init__(self):
        self.kpc = 3.086e21
        self.pc = 3.086e18
        self.c = 2.997925e10
        self.yr = 3.1556925e7
        self.myr = 3.1556925e13
        self.kyr = 3.1556925e10
        self.radian = 57.29577951308232
        self.msol = 1.989e33
        self.mprot = 1.672661e-24
        self.melec = 9.10956e-28
        self.melec_csq = self.melec * self.c * self.c
        self.mprot_csq = self.mprot * self.c * self.c
        self.e = 4.8035e-10     # fundamental charge 
        self.ev = 1.602192e-12  # electron volts in CGS
        self.kb = 1.38062e-16   # boltzmann 
        self.h = 6.6262e-27     # plank 
        self.hbar = self.h / 2.0 / np.pi      
        self.g = 6.670e-8       # gravitational 
        self.hbar_c = self.hbar * self.c
        self.alpha = self.e * self.e / self.hbar_c
        self.thomson = 0.66524e-24

# import this with 
# from msynchro.units import unit 
unit = units() 