Examples
--------------------------------------


Synchrotron Calculation
================================================
A synchrotron calculation can be carried out using a few simple commands 

.. code:: python

    import numpy as np
    import msynchro
    frequencies = np.logspace(8,12,1000)
    energies = np.logspace(8,12,1000)
    ne = energies ** -2.0 
    B = 1e-6
    spec = msynchro.Ptot(frequencies, energies, ne, B)

Particle evolution calculation
================================================
A population of particles can be evolved using the tridiagonal matrix algorithm by specifying arrays holding the initial particle states, and the cooling rates. 