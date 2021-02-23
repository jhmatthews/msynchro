Description
--------------------------------------

msynchro is a simple python package that allows calculation of optically thin synchrotron emission from arbitrary electron distributions. It also include a tridiagonal matrix algorithm for evolving distributions of particles subject to cooling. This page contains a brief description of the main methods.

Synchrotron Calculation
================================================

.. todo:: write description.


Tridiagonal matrix algorithm
================================================
The TDM algorithm is used to solve a continuity equation of the form. This discussion broadly follows that given by `Chiaberge and Ghisellini (1999) <https://ui.adsabs.harvard.edu/abs/1999MNRAS.306..551C/abstract>`_.  The algorithm is intended to solve equations of the following form

.. math::

	\frac{dn(E,t)}{dE} = \frac{d}{dE}\left[\frac{dE}{dt} n(E,t) \right] + S(E,t) - \frac{n}{\tau_{\rm loss}(E,t)}

which can be written in discrete form as 

.. math::

	\frac{n_j^{i+1} - n_j^i}{\Delta t} = \frac{F^{i+1}_{j+1/2} - F^{i+1}_{j-1/2}}{\Delta E} + Q^i_j - \frac{n_j^{i+1}}{\tau_{j,{\rm loss}}}

where the indices :math:`j` and :math:`i` denote the energy bin and time step, respectively. 

.. todo:: finish this discussion.