Description
--------------------------------------

msynchro is a simple python package that allows calculation of optically thin synchrotron emission from arbitrary electron distributions. It also include a tridiagonal matrix algorithm for evolving distributions of particles subject to cooling. This page contains a brief description of the main methods.

Tridiagonal matrix algorithm
================================================
The TDM algorithm is used to solve a continuity equation of the form. This discussion broadly follows that given by `Chiaberge and Ghisellini (1999) <https://ui.adsabs.harvard.edu/abs/1999MNRAS.306..551C/abstract>`_.  The algorithm is intended to solve equations of the following form

.. math::

	\frac{dn(E,t)}{dE} = \frac{d}{dE}\left[\dot{E} n(E,t) \right] + Q(E,t) - \frac{n}{\tau_{\rm loss}(E,t)},

where  :math:`n(E,t)` is the differential spectrum at energy :math:`E` and time :math:`t`. :math:`\tau_{\rm loss}` is an escape or loss timescale, :math:`Q(E,t)` is the source term, and :math:`\dot{E}` is the cooling rate. This differential equation can be written in discrete form as 

.. math::

	\frac{n_j^{i+1} - n_j^i}{\Delta t} = \frac{F^{i+1}_{j+1/2} - F^{i+1}_{j-1/2}}{\Delta E} + Q^i_j - \frac{n_j^{i+1}}{\tau_{j,{\rm loss}}}

where the indices :math:`j` and :math:`i` denote the energy bin and time step, respectively. This equation can be written in tridiagonal matrix form 

.. math::

	c_j n_{j+1}^{i+1} + b_j n_{j}^{i+1} + a_j n_{j-1}^{i+1} = S^i_j 

where :math:`S^i_j = N^i_j + Q^i_j \Delta t` and 

.. math::

	a_j & = & 0 \nonumber \\
	b_j & = & 1+ \frac{\Delta t}{\tau_{\rm loss}} + \frac{\Delta t \, \dot{E}_{j-1/2}}
	{\Delta E_j} \\
	c_j & = & -\frac{\Delta t \, \dot{E}_{j+1/2}}{\Delta E_j}. \nonumber

This equation can then be written at each time step as an :math:`(m \times m)` tridiagonal matrix, with the matrix equation to solve as follows (dropping the :math:`i` superscript):

.. math::
	\begin{bmatrix}
	   {b_1} & {c_1} & {   } & {   } & { 0 } \\
	   {a_2} & {b_2} & {c_2} & {   } & {   } \\
	   {   } & {a_3} & {b_3} & \ddots & {   } \\
	   {   } & {   } & \ddots & \ddots & {c_{m-1}}\\
	   { 0 } & {   } & {   } & {a_m} & {b_m}\\
	\end{bmatrix}
	\begin{bmatrix}
	   {n_1 }  \\
	   {n_2 }  \\
	   {n_3 }  \\
	   \vdots   \\
	   {n_m }  \\
	\end{bmatrix}
	=
	\begin{bmatrix}
	   {S_1 }  \\
	   {S_2 }  \\
	   {S_3 }  \\
	   \vdots   \\
	   {S_m }  \\
	\end{bmatrix}

which is solved using the `tridiagonal matrix algorithm <https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm>`_. The actual algorithm is written as a C extension and an example of it's use is found under :doc:`examples`.

Synchrotron Calculation
================================================

We calculate the synchrotron emissivity :math:`\epsilon_\nu` of an electron distribution :math:`n(E)` with 

.. math::
	\epsilon_\nu = \int^{E_{max}}_{E_{min}} n(E) \, P_\nu(E) \, dE 

Here :math:`P_s(\nu,\gamma)`  is the single particle synchrotron 
emissivity averaged over an isotropic distribution of pitch angles, given by 
(`Crusius & Schlickeiser 1986 <https://ui.adsabs.harvard.edu/abs/1986A%26A...164L..16C/abstract>`_ )

.. math:: 
	P_\nu(\gamma) & = & \frac{3 \sqrt{3}}{\pi} \frac{\sigma_{T} c 
	U_{B}}{\nu_{B}} t^{2} \\ 
	&  \times & \left\{K_{4/3}(t) \, K_{1/3}(t)- \frac{3}{5} t \, 
	 [K_{4/3}^{2}(t)-K_{1/3}^{2}(t)] \right\}

where :math:`t=\nu/(3\gamma^2\nu_B)`, :math:`\nu_B = e B/(2\pi m_e c)`, and the :math:`K` quantities are the usual modified Bessel functions.

.. todo:: write description.