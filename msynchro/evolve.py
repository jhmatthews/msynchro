from msynchro.units import unit
import msynchro
import numpy as np 

def particle_evolve(energy_edges, energy_loss_rate, tloss_discrete, source, n_i, dt):
	'''
	Evolve a particle distribution for one time step. any units allowed as
	long as all energy and time units are consistent. Using a TDMA solver to evolve
	the distribution following the method of Chiaberge & Chisellini 1999.

	Parameters:
		energy_edges 		array-like
							edges of energy bins, len(N)

		energy_loss_rate 	array-like or float
							dE/dt at edges of energy bins, len (N)

		tloss_discrete 		array-like or float
							tau_loss in each bin, len (N-1)

		source 				array-like
							source term, len (N-1)

		n_i 				array-like
							len (N-1) array holding initial state

		dt 					float
							time-step, needs to have consistent units with
							tloss_discrete and energy_loss_rate

	Returns:
		n_iplusone 			array-like
							array holding the distribution at the next time step
	'''


	# find the lower and upper bin boundaries, bin centres and bin sizes
	E1 = energy_edges[:-1]
	E2 = energy_edges[1:]
	Ecen = 0.5 * (E1 + E2)
	Ebins = E2 - E1

	# check if energy_loss_rate is array like or scalar
	if np.isscalar(energy_loss_rate):
		if energy_loss_rate == None or energy_loss_rate == 0.0:
			loss_term = 0.0
			v3_loss_term = 0.0
		else:
			# allow for constant loss rate
			loss_term = np.ones_like(Ecen) * energy_loss_rate * dt
			v3_loss_term = -np.ones_like(Ecen) * energy_loss_rate * dt
	else:
		# get the j+1/2 and j-1/2 energy loss rates 
		Edot1 = energy_loss_rate[:-1]
		Edot2 = energy_loss_rate[1:]
		loss_term = (dt * Edot1 / Ebins)
		v3_loss_term = -dt * Edot2 / Ebins

	# check if tloss_discrete is array like or scalar
	if np.isscalar(tloss_discrete):
		if tloss_discrete == None or tloss_discrete == 0.0:
			discrete_term = 0.0
		else:
			# allow for single timescale
			discrete_term = (dt / tloss_discrete)
	else:
		discrete_term = (dt / tloss_discrete)

	# set up the four terms to pass to the TDMA solver 
	a = np.zeros_like(Ecen)
	b = 1.0 + loss_term + discrete_term
	c = v3_loss_term * np.ones_like(Ecen)
	d = n_i + (source * dt)

	# run the TDMA solver 
	n_iplusone = msynchro.tdma.TDMASolver(a, b, c, d)

	return (n_iplusone)


def get_dt(energies, loss_rate_cen, tloss_discrete, source, n_i, relative_threshold=1e-15):
	'''
	Calculate a time step. Units consistent with
	tloss_discrete and energy_loss_rate

	Parameters:
		energy_edges 		array-like
							edges of energy bins, len(N)

		energy_loss_rate 	array-like or float
							dE/dt at edges of energy bins, len (N)

		tloss_discrete 		array-like or float
							tau_loss in each bin, len (N-1)

		source 				array-like
							source term, len (N-1)

		n_i 				array-like
							len (N-1) array holding initial state


	Returns:
		n_iplusone 			array-like
							array holding the distribution at the next time step
	'''
	
	if np.all(n_i == 0.0):
		return (0.0)

	minimum_n_i = np.max(n_i) * relative_threshold
	#print (np.max(n_i), minimum_n_i)
	select = (n_i > minimum_n_i)

	if np.isscalar(loss_rate_cen):
		if loss_rate_cen == 0.0:
			dndt1 = 0.0
		else:
			dndt1 = n_i[select] / energies[select] * loss_rate_cen
	else:
		dndt1 = n_i[select] / energies[select] * loss_rate_cen[select]

	if np.isscalar(tloss_discrete):
		if tloss_discrete == 0.0:
			dndt2 = 0.0
		else:
			dndt2 = n_i[select] / tloss_discrete
	else:
		dndt2 = n_i[select] / tloss_discrete[select]

	dt1 = np.nanmin(n_i[select] / dndt1)
	dt2 = np.nanmin(n_i[select] / dndt2)


	delta_t = 0.4 * np.nanmin( (dt1, dt2) )

	return (delta_t)