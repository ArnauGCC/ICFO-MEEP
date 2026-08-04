import meep as mp
import numpy as np


"""
Computes if the region given has reached a stable (converged) state, 
this is, there are no significant fluctuations of the energy inside the region 
in some time (window)

There must be some fields to reach a converged state
"""
def make_stop_when_converged(center, size,	# Energy object
							 window=10,
							 tolerance=1e-6,
							 min_field=1e-6, print_err=False, err_rate=100):

	energy_history = []
	i = 0

	def stop_when_field_const(sim):
		Ez = sim.get_array(center=center, size=size, component=mp.Ez)
		energy = np.sum(np.abs(Ez)**2)

		energy_history.append(energy)

		if len(energy_history) < 2*window:
			return False

		avg1 = np.mean(energy_history[-2*window:-window], axis=0)
		avg2 = np.mean(energy_history[-window:], axis=0)

		# Do not check convergence before energy has arrived
		if np.max(np.abs(avg2)) < min_field:
			return False

		error = np.mean(np.abs(avg2 - avg1) / (np.abs(avg1) + 1e-20))
		if print_err and i%err_rate==0:	print("ERROR: ", error)
			
		return error < tolerance

	return stop_when_field_const


"""
Writes output to the file output.txt
If file not exists is created in working directory
"""
def write_output(text):
	if mp.am_master():
		with open("output.txt", "a") as file:
			file.write(str(text) + "\n")


"""
This function doesn't have to be called
Computes the maximum efficiency of a given frequency between two values of a parameter in eff_parameter
with known source time and flux 
"""
def find_max_efficiency_(coupler, param, min_s, max_s, n_steps, freq, stage, src_time, src_flux):
	write_output("")
	write_output("STAGE:  " + str(stage))

	eff = []
	steps = np.linspace(min_s, max_s, n_steps)
	write_output("STEPS:")
	write_output(steps)

	for s in steps:
		print("EXECUTING MEEP WITH PARAM = ", s)
		coupler.set_global_param(param, s)
		#eff_params[param] = s
		eff.append(coupler.main(freq, src_time, src_flux))

	write_output("EFF:")
	write_output(eff)

	if n_steps == 1:
		return eff[0], steps[0]

	max_v= max(eff)
	max_i = eff.index(max_v)

	second_v = max(n for n in eff if n != max_v)
	second_i = eff.index(second_v)

	step_size = (max_s - min_s) / n_steps

	# The second max is not a neighbour
	if abs(max_i - second_i) > 1:
		if max_i + 1 < n_steps and abs(max_v - eff[max_i + 1]) < 0.01 and stage > 0:
					return max_v, steps[max_i]
		elif max_i - 1 >= 0 and abs(max_v - eff[max_i - 1]) < 0.01 and stage > 0:
			return max_v, steps[max_i]
		
		return find_max_efficiency_(coupler, param, min_s, max_s, 2*n_steps, freq, stage+1, src_time, src_flux)

	elif max_i == 0:
		return find_max_efficiency_(coupler, param, min_s - n_steps*step_size/2, max_s - n_steps*step_size/2, n_steps, freq, stage+1, src_time, src_flux)

	elif max_i == n_steps-1:
		return find_max_efficiency_(coupler, param, min_s + n_steps*step_size/2, max_s + n_steps*step_size/2, n_steps, freq, stage+1, src_time, src_flux)

	# Differnce is lower than 1%
	if stage > 0 and abs(max_v - second_v) < 0.01:
		return max_v, steps[max_i]

	if max_i == n_steps - 1:
		return find_max_efficiency_(coupler, param, steps[max_i-1], steps[max_i] + step_size, n_steps, freq, stage+1, src_time, src_flux)

	elif max_i == 0:
		return find_max_efficiency_(coupler, param, steps[0] - step_size, steps[1], n_steps, freq, stage+1, src_time, src_flux)
		
	elif max_i - second_i > 0:
		return find_max_efficiency_(coupler, param, steps[second_i], steps[max_i] + step_size,
							n_steps if abs(max_v - second_v) > 0.015 or int(n_steps/2) <= 2 
							else int(n_steps/2), 
							freq, stage+1, src_time, src_flux)
	else:
		return find_max_efficiency_(coupler, param, steps[max_i] - step_size, steps[second_i], 
									n_steps if abs(max_v - second_v) > 0.015 or int(n_steps/2) <= 2 
									else int(n_steps/2), 
									freq, stage+1, src_time, src_flux)
