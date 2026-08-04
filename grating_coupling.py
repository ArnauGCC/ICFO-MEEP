import meep as mp
import math
import matplotlib.pyplot as plt
import numpy as np
from packages.figures import *
from packages.utils import *

class GratingCoupler:
	
	""" EXECUTION (Parameters for main() function) """
	run_meep = True							# Runs the simulation, False to check geometry
	do_plots = False						# True:	Plots the geometry (and fields if run_meep)
	end_src = False							# Computes efficencies turning on and off the source
	compute_eff = True						# Computes the efficiency of the coupled waveguide
	compute_src_time = True					# True: src_time is computed in main(), else: src_time is set manual
	compute_src_power = True				# True: main() computes the source power, else: src_power is set manual
	show_region_converged_state = False		# Shows the region used to determinate a converged state (if do_plots)


	""" PARAMETERS """
	freq = 1
	eff_params = {							# Variables that can be used find max efficiency 
		"theta_deg":            16, 		# Degrees of inclination from normal incidence to the waveguide
		"grat_period_factor":   1, 			# The grating period is gp = eff_params['grat_period_factor'] / freq
		"wg_width_factor":      1,			# The waveguide widyh is wg_width = eff_params["wg_width_factor"] / (2*freq)
		"grat_height_factor":   1/2,		# The grating height is gh = eff_params['grat_height_factor'] * wg_width
		"grat_duty_cycle":		0.5,		# Sets the grating duty cycle
	}
	n_wg = 1.5								# waveguide refraction index
	n_wlengths_power_measure = 30			# Number of wavelengths spent to measure the efficiency: self.compute_power or self.compute_modes_coeff
	n_cells = 40							# Number of cells of the grating

	res_factor = 20							# Number of pixels for wavelength in the highest refraction index 
	src_time = -1							# if not compute_src_time, is assumed a converged state after src_time time units to compute efficiency
	src_power = -1							# if not compute_src_power, then must specify src_flux to compute efficiency
	h5_name = ""

	def compute_initial_parameters(self, freq, gp, wg_width, gh, gdc):
		pad_src_wg = 5/freq
		pad_inf = 3/freq
		dpml = 1/freq

		sx = int(self.n_cells*gp*2.25)
		sy = int(wg_width+2*dpml+pad_src_wg+pad_inf)

		wg_y = -sy/2 + dpml + pad_inf + wg_width/2
		geometry = create_h_grating(gp, gh, gdc, self.n_cells, wg_width, self.n_wg, mp.Vector3(-sx/4, wg_y))
		
		src_size = mp.Vector3(x = sx/2)
		src_center = mp.Vector3(-sx/2 + dpml + src_size.x/2, sy/2 - dpml)
		theta = math.radians(90 - self.eff_params["theta_deg"])
		k = mp.Vector3(1).rotate(mp.Vector3(0, 0, -1), theta)
		beam_w0 = sx/2

		src = [mp.GaussianBeam2DSource(
				src=mp.ContinuousSource(freq),
				center=src_center,
				size=src_size,
				beam_x0=sy*k/4,                 # relatiu al centre de la font
				beam_kdir=k,
				beam_w0=beam_w0,                # beam waist
				beam_E0=mp.Vector3(0, 0, 1),
				)]

		return dpml, sx, sy, wg_y, geometry, src_size, src_center, k, beam_w0, src


	def main(self, freq=None, src_time=None, src_power=None):
		if freq is None:
			freq = self.freq

		if not self.compute_src_power and src_power is None:
			src_power = self.src_power

		if not self.compute_src_time and src_time is None:
			src_time = self.src_time


		wg_width = 	self.eff_params["wg_width_factor"] / (2*freq)
		gp = 		self.eff_params['grat_period_factor'] / freq
		gdc = 		self.eff_params['grat_duty_cycle']
		gh = 		self.eff_params['grat_height_factor'] * wg_width

		df = 0
		nfreq =  1

		dpml, sx, sy, wg_y, geometry, src_size, src_center, k, beam_w0, src = self.compute_initial_parameters(freq, gp, wg_width, gh, gdc)
		resolution = int(self.res_factor * self.n_wg * freq)

		src_fr = mp.FluxRegion(center=src_center - mp.Vector3(y=1), size=src_size*1.2)

		"""Compute source power to get efficiency (only in stable state)"""
		if self.compute_src_power and not self.end_src:
			sim_src = mp.Simulation(cell_size=mp.Vector3(sx, sy),
									sources=src,
									resolution=resolution,
									boundary_layers=[mp.PML(dpml)]
									)

			sim_src.run(until=10)

			src_region = sim_src.add_flux(freq, df, nfreq, src_fr)
			sim_src.run(until= self.n_wlengths_power_measure/freq)
			src_power = mp.get_fluxes(src_region)[0]


		sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
							geometry=geometry,
							sources=src,
							resolution=resolution,
							boundary_layers=[mp.PML(dpml)]
							)

		conv_region_size = mp.Vector3(sx/4, wg_width)
		conv_region_center = mp.Vector3((sx - conv_region_size.x)/2, wg_y)

		wg_fr = mp.FluxRegion(center=mp.Vector3(sx/2 - 2*dpml, wg_y), size=mp.Vector3(y=wg_width*2))
		
		if self.run_meep:
			stop_cond = make_stop_when_converged(center=conv_region_center, size=conv_region_size,
												window=int(2*10*resolution/(freq*self.n_wg)), print_err=False, 
												tolerance=1e-5, err_rate=5000)

			if self.end_src and self.compute_eff:
				wg_region = sim.add_flux(freq, df, nfreq, wg_fr)

			sim.run(
#					mp.at_beginning(mp.with_prefix(f"{self.h5_name}-", mp.output_epsilon)),
#					mp.at_every(1, mp.to_appended(f"{self.h5_name}-ez", mp.output_efield_z)),
					until=stop_cond if self.compute_src_time else src_time)


			if self.end_src:
				time = sim.meep_time()
				sim.change_sources([])

			if self.compute_eff:
				
				if not self.end_src:
					wg_region = sim.add_flux(freq, df, nfreq, wg_fr)

				if self.h5_name != "":	self.h5_name = self.h5_name + '-'
				sim.run(
	                    mp.at_beginning(mp.with_prefix(f"{self.h5_name}", mp.output_epsilon)),
	                    mp.at_every(1, mp.to_appended(f"{self.h5_name}ez", mp.output_efield_z)),
						until= mp.stop_when_energy_decayed(dt=int(5/freq), decay_by=1e-4) if self.end_src 
						else self.n_wlengths_power_measure/freq)

				"""Compute source power to get efficiency (turning on and off)"""
				if self.compute_src_power and self.end_src:
					src = [mp.GaussianBeam2DSource(
								src=mp.ContinuousSource(freq, end_time=time),
								center=src_center,
								size=src_size,
								beam_x0=sy*k/4,                 # relatiu al centre de la font
								beam_kdir=k,
								beam_w0=beam_w0,                # beam waist
								beam_E0=mp.Vector3(0, 0, 1),
								)]

					sim_src = mp.Simulation(cell_size=mp.Vector3(sx, sy),
													sources=src,
													resolution=resolution,
													boundary_layers=[mp.PML(dpml)]
													)
					
					src_region = sim_src.add_flux(freq, df, nfreq, src_fr)
					sim_src.run(
#						mp.at_every(1, mp.to_appended(f"src-ez", mp.output_efield_z)),
						until=mp.stop_when_energy_decayed(dt=int(5/freq), decay_by=1e-2))
					src_power = mp.get_fluxes(src_region)[0]

				eff = -mp.get_fluxes(wg_region)[0]/src_power
				return eff

				

		if self.do_plots:
			if self.show_region_converged_state:
				fr = mp.FluxRegion(center=conv_region_center, size=conv_region_size)
				sim.add_flux(freq, 0, 1, fr)

			sim.plot2D(fields=mp.Ez,
						eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none'},
						field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
						boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
						output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)))

			plt.show()


	def set_global_param(self, param, value):
		self.eff_params[param] = value


	def get_variables(self):
		return list(self.eff_params.keys())


	def find_max_efficiency(self, param, min_s, max_s, n_steps, freq):
		"""
		Computes the maximum efficiency of a given frequency (global parameter) between two values of a 
		variable in self.eff_parameter (to get the list of variables use: get_variables())
		Partial results are printed in the created file output.txt
		Returns the maximum efficieny with the corresponding value of the parameter  
		"""

		if param not in self.eff_params:	raise Exception(f"The variable ({param}) must be one of the "\
			"self.eff_params keys: {self.get_variables()}. To get the list use: get_variables()")

		self.run_meep = True
		self.do_plots = False
		self.end_src = False
		self.compute_eff = True
		self.compute_src_power = param == "theta_deg"
		self.compute_src_time = True
		self.show_region_converged_state = False

		src_power = None
		"""Compute Source Power"""
		if not self.compute_src_power:
			wg_width = 	self.eff_params["wg_width_factor"] / (2*freq)
			gp = 		self.eff_params['grat_period_factor'] / freq
			gdc = 		self.eff_params['grat_duty_cycle']
			gh = 		self.eff_params['grat_height_factor'] * wg_width
			resolution = int(self.res_factor * self.n_wg * freq)
	
			dpml, sx, sy, wg_y, geometry, src_size, src_center, k, beam_w0, src = self.compute_initial_parameters(freq, gp, wg_width, gh, gdc)

			sim_src = mp.Simulation(cell_size=mp.Vector3(sx, sy),
									sources=src,
									resolution=resolution,
									boundary_layers=[mp.PML(dpml)]
									)

			sim_src.run(until=10)
			src_fr = mp.FluxRegion(center=src_center - mp.Vector3(y=1), size=src_size*1.2)

			src_region = sim_src.add_flux(freq, 0, 1, src_fr)
			sim_src.run(until= self.n_wlengths_power_measure/freq)
			src_power = mp.get_fluxes(src_region)[0]

		time = None
		if param not in ["theta_deg", "grat_period_factor", "grat_height_factor", "grat_duty_cycle"]:
			self.compute_src_time = False

			"""Time for converged state is maximum for larger wg_widths"""
			if param == "wg_width_factor":
				self.set_global_param(param, max_s)

			sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
										geometry=geometry,
										sources=src,
										resolution=resolution,
										boundary_layers=[mp.PML(dpml)]
										)

			conv_region_size = mp.Vector3(sx/4, wg_width)
			conv_region_center = mp.Vector3((sx - conv_region_size.x)/2, wg_y)

			stop_cond = make_stop_when_converged(center=conv_region_center, size=conv_region_size,
													window=int(2*10*resolution/(freq*self.n_wg)), tolerance=1e-5)

			sim.run(until=stop_cond)
			time = sim.meep_time()

		return find_max_efficiency_(self, param, min_s, max_s, n_steps, freq, 0, time, src_power)


if __name__ == "__main__":

	"""
	grat_period_factor: 0.5, 2
	wg_width_factor:	0.5, 1.5
	grat_height_factor:	0.2, 0.8
	grat_duty_cycle:	0.2, 0.8
	"""



	"""
	ef, fact = find_max_efficiency_("theta_deg", 16, 16.5, 2, freq, 0, -1, -1)

    write_output(str(ef))
    write_output(str(fact))
"""
#    eff_params['theta_deg'] = 16.12142857142857
#    print(main())