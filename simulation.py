from grating_coupling import GratingCoupler
from packages.figures import Layer
from packages.utils import *
import meep as mp
#import numpy as np

gc = GratingCoupler()
vars = gc.get_variables()           # ['theta_deg', 'n_cells', 'grat_period', 'grat_height', 'grat_depth_factor', 'grat_duty_cycle']


gc.set_global_param(vars[0], 45)
gc.set_global_param(vars[1], 20)
gc.set_global_param(vars[2], 0.41)
gc.set_global_param(vars[3], 0.7)
gc.set_global_param(vars[4], 0.5)
gc.set_global_param(vars[5],0.5)

gc.src_freq = 1/0.6321

gc.n_grat = 2.30245
#gc.bottom_layers.append(Layer(height=0.7, index=2.30245))
gc.bottom_layers.append(Layer(height=mp.inf, index=2.2024))


#gc.manual_res = int(25 * 3.47 * 0.975)
#gc.set_manual_resolution = True
#gc.h5_frames_per_wlength = 20
#gc.n_wlengths_power_measure = 5
gc.width_sim_scale = 1.4

gc.compute_eff_by_modes = True
gc.n_modes_to_compute = 8


#gc.run_meep = False
#gc.h5_frames_per_wlength=20
#gc.h5_file_transistent = True
#gc.do_plots = True
#gc.show_region_converged_state = True

#write_output(gc.main())


"""
write_output("\n\n\n*******************************************************************")
write_output("*******************************************************************")
write_output("SIM LAB Grating n_cells = 20 changing depth - duty cycle")

parameters = gc.EffParams(theta_deg=45, 
                          n_cells=20, 
                          grat_period=0.41, 
                          grat_height=0.7,          #gc.Interval(50e-3, 1.5, res_s=50e-3) 
                          grat_depth_factor=gc.Interval(0, 1, res_s=0.025), 
                          grat_duty_cycle=gc.Interval(0, 1, res_s=0.025))

result, steps = gc.scan_all_parameters(gc.src_freq, parameters)

write_output("\n\n")
write_output(result.dims)
write_output(result.data)
write_output("\n\n")
write_output(steps.data)

"""



write_output("\n\n\n*******************************************************************")
write_output("******************************************************************")
write_output("SIM LAB Grating n_cells = 20 changing height - depth")

parameters = gc.EffParams(theta_deg=45, 
                          n_cells=20, 
                          grat_period=0.41, 
                          grat_height=gc.Interval(0.2, 1.2, res_s=0.025), 
                          grat_depth_factor=gc.Interval(0, 1, res_s=0.025), 
                          grat_duty_cycle=0.5)

result, steps = gc.scan_all_parameters(gc.src_freq, parameters)

write_output("\n\n")
write_output(result.dims)
write_output(result.data)
write_output("\n\n")
write_output(steps.data)


write_output("\n\n\n*******************************************************************")
write_output("******************************************************************")
write_output("SIM LAB Grating n_cells = 20 changing period - height")

parameters = gc.EffParams(theta_deg=45, 
                          n_cells=20, 
                          grat_period=gc.Interval(0.3, 1.3, res_s=0.025), 
                          grat_height=gc.Interval(0.2, 1.2, res_s=0.025), 
                          grat_depth_factor=0.5, 
                          grat_duty_cycle=0.5)

result, steps = gc.scan_all_parameters(gc.src_freq, parameters)

write_output("\n\n")
write_output(result.dims)
write_output(result.data)
write_output("\n\n")
write_output(steps.data)