from grating_coupling import GratingCoupler
from packages.figures import Layer
from packages.utils import write_output
import meep as mp
#import numpy as np

gc = GratingCoupler()
vars = gc.get_variables()           # ['theta_deg', 'n_cells', 'grat_period_factor', 'wg_width_factor', 'grat_height_factor', 'grat_duty_cycle']


#gc.set_global_param(vars[0], 10)
gc.set_global_param(vars[1], 12)
gc.set_global_param(vars[2], 0.5)
gc.set_global_param(vars[3], 1.5)
gc.set_global_param(vars[4], 0.5)
gc.set_global_param(vars[5],0.5)


gc.n_wg = 3.47
gc.set_global_param(vars[0], 25)
gc.set_global_param(vars[3], 1.1)


gc.src_freq = 1.225
gc.res_factor = 20


gc.bottom_layers.append(Layer(2.27, 1.44))
gc.bottom_layers.append(Layer(mp.inf, 3.47))

#gc.manual_res = int(25 * 3.47 * 0.975)
#gc.set_manual_resolution = True
#gc.h5_frames_per_wlength = 20
#gc.n_wlengths_power_measure = 5
gc.width_sim_scale = 1.5

#gc.compute_src_time = False
#gc.src_time = 290/4

gc.compute_eff_by_modes = True
gc.n_modes_to_compute = 5
mp.verbosity(1)

#gc.h5_file_transistent = True

gc.do_plots = True
write_output(gc.main())

