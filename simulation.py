from grating_coupling import GratingCoupler
from packages.figures import Layer
from packages.utils import write_output
import meep as mp
#import numpy as np

gc = GratingCoupler()
vars = gc.get_variables()           # ['theta_deg', 'n_cells', 'grat_period_factor', 'wg_width_factor', 'grat_height_factor', 'grat_duty_cycle']

"""
gc.set_global_param(vars[0], 37.325)
gc.set_global_param(vars[1], 15)
gc.set_global_param(vars[2], 1)
gc.set_global_param(vars[3], 1.457273)
gc.set_global_param(vars[4], 0.5354)
gc.set_global_param(vars[5],0.4906)
"""

"""
gc.n_wg = 2
gc.res_factor = 30

vals = [37.325, np.int64(14), np.float64(0.9984553850156087), np.float64(1.4572727272727273), np.float64(0.5354371378180902), np.float64(0.4906462585034014)]
for i in range(len(vals)):
    gc.set_global_param(vars[i], vals[i])
"""

gc.set_global_param(vars[0], 12)
gc.set_global_param(vars[1], 8)
gc.set_global_param(vars[2], 0.5)
gc.set_global_param(vars[3], 1.5)
gc.set_global_param(vars[4], 0.5)
gc.set_global_param(vars[5],0.5)

gc.n_wg = 3.47
gc.src_freq = 0.834

#write_output("")
#write_output("SIM WITH ONLY THE WG n=3.47")
#eff, n_cells = gc.find_max_efficiency(vars[1], 5, 40, 8, 1)

#gc.set_global_param(vars[1], n_cells)
#write_output("EXECUTING SIMULATION WITH N_CELLS = " + str(n_cells))


#freqs, effs = gc.get_eff_in_freq_range(0.8, 0.95, f_res=0.002)
#write_output(freqs)
#write_output(effs)

gc.res_factor = 25
gc.h5_file_transistent = True
gc.compute_only_1st_mode = True


gc.set_global_param(vars[0], 25)
gc.set_global_param(vars[1], 12)
gc.set_global_param(vars[2], 0.5)
gc.set_global_param(vars[3], 1)
gc.set_global_param(vars[4], 0.5)
gc.set_global_param(vars[5],0.5)

#gc.n_default = 1.46
#gc.bottom_layers.append(Layer(2.27, 1.44))
#gc.bottom_layers.append(Layer(mp.inf, 3.47))

gc.manual_res = int(25 * 3.47 * 0.95)
gc.set_manual_resolution = True
gc.h5_frames_per_wlength = 20
gc.n_wlengths_power_measure = 5
gc.width_sim_scale = 1.8

gc.src_freq = 0.914
print(gc.main())
