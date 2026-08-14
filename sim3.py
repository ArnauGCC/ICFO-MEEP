from grating_coupling import GratingCoupler
from packages.figures import Layer
from packages.utils import write_output
import numpy as np
import meep as mp


gc = GratingCoupler()
vars = gc.get_variables()           # ['theta_deg', 'n_cells', 'grat_period_factor', 'wg_width_factor', 'grat_height_factor', 'grat_duty_cycle']

gc.set_global_param(vars[0], 12)
gc.set_global_param(vars[1], 6)
gc.set_global_param(vars[2], 0.5)
gc.set_global_param(vars[3], 1.5)
gc.set_global_param(vars[4], 0.5)
gc.set_global_param(vars[5],0.5)

# Optimal values for src_freq = freq (= 1)
#vals = [np.float64(10), np.int64(15), np.float64(0.5215107960457857), np.float64(1.3727272727272728), np.float64(0.430952380952381), np.float64(0.447215536902577)]
#for i in range(len(vals)):
#    gc.set_global_param(vars[i], vals[i])

gc.n_wg = 3.47
gc.res_factor = 20
gc.bottom_layers.append(Layer(2.27, 1.44))
gc.bottom_layers.append(Layer(mp.inf, 3.47))
gc.n_default = 1.46

gc.compute_only_1st_mode = True


"""
write_output("")
write_output("SIM N_EXT = 1.46 ONLY 1st MODE theta = 12º")

freqs, effs = gc.get_eff_in_freq_range(0.8, 1, f_res=0.002)
write_output(freqs)
write_output(effs)

"""
gc.set_global_param(vars[0], 10)
#write_output("")
#write_output("SIM N_EXT = 1.46 ONLY 1st MODE theta = 10º")
gc.src_freq = 0.852
gc.set_global_param(vars[4], 0.1)
gc.optimize_all_parameters(1, 10, 6, grat_height_factor=0.1)

#gc.h5_file_transistent = True
#gc.src_freq = 0.852
#print(gc.main())

