from grating_coupling import GratingCoupler
from packages.figures import Layer
from packages.utils import write_output
import meep as mp
import numpy as np

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

gc.n_wg = 2
gc.res_factor = 40

vals = [37.325, np.int64(14), np.float64(0.9984553850156087), np.float64(1.4572727272727273), np.float64(0.5354371378180902), np.float64(0.4906462585034014)]
for i in range(len(vals)):
    gc.set_global_param(vars[i], vals[i])


#gc.bottom_layers.append(Layer(2.27, 1.5))
#gc.bottom_layers.append(Layer(mp.inf, 2))
#gc.run_meep = False
#gc.do_plots = True

gc.h5_file_transistent = True
print(gc.main())
#gc.find_max_efficiency(vars[1], 5, 40, 8, 1)

#gc.optimize_all_parameters(1, theta_deg=37.325, n_cells=GratingCoupler.Interval(10.625, 20, 8))
#
#write_output("###########################################################################")
#write_output("###########################################################################")
#write_output("")
#write_output("")




