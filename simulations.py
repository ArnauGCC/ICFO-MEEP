from grating_coupling import GratingCoupler
from packages.figures import Layer
from packages.utils import write_output
import meep as mp

gc = GratingCoupler()
gc.res_factor = 15
gc.h5_file_transistent = False
gc.compute_only_1st_mode = True
gc.n_wg = 3.47


vars = gc.get_variables()           # ['theta_deg', 'n_cells', 'grat_period_factor', 'wg_width_factor', 'grat_height_factor', 'grat_duty_cycle']
gc.set_global_param(vars[0], 10)
gc.set_global_param(vars[1], 12)
gc.set_global_param(vars[2], 0.5)
gc.set_global_param(vars[3], 1.5)
gc.set_global_param(vars[4], 0.5)
gc.set_global_param(vars[5],0.5)


write_output("")
write_output("SIM N_WG=3.47 ONLY 1st MODE theta = 10º N_CELLS = 12 [Only WG]")
freqs, effs = gc.get_eff_in_freq_range(0.8, 2, f_res=0.05)
write_output("")


write_output("")
write_output("SIM N_WG=3.47 ONLY 1st MODE theta = 10º [Only WG, Main lobe]")
max_eff = max(effs)
eff_i = effs.index(max_eff)
max_freq = freqs[eff_i]

freqs, effs = gc.get_eff_in_freq_range(max_freq - 0.05, max_freq + 0.05, f_res=0.002)
write_output("")



gc.set_global_param(vars[0], 25)

write_output("")
write_output("SIM N_WG=3.47 ONLY 1st MODE theta = 25º N_CELLS = 12 [Only WG]")
freqs, effs = gc.get_eff_in_freq_range(0.8, 2, f_res=0.05)
write_output("")


write_output("")
write_output("SIM N_WG=3.47 ONLY 1st MODE theta = 25º [Only WG, Main lobe]")
max_eff = max(effs)
eff_i = effs.index(max_eff)
max_freq = freqs[eff_i]

freqs, effs = gc.get_eff_in_freq_range(max_freq - 0.05, max_freq + 0.05, f_res=0.002)
write_output("")


gc.set_global_param(vars[0], 10)
gc.bottom_layers.append(Layer(2.27, 1.44))
gc.bottom_layers.append(Layer(mp.inf, 3.47))

write_output("")
write_output("###########################################################")
write_output("###########################################################")
write_output("")

write_output("")
write_output("SIM N_WG=3.47 ONLY 1st MODE N_EXT = 1 theta = 10º N_CELLS = 12 [Layer Stack]")
freqs, effs = gc.get_eff_in_freq_range(0.8, 2, f_res=0.05)
write_output("")


write_output("")
write_output("SIM N_WG=3.47 ONLY 1st MODE N_EXT = 1 theta = 10º N_CELLS = 12 [Layer Stack]")
max_eff = max(effs)
eff_i = effs.index(max_eff)
max_freq = freqs[eff_i]

freqs, effs = gc.get_eff_in_freq_range(max_freq - 0.05, max_freq + 0.05, f_res=0.002)
write_output("")



gc.n_default = 1.46
write_output("")
write_output("###########################################################")
write_output("###########################################################")
write_output("")

write_output("")
write_output("SIM N_WG=3.47 ONLY 1st MODE N_EXT = 1.46 theta = 10º N_CELLS = 12 [Layer Stack]")
freqs, effs = gc.get_eff_in_freq_range(0.8, 2, f_res=0.05)
write_output("")


write_output("")
write_output("SIM N_WG=3.47 ONLY 1st MODE N_EXT = 1.46 theta = 10º N_CELLS = 12 [Layer Stack]")
max_eff = max(effs)
eff_i = effs.index(max_eff)
max_freq = freqs[eff_i]

freqs, effs = gc.get_eff_in_freq_range(max_freq - 0.05, max_freq + 0.05, f_res=0.002)
write_output("")



gc.set_global_param(vars[0], 25)

write_output("")
write_output("SIM N_WG=3.47 ONLY 1st MODE N_EXT = 1.46 theta = 25º N_CELLS = 12 [Layer Stack]")
freqs, effs = gc.get_eff_in_freq_range(0.8, 2, f_res=0.05)
write_output("")


write_output("")
write_output("SIM N_WG=3.47 ONLY 1st MODE N_EXT = 1.46 theta = 25º N_CELLS = 12 [Layer Stack]")
max_eff = max(effs)
eff_i = effs.index(max_eff)
max_freq = freqs[eff_i]

freqs, effs = gc.get_eff_in_freq_range(max_freq - 0.05, max_freq + 0.05, f_res=0.002)
write_output("")