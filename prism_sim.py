from prism_coupling import PrismCoupler
from packages.utils import *

pc = PrismCoupler()
vars = pc.get_variables()           # ['wg_width', 'pad', 'offs_deg']

pc.set_global_param(vars[0], 0.55)
pc.set_global_param(vars[1], 0.46)
pc.set_global_param(vars[2], 0)



pc.prism_length = 25
pc.n = 2
pc.compute_power = False
pc.compute_modes_coeff = True

#pc.set_same_pad = True
pc.run_meep = True
pc.do_plots = True
pc.main(0.75)

"""
write_output("WG_HEIGHT: 0.55, INDEX_DIFF: 0.275")
iss = compute_freqs(0.435, 0.525, f_res=0.005)
for i in iss:
    pc.set_global_param(vars[1], i)
    write_output("PAD: ", i, "EFF: ", pc.main(0.75))
"""
"""
wgws = compute_freqs(0.1, 0.35, f_res=0.0125)

effs = []
for wgw in wgws:
    write_output("Executing with wg_width = ", wgw)
    pc.set_global_param(vars[0], wgw)
    eff = pc.get_eff_in_freq_range(0.7, 3, f_res=0.05, set_same_resolution=True, set_same_pad=True)

write_output("\n\nFINAL RESULT:")    
write_output(effs)
"""