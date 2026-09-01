import matplotlib.pyplot as plt
import h5py


DO_PLOTS_AT_END = True

plot_only_waveguide = True
plot_stack_n_ext_1_00 = True
plot_stack_n_ext_1_46 = True
plot_sims_theta_10 = True
plot_sims_theta_25 = True

h5_sims = {
    "only_wg_10":       "/freq_range_only_grating_theta10",
    "only_wg_25":       "/freq_range_only_grating_theta25",
    "stack_1_10":       "/freq_range_layer_stack_theta10_next1",
    "stack_1_25":       "/freq_range_layer_stack_theta25_next1",
    "stack_1-46_10":    "/freq_range_layer_stack_theta10_next1-46",
    "stack_1-46_25":    "/freq_range_layer_stack_theta25_next1-46",
}


def do_plot_(file, simulation, title, large=False):
    try:
        if large:
            freqs = file[simulation]["freqs_large"][:]
            effs = file[simulation]["fund_mode_eff_large"][:]
        else:
            freqs = file[simulation]["freqs_main_lobe"][:]
            effs = file[simulation]["fund_mode_eff_main_lobe"][:]        
                    
        plt.figure()
        plt.plot(freqs, effs)
        plt.xlabel(r'Frequency')
        plt.ylabel(r'Efficiency (Fundamental Mode Power)')
        plt.title(title)
        plt.grid(True, which="both", alpha=0.3)
    except Exception as e:
        print("-------------DATA NOT FOUND-------------")
        print(e)
        print("Corresponding with plot:")
        print(title)
        print("----------------------------------------")


def do_plot(file, sim, text):
    if plot_sims_theta_10:
        ### Theta = 10º ###
        simulation = h5_sims[sim+"_10"]
        if not DO_PLOTS_AT_END:
            do_plot_(file, simulation, fr"Grating Coupler N_WG=3.47 $\theta$ = 10º N_CELLS = 12 [{text}]", True)
    
        do_plot_(file, simulation, fr"Grating Coupler Main Lobe N_WG=3.47 $\theta$ = 10º N_CELLS = 12 [{text}]")
       
    if plot_sims_theta_25:
        ### Theta = 25º ###
        simulation = h5_sims[sim+"_25"]
        if not DO_PLOTS_AT_END:
            do_plot_(file, simulation, fr"Grating Coupler N_WG=3.47 $\theta$ = 25º N_CELLS = 12 [{text}]", True)

        do_plot_(file, simulation, fr"Grating Coupler Main Lobe N_WG=3.47 $\theta$ = 25º N_CELLS = 12 [{text}]")

    if not DO_PLOTS_AT_END:     plt.show()


def main():
    with h5py.File("grating_coupler_results.h5", "r") as f: 
        if plot_only_waveguide:         do_plot(f, "only_wg",       "Only WG")
        if plot_stack_n_ext_1_00:       do_plot(f, "stack_1",       "Layer Stack, n_ext=1")
        if plot_stack_n_ext_1_46:       do_plot(f, "stack_1-46",    "Layer Stack, n_ext=1.46")
        if DO_PLOTS_AT_END:             plt.show()


if __name__ == "__main__":
    main()