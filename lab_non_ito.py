import numpy as np
from packages.utils import *
import matplotlib.pyplot as plt
import h5py


def do_plot(X, Y, matrix, title, xlabel='Duty Cycle', ylabel='Depth Factor'):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X, Y, matrix, cmap='viridis')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel('Efficiency')
    plt.title(f"lab grating without ITO with n_cells=20 grat_height (PE layer) = 0.7 {title}")


def main():
    sim = "/lab_non_ito_chng_depth-duty"
    with h5py.File("grating_coupler_results.h5", "r") as f: 
        X, Y = np.meshgrid(f[sim]['x_values'][:], f[sim]['y_values'][:])

        do_plot(X, Y, f[sim]["total_power_R"], "TOTAL POWER Right")
        do_plot(X, Y, f[sim]["mode0R"], "MODE 0R")
        do_plot(X, Y, f[sim]["mode1R"], "MODE 1R")
        do_plot(X, Y, f[sim]["mode2R"], "MODE 2R")
        do_plot(X, Y, f[sim]["total_power_L"], "TOTAL POWER Left")

    plt.show()


if __name__ == "__main__":
    main()