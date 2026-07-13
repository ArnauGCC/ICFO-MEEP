from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import numpy as np


file_name = "/home/arnaugc/ICFO/MEEP/prismes/dft-gam0.65_al20_n2.0_pad0.5_w1.npz"
data = np.load(file_name)

freqs = data["freqs"]
power = data["power"]

plt.plot(freqs, power)
plt.gca().xaxis.set_major_locator(MultipleLocator(0.75*0.8*0.10))
plt.xlabel(r'frequency $f (kHz)$')
plt.ylabel("Fourier transform (PDF)")
plt.grid(True, which="both", alpha=0.3)
plt.show()