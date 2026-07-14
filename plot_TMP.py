from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import numpy as np


file_name = "/home/arnaugc/ICFO/MEEP/TMP-dft-off0_fcen0.25_w1.7777777777777777_al45_n2.25_pad0.5.npz"
data = np.load(file_name)

src_SPD = data["src_power"]
wg_SPD = data["wg_power"]
refl_SPD = data["refl_power"]
freqs = data["freqs"]
df = data["df"]

plt.figure()
plt.plot(freqs,src_SPD)
plt.xlabel(r'frequency $f (kHz)$')
plt.ylabel("Fourier transform (SPD)")
plt.title("Source SPD")
plt.grid(True, which="both", alpha=0.3)
plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

plt.figure()
plt.plot(freqs,wg_SPD)
plt.xlabel(r'frequency $f (kHz)$')
plt.ylabel("Fourier transform (SPD)")
plt.title("Waveguide SPD")
plt.grid(True, which="both", alpha=0.3)
plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

plt.figure()
plt.plot(freqs,refl_SPD)
plt.xlabel(r'frequency $f (kHz)$')
plt.ylabel("Fourier transform (SPD)")
plt.title("Reflction SPD")
plt.grid(True, which="both", alpha=0.3)
plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

plt.show()