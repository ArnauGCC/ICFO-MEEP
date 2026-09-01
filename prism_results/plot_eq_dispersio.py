import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import brentq

def f(x, n, d):
    return (
        np.pi * d * np.sqrt(n**2 - x**2)
        - np.arctan(
            np.sqrt((x**2 - 1)/(n**2 - x**2))
        )
    )



f_min = 0.7
f_max = 3.0
f_res = 0.1
f_steps = round((f_max - f_min)/f_res) +1

do_prev_plot = False

n = 1.50                                # index waveguide

freqs = np.linspace(f_min, f_max, f_steps)
wg_factors = [np.float64(0.7045454545454546), np.float64(0.7522727272727273), np.float64(0.7763636363636364), np.float64(0.8045454545454546), np.float64(0.8045454545454546), np.float64(0.8477272727272728), np.float64(0.8763636363636365), np.float64(0.8954545454545455), np.float64(0.9045454545454545), np.float64(0.9045454545454545), np.float64(0.9236363636363636), np.float64(0.9310227272727273), np.float64(0.9524999999999999), np.float64(0.9572727272727273), np.float64(0.9763636363636363), np.float64(0.9763636363636363), np.float64(0.9954545454545454), np.float64(0.9954545454545454), np.float64(0.9854545454545455), np.float64(1.0045454545454546), np.float64(1.0045454545454546), np.float64(1.0045454545454546), np.float64(1.0045454545454546), np.float64(1.0045454545454546)]

wg_widths = wg_factors/(2*freqs*n)
eff_wg_width = [0.7029030933779804, 0.7004010622479051, 0.7050796035058626, 0.7106936385651379, 0.711708101020767, 0.71733140326665, 0.7203518202258516, 0.7196660928141703, 0.7139453022686275, 0.7155084474443585, 0.7142760857853228, 0.7037149237638584, 0.6973642456025235, 0.6929148311538778, 0.6903930805539424, 0.6858767355867705, 0.6807204156015279, 0.6761178986136613, 0.6691303211446368, 0.6638359502025389, 0.6564535716529222, 0.647508085528874, 0.6361776335178501, 0.6295057526726097]

if do_prev_plot:
    f_max = 2.5
    f_steps = round((f_max - f_min)/f_res) +1
    freqs_prev = np.linspace(f_min, f_max, f_steps)

    wg_factors_prev = [np.float64(0.7857142857142857), np.float64(0.7576530612244897), np.float64(0.7755102040816326), np.float64(0.8137755102040817), np.float64(0.8137755102040817), np.float64(0.8714285714285714), np.float64(0.8622448979591837), np.float64(0.9005102040816326), np.float64(0.90625), np.float64(0.9183673469387755), np.float64(0.9183673469387755), np.float64(0.9566326530612245), np.float64(0.9566326530612245), np.float64(0.9566326530612245), np.float64(0.9948979591836735), np.float64(0.9948979591836735), np.float64(0.9964285714285714), np.float64(0.9964285714285714), np.float64(0.9955357142857143)]
    eff_wg_width_prev = [0.6882152669262575, 0.6974634137839019, 0.7014771913905148, 0.7049594411871288, 0.7075875802012354, 0.709520078659464, 0.7107515556515785, 0.7116505204116244, 0.7064925173008624, 0.7057978774299665, 0.7017332109972281, 0.694126457986922, 0.6866706801799823, 0.679922625277651, 0.6762803267971592, 0.6723638441510565, 0.6667557245319612, 0.6611362236561644, 0.6544778701219338]
    wg_widths_prev = wg_factors_prev/(2*freqs_prev*n)


    plt.figure()
    plt.plot(wg_widths_prev, freqs_prev)
    plt.xlabel(r'Waveguide width')
    plt.ylabel(r'Frequency')
    plt.title("Equació de dispersió")
    plt.grid(True, which="both", alpha=0.3)

    plt.figure()
    plt.plot(wg_widths_prev, eff_wg_width_prev)
    plt.xlabel(r'Waveguide width')
    plt.ylabel(r'Efficiency')
    plt.title("Waveguide Efficiency")
    plt.grid(True, which="both", alpha=0.3)

    plt.figure()
    plt.plot(freqs_prev, eff_wg_width_prev)
    plt.xlabel(r'Frequency')
    plt.ylabel(r'Efficiency')
    plt.title("Waveguide Efficiency")
    plt.grid(True, which="both", alpha=0.3)

    plt.show()




plt.figure()
plt.plot(wg_widths, freqs)
plt.xlabel(r'Waveguide width')
plt.ylabel(r'Frequency')
plt.title("Equació de dispersió")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(wg_widths, eff_wg_width)
plt.xlabel(r'Waveguide width')
plt.ylabel(r'Efficiency')
plt.title("Waveguide Efficiency")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(freqs, eff_wg_width)
plt.xlabel(r'Frequency')
plt.ylabel(r'Efficiency')
plt.title("Waveguide Efficiency")
plt.grid(True, which="both", alpha=0.3)

plt.show()


eff_pad_wg = [0.1262279008396752, 0.14357630364251836, 0.14249090803693473, 0.15504475336603185, 0.15533905500050416, 0.17835075506292258, 0.18302295326674037, 0.18225209530753983, 0.1849700021750237, 0.1904163886201993, 0.19938226163548484, 0.20782511561958758, 0.21819503961304942, 0.21705976454794468, 0.2148379272276824, 0.2154668058175726, 0.21966784877045087, 0.21827774947324852, 0.2226166651340914, 0.22117428800286626, 0.22284505876498695, 0.22632136014070403, 0.2209122238966025, 0.22325688647184694]
pad_wg = [np.float64(0.42363636363636364), np.float64(0.36363636363636365), np.float64(0.34272727272727277), np.float64(0.2854545454545455), np.float64(0.2572727272727273), np.float64(0.24938579545454548), np.float64(0.22276859504132235), np.float64(0.20454545454545459), np.float64(0.20454545454545459), np.float64(0.18545454545454548), np.float64(0.18545454545454548), np.float64(0.16636363636363638), np.float64(0.1572727272727273), np.float64(0.14902892561983472), np.float64(0.13818181818181818), np.float64(0.1525), np.float64(0.14272727272727273), np.float64(0.12363636363636366), np.float64(0.12363636363636366), np.float64(0.12450413223140497), np.float64(0.11912396694214877), np.float64(0.11365702479338843), np.float64(0.12363636363636366), np.float64(0.10454545454545455)]

plt.figure()
plt.plot(pad_wg, freqs)
plt.xlabel(r'Pad prism - waveguide')
plt.ylabel(r'Frequency')
plt.title("Max eff of waveguide power")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(pad_wg, eff_pad_wg, marker='o', linestyle='None')
plt.xlabel(r'Pad prism - waveguide')
plt.ylabel(r'Efficiency')
plt.title("Waveguide Power Efficiency")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(freqs, eff_pad_wg)
plt.xlabel(r'Frequency')
plt.ylabel(r'Efficiency')
plt.title("Waveguide Power Efficiency")
plt.grid(True, which="both", alpha=0.3)

plt.show()

f_max = 2.8
f_steps = round((f_max - f_min)/f_res) +1
freqs = np.linspace(f_min, f_max, f_steps)

n = 1.5
n_p = n + 0.3
alpha = math.radians(45)

thetas_teo = []

for w in wg_widths[0:-2]:
    n_eff = brentq(f, 1.01, n-0.01, args=(n, w))
    theta_inc = math.asin( math.sin(math.asin(n_eff/n_p) - alpha) * n_p) + alpha

    thetas_teo.append(math.degrees(theta_inc))


eff_off_deg = [np.float64(0.14642618988402584), np.float64(0.16529916324760222), np.float64(0.16421937203871956), np.float64(0.17383908310210705), np.float64(0.1781161967516983), np.float64(0.19506078318882566), np.float64(0.19537600262507623), np.float64(0.19088468227754118), np.float64(0.19636410622106712), np.float64(0.19835222357310572), np.float64(0.21862520012315795), np.float64(0.2208027528573794), np.float64(0.23329576840180066), np.float64(0.22817034070673806), np.float64(0.22781066840028008), np.float64(0.21952860840736144), np.float64(0.22908739299616193), np.float64(0.2268568507816542), np.float64(0.22416664025401345), np.float64(0.23106457480425205), np.float64(0.22798977234315562), np.float64(0.22996461565341952)]
off_deg = [np.float64(3.04545455), np.float64(2.95454545), np.float64(5.04545455), np.float64(4.95454545), np.float64(7.04545455), np.float64(6.09090909), np.float64(7.04545455), np.float64(7.90909091), np.float64(9.04545455), np.float64(9.6020057), np.float64(9.60621236), np.float64(9.6039663), np.float64(9.60845528), np.float64(10.94640562), np.float64(10.94294029), np.float64(12.2844333), np.float64(10.93336416), np.float64(12.21915326), np.float64(13.00365986), np.float64(12.27644804), np.float64(13.38564468), np.float64(13.87262521)]

thetas = np.add(thetas_teo, off_deg)


plt.figure()
plt.plot(freqs, eff_off_deg)
plt.xlabel(r'Frequency')
plt.ylabel(r'Efficiency')
plt.title("Waveguide Power Efficiency - Incident angle ajusted")
plt.grid(True, which="both", alpha=0.3)

plt.show()


"""
max_effs = []
max_pads = []
for i in range(0, len(freqs)):
    if eff_pad_wg[i] > eff_pad_coupled[i]:
        max_effs.append(eff_pad_wg[i])
        max_pads.append(pad_wg[i])
    else:
        max_effs.append(eff_pad_coupled[i])
        max_pads.append(pad_coupled[i])

print(str(max_effs))
print("\n \n")
print(str(max_pads))"""