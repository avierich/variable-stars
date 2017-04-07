# periodFit.py
#
# David Lister
# March, 2017
#

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
from phys_toolkit import *

def fftMaxima(x, y, minDif = 0.05):
    if not evenSpacing(x):
        x, y = spaceEvenly(x, y)
    fft = np.fft.rfft(y)
    fft_mag = np.log(abs(fft))
    fft_phi = np.angle(fft)
    fft_phi_dict = {fft_mag[i]:fft_phi[i] for i in range(len(fft_mag))}
    freq = np.fft.rfftfreq(len(y), d = spacing(x))
    maxima = argrelextrema(fft_mag, np.greater)[0]
    dif = (2.0 * fft_mag[maxima] - (fft_mag[maxima-1] + fft_mag[maxima+1])) / (2.0 * fft_mag[maxima])
    pts = [maxima[i] for i in range(len(maxima)) if dif[i] >= minDif ]
    return freq[pts], fft_phi[pts]

def evenSpacing(lst):
    test = len(set(np.round([lst[i] - lst[i-1] for i in range(1, len(lst))], decimals = 8)))
    if test > 1:
        return False
    return True

def spaceEvenly(x, y):
    meanSpacing = sum([x[i] - x[i-1] for i in range(1, len(x))]) / (len(x) - 1.0)
    spline = interp1d(x, y, kind='cubic')
    newX = np.arange(x[0], x[-1], meanSpacing)
    newY = spline(newX)

    return newX, newY
        
def spacing(x):
    return x[1] - x[0]

def fitSine(x, y, x_err, y_err):
    freq, angle = fftMaxima(x, y)
    print(freq)
    print(angle)
    init = [max(y) - min(y), freq[0], angle[0], 0]
    odr, out, fit = fit_data_to_model(x, y, x_err, y_err, sineModel, init)
    return out

def getPeriod(x, y):
    x_val = [ i.nominal_value for i in x ]
    x_err = [ i.std_dev for i in x ]
    y_val = [ i.nominal_value for i in y ]
    y_err = [ i.std_dev for i in y ]

    out = fitSine(x_val, y_val, x_err, y_err)

    return out

sine = lambda x, A, f, p: A * np.sin(x * 2.0 * np.pi * f + p)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    print("Testing evenSpacing")
    import random
    x = [i + random.gauss(0, 0.1) for i in range(10)]
    if evenSpacing(x) == False:
        print("Success - Correctly identified uneven list")
    else:
        print("Failure - Uneven list incorrectly identified")

    if evenSpacing(range(10)) == True:
        print("Success - Correctly identified even list")
    else:
        print("Failure - Even list incorrectly identified")

    print("Testing spaceEvenly")

    y = range(10)

    newX, newY = spaceEvenly(x, y)

    if evenSpacing(newX):
        print("Success - spaceEvenly appears to work")

    else:
        print("Failure - spaceEvenly failed")

    print("Testing getPeriod")

    F1 = 0.832
    F2 = 2.34
    A1 = 4.65
    A2 = 2.76
    P1 = 0.12
    P2 = -0.78

    x = np.linspace(0, 2, 100)
    y1 = sine(x, A1, F1, P1)
    y2 = sine(x, A2, F2, P2)
    y = y1 + y2
##    plt.plot(x, y1)
##    plt.plot(x, y2)
##    plt.plot(x, y)
##    plt.show()
    
  

    V = [ 10.32577373,  10.37217067,  10.38563709,  10.40279055,  10.41469995,
  10.40038288,  10.34226395,  10.23093946,  10.05760734,   9.89024585,
   9.84437577,   9.89251891,  10.12501703,  10.1971078,   10.28285867,
  10.3169841,   10.35169678,  10.37890353,  10.39951768,  10.4154415,
  10.39797001,  10.35791945,  10.23844714,  10.06747097,   9.93123516,
   9.86515774]

    V_err = [0.17353996] * len(V)

    x = [ 0.0 ,0.00368056,  0.00762731,  0.01152778,  0.01873843,  0.02655093,
  0.03042824,  0.03434028,  0.0380787,   0.04188657,  0.04587963,  0.04978009,
  0.06071759,  0.06445602,  0.07237269,  0.07619213,  0.08002315,  0.0840162,
  0.08780093,  0.09184028,  0.09570602,  0.1030787,   0.10702546,  0.11078704,
  0.1146412,   0.1183912 ]
    x_err = [0.001] * len(V)

    Vu = [ufloat(V[i], V_err[i]) for i in range(len(V))]
    Xu = [ufloat(x[i], x_err[i]) for i in range(len(x))]

    print(Vu)
    print(Xu)

##    freq, angle = fftMaxima(x, V)
##    init = [max(V) - min(V), freq[1], angle[1], 0]
##    odr, out, fit = fit_data_to_model(x, V, x_err, V_err, sineModel, init)
##
##    plt.plot(x, V)
##    plt.show()

    fit = getPeriod(Xu, Vu)
    print(fit)

    if abs(fit.nominal_value + 0.07317248932232885)< 0.001:
        print("Success - getPeriod appears to work")
    else:
        print("Failed - getPeriod does not appears to work")
