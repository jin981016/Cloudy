import numpy as np 
import matplotlib.pyplot as plt
import matplotlib 
from scipy.optimize import curve_fit
import pandas as pd 

df1 = np.genfromtxt('cngc6741_3600s.0030.txt', names=["wavelength","flux"]) 
df2 = np.genfromtxt('cngc6741_3600s.0055.txt', names=["wavelength","flux"])

ww1 = df1['wavelength']
ff1 = df1['flux']*10**(11)

ww2 = df2['wavelength']
ff2 = df2['flux']*10**(11)

print(ww1)

# Gaussian fitting
def gauss(A,x,mu,sigma):
    f = A*np.exp(-(x-mu)**2/(2.*sigma**2))
    return f

    

#popt, pcov = curve_fit(gauss, )






# # Plotting data and fitting

# plt.subplot(2,2,1)
# plt.plot(ww1, ff1, label = 'Data', color = 'black')
# # plt.plot(ww, gg1, label = 'Gaussian Fitting', 'r-')
# plt.xlabel(r'Wavelength[$\AA$]')
# plt.ylabel(r'Flux($10^{-11}$erg$s$$^{-1}$$cm^{-2}$$\AA^{-1}$)')
# plt.xlim(4840,4865)
# plt.grid()
# plt.legend()

# plt.subplot(2,2,3)
# plt.plot(ww1,ff1, label = 'Data', color = 'black')
# #plt.plot(ww, gg2, label = 'Gaussian Fitting','r-')
# plt.xlabel(r'Wavelength[$\AA$]')
# plt.ylabel(r'Flux($10^{-11}$erg$s$$^{-1}$$cm^{-2}$$\AA^{-1}$)')
# plt.xlim(4845,4865)
# plt.ylim(0,0.1)
# plt.grid()
# plt.legend()

# plt.subplot(2,2,2)
# plt.plot(ww2, ff2, label = 'Data', color = 'black')
# # plt.plot(ww, gg1, label = 'Gaussian Fitting','r-')
# plt.xlabel(r'Wavelength[$\AA$]')
# plt.ylabel(r'Flux($10^{-11}$erg$s$$^{-1}$$cm^{-2}$$\AA^{-1}$)')
# plt.xlim(6520,6600)
# plt.grid()
# plt.legend()

# plt.subplot(2,2,4)
# plt.plot(ww2,ff2, label = 'Data', color = 'black')
# # plt.plot(ww, gg2, label = 'Gaussian Fitting','r-')
# plt.xlabel(r'Wavelength[$\AA$]')
# plt.ylabel(r'Flux($10^{-11}$erg$s$$^{-1}$$cm^{-2}$$\AA^{-1}$)')
# plt.xlim(6520,6570)
# plt.ylim(0,0.01)
# plt.grid()
# plt.legend()


# plt.show()
