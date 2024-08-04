import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy



df1 = pd.read_csv('cngc6884_3600s.0030.txt', sep='  ', names = ['wavelength','flux'])


ww1 = np.array(df1['wavelength'])
ff1 = np.array(df1['flux']*10**12)


# Data

def gauss(x, A, mu, sigma):
    f = A*np.exp(-(x-mu)**2/(2*sigma**2))
    return f



wl_r,wl_he,wl_hb,flux_r,flux_he,flux_hb = [],[],[],[],[],[]
for i in range(len(df1)):
    if 4848. < ww1[i] and ww1[i] < 4851.5:
        wl_r.append(ww1[i])
        flux_r.append(ff1[i])
    elif 4858.< ww1[i] and ww1[i] <= 4859.3:
        wl_he.append(ww1[i])
        flux_he.append(ff1[i])    
    elif 4860. <= ww1[i] and ww1[i] <4861.1:
        wl_hb.append(ww1[i])
        flux_hb.append(ff1[i])  



# H beta line 4861
p_hb = np.max(flux_hb)
cwl_hb = np.mean(wl_hb)
sigma_hb = np.std(wl_hb)

x1 = np.linspace(np.min(wl_hb),np.max(wl_hb),1000)
y1 = gauss(x1,p_hb,cwl_hb,sigma_hb)

popt1,pcov1 = curve_fit(gauss,x1,y1, p0=[p_hb,cwl_hb,sigma_hb])

A1 = popt1[0]
mu1 = popt1[1]
sigma1 = popt1[2]

yy1 = gauss(x1,*popt1)

# plt.plot(ww1,ff1,label = 'Original Data', color ='black',alpha = 0.6,lw =1)
# plt.plot(x1, yy1,color = 'red', linestyle='--', label = 'Gaussian Fitting')
# plt.xlabel(r'Wavelength[$\AA$]')
# plt.ylabel(r'Flux[$10^{-12}$erg$s^{-1}$$cm^{2}$$\AA^{-1}$]')
# plt.text(4843,0.13,'Raman HeII 4851',fontdict={'fontsize':'10','fontweight':'bold'})
# plt.annotate('',xy=(4848.5,0.09),xytext=(4847.3,0.12),arrowprops=dict(arrowstyle='->',color='black',lw=2))
# plt.text(4853,0.15,'HeII 4859',fontdict={'fontsize':'10','fontweight':'bold'})
# plt.annotate('',xy=(4857.3,0.12),xytext=(4856,0.14),arrowprops=dict(arrowstyle='->',color='black',lw=2))
# plt.text(4862.5,0.2, r'H$\beta$',fontdict={'fontsize':'13','fontweight':'bold'})
# plt.grid(linestyle=':')
# plt.xlim(4840,4870)
# #plt.ylim(0,0.6)
# plt.legend() 
# plt.show()

# He II 4859
p_he = np.max(flux_he)
cwl_he = np.mean(wl_he)
sigma_he = np.std(wl_he)

x2 = np.linspace(np.min(wl_he),np.max(wl_he),1000)
y2 = gauss(x2,p_he,cwl_he,sigma_he)

popt2,pcov2 = curve_fit(gauss,x2,y2,p0=[p_he,cwl_he,sigma_he])

A2 = popt2[0]
mu2 = popt2[1]
sigma2 = popt2[2]



# plt.plot(ww1,ff1,label = 'Original Data', color ='black',alpha = 0.6,lw =1)
# plt.plot(x2, y2,color = 'red', linestyle='--', label = 'Gaussian Fitting')
# plt.xlabel(r'Wavelength[$\AA$]')
# plt.ylabel(r'Flux[$10^{-12}$erg$s^{-1}$$cm^{2}$$\AA^{-1}$]')
# plt.text(4843,0.13,'Raman HeII 4851',fontdict={'fontsize':'10','fontweight':'bold'})
# plt.annotate('',xy=(4848.5,0.09),xytext=(4847.3,0.12),arrowprops=dict(arrowstyle='->',color='black',lw=2))
# plt.text(4853,0.15,'HeII 4859',fontdict={'fontsize':'10','fontweight':'bold'})
# plt.annotate('',xy=(4857.3,0.12),xytext=(4856,0.14),arrowprops=dict(arrowstyle='->',color='black',lw=2))
# plt.text(4862.5,0.2, r'H$\beta$',fontdict={'fontsize':'13','fontweight':'bold'})
# plt.grid(linestyle=':')
# plt.xlim(4840,4870)
# plt.ylim(0,0.6)
# plt.legend() 
# plt.show()

# Raman He II 4851
p_r = np.max(flux_r)
cwl_r = np.mean(wl_r)
sigma_r = np.std(wl_r)

x3 = np.linspace(np.min(wl_r),np.max(wl_r),1000)
y3 = gauss(x3,p_r,cwl_r,sigma_r)

popt3,pcov3 = curve_fit(gauss,x3,y3,p0=[p_r,cwl_r,sigma_r])

A3 = popt3[0]
mu3 = popt3[1]
sigma3 = popt3[2]



#Guassian Fitthing

def gauss_fit(x,A1,mu1,sigma1,A2,mu2,sigma2,A3,mu3,sigma3):
    ff = gauss(x,A1,mu1,sigma1)+gauss(x,A2,mu2,sigma2)+gauss(x,A3,mu3,sigma3)
    return ff

xx = np.linspace(4800,4900,10000)
yy = gauss_fit(xx,A1,mu1,sigma1,A2,mu2,sigma2,A3,mu3,sigma3)

offset = 0.03

xlist,ylist =[],[]
for ii in range(len(yy)):
    if yy[ii] <= offset:
        yy[ii] = offset
        xlist.append(xx[ii])
        ylist.append(yy[ii])
    elif yy[ii] >= offset:
        yy[ii] = yy[ii]
        xlist.append(xx[ii])
        ylist.append(yy[ii])            

xx_f = np.array(xlist)
yy_f = np.array(ylist)


plt.plot(ww1,ff1,label = 'Original Data', color ='black',alpha = 0.6,lw =1)
plt.plot(xx_f, yy_f,color = 'red', linestyle='--', label = 'Gaussian Fitting')
plt.xlabel(r'Wavelength[$\AA$]')
plt.ylabel(r'Flux[$10^{-12}$erg$s^{-1}$$cm^{2}$$\AA^{-1}$]')
plt.text(4843,0.13,'Raman HeII 4851',fontdict={'fontsize':'10','fontweight':'bold'})
plt.annotate('',xy=(4848.5,0.09),xytext=(4847.3,0.12),arrowprops=dict(arrowstyle='->',color='black',lw=2))
plt.text(4853,0.15,'HeII 4859',fontdict={'fontsize':'10','fontweight':'bold'})
plt.annotate('',xy=(4857.3,0.12),xytext=(4856,0.14),arrowprops=dict(arrowstyle='->',color='black',lw=2))
plt.text(4862.5,0.2, r'H$\beta$',fontdict={'fontsize':'13','fontweight':'bold'})
plt.grid(linestyle=':')
plt.xlim(4840,4870)
plt.ylim(0,0.4)
plt.title('Raman HeII 4851 of NGC 6884')
plt.legend(loc='upper right')
plt.savefig('4851_NGC6884.png') 
#plt.show()

