import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy

#===================================================================== Data file

df = pd.read_csv('ngc6886_2400s.0051.txt', sep='  ',names = ['wavelength','flux']) #alpha

ww = df['wavelength']
ff = df['flux']*10**12
nn = len(df)


#==================================================================== Data fittering

wl_he27,fl_he27=[],[]
wl_he60,fl_he60=[],[]
wl_ha,fl_ha=[],[]
wl_d,fl_d=[],[]
wl,fl=[],[]
for i in range(nn):
    if 6526 < ww[i] and ww[i] < 6527.5:
        wl_he27.append(ww[i])
        fl_he27.append(ff[i])    
    elif 6558 < ww[i] and ww[i] <= 6560.5:
        wl_he60.append(ww[i])
        fl_he60.append(ff[i])        
    elif 6562.3 <= ww[i] and ww[i]< 6564.7:
        wl_ha.append(ww[i])
        fl_ha.append(ff[i])
    elif 6546.5 < ww[i] and ww[i] <6550:
        wl_d.append(ww[i])
        fl_d.append(ff[i])
    else:
        wl.append(ww[i])
        fl.append(ff[i])



wl_w,fl_w=[],[]
wl_r,fl_r=[],[]
for j in range(len(wl)):
    if 6542<wl[j] and wl[j]<6551:
        wl_r.append(wl[j])
        fl_r.append(fl[j])
    elif 6558.<wl[j] and wl[j]<6559.:
        wl_w.append(wl[j])
        fl_w.append(fl[j])
    elif 6564.8<wl[j] and wl[j]<6567.6:
        wl_w.append(wl[j])
        fl_w.append(fl[j])

wl_w = np.array(wl_w)
fl_w = np.array(fl_w)

wing = np.stack((wl_w,fl_w),axis=1)

wl_r = np.array(wl_r)
fl_r = np.array(fl_r)

raman = np.stack((wl_r,fl_r),axis=1)



#===================================================================== Gaussian Fitting

def gauss(x,A, mu, sigma):
    f1 = A*np.exp(-(x-mu)**2/(2*sigma**2))
    return f1


# He II 6527
p_he27 = np.max(fl_he27)
cwl_he27 = np.mean(wl_he27)
std_he27 = np.std(wl_he27)


x1 = np.linspace(6500,6570,100)

popt1,pcov1 = curve_fit(gauss,wl_he27,fl_he27,p0=[p_he27,cwl_he27,std_he27],maxfev = 5000)

A1 = popt1[0] ; mu1 = popt1[1] ; sigma1 = popt1[2]

# He II 6560
p_he60 = np.max(fl_he60)
cwl_he60 = np.mean(wl_he60)
std_he60 = np.std(wl_he60)

x2 = np.linspace(6558.8,np.max(wl_he60),100)

popt2,pcov2 = curve_fit(gauss,wl_he60,fl_he60,p0=[p_he60,cwl_he60,std_he60],maxfev = 5000)
A2 = popt2[0] ; mu2 = popt2[1] ; sigma2 = popt2[2] 

# H alpha line
p_ha = np.max(fl_ha)
cwl_ha = np.mean(wl_ha)
std_ha = np.std(wl_ha)

x3 = np.linspace(np.min(wl_ha),6565,100)

popt3,pcov3 = curve_fit(gauss,wl_ha,fl_ha,p0=[p_ha,cwl_ha,std_ha],maxfev = 5000)

A3 = popt3[0] ; mu3 = popt3[1] ; sigma3 = popt3[2]


# H alpha wing
p_w = np.max(fl_w)
cwl_w = np.mean(wl_w)
std_w = np.std(wl_w)

x4 = np.linspace(np.min(wl_w),np.max(wl_w),100)

popt4,pcov4 = curve_fit(gauss,wl_w,fl_w,p0=[p_w,cwl_w,std_w],maxfev = 5000)

A4 = popt4[0] ; mu4 = popt4[1] ; sigma4 = popt4[2]

# Raman He II 6545
p_r = np.max(fl_r)
cwl_r = np.mean(wl_r)
std_r = np.std(wl_r)

x5 = np.linspace(np.min(wl_r),np.max(wl_r),100)

popt5,pcov5 = curve_fit(gauss,wl_r,fl_r,p0=[p_r,cwl_r,std_r],maxfev = 5000)

A5 = popt5[0] ; mu5 = popt5[1] ; sigma5 = popt5[2] 


#=======================================================================Fitting

# Method1

def Alpha1(x, A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3, A4, mu4, sigma4, A5, mu5, sigma5):
    f1 = A2*np.exp(-(x-mu2)**2/(2*sigma2**2))+A3*np.exp(-(x-mu3)**2/(2*sigma3**2))+\
        A4*np.exp(-(x-mu4)**2/(2*sigma4**2))+A1*np.exp(-(x-mu1)**2/(2*sigma1**2))+\
        A5*np.exp(-(x-mu5)**2/(2*sigma5**2))
    return f1


xx = np.linspace(6520,6600,10000)
yy = Alpha1(xx,A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3, A4, mu4, sigma4, A5, mu5, sigma5)

yy = np.array(yy)
offset = 0.0055

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


#========================================================================= Plotting

# All spectrum
plt.plot(ww,ff, label = 'Original Data',color ='black',alpha = 0.6,lw =1)
#plt.plot(xx_f,yy_f,color='red', linestyle='--', label = 'Gaussian fitting')
#plt.xlim(6520,6570)
#plt.ylim(0,0.06) 
plt.xlabel(r'Wavelength[$\AA$]')
plt.ylabel(r'Flux[$10^{-12}$erg$s^{-1}$$cm^{2}$$\AA^{-1}$]')
# plt.text(6521,0.028,'HeII 6527',fontdict={'fontsize':'10','fontweight':'bold'})
# plt.annotate('',xy=(6525.7,0.02),xytext=(6524.8,0.027),arrowprops=dict(arrowstyle='->',color='black',lw=2))
# plt.text(6532,0.018,'Raman HeII 6545',fontdict={'fontsize':'10','fontweight':'bold'})
# plt.annotate('',xy=(6543,0.01),xytext=(6539.5,.017),arrowprops=dict(arrowstyle='->',color='black',lw=2))
# plt.text(6550,0.031,'HeII 6560',fontdict={'fontsize':'10','fontweight':'bold'})
# plt.annotate('',xy=(6557.3,0.025),xytext=(6553.3,0.03),arrowprops=dict(arrowstyle='->',color='black',lw=2))
# plt.text(6566,0.04, r'H$\alpha$',fontdict={'fontsize':'13','fontweight':'bold'})
plt.title('NGC 6886')
plt.grid(linestyle=':')
plt.legend(loc = 'upper right')
plt.show()



#Gaussian fitting(all)

plt.plot(ww,ff, label = 'Original Data',color ='black',alpha = 0.6,lw =1)
plt.plot(xx_f,yy_f,color='red', linestyle='--', label = 'Gaussian fitting')
# plt.xlim(6520,6570)
# plt.ylim(0,0.06) 
plt.xlabel(r'Wavelength[$\AA$]')
plt.ylabel(r'Flux[$10^{-12}$erg$s^{-1}$$cm^{2}$$\AA^{-1}$]')
plt.text(6521,0.028,'HeII 6527',fontdict={'fontsize':'10','fontweight':'bold'})
plt.annotate('',xy=(6525.7,0.02),xytext=(6524.8,0.027),arrowprops=dict(arrowstyle='->',color='black',lw=2))
plt.text(6532,0.018,'Raman HeII 6545',fontdict={'fontsize':'10','fontweight':'bold'})
plt.annotate('',xy=(6543,0.01),xytext=(6539.5,.017),arrowprops=dict(arrowstyle='->',color='black',lw=2))
plt.text(6550,0.031,'HeII 6560',fontdict={'fontsize':'10','fontweight':'bold'})
plt.annotate('',xy=(6557.3,0.025),xytext=(6553.3,0.03),arrowprops=dict(arrowstyle='->',color='black',lw=2))
plt.text(6566,0.04, r'H$\alpha$',fontdict={'fontsize':'13','fontweight':'bold'})
plt.title('NGC 6886')
plt.grid(linestyle=':')
plt.legend(loc = 'upper right')
plt.show()


#Gaussian fitting(xlim)

plt.plot(ww,ff, label = 'Original Data',color ='black',alpha = 0.6,lw =1)
plt.plot(xx_f,yy_f,color='red', linestyle='--', label = 'Gaussian fitting')
plt.xlim(6520,6570)
#plt.ylim(0,0.06) 
plt.xlabel(r'Wavelength[$\AA$]')
plt.ylabel(r'Flux[$10^{-12}$erg$s^{-1}$$cm^{2}$$\AA^{-1}$]')
plt.text(6521,0.028,'HeII 6527',fontdict={'fontsize':'10','fontweight':'bold'})
plt.annotate('',xy=(6525.7,0.02),xytext=(6524.8,0.027),arrowprops=dict(arrowstyle='->',color='black',lw=2))
plt.text(6532,0.018,'Raman HeII 6545',fontdict={'fontsize':'10','fontweight':'bold'})
plt.annotate('',xy=(6543,0.01),xytext=(6539.5,.017),arrowprops=dict(arrowstyle='->',color='black',lw=2))
plt.text(6550,0.031,'HeII 6560',fontdict={'fontsize':'10','fontweight':'bold'})
plt.annotate('',xy=(6557.3,0.025),xytext=(6553.3,0.03),arrowprops=dict(arrowstyle='->',color='black',lw=2))
plt.text(6566,0.04, r'H$\alpha$',fontdict={'fontsize':'13','fontweight':'bold'})
plt.title('Raman HeII 6545 in NGC 6886')
plt.grid(linestyle=':')
plt.legend(loc = 'upper right')
plt.show()


#Gaussian fitting(xlim, ylim)

plt.plot(ww,ff, label = 'Original Data',color ='black',alpha = 0.6,lw =1)
plt.plot(xx_f,yy_f,color='red', linestyle='--', label = 'Gaussian fitting')
plt.xlim(6520,6570)
plt.ylim(0,0.06) 
plt.xlabel(r'Wavelength[$\AA$]')
plt.ylabel(r'Flux[$10^{-12}$erg$s^{-1}$$cm^{2}$$\AA^{-1}$]')
plt.text(6521,0.028,'HeII 6527',fontdict={'fontsize':'10','fontweight':'bold'})
plt.annotate('',xy=(6525.7,0.02),xytext=(6524.8,0.027),arrowprops=dict(arrowstyle='->',color='black',lw=2))
plt.text(6532,0.018,'Raman HeII 6545',fontdict={'fontsize':'10','fontweight':'bold'})
plt.annotate('',xy=(6543,0.01),xytext=(6539.5,.017),arrowprops=dict(arrowstyle='->',color='black',lw=2))
plt.text(6550,0.031,'HeII 6560',fontdict={'fontsize':'10','fontweight':'bold'})
plt.annotate('',xy=(6557.3,0.025),xytext=(6553.3,0.03),arrowprops=dict(arrowstyle='->',color='black',lw=2))
plt.text(6566,0.04, r'H$\alpha$',fontdict={'fontsize':'13','fontweight':'bold'})
plt.title('Raman HeII 6545 in NGC 6886')
plt.grid(linestyle=':')
plt.legend(loc = 'upper right')
plt.show()




# Method2

# xx = np.linspace(6500,6600,10000)
# xlist,ylist = [],[]
# for k in range(len(xx)):
#     xnew = xx[k]
#     if xnew >= np.min(wl_he27) and xnew <= np.max(wl_he27):
#         ynew = gauss(xnew,*popt1)
#         xlist.append(xnew)
#         ylist.append(ynew)
#     elif xnew >= np.min(wl_he60) and xnew <= np.max(wl_he60):
#         ynew = gauss(xnew,*popt2)
#         xlist.append(xnew)
#         ylist.append(ynew)
#     elif xnew >= np.min(wl_ha) and xnew <= 6565:
#         ynew = gauss(xnew,*popt3)
#         xlist.append(xnew)
#         ylist.append(ynew)
#     elif xnew >= np.min(wl_r) and xnew <= np.max(wl_r):
#         ynew = gauss(xnew,*popt5)
#         xlist.append(xnew)
#         ylist.append(ynew)
#     elif xnew >= np.min(wl_r) and xnew <= np.min(wl_he60):
#         ynew = gauss(xnew,*popt4)
#         xlist.append(xnew)
#         ylist.append(ynew)
#     elif xnew>= np.max(wl_ha) and xnew <= np.max(wl_w):
#         ynew = gauss(xnew,*popt4)
#         xlist.append(xnew)
#         ylist.append(ynew)
#     # else:
#     #     ynew = 0.005
#     #     xlist.append(xnew)
#     #     ylist.append(ynew)


# xx = np.array(xlist)
# yy = np.array(ylist)


# plt.plot(ww,ff)
# plt.plot(xx,yy)
# plt.xlim(6520,6570)
# plt.ylim(0,0.06)
# plt.show()