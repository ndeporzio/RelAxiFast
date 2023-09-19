import matplotlib.font_manager
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
import seaborn as sns
import subprocess
from matplotlib.lines import Line2D
from matplotlib import rc
import pylab as pl 

omega_cdm_LCDM = 0.1127
meq = -27.44

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelAxiFast.nosync"

c = 3.0e8

Hdimensionless = np.loadtxt(rfpath+"/axion_cad2.dat", skiprows=0)[:, [0, 2]]
H = np.matmul(Hdimensionless, [[1, 0],[0, 100.*1000./c]]) #rescale second column of Hdimensionless 
Xcomoving = test=np.loadtxt("/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelAxiFast.nosync/Comoving_Horizon.txt", skiprows=0)
kcrossing = np.array(Xcomoving)
kcrossing[:,1] = 1./kcrossing[:,1]

m_ax = np.logspace(-32, -22, 6)

m_ax_osc = np.loadtxt(rfpath+"/plots/Figure_1_m_ax.txt")*(1.5e29)
a_osc = np.loadtxt(rfpath+"/plots/Figure_1_axion_aosc.txt")[:,2]


H_interp = scipy.interpolate.interp1d(H[:,0], H[:,1])
kcrossing_interp = scipy.interpolate.interp1d(kcrossing[:,0], kcrossing[:,1], bounds_error=False, fill_value=((kcrossing[0,1],kcrossing[-1,1])))
mosc_interp = scipy.interpolate.interp1d(a_osc, m_ax_osc, bounds_error=False, fill_value=((np.max(m_ax_osc), np.min(m_ax_osc))))
aosc_interp = scipy.interpolate.interp1d(m_ax_osc, a_osc)



plt.figure()
plt.plot(H[:,0], H[:,1], label=r"$H(a)$")
plt.plot(a_osc, m_ax_osc, label=r"$m_{\rm osc}$", color='grey')
plt.scatter(aosc_interp(m_ax*(1.5e29)), m_ax*(1.5e29), color="black")
for idx in range(len(m_ax)): 
    plt.text(0.4*aosc_interp(m_ax[idx]*(1.5e29)), 0.4*m_ax[idx]*(1.5e29), f"{m_ax[idx]:.0e}", fontsize=30)

plt.plot(kcrossing[:,0], kcrossing[:,1], color='green', label=r"$k_{\rm crossing}$ (comoving)")
kJ=[]
for idx in range(len(m_ax)):
    kJ.append(np.pi * np.sqrt(m_ax[idx] * H[:, 1] * 1.5e29))
    line = plt.plot(H[:,0], kJ[idx], color='red', alpha=(idx+1)/(len(m_ax)+1), label=r"$k_{\rm Jeans}(m_\phi=$"+f"{m_ax[idx]:.0e}"+r" eV)")
plt.fill_between(H[:,0], kJ[0], 1.0e8, color='red', alpha=0.05)
#plt.fill_between(H[:,0], kJ, kcrossing_interp(H[:,0]), color='red', alpha=0.05)
plt.fill_between(kcrossing[:,0], kcrossing[:,1], 0, color='green', alpha=.1)
plt.fill_between(kcrossing[:,0], mosc_interp(kcrossing[:,0]), 0., color='grey', alpha=.1)
#plt.fill_between(kcrossing[:,0], mosc_interp(kcrossing[:,0]), np.min(m_ax_osc), color='grey', alpha=.1)
#plt.fill_between(kcrossing[:,0], mosc_interp(kcrossing[:,0]), 1.0e8, color='grey', alpha=.1)

plt.plot([np.min(kcrossing[:,0]), np.max(a_osc)], [np.power(10.,-0.5), np.power(10.,-0.5)], linestyle='dashdot', color='black', linewidth=5)
plt.text(1.0e-7, 0.5e0, r"$k_{\rm *}$", fontsize=40)

plt.plot([np.min(kcrossing[:,0]), np.max(a_osc)], [1.0e-4, 1.0e-4], linestyle='dashed', color='black', linewidth=5)
plt.text(1.0e-7, 2.0e-4, r"$k_{\rm ref}$", fontsize=40)

plt.plot([1./201, 1./201], [0., mosc_interp(1./201)], linestyle='dotted', color='black', linewidth=5)
plt.text(2.8e-3, 3.0e-7, r"$z_{\rm ini}=200$", rotation='vertical', fontsize=40)

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$a$')
plt.ylabel(r'${\rm [Mpc}^{-1}{\rm ]}$')
plt.legend(fontsize=30, loc="upper right")
plt.text(1.0e-5, 1.0e6, "Jeans Suppressed", color='red', alpha=0.5)
plt.text(3.0e-10, 0.3e-7, "Horizon Suppressed", color='green', alpha=0.5)
plt.text(3.0e-10, 1.0e-2, "Oscillation Suppressed", color='grey', alpha=0.5)
plt.ylim((1.0e-8, np.max(m_ax_osc)))
plt.xlim((np.min(kcrossing[:,0]), np.max(a_osc)))
plt.show()
