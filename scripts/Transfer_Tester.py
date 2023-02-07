import matplotlib.font_manager
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
import seaborn as sns
import subprocess
from matplotlib.lines import Line2D
from matplotlib import rc


transferprepath = "/n/holyscratch01/dvorkin_lab/ndeporzio/RelAxiFast_"
transferpostpath = "/Boltzmann_2/transfer_files_0/"

omega_cdm_LCDM = 0.1127
omega_ax = np.linspace(0.06, 0.10, 17)*omega_cdm_LCDM
m_ax = np.logspace(-32., -22., 51)

Nruns = int(len(omega_ax)*len(m_ax))

Tcdm_zk__m = len(m_ax)*[0]
Tb_zk__m = len(m_ax)*[0]
Tphi_zk__m = len(m_ax)*[0]
klists = len(m_ax)*[0]
zlists = len(m_ax)*[0]

Tcdminterp_zk__m = len(m_ax)*[0]
Tbinterp_zk__m = len(m_ax)*[0]
Tphiinterp_zk__m = len(m_ax)*[0]

for idx in range(Nruns):
    m_idx = np.mod(idx, len(m_ax))
    o_idx = np.floor_divide(idx, len(m_ax))

    zlist = []
    Tcdm_zk = []
    Tb_zk = []
    Tphi_zk = []

    if (o_idx==12): 
        #The 9% CDM case... 

        #Collect all the redshift values 
        filenames = os.listdir(transferprepath+str(idx+1)+transferpostpath)
        for idx, val in enumerate(filenames): 
            if (val[0:15]=='_transfer_out_z'): 
                zlist.append(float(val[15:]))
        zlist = np.sort(zlist)     
        klist = np.loadtxt(transferprepath+str(idx+1)+transferpostpath+'_transfer_out_z'+f"{zlist[0]:.3f}")[:,0]
   
        for zidx, zval in enumerate(zlist):
            transfers=np.loadtxt(transferprepath+str(idx+1)+transferpostpath+'_transfer_out_z'+f"{zlist[0]:.3f}")
            Tcdm_zk.append(transfers[:,1]) 
            Tb_zk.append(transfers[:,2]) 
            Tphi_zk.append(transfers[:,6]) 

        Tcdm_zk__m.append(Tcdm_zk)        
        Tb_zk__m.append(Tb_zk)        
        Tphi_zk__m.append(Tphi_zk)        
        klists.append(klist)
        zlists.append(zlist)

for idx in range(len(m_ax)): 
    Tcdminterp_zk__m[idx] = scipy.interpolate.interp2d(klists[idx], zlists[idx], Tcdm_zk__m[idx], bounds_error=True)    
    Tbinterp_zk__m[idx] = scipy.interpolate.interp2d(klists[idx], zlists[idx], Tb_zk__m[idx], bounds_error=True)    
    Tphiinterp_zk__m[idx] = scipy.interpolate.interp2d(klists[idx], zlists[idx], Tphi_zk__m[idx], bounds_error=True)    


zplot = zlists[0]

np.savetxt("TEST_zlists.txt", zlists)
np.savetxt("TEST_zplot.txt", zplot)
np.savetxt("TEST_klists.txt", klists)
np.savetxt("TEST_mplot.txt", m_ax) 

Tcdm_mz__k = np.zeros((len(zinterp), len(m_ax)))
Tb_mz__k = np.zeros((len(zinterp), len(m_ax)))
Tphi_mz__k = np.zeros((len(zinterp), len(m_ax)))

# First for 10^{-4} [Mpc]^{-1}
for m_idx, m_val in enumerate(m_ax): 
    for z_idx, z_val in enumerate(z_plot): 
        Tcdm_mz__k[z_idx, m_idx] = Tcdminterp_zk__m[m_idx](np.power(10., -4.)/0.70148, z_val)
        Tb_mz__k[z_idx, m_idx] = Tbinterp_zk__m[m_idx](np.power(10., -4.)/0.70148, z_val)
        Tphi_mz__k[z_idx, m_idx] = Tphiinterp_zk__m[m_idx](np.power(10., -4.)/0.70148, z_val)
np.savetxt("TEST_Tcdm_k1.0e-4.txt", Tcdm_mz__k) 
np.savetxt("TEST_Tb_k1.0e-4.txt", Tb_mz__k) 
np.savetxt("TEST_Tphi_k1.0e-4.txt", Tphi_mz__k) 

# Then for 10^{-0.5} [Mpc]^{-1}
for m_idx, m_val in enumerate(m_ax): 
    for z_idx, z_val in enumerate(z_plot): 
        Tcdm_mz__k[z_idx, m_idx] = Tcdminterp_zk__m[m_idx](np.power(10., -0.5)/0.70148, z_val)
        Tb_mz__k[z_idx, m_idx] = Tbinterp_zk__m[m_idx](np.power(10., -0.5)/0.70148, z_val)
        Tphi_mz__k[z_idx, m_idx] = Tphiinterp_zk__m[m_idx](np.power(10., -0.5)/0.70148, z_val)
np.savetxt("TEST_Tcdm_k1.0e-0.5.txt", Tcdm_mz__k) 
np.savetxt("TEST_Tb_k1.0e-0.5.txt", Tb_mz__k) 
np.savetxt("TEST_Tphi_k1.0e-0.5.txt", Tphi_mz__k) 

