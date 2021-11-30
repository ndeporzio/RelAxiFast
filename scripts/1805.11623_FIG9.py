import matplotlib.pyplot as plt
import numpy as np
import os
import scipy 
import seaborn as sns
import subprocess

from matplotlib.lines import Line2D

sns.set()

axioncamb_data_eulbias = []
axioncamb_data_lagbias = []
axioncamb_data_tf = []

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast/"
rfpath_outputsuffix = "output/result-0/"
outpath = "/Users/nicholasdeporzio/Desktop/"

Mnu = np.array([0.0, 90.0e-3]) # Units: eV
redshift = 0.7
kmin = 5.0e-5
kmax = 1.5
Nk = 50

h_lcdm = 0.67
omega_nu = Mnu/93.2
omega_cdm_LCDM = 0.12
omega_b_LCDM = 0.022
#Omega_M_LCDM = 0.27464
Omega_M_LCDM = (omega_b_LCDM + omega_cdm_LCDM)/np.power(h_lcdm, 2.)

f_cdm_LCDM = omega_cdm_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
f_b_LCDM = omega_b_LCDM/(omega_cdm_LCDM + omega_b_LCDM)

#h = np.sqrt((omega_b_LCDM+omega_cdm_LCDM+omega_nu)/Omega_M_LCDM)
h = np.array(len(Mnu)*[h_lcdm])

omega_cdm = omega_cdm_LCDM - omega_nu
f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)

# Set solver to axionCAMB
os.chdir(rfpath+'/include/')
reading_file = open("common.h", "r")
new_file_content = ""
for line in reading_file:
    stripped_line = line.strip()
    new_line = stripped_line.replace(
        "define boltzmann_tag  _CLASS_", "define boltzmann_tag  _AXIONCAMB_"
    ).replace(
        "define boltzmann_tag  _CAMB_", "define boltzmann_tag  _AXIONCAMB_"
    )
    new_file_content += "    " + new_line +"\n"
reading_file.close()
writing_file = open("common.h", "w")
writing_file.write(new_file_content)
writing_file.close()
os.chdir(rfpath)
os.system('make')

# Calculate Bias using axionCAMB
for m_idx, m_val in enumerate(Mnu):  
    print("Running RelicFast + axionCAMB for \Sigma m_nu = ", Mnu[m_idx], ", z = ", redshift)
    
    os.system('rm -rf '+rfpath_outputsuffix)
    os.system('rm -rf '+'Boltzmann_2')
    os.system('rm ./run.ini')
    
    # RUN RelicFast + axionCAMB for each neutrino mass choice 
    reading_file = open("1805.11623.ini", "r")
    new_file_content = ""
    for line in reading_file:
        stripped_line = line.strip()
        new_line = stripped_line.replace(
            "mnu1 = 0.0", "mnu1 = "+f'{Mnu[m_idx]/3.:.2e}'
        ).replace(
            "mnu2 = 0.0", "mnu2 = "+f'{Mnu[m_idx]/3.:.2e}'
        ).replace(
            "m_SN = 0.02", "m_SN = "+f'{Mnu[m_idx]/3.:.2e}'
        ).replace(
            "hubble = 0.701", "hubble = "+f'{h[m_idx]:.4f}'
        ).replace(
            "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{redshift:.3f}'
        ).replace(
            "kbot = 1e-4", "kbot = "+f'{kmin:.3e}'
        ).replace(
            "ktop = 0.7", "ktop = "+f'{kmax:.3e}'
        ).replace(
            "N_klong = 1", "N_klong = "+f'{Nk:d}'
        ).replace(
            "omegab = 0.02226", "omegab = "+f'{omega_b_LCDM:.3f}'
        ).replace(
            "omegac = 0.11271", "omegac = "+f'{omega_cdm_LCDM:.3f}'
        ).replace(
            "Mhalo_min = 1e13", "Mhalo_min = "+f'{1e13/h_lcdm:.3e}'
        )
        if (m_val!=0.0): 
            new_line = new_line.replace(
                "tag_sterile_nu = 0", "tag_sterile_nu = 1"
            )
        new_file_content += "    " + new_line +"\n"
    reading_file.close()
    writing_file = open("run.ini", "w")
    writing_file.write(new_file_content)
    writing_file.close()
    
    os.system('./relicfast run.ini')
    
    # Collect bias for requested redshift 
#    axioncamb_data_eulbias.append(
#        np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{redshift:.2f}'+'_M13.48_Nk50.dat', skiprows=1)
#    )
#    axioncamb_data_lagbias.append(
#        np.loadtxt(rfpath_outputsuffix+'bias_Lagrangian_z'+f'{redshift:.2f}'+'_M13.48_Nk50.dat', skiprows=1)
#    )
    axioncamb_data_eulbias.append(
        np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{redshift:.2f}'+'_M'
                   +f'{13-np.log10(h_lcdm):.2f}'
                   +'_Nk50.dat', skiprows=1)
    )
    axioncamb_data_lagbias.append(
        np.loadtxt(rfpath_outputsuffix+'bias_Lagrangian_z'+f'{redshift:.2f}'+'_M'
                   +f'{13-np.log10(h_lcdm):.2f}'
                   +'_Nk50.dat', skiprows=1)
    )
    axioncamb_data_tf.append(
        np.loadtxt(rfpath_outputsuffix+'z'+f'{redshift:.2f}'+'TF_CAMB.dat', skiprows=1)
    )
    
    os.system('mv ./run.ini ./run_axioncamb_'+str(m_idx)+'.ini')

kplot = np.geomspace(10**-4, 10**0, 41) #Units: h Mpc^-1
kfs = np.array([0.08*(mval/3./0.1)*h[m_idx]/np.sqrt(1+redshift) 
    for midx, mval in enumerate(Mnu)])

colors = sns.color_palette("seismic", len(Mnu))

plt.figure(figsize=(15, 7.5))
plt.xscale('log')

for m_idx, m_val in enumerate(Mnu): 
    ref = np.loadtxt(rfpath+"/scripts/1805.11623_FIG9_REF_"+str(int(1000.*m_val))+"meV.csv", delimiter=',')
    refinterp = scipy.interpolate.interp1d(ref[:,0], ref[:,1])
    ref_refactor = refinterp(10**-4)
    
    kvals = axioncamb_data_lagbias[m_idx][:, 0]/h_lcdm # h/Mpc
    lagbiasvals = axioncamb_data_lagbias[m_idx][:, 1]

    lagbiasinterp = scipy.interpolate.interp1d(kvals, lagbiasvals)
    lagbiasplot = lagbiasinterp(kplot) # b(k[1/Mpc])
    
    plt.plot(
        kplot, 
        lagbiasplot/lagbiasinterp(10**-4), 
        label=r'Nick Data, $\Sigma m_\nu$ = '+f'{Mnu[m_idx]:.2f}', 
        color=colors[m_idx]
    )

    plt.plot(
        ref[:,0], 
        ref[:, 1]/ref_refactor, 
        color=colors[m_idx],  
        marker="o", 
        linestyle="dashed"
    )
    
plt.xlabel(r'$k ~[h ~{\rm Mpc}^{-1}]$', fontsize=15)
plt.ylabel(r'$b_1^L(k)/b_1^L(k_{\rm ref})$', fontsize=15)
plt.title(r'1805.11623 Fig 9, axionCAMB', fontsize=15)
plt.legend(fontsize=15)
plt.grid(True, which='both', axis='both')
plt.savefig("1805.11623_FIG9.png")    
