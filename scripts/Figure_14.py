######################################################
####         USER INPUTS
######################################################

import matplotlib.font_manager
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy 
import seaborn as sns
import subprocess
from matplotlib.lines import Line2D
from matplotlib import rc

#sns.set()
#sns.set_style(style='white')
rc('font', **{'serif': ['Computer Modern']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'
matplotlib.rcParams.update({
    "font.weight" : "bold",
    "font.size" : 110,
    "axes.labelsize" : 110,
    "axes.labelpad" : 8.0,  
    "xtick.labelsize" : 60, 
    "ytick.labelsize" : 60, 
    "xtick.major.size" : 30,
    "xtick.major.width" : 5,
    "xtick.minor.size" : 20,
    "xtick.minor.width" : 3,
    "ytick.major.size" : 30,
    "ytick.major.width" : 5,
    "ytick.minor.size" : 20,
    "ytick.minor.width" : 3,
    "legend.fontsize" : 60, 
    "figure.dpi" : 100, 
    "figure.figsize" : [30, 30],
    "figure.constrained_layout.use" : True, 
    "figure.constrained_layout.wspace": 0.1,
    "savefig.pad_inches" : 0.1

})

######################################################

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast.nosync/"
rfpath_outputsuffix = "output/result-0/"
outpath = "/Users/nicholasdeporzio/Desktop/"

sum_massive_nu = 0.
redshift = 0.65
kmin = 1.0e-5
kmax = 1.1
Nk = 50

omega_cdm_LCDM = 0.1127
omega_b_LCDM = 0.02226
Omega_M_LCDM = 0.27464

######################################################

omega_nu = sum_massive_nu/93.2

kref = np.power(10., -4.) 
omega_ax = np.array([omega_cdm_LCDM*1.0e-12, omega_cdm_LCDM*0.022])
m_ax = np.array([
    np.power(10., -29.),
    np.power(10., -28.),
    np.power(10., -27.)
])

f_cdm_LCDM = omega_cdm_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
f_b_LCDM = omega_b_LCDM/(omega_cdm_LCDM + omega_b_LCDM)

omega_nu = sum_massive_nu/93.14

omega_cdm = omega_cdm_LCDM - omega_ax
f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)

h = 0.70148
#h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm_LCDM+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))
#h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))

######################################################
####         INTERNAL 
######################################################

rfpath_boltzmannsuffix = "Boltzmann_2/transfer_files_0/"

data_pm = []
data_tf = []
data_eulbias = []
data_lagbias = []
axion_background = []
relicfast_pss = []

os.chdir(rfpath+'/include/')
reading_file = open("common.h", "r")
new_file_content = ""
for line in reading_file:
    stripped_line = line.strip()
    if "#define length_transfer_axioncamb " in stripped_line:
        expected_axioncamb_output_lines = int(''.join(filter(str.isdigit, stripped_line)))
    new_line = stripped_line.replace(
            "#define boltzmann_tag  _CLASS_", "#define boltzmann_tag _AXIONCAMB_"
    ).replace(
        "#define boltzmann_tag  _CAMB_", "#define boltzmann_tag _AXIONCAMB_"
    )
        
    new_file_content += "    " + new_line +"\n"
reading_file.close()
writing_file = open("common.h", "w")
writing_file.write(new_file_content)
writing_file.close()
os.chdir(rfpath)
os.system('make')

os.chdir(rfpath)

# Run RelicFast for each axion abundance 
for m_idx, m_val in enumerate(m_ax): 
    for oax_idx, oax_val in enumerate(omega_ax): 
        print("Running RelicFast + axionCAMB for omega_ax = ", oax_val)
        
        # Clear old data 
        os.system('rm -r '+rfpath_outputsuffix)
        os.system('rm -r '+'Boltzmann_0')
        os.system('rm -r '+'Boltzmann_1')
        os.system('rm -r '+'Boltzmann_2')
        os.system('rm ./run.ini')
        
        # RUN RelicFast + solver for each neutrino mass choice 
        reading_file = open("2006.09395.ini", "r")
        new_file_content = ""
        for line in reading_file:
            stripped_line = line.strip()
            new_line = stripped_line.replace(
                "mnu1 = 0.0", "mnu1 = "+f'{sum_massive_nu/3.:.2e}'
            ).replace(
                "mnu2 = 0.0", "mnu2 = "+f'{sum_massive_nu/3.:.2e}'
            ).replace(
                "m_SN = 0.02", "m_SN = "+f'{sum_massive_nu/3.:.2e}'
            ).replace(
                "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{redshift:.3f}'
            ).replace(
                "z_collapse_top = 1.4", "z_collapse_top = 1.65"
            ).replace(
                "N_zcoll = 1", "N_zcoll = 1"
            ).replace(
                "hubble = 0.701", "hubble = "+f'{h:.6f}'
            ).replace(
                "omega_ax = 1.0e-9", "omega_ax = "+f'{oax_val:.6e}'
            ).replace(
                "N_klong = 1", "N_klong = "+str(Nk)
            ).replace(
                "omegac = 0.11271", "omegac = "+f'{omega_cdm[oax_idx]:.6e}'
            ).replace(
                "m_ax = 1.0e-22", "m_ax = "+f'{m_val:.6e}'
            ).replace(
                "kbot = 1e-4", "kbot = "+f'{kmin:.3e}'
            ).replace(
                "ktop = 0.7", "ktop = "+f'{kmax:.3e}'
            )
            
            
            if (sum_massive_nu!=0.): 
                new_line = new_line.replace(
                    "tag_sterile_nu = 0", "tag_sterile_nu = 1"
                )
                
            new_file_content += "    " + new_line +"\n"
        reading_file.close()
        writing_file = open("run.ini", "w")
        writing_file.write(new_file_content)
        writing_file.close()
        
        lines_match_flag=False
        while lines_match_flag==False:
            os.system('./relicfast run.ini')
            with open(rfpath+"/Boltzmann_2/transfer_files_0/_transfer_out_z200.000", 'r') as fp:
                ac_lines_out = sum(1 for line in fp)
            fp.close()
    
            if ac_lines_out==expected_axioncamb_output_lines:
                lines_match_flag=True
            else:
                print("axionCAMB output doesn't match RelicFast compile parameters. Recompiling: "
                    +str(expected_axioncamb_output_lines)+"-->"+str(ac_lines_out))
                os.chdir(rfpath+'/include/')
                reading_file = open("common.h", "r")
                new_file_content = ""
                for line in reading_file:
                    stripped_line = line.strip()
                    new_line = stripped_line.replace(
                            "#define length_transfer_axioncamb "+str(expected_axioncamb_output_lines),
                            "#define length_transfer_axioncamb "+str(ac_lines_out)
                    )
                    new_file_content += "    " + new_line +"\n"
                reading_file.close()
                writing_file = open("common.h", "w")
                writing_file.write(new_file_content)
                writing_file.close()
                os.chdir(rfpath)
                os.system('make')
    
                expected_axioncamb_output_lines = int(ac_lines_out)
        
        # Collect axion background evolution
        axion_background.append(np.loadtxt("/Users/nicholasdeporzio/Downloads/axion_background.dat"))    
        
        # Collect power spectrum for requested redshift 
        output_dirs = os.listdir(rfpath_boltzmannsuffix)
        output_dirs.remove('_params.ini')
        
        z_vals = np.array([])
        
        for str_idx, str_val in enumerate(output_dirs):
            if (str_val[0:15]=='_transfer_out_z'): 
                z_vals = np.append(
                    z_vals, 
                    np.float(str_val.split('_transfer_out_z')[1])
                )
                    
        z_vals = np.sort(z_vals)
            
    
        pm_idx = np.argmin(np.abs(z_vals - redshift))
        print('Requested/found redshift: ', redshift, z_vals[pm_idx])
        print('\t '+rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
        print('\t '+rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
        
        data_pm.append(
            np.loadtxt(rfpath_boltzmannsuffix+'_matterpower_'+str(pm_idx+1)+'.dat')
        )
        
        data_tf.append(
            np.loadtxt(rfpath_boltzmannsuffix+'_transfer_out_z'+f'{z_vals[pm_idx]:.3f}')
        )
        
        data_eulbias.append(
            np.loadtxt(rfpath_outputsuffix+'bias_Euler_z'+f'{z_vals[pm_idx]:.2f}'+'_M13.00_Nk'+str(Nk)+'.dat', skiprows=1)
        )
        data_lagbias.append(
            np.loadtxt(rfpath_outputsuffix+'bias_Lagrangian_z'+f'{z_vals[pm_idx]:.2f}'+'_M13.00_Nk'+str(Nk)+'.dat', skiprows=1)
        )
    
        relicfast_pss.append(np.loadtxt(
            rfpath_outputsuffix+'power_spectra_z'+f"{redshift:.2f}"+"_M13.00_Nk"+f"{Nk:d}"+".dat"
            , skiprows=1)
        )
    
        os.system('mv ./run.ini ./run_'+str(oax_idx)+'.ini')
        os.system('mv ./axionCAMB_Current/params_collapse.ini ./axionCAMB_Current/params_collapse_'+str(oax_idx)+'.ini')    

#####################################

colors = sns.color_palette('magma', len(m_ax))
kplot = np.geomspace(10**-3.0, 1.0, 100)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1, 1]})
for m_idx, m_val in enumerate(m_ax): 
    for oax_idx, oax_val in enumerate(omega_ax): 
        k_vals = relicfast_pss[m_idx*len(omega_ax)+oax_idx][:, 0]
        pmm_vals = relicfast_pss[m_idx*len(omega_ax)+oax_idx][:, 1]
        phh_vals = relicfast_pss[m_idx*len(omega_ax)+oax_idx][:, 3]
        
        pmm_interp = scipy.interpolate.interp1d(k_vals, pmm_vals)
        phh_interp = scipy.interpolate.interp1d(k_vals, phh_vals)
    
        if oax_idx==0:
            pmm_LCDM = pmm_interp(kplot)
            phh_LCDM = phh_interp(kplot)
            pmm_LCDM_ref = pmm_interp(kref)
            phh_LCDM_ref = phh_interp(kref) 
            print(pmm_LCDM)
            print(phh_LCDM)
    
        else:
            pmm = pmm_interp(kplot)
            phh = phh_interp(kplot) 
            Rmm = pmm/pmm_LCDM
            Rhh = phh/phh_LCDM
            Rmm_ref = pmm_interp(kref)/pmm_LCDM_ref
            Rhh_ref = phh_interp(kref)/phh_LCDM_ref
    
            yplot1 = Rmm/Rmm_ref
            yplot2 = Rhh/Rhh_ref    
     
            ax1.plot(#Matter
                kplot,
                yplot1, 
                color=colors[m_idx], 
                linewidth=5., 
                linestyle='dashed'
            )
            ax1.plot(#Halo 
                kplot,
                yplot2, 
                label=(r"$m_\phi = 10^{"+f"{np.log10(m_val):.0f}"+r"} {\rm ~eV}$"), 
                color=colors[m_idx], 
                linewidth=5.,
                linestyle='solid'
            )
            ax2.plot(
                kplot, 
                yplot2-yplot1, 
                color=colors[m_idx], 
                linewidth=5., 
                linestyle='dotted'
            )
            ax3.plot(
                kplot, 
                (yplot2-yplot1)/yplot1, 
                color=colors[m_idx], 
                linewidth=5., 
                linestyle='dotted'
            )
ax1.axvspan(0.025*h, 0.23*h, alpha=0.5, color='grey')
#ax.plot([min(kplot), max(kplot)], [1., 1.], color='black', linewidth=5.)
ax1.set_xlim((min(kplot), max(kplot)))
#ax.plot([0.7*0.015, 0.7*0.015], [0.999, 1.008], color='red', label=r'$k_{eq}$')
#ax.plot([0.024, 0.024], [0.999, 1.008], color='blue', label=r'$k_{*}$')
ax1.set_xscale('log')
ax1.set_ylabel(r'$R(k)/R(k_{\rm ref})$')
ax1.set_ylim((0.69, 1.01))
ax1.tick_params(axis='both')
ax1.legend()
ax1.grid(False)
ax2.axvspan(0.025*h, 0.23*h, alpha=0.5, color='grey')
ax2.set_ylabel(r'$\Delta_{\rm H-M}$', fontsize=30)
ax3.axvspan(0.025*h, 0.23*h, alpha=0.5, color='grey')
ax3.set_ylabel(r'$\bar{\Delta}_{\rm (H-M)/M}$', fontsize=30)
ax3.set_xlabel(r'$k ~[{\rm Mpc}^{-1}]$')
plt.savefig(rfpath+"plots/Figure_14.png")

#fig, ax = plt.subplots(1, 1)
#for oax_idx, oax_val in enumerate(omega_ax): 
#    k_vals = relicfast_pss[oax_idx][:, 0]
#    pmm_vals = relicfast_pss[oax_idx][:, 1]
#    phh_vals = relicfast_pss[oax_idx][:, 3]
#    
#    pmm_interp = scipy.interpolate.interp1d(k_vals, pmm_vals)
#    phh_interp = scipy.interpolate.interp1d(k_vals, phh_vals)
#
#    if oax_idx==0:
#        pmm_LCDM = pmm_interp(kplot)
#        phh_LCDM = phh_interp(kplot) 
#        print(pmm_LCDM)
#        print(phh_LCDM)
#
#    else:
#        pmm = pmm_interp(kplot)
#        phh = phh_interp(kplot) 
#        print(pmm)
#        print(phh)
#        yplot1 = pmm/pmm_LCDM
#        yplot2 = phh/phh_LCDM    
#
#        ax.plot(
#            kplot,
#            yplot2, 
#            label=r"$R_h(k)$", 
#            color=colors[oax_idx-1], 
#            linewidth=5.,
#            linestyle='dashed'
#        )
##ax.plot([min(kplot), max(kplot)], [1., 1.], color='black', linewidth=5.)
#ax.set_xlim((min(kplot), max(kplot)))
##ax.plot([0.7*0.015, 0.7*0.015], [0.999, 1.008], color='red', label=r'$k_{eq}$')
##ax.plot([0.024, 0.024], [0.999, 1.008], color='blue', label=r'$k_{*}$')
#ax.set_xscale('log')
#ax.set_xlabel(r'$k ~[{\rm Mpc}^{-1}]$')
#ax.set_ylabel(r'$R(k)$')
#ax.tick_params(axis='both')
#ax.legend()
#ax.grid(False)
#ax.set_title(r"$m_\phi = 10^{"+f"{np.log10(m_ax):.1f}"+r"}$ eV")
#plt.savefig(rfpath+"plots/Figure_8b_logmax"+f"{np.log10(m_ax):.1f}"+".png")
