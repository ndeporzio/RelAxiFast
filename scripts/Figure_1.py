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
    "axes.labelpad" : 10.0,  
    "xtick.labelsize" : 60, 
    "xtick.major.size" : 30,
    "xtick.major.width" : 5,
    "xtick.major.pad" : 15, 
    "xtick.minor.size" : 20,
    "xtick.minor.width" : 3,
    "xtick.direction" : "out", 
    "ytick.labelsize" : 60,
    "ytick.major.size" : 30,
    "ytick.major.width" : 5,
    "ytick.major.pad" : 10,   
    "ytick.minor.size" : 20,
    "ytick.minor.width" : 3,
    "ytick.direction" : "out", 
    "legend.fontsize" : 60, 
    "figure.dpi" : 100, 
    "figure.figsize" : [30, 30],
    "figure.constrained_layout.use" : True, 
    "figure.constrained_layout.wspace": 0.1,
    "savefig.pad_inches" : 0.1

})
use_existing_data=True
data_save_level=2

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelicFast.nosync/"
acpath = (rfpath+"axionCAMB_Current/")
acpath_output = "Boltzmann_2/transfer_files_0/"

M_ax = np.logspace(-32.5, -22.0, 106)
M_nu = 0. # Units: eV
redshift = 0.7
kmin = 5.0e-5
kmax = 1.5
Nk = 100
M_halo = 1.0e13

h_lcdm = 0.67
omega_nu = M_nu/93.2
omega_cdm_LCDM = 0.12
omega_b_LCDM = 0.022
Omega_M_LCDM = (omega_b_LCDM + omega_cdm_LCDM)/np.power(h_lcdm, 2.)
print(Omega_M_LCDM)
Omega_ax = (
    np.array([1.0e-9, 0.01, 0.05, 0.1])
    *omega_cdm_LCDM
    /np.power(h_lcdm, 2.)
)

######################################################
####         INTERNAL
######################################################

omega_ax = Omega_ax*np.power(h_lcdm, 2.)
omega_cdm = omega_cdm_LCDM - omega_ax

f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)

axion_aosc = np.zeros((len(M_ax), len(omega_ax)))
axion_kfs = np.zeros((len(M_ax), len(omega_ax)))

# Try to load data 
data_loaded=False
if (use_existing_data==True): 
    m_path = (rfpath+"plots/Figure_1_m_ax.txt")
    omega_path = (rfpath+"plots/Figure_1_omega_ax.txt")
    aosc_path = (rfpath+"plots/Figure_1_axion_aosc.txt")
    zosc_path = (rfpath+"plots/Figure_1_axion_zosc.txt")
    kfs_path = (rfpath+"plots/Figure_1_axion_kfs.txt")
    if (
        os.path.exists(m_path) 
        and os.path.exists(omega_path) 
        and os.path.exists(aosc_path)
        and os.path.exists(zosc_path) 
        and os.path.exists(kfs_path)
    ):
        print("Loading existing data...")
        M_ax = np.loadtxt(m_path)
        omega_ax = np.loadtxt(omega_path)
        axion_aosc = np.loadtxt(aosc_path)
        axion_zosc = np.loadtxt(zosc_path)
        axion_kfs = np.loadtxt(kfs_path)
        
        data_loaded=True

# Calculate Pmm using axionCAMB
if data_loaded==False: 
    os.chdir(acpath)
    for m_idx, m_val in enumerate(M_ax):
        print("Running axionCAMB for m_ax = ", m_val)
        for o_idx, o_val in enumerate(omega_ax):
            print("\t For omega_ax = ", o_val)
    
            axion_kfs[m_idx, o_idx] = (
                np.pi 
                * np.sqrt(m_val*1.56*np.power(10., 29)) 
                * np.power((h_lcdm/2997.)*np.power(1.+redshift, 3.), 0.5)
            ) 
            # axionCAMB for each mass/abundance choice
            reading_file = open(acpath+"params_base.ini", "r")
            new_file_content = ""
            for line in reading_file:
                stripped_line = line.strip()
                new_line = str(stripped_line)
                #new_line = stripped_line.replace("mnu1 = 0.0", "mnu1 = "+f'{M_nu:.2e}')
                #if (M_nu!=0.0):
                #    new_line = new_line.replace(
                #        "tag_sterile_nu = 0", "tag_sterile_nu = 1"
                #    )
                new_file_content += "    " + new_line +"\n"
            new_file_content += "    " + "ombh2 = " + f"{omega_b_LCDM:.4f}" +"\n"
            new_file_content += "    " + "omch2 = " + f"{omega_cdm[o_idx]:.4f}" +"\n"
            new_file_content += "    " + "omnuh2 = " + f"{omega_nu:.4f}" +"\n"
            new_file_content += "    " + "hubble = " + f"{100.*h_lcdm:.4f}" +"\n"
            new_file_content += "    " + "omaxh2 = " + f"{o_val:.4e}" +"\n"
            new_file_content += "    " + "m_ax = " + f"{m_val:.4e}" +"\n"
            new_file_content += "    " + "massless_neutrinos = 3.046"+"\n"
            new_file_content += "    " + "nu_mass_eigenstates = 1"+"\n"
            new_file_content += "    " + "massive_neutrinos = 0"+"\n"
            new_file_content += "    " + "scalar_amp(1) = 2.20e-09"+"\n"
            new_file_content += "    " + "scalar_spectral_index(1) = 0.9655"+"\n"
            new_file_content += "    " + "transfer_num_redshifts = 1"+"\n"
            new_file_content += "    " + "transfer_filename(1) = transfer_out_z" + f"{redshift:.3f}" +"\n"
            new_file_content += "    " + "transfer_redshift(1) = " + f"{redshift:.3f}" +"\n"
            new_file_content += "    " + "output_root = Boltzmann_2/transfer_files_0/"+"\n"
    
            reading_file.close()
            writing_file = open(acpath+"params_collapse_"+f"{m_idx+o_idx:d}"+".ini", "w")
            writing_file.write(new_file_content)
            writing_file.close()
    
            os.system('./camb '+"params_collapse_"+f"{m_idx+o_idx:d}"+".ini")
    
            axion_aosc[m_idx, o_idx] = np.loadtxt("/Users/nicholasdeporzio/Downloads/axion_aosc.dat")
    
    axion_zosc = ((1./axion_aosc)-1.)

    if data_save_level>0:     
        np.savetxt(rfpath+"plots/Figure_1_m_ax.txt", M_ax)
        np.savetxt(rfpath+"plots/Figure_1_omega_ax.txt", omega_ax)
        np.savetxt(rfpath+"plots/Figure_1_axion_aosc.txt", axion_aosc)
        np.savetxt(rfpath+"plots/Figure_1_axion_zosc.txt", axion_zosc)
        np.savetxt(rfpath+"plots/Figure_1_axion_kfs.txt", axion_kfs)


plot_x = np.geomspace(np.min(M_ax), np.max(M_ax), 100) #Units: [h Mpc^-1]
colors = sns.color_palette("magma", len(omega_ax))
fig, ax = plt.subplots(1, 1)
for o_idx, o_val in enumerate(omega_ax[0:1]):
    osc_interp = scipy.interpolate.interp1d(np.log10(M_ax), np.log10(axion_zosc[:,o_idx]))
    plot_y = np.power(10., osc_interp(np.log10(plot_x)))
    ax.plot(
        plot_x,
        plot_y,
        #M_ax, 
        #axion_zosc[:,o_idx],
        label=r'$\Omega_{\phi}/\Omega_{d} = '+f'{omega_ax[o_idx]/omega_cdm_LCDM:.2f}'+r"$",
        #color=colors[o_idx],
        color="black",
        linewidth=5.0
    )
dlogydlogx1 = 0.66
b1 = osc_interp(-30.)-dlogydlogx1*(-30.)  
dlogydlogx2 = 0.50
b2 = osc_interp(-22.)-dlogydlogx2*(-22.) 
ax.plot(plot_x, np.power(10., np.log10(plot_x)*dlogydlogx1 + b1), linestyle='dashed', color='tab:orange', linewidth=5)
ax.plot(plot_x, np.power(10., np.log10(plot_x)*dlogydlogx2 + b2), linestyle='dashed', color='tab:cyan', linewidth=5)
ax.set_yscale('log')
ax.set_ylabel(r'$z_{\rm osc}$')
ax.set_yticks(np.logspace(-1, 6, 8))
ax.set_xscale('log')
ax.set_xlabel(r'$m_\phi ~{\rm [eV]}$')
ax.set_xticks(np.logspace(-32, -22, 6))
#ax.grid(False)
ax.set_ylim((0.09, 2.0e6))
ax.set_xlim((1.0e-33, 1.0e-21))

ax.fill_between([1.0e-33, 1.0e-21], [0.09, 0.09], [200., 200.], color='grey', alpha=0.5)
ax.fill_between([1.0e-33, 1.0e-21], [1100, 1100], [1.3*np.power(10., 6.), 1.3*np.power(10., 6.)], color='red', alpha=0.1)
#ax.plot([1.0e-33, 1.0e-21], [1100, 1100], color='red', alpha=0.5, linewidth=5)
#ax.plot([1.0e-33, 1.0e-21], [1.3*np.power(10., 6.), 1.3*np.power(10., 6.)], color='red', alpha=0.5, linewidth=5, linestyle='dashed')

ax_twin = ax.twinx()

#ax_twin.plot((plot_x), (1./(plot_y+1.)), color='cyan', linestyle='dashed')
#ax_twin.set_xscale('log')
#ax_twin.set_xticks(np.logspace(-32, -22, 6))
ax_twin.set_yscale('log')
#ax_twin.grid(False)
#ax_twin.tick_params(axis='both')
#ax_twin.tick_params(axis='x', direction="in", length=50, width=5)
ax_twin.set_ylabel(r"${a}_{\rm osc}$")
ax_twin.set_ylim(ax.get_ylim())
ax_twin.set_yticks(np.concatenate((
    [np.power(10., -1.*np.log10(0.9)), np.power(10., -1.*np.log10(0.5))], np.logspace(1., 6., 6)))-1.)
ax_twin.set_yticklabels(
    [r'$0.9$', r'$0.5$', r'$10^{-1}$', r'$10^{-2}$', r'$10^{-3}$', r'$10^{-4}$', r'$10^{-5}$', r'$10^{-6}$'])
#ax_twin.set_yticks((1./np.logspace(-1., 6., 8))-1.)
#ax_twin.set_yticklabels(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])
#ax.legend(fontsize=25)
plt.savefig(rfpath+"plots/Figure_1.png")

######################################################
####         EXTRA PLOTS
######################################################
plot_x = np.geomspace(np.min(M_ax), np.max(M_ax), 100) #Units: [h Mpc^-1]
colors = sns.color_palette("magma", len(omega_ax))
fig, ax = plt.subplots(1, 1)
for o_idx, o_val in enumerate(omega_ax):
    osc_interp = scipy.interpolate.interp1d(np.log10(M_ax), np.log10(axion_zosc[:,o_idx]))
    plot_y = np.power(10., osc_interp(np.log10(plot_x)))
    if o_idx==0: 
        norm = np.array(plot_y)
        continue
    ax.plot(
        plot_x,
        plot_y/norm-1.,
        #M_ax, 
        #axion_zosc[:,o_idx],
        label=r'$\Omega_{\phi}/\Omega_{d} = '+f'{omega_ax[o_idx]/omega_cdm_LCDM:.2f}'+r"$",
        color=colors[o_idx],
        linewidth=1.0
    )
ax.set_xscale('log')
ax.set_xlabel(r'$M_\phi {\rm ~[eV]}$')
ax.set_ylabel(r'$z_{\rm osc}$')
ax.grid(False)
ax.tick_params(axis='both')
plt.savefig(rfpath+"plots/Figure_1b.png")


fig, ax = plt.subplots(1, 1)
for o_idx, o_val in enumerate(omega_ax[0:1]):
    osc_interp = scipy.interpolate.interp1d(np.log10(M_ax), np.log10(axion_zosc[:,o_idx]))
    plot_y = np.diff(osc_interp(np.log10(plot_x)))/np.diff(np.log10(plot_x))
    ax.plot(
        plot_x[0:-1],
        plot_y,
        color=colors[o_idx],
        linewidth=5.0
    )
ax.set_xlabel(r'$\Delta \log{(M_\phi/[{\rm eV}])}$')
ax.set_xscale('log') 
ax.set_ylabel(r'$\Delta \log{(z_{\rm osc})}$')
ax.grid(True, which='both', axis='both')
ax.tick_params(axis='both')
plt.savefig(rfpath+"plots/Figure_1c.png")
