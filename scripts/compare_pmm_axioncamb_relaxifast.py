import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
import seaborn as sns
import subprocess

from matplotlib.lines import Line2D

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

omega_ax = Omega_ax*np.power(h_lcdm, 2.)
omega_cdm = omega_cdm_LCDM - omega_ax

f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)

axion_aosc = np.zeros((len(M_ax), len(omega_ax)))
axion_kfs = np.zeros((len(M_ax), len(omega_ax)))

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

