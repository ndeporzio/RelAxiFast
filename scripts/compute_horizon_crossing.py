import matplotlib.font_manager
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
import seaborn as sns
import subprocess
from matplotlib.lines import Line2D
from matplotlib import rc

rc('font', **{'serif': ['Computer Modern']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'
matplotlib.rcParams.update({
    "font.weight" : "bold",
    "font.size" : 60,
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
    "figure.figsize" : [30, 20],
    "figure.constrained_layout.use" : True,
    "figure.constrained_layout.wspace": 0.1,
    "savefig.pad_inches" : 0.1

})

rfpath = "/Users/nicholasdeporzio/Documents/Academic/Projects/P005_FuzzyCdmBias/RelAxiFast.nosync"

rfpath_outputsuffix = "output/result-0/"

outpath = "/Users/nicholasdeporzio/Desktop/"

omega_cdm_LCDM = 0.1127
omega_b_LCDM = 0.02226
Omega_M_LCDM = 0.27464

m_ax = np.power(10., -32.5)
omega_ax = 0.1*omega_cdm_LCDM

sum_massive_nu = 0.
redshift = 0.65
kmin = 5.0e-5
kmax = 1.5
Nk = 50

omega_nu = sum_massive_nu/93.2
Mnu = sum_massive_nu # Units: eV

f_cdm_LCDM = omega_cdm_LCDM/(omega_cdm_LCDM + omega_b_LCDM)
f_b_LCDM = omega_b_LCDM/(omega_cdm_LCDM + omega_b_LCDM)

omega_nu = Mnu/93.14

omega_cdm = omega_cdm_LCDM - omega_ax
f_cdm = omega_cdm/(omega_cdm + omega_b_LCDM)
f_b = omega_b_LCDM/(omega_cdm + omega_b_LCDM)

h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm_LCDM+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))
#h = (0.70148)*np.sqrt((omega_b_LCDM+omega_cdm+omega_nu+omega_ax)/(omega_b_LCDM+omega_cdm_LCDM))

kfs = np.pi * np.sqrt(m_ax*1.56*np.power(10., 29)) * np.power((h/2997.)*np.power(1.+redshift, 3.), 0.5)
keq = 0.70148*0.015

rfpath_boltzmannsuffix = "Boltzmann_2/transfer_files_0/"

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


print("Running RelicFast + axionCAMB for m_ax = ", m_ax)

# Clear old data 
os.system('rm -r '+rfpath_outputsuffix)
os.system('rm -r '+'Boltzmann_0')
os.system('rm -r '+'Boltzmann_1')
os.system('rm -r '+'Boltzmann_2')
os.system('rm ./run.ini')

reading_file = open("2006.09395.ini", "r")
new_file_content = ""
for line in reading_file:
    stripped_line = line.strip()
    new_line = stripped_line.replace(
        "mnu1 = 0.0", "mnu1 = "+f'{Mnu:.2e}'
    ).replace(
        "mnu2 = 0.0", "mnu2 = "+f'{Mnu:.2e}'
    ).replace(
        "m_SN = 0.02", "m_SN = "+f'{Mnu:.2e}'
    ).replace(
        "z_collapse_bot = 0.7", "z_collapse_bot = "+f'{redshift:.3f}'
    ).replace(
        "z_collapse_top = 1.4", "z_collapse_top = 1.65"
    ).replace(
        "N_zcoll = 1", "N_zcoll = 1"
    ).replace(
        "kbot = 1e-4", "kbot = "+f'{kmin:.3e}'
    ).replace(
        "ktop = 0.7", "ktop = "+f'{kmax:.3e}'
    ).replace(
        "hubble = 0.701", "hubble = "+f'{h:.6f}'
    ).replace(
        "omega_ax = 1.0e-9", "omega_ax = "+f'{omega_ax:.6e}'
    ).replace(
        "N_klong = 1", "N_klong = "+str(Nk)
    ).replace(
        "omegac = 0.11271", "omegac = "+f'{omega_cdm:.6e}'
    ).replace(
        "m_ax = 1.0e-22", "m_ax = "+f'{m_ax:.6e}'
    )

    if (Mnu!=0.):
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
    with open(rfpath+"/Boltzmann_2/transfer_files_0/_transfer_out_z0.000", 'r') as fp:
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

#Collect dimensionless H defined as H/(100 km/s/Mpc)....
Hdimensionless = np.loadtxt(rfpath+"/axion_cad2.dat", skiprows=0)[:, [0, 2]] 
m_ax_osc = np.loadtxt(rfpath+"/plots/Figure_1_m_ax.txt", skiprows=0)
a_osc = np.loadtxt(rfpath+"/plots/Figure_1_axion_aosc.txt", skiprows=0)[:,2] #5% abundance (though all similar)
a_eq = np.float(np.loadtxt(rfpath+"/axion_aeq.dat")) 
c = 3.0e8 #Units of m/s 

#Convert H to units of (m/s/Mpc) 
H = np.matmul(Hdimensionless, [[1, 0],[0, 100.*1000.]]) #rescale second column of Hdimensionless 

#Compute comoving horizon in units of Mpc
amin = np.min(H[:,0]) 
amax = np.max(H[:,0]) 
N_a = 1000

Xcomoving = np.zeros((N_a, 2))

Xcomoving_ini = (c*amin)/(100.*1000.*h*np.sqrt(Omega_M_LCDM*a_eq))
print("Initial comoving horizon at a_ini=", amin, ": ", Xcomoving_ini) 

aplot = np.geomspace(amin, amax, N_a)
Xintegrand = c/(np.power(H[:,0], 2.)*H[:,1])
Xintegrand_interp = scipy.interpolate.interp1d(H[:,0], Xintegrand, kind="linear")
integrand = lambda a: Xintegrand_interp(a)
for idx in range(N_a): 
    print("Evaluate comoving distance #", idx+1, " of ", N_a) 
    a = aplot[idx]
    Xcomoving[idx, 0] = a
    Xcomoving[idx, 1] = scipy.integrate.quad(integrand, amin, a)[0] + Xcomoving_ini

np.savetxt("Comoving_Horizon.txt", Xcomoving) 

 


 
 



