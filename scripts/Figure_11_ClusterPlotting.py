#############################

# TO BE RUN FROM DIRECTORY CONTAINING CLUSTER OUTPUTS
# OF ONE DIRECTORY FOR EACH MASS/ABUNDANCE CHOICE


import matplotlib.font_manager
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
import seaborn as sns
import subprocess
from matplotlib.lines import Line2D
from matplotlib import rc

omega_cdm_LCDM = 0.1127
m_ax = np.loadtxt("./m_ax.txt")
omega_ax = np.loadtxt("./omega_ax.txt") 
nruns = len(m_ax)*len(omega_ax) 

#sns.set()
#sns.set_style(style='white')
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
    "figure.figsize" : [30, 30],
    "figure.constrained_layout.use" : True,
    "figure.constrained_layout.wspace": 0.1,
    "savefig.pad_inches" : 0.1

})

# Scales to compute bias values at 
kref = np.power(10., -0.5) # Units: 

b1e = np.zeros((len(m_ax), len(omega_ax))) # Units: none
b1l = np.zeros((len(m_ax), len(omega_ax))) # Units: none
deltacritlow = np.zeros((len(m_ax), len(omega_ax))) # Units: none
deltainitlow = np.zeros((len(m_ax), len(omega_ax))) # Units: none
deltacrithigh = np.zeros((len(m_ax), len(omega_ax))) # Units: none
deltainithigh = np.zeros((len(m_ax), len(omega_ax))) # Units: none
deltalonglow = np.zeros((len(m_ax), len(omega_ax))) # Units: none
deltalonghigh = np.zeros((len(m_ax), len(omega_ax))) # Units: none

b1e_ref = np.zeros((len(m_ax), len(omega_ax))) # Units: none
b1l_ref = np.zeros((len(m_ax), len(omega_ax))) # Units: none
deltacritlow_ref = np.zeros((len(m_ax), len(omega_ax))) # Units: none
deltainitlow_ref = np.zeros((len(m_ax), len(omega_ax))) # Units: none
deltacrithigh_ref = np.zeros((len(m_ax), len(omega_ax))) # Units: none
deltainithigh_ref = np.zeros((len(m_ax), len(omega_ax))) # Units: none
deltalonglow_ref = np.zeros((len(m_ax), len(omega_ax))) # Units: none
deltalonghigh_ref = np.zeros((len(m_ax), len(omega_ax))) # Units: none

sigmaMinit = np.zeros((len(m_ax), len(omega_ax))) # Units: none
sigmaMcoll = np.zeros((len(m_ax), len(omega_ax))) # Units: none
dsigmaMinit =np.zeros((len(m_ax), len(omega_ax))) # Units: none
deltacritsim = np.zeros((len(m_ax), len(omega_ax))) # Units: none
b1e_step = np.zeros((len(m_ax), len(omega_ax))) # Units: none
b1l_step = np.zeros((len(m_ax), len(omega_ax))) # Units: none
hmf = np.zeros((len(m_ax), len(omega_ax))) # Units: none

def fmt(x):
    s = f"{(x-1.)*100.:.0f}"
    return rf"${s} \%$" if plt.rcParams["text.usetex"] else f"{s} %"

for idx in range(nruns):
    m_idx = np.mod(idx, len(m_ax)) 
    o_idx = np.floor_divide(idx, len(m_ax))

    b1ldata = np.loadtxt("./FIG11_"+str(idx+1)+"/bias_Lagrangian_z0.65_M13.00_Nk50.dat", skiprows=1)
    b1linterp = scipy.interpolate.interp1d(np.log10(b1ldata[:,0]), b1ldata[:,1])
    b1l[m_idx, o_idx] = b1linterp(np.log10(kref))
    b1l_ref[m_idx, o_idx] = b1linterp(-4.) 
    b1l_step[m_idx, o_idx] = b1l[m_idx, o_idx]/b1l_ref[m_idx, o_idx] 

    b1edata = np.loadtxt("./FIG11_"+str(idx+1)+"/bias_Euler_z0.65_M13.00_Nk50.dat", skiprows=1)
    b1einterp = scipy.interpolate.interp1d(np.log10(b1edata[:,0]), b1edata[:,1])
    b1e[m_idx, o_idx] = b1einterp(np.log10(kref))
    b1e_ref[m_idx, o_idx] = b1einterp(-4.) 
    b1e_step[m_idx, o_idx] = b1e[m_idx, o_idx]/b1e_ref[m_idx, o_idx] 

    deltacritdatalow = np.loadtxt("./FIG11_"+str(idx+1)+"/delta_crit_z0.65_M13.00_Nk50.dat", skiprows=1)[0:50,[1,2]]
    deltacritinterplow = scipy.interpolate.interp1d(np.log10(deltacritdatalow[:,0]), deltacritdatalow[:,1])
    deltacritlow[m_idx, o_idx] = deltacritinterplow(np.log10(kref))
    deltacritlow_ref[m_idx, o_idx] = deltacritinterplow(-4.)

    deltainitdatalow = np.loadtxt("./FIG11_"+str(idx+1)+"/delta_initial_z0.65_M13.00_Nk50.dat", skiprows=1)[0:50,[1,2]]
    deltainitinterplow = scipy.interpolate.interp1d(np.log10(deltainitdatalow[:,0]), deltainitdatalow[:,1])
    deltainitlow[m_idx, o_idx] = deltainitinterplow(np.log10(kref)) 
    deltainitlow_ref[m_idx, o_idx] = deltainitinterplow(-4.) 

    deltacritdatahigh = np.loadtxt("./FIG11_"+str(idx+1)+"/delta_crit_z0.65_M13.00_Nk50.dat", skiprows=1)[50:100,[1,2]]
    deltacritinterphigh = scipy.interpolate.interp1d(np.log10(deltacritdatahigh[:,0]), deltacritdatahigh[:,1])
    deltacrithigh[m_idx, o_idx] = deltacritinterphigh(np.log10(kref))
    deltacrithigh_ref[m_idx, o_idx] = deltacritinterphigh(-4.)

    deltainitdatahigh = np.loadtxt("./FIG11_"+str(idx+1)+"/delta_initial_z0.65_M13.00_Nk50.dat", skiprows=1)[50:100,[1,2]]
    deltainitinterphigh = scipy.interpolate.interp1d(np.log10(deltainitdatahigh[:,0]), deltainitdatahigh[:,1])
    deltainithigh[m_idx, o_idx] = deltainitinterphigh(np.log10(kref)) 
    deltainithigh_ref[m_idx, o_idx] = deltainitinterphigh(-4.) 

    deltalonglow[m_idx, o_idx] = 0.

    deltalongcollapsedata = np.loadtxt("./FIG11_"+str(idx+1)+"/delta_crit_z0.65_M13.00_Nk50.dat", skiprows=1)[50:100,[0,1]]
    deltalongcollapseinterp = scipy.interpolate.interp1d(np.log10(deltalongcollapsedata[:,1]), deltalongcollapsedata[:,0])
    deltalonghigh[m_idx, o_idx] = deltalongcollapseinterp(np.log10(kref)) 
    deltalonghigh_ref[m_idx, o_idx] = deltalongcollapseinterp(-4.) 

    sigmaM = np.loadtxt("./FIG11_"+str(idx+1)+"/sigmaM_z0.65_M13.00_Nk50.dat", skiprows=1); 
    sigmaMinit[m_idx, o_idx] = sigmaM[0]         
    sigmaMcoll[m_idx, o_idx] = sigmaM[1]
    dsigmaMinit[m_idx, o_idx] = sigmaM[4]

    deltacritsim[m_idx, o_idx] = np.float(np.loadtxt("./FIG11_"+str(idx+1)+"/delta_crit_sim_z0.65_M13.00_Nk50.dat", skiprows=0))
    hmf[m_idx, o_idx] = np.float(np.loadtxt("./FIG11_"+str(idx+1)+"/HMF_z0.65_M13.00_Nk50.dat", skiprows=0))


#########################################
#########################################
#########################################


Z = np.transpose(np.nan_to_num(b1l))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), shading="auto", cmap='magma')
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{m_\phi / {\rm [eV]}}$")
ax.set_ylabel("$\omega_\phi / \omega_d$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$b^1_L(k)$")#, size=50
)
plt.savefig("./Figure_11_b1L.png")
np.savetxt("Figure_11_b1L.txt", Z)

Z = np.transpose(np.nan_to_num(b1l_ref))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), shading="auto", cmap='magma')
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{m_\phi / {\rm [eV]}}$")
ax.set_ylabel("$\omega_\phi / \omega_d$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$b^1_L(k)$")#, size=50
)
plt.savefig("./Figure_11_b1L_ref.png")
np.savetxt("Figure_11_b1L_ref.txt", Z)

Z = np.transpose(np.nan_to_num(b1l_step))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=1.0, vmax=np.max(b1l_step), shading="auto", cmap='magma')
if (np.max(Z)>1.01):
    CS = ax.contour(X, Y, Z, np.linspace(1.01, 1.05, 5), colors='white')
    ax.clabel(CS, CS.levels, inline=False, fmt=fmt)
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$b^1_L(k)~/~b^1_L(k_{\rm ref})$")#, size=50
)
plt.savefig("Figure_11_b1Lstep.png")
np.savetxt("Figure_11_b1Lstep.txt", Z)

Z = np.transpose(np.nan_to_num(deltacritlow))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(deltacritlow), vmax=np.max(deltacritlow), shading="auto", cmap='magma')
if (np.max(Z)>1.01):
    CS = ax.contour(X, Y, Z, np.linspace(1.01, 1.05, 5), colors='white')
    ax.clabel(CS, CS.levels, inline=False, fmt=fmt)
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\delta_{crit, low}(k)$")#, size=50
)
plt.savefig("./Figure_11_deltacritlow.png")
np.savetxt("Figure_11_deltacritlow.txt", Z)

Z = np.transpose(np.nan_to_num(deltacritlow_ref))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(deltacritlow), vmax=np.max(deltacritlow), shading="auto", cmap='magma')
if (np.max(Z)>1.01):
    CS = ax.contour(X, Y, Z, np.linspace(1.01, 1.05, 5), colors='white')
    ax.clabel(CS, CS.levels, inline=False, fmt=fmt)
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\delta_{crit, low}(k)$")#, size=50
)
plt.savefig("./Figure_11_deltacritlow_ref.png")
np.savetxt("Figure_11_deltacritlow_ref.txt", Z)

Z = np.transpose(np.nan_to_num(deltacrithigh))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(deltacrithigh), vmax=np.max(deltacrithigh), shading="auto", cmap='magma')
if (np.max(Z)>1.01):
    CS = ax.contour(X, Y, Z, np.linspace(1.01, 1.05, 5), colors='white')
    ax.clabel(CS, CS.levels, inline=False, fmt=fmt)
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\delta_{crit, high}(k)$")#, size=50
)
plt.savefig("./Figure_11_deltacrithigh.png")
np.savetxt("Figure_11_deltacrithigh.txt", Z)

Z = np.transpose(np.nan_to_num(deltacrithigh_ref))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(deltacrithigh), vmax=np.max(deltacrithigh), shading="auto", cmap='magma')
if (np.max(Z)>1.01):
    CS = ax.contour(X, Y, Z, np.linspace(1.01, 1.05, 5), colors='white')
    ax.clabel(CS, CS.levels, inline=False, fmt=fmt)
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\delta_{crit, high}(k)$")#, size=50
)
plt.savefig("./Figure_11_deltacrithigh_ref.png")
np.savetxt("Figure_11_deltacrithigh_ref.txt", Z)

Z = np.transpose(np.nan_to_num(deltainitlow))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(deltainitlow), vmax=np.max(deltainitlow), shading="auto", cmap='magma')
if (np.max(Z)>1.01):
    CS = ax.contour(X, Y, Z, np.linspace(1.01, 1.05, 5), colors='white')
    ax.clabel(CS, CS.levels, inline=False, fmt=fmt)
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\delta_{init, low}(k)$")#, size=50
)
plt.savefig("./Figure_11_deltainitlow.png")
np.savetxt("Figure_11_deltainitlow.txt", Z)  

Z = np.transpose(np.nan_to_num(deltainitlow_ref))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(deltainitlow), vmax=np.max(deltainitlow), shading="auto", cmap='magma')
if (np.max(Z)>1.01):
    CS = ax.contour(X, Y, Z, np.linspace(1.01, 1.05, 5), colors='white')
    ax.clabel(CS, CS.levels, inline=False, fmt=fmt)
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\delta_{init, low}(k)$")#, size=50
)
plt.savefig("./Figure_11_deltainitlow_ref.png")
np.savetxt("Figure_11_deltainitlow_ref.txt", Z)  

Z = np.transpose(np.nan_to_num(deltainithigh))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(deltainithigh), vmax=np.max(deltainithigh), shading="auto", cmap='magma')
if (np.max(Z)>1.01):
    CS = ax.contour(X, Y, Z, np.linspace(1.01, 1.05, 5), colors='white')
    ax.clabel(CS, CS.levels, inline=False, fmt=fmt)
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\delta_{init, high}(k)$")#, size=50
)
plt.savefig("./Figure_11_deltainithigh.png")
np.savetxt("Figure_11_deltainithigh.txt", Z)  

Z = np.transpose(np.nan_to_num(deltainithigh_ref))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(deltainithigh), vmax=np.max(deltainithigh), shading="auto", cmap='magma')
if (np.max(Z)>1.01):
    CS = ax.contour(X, Y, Z, np.linspace(1.01, 1.05, 5), colors='white')
    ax.clabel(CS, CS.levels, inline=False, fmt=fmt)
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\delta_{init, high}(k)$")#, size=50
)
plt.savefig("./Figure_11_deltainithigh_ref.png")
np.savetxt("Figure_11_deltainithigh_ref.txt", Z)  

Z = np.transpose(np.nan_to_num(sigmaMinit))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), shading="auto", cmap='magma')
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\sigma_M(z_{ini})$")#, size=50
)
plt.savefig("./Figure_11_sigmaMinit.png")
np.savetxt("Figure_11_sigmaMinit.txt", Z)

Z = np.transpose(np.nan_to_num(sigmaMcoll))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), shading="auto", cmap='magma')
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\sigma_M(z_{ini})$")#, size=50
)
plt.savefig("./Figure_11_sigmaMcoll.png")
np.savetxt("Figure_11_sigmaMcoll.txt", Z)

Z = np.transpose(np.nan_to_num(deltalonglow))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), shading="auto", cmap='magma')
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\delta_{long, low}$")#, size=50
)
plt.savefig("./Figure_11_deltalonglow.png")
np.savetxt("Figure_11_deltalonglow.txt", Z)

Z = np.transpose(np.nan_to_num(deltalonghigh))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), shading="auto", cmap='magma')
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\delta_{long, high}$")#, size=50
)
plt.savefig("./Figure_11_deltalonghigh.png")
np.savetxt("Figure_11_deltalonghigh.txt", Z)

Z = np.transpose(np.nan_to_num(deltalonghigh_ref))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), shading="auto", cmap='magma')
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\delta_{long, high}$")#, size=50
)
plt.savefig("./Figure_11_deltalonghigh_ref.png")
np.savetxt("Figure_11_deltalonghigh_ref.txt", Z)

Z = np.transpose(np.nan_to_num(deltacritsim))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), shading="auto", cmap='magma')
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"$\delta_{crit, sim}$")#, size=50
)
plt.savefig("./Figure_11_deltacritsim.png")
np.savetxt("Figure_11_deltacritsim.txt", Z)

Z = np.transpose(np.nan_to_num(hmf))
X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
fig, ax = plt.subplots(1,1)
hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), shading="auto", cmap='magma')
ax.tick_params(axis='both')
ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
ax.set_ylim((0.,0.1))
cbar = plt.colorbar(hmap)
cbar.set_label(label=(
    r"HMF")#, size=50
)
plt.savefig("./Figure_11_hmf.png")
np.savetxt("Figure_11_hmf.txt", Z)

