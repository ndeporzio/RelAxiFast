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
krefs = np.logspace(-3.5, -0.5, 4) # Units: 

b1e = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
b1l = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
b1e_step = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
b1l_step = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none

#b1e = [0]*len(krefs)
#b1l = [0]*len(krefs)
#b1e_step = [0]*len(krefs)
#b1l_step = [0]*len(krefs)

def fmt(x):
    s = f"{(x-1.)*100.:.0f}"
    return rf"${s} \%$" if plt.rcParams["text.usetex"] else f"{s} %"

for kidx, kval in enumerate(krefs):
    for idx in range(nruns):
        m_idx = np.mod(idx, len(m_ax)) 
        o_idx = np.floor_divide(idx, len(m_ax))

        b1e[kidx, m_idx, o_idx] = np.float(np.loadtxt("./FIG11_"+str(idx+1)+"/Figure_11_b1e_logk"+f"{np.log10(kval):.3f}.txt")) 
        b1l[kidx, m_idx, o_idx] = np.float(np.loadtxt("./FIG11_"+str(idx+1)+"/Figure_11_b1l_logk"+f"{np.log10(kval):.3f}.txt")) 
        b1e_step[kidx, m_idx, o_idx] = np.float(np.loadtxt("./FIG11_"+str(idx+1)+"/Figure_11_b1estep_logk"+f"{np.log10(kval):.3f}.txt")) 
        b1l_step[kidx, m_idx, o_idx] = np.float(np.loadtxt("./FIG11_"+str(idx+1)+"/Figure_11_b1lstep_logk"+f"{np.log10(kval):.3f}.txt")) 

    Z = np.transpose(np.nan_to_num(b1l[kidx]))
    X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
    #interp = scipy.interpolate.interp2d(X, Y, Z, kind='linear') 
    #xn = np.arange(np.min(np.log10(m_ax)), np.max(np.log10(m_ax)), 0.01)
    #yn = np.arange(np.min(omega_ax/omega_cdm_LCDM), np.max(omega_ax/omega_cdm_LCDM), .001)
    #Z = interp(xn,yn)
    #X, Y = np.meshgrid(xn, yn) 
    fig, ax = plt.subplots(1,1)
    hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), shading="auto")
    ax.tick_params(axis='both')
    ax.set_xlabel(r"$\log{m_\phi / {\rm [eV]}}$")
    ax.set_ylabel("$\omega_\phi / \omega_d$")
    ax.set_ylim((0.,0.1))
    cbar = plt.colorbar(hmap)
    plt.savefig("./Figure_11_b1l_logk"+f"{np.log10(kval):.3f}.png")

    Z = np.transpose(np.nan_to_num(b1l_step[kidx]))
    X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
    #interp = scipy.interpolate.interp2d(X, Y, Z, kind='linear') 
    #xn = np.arange(np.min(np.log10(m_ax)), np.max(np.log10(m_ax)), 0.01)
    #yn = np.arange(np.min(omega_ax/omega_cdm_LCDM), np.max(omega_ax/omega_cdm_LCDM), .001)
    #Z = interp(xn,yn)
    #X, Y = np.meshgrid(xn, yn) 
    fig, ax = plt.subplots(1,1)
    hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(b1l_step), vmax=np.max(b1l_step), shading="auto", cmap='magma')
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
    plt.savefig("./Figure_11_b1lstep_logk"+f"{np.log10(kval):.3f}.png")
    if (kidx==(len(krefs)-1)):
        plt.savefig("Figure_11.png")
        np.savetxt("Figure_11_zdata.txt", Z) 
