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
deltacritlow = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
deltainitlow = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
deltacritsteplow = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
deltainitsteplow = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
deltacrithigh = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
deltainithigh = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
deltacritstephigh = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
deltainitstephigh = np.zeros((len(krefs), len(m_ax), len(omega_ax))) # Units: none
sigmaMinit = np.zeros((len(m_ax), len(omega_ax))) # Units: none
sigmaMcoll = np.zeros((len(m_ax), len(omega_ax))) # Units: none
dsigmaMinit =np.zeros((len(m_ax), len(omega_ax))) # Units: none

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

        deltacritdatalow = np.loadtxt("./FIG11_"+str(idx+1)+"/delta_crit_z0.65_M13.00_Nk50.dat", skiprows=1)[0:50,[1,2]]
        deltainitdatalow = np.loadtxt("./FIG11_"+str(idx+1)+"/delta_initial_z0.65_M13.00_Nk50.dat", skiprows=1)[0:50,[1,2]]
        deltacritinterplow = scipy.interpolate.interp1d(deltacritdatalow[:,0], deltacritdatalow[:,1])
        deltainitinterplow = scipy.interpolate.interp1d(deltainitdatalow[:,0], deltainitdatalow[:,1])
        deltacritlow[kidx, m_idx, o_idx] = deltacritinterplow(kval)
        deltainitlow[kidx, m_idx, o_idx] = deltainitinterplow(kval) 
        deltacritsteplow[kidx, m_idx, o_idx] = deltacritinterplow(kval)/deltacritinterplow(1.0e-4)
        deltainitsteplow[kidx, m_idx, o_idx] = deltainitinterplow(kval)/deltacritinterplow(1.0e-4)

        deltacritdatahigh = np.loadtxt("./FIG11_"+str(idx+1)+"/delta_crit_z0.65_M13.00_Nk50.dat", skiprows=1)[50:100,[1,2]]
        deltainitdatahigh = np.loadtxt("./FIG11_"+str(idx+1)+"/delta_initial_z0.65_M13.00_Nk50.dat", skiprows=1)[50:100,[1,2]]
        deltacritinterphigh = scipy.interpolate.interp1d(deltacritdatahigh[:,0], deltacritdatahigh[:,1])
        deltainitinterphigh = scipy.interpolate.interp1d(deltainitdatahigh[:,0], deltainitdatahigh[:,1])
        deltacrithigh[kidx, m_idx, o_idx] = deltacritinterphigh(kval)
        deltainithigh[kidx, m_idx, o_idx] = deltainitinterphigh(kval) 
        deltacritstephigh[kidx, m_idx, o_idx] = deltacritinterphigh(kval)/deltacritinterphigh(1.0e-4)
        deltainitstephigh[kidx, m_idx, o_idx] = deltainitinterphigh(kval)/deltacritinterphigh(1.0e-4)

        if kidx==0: 
            sigmaM = np.loadtxt("./FIG11_"+str(idx+1)+"/sigmaM_z0.65_M13.00_Nk1.dat", skiprows=1); 
            sigmaMinit[m_idx, o_idx] = sigmaM[0]         
            sigmaMcoll[m_idx, o_idx] = sigmaM[1]
            dsigmaMinit[m_idx, o_idx] = sigmaM[4]

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
    plt.savefig("./Figure_11_b1lstep_logk"+f"{np.log10(kval):.3f}.png")
    if (kidx==(len(krefs)-1)):
        plt.savefig("Figure_11.png")
        np.savetxt("Figure_11_bLstep.txt", Z)

    Z = np.transpose(np.nan_to_num(deltacritlow[kidx]))
    X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
    #interp = scipy.interpolate.interp2d(X, Y, Z, kind='linear') 
    #xn = np.arange(np.min(np.log10(m_ax)), np.max(np.log10(m_ax)), 0.01)
    #yn = np.arange(np.min(omega_ax/omega_cdm_LCDM), np.max(omega_ax/omega_cdm_LCDM), .001)
    #Z = interp(xn,yn)
    #X, Y = np.meshgrid(xn, yn) 
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
        r"$\delta_{crit}(k)~/~\delta_{crit}(k_{\rm ref})$")#, size=50
    )
    plt.savefig("./Figure_11_deltacrit_logk"+f"{np.log10(kval):.3f}.png")
    if (kidx==(len(krefs)-1)):
        np.savetxt("Figure_11_deltacritlow.txt", Z)

    Z = np.transpose(np.nan_to_num(deltainitlow[kidx]))
    X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
    #interp = scipy.interpolate.interp2d(X, Y, Z, kind='linear') 
    #xn = np.arange(np.min(np.log10(m_ax)), np.max(np.log10(m_ax)), 0.01)
    #yn = np.arange(np.min(omega_ax/omega_cdm_LCDM), np.max(omega_ax/omega_cdm_LCDM), .001)
    #Z = interp(xn,yn)
    #X, Y = np.meshgrid(xn, yn) 
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
        r"$\delta_{crit}(k)~/~\delta_{crit}(k_{\rm ref})$")#, size=50
    )
    plt.savefig("./Figure_11_deltainitlow_logk"+f"{np.log10(kval):.3f}.png")
    if (kidx==(len(krefs)-1)):
        np.savetxt("Figure_11_deltainitlow.txt", Z)  

    Z = np.transpose(np.nan_to_num(deltacrithigh[kidx]))
    X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
    #interp = scipy.interpolate.interp2d(X, Y, Z, kind='linear') 
    #xn = np.arange(np.min(np.log10(m_ax)), np.max(np.log10(m_ax)), 0.01)
    #yn = np.arange(np.min(omega_ax/omega_cdm_LCDM), np.max(omega_ax/omega_cdm_LCDM), .001)
    #Z = interp(xn,yn)
    #X, Y = np.meshgrid(xn, yn) 
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
        r"$\delta_{crit}(k)~/~\delta_{crit}(k_{\rm ref})$")#, size=50
    )
    plt.savefig("./Figure_11_deltacrithigh_logk"+f"{np.log10(kval):.3f}.png")
    if (kidx==(len(krefs)-1)):
        np.savetxt("Figure_11_deltacrithigh.txt", Z)

    Z = np.transpose(np.nan_to_num(deltainithigh[kidx]))
    X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
    #interp = scipy.interpolate.interp2d(X, Y, Z, kind='linear') 
    #xn = np.arange(np.min(np.log10(m_ax)), np.max(np.log10(m_ax)), 0.01)
    #yn = np.arange(np.min(omega_ax/omega_cdm_LCDM), np.max(omega_ax/omega_cdm_LCDM), .001)
    #Z = interp(xn,yn)
    #X, Y = np.meshgrid(xn, yn) 
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
        r"$\delta_{crit}(k)~/~\delta_{crit}(k_{\rm ref})$")#, size=50
    )
    plt.savefig("./Figure_11_deltainithigh_logk"+f"{np.log10(kval):.3f}.png")
    if (kidx==(len(krefs)-1)):
        np.savetxt("Figure_11_deltainithigh.txt", Z)  


    #Z = np.transpose(np.nan_to_num(deltacritstep[kidx]))
    #X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
    ##interp = scipy.interpolate.interp2d(X, Y, Z, kind='linear') 
    ##xn = np.arange(np.min(np.log10(m_ax)), np.max(np.log10(m_ax)), 0.01)
    ##yn = np.arange(np.min(omega_ax/omega_cdm_LCDM), np.max(omega_ax/omega_cdm_LCDM), .001)
    ##Z = interp(xn,yn)
    ##X, Y = np.meshgrid(xn, yn) 
    #fig, ax = plt.subplots(1,1)
    #hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(deltacritstep), vmax=np.max(deltacritstep), shading="auto", cmap='magma')
    #if (np.max(Z)>1.01):
    #    CS = ax.contour(X, Y, Z, np.linspace(1.01, 1.05, 5), colors='white')
    #    ax.clabel(CS, CS.levels, inline=False, fmt=fmt)
    #ax.tick_params(axis='both')
    #ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
    #ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
    #ax.set_ylim((0.,0.1))
    #cbar = plt.colorbar(hmap)
    #cbar.set_label(label=(
    #    r"$\delta_{crit}(k)~/~\delta_{crit}(k_{\rm ref})$")#, size=50
    #)
    #plt.savefig("./Figure_11_deltacritstep_logk"+f"{np.log10(kval):.3f}.png")
    #if (kidx==(len(krefs)-1)):
    #    np.savetxt("Figure_11_deltacritstep.txt", Z)

    #Z = np.transpose(np.nan_to_num(deltainitstep[kidx]))
    #X, Y = np.meshgrid(np.log10(m_ax), omega_ax/omega_cdm_LCDM)
    ##interp = scipy.interpolate.interp2d(X, Y, Z, kind='linear') 
    ##xn = np.arange(np.min(np.log10(m_ax)), np.max(np.log10(m_ax)), 0.01)
    ##yn = np.arange(np.min(omega_ax/omega_cdm_LCDM), np.max(omega_ax/omega_cdm_LCDM), .001)
    ##Z = interp(xn,yn)
    ##X, Y = np.meshgrid(xn, yn) 
    #fig, ax = plt.subplots(1,1)
    #hmap = ax.pcolormesh(X, Y, Z, vmin=np.min(deltainitstep), vmax=np.max(deltainitstep), shading="auto", cmap='magma')
    #if (np.max(Z)>1.01):
    #    CS = ax.contour(X, Y, Z, np.linspace(1.01, 1.05, 5), colors='white')
    #    ax.clabel(CS, CS.levels, inline=False, fmt=fmt)
    #ax.tick_params(axis='both')
    #ax.set_xlabel(r"$\log{\left(m_\phi ~/~ {\rm [eV]}\right)}$")
    #ax.set_ylabel(r"$\omega_{\phi,0} ~/~ \omega_{{\rm d}, 0}$")
    #ax.set_ylim((0.,0.1))
    #cbar = plt.colorbar(hmap)
    #cbar.set_label(label=(
    #    r"$\delta_{crit}(k)~/~\delta_{crit}(k_{\rm ref})$")#, size=50
    #)
    #plt.savefig("./Figure_11_deltainitstep_logk"+f"{np.log10(kval):.3f}.png")
    #if (kidx==(len(krefs)-1)):
    #    np.savetxt("Figure_11_deltainitstep.txt", Z)  


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
plt.savefig("./sigmaMinit.png")
np.savetxt("sigmaMinit.txt", Z)

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
plt.savefig("./sigmaMcoll.png")
np.savetxt("sigmaMcoll.txt", Z)

Z = np.transpose(np.nan_to_num(dsigmaMinit))
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
plt.savefig("./dsigmaMinit.png")
np.savetxt("dsigmaMinit.txt", Z)

