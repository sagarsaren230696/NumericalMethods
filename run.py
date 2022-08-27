from cProfile import label
import numpy as np
import pandas as pd
import subprocess as sp

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
import matplotlib as mpl

import os
import CoolProp.CoolProp as CP
from preos import Molecule,preos


kb = 1.3806448e-23 # Boltzmann constant J/K
N = 6.023e23 # Avogadro's number

data = {"Energy":[],"Pressure(bar)":[],"Loading(vSTP/v)":[]}
df = pd.DataFrame(data)

Q0 = np.linspace(4,20,9) # Heat of adsorption in kJ/mol
U0 = -Q0/kb/N*1000

plt.style.use('classic')
parameters = {'axes.labelsize': 24,
            'axes.titlesize': 18,
            "legend.fontsize" : 16,
            'figure.autolayout': True,
            'axes.edgecolor':'black',
            'axes.linewidth':2,
            'axes.titlecolor': "black",
            'legend.fontsize': 16,
            'lines.linewidth':3,
            'mathtext.fontset':'stix',
            'font.family':'STIXGeneral'}
plt.rcParams.update(parameters)
colors = plt.cm.autumn(np.linspace(0,1,len(Q0)))

fig,ax = plt.subplots()
ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
ax.tick_params(direction='in', which="major", length=6, width=2, labelsize=18)
ax.tick_params(direction='in', which="minor", length=4, width=2, labelsize=18)
# pressures = np.concatenate((np.linspace(0.1,0.9,9)*1e6,np.linspace(1,10,10)*1e6),axis=None)
pressures = np.linspace(5000,80000,16)
print(pressures)
uptakes = []
# methane = Molecule("methane", -82.59 + 273.15, 45.99, 0.011)
argon = Molecule("argon", 150.69 + 273.15, 48.63, -0.00219)
for i,p in enumerate(pressures):
    print("Simulating empty_GCMC at pressure {}".format(p))
    fugacity = preos(argon,87,p*1e-5,plotcubic=False,printresults=False)['fugacity(bar)']*1e5
    print(fugacity)
    if os.path.exists(r"./GrapheneArgonResults/SimulationResults/SimulationResults_{:.2f}.txt".format(p)) == False:
        sp.call([r"./empty_GCMC.exe","20000",f"{fugacity}",f"{p}"])
    dataFile = pd.read_csv(r"./GrapheneArgonResults/SimulationResults/SimulationResults_{:.2f}.txt".format(p),sep=',',header=2)
    uptakes.append(dataFile["Loading(vSTP/v)"].values[0])
ax.plot(pressures/1e6,uptakes,label="Steele 10-4-3 potential")
LJresults = pd.read_csv(r"./GrapheneArgonResults/OutputDataNew_CH4_graphite-sheet_3-layers_40A.csv")
ax.plot(LJresults['Pressure'].values/1e6,LJresults['uptakeVol'].values,label="LJ potential")
ax.set_xlabel("Pressure (MPa)")
ax.set_ylabel("argon uptake (cc(STP)/cc)")
ax.legend(loc="best")
plt.savefig(r"./GrapheneArgonResults/uptake_profile.png",dpi=300)


fig,ax = plt.subplots()
ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
ax.tick_params(direction='in', which="major", length=6, width=2, labelsize=18)
ax.tick_params(direction='in', which="minor", length=4, width=2, labelsize=18)
colors_pressure = plt.cm.rainbow(np.linspace(0,1,len(pressures)))
for i,p in enumerate(pressures):
    dataFile = pd.read_csv(f"./GrapheneArgonResults/NumberDensity/NumberDensity_{p:.2f}.csv")
    dataFile.columns = ["zPositions","avgNumOfMethanesAtZ"]
    ax.plot(dataFile["zPositions"].values,dataFile["avgNumOfMethanesAtZ"].values*160/6.023,color=colors_pressure[i],label=f"{p*1e-6} MPa",linewidth=0.5)
    # ax.hlines(CP.PropsSI('D','P',p,'T',298.15,'methane')*1e-3,0,40,colors=colors_pressure[i],linestyles='dashed',linewidth=1)
ax.set_xlabel("Z distance (A)")
ax.set_ylabel("Methane Number Density (molec./$A^3$)")
# ax.legend(loc="best")
cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=round(pressures[0]*1e-6,2),vmax=round(pressures[-1]*1e-6,2)), 
                    cmap=mpl.cm.rainbow), orientation='vertical', label='Pressure (MPa)')
cbar.ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
cbar.ax.tick_params(direction='in', which="major", length=6, width=2, labelsize=18)
cbar.ax.tick_params(direction='in', which="minor", length=4, width=2, labelsize=18)
cbar.ax.set_yscale('log')
plt.savefig(r"./GrapheneArgonResults/NumberDensity.png",dpi=300)

# fig,ax = plt.subplots()
# ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
# ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
# ax.tick_params(direction='in', which="major", length=6, width=2, labelsize=18)
# ax.tick_params(direction='in', which="minor", length=4, width=2, labelsize=18)
# for i,u in enumerate(U0):
#     print("Simulating empty_GCMC with background energy {}".format(u))
#     # sp.call([r"./empty_GCMC.exe","{}".format(u),"{}".format(20000)])
#     dataFile = pd.read_csv(r"./results_empty_box/SIM_U{:.2f}.txt".format(u),sep=',',header=2)
#     ax.plot(dataFile["Pressure(bar)"].values/dataFile["Pressure(bar)"].values[-1],dataFile["Loading(vSTP/v)"].values,color=colors[i],label="{} kJ/mol".format(Q0[i]))
# ax.set_xlabel("Relative pressure (-)")
# ax.set_ylabel("Methane uptake (cc(STP)/cc)")
# ax.legend(loc="best")
# plt.savefig(r"./results_empty_box/uptake_profile.png",dpi=300)


# fig,ax = plt.subplots()
# ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
# ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
# ax.tick_params(direction='in', which="major", length=6, width=2, labelsize=18)
# ax.tick_params(direction='in', which="minor", length=4, width=2, labelsize=18)
# xArr = []
# yArr = []
# for i,u in enumerate(U0):
#     dataFile = pd.read_csv(r"./results_empty_box/SIM_U{:.2f}.txt".format(u),sep=',',header=2)
#     xArr.append(Q0[i])
#     yArr.append(dataFile["Loading(vSTP/v)"].values[3]-dataFile["Loading(vSTP/v)"].values[1])
# ax.plot(xArr,yArr)
# ax.set_xlabel("$\Delta H_{ads}$ (kJ/mol)")
# ax.set_ylabel("Deliverable capacity (cc(STP)/cc)")
# plt.savefig(r"./results_empty_box/deliverable capacity.png",dpi=300)


# listOfSiteToSiteDistance = np.linspace(4,10,7)
# listOfOptimumDC = []
# listOfOptimumDC2 = []
# listOfOptimumEnergies = []
# listOfOptimumEnergies2 = []
# fig,ax = plt.subplots()
# ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
# ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
# ax.tick_params(direction='in', which="major", length=6, width=2, labelsize=18)
# ax.tick_params(direction='in', which="minor", length=4, width=2, labelsize=18)
# xArr = []
# yArr = []
# for d in listOfSiteToSiteDistance:
#     # sp.call(["python",r"./scale_and_replicate_fcc.py","{}".format(d)])
#     maxV = -1e10
#     maxV2 = -1e10
#     for i,u in enumerate(U0):
#         print("Simulating lattice_GCMC with background energy {}".format(u))
#         # sp.call([r"./lattice_GCMC.exe","{}".format(u),"200","{}".format(d),"0"])
#         dataFile = pd.read_csv(r"./results_lattice_model/SIM_UnoGG{:.2f}_d{:.2f}.txt".format(u,d),sep=',',header=3)
#         # sp.call([r"./lattice_GCMC.exe","{}".format(u),"200","{}".format(d),"1"])
#         dataFile2 = pd.read_csv(r"./results_lattice_model/SIM_U{:.2f}_d{:.2f}.txt".format(u,d),sep=',',header=3)
#         DC = dataFile["Loading(vSTP/v)"].values[3]-dataFile["Loading(vSTP/v)"].values[1]
#         DC2 = dataFile2["Loading(vSTP/v)"].values[3]-dataFile2["Loading(vSTP/v)"].values[1]
#         if maxV <= DC:
#             maxV = DC
#             iMax = i
#         if maxV2 <= DC2:
#             maxV2 = DC2
#             iMax2 = i
#     listOfOptimumDC.append(maxV)
#     listOfOptimumDC2.append(maxV2)
#     listOfOptimumEnergies.append(Q0[iMax])
#     listOfOptimumEnergies2.append(Q0[iMax2])
#     minC,maxC = min(listOfOptimumEnergies+listOfOptimumEnergies2),max(listOfOptimumEnergies+listOfOptimumEnergies2)
# # print(listOfOptimumEnergies)
# # print(listOfOptimumEnergies2)
# ax.plot(listOfSiteToSiteDistance,listOfOptimumDC,label="w/o CH4-CH4 interaction")
# ax.plot(listOfSiteToSiteDistance,listOfOptimumDC2,'--',label="with CH4-CH4 interaction")
# ax.set_ylim(bottom=0.0)
# plt.scatter(listOfSiteToSiteDistance,listOfOptimumDC,c=listOfOptimumEnergies,cmap="autumn",s=200)
# plt.clim(minC,maxC)
# plt.scatter(listOfSiteToSiteDistance,listOfOptimumDC2,c=listOfOptimumEnergies2,cmap="autumn",s=200)
# plt.clim(minC,maxC)
# plt.colorbar()
# ax.set_xlabel("d (A)")
# ax.set_ylabel("Maximum DC (cc(STP)/cc)")
# ax.legend(loc='best')
# plt.savefig(r"./results_lattice_model/max DC.png",dpi=300)

