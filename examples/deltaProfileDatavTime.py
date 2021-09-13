import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from ahfhalotools.objects import Cluster
import ahfhalotools.filetools as ft
import ahfhalotools.analysis as analysis

## -----------Load Cluster Instances--------------##

#define base file name (these are for full files)
fileNameBaseGX = "GadgetX-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseGiz = "GIZMO-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseMus = "GadgetMUSIC-NewMDCLUSTER_0001.z{z:.3f}"
#define directory
gxdir = "GadgetX\\NewMDCLUSTER_0001\\"
gizdir = "gizmo\\"
musdir = "music\\"
#define snapshots to load
snapNosGX = np.arange(97,129)
snapNosGiz = np.arange(97,129)
#n.b. music doesnt have snapshot numbers (bruh)

#get snap num to redshift map
zsMus = ft.getMusZs(musdir)
rsmapGiz = ft.getSnapNumToZMapGiz(gizdir)
rsmapGX = ft.getSnapNumToZMapGX(gxdir)

snapNosMus = np.arange(129-len(zsMus),129)
#get redshifts from snapNos
zsGX = np.array([rsmapGX[num] for num in snapNosGX])
zsGiz = np.array([rsmapGiz[num] for num in snapNosGiz])

#truncated file names
truncFileGX = "first10\\GadgetX\\GX_first10_snap{snap:0=3d}.z{z:.3f}"
truncFileGiz = "first10\\Gizmo\\giz_first10_snap{snap:0=3d}.z{z:.3f}"
truncFileMus = "first10\\Music\\mus_first10_snap.z{z:.3f}"

#CLUSTER OBJECTS
gxCluster = Cluster(truncFileGX, snapNosGX, zsGX)
gizCluster = Cluster(truncFileGiz, snapNosGiz, zsGiz)
musCluster = Cluster(truncFileMus, snapNosMus, zsMus)
clusters = [gxCluster,gizCluster,musCluster]
clusterNames = ["GadgetX","Gizmo","GadgetMUSIC"]
clusterColors = ["C0","C1","C2"]

##----------PLOT-----------##

## comparing sims on deltaM/M over time for total mass, gas mass and Stellar
# mass

fig, axes = plt.subplots(1,3,figsize=(15,6),sharey='row')
haloID = 128000000000001
radius = 100

axisTitles = ["Local Total Mass Density", "Local Gas Density", "Local Stellar Density"]
quantities = ["locDens", "gasLocDens", "starLocDens"]

#plot on each axis
for i in range(3):
    ax = axes[i]

    #loop through each simulation
    for j in range(3):
        color = clusterColors[j]
        cluster = clusters[j]

        ## can get delta values in two ways:
        #age, deltaM = cluster.funcOfAgeDeltaProfileData(haloID,quantities[i],radius)
        age, deltaM = cluster.funcOfAgeProfileData(haloID,"delta"+quantities[i],radius)

        _, M = cluster.funcOfAgeProfileData(haloID,quantities[i],radius)
        #M is going to be 1 element longer than deltaM, as we cannot get delta
        # for a quantity with unknown previous value
        # note that delta M should be divided by M of halo of *earlier* snapshot
        M = M[:-1]
        deltaMoverM = deltaM/M
        #deltaMoverM = deltaM

        #plot
        ax.plot(age,deltaMoverM, label=clusterNames[j],c=color)

        #add in line to represent time of largest merge
        mergeZ, mergeSize = cluster.getLargestMergeZInRange(haloID,0,1,scheme="halodata",fractional=False)
        mergeTime = analysis.tfromz(mergeZ)
        ls = ["--","-.",":"][j]
        ax.axvline(mergeTime,c=color,alpha=0.7,lw=(4-j),ls=ls)

    #set labels
    ax.set_title(axisTitles[i])
    ax.set_xlabel("age, $t$ / Gyr")
    ax.set_ylabel("$\\frac{\\Delta \\rho}{\\rho} \\mid r="+str(radius)+"\\; kpc \\: h^{-1}$")
    ax.legend()

    #annotate merge lines
    pad=15
    topCoord = ax.get_ylim()[1]
    toppad = 0.02*(topCoord-ax.get_ylim()[0])
    if i == 2: ax.annotate('Largest merger events\nusing halo data',xy=(mergeTime,topCoord-toppad),xytext=(-pad,0),xycoords="data",
                   textcoords="offset points", ha="right", va="top", c = 'dimgrey')


#set figure title
fig.suptitle("$\\frac{\\Delta \\rho}{\\rho}$ as a function of age at $r="+str(radius)+"\\; kpc \\: h^{-1}$")

plt.tight_layout()
plt.show()
