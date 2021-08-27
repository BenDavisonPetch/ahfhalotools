import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.tri as tri
import ahfhalotools.filetools as ft
from ahfhalotools.objects import Cluster

## -----------Load Cluster Instances--------------##

#define base file name (these are for full files)
fileNameBaseGX = "GadgetX-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseGiz = "GIZMO-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseMus = "GadgetMUSIC-NewMDCLUSTER_0001.z{z:.3f}"
#define directory
gxdir = "gadgetX\\"
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

##-----------Plot-------------##

## 2D plots of profile data
fig,ax = plt.subplots(2,3,figsize=(12,7))
quantities = ["starLocDens","gasLocDens"]
qnames = ["Stellar","Gas"]
haloID = 128000000000001

#normalize = mcolors.Normalize(vmin=4.5,vmax=10.3)
colmap = cm.plasma

for colIndex in range(3):
    cluster = clusters[colIndex]
    for rowIndex in range(2):
        ax1 = ax[rowIndex,colIndex]
        quantity = quantities[rowIndex]
        # x and y coordinates of points are r and z
        r,z,vals = cluster.funcOfRandZProfileDataPoints(haloID,quantity,removeZeroes = True)
        vals = vals * (r ** 2)

        #remove negative radii (test)
        #indexes = np.argwhere(r<0)
        #r = np.delete(r, indexes)
        #z = np.delete(z, indexes)
        #vals = np.delete(vals,indexes)

        r = np.log10(abs(r))
        vals = np.log10(abs(vals))

        #create triangulation
        triang = tri.Triangulation(r,z)

        #tripcolor plot
        tpc = ax1.tripcolor(triang,vals,cmap=colmap)
        cbar = fig.colorbar(tpc,ax=ax1)

        #add colorbar to last column
        if colIndex == 2:
            cbar.set_label("$\\log_{10}(\\; \\rho_{loc} \\; r^2 \\;)$")

        #add merge line to plot
        mergeZ,mergeSize = cluster.getLargestMergeZInRange(haloID,0,1)
        #print(mergeSize)
        mergeColor = 'k'
        ax1.add_line(Line2D(ax1.get_xlim(),[mergeZ,mergeZ],color=mergeColor))
        mergepad=5
        labelx = ax1.get_xlim()[0] + 0.1
        ax1.annotate('Largest merger event',xy=(labelx,mergeZ),xytext=(0,mergepad),xycoords="data",
                       textcoords="offset points", ha="left", va="bottom", c = mergeColor,
                       size='large')
        ax1.annotate('Size = {0}'.format(mergeSize), xy=(labelx,mergeZ),
                       xytext=(0,-mergepad), xycoords="data",
                       textcoords="offset points", ha="left", va="top", c = mergeColor,
                       size='large')


#label axes
ax[0,0].set_title("GadgetX")
ax[0,1].set_title("Gizmo")
ax[0,2].set_title("GadgetMUSIC")
[ax[1,i].set_xlabel("$\\log_{10}(r/kpc \\; h^{-1})$") for i in range(3)]
[ax[i,0].set_ylabel("Redshift, $z$") for i in range(2)]
fig.suptitle("Local Density x Radius Squared as a function of z and r")

#label rows
pad = 5
for i in range(2):
    ax[i,0].annotate(qnames[i],xy=(0,0.5),xytext=(-ax[i,0].yaxis.labelpad-pad,0),
            xycoords=ax[i,0].yaxis.label, textcoords='offset points',
            ha = 'right', va='center',size='large',rotation=90)

plt.tight_layout()
plt.show()
