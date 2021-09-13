import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from ahfhalotools.objects import Cluster
import ahfhalotools.filetools as ft


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

## 3D plots of profile data
fig, (ax1, ax2) = plt.subplots(ncols=2, subplot_kw={'projection':'3d'})
ax1.view_init(50,-73)
ax2.view_init(50,-73)
#first plot
rr,zz,vals = musCluster.funcOfRandZProfileDataPoints(128000000000001,"starLocDens",removeZeroes = True)
vals = vals * (rr**2)

#might be able to get custom coloring using:
#https://stackoverflow.com/questions/63298864/plot-trisurface-with-custom-color-array

#ax1.plot_trisurf(np.log10(abs(rr)),zz,np.log10(abs(vals)), cmap=cm.plasma,alpha = 1,label="Stellar")
ax1.plot_trisurf(np.log10(abs(rr)),zz,np.log10(abs(vals)), cmap=cm.plasma,alpha = 1,label="Stellar")
#second plot
rr,zz,vals = musCluster.funcOfRandZProfileDataPoints(128000000000001,"gasLocDens",removeZeroes = True)
vals = vals * (rr**2)
ax2.plot_trisurf(np.log10(abs(rr)),zz,np.log10(abs(vals)), cmap=cm.plasma,alpha = 1,label="Gas")

#setting labels
ax1.set_xlabel("$\\log_{10}(\\; r \\; / \\; kpc \\; h^{-1}\\;)$")
ax2.set_xlabel("$\\log_{10}(\\; r \\; / \\; kpc \\; h^{-1}\\;)$")
ax1.set_ylabel("Redshift, $z$")
ax2.set_ylabel("Redshift, $z$")
ax1.set_zlabel("$\\log_{10}(\\; \\rho_{loc} \\: r^2 \\; )$")
ax2.set_zlabel("$\\log_{10}(\\; \\rho_{loc} \\: r^2 \\; )$")
ax1.set_title("Local Stellar Density")
ax2.set_title("Local Gas Density")
#ax.legend()
plt.show()
