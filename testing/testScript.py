import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from ahfhalotools.objects import Cluster
import ahfhalotools.filetools as ft

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

haloID = 128000000000001
zs,values = gxCluster.funcOfZHaloData(haloID, "Mvir")
print("#HALO DATA#\n\nzs:\n{0}\nMvir:\n{1}\n\n".format(zs,values))
dzs,deltaValues = gxCluster.funcOfZDeltaHaloData(haloID, "Mvir")
print("#DELTA HALO DATA#\n\nzs:\n{0}\ndeltaMvir:\n{1}\n\n".format(dzs,deltaValues))

print("\n\n## V TIME ##\n\n")

ages,avalues = gxCluster.funcOfAgeHaloData(haloID, "Mvir")
print("#HALO DATA#\n\nages:\n{0}\nMvir:\n{1}\n\n".format(ages,avalues))
dages,adeltaValues = gxCluster.funcOfAgeDeltaHaloData(haloID, "Mvir")
print("#DELTA HALO DATA#\n\nages:\n{0}\ndeltaMvir:\n{1}\n\n".format(dages,adeltaValues))
