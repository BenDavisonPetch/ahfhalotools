import numpy as np
import ahfhalotools.filetools as ft
from ahfhalotools.objects import Cluster

#define base file name for full, untruncated files
fileNameBaseGX = "GadgetX-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseGiz = "GIZMO-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseMus = "GadgetMUSIC-NewMDCLUSTER_0001.z{z:.3f}"
#define directory for untruncated files
gxdir = "gadgetX\\"
gizdir = "gizmo\\"
musdir = "music\\"
#define snapshots to load
snapNosGX = np.arange(97,129)
snapNosGiz = np.arange(97,129)
#n.b. music doesnt have snapshot numbers (bruh)

#get snap num to redshift map (ft.getMusZs will return all redshifts present in
# music's directory)
zsMus = ft.getMusZs(musdir)
rsmapGiz = ft.getSnapNumToZMapGiz(gizdir)
rsmapGX = ft.getSnapNumToZMapGX(gxdir)

snapNosMus = np.arange(129-len(zsMus),129)
#get redshifts from snapNos
zsGX = np.array([rsmapGX[num] for num in snapNosGX])
zsGiz = np.array([rsmapGiz[num] for num in snapNosGiz])

#define where truncated files will go
truncFileGiz = "data\\cluster1trunc\\GIZMO\\GIZMO_C1_snap{snap:0=3d}.z{z:.3f}"
truncFileMus = "data\\cluster1trunc\\GadgetMUSIC\\GadgetMUSIC_C1_snap.z{z:.3f}"
truncFileGX = "data\\cluster1trunc\\GadgetX\\GadgetXX_C1_snap{snap:0=3d}.z{z:.3f}"

## Truncation ##

#define how many halos truncated files should contain information on
truncSize = 5

ft.truncateFiles(gxdir + fileNameBaseGX,snapNosGX,zsGX,truncFileGX,truncSize)
ft.truncateFiles(gizdir + fileNameBaseGiz,snapNosGiz,zsGiz,truncFileGiz,truncSize)
ft.truncateFiles(musdir + fileNameBaseMus,snapNosMus,zsMus,truncFileMus,truncSize)

#load truncated files once we're done (note, you only need to run truncateFiles
# ONCE for each set of data)
gxCluster = Cluster(truncFileGX, snapNosGX, zsGX)
gizCluster = Cluster(truncFileGiz, snapNosGiz, zsGiz)
musCluster = Cluster(truncFileMus, snapNosMus, zsMus)
