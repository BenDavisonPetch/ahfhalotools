import numpy as np
import ahfhalotools.filetools as ft
from ahfhalotools.objects import Cluster

#define base file name for full, untruncated files
fileNameBaseGX = "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.snap_{snap:0=3d}.z{z:.3f}"
#define directory for untruncated files
directory = "{simName}/"
#define snapshots to load
snapNosGX = np.arange(97,129)

#get snap num to redshift map (ft.getMusZs will return all redshifts present in
# music's directory)
rsmapGX = ft.getSnapNumToZMapGX("GadgetX/NewMDCLUSTER_0001")

snapNosMus = np.arange(129-len(zsMus),129)
#get redshifts from snapNos
zsGX = np.array([rsmapGX[num] for num in snapNosGX])

#define where truncated files will go
outputDir = "GadgetXTrunc/"

## Truncation ##

#define how many halos truncated files should contain information on
truncSize = 1

#cluster numbers
clusterNums = [1,2]
simName = "GadgetX"

ft.truncateClusters(clusterNums, snapNums, zs, simName, truncSize, outputDir,
                     directory = directory)

#load clusters once we're done
clusters = ft.loadClusters(clusterNums, snapNums, zs, simName, directory = outputDir)
print(clusters)
