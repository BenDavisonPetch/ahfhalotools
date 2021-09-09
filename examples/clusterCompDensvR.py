import numpy as np
import ahfhalotools.filetools as ft
from ahfhalotools.objects import Cluster

#define base file name for full, untruncated files
fileNameBaseGX = "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseGiz = fileNameBaseGX
fileNameBaseMus = "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.z{z:.3f}"
#define directory to read from
directory = "TruncData/{simName}/"

#define snapshots to use
snapNos = np.arange(97,129)

## Truncation ##

#define how many halos truncated files should contain information on
truncSize = 1

#cluster numbers
clusterNums = np.arange(1,325)

# N.B. for GadgetMUSIC data, for some reason, for cluster numbers 22 and above,
# there are only 17 snapshots. The snapshots with redshift less than one are
# 7-17, so we truncate these separately.

# #up to 128
# musClusters1 = ft.loadClusters(np.arange(1,22), snapNos, "GadgetMUSIC",
#                                directory = directory,
#                                fileBaseFmt = fileNameBaseMus,
#                                printProgress = True, skipmtree = True)
# #up to 17
# musClusters2 = ft.loadClusters(np.arange(22,325), np.arange(7,18), "GadgetMUSIC",
#                                directory = directory,
#                                fileBaseFmt = fileNameBaseMus,
#                                printProgress = True, skipmtree = True)
#
# gizClusters = ft.loadClusters(clusterNums, snapNos, "GIZMO",
#                               directory = directory, printProgress = True,
#                               skipmtree = True)

gxClusters = ft.loadClusters(clusterNums, snapNos, "GadgetX",
                             directory = directory, printProgress = True,
                             skipmtree = True)
print(gxClusters)
