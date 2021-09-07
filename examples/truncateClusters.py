#!/usr/bin/env python3.8

"""
FILE TRUNCATION SCRIPT
This script was used on castor to truncate all of the clusters at once before
downloading. It takes a **very** long time to run.
"""

import numpy as np
import ahfhalotools.filetools as ft
from ahfhalotools.objects import Cluster

#define base file name for full, untruncated files
fileNameBaseGX = "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseGiz = fileNameBaseGX
fileNameBaseMus = "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.z{z:.3f}"
#define directory for untruncated files
directory = "/home/nifty2014/TheThreeHundred/catalogues/AHF/{simName}/"
#define snapshots to load
snapNos = np.arange(97,129)

#get snap num to redshift map
rsmapGX = ft.getSnapNumToZMapGX("/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_0001")
rsmapGiz = ft.getSnapNumToZMapGiz("/home/nifty2014/TheThreeHundred/catalogues/AHF/GIZMO/NewMDCLUSTER_0001")

#get redshifts from snapNos
zsGX = np.array([rsmapGX[num] for num in snapNos])
zsGiz = np.array([rsmapGiz[num] for num in snapNos])
#have to get redshifts of music snaps directly from directory. We only want the
# last len(snapNos) entries
zsMus = ft.getMusZs("/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetMUSIC/NewMDCLUSTER_0001")[-len(snapNos):]
assert(len(zsGX)==len(zsGiz)==len(zsMus))

#define where truncated files will go
outputDir = "/home/nifty2014/TheThreeHundred/playground/BDP/TruncData/{simName}"

## Truncation ##

#define how many halos truncated files should contain information on
truncSize = 1

#cluster numbers
clusterNums = np.arange(1,325)

print("====================================")
print("====== GadgetMUSIC Truncation ======")
print("====================================\n")

ft.truncateClusters(clusterNums, snapNos, zsMus, "GadgetMUSIC", truncSize, outputDir,
                     directory = directory, skipmtree = True,
                     fileBaseFmt = fileNameBaseMus)

print("\n\n")
print("==============================")
print("====== GIZMO Truncation ======")
print("==============================\n")

ft.truncateClusters(clusterNums, snapNos, zsGiz, "GIZMO", truncSize, outputDir,
                     directory = directory, skipmtree = True)

print("\n\n")
print("================================")
print("====== GadgetX Truncation ======")
print("================================\n")
ft.truncateClusters(clusterNums, snapNos, zsGX, "GadgetX", truncSize, outputDir,
                     directory = directory, skipmtree = True)


#load gadget X clusters once we're done to check it works
clusters = ft.loadClusters(clusterNums, snapNos, zsGX, "GadgetX", directory = outputDir)
print(clusters)
