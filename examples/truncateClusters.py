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

#define snapshots to use
snapNos = np.arange(97,129)

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
# N.B. for GadgetMUSIC data, for some reason, for cluster numbers 22 and above,
# there are only 17 snapshots. The snapshots with redshift less than one are
# 7-17, so we truncate these separately.
ft.truncateClusters(np.arange(1,22), snapNos, "GadgetMUSIC", truncSize, outputDir,
                     directory = directory, skipmtree = True,
                     fileBaseFmt = fileNameBaseMus)
ft.truncateClusters(np.arange(22,325), np.arange(7,18), "GadgetMUSIC", truncSize, outputDir,
                     directory = directory, skipmtree = True,
                     fileBaseFmt = fileNameBaseMus)

print("\n\n")
print("==============================")
print("====== GIZMO Truncation ======")
print("==============================\n")

ft.truncateClusters(clusterNums, snapNos, "GIZMO", truncSize, outputDir,
                     directory = directory, skipmtree = True)

print("\n\n")
print("================================")
print("====== GadgetX Truncation ======")
print("================================\n")
# cluster 10 for GadgetX appears to be missing snapshot 100 so will be absent
# cluster 86 also has something weird in the data so will be skipped too
ft.truncateClusters(clusterNums, snapNos, "GadgetX", truncSize, outputDir,
                    directory = directory, skipmtree = True)


#load gadget X clusters once we're done to check it works
clusters = ft.loadClusters(clusterNums, snapNos, "GadgetX", directory = outputDir)
print(clusters)
