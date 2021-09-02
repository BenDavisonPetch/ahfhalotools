# Merger Detection Scheme Guide

In AHFHaloTools, each of the functions that detect mergers allow for two optional arguments that can be used to specify how mergers are quantified. These arguments are `fractional : bool`, and `scheme : str`. This document serves as an explanation of the differences between schemes and the usage of these arguments.

## Fractional vs Absolute Merger Sizes
`fractional` is a boolean argument that specifies whether calculated merger sizes should be fractional or in absolute size. When fractional merger sizes are calculated, they are relative to the size of the __father__ halo.

## Schemes
There are two distinct type of schemes:
* The Merger Tree (mtree) schemes, and
* The halo data scheme

### Merger Tree (mtree) Schemes
The mtree schemes use `.AHF_mtree` and `.BDP_enchalos` files to determine the size of mergers using the number of shared particles from progenitors. In each mtree scheme, the merger tree entry of the halo is scanned, and any progenitors that were subhalos of the father halo are discarded, as is the entry for the father halo itself. This leaves a list of "new" particle contributions to the halo (ie particles that weren't already in the halo in the last snapshot)[^1].

__mtree-sum__ is the recommended scheme out of the mtree schemes. It sums each of the particle contributions to determine the merger size.

__mtree-largest__ uses _only one_ progenitor's particle contribution, specifically the largest.

__mtree-first__ also uses _only one_ progenitor's particle contribution, specifically the progenitor (that wasn't already a sub-halo or the father) farthest up in the merger tree. The order of entries in the merger tree is slightly ambiguous as per [AHF's documentation](http://popia.ft.uam.es/AHF/files/AHF.pdf), therefore this scheme is not recommended unless under specific use-cases.

If `fractional` is False, these sizes will just be the number of particles added to the halo.
If `fractional` is True, these sizes will be expressed as fractions of the number of particles of the father halo.
Please note that `mtree` fractional schemes generally result in ___very low merge sizes___ (of the order of 0.001)

### The Halo Data Scheme (default)
The halodata scheme determines merger size as the change in total virial mass between snapshots, using data from the `.AHF_halos` files.
If `fractional` is False, the merger size will be in Msol/h
If `fractional` is True, the merger size will be given as deltaM / M

## Scheme Choice
Because AHF ignores halos below a certain mass threshold (both in the halo finder itself and when generating merger trees), the actual number of particles introduced into a halo between snapshots can be larger than the sum of contributions in the .AHF_mtree file. As such, the `halodata` scheme is best for detecting points at which large amounts of mass have been added to the halo as a result of mergers with many small structures, and the `mtree` schemes are best for detecting mergers with only large structures.

![Comparison of Merger Detection Schemes](https://raw.githubusercontent.com/BenDavisonPetch/ahfhalotools/main/gallery/mergerSchemeComp.png)

[^1]: If this wasn't done, halos that were already sub-halos of the father halo, and so already contributed to the number of particles in the father, would be counted again as "new particles"
