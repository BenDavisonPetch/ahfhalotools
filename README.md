# AHFHaloTools

AHFHaloTools is a Python 3 library for the analysis of data produced by AMIGA's Halo Finder (AHF).

## Features
* File parsing for `.AHF_halos`, `.AHF_profiles`, `.AHF_mtree(_idx)`
* File truncation
* Analysis of radial profile data
* Analysis of integral properties (halo data)
* Proper tracking of halos through time
* Merger detection
* Analysis of enclosed halos as indicator of host halo kinematics

## Dependencies
* [NumPy](https://numpy.org/)
* [matplotlib](https://matplotlib.org/)
* [AstroPy](https://www.astropy.org/)

## Installation
`pip install ahfhalotools`

## Usage and Documentation
The majority of analysis is enabled by the `ahfhalotools.objects.Cluster` object.
Example analysis and truncation scripts are available in `/examples`.

Documentation is available in two places: the first is as docstrings within the
code, which can be viewed using the built-in ``help`` function:

```python
    >>> from ahfhalotools.objects import Cluster
    >>> help(Cluster)
    Help on class Cluster in module ahfhalotools.objects: ...
```

The other way to view documentation is in a web browser [here](https://htmlpreview.github.io/?https://raw.githubusercontent.com/BenDavisonPetch/ahfhalotools/main/docs/ahfhalotools/index.html)

## AHF documentation
For information on AMIGA's Halo Finder, including documentation and output file formats, you can visit their website [here](http://popia.ft.uam.es/AHF/Download.html).

## Gallery
![Velocity Distribution of Enclosed Halo as a Function of Redshift](https://raw.githubusercontent.com/BenDavisonPetch/ahfhalotools/main/gallery/enclosedHalovDist.png)
![Velocity Distribution of Enclosed Halo as a Function of Redshift with Skew Plots](https://raw.githubusercontent.com/BenDavisonPetch/ahfhalotools/main/gallery/enclosedHalovDistwSkewPlots.png)
![Local Density times Radius Squared as a function of z and r](https://raw.githubusercontent.com/BenDavisonPetch/ahfhalotools/main/gallery/localDensityZR.png)
![Local Density times Radius Squared as a function of z and r for select redshifts](https://raw.githubusercontent.com/BenDavisonPetch/ahfhalotools/main/gallery/selectzs.png)
