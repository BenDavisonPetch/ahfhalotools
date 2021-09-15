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

## Documentation
View HTML documentation __[here](https://htmlpreview.github.io/?https://raw.githubusercontent.com/BenDavisonPetch/ahfhalotools/main/docs/ahfhalotools/index.html)__

Documentation is also available in docstrings within the
code, which can be viewed using the built-in ``help`` function:

```python
    >>> from ahfhalotools.objects import Cluster
    >>> help(Cluster)
    Help on class Cluster in module ahfhalotools.objects: ...
```

The majority of analysis is enabled by the `ahfhalotools.objects.Cluster` object.

## Examples
Example scripts are available in `/examples`. These scripts are provided without data, as the data files are large and can't go on GitHub. To run them on a local machine, data must be downloaded and truncated, and the paths to the data in the scripts should be updated to reflect the location of the files.
Alternatively the scripts could be deployed to popia/castor and run there, after updating the directory paths in the code. Data should still be truncated before running, otherwise the scripts will execute very slowly.

## AHF documentation
For information on AMIGA's Halo Finder, including documentation and output file formats, you can visit their website [here](http://popia.ft.uam.es/AHF/Download.html).

## Gallery
![Velocity Distribution of Enclosed Halo as a Function of Redshift](https://raw.githubusercontent.com/BenDavisonPetch/ahfhalotools/main/gallery/enclosedHalovDist.png)
![Velocity Distribution of Enclosed Halo as a Function of Redshift with Skew Plots](https://raw.githubusercontent.com/BenDavisonPetch/ahfhalotools/main/gallery/enclosedHalovDistwSkewPlots.png)
![Local Density times Radius Squared as a function of z and r](https://raw.githubusercontent.com/BenDavisonPetch/ahfhalotools/main/gallery/localDensityZR.png)
![Local Density times Radius Squared as a function of z and r for select redshifts](https://raw.githubusercontent.com/BenDavisonPetch/ahfhalotools/main/gallery/selectzs.png)
![Local Gas Density as a function of Radius, with quartiles of dynamical state highlighted](https://raw.githubusercontent.com/BenDavisonPetch/ahfhalotools/main/clusterCompDensvR/quartile%20selection/gasDensMBP.png)
