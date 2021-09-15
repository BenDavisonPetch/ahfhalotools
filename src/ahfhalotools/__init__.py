"""
A Python 3 library for the analysis of data produced by AMIGA's Halo Finder (AHF).

For more information visit https://github.com/BenDavisonPetch/ahfhalotools

Usage and Examples
------------------
The majority of analysis is done via the ahfhalotools.objects.Cluster class. For
examples on usage, there are example scripts available at:
https://github.com/BenDavisonPetch/ahfhalotools/tree/main/examples
These scripts are provided without data, as the data files are large and can't
go on GitHub. To run them on a local machine, data must be downloaded and
truncated, and the paths to the data in the scripts should be updated to reflect
the location of the files.
Alternatively the scripts could be deployed to popia/castor and run there, after
updating the directory paths in the code. Data should still be truncated before
running, otherwise the scripts will execute very slowly.

Documentation
-------------
Documentation is available in two places: the first is as docstrings within the
code, which can be viewed using the built-in ``help`` function:

    >>> from ahfhalotools.objects import Cluster
    >>> help(Cluster)
    Help on class Cluster in module ahfhalotools.objects: ...

The other way to view documentation is in a web browser at:
https://htmlpreview.github.io/?https://raw.githubusercontent.com/BenDavisonPetch/ahfhalotools/main/docs/ahfhalotools/index.html
"""

#__all__ = ["analysis","filetools","objects"]
