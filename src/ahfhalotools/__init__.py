"""
A Python 3 library for the analysis of data produced by AMIGA's Halo Finder (AHF).

For more information visit https://github.com/BenDavisonPetch/ahfhalotools

The majority of analysis is done via the ahfhalotools.objects.Cluster class. For
examples on usage, there are example scripts available at:
https://github.com/BenDavisonPetch/ahfhalotools/tree/main/examples

Documentation is available as docstrings within the code, which can be viewed
using the built-in ``help`` function:

    >>> from ahfhalotools.objects import Cluster
    >>> help(Cluster)
    Help on class Cluster in module ahfhalotools.objects: ...
"""

__all__ = ["analysis","filetools","objects"]
