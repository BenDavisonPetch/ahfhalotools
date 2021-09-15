"""
Some useful methods for making 3D plots of simulations
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib

'''========================================================='''
'''                       3D PLOTTING                       '''
'''========================================================='''
class Sphere:
    """
    An object that allows for 3D plotting of a sphere. Includes
    functionality to plot sphere representing a halo in 3D

    Attributes
    ----------
    x : np.array of dtype float
        2D Grid of x coordinates of vertexes
    y : np.array of dtype float
        2D Grid of y coordinates of vertexes
    z : np.array of dtype float
        2D Grid of z coordinates of vertexes
    color : str
    alpha : float
    label : str
    """
    def __init__(self,xc,yc,zc,r, ures=20 , vres=10 ,color='',alpha=1,label=""):
        """
        Initialises Sphere instance

        Parameters
        ----------
        xc : float
            x coordinate of centre of sphere
        yc : float
            y coordinate of centre of sphere
        zc : float
            z coordinate of centre of sphere
        r : float
            radius of sphere

        Other parameters
        ----------------
        ures : int
            Sets resolution of mesh. A higher ures -> more vertexes.
            Defaults to 20
        vres : int
            Sets resolution of mesh. A higher ures -> more vertexes.
            Defaults to 10
        color : str, optional
        alpha : float
            Defaults to 1
        label : str, optional
            Label to use for legends
        """
        #stores list of vertexes in self.x,self.y,self.z
        self.color=color
        self.alpha=alpha
        self.label=label

        #code provided from
        # https://stackoverflow.com/questions/24123659/scatter-plot-3d-with-labels-and-spheres

        #create matrix of points for unit sphere centred at origin
        u,v = np.mgrid[0:2*np.pi:ures*1j,0:np.pi:vres*1j]
        x=np.cos(u)*np.sin(v)
        y=np.sin(u)*np.sin(v)
        z=np.cos(v)

        #scale
        self.x = r*x + xc
        self.y = r*y + yc
        self.z = r*z + zc

    @classmethod
    def fromHalo(cls, halo,ures=20,vres=10,color='',alpha=1,label=""):
        """
        Initialises a Sphere object from a Halo object

        Parameters
        ----------
        halo : Halo instance

        Other Parameters
        ----------------
        ures : int
            Defaults to 20
            Specifies resolution of sphere. Higher resolution -> more
            vertexes
        vres : int
            Defaults to 10
        color : str, optional
        alpha : float
            Defaults to 1
        label : str, optional

        Returns
        -------
        sphere : Sphere instance
        """
        xc = halo.Xc
        yc = halo.Yc
        zc = halo.Zc
        Rvir = halo.Rvir
        return cls(xc,yc,zc,Rvir,
                ures=ures,vres=vres,color=color,alpha=alpha,label=label)


def plotSpheres(ax,spheres,style='surface'):
    """
    Takes array of Sphere objects ands plots them on a 3D axis

    Parameters
    ----------
    ax : matplotlib.Axes.axes instance
        The 3D axes instance to plot the sphere onto
    sphere : list of Sphere instances
        List of Spheres to plot
    style : str
        Can either specify style as "surface" or "wireframe"
        Defaults to "surface"

    Returns
    -------
    children : list of mpl_toolkits.mplot3d.art3d.Poly3DCollection objects
    """
    children = []

    if style=='surface':
        for sphere in spheres:
            artist=ax.plot_surface(sphere.x,sphere.y,sphere.z,color=sphere.color,
                    alpha=sphere.alpha,label=sphere.label)

            #put in some code to fix error described by
            #https://stackoverflow.com/questions/55531760/is-there-a-way-to-label-multiple-3d-surfaces-in-matplotlib/55534939
            artist._facecolors2d = artist._facecolor3d
            artist._edgecolors2d = artist._edgecolor3d

            #add artist to list of children to return
            children.append(artist)
        return children

    elif style=='wireframe':
        for sphere in spheres:
            artist=ax.plot_wireframe(sphere.x,sphere.y,sphere.z,color=sphere.color,
                    alpha=sphere.alpha,label=sphere.label)

            #put in some code to fix error described by
            #https://stackoverflow.com/questions/55531760/is-there-a-way-to-label-multiple-3d-surfaces-in-matplotlib/55534939
            artist._facecolors2d = artist._facecolor3d
            artist._edgecolors2d = artist._edgecolor3d

            #add artist to list of children to return
            children.append(artist)
        return children

    else:
        print("WARNING: invalid style for plotSpheres(), valid styles are\
            'surface' or 'wireframe'")
        return None
