import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import WMAP9
from .. import analysis

'''========================================================='''
'''                     PLOT FUNCTIONS                      '''
'''========================================================='''

def plot4halocb(x_func,y_func,snaps,zs,title,x_label,y_label, debug = False):
    """
    Plots a 2x2 grid of profile data vs profile data for 4 most massive halos
    with a colorbar for z.

    .. deprecated:: 0.0.0
        plot4halocb uses outdated data format of a list of Snapshot instances.
        This results in inaccurate tracking of halos through time.
        Instead, plotting should be done manually, with data obtained via the
        ahfhalotools.objects.Cluster object

    Parameters
    ----------
    x_func : lambda
        lambda that should take Halo object and return array representing x
        values.
        e.g:
        lambda halo : halo.radii()
    y_func : lambda
        lambda that should take Halo object and return array representing y
        values.
    snaps : list of Snapshot instances
    zs : list of floats
        List of redshifts that correspond to snapshots in snaps
    title : str
        Title of plot
    x_label : str
        x axis label of plot
    y_label : str
        y axis label of plot
    """
    widths = [1,1,0.1]
    gs_kw = dict(width_ratios=widths)
    fig,ax = plt.subplots(2,3,figsize=(9,6),gridspec_kw=gs_kw)
    fig.suptitle(title,wrap=True)

    #setup color stuff
    normalize = mcolors.Normalize(vmin=zs.min(),vmax=zs.max())
    colormap = cm.plasma


    for i in range(4):
        #row and column indexes
        row = i//2
        col = i%2
        ax[row,col].set_title("Halo {0}".format(i+1))
        ax[row,col].set_xscale("log")
        ax[row,col].set_yscale("log")
        if row == 1: ax[row,col].set_xlabel(x_label)
        if col == 0: ax[row,col].set_ylabel(y_label)

        for snap in snaps:
            halo = snap.halos[i]
            x = x_func(halo)
            y = y_func(halo)

            #DEBUG
            if y.min() < 100 and debug:
                print("WEIRD DATA")
                print("y min is {0} for snapshot {1}".format(y.min(),snap.snapNo))
                print(x)
                print(y)

            #r,dens = removeNegatives(r,dens)

            ax[row,col].plot(x,y,label="$z={z:.3f}$".format(z=snap.z),color=colormap(normalize(snap.z)))

        #ax[row,col].set_xlim(0,5)
        #ax[row,col].set_ylim(0,300000000)

        #ax[row,col].legend()


    #---- colorbar legend ----
    scalarmappable = cm.ScalarMappable(norm=normalize, cmap=colormap)
    scalarmappable.set_array(zs)

    gs = ax[0,2].get_gridspec()
    #remove plots on rightmost column
    ax[0,2].remove()
    ax[1,2].remove()
    #create big axis for color bar
    cax = fig.add_subplot(gs[:,2])
    fig.colorbar(scalarmappable, cax=cax)
    cax.set_ylabel("redshift, $z$")


    plt.tight_layout()
    plt.show()

def plot4halovtimeradinterp(x_func,v_func,snaps,zs,title,x_label,y_label,radii):
    """
    Plots a 2x2 grid of profile data vs time for 4 most massive halos.

    .. deprecated:: 0.0.0
        plot4halovtimeradinterp uses the outdated data format of a list of
        Snapshot instances.
        This results in inaccurate tracking of halos through time.
        Instead, plotting should be done manually, with data obtained via the
        ahfhalotools.objects.Cluster object

    Parameters
    ----------
    x_func : lambda
        lambda that should take list of redshifts (floats) and return array
        representing x values.
        e.g:
        lambda zs : analysis.tfromz(z)
    v_func : lambda
        lambda that should take Halo object and return array representing
        profile data in question as a *function of radius*
    snaps : list of Snapshot instances
    zs : list of floats
        List of redshifts that correspond to snapshots in snaps
    title : str
        Title of plot
    x_label : str
        x axis label of plot
    y_label : str
        y axis label of plot
    radii : list of floats
        List of radii to interpolate at
    """
    fig,ax = plt.subplots(2,2,figsize=(9,6))
    fig.suptitle(title,wrap=True)

    for i in range(4):
        #row and column indexes
        row = i//2
        col = i%2
        #calculate expected timescale from halo at z=0
        halo = snaps[0].halos[i]
        sigV = halo.sigV * 1000 #km/s to m/s
        Rvir = halo.Rvir * 3.086e+19 / (WMAP9.H(0).value/100) #kpc/h to m
        timescale = Rvir/sigV / (60**2 * 24 * 365.26 * 1e9)

        ax[row,col].set_title("Halo {hid}, expected timescale {t:.3f}Gyr".format(hid=i+1,t=timescale))
        #ax[row,col].set_xscale("log")
        ax[row,col].set_yscale("log")
        if row == 1: ax[row,col].set_xlabel(x_label)
        if col == 0: ax[row,col].set_ylabel(y_label)

        for radius in radii:
            x = x_func(zs)
            y = analysis.getAsFuncOfTimeAtRadius(snaps,i,v_func,radius)
            ax[row,col].plot(x,y, label="r = {0}".format(radius))


        #ax[row,col].set_xlim(0,5)
        #ax[row,col].set_ylim(0,300000000)
        ax[row,col].legend()


    plt.tight_layout()
    plt.show()

def plot4halovtime(x_func,v_funcs,v_names,snaps,zs,title,x_label,y_label, yscale = "log",
                   ls = None):
    """
    Plots a 2x2 grid of halo data vs time for 4 most massive halos.

    .. deprecated:: 0.0.0
        plot4halovtime uses the outdated data format of a list of
        Snapshot instances.
        This results in inaccurate tracking of halos through time.
        Instead, plotting should be done manually, with data obtained via the
        ahfhalotools.objects.Cluster object
    """
    #plots a halo data piece against time for 4 different halos
    #x_func is lambda that takes array of zs and returns x-values to plot against
    #v_funcs is array of lambdas that takes halo and returns halo data value being plotted on y
    #v_names is names of values from v_funcs
    if not ls: ls = [""]*len(v_funcs)

    fig,ax = plt.subplots(2,2,figsize=(9,6))
    fig.suptitle(title,wrap=True)

    for i in range(4):
        #row and column indexes
        row = i//2
        col = i%2

        ax[row,col].set_title("Halo {hid}".format(hid=i+1))
        #ax[row,col].set_xscale("log")
        ax[row,col].set_yscale(yscale)
        if row == 1: ax[row,col].set_xlabel(x_label)
        if col == 0: ax[row,col].set_ylabel(y_label)

        x = x_func(zs)
        for j in range(len(v_funcs)):
            v_func = v_funcs[j]
            y = [v_func(snap.halos[i]) for snap in snaps]
            ax[row,col].plot(x,y,ls[j],label=v_names[j])


        #ax[row,col].set_xlim(0,5)
        #ax[row,col].set_ylim(0,300000000)
        if len(v_funcs) > 1: ax[row,col].legend(loc="upper right")

    plt.tight_layout()
    plt.show()

def plot1halovtimecompfuncs(haloid,x_func,v_funcs,value_names,snaps,zs,title,x_label,y_label,radii):
    """
    Plots a 2x2 grid of different pieces of halo data vs time for a particular
    halo.

    .. deprecated:: 0.0.0
        plot1halovtimecompfuncs uses the outdated data format of a list of
        Snapshot instances.
        This results in inaccurate tracking of halos through time.
        Instead, plotting should be done manually, with data obtained via the
        ahfhalotools.objects.Cluster object
    """
    #radii is array of radii to test at
    #x_func is lambda that takes array of zs and returns x-values to plot against
    #v_funcs is array of lambdas that takes halo and returns array of value in question as a
    #  function of **radius**
    fig,ax = plt.subplots(2,2,figsize=(9,6))
    fig.suptitle(title,wrap=True)

    halo = snaps[0].halos[haloid]
    #timescale calc
    sigV = halo.sigV * 1000 #km/s to m/s
    Rvir = halo.Rvir * 3.086e+19 / (WMAP9.H(0).value/100) #kpc/h to m
    timescale = Rvir/sigV / (60**2 * 24 * 365.26 * 1e9)
    x = x_func(zs)

    for i in range(4):
        #row and column indexes
        row = i//2
        col = i%2
        radius = radii[i]

        ax[row,col].set_title("r = {radius}, expected timescale {t:.3f}Gyr".format(radius=radius,t=timescale))
        #ax[row,col].set_xscale("log")
        ax[row,col].set_yscale("log")
        if row == 1: ax[row,col].set_xlabel(x_label)
        if col == 0: ax[row,col].set_ylabel(y_label)

        for j in range(len(v_funcs)):
            y = analysis.getAsFuncOfTimeAtRadius(snaps,haloid,v_funcs[j],radius)
            ax[row,col].plot(x,y, label=value_names[j])


        #ax[row,col].set_xlim(0,5)
        #ax[row,col].set_ylim(0,300000000)
        ax[row,col].legend()


    plt.tight_layout()
    plt.show()

def plothalovtimecompsims_radinterp(ax, haloid,x_func,v_func,sims,sim_names,title,x_label,y_label,radius):
    """
    Takes a matplotlib.Axes.axes object and plots onto it a piece of profile
    data, interpolated at a specified radius, as a function of time, for each
    simulation.

    .. deprecated:: 0.0.0
        plothalovtimecompsims_radinterp uses the outdated data format of a list
        of Snapshot instances.
        This results in inaccurate tracking of halos through time.
        Instead, plotting should be done manually, with data obtained via the
        ahfhalotools.objects.Cluster object
    """
    #ax is axis object that is plotted onto
    #x_func is lambda that takes array of zs and returns x-values to plot against
    #v_func is lambda that takes halo and returns array of value in question as a
    #  function of **radius**
    #sims is an array of arrays of snaps corresponding to each simulation
    ax.set_yscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

    for i in range(len(sims)):
        snapshots = sims[i]
        zs = np.array([snapshot.z for snapshot in snapshots])
        #NOTE: above line is not very efficient, if scaled to more data just
        #    give z arrays in arguments of function
        x = x_func(zs)
        y = analysis.getAsFuncOfTimeAtRadius(snapshots,haloid,v_func,radius)
        ax.plot(x,y, label=sim_names[i])

        ax.legend()
    #plt.tight_layout()
    #plt.show()

def plothalovtimecompsims(haloid,x_func,v_func,sims,sim_names,title,x_label,y_label,yscale="log"):
    """
    Plots a graph of halo data as a function of time for each simulation.

    .. deprecated:: 0.0.0
        plothalovtimecompsims uses the outdated data format of a list of
        Snapshot instances.
        This results in inaccurate tracking of halos through time.
        Instead, plotting should be done manually, with data obtained via the
        ahfhalotools.objects.Cluster object
    """
    #x_func is lambda that takes array of zs and returns x-values to plot against
    #v_func is lambda that takes halo and returns halodata being plotted
    #sims is an array of arrays of snaps corresponding to each simulation
    fig,ax = plt.subplots(1,1,figsize=(9,6))
    fig.suptitle(title,wrap=True)

    #ax[row,col].set_xscale("log")
    ax.set_yscale(yscale)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    for i in range(len(sims)):
        snapshots = sims[i]
        zs = np.array([snapshot.z for snapshot in snapshots])
        #NOTE: above line is not very efficient, if scaled to more data just
        #    give z arrays in arguments of function
        x = x_func(zs)
        y = [v_func(snap.halos[haloid]) for snap in snapshots]
        ax.plot(x,y, label=sim_names[i])

        ax.legend()
    plt.tight_layout()
    plt.show()
