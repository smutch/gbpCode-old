"""Usage: view_grid_slice.py <grid_file> <slice_axis> (<fps> [--bitrate=<val>] | --single=<val>) [--novel]

-h, --help        show this help
--single=<val>    single slice at `val` h^-1 Mpc
--bitrate=<val>   bitrate of movie [default: 3600]
--novel           don't plot velocity quivers

Note that the movie generation may not work for some matplotlib backends.  On a
mac try qt4 for the best performance.
"""


import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import LogNorm

try:
    from astropy.utils.console import ProgressBar
except:
    print "The astropy package is not available."
    print "Please `pip install astropy` to have a progress bar showing status of movie creation."

try:
    from docopt import docopt
except:
    print "docopt is required to parse command line inputs..."
    print "Please `pip install docopt`"
    sys.exit(1)

# parse the command line options
arg = docopt(__doc__, help=True, version=1.0)
LABELS = np.array(['x','y','z'])
GRID_FNAME = arg['<grid_file>']
AXIS = np.where(LABELS==arg['<slice_axis>'])[0][0]
if arg['--single']:
    MOVIE = False
    INDEX = int(arg['--single']) # will be converted to a real index below
else:
    MOVIE = True
    FPS = int(arg['<fps>'])
    INDEX = 0
    BITRATE = int(arg['--bitrate'])
if arg['--novel']: 
    PLOTVEL=False
else:
    PLOTVEL=True

def read_grid_header(fin):
    """Read in the grid file header from the open file handle `fin`."""
    
    n_cell = np.fromfile(fin, 'i4', 3)
    box_size = np.fromfile(fin, 'f8', 3)
    n_grids = np.fromfile(fin, 'i4', 1)[0]
    ma_scheme = np.fromfile(fin, 'i4', 1)[0]
    
    print "-"*20
    print "n_cell    = ",n_cell
    print "box_size  = ",box_size
    print "n_grids   = ",n_grids
    print "ma_scheme = ", ma_scheme
    print "-"*20
    
    return n_cell, box_size, n_grids, ma_scheme

def read_grids(fname):
    """Read in a density grid from file `fname`."""
    
    grid = {}
    
    with open(fname, "rb") as fin:
        n_cell, box_size, n_grids, ma_scheme = read_grid_header(fin)
        n_elem = n_cell.cumprod()[-1]
        if not PLOTVEL:
            n_grids = 1
        for i_grid in xrange(n_grids):
            ident = np.fromfile(fin, 'S32', 1)[0]
            print("Reading grid "+ident)
            # For some reason numpy freaks out with a (1024,1024,1024)f4 read
            # so we need to just read in raw elements and reshape the result...
            grid[ident] = np.fromfile(fin, 'f4', n_elem)
            grid[ident].shape = n_cell
    
    return grid, n_cell, box_size


def animate(index):
    """This function is called for each new frame and updates the data being
    plotted."""

    # create a new slice
    s = list(np.s_[:,:,:])
    s[AXIS] = index

    # update the data in the density image and quivers
    cax.set_data(grid['rho_r_dark'][s])
    if PLOTVEL:
        quiver.set_UVC(grid['v_'+LABELS[(AXIS-1)%3]+'_r_dark'][s][quiver_s],
                       grid['v_'+LABELS[(AXIS+1)%3]+'_r_dark'][s][quiver_s])

    # update the figure title
    slice_pos = index*(cell_half_width[AXIS]*2.)+cell_half_width[AXIS]
    plt.title(LABELS[AXIS]+" = {:03.2f}".format(slice_pos)+r"h$^{-1}$ Mpc")

    # redraw
    fig.canvas.draw()
    
    # update the progress bar
    try:
        progressbar.update()
    except:
        pass



# read in the grids
grid, n_cell, box_size = read_grids(GRID_FNAME)

# We will use a log10 normalisation for the density field so multiply the
# densities by 1e10 to get h^-1 Msol units.
grid['rho_r_dark'] = grid['rho_r_dark']*1.e10
min_density, max_density = grid['rho_r_dark'][grid['rho_r_dark']>0].min(), grid['rho_r_dark'].max()

# mask the velocity fields to ignore zero density cells
if PLOTVEL:
    print "Masking the velocity fields to remove zero density cells..."
    bool_arr = grid['rho_r_dark'] < min_density
    for k in grid.iterkeys():
        if k is not 'rho_r_dark':
            grid[k] = np.ma.masked_where(bool_arr, grid[k], copy=False)

print "Plotting..."
# set up the figure
fig = plt.figure(0, figsize=(16, 12))
ax = plt.subplot(111)

# calculcate useful quantities
cell_half_width = box_size/n_cell/2.
plotted_axis = np.array(range(3), 'i8')
plotted_axis = plotted_axis[plotted_axis!=AXIS]
extent = (cell_half_width[plotted_axis[0]],
          box_size[plotted_axis[0]]-cell_half_width[plotted_axis[0]],
          cell_half_width[plotted_axis[1]],
          box_size[plotted_axis[1]]-cell_half_width[plotted_axis[1]])
sliced_grid = {}

if not MOVIE:
    INDEX = min(int(INDEX/(cell_half_width[AXIS]*2.)), n_cell[AXIS]-1)

# create the requested slice
s = list(np.s_[:,:,:])
s[AXIS] = INDEX

# plot the density grid
palette = plt.cm.jet
palette.set_bad('k')
cax = plt.imshow(grid['rho_r_dark'][s], origin='lower', extent=extent,
                 interpolation='bicubic',
                 vmin=min_density, vmax=max_density,
                 norm=LogNorm(),
                 cmap=palette)

# add a colorbar
cb = plt.colorbar(cax)

# add quivers for the velocity field
if PLOTVEL:
    X,Y = np.meshgrid(np.linspace(extent[0]+cell_half_width[plotted_axis[0]],
                                  extent[1]-cell_half_width[plotted_axis[0]],
                                  n_cell[plotted_axis[0]]),
                      np.linspace(extent[2]+cell_half_width[plotted_axis[1]],
                                  extent[3]-cell_half_width[plotted_axis[1]],
                                  n_cell[plotted_axis[1]]))
    quiver_s = list(np.s_[:,:])
    for i,a in enumerate(plotted_axis):
        if n_cell[a]>64:
            quiver_s[i] = np.s_[::int(n_cell[a]/64)]
    quiver = ax.quiver(X[quiver_s],Y[quiver_s],
                       grid['v_'+LABELS[(AXIS-1)%3]+'_r_dark'][s][quiver_s],
                       grid['v_'+LABELS[(AXIS+1)%3]+'_r_dark'][s][quiver_s],
                       units='xy', scale=500,
                       color='w', alpha=0.6,
                       pivot='tail',
                       headaxislength=5,
                       headwidth=2.5)

    # add a key for the quivers that indicates the normailsation of the lengths
    quiverkey = ax.quiverkey(quiver, 0.9, 0.93, 500, r'500 ks$^{-1}$', 
                             fontproperties={'weight': 'bold'}, labelcolor='w', alpha=1, coordinates='axes')

# set the axis labels and figure title
units = r"h$^{-1}$ [Mpc]"
plt.xlabel(LABELS[(AXIS+1)%3]+units)
plt.ylabel(LABELS[(AXIS-1)%3]+units)
cb.set_label(r"h$^{-1}$ [M$_{\odot}$ Mpc$^{-3}$]")
slice_pos = INDEX*(cell_half_width[AXIS]*2.)+cell_half_width[AXIS]
plt.title(LABELS[AXIS]+" = {:03.2f}".format(slice_pos)+r"h$^{-1}$ Mpc")


# save the image or create a movie by recursively calling the `animate` function
if not MOVIE:
    plt.draw()
    plt.savefig("slice_"+LABELS[AXIS]+"_{:04d}.png".format(INDEX))
else:
    print "Generating movie..."
    try:
        progressbar = ProgressBar(n_cell[AXIS]+1)
    except:
        pass

    anim = animation.FuncAnimation(fig, animate, frames=n_cell[AXIS], blit=True)
    # Save as mp4. This requires ffmpeg to be installed.
    anim.save('slice_'+LABELS[AXIS]+'.mp4', fps=FPS, bitrate=BITRATE)

    try:
        progressbar.__exit__(None,None,None)
    except:
        pass

