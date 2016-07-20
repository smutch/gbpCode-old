#! /bin/env python
""" This script generates a series of plots for analysing the results of an MCMC run."""
import sys
import os 
import numpy as np
from numpy import log10
import matplotlib
matplotlib.use('Agg')
import pylab as plt
from matplotlib.patches import Rectangle
 
sys.path.append(os.path.expandvars('$GBP_SRC/gbpMath/gbpMCMC'))
from MCMC_pyTools import MCMC

__author__    = "Simon Mutch"
__date__      = "25/10/11"
__copyright__ = "Simon Mutch"
__version__   = "1.0"
__email__     = "smutch@astro.swin.edu.au"
 
def MCMC_prob_plot(ax, ds_x, ds_y, run, n_integrate, 
        prob_type='posterior',
        swap=False, 
        legend=True,
        contours=True,
        scatter_points=True):
 
    # Read in the coverage map
    coverage = run.read_coverage_map(ds_x,ds_y, swap=swap)
    conf_68pc = coverage['conf_68pc']
    conf_95pc = coverage['conf_95pc']
    coverage_P = coverage['coverage_P']
    max_L = coverage['max_L']
    
    coverage_lo  =0
    coverage_hi  =10*n_integrate/(run.coverage_size*run.coverage_size)
 
    if swap==True:
        tmp=ds_x; ds_x=ds_y; ds_y=tmp
 
    coverage_P[np.where(coverage_P==0)] = -999
 
    # Plot the probabilities
    if prob_type=='posterior':
        im=ax.imshow(coverage_P,origin='lower',interpolation='bicubic',aspect='auto',
                        extent=[run.coverage_min[ds_x],run.coverage_max[ds_x],
                                run.coverage_min[ds_y],run.coverage_max[ds_y]], label='_none_', cmap=plt.cm.Blues)
        im.set_clim(coverage_lo,coverage_hi)
    elif prob_type=='max_L':
        im=ax.imshow(max_L,origin='lower',interpolation='bicubic',aspect='auto',
                        extent=[run.coverage_min[ds_x],run.coverage_max[ds_x],
                                run.coverage_min[ds_y],run.coverage_max[ds_y]], label='_none_', cmap=plt.cm.Blues)
        im.set_clim(0,1)
 
    im.get_cmap().set_under((1,1,1))
 
    # Plot the contours
    if contours:
        x=np.linspace(run.coverage_min[ds_x],run.coverage_max[ds_x],run.coverage_size)
        y=np.linspace(run.coverage_min[ds_y],run.coverage_max[ds_y],run.coverage_size)
 
        ax.contour(x,y,coverage_P,[conf_68pc,conf_95pc],colors='black', label='_none_')
 
    # Plot the starting and best fit parameter positions
    if scatter_points:
        ax.scatter(run.P_init[ds_x],run.P_init[ds_y],s=100,linewidth=3,c='green',
                    edgecolors='green', marker='+', zorder=100,
                    label="(%1.2e , %1.2e)"%(10.0**run.P_init[ds_x],10.0**run.P_init[ds_y]))
        ax.scatter(run.P_best[ds_x],run.P_best[ds_y],s=100,linewidth=1,c='red',
                    edgecolors='red', marker='o', zorder=100,
                    label="(%1.2e , %1.2e)"%(10.0**run.P_best[ds_x],10.0**run.P_best[ds_y]))
        ax.scatter(run.P_max_l[ds_x],run.P_max_l[ds_y],s=100,linewidth=1,c='orange',
                    edgecolors='red', marker='d', zorder=100,
                    label="(%1.2e , %1.2e)"%(10.0**run.P_max_l[ds_x],10.0**run.P_max_l[ds_y]))
 
    if legend==True:
        leg = ax.legend(loc='upper right', labelspacing=0.2,
                        scatterpoints=1, frameon=True)
        leg.legendPatch.set_alpha(0.2)
        leg.legendPatch.set_color('black')
        for t in leg.get_texts():
            t.set_size('small')
            t.set_color('black')
 
    return im
 
 
def MCMC_hist_plot(ax, i_DS, run):
 
    # Read in the histogram data
    (x, hist_p, hist_mean, hist_max, mask_68, 
            mask_95, best_val) = run.read_histogram(i_DS)
 
    # Plot the histogram
    ax.plot(x, hist_p, 'k-', linewidth=2)
    zero_line = np.linspace(0.,0.,run.coverage_size)
    ax.fill_between(x, zero_line, hist_p, where=mask_95>0, facecolor="#E0CA99",label=r'95% Confidence')
    ax.fill_between(x, zero_line, hist_p, where=mask_68>0, facecolor="#E39F0E",label=r'68% Confidence')
    ax.plot(x,hist_max, c='red', linewidth=1)
    ax.plot([best_val,best_val],[0.,1.],c='black',linewidth=2)
 
 
 
def plot_prob_matrix(mcmc_dir, limits='auto', prob_type='posterior',
                     scatter_points=None):
 
    mcmc_run = MCMC.MCMCrun(mcmc_dir)
    mcmc_run.read_best_fit_params()
    chain0 = MCMC.Chain(mcmc_run, 0, quiet=True)
 
    NPARAM = mcmc_run.n_P
 
    # We want a NPARAMxNPARAM grid of plots
    fig, fig_ax_arr = plt.subplots(NPARAM,NPARAM,figsize=(16,12))
    # Increase the space between plots to fit the labels in
    fig.subplots_adjust(hspace=0.0, wspace=0.0)
 
    # Set the background color to the same as the base color of the map
    plt.setp([ax for ax in fig.axes], axis_bgcolor="white")
    [ax.tick_params(labelsize='small') for ax in fig.axes]
 
    if limits=='priors':
        # Set the axis limits
        for j in xrange(NPARAM):
            for i in xrange(NPARAM):
                fig.axes[(j*NPARAM)+i].set_xlim((
                    mcmc_run.P_limit_min[i], mcmc_run.P_limit_max[i]))
                fig.axes[(j*NPARAM)+i].set_ylim((
                    mcmc_run.P_limit_min[j], mcmc_run.P_limit_max[j]))
 
    # Set the color of the ticks
    # [ax.tick_params(color="white", which="both") for ax in fig.axes]
 
    # Set the axis labels
    for i in xrange(NPARAM):
        fig.axes[((NPARAM-1)*NPARAM)+i].set_xlabel(mcmc_run.P_name[i])
        fig.axes[i*NPARAM].set_ylabel(mcmc_run.P_name[i])
 
    # Read the in the coverage map for each axis and plot the probability map
    for j in xrange(1,NPARAM):
        for i in xrange(j):
            im = MCMC_prob_plot(fig.axes[(j*NPARAM)+i],  
                    i, j, mcmc_run, chain0.n_integrate, legend=False,
                    prob_type=prob_type)
            if i>0 :
                plt.setp(fig.axes[(j*NPARAM)+i].get_yticklabels(), visible=False)
 
    # Plot the histograms
    for i in xrange(NPARAM):
        MCMC_hist_plot(fig.axes[(i*NPARAM)+i], i, mcmc_run)
 
    # If there are any scatter points to plot then do so...
    if scatter_points != None:
        for j in xrange(1,NPARAM):
            for i in xrange(j):
                fig.axes[(j*NPARAM)+i].scatter(scatter_points[:,i],
                                                scatter_points[:,j],
                                                color='purple',
                                                marker='x',)
        for i in xrange(NPARAM):
            for k in xrange(scatter_points.shape[0]):
                fig.axes[(i*NPARAM)+i].axvline(scatter_points[k,i],
                                               color='purple',
                                               lw=0.5)
 
    # If plot limits = priors was not requested then use the
    # histogram ranges to set all of the plot ranges.
    if limits=='auto':
        for j in xrange(NPARAM):
            for i in xrange(NPARAM):
                fig.axes[(j*NPARAM)+i].set_xlim(
                    fig.axes[(i*NPARAM)+i].get_xlim())
                fig.axes[(j*NPARAM)+i].set_ylim(
                    fig.axes[(j*NPARAM)+j].get_xlim())
 
    # Remove the axes from the plots we are not using
    for j in xrange(NPARAM):
        for i in xrange(j+1,NPARAM):
            fig.axes[(j*NPARAM)+i].set_axis_off()
 
    # Add a second y-axis to all of the probability plots
    for i in xrange(NPARAM):
        current_ax=fig.axes[(i*NPARAM)+i]
        current_ax.yaxis.set_ticks_position('right')
        current_ax.yaxis.set_label_position('right')
        current_ax.set_ylabel(r"$\mathrm{Prob}$")
        current_ax.yaxis.label.set_rotation(270)
        current_ax.set_ylim((0,1))
        current_ax.get_yticklabels()[0].set_visible(False)
 
    # Remove the first value from some of the y-axes
    for i in xrange(NPARAM-1):
        fig.axes[(i*NPARAM)].get_yticklabels()[0].set_visible(False)
 
    # Rotate the tick labels on the x-axes and remove some
    for i in xrange(NPARAM):
        labels = fig.axes[((NPARAM-1)*NPARAM)+i].get_xticklabels() 
        for label in labels: 
            label.set_rotation(90) 
        if i>0:
            labels[0].set_visible(False)
    for j in xrange(NPARAM-1):
        for i in xrange(NPARAM):
            plt.setp(fig.axes[(j*NPARAM)+i].get_xticklabels(), visible=False)
 
    # Plot the colour bar
    plt.subplots_adjust(bottom=0.2, top=0.95)
    cax = fig.add_subplot(111)
    cax.set_position([0.1, 0.05, 0.8, 0.05])
    plt.colorbar(im,cax=cax, orientation='horizontal')
    cax.text(0.5,1.2,"# of propositions", transform=cax.transAxes, horizontalalignment='center', verticalalignment='baseline', size='x-large')
 
    # Save the plot
    fig_fout = mcmc_run.filename_root+'./plots/prob_matrix.png'
    plt.savefig(fig_fout, transparent=False)
    print 'Saved file to', fig_fout
 
 
 
if __name__ == '__main__':
 
    mcmc_dir = sys.argv[1]+'/'
    #plot_prob_matrix(mcmc_dir, limits='priors', prob_type='posterior')
    plot_prob_matrix(mcmc_dir, limits='auto', prob_type='posterior')

