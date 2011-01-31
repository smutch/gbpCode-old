#! /usr/bin/env python2.6
import sys
from pylab import *
from matplotlib.font_manager import FontProperties
from string import *

# Todo: Add/read array info to/from run.dat
#       Add array names to column headers of results files
#       Fix position of max likeliood lines in plots
#       Clean-up prop-report of compute_MCMC and add prop_report switch
#       Add loops over n_chains in this code

# Check syntax
if(len(sys.argv)!=2):
	print
	print "Syntax:",sys.argv[0],"MCMC_directory"
	print "-------"
	sys.exit(101)
filename_root = str(sys.argv[1])+'/'

# User editable stuff
flag_results_on    =1
flag_datasets_on   =1
flag_histograms_on =1
flag_coverage_on   =1
flag_trace_on      =1
my_chain           =0
P_results_ordinate =0
DS_results_ordinate=0
color_95pc='#90EE90'
color_68pc='#226622'

# These should not change
MCMC_NAME_SIZE    =256
filename_plot_root=filename_root+'/plots/'

# Read (and print) header info 
fd                  =open(filename_root+'/run.dat','rb')
problem_name        =fromfile(file=fd, dtype=(str_, MCMC_NAME_SIZE), count=1)[0].split('\x00')[0]
n_avg               =fromfile(file=fd,dtype='i',count=1)
n_avg_covariance    =fromfile(file=fd,dtype='i',count=1)
flag_autocor_on_file=fromfile(file=fd,dtype='i',count=1)
print 'Problem name=',problem_name
print 'n_avg       =',n_avg
print 'flag_autocor=',flag_autocor_on_file

# Read (and print) parameter names to stdout
print 'MCMC Parameters:'
n_P_in     =fromfile(file=fd,dtype='i',count=1)
n_P        =n_P_in[0]
P_name     = []
P_init     = []
P_limit_min= []
P_limit_max= []
for i_P in xrange(0,n_P):
	P_name.append(fromfile(file=fd, dtype=(str_, MCMC_NAME_SIZE), count=1)[0].split('\x00')[0])
	P_init.append(fromfile(file=fd,dtype='d',count=1))
	P_limit_min.append(fromfile(file=fd,dtype='d',count=1))
	P_limit_max.append(fromfile(file=fd,dtype='d',count=1))
	print '  ',P_name[i_P]
print

# Read (and print) parameter array names
n_P_arrays_in=fromfile(file=fd,dtype='i',count=1)
n_P_arrays   =n_P_arrays_in[0]
if(n_P_arrays>0):
	print 'MCMC Parameter array(s):'
	P_array_name= []
	for i_array in xrange(0,n_P_arrays):
		P_array_name.append(fromfile(file=fd, dtype=(str_, MCMC_NAME_SIZE), count=1)[0].split('\x00')[0])
		junk=fromfile(file=fd,dtype='d',count=n_P)
		print '  ',P_array_name[i_array]
	print


# Read (and print) the dataset names
print 'MCMC Datasets:'
n_DS_in  =fromfile(file=fd,dtype='i',count=1)
n_DS     =n_DS_in[0]
n_M_total=sum(n_M)
DS_name  =[]
M_target =[]
dM_target=[]
n_M      =[]
DS_array_name    = []
n_DS_arrays      = []
DS_arrays_offset = []
n_DS_arrays_total= 0
for i_DS in xrange(0,n_DS):
	DS_name.append(fromfile(file=fd, dtype=(str_, MCMC_NAME_SIZE), count=1)[0].split('\x00')[0])
	n_M.append(fromfile(file=fd,dtype='i',count=1))
	M_target.append(fromfile(file=fd,dtype='d',count=n_M[i_DS]))
	dM_target.append(fromfile(file=fd,dtype='d',count=n_M[i_DS]))
	n_DS_arrays.append(fromfile(file=fd,dtype='i',count=1))
	print '  ',DS_name[i_DS]
	if(n_DS_arrays[i_DS]>0):
        	for i_array in xrange(0,n_DS_arrays[i_DS]):
	                DS_array_name.append(fromfile(file=fd, dtype=(str_, MCMC_NAME_SIZE), count=1)[0].split('\x00')[0])
			junk=fromfile(file=fd,dtype='d',count=n_M[i_DS])
        		print '     array #'+str(i_array+1).zfill(2)+': '+DS_array_name[n_DS_arrays_total+i_array]
        DS_arrays_offset.append(n_DS_arrays_total)
	n_DS_arrays_total=n_DS_arrays_total+n_DS_arrays[i_DS]
print
fd.close()

# Read (and print) the number of iterations
fd                    =open(filename_root+'/chains/chain_iterations_'+str(my_chain).zfill(6)+'.dat','rb')
n_iterations          =fromfile(file=fd,dtype='i',count=1)
n_iterations_burn     =fromfile(file=fd,dtype='i',count=1)
n_iterations_integrate=n_iterations-n_iterations_burn
n_burn                =n_avg*n_iterations_burn
n_integrate           =n_avg*n_iterations_integrate
n_total               =n_burn+n_integrate
fd.close()
print 'n_burn       =',n_burn
print 'n_integrate  =',n_integrate
print

# Read coverage information
fd=open(filename_root+'/results/coverage.dat','rb')
n_coverage   =fromfile(file=fd,dtype='i',count=1)
coverage_size=fromfile(file=fd,dtype='i',count=1)
coverage_min =fromfile(file=fd,dtype='d',count=n_P)
coverage_max =fromfile(file=fd,dtype='d',count=n_P)
fd.close()

# Plot results
if (flag_results_on==1 and n_P_arrays>0):
	print "Creating results plots:"
	lx       = []
	lP_init  = []
	lP_best  = []
	lP_lo_68 = []
	lP_hi_68 = []
	lP_lo_95 = []
	lP_hi_95 = []
	for line in file(filename_root+'/results/fit_for_parameters.dat'):
		line = line.split()
		if(line[0][0]!='#'):
			lx.append(line[1])
			lP_init.append(line[2])
			lP_best.append(line[5])
			lP_lo_68.append(line[6])
			lP_hi_68.append(line[7])
			lP_lo_95.append(line[8])
			lP_hi_95.append(line[9])
	x      =array(lx,dtype='d')
        P_init =array(lP_init,dtype='d')
        P_best =array(lP_best,dtype='d')
        P_lo_68=array(lP_lo_68,dtype='d')
        P_hi_68=array(lP_hi_68,dtype='d')
        P_lo_95=array(lP_lo_95,dtype='d')
        P_hi_95=array(lP_hi_95,dtype='d')
	# Create plot here
        plt.figure()
        plt.cla()
        ax = plt.subplot(111)
	plt.axis([4.0e-3,5.0e-1,1e2,2e5],linewidth=2,fontsize='large')
	plt.ylabel(problem_name,fontsize='large')
	plt.xlabel(P_array_name[P_results_ordinate],fontsize='large')
	gca().yaxis.get_major_formatter().set_powerlimits((-1,1))
	ax.set_xscale('log')
	ax.set_yscale('log')
	p0=ax.fill_between(x, P_lo_95, P_hi_95, where=P_hi_95>P_lo_95, facecolor=color_95pc,label=r'95% Confidence')
	p0=ax.fill_between(x, P_lo_68, P_hi_68, where=P_hi_68>P_lo_68, facecolor=color_68pc,label=r'68% Confidence')
	p1=ax.plot(x,P_best,c='black',linewidth=2,label=r'Best fit')
	p2=ax.plot(x,P_init,c='blue',  linewidth=2,label=r'Smith et al z=0.6, WMAP-5')
	#plt.legend((p1,p2),(problem_name,'Initial values'),shadow = True)
        filename_out = 'fit_for_parameters.png'
        print '  Writing',filename_out
        filename_out =filename_plot_root+filename_out
        plt.savefig(filename_out)
	print

# Plot dataset fits
if (flag_datasets_on==1 and n_DS_arrays_total>0):
	print "Creating dataset-fit plots:"
	for i_DS in xrange(0,n_DS):
		if(n_DS_arrays[i_DS]>0):
		        lx       = []
        		lM_DS    = []
    	    		ldM_DS   = []
	       	 	lM_best  = []
       		 	lM_lo_68 = []
       		 	lM_hi_68 = []
       		 	lM_lo_95 = []
       		 	lM_hi_95 = []
			if(n_DS>1):
				filename_read=filename_root+'/results/fit_for_dataset_'+str(i_DS).zfill(5)+'.dat'
			else:
				filename_read=filename_root+'/results/fit_for_dataset.dat'
        		for line in file(filename_read):
                		line = line.split()
				if(line[0][0]!='#'):
        	        		lx.append(line[0])
                			lM_DS.append(line[1])
                			ldM_DS.append(line[2])
                			lM_best.append(line[5])
	                		lM_lo_68.append(line[6])
        	        		lM_hi_68.append(line[7])
                			lM_lo_95.append(line[8])
					lM_hi_95.append(line[9])
	                x      =array(lx,dtype='d')
        	        M_DS   =array(lM_DS,dtype='d')
                	dM_DS  =array(ldM_DS,dtype='d')
	                M_best =array(lM_best,dtype='d')
        	        M_lo_68=array(lM_lo_68,dtype='d')
                	M_hi_68=array(lM_hi_68,dtype='d')
	                M_lo_95=array(lM_lo_95,dtype='d')
        	        M_hi_95=array(lM_hi_95,dtype='d')
			dM_hi_DS=dM_DS
			dM_lo_DS=dM_DS
			for i_M in xrange(0,n_M[i_DS]):
				if((M_DS[i_M]-dM_lo_DS[i_M])<0):
					dM_lo_DS[i_M]=0.99999*M_DS[i_M]
	        	fd.close()
        		plt.figure()
        		# Create plot here
	        	plt.cla()
	        	ax = plt.subplot(111)
			plt.axis([4.0e-3,5.0e-1,1e2,2e5],linewidth=2,fontsize='large')
			plt.ylabel(DS_name[i_DS],fontsize='large')
			plt.xlabel(DS_array_name[DS_arrays_offset[i_DS]+DS_results_ordinate],fontsize='large')
			gca().yaxis.get_major_formatter().set_powerlimits((-1,1))
			ax.set_xscale('log')
		        ax.set_yscale('log')
	       	 	p0=ax.fill_between(x, M_lo_95, M_hi_95, where=M_hi_95>M_lo_95, facecolor=color_95pc,label=r'95% Confidence')
       		 	p0=ax.fill_between(x, M_lo_68, M_hi_68, where=M_hi_68>M_lo_68, facecolor=color_68pc,label=r'68% Confidence')
       		 	d1=errorbar(x,M_DS,[dM_lo_DS,dM_hi_DS],ecolor='red',fmt=None,linewidth=2)
			dummy1=[1e-20,2e-20]
			dummy2=[1e-20,2e-20]
			d1a   =plot(dummy1,dummy2,c='red')
        		p1    =ax.plot(x,M_best,c='black',linewidth=2,label=r'Best fit')
			# Create plot here
			if(n_DS>1):
	   	    	 	filename_out = 'fit_for_dataset_'+str(i_DS)+'.png'
			else:
				filename_out = 'fit_for_dataset.png'
                	print '  Writing',filename_out
                	filename_out =filename_plot_root+filename_out
                	plt.savefig(filename_out)	
        print

# Create histogram plots
if (flag_histograms_on==1):
        print "Creating histogram plots:"
        fd=open(filename_root+'/results/histograms.dat','rb')
        ii=0
        for jj in xrange(0,n_P):
                  best_val  =fromfile(file=fd,dtype='d',count=1)
                  lo_68pc   =fromfile(file=fd,dtype='d',count=1)
                  hi_68pc   =fromfile(file=fd,dtype='d',count=1)
                  lo_95pc   =fromfile(file=fd,dtype='d',count=1)
                  hi_95pc   =fromfile(file=fd,dtype='d',count=1)
                  hist_P    =fromfile(file=fd,dtype='int64',count=coverage_size[0])
		  x         =linspace(coverage_min[jj],coverage_max[jj],coverage_size)
		  mask_68   =zeros(coverage_size[0])
		  mask_95   =zeros(coverage_size[0])
		  for kk in xrange(0,coverage_size):
			if(x[kk]>=lo_68pc and x[kk]<=hi_68pc):
				mask_68[kk]=1
                        if(x[kk]>=lo_95pc and x[kk]<=hi_95pc):
                                mask_95[kk]=1
		  zero_line =linspace(0.,0.,coverage_size)
		  hist_P=double(hist_P)/double(hist_P.max())

                  # Create plot
                  plt.figure()
                  plt.cla()
                  ax = plt.subplot(111)
		  p0=ax.plot(x,hist_P,c='black',linewidth=2)
                  p2=ax.fill_between(x, zero_line, hist_P, where=mask_95>0, facecolor=color_95pc,label=r'95% Confidence')
		  p1=ax.fill_between(x, zero_line, hist_P, where=mask_68>0, facecolor=color_68pc,label=r'68% Confidence')
		  plot([best_val[0],best_val[0]],[0.,1.],c='black',linewidth=2)
                  ax.set_xlim((coverage_min[jj],coverage_max[jj]))
                  ax.set_ylim((0,1.1))
                  ax.set_xlabel(P_name[jj])
                  ax.set_ylabel('Probability')
                  ax.set_title(problem_name)
                  filename_out = 'histogram_'+str(jj).zfill(5)+'.png'
                  print '  Writing',filename_out
                  filename_out =filename_plot_root+filename_out
                  plt.savefig(filename_out)
	fd.close()
	print

# Create coverage maps
if (flag_coverage_on==1):
	print "Creating",n_coverage,"coverage plots:"
	fd=open(filename_root+'/results/coverage.dat','rb')
	n_coverage   =fromfile(file=fd,dtype='i',count=1)
	coverage_size=fromfile(file=fd,dtype='i',count=1)
	coverage_min =fromfile(file=fd,dtype='d',count=n_P)
	coverage_max =fromfile(file=fd,dtype='d',count=n_P)
	coverage_lo  =0
	coverage_hi  =10*n_integrate/(coverage_size[0]*coverage_size[0])
        ii=0
        for jj in xrange(0,n_P):
                for kk in xrange(jj+1,n_P):
		  x=linspace(coverage_min[jj],coverage_max[jj],coverage_size)
		  y=linspace(coverage_min[kk],coverage_max[kk],coverage_size)
		  conf_68pc =fromfile(file=fd,dtype='d',count=1)
		  conf_95pc =fromfile(file=fd,dtype='d',count=1)
                  coverage_1=fromfile(file=fd,dtype='int64',count=coverage_size[0]*coverage_size[0])
                  coverage_0=fromfile(file=fd,dtype='int64',count=coverage_size[0]*coverage_size[0])
                  coverage_P=fromfile(file=fd,dtype='int64',count=coverage_size[0]*coverage_size[0])
                  coverage_1=transpose(reshape(coverage_1,[coverage_size[0],coverage_size[0]]))
                  coverage_0=transpose(reshape(coverage_0,[coverage_size[0],coverage_size[0]]))
                  coverage_P=transpose(reshape(coverage_P,[coverage_size[0],coverage_size[0]]))

		  # Create probability maps
                  plt.figure(facecolor='white')
                  plt.cla()
                  ax = plt.subplot(111)
                  im=plt.imshow(coverage_P,origin='lower',interpolation='bicubic',aspect='auto',extent=[coverage_min[jj],coverage_max[jj],coverage_min[kk],coverage_max[kk]])
                  im.set_clim(coverage_lo,coverage_hi)
                  cont2 = plt.contour(x,y,coverage_P,[conf_68pc,conf_95pc],colors='red')
		  pp1=ax.scatter(P_init[jj],P_init[kk],s=250,linewidth=3,c='lawngreen',marker='+')
		  ax.set_xlim((coverage_min[jj],coverage_max[jj]))
		  ax.set_ylim((coverage_min[kk],coverage_max[kk]))
                  ax.set_xlabel(P_name[jj])
                  ax.set_ylabel(P_name[kk])
                  ax.set_title('Probability Map')
                  cb=plt.colorbar(im)
		  cb.ax.set_ylabel('# of propositions')
                  filename_out = 'probable_'+str(ii).zfill(5)+'.png'
                  print '  Writing',filename_out
                  filename_out =filename_plot_root+filename_out
                  plt.savefig(filename_out)

		  # Create coverage maps
                  plt.figure(figsize=(8,5))
                  plt.clf()
                  ax = axes([0.125,0.225,0.4,0.7],frameon=True)
		  ax.set_xlim((coverage_min[jj],coverage_max[jj]))
		  ax.set_ylim((coverage_min[kk],coverage_max[kk]))
		  ax.xaxis.set_major_locator(MaxNLocator(4))
		  ax.yaxis.set_major_locator(MaxNLocator(4))
                  im=plt.imshow(coverage_1,origin='lower',interpolation='bicubic',aspect='auto',extent=[coverage_min[jj],coverage_max[jj],coverage_min[kk],coverage_max[kk]])
                  im.set_clim(coverage_lo,coverage_hi)
                  ax.images.append(im)
		  cont2 = plt.contour(x,y,coverage_P,[conf_68pc,conf_95pc],colors='red')
                  plt.rcParams['font.size'] = 12.
                  ax.set_xlabel(P_name[jj])
                  ax.set_ylabel(P_name[kk])
                  ax.set_title('Accepted Propositions')
                  ax = axes([0.525,0.225,0.4,0.7],frameon=True)
		  ax.set_xlim((coverage_min[jj],coverage_max[jj]))
		  ax.set_ylim((coverage_min[kk],coverage_max[kk]))
		  ax.xaxis.set_major_locator(MaxNLocator(4))
		  ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useOffset=False))
		  setp(ax.get_yticklabels(),visible=False)
	 	  im=plt.imshow(coverage_0,origin='lower',interpolation='bicubic',aspect='auto',extent=[coverage_min[jj],coverage_max[jj],coverage_min[kk],coverage_max[kk]])
                  im.set_clim(coverage_lo,coverage_hi)
		  cont2 = plt.contour(x,y,coverage_P,[conf_68pc,conf_95pc],colors='red')
                  ax.images.append(im)
                  plt.rcParams['font.size'] = 12.
                  ax.set_xlabel(P_name[jj])
                  ax.set_title('Rejected Propositions')
                  cax = axes([0.125,0.085,0.8,0.025],frameon=False)
                  cb=plt.colorbar(im,cax=cax,orientation='horizontal')
                  cb.ax.set_xlabel('# of propositions')
                  filename_out = 'coverage_'+str(ii).zfill(5)+'.png'
                  print '  Writing',filename_out
                  filename_out =filename_plot_root+filename_out
                  plt.savefig(filename_out,transparent=False)
                  ii=ii+1
	fd.close()
        print

#Create trace plots
if (flag_trace_on==1):
        print "Creating trace plots:"
        fd2=open(filename_root+'/results/histograms.dat','rb')
        ii=0
	nullfmt        = NullFormatter()
	left, width    = 0.15, 0.65
	bottom, height = 0.15, 0.65
	bottom_h = left_h = left+width+0.02
	rect_chi2  = [left, bottom, 0.8, height]
	rect_trace = [left, bottom, width, height]
	rect_hist  = [left_h, bottom, 0.15, height]
        fd=open(filename_root+'/chains/chain_stats_000000.dat','rb')
	P_min_p    = zeros((n_iterations,n_P))
	P_avg_p    = zeros((n_iterations,n_P))
	P_max_p    = zeros((n_iterations,n_P))
	dP_avg_p   = zeros((n_iterations,n_P))
	dP_sub_p   = zeros((n_iterations,n_P))
	ln_Pr_min_p= zeros((n_iterations))
	ln_Pr_avg_p= zeros((n_iterations))
	ln_Pr_max_p= zeros((n_iterations))
        P_min_c    = zeros((n_iterations,n_P))
        P_avg_c    = zeros((n_iterations,n_P))
        P_max_c    = zeros((n_iterations,n_P))
        dP_avg_c   = zeros((n_iterations,n_P))
        dP_sub_c   = zeros((n_iterations,n_P))
        ln_Pr_min_c= zeros((n_iterations))
        ln_Pr_avg_c= zeros((n_iterations))
        ln_Pr_max_c= zeros((n_iterations))
	slopes     = zeros((n_iterations,n_P))
	drift      = zeros((n_iterations,n_P))
	x          = linspace(1,n_iterations,n_iterations)
	for i in xrange(0,n_iterations):
		P_min_p[i,:]   =fromfile(file=fd,dtype='d',count=n_P)
		P_avg_p[i,:]   =fromfile(file=fd,dtype='d',count=n_P)
		P_max_p[i,:]   =fromfile(file=fd,dtype='d',count=n_P)
		dP_avg_p[i,:]  =fromfile(file=fd,dtype='d',count=n_P)
		dP_sub_p[i,:]  =fromfile(file=fd,dtype='d',count=n_P)
		ln_Pr_min_p[i] =fromfile(file=fd,dtype='d',count=1)
		ln_Pr_avg_p[i] =fromfile(file=fd,dtype='d',count=1)
		ln_Pr_max_p[i] =fromfile(file=fd,dtype='d',count=1)
		P_min_c[i,:]   =fromfile(file=fd,dtype='d',count=n_P)
		P_avg_c[i,:]   =fromfile(file=fd,dtype='d',count=n_P)
		P_max_c[i,:]   =fromfile(file=fd,dtype='d',count=n_P)
		dP_avg_c[i,:]  =fromfile(file=fd,dtype='d',count=n_P)
		dP_sub_c[i,:]  =fromfile(file=fd,dtype='d',count=n_P)
		ln_Pr_min_c[i] =fromfile(file=fd,dtype='d',count=1)
		ln_Pr_avg_c[i] =fromfile(file=fd,dtype='d',count=1)
		ln_Pr_max_c[i] =fromfile(file=fd,dtype='d',count=1)
		slopes[i,:]    =fromfile(file=fd,dtype='d',count=n_P)
		drift[i,:]     =fromfile(file=fd,dtype='d',count=n_P)
	ln_Pr_min_p/=n_M_total
	ln_Pr_avg_p/=n_M_total
	ln_Pr_max_p/=n_M_total
	ln_Pr_min_c/=n_M_total
	ln_Pr_avg_c/=n_M_total
	ln_Pr_max_c/=n_M_total
	#ln_Pr_limit_min=0.9*min(ln_Pr_min_p[n_iterations_burn:])
	#ln_Pr_limit_max=1.1*min(ln_Pr_min_p[n_iterations_burn:])
	ln_Pr_limit_min=-10.00
	ln_Pr_limit_max=  0.05

	# ln_Pr trace
	plt.figure(figsize=(8,6))
	plt.clf()
	axTrace = plt.axes(rect_chi2)
	title(r'$\ln(Pr)$ trace plot for '+problem_name)
        min_line =linspace(ln_Pr_limit_min,ln_Pr_limit_min,n_iterations)
        max_line =linspace(ln_Pr_limit_max,ln_Pr_limit_max,n_iterations)
        axTrace.fill_between(x,ln_Pr_limit_min,ln_Pr_limit_max, where=x<=n_iterations_burn,facecolor='lightgrey',label=r'Burn Interval')
	axTrace.plot(x,ln_Pr_avg_c,c='black')
        axTrace.fill_between(x,ln_Pr_min_p,ln_Pr_max_p,facecolor='orange',edgecolor='none')
        axTrace.fill_between(x,ln_Pr_min_c,ln_Pr_max_c,facecolor='yellow',edgecolor='none')
        axTrace.set_xlim((1,n_iterations))
        axTrace.set_ylim((ln_Pr_limit_min,ln_Pr_limit_max))
	axTrace.set_xlabel("Averaging Interval")
        axTrace.set_ylabel(r'$\ln(Pr)$')
	#axHist  = plt.axes(rect_hist)
	#axHist.yaxis.set_major_formatter(nullfmt)
	#axHist.xaxis.set_major_formatter(nullfmt)
        #axHist.plot(hist_P,xx,c='black',linewidth=2)
        #axHist.fill_betweenx(xx,zero_line,hist_P, where=mask_95>0, facecolor=color_95pc,label=r'95% Confidence')
        #axHist.fill_betweenx(xx,zero_line,hist_P, where=mask_68>0, facecolor=color_68pc,label=r'68% Confidence')
        #axHist.plot([0.,1.],[best_val[0],best_val[0]],c='black',linewidth=2)
        #axHist.set_ylim((coverage_min[i],coverage_max[i]))
        #axHist.set_xlim((0,1.1))

	# Write plot
        filename_out = 'ln_Pr_trace.png'
        print '  Writing',filename_out
        filename_out =filename_plot_root+filename_out
        plt.savefig(filename_out,transparent=False)

	for i in xrange(0,n_P):
		# Create trace
                min_line =linspace(coverage_min[i],coverage_min[i],n_iterations)
                max_line =linspace(coverage_max[i],coverage_max[i],n_iterations)
		lo_line_1=P_avg_p[:,i]-dP_avg_p[:,i]
                hi_line_1=P_avg_p[:,i]+dP_avg_p[:,i]
                lo_line_2=P_avg_p[:,i]-dP_sub_p[:,i]
                hi_line_2=P_avg_p[:,i]+dP_sub_p[:,i]
		lo_line_3=P_avg_c[:,i]-dP_avg_c[:,i]
                hi_line_3=P_avg_c[:,i]+dP_avg_c[:,i]
                lo_line_4=P_avg_c[:,i]-dP_sub_c[:,i]
                hi_line_4=P_avg_c[:,i]+dP_sub_c[:,i]

		P_limit_min=0.9*min(P_min_p[n_iterations_burn:,i])
		P_limit_max=1.1*max(P_max_p[n_iterations_burn:,i])
		plt.figure(figsize=(8,6))
		plt.clf()
		axTrace = plt.axes(rect_trace)
		title('Trace plot for '+problem_name)
                axTrace.fill_between(x,min_line,max_line, where=x<=n_iterations_burn,facecolor='lightgrey',label=r'Burn Interval')
		axTrace.plot(x,P_avg_c[:,i],c='black')
                axTrace.fill_between(x,lo_line_1,hi_line_1,facecolor='orange',edgecolor='none')
                axTrace.fill_between(x,lo_line_2,hi_line_2,facecolor='yellow',edgecolor='none')
                axTrace.fill_between(x,lo_line_3,hi_line_3,facecolor='red',edgecolor='none')
                axTrace.fill_between(x,lo_line_4,hi_line_4,facecolor='blue',edgecolor='none')
                axTrace.set_xlim((1,n_iterations))
                axTrace.set_ylim((coverage_min[i],coverage_max[i]))
                axTrace.set_ylabel(P_name[i])
	  	axTrace.set_xlabel("Averaging Interval")

                # Create histogram
                best_val  =fromfile(file=fd2,dtype='d',count=1)
                lo_68pc   =fromfile(file=fd2,dtype='d',count=1)
                hi_68pc   =fromfile(file=fd2,dtype='d',count=1)
                lo_95pc   =fromfile(file=fd2,dtype='d',count=1)
                hi_95pc   =fromfile(file=fd2,dtype='d',count=1)
                hist_P    =fromfile(file=fd2,dtype='int64',count=coverage_size[0])
                xx        =linspace(coverage_min[i],coverage_max[i],coverage_size[0])
                hist_P    =double(hist_P)/double(hist_P.max())
                zero_line =linspace(0.,0.,coverage_size)
                mask_68   =zeros(coverage_size[0])
                mask_95   =zeros(coverage_size[0])
                for kk in xrange(0,coverage_size):
                        if(xx[kk]>=lo_68pc and xx[kk]<=hi_68pc):
                                mask_68[kk]=1
                        if(xx[kk]>=lo_95pc and xx[kk]<=hi_95pc):
                                mask_95[kk]=1
		axHist  = plt.axes(rect_hist)
		axHist.yaxis.set_major_formatter(nullfmt)
		axHist.xaxis.set_major_formatter(nullfmt)
                axHist.plot(hist_P,xx,c='black',linewidth=2)
                axHist.fill_betweenx(xx,zero_line,hist_P, where=mask_95>0, facecolor=color_95pc,label=r'95% Confidence')
                axHist.fill_betweenx(xx,zero_line,hist_P, where=mask_68>0, facecolor=color_68pc,label=r'68% Confidence')
                axHist.plot([0.,1.],[best_val[0],best_val[0]],c='black',linewidth=2)
                axHist.set_ylim((coverage_min[i],coverage_max[i]))
                axHist.set_xlim((0,1.1))

		# Write plot
                filename_out = 'chain_trace_000_'+str(i).zfill(5)+'.png'
                print '  Writing',filename_out
                filename_out =filename_plot_root+filename_out
                plt.savefig(filename_out,transparent=False)
	fd.close()
	fd2.close()
        print

# Finished
print 'Done'
sys.exit()

