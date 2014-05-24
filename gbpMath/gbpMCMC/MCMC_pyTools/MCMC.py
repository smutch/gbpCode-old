#!/usr/bin/env python

import numpy as np
from utility import line_break, mkdir
import random

class MCMCrun(object):

    """ MCMC run class."""

    def __init__(self, filename_root, quiet=False):

        """ Initialise the MCMC run.

        Input:
        filename_root (str)  -> the root directory of the gbpMCMC output

        """

        # These should not change
        self.__MCMC_NAME_SIZE__    =256
        self.__filename_plot_root__=filename_root+'/plots/'

        self.filename_root = filename_root

        # Read (and print) header info.
        fd                       =open(filename_root+'/run.dat','rb')
        self.problem_name        =np.fromfile(file=fd, dtype=(np.str_, self.__MCMC_NAME_SIZE__), count=1)[0].split('\x00')[0]
        self.n_chains            =np.fromfile(file=fd,dtype='i',count=1)[0]
        self.n_avg               =np.fromfile(file=fd,dtype='i',count=1)[0]
        self.flag_autocor_on_file=np.fromfile(file=fd,dtype='i',count=1)[0]
        self.flag_no_map_write   =np.fromfile(file=fd,dtype='i',count=1)[0]

        if not quiet:
            line_break()
            print 'Problem name       =',self.problem_name
            print 'n_chains           =',self.n_chains
            print 'n_avg              =',self.n_avg
            print 'flag_autocor       =',self.flag_autocor_on_file
            print 'flag_no_map_write  =',self.flag_no_map_write

        # Read (and print) parameter names.
        if not quiet:
            print '\nMCMC Parameters:'
        self.n_P=np.fromfile(file=fd,dtype='i',count=1)[0]
        self.P_name= []
        self.P_init= []
        self.P_limit_min= []
        self.P_limit_max= []
        for i_P in xrange(0,self.n_P):
            self.P_name.append(np.fromfile(file=fd, dtype=(np.str_, self.__MCMC_NAME_SIZE__), count=1)[0].split('\x00')[0])
            self.P_init.append(np.fromfile(file=fd,dtype='d',count=1)[0])
            self.P_limit_min.append(np.fromfile(file=fd,dtype='d',count=1)[0])
            self.P_limit_max.append(np.fromfile(file=fd,dtype='d',count=1)[0])
            if not quiet:
                print '  ',self.P_name[i_P]
        
        if not quiet:
            print

        # Read (and print) parameter array names
        self.n_P_arrays=np.fromfile(file=fd,dtype='i',count=1)[0]
        if(self.n_P_arrays>0):
            if not quiet:
                print 'MCMC Parameter array(s):'
            self.P_array_name= []
            for i_array in xrange(0,self.n_P_arrays):
                self.P_array_name.append(np.fromfile(file=fd, dtype=(np.str_, self.__MCMC_NAME_SIZE__)
                    , count=1)[0].split('\x00')[0])
                np.fromfile(file=fd,dtype='d',count=self.n_P)
                if not quiet:
                    print '    ',self.P_array_name[i_array]
            if not quiet:
                print

        # Read (and print) the dataset names
        if not quiet:
            print 'MCMC Datasets:'
        self.n_DS  =np.fromfile(file=fd,dtype='i',count=1)[0]
        self.n_M      =[]
        self.DS_name  =[]
        self.M_target =[]
        self.dM_target=[]
        self.DS_array_name    = []
        self.n_DS_arrays      = []
        self.DS_arrays_offset = []
        self.n_DS_arrays_total= 0
        for i_DS in xrange(0,self.n_DS):
            self.DS_name.append(np.fromfile(file=fd, dtype=(np.str_, self.__MCMC_NAME_SIZE__), count=1)[0].split('\x00')[0])
            self.n_M.append(np.fromfile(file=fd,dtype='i',count=1)[0])
            self.M_target.append(np.fromfile(file=fd,dtype='d',count=self.n_M[i_DS]))
            self.dM_target.append(np.fromfile(file=fd,dtype='d',count=self.n_M[i_DS]))
            self.n_DS_arrays.append(np.fromfile(file=fd,dtype='i',count=1)[0])
            if not quiet:
                print '    ',self.DS_name[i_DS]
            if(self.n_DS_arrays[i_DS]>0):
                for i_array in xrange(0,self.n_DS_arrays[i_DS]):
                    self.DS_array_name.append(np.fromfile(file=fd, dtype=(np.str_, self.__MCMC_NAME_SIZE__), 
                        count=1)[0].split('\x00')[0])
                    np.fromfile(file=fd,dtype='d',count=self.n_M[i_DS])
                    if not quiet:
                        print '        array #'+str(i_array+1).zfill(2)+': '+self.DS_array_name[self.n_DS_arrays_total+i_array]
                self.DS_arrays_offset.append(self.n_DS_arrays_total)
            else:
                if not quiet:
                    print "        No associated arrays."
            self.n_DS_arrays_total=self.n_DS_arrays_total+self.n_DS_arrays[i_DS]
        self.n_M_total=sum(self.n_M)
        if not quiet:
            print

        # Close the file
        fd.close()

        if not quiet:
            line_break()


    def read_param_results(self, quiet=False):

        """ Read in and return the best fit parameters. 
        
        Input:
            quiet - print the param results file? (true/false)
        """

        self.P_best    = np.zeros(self.n_P, float)
        self.P_upper95 = np.zeros(self.n_P, float)
        self.P_lower95 = np.zeros(self.n_P, float)
        self.P_upper68 = np.zeros(self.n_P, float)
        self.P_lower68 = np.zeros(self.n_P, float)
        self.P_max_l   = np.zeros(self.n_P, float)
        self.P_std     = np.zeros(self.n_P, float)
        self.P_av      = np.zeros(self.n_P, float)
        self.P_initial = np.zeros(self.n_P, float)

        if quiet==False:
            line_break()
        
        i_P = 0
        for line in file(self.filename_root+'/results/fit_for_parameters.dat'):
            if quiet==False:
                print line.rstrip('\n')
            line = line.split()
            if(line[0][0]!='#'):
                self.P_max_l[i_P]   = np.float(line[-1])
                self.P_upper95[i_P] = np.float(line[-2])
                self.P_lower95[i_P] = np.float(line[-3])
                self.P_upper68[i_P] = np.float(line[-4])
                self.P_lower68[i_P] = np.float(line[-5])
                self.P_best[i_P]    = np.float(line[-6])
                self.P_std[i_P]     = np.float(line[-7])
                self.P_av[i_P]      = np.float(line[-8])
                self.P_initial[i_P] = np.float(line[-9])
                i_P+=1

        if quiet==False:
            line_break()


    def read_coverage_info(self):

        """Read in the coverage info."""

        # Read coverage information
        fd=open(self.filename_root+'/results/coverage.dat','rb')
        self.n_coverage   =np.fromfile(file=fd,dtype='i',count=1)[0]
        self.coverage_size=np.fromfile(file=fd,dtype='i',count=1)[0]
        self.coverage_min =np.fromfile(file=fd,dtype='d',count=self.n_P)
        self.coverage_max =np.fromfile(file=fd,dtype='d',count=self.n_P)
        fd.close()      # close the file



    def read_coverage_map(self, p_x, p_y, swap=False):

        """Read the coverage data of a given dataset. 
        
        Input:
            p_x  - Index of the x-axis parameter 
            p_y  - Index of the y-axis parameter
            swap - Swap the axis round (true/false)

        Returns a dict with the following entries:
            conf_68pc, conf_95pc, coverage_1, coverage_0, coverage_P, 
            max_L, mean_L
        """

        fd=open(self.filename_root+'/results/coverage.dat','rb')
        self.n_coverage   =np.fromfile(file=fd,dtype='i',count=1)[0]
        self.coverage_size=np.fromfile(file=fd,dtype='i',count=1)[0]
        self.coverage_min =np.fromfile(file=fd,dtype='d',count=self.n_P)
        self.coverage_max =np.fromfile(file=fd,dtype='d',count=self.n_P)
        coverage_size = self.coverage_size

        jj=0
        for jj in xrange(0,p_x):
            for kk in xrange(jj+1,self.n_P):
                conf_68pc =np.fromfile(file=fd,dtype='d',count=1)[0]
                conf_95pc =np.fromfile(file=fd,dtype='d',count=1)[0]
                coverage_1=np.fromfile(file=fd,dtype='int64',count=coverage_size*coverage_size)
                coverage_0=np.fromfile(file=fd,dtype='int64',count=coverage_size*coverage_size)
                coverage_P=np.fromfile(file=fd,dtype='int64',count=coverage_size*coverage_size)
                max_L     =np.fromfile(file=fd,dtype='d',count=coverage_size*coverage_size)
                mean_L    =np.fromfile(file=fd,dtype='d',count=coverage_size*coverage_size)

        if (p_x>0): jj+=1
        for kk in xrange(jj+1,p_y+1):
            conf_68pc =np.fromfile(file=fd,dtype='d',count=1)[0]
            conf_95pc =np.fromfile(file=fd,dtype='d',count=1)[0]
            coverage_1=np.fromfile(file=fd,dtype='int64',count=coverage_size*coverage_size)
            coverage_0=np.fromfile(file=fd,dtype='int64',count=coverage_size*coverage_size)
            coverage_P=np.fromfile(file=fd,dtype='int64',count=coverage_size*coverage_size)
            max_L     =np.fromfile(file=fd,dtype='d',count=coverage_size*coverage_size)
            mean_L    =np.fromfile(file=fd,dtype='d',count=coverage_size*coverage_size)

        max_L_norm =np.double(max_L.max())
        mean_L_norm=np.double(mean_L.max())

        if swap==False:
            coverage_1=np.transpose(np.reshape(coverage_1,[coverage_size,coverage_size]))
            coverage_0=np.transpose(np.reshape(coverage_0,[coverage_size,coverage_size]))
            coverage_P=np.transpose(np.reshape(coverage_P,[coverage_size,coverage_size]))
            max_L      =10**(np.transpose(np.reshape(max_L,[coverage_size,coverage_size]))-max_L_norm)
            mean_L     =10**(np.transpose(np.reshape(mean_L,[coverage_size,coverage_size]))-mean_L_norm)

        else:
            coverage_1=np.reshape(coverage_1,[coverage_size,coverage_size])
            coverage_0=np.reshape(coverage_0,[coverage_size,coverage_size])
            coverage_P=np.reshape(coverage_P,[coverage_size,coverage_size])
            max_L      =10**(np.reshape(max_L,[coverage_size,coverage_size])-max_L_norm)
            mean_L     =10**(np.reshape(mean_L,[coverage_size,coverage_size])-mean_L_norm)

        fd.close()

        coverage_dict = {
                'conf_68pc':conf_68pc,
                'conf_95pc':conf_95pc,
                'coverage_1':coverage_1,
                'coverage_0':coverage_0,
                'coverage_P':coverage_P,
                'max_L':max_L,
                'mean_L':mean_L
                }

        return coverage_dict


    def read_dataset_fit(self, i_DS=0):

        """Read the fit for a given dataset. 
        
        Input:
            i_DS - Index of the requested dataset 

        Returns a dict with the following entries:
           x, M_DS, dM_DS, M_best, M_lo_68, M_hi_68, M_lo_95, M_hi_95,
           dM_hi_DS, dM_lo_DS
        """


        if(self.n_DS_arrays[i_DS]>0):
            lx       = []
            lM_DS    = []
            ldM_DS   = []
            lM_best  = []
            lM_lo_68 = []
            lM_hi_68 = []
            lM_lo_95 = []
            lM_hi_95 = []
            if(self.n_DS>1):
                filename_read=self.filename_root+'/results/fit_for_dataset_'+str(i_DS).zfill(5)+'.dat'
            else:
                filename_read=self.filename_root+'/results/fit_for_dataset.dat'

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

            x      =np.array(lx,dtype=np.float64)
            M_DS   =np.array(lM_DS,dtype=np.float64)
            dM_DS  =np.array(ldM_DS,dtype=np.float64)
            M_best =np.array(lM_best,dtype=np.float64)
            M_lo_68=np.array(lM_lo_68,dtype=np.float64)
            M_hi_68=np.array(lM_hi_68,dtype=np.float64)
            M_lo_95=np.array(lM_lo_95,dtype=np.float64)
            M_hi_95=np.array(lM_hi_95,dtype=np.float64)
            dM_hi_DS=dM_DS
            dM_lo_DS=dM_DS
            for i_M in xrange(0,self.n_M[i_DS]):
                if((M_DS[i_M]-dM_lo_DS[i_M])<0):
                    dM_lo_DS[i_M]=0.99999*M_DS[i_M]

            return {'x'         : x,
                    'M_DS'      : M_DS,
                    'dM_DS'     : dM_DS,
                    'M_best'    : M_best,
                    'M_lo_68'   : M_lo_68,
                    'M_hi_68'   : M_hi_68,
                    'M_lo_95'   : M_lo_95,
                    'M_hi_95'   : M_hi_95,
                    'dM_hi_DS'  : dM_hi_DS,
                    'dM_lo_DS'  : dM_lo_DS}

        else:
            raise ValueError("There are no arrays associated with this\
                             dataset.")



    def read_best_fit_params(self):

        """ Read in and return the best fit parameters. """

        lP_best  = []
        lP_max_l = []

        for line in file(self.filename_root+'/results/fit_for_parameters.dat'):
            line = line.split()
            if(line[0][0]!='#'):
                lP_best.append(line[4])
                lP_max_l.append(line[9])
            self.P_best =np.array(lP_best,dtype='d')
            self.P_max_l =np.array(lP_max_l,dtype='d')

        return self.P_best, self.P_max_l


    def read_histogram(self, i_DS):

        """ Read in the histogram for a single parameter. 
        
        Args:
            i_DS   :    The index of the parameter.
        
        Returns:
            x :         Parameters values
            hist_P :    Histogram values
            hist_mean : ??
            hist_max :  ??
            mask_68 :   ??
            mask_95 :   ??
            best_val :  ??
        """

        fd=open(self.filename_root+'/results/histograms.dat','rb')
        for jj in xrange(0,i_DS+1):
            best_val  = np.fromfile(file=fd,dtype='d',count=1)[0]
            lo_68pc   = np.fromfile(file=fd,dtype='d',count=1)[0]
            hi_68pc   = np.fromfile(file=fd,dtype='d',count=1)[0]
            lo_95pc   = np.fromfile(file=fd,dtype='d',count=1)[0]
            hi_95pc   = np.fromfile(file=fd,dtype='d',count=1)[0]
            peak_val  = np.fromfile(file=fd,dtype='d',count=1)[0]
            hist_P    = np.fromfile(file=fd,dtype='uint64',count=self.coverage_size)
            hist_mean = np.fromfile(file=fd,dtype='d',count=self.coverage_size)
            hist_max  = np.fromfile(file=fd,dtype='d',count=self.coverage_size)

        x         =np.linspace(self.coverage_min[jj],self.coverage_max[jj],self.coverage_size)
        mask_68   =np.zeros(self.coverage_size)
        mask_95   =np.zeros(self.coverage_size)
        for kk in xrange(0,self.coverage_size):
           if(x[kk]>=lo_68pc and x[kk]<=hi_68pc):
              mask_68[kk]=1
           if(x[kk]>=lo_95pc and x[kk]<=hi_95pc):
              mask_95[kk]=1
        hist_P        =np.double(hist_P)/np.double(hist_P.max())
        hist_max_norm =np.double(hist_max.max())
        hist_mean_norm=np.double(hist_mean.max())
        hist_max      =10**(np.double(hist_max)-hist_max_norm)
        hist_mean     =10**(np.double(hist_mean)-hist_mean_norm)

        fd.close()

        return x, hist_P, hist_mean, hist_max, mask_68, mask_95, best_val




    def print_most_likely(self, i_DS=None):

        """ Print the most likely values for 1 or all parameters. """

        print
        print "-----------------------------"
        print "Most likely parameter values:"
        print "-----------------------------"

        if i_DS == None:
            for i, best in enumerate(self.P_best):
                print "%-25s\t=\t%-1.3e" % (self.P_name[i], 10.0**best)
                print "\t\t68%%: -%.3e  +%.3e" % (10.**best - 10.**self.P_lower68[i],
                                             10.**self.P_upper68[i] - 10.**best)
                print "\t\t95%%: -%.3e  +%.3e" % (10.**best - 10.**self.P_lower95[i],
                                             10.**self.P_upper95[i] - 10.**best)
        else:
            if type(i_DS) == type(list()):
                for i in i_DS: 
                    best = self.P_best[i]
                    print "%-25s\t=\t%-1.3e" % (self.P_name[i], 10.0**(self.P_max_l[i]))
                    print "\t\t68%%: -%.3e  +%.3e" % (10.**best - 10.**self.P_lower68[i],
                                                 10.**self.P_upper68[i] - 10.**best)
                    print "\t\t95%%: -%.3e  +%.3e" % (10.**best - 10.**self.P_lower95[i],
                                                 10.**self.P_upper95[i] - 10.**best)
            else:
                assert type(i_DS) == type(int())
                best = self.P_best[i_DS]
                print "%-25s\t=\t%-1.3e" % (self.P_name[i_DS], 10.0**(self.P_max_l[i_DS]))
                print "\t\t68%%: -%.3e  +%.3e" % (10.**best - 10.**self.P_lower68[i_DS],
                                             10.**self.P_upper68[i_DS] - 10.**best)
                print "\t\t95%%: -%.3e  +%.3e" % (10.**best - 10.**self.P_lower95[i_DS],
                                             10.**self.P_upper95[i_DS] - 10.**best)
        print "------------------------"
        print


class Chain(object):

    """Class representing an MCMC chain."""

    def __init__(self, run, i_chain, quiet=False):

        if quiet == False:
            print
            line_break()
            print "Chain %03d:" %(i_chain)
            print "----------"

        self.i_chain = i_chain
        self.run = run

        # Read (and print) the number of iterations
        fd                         = open(run.filename_root+'/chains/chain_config_'+str(i_chain).zfill(6)+'.dat','rb')
        self.n_iterations          = np.fromfile(file=fd,dtype='i',count=1)[0]
        self.n_iterations_burn     = np.fromfile(file=fd,dtype='i',count=1)[0]
        self.temp                  = np.fromfile(file=fd,dtype=np.float64,count=1)[0]
        self.n_iterations_integrate= self.n_iterations-self.n_iterations_burn
        self.n_burn                = run.n_avg*self.n_iterations_burn
        self.n_integrate           = run.n_avg*self.n_iterations_integrate
        self.n_total               = self.n_burn+self.n_integrate
        self.covariance_matrix     = np.fromfile(file=fd,dtype=np.float64, count=run.n_P*run.n_P)
        self.covariance_matrix     = self.covariance_matrix.reshape((run.n_P, run.n_P))
        fd.close()  # close the file
        if quiet == False:
            print 'n_iterations           =',self.n_iterations
            print 'n_iterations_burn      =',self.n_iterations_burn
            print 'n_iterations_integrate =',self.n_iterations_integrate
            print 'n_burn                 =',self.n_burn
            print 'n_integrate            =',self.n_integrate
            print 'n_total                =',self.n_total
            print 'temperature            =',self.temp
            line_break()

        def __repr__(object):
            return "Chain"

    def read_trace(self, phase):
        """Read in the chain trace parameter values.
        
        Args:
            phase   :   'integration' OR 'burnin'
        
        Returns:
            flag_success   : char value (!=0 if proposal was successful)
            ln_likelihood  : natural log of likelihood value of propositions
            param_vals     : parameter values for each proposition

        TODO - Read in the dataset results as well if
               run.flag_autocor_on_file==0
        """
        
        run = self.run
        fin = open(run.filename_root+'/chains/chain_trace_{:06d}'.format(self.i_chain)+'.dat', 'rb')

        if(phase)=='integration':
            # Seek past the burn
            byte_seek = 1+8+(8*run.n_P)
            if run.flag_no_map_write==0:
                for i_DS in xrange(run.n_DS):
                    byte_seek += 8*run.n_M[i_DS]
            fin.seek(byte_seek*self.n_burn,1)
            n_links = self.n_integrate
        elif(phase)=='burnin':
            n_links = self.n_burn
        
        # Allocate the storage arrays
        n_P = run.n_P
        flag_success = np.zeros(n_links, dtype='S1')
        ln_likelihood = np.zeros(n_links, np.double)
        p_vals_dtype = np.dtype((np.float64, n_P))
        param_vals = np.zeros(n_links, dtype=p_vals_dtype)
        
        # Read the values
        for l in xrange(n_links):
            flag_success[l] = np.fromfile(fin, dtype='S1', count=1)[0]
            ln_likelihood[l] = np.fromfile(fin, dtype=np.float64, count=1)[0]
            param_vals[l] = np.fromfile(fin, dtype=p_vals_dtype, count=1)[0]
            if run.flag_no_map_write==0:
                for i_DS in xrange(run.n_DS):
                    fin.seek(8*run.n_M[i_DS], 1)

        fin.close()
        return flag_success, ln_likelihood, param_vals 

    def _chain_(self,success,props):
        """Construct the chain using the trace."""
        pcur = self.run.P_init
        for i in xrange(len(success)):
            if success[i]=='':
                props[i] = pcur
            else:
                pcur = props[i]
        return props


    def PCA(self, quiet=False):
        """Carry out a Principle Component Analysis (PCA) of the successful
        integration phase propositions.

        Args:
            verbose   :  If True then print the energy fraction and eigenvector
                         summary for the principal components containing >=90%
                         of the total energy.

        Returns:
            eigenval  :  The eigenvalues of the PCA
            cumenergy :  The cumulative energy fraction of the components
            eigenvec  :  The eigenvectors (principle components)
        """

        # Read in the trace
        success, ln_likelihood, props = self.read_trace('integration')

        # Construct the chain using the trace
        props = self._chain_(success, props)

        # These now do not match props and should not be used...
        del(success, ln_likelihood)
        
        # subtract the mean (along columns)
        M = (props-np.mean(props.T,axis=1)).T 
        # WARNING - ^!OR!v
        # subtract the marginalised best parameters along the columns
        # self.run.read_param_results()
        # M = (props-self.run.P_best).T 
        
        # Calculate the eigenvectors and eigenvalues
        eigenval,eigenvec = np.linalg.eig(np.cov(M))

        # Sort the eigenvectors and values according to the eigenvalues
        w = np.argsort(eigenval)[::-1]
        eigenval = eigenval[w]
        eigenvec = eigenvec[w]

        # Calculate the cumulative energy fraction
        eigenval_sum = eigenval.sum() 
        cumenergy = eigenval.cumsum()/eigenval_sum

        if not quiet:
            line_break()
            print "PCA\n----------"
            pc_90 = np.where(cumenergy>=0.9)[0][0]
            strlen = np.array([ len(self.run.P_name[i]) for i in xrange(self.run.n_P) ]
                              , int).max() + 3
            for pc in xrange(eigenval.size):
                print
                print "Principle component #{:d}:".format(pc)
                print "Energy fraction = {:.3f} (running total = {:.3f})"\
                        .format(eigenval[pc]/eigenval_sum, \
                                eigenval[:pc+1].sum()/eigenval_sum)
                for p in np.argsort(np.abs(eigenvec[pc]))[::-1]:
                    print "{:s} : {:-.3f}  (^2 = {:-.3f})".format(
                        self.run.P_name[p].ljust(strlen),
                        eigenvec[pc][p],
                        eigenvec[pc][p]**2.)
                if pc == pc_90:
                    print
                    print "^^^^^ 90% ^^^^^"
            print


        return eigenval, cumenergy, eigenvec


    def sample_chain(self, n_samples, seed=None):
        """Construct the chain from the trace and then randomly sample from it.

        Args:
            n_samples   :  The number of samples to draw
            seed        :  Optional seed for random number generator

        Returns:
            samples     :  The sampled chain
        """

        # Read in the trace
        success, ln_likelihood, props = self.read_trace('integration')

        # Construct the chain using the trace
        props = self._chain_(success, props)

        # These now do not match props and should not be used...
        del(success, ln_likelihood)

        # Set up the RNG
        random.seed(seed) 

        # Sample the chain and return the values
        return random.sample(props, n_samples)



def check_param_compatibility(run_list, param):
    
    param1 = run_list[0].__getattribute__(param)
    if any([run.__getattribute__(param) !=param1 for run in run_list[1:]]):
        raise ValueError('Run parameters `%s` do not match!...')
    else:
        return True


def join_runs(run_list, joined_fname_root):

    try:
        shutil.__name__
        deque.__name__
    except NameError:
        import shutil
        from collections import deque

    # Create the directory structure
    mkdir(joined_fname_root+'/')
    mkdir(joined_fname_root+'/chains/')
    mkdir(joined_fname_root+'/results/')
    mkdir(joined_fname_root+'/plots/')

    # Do some checks to ensure that the runs we are joining are compatible with
    # each other...
    for param in ['problem_name', 'n_avg', 'flag_autocor_on_file',
            'flag_no_map_write', 'n_P', 'P_name', 'P_limit_min', 
            'P_limit_max', 'n_DS_arrays_total', 'n_DS']:
        check_param_compatibility(run_list, param)
    
    # Copy the run.dat file for the first run as our new run.dat file
    shutil.copy(run_list[0].filename_root+'/run.dat',
            joined_fname_root)


    # Read in the chain config files and work out the values for our new file
    n_iterations = 0
    n_iterations_burn = None
    chain_list = deque()
    for run in run_list:
        for i_chain in xrange(run.n_chains):
            chain = Chain(run, i_chain)
            chain_list.append(chain)
            if n_iterations_burn==None: 
                n_iterations_burn=chain.n_iterations_burn
                n_iterations += chain.n_iterations
            else:
                n_iterations += chain.n_iterations-chain.n_iterations_burn
    covariance_matrix = chain_list[0].covariance_matrix
    temp = chain_list[0].temp

    # Write the new chain config file 
    fout = open(joined_fname_root+'/chains/chain_config_%06d.dat'%(0), 'wb')
    np.array([n_iterations], dtype=np.int32).tofile(fout)
    np.array([n_iterations_burn], dtype=np.int32).tofile(fout)
    np.array([temp], dtype=np.float64).tofile(fout)
    covariance_matrix.flatten().tofile(fout)
    fout.close()


    # Concatenate all of the chain stats files together
    fout = open(joined_fname_root+'/chains/chain_stats_%06d.dat'%(0), 'wb')
    first_burn = True
    for i_run, run in enumerate(run_list):
        for i_chain in xrange(run.n_chains):
            fin = open(run.filename_root+'/chains/chain_stats_%06d.dat'%(i_chain), 'rb')
            if first_burn==False:
                # Seek past the burn of this chain
                n_iterations_burn = chain_list[i_run+i_chain].n_iterations_burn
                byte_seek = (8*run.n_P*12)+(8*6)
                if run.flag_autocor_on_file==1:
                    byte_seek+=(run.n_avg-1)*8
                fin.seek(byte_seek*n_iterations_burn,1)
            else:
                first_burn = False
            shutil.copyfileobj(fin, fout)
            fin.close()
    fout.close()


    # Concatenate all of the chain trace files together
    fout = open(joined_fname_root+'/chains/chain_trace_%06d.dat'%(0), 'wb')
    first_burn = True
    for run in run_list:
        for i_chain in xrange(run.n_chains):
            fin = open(run.filename_root+'/chains/chain_trace_%06d.dat'%(i_chain), 'rb')
            if first_burn == False:
                # Seek past the burn of this chain
                chain = chain_list[i_run+i_chain]
                n_iterations_burn = chain.n_iterations_burn
                n_avg = run.n_avg
                byte_seek = 1+8+(8*run.n_P)
                if run.flag_no_map_write==0:
                    for i_DS in xrange(run.n_DS):
                        byte_seek += 8*run.n_M[i_DS]
                fin.seek(byte_seek*n_iterations_burn*n_avg,1)
            else:
                first_burn = False
            shutil.copyfileobj(fin, fout)
            fin.close()
    fout.close()


def restore_failed_chain(run, i_chain, n_iterations_burn, n_iterations, n_avg):
    '''Try and restore a chain file to the last completed integration
    iteration.  This is used to recover from a crashed run.'''

    import shutil
    import os

    print "i_chain = %d" % i_chain
    print "n_itertions_burn = %d" % n_iterations_burn
    print "n_iterations = %d" % n_iterations
    print "n_avg = %d" % n_avg
    
    # Read the chain trace and see how far we got before we failed
    fname = run.filename_root+'/chains/chain_trace_%06d.dat'%(i_chain)

    # Work out the number of bytes in the burnin
    n_avg = run.n_avg
    byte_seek = 1+8+(8*run.n_P)
    if run.flag_no_map_write==0:
        for i_DS in xrange(run.n_DS):
            byte_seek += 8*run.n_M[i_DS]
    burn_bytes = long(byte_seek*n_iterations_burn*n_avg)

    # Now work out the number of successfully completed interations in the
    # chain trace.
    bytes_left = os.path.getsize(fname) - burn_bytes 
    print "burn_bytes = %ld" % burn_bytes
    print "bytes_left = %ld" % bytes_left
    print "file_size = %ld" % os.path.getsize(fname)
    print "byte_seek = %ld" % byte_seek
    completed_iterations = int(bytes_left / (byte_seek*n_avg))
    bytes_requested = long(completed_iterations*n_iterations*n_avg) + burn_bytes
    print "Found %d successfully completed integration phase iterations..." % completed_iterations
    print "Creating new file with %ld bytes..." % bytes_requested

    # Copy the old chain trace file for backup
    fname_bu = fname+'.backup'
    shutil.move(fname, fname_bu)

    # Write the new trace file
    fout = open(fname, 'wb')
    fin = open(fname_bu, 'rb')
    fout.write(fin.read(bytes_requested))
    fin.close()
    fout.close()

    print "Created new chain file:",fname
    print "Original chain file copied to:",fname_bu
    print


class ConvergenceTests(object):
    """A collection of convergence tests for MCMC chains.

    Currently implemented:
        Rubin-Gelman statistic (Requires >= 2 chains)
    """

    def __init__(self,chains):
        # Make sure that chains is a list
        if type(chains)!=list:
            self.chains = [chains,]
        else:
            self.chains = chains
        assert type(self.chains[0]) == Chain,\
                "Valid chain object not passed..."
        self.nchains = len(chains)
    
    def rubin_gelman(self):
        """Rubin-Gelman statistic for chain convergence.

        This test requires >=2 chains.

        See
        http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_introbayes_sect008.htm
        for details on the implementation.
        
        Returns:
            PSRF - The potential square reduction factor
                   (target = 1.0 for converged chains)
        """

        if self.nchains<2:
            raise ValueError('nchains is less than 2')   

        nchains = float(self.nchains) 
        chain_list = self.chains
        nparams =  chain_list[0].run.n_P
        param_names = chain_list[0].run.P_name
        chain_mean = np.zeros((nchains, nparams), float)
        chain_var  = np.zeros((nchains, nparams), float)
        
        # We need the same length for each chain hence we are restricted to the
        # length of the shortest chain.
        lengths = np.array([chain.n_integrate for chain in chain_list], float)
        if any(lengths!=lengths[0]):
            length = lengths.min()
        else:
            length = lengths[0]

        for i_chain,chain in enumerate(chain_list):
            success, ln_likelihood, props = chain.read_trace('integration')
            props = chain._chain_(success, props)  # Reconstruct chain from trace
            del(success, ln_likelihood)
            if props.shape[0]>length:
                props = props[:length]  # Trim the chain if necessary
            chain_mean[i_chain] = props.mean(axis=0)
            chain_var[i_chain]  = props.var(ddof=1, axis=0)

        print
        print "Rubin-Gelman Convergence Test:"
        print "------------------------------"
        print "Potential square reduction factors (target = 1.0):"

        PSRF = np.zeros(nparams, float)
        for i_param, param in enumerate(param_names):
            B = chain_mean[:,i_param].var(ddof=1)*length  # inter-chain var
            W = chain_var[:,i_param].mean()  # intra-chain var

            # Posterior marginal variance
            V = ((length-1)/length*W) + ((nchains+1)/(length*nchains)*B)

            # potential square reduction factor
            PSRF[i_param] = np.sqrt(V/W)

            # print result
            print "\t%-30s -> %.2f"%(param,PSRF[i_param])

        return PSRF

