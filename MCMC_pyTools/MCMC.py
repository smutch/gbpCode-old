#!/usr/bin/env python

import numpy as np
from utility import line_break, mkdir

class MCMCrun(object):

    """ MCMC run class."""

    def __init__(self, filename_root):

        """ Initialise the MCMC run.

        Input:
        filename_root (str)  -> the root directory of the gbpMCMC output

        """

        # These should not change
        self.__MCMC_NAME_SIZE__    =256
        self.__filename_plot_root__=filename_root+'/plots/'

        line_break()

        self.filename_root = filename_root

        # Read (and print) header info.
        fd                       =open(filename_root+'/run.dat','rb')
        self.problem_name        =np.fromfile(file=fd, dtype=(np.str_, self.__MCMC_NAME_SIZE__), count=1)[0].split('\x00')[0]
        self.n_chains            =np.fromfile(file=fd,dtype='i',count=1)[0]
        self.n_avg               =np.fromfile(file=fd,dtype='i',count=1)[0]
        self.flag_autocor_on_file=np.fromfile(file=fd,dtype='i',count=1)[0]
        self.flag_no_map_write   =np.fromfile(file=fd,dtype='i',count=1)[0]

        print 'Problem name       =',self.problem_name
        print 'n_chains           =',self.n_chains
        print 'n_avg              =',self.n_avg
        print 'flag_autocor       =',self.flag_autocor_on_file
        print 'flag_no_map_write  =',self.flag_no_map_write

        # Read (and print) parameter names.
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
            print '  ',self.P_name[i_P]
        print

        # Read (and print) parameter array names
        self.n_P_arrays=np.fromfile(file=fd,dtype='i',count=1)[0]
        if(self.n_P_arrays>0):
            print 'MCMC Parameter array(s):'
            self.P_array_name= []
            for i_array in xrange(0,self.n_P_arrays):
                self.P_array_name.append(np.fromfile(file=fd, dtype=(np.str_, self.__MCMC_NAME_SIZE__)
                    , count=1)[0].split('\x00')[0])
                junk=np.fromfile(file=fd,dtype='d',count=self.n_P)
                print '    ',self.P_array_name[i_array]
            print

        # Read (and print) the dataset names
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
            print '    ',self.DS_name[i_DS]
            if(self.n_DS_arrays[i_DS]>0):
                for i_array in xrange(0,self.n_DS_arrays[i_DS]):
                    self.DS_array_name.append(np.fromfile(file=fd, dtype=(np.str_, self.__MCMC_NAME_SIZE__), 
                        count=1)[0].split('\x00')[0])
                    junk=np.fromfile(file=fd,dtype='d',count=self.n_M[i_DS])
                    print '        array #'+str(i_array+1).zfill(2)+': '+self.DS_array_name[self.n_DS_arrays_total+i_array]
                self.DS_arrays_offset.append(self.n_DS_arrays_total)
            else:
                print "        No associated arrays."
            self.n_DS_arrays_total=self.n_DS_arrays_total+self.n_DS_arrays[i_DS]
        self.n_M_total=sum(self.n_M)
        print

        # Close the file
        fd.close()

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
                self.P_best[i_P]    = np.float(line[-1])
                self.P_upper95[i_P] = np.float(line[-2])
                self.P_lower95[i_P] = np.float(line[-3])
                self.P_upper68[i_P] = np.float(line[-4])
                self.P_lower68[i_P] = np.float(line[-5])
                self.P_max_l[i_P]   = np.float(line[-6])
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


    def read_best_fit_params(self):

        """ Read in and return the best fit parameters. """

        lP_best  = []

        for line in file(self.filename_root+'/results/fit_for_parameters.dat'):
            line = line.split()
            if(line[0][0]!='#'):
                lP_best.append(line[4])
            self.P_best =np.array(lP_best,dtype='d')


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
        else:
            if type(i_DS) == type(list()):
                for i in i_DS: print "%15s\t=\t%-1.3e" % (self.P_name[i], 10.0**(self.P_best[i]))
            else:
                assert type(i_DS) == type(int())
                print "%-25s\t=\t%-1.3e" % (self.P_name[i_DS], 10.0**(self.P_best[i_DS]))
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
        self.covariance_matrix.reshape((run.n_P, run.n_P))
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



