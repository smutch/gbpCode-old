#!/usr/bin/env python

import numpy as np
from utility import line_break

class MCMCrun:

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

        Returns:
            coverage_lo, coverage_hi, conf_68pc, conf_95pc, coverage_1, coverage_0, coverage_P
        """

        fd=open(self.filename_root+'/results/coverage.dat','rb')
        self.n_coverage   =np.fromfile(file=fd,dtype='i',count=1)[0]
        self.coverage_size=np.fromfile(file=fd,dtype='i',count=1)[0]
        self.coverage_min =np.fromfile(file=fd,dtype='d',count=self.n_P)
        self.coverage_max =np.fromfile(file=fd,dtype='d',count=self.n_P)
        coverage_size = self.coverage_size
        coverage_lo  =0
        coverage_hi  =10*self.n_integrate/(coverage_size*coverage_size)

        jj=0
        for jj in xrange(0,p_x):
            for kk in xrange(jj+1,self.n_P):
                conf_68pc =np.fromfile(file=fd,dtype='d',count=1)[0]
                conf_95pc =np.fromfile(file=fd,dtype='d',count=1)[0]
                coverage_1=np.fromfile(file=fd,dtype='int64',count=coverage_size*coverage_size)
                coverage_0=np.fromfile(file=fd,dtype='int64',count=coverage_size*coverage_size)
                coverage_P=np.fromfile(file=fd,dtype='int64',count=coverage_size*coverage_size)

        if (p_x>0): jj+=1
        for kk in xrange(jj+1,p_y+1):
            conf_68pc =np.fromfile(file=fd,dtype='d',count=1)[0]
            conf_95pc =np.fromfile(file=fd,dtype='d',count=1)[0]
            coverage_1=np.fromfile(file=fd,dtype='int64',count=coverage_size*coverage_size)
            coverage_0=np.fromfile(file=fd,dtype='int64',count=coverage_size*coverage_size)
            coverage_P=np.fromfile(file=fd,dtype='int64',count=coverage_size*coverage_size)

        if swap==False:
            coverage_1=np.transpose(np.reshape(coverage_1,[coverage_size,coverage_size]))
            coverage_0=np.transpose(np.reshape(coverage_0,[coverage_size,coverage_size]))
            coverage_P=np.transpose(np.reshape(coverage_P,[coverage_size,coverage_size]))
        else:
            coverage_1=np.reshape(coverage_1,[coverage_size,coverage_size])
            coverage_0=np.reshape(coverage_0,[coverage_size,coverage_size])
            coverage_P=np.reshape(coverage_P,[coverage_size,coverage_size])

        fd.close()

        return coverage_lo, coverage_hi, conf_68pc, conf_95pc, coverage_1, coverage_0, coverage_P



class Chain:

    """Class representing an MCMC chain."""

    def __init__(self, run, my_chain):

        print
        line_break()
        print "Chain %03d:" %(my_chain)

        self.my_chain = my_chain

        # Read (and print) the number of iterations
        fd                         =open(run.filename_root+'/chains/chain_config_'+str(my_chain).zfill(6)+'.dat','rb')
        self.n_iterations          =np.fromfile(file=fd,dtype='i',count=1)[0]
        self.n_iterations_burn     =np.fromfile(file=fd,dtype='i',count=1)[0]
        self.n_iterations_integrate=self.n_iterations-self.n_iterations_burn
        self.n_burn                =run.n_avg*self.n_iterations_burn
        self.n_integrate           =run.n_avg*self.n_iterations_integrate
        self.n_total               =self.n_burn+self.n_integrate
        fd.close()  # close the file
        print 'n_iterations           =',self.n_iterations
        print 'n_iterations_burn      =',self.n_iterations_burn
        print 'n_iterations_integrate =',self.n_iterations_integrate
        print 'n_burn                 =',self.n_burn
        print 'n_integrate            =',self.n_integrate
        print 'n_total                =',self.n_total

        line_break()


