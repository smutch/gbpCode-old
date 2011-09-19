from MCMC import MCMCrun, Chain

if __name__ == '__main__':
    test_run = MCMCrun('/home/ssi/smutch/Work/Models/MCMC_SAGE/gigglez_runs/mcmc_runs/gigglez_mr_01_MCMC')
    # Read in the fit parameter values
    test_run.read_param_results()

    chain0 = Chain(test_run, 0)


