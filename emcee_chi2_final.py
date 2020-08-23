#Created on Thu Aug 13 13:46:51 2020
#
#####################  emcee_chi2_final.py  #####################
#
#
### The following script runs the "emcee" MCMC sampler, using a chi-squared metric as the cost function: chi^2 = -0.5 * ( (data - model)^2 / error^2 ) to fit a SCHECHTER function to scatter plot
#
#
### Section summary:
#
### PROGRAM START
#
### (0)    import modules, define functions; FLAGS
### (1)    setup lists and DEFINE #DIM, #WALKERS, #STEPS
### (2)    initialize MCMC guess w/ LEAST-SQUARES estimate
### (3)    MCMC SIMULATION
### (4)    Visualize
### (5)    print to OUTPUT
#
### PROGRAM END
#
#
#
#
## SECTION (0): import MODULES, DEFINE functions, FLAGS
#
## MODULES
import pandas as pd
import emcee
import time
import scipy.optimize as op
import scipy.integrate as integrate
import corner
#
#
#
# FLAGS: 0 = off (i.e. skip), 1 = on (i.e. execute code)
#
#
#
SF_flag = 0          # star-forming population
Q_flag = 1           # Q pop
T_flag = 0           # total pop
#
#
#
######
###### MAY NEED TO EDIT ###### 
######
#
# initialize walkers
ndim, nwalkers, nsteps = 3, 1000, 100000    # (# of parameters), (#walkers), (#steps/walker)
#
#
## SECTION (1): setup lists
#
### SETUP appropriate lists and convert them to arrays for use in MCMC code:
## this first bit replaces 'NaN' entries in SF_error (errors to the SMF) with 1e5. This only applies to pixels for which the SMF count = 0, hence the error (given by sqrt(smf) is NaN. It is done for the purpose of computing chi_squared later on
#for ii in range(len(SF_error)):
#    if np.isnan(SF_error[ii]) == 1:
#        SF_error[ii] = 1e5
#
## this first bit drops all 'NaN' entries in SF_error, shortening the length of the array
#SF_midbins_mcmc = []
#SF_smf_mcmc = []
#SF_error_mcmc = []
#for ii in range(len(SF_error)):
#    if np.isnan(SF_error[ii]) == 0:
#        SF_midbins_mcmc.append(SF_midbins[ii])
#        SF_smf_mcmc.append(SF_smf[ii])
#        SF_error_mcmc.append(SF_error[ii])
#
## this does the exact same as above, but drops all zero-values in SF_smf
SF_midbins_mcmc = []
SF_smf_mcmc = []
SF_error_mcmc = []
for ii in range(len(SF_smf)):
    if SF_smf[ii] != 0:
        SF_midbins_mcmc.append(SF_midbins[ii])
        SF_smf_mcmc.append(SF_smf[ii])
        SF_error_mcmc.append(SF_error[ii])
#
## this does the exact same as above, but for Q_smf
Q_midbins_mcmc = []
Q_smf_mcmc = []
Q_error_mcmc = []
for ii in range(len(Q_smf)):
    if Q_smf[ii] != 0:
        Q_midbins_mcmc.append(SF_midbins[ii])
        Q_smf_mcmc.append(Q_smf[ii])
        Q_error_mcmc.append(Q_error[ii])
#
#
# convert lists to arrays
SF_smf_mcmc = np.array(SF_smf_mcmc)
SF_error_mcmc = np.array(SF_error_mcmc)
Q_smf_mcmc = np.array(Q_smf_mcmc)
Q_error_mcmc = np.array(Q_error_mcmc)
#
#
## define x array to generate points to plot Schechter fit
x = np.linspace(SF_midbins[0],SF_midbins[len(SF_midbins)-1],num=1000)#
#
#
### MCMC CODE:
#
#
#
## SECTION (2): inialize MCMC guess w/ least-squares estimate
#
# define guesses for initial guesses - SINGLE schechter
M_star_guess = 10.55
phi_guess = 0.03
alpha_guess = -1
#
# define guesses for initial guesses - DOUBLE schechter
#M_star_guess = 10.5
#phi_guess1 = 1
#alpha_guess1 = -1
#hi_guess2 = 1
#alpha_guess2 = -1
#
#
#
## open a file to print to
#
#
## check it directories to store outputs exist. if not, create them
output_dir = '/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/diagnostic_outputs/smf_fits/binning_method_%i'%membership_correction_binning_flag
check_folder = os.path.isdir(output_dir)
#
## If folder doesn't exist, then create it.
if not check_folder:
    os.makedirs(output_dir)
    if project_diagnostic_flag == 1 or diag_flag_master == 2:
        print("\nCreated folder : "+output_dir)
else:
    pass#print(output_dir, "\nfolder already exists.")
#
f = open('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/diagnostic_outputs/smf_fits/binning_method_%i'%membership_correction_binning_flag+'/z_spec_%.3f'%z_cutoff[0]+'_phot_cutoff_%.3f'%z_cutoff[1]+'_phot_cutoff_smf_fits.txt','w+')
#
#
## to be used in building strings throughout program
delim = ','   
## write a header for the file, start with hashtag to identify comment
header1 = 'z_spec_cutoff'+delim+'z_phot_cutoff'+delim+'type'+delim+'M_star'+delim+'M_star_sgima'+delim+'phi'+delim+'phi_sigma'+delim+'alpha'+delim+'alpha_sigma\n'
#
f.write(header1)
#
#
# define likelihood function for SINGLE schechter function; based on Eq.(2) of VDB 2013 et al
def lnlike(theta, midbins, smf, smf_error):
    M_star, phi, alpha = theta
    model = np.log(10)*phi*(10**((midbins-M_star)*(1+alpha)))*np.exp(-10**(midbins-M_star))
    return -0.5*np.sum((smf - model)**2 / smf_error**2)
#
# define likelihood function for DOUBLE schechter function; based on Eq.(11) of Weigel et al 2016
def lnlike2(theta, midbins, smf, smf_error):
    M_star, phi1, alpha1, phi2, alpha2 = theta
    model = np.log(10)*np.exp(-10**(midbins-M_star))*((phi1*(10**((midbins-M_star)*(1+alpha1))))+(phi2*(10**((midbins-M_star)*(1+alpha2)))))
    return -0.5*np.sum((smf - model)**2 / smf_error**2)

#
#
if SF_flag == 1:
    # optimize SF
    SF_nll = lambda *args: -lnlike(*args)       # single schechter fit - SF
    #
    ## fits to a single curve
    SFresult = op.minimize(SF_nll,[M_star_guess, phi_guess, alpha_guess], args=(SF_midbins_mcmc, SF_smf_mcmc, SF_error_mcmc), method='Nelder-Mead')
    SFM_ml, SFphi_ml, SFalpha_ml = SFresult['x']
    #
    # full fits
    SF_model_ml = np.log(10)*SFphi_ml*(10**((x-SFM_ml)*(1+SFalpha_ml)))*np.exp(-10**(x-SFM_ml))
#
if Q_flag ==1:
    # optimize Q
    Q_nll = lambda *args: -lnlike(*args)       # single schechter fit - Q
    #
    ## fits to a single curve
    Qresult = op.minimize(Q_nll,[M_star_guess, phi_guess, alpha_guess], args=(Q_midbins_mcmc, Q_smf_mcmc, Q_error_mcmc), method='Nelder-Mead')
    QM_ml, Qphi_ml, Qalpha_ml = Qresult['x']
    #
    # full fits
    Q_model_ml = np.log(10)*Qphi_ml*(10**((x-QM_ml)*(1+Qalpha_ml)))*np.exp(-10**(x-QM_ml))
#
#
if T_flag ==1:
    # optimize Total
    T_nll = lambda *args: -lnlike(*args)       # single schechter fit - total
    #
    ## fits to a single curve
    Tresult = op.minimize(T_nll,[M_star_guess, phi_guess, alpha_guess], args=(Q_midbins_mcmc, total_smf, total_error), method='Nelder-Mead')
    TM_ml, Tphi_ml, Talpha_ml = Tresult['x']
    #
    # full fits
    T_model_ml = np.log(10)*Tphi_ml*(10**((x-TM_ml)*(1+Talpha_ml)))*np.exp(-10**(x-TM_ml))
#
#
## generate points to plot Schechter fit
x = np.linspace(SF_midbins[0],SF_midbins[len(SF_midbins)-1],num=1000)#
#
# 
#
###################
###################
#
## Uncertainty estimation
#
# define "prior" function - SINGLE schechter
def lnprior(theta):
    M_star, phi, alpha = theta
    if 9.5 < M_star < 15 and -10 < phi < 30 and -2 < alpha < -0.5:
        return 0.0
    else:
        return -np.inf
#
# define "prior" function - DOUBLE schechter
def lnprior2(theta):
    M_star, phi1, alpha1, phi2, alpha2 = theta
    if 9.5 < M_star < 15 and -10 < phi1 < 30 and -2 < alpha1 < -0.5 and -10 < phi2 < 30 and -2 < alpha2 < -0.5:
        return 0.0
    else:
        return -np.inf
#
# combine with lnlike to write a function for the probability function - SINGLE schchter:
def lnprob(theta, midbins, smf, error):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    else:
        return lp + lnlike(theta, midbins, smf,error)
#
###### combine with lnlike to write a function for the probability function - DOUBLE schechter:
#####def lnprob2(theta, midbins, smf, error):
#####    lp = lnprior2(theta)
#####    if not np.isfinite(lp):
#####        return -np.inf
#####    else:
#####        return lp + lnlike2(theta, midbins, smf,error)
#
#
## STAR-FORMING LOOP
if SF_flag ==1:
    #
    print('Starting SF loop for ',nwalkers,' walkers and ',nsteps,' steps...')
#
    t0 = time.time()       # start timer to run emcee.EnsembleSampler()
    ## setup initial positions of walkers in a Guassian ball around the least-squares position    
    SFpos = SFresult['x'] + 1e-4*np.random.randn(nwalkers, ndim)
    # setup the emcee sampler
    SFsampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(SF_midbins_mcmc, SF_smf_mcmc,SF_error_mcmc))
    SFsampler.run_mcmc(SFpos, nsteps)  #run sampler from initial position "pos" for "nsteps" steps
    #
    t1 = time.time()       # stop timer
    total_time1 = t1-t0     # calculate time elapsed
    hours = total_time1 / 3600
    print('SF loop done.')
    print('Time elapsed for SF: ', total_time1, ' seconds, or: ',hours, ' hours.') # display time to run emcee Sampler   
    #
    ### delete burn-in (by specifying which cells to keep from *sampler.chain) and reshape for plotting corner plots
    print('Starting chain flattening...')
    t0 = time.time()       # start timer to flatten chain
    SFsamples = SFsampler.chain[:, :, :].reshape((-1, ndim))
    t1 = time.time()       # stop timer
    total_time2 = t1-t0     # calculate time elapsed
    hours = total_time2 / 3600
    print('Chain flatenning done.')
    print('Time elapsed: ', total_time2, ' seconds, or: ',hours, ' hours.') # display time to run
    #
    print('Creating corner plot & displaying results...')
    t0 = time.time()       # start timer to create corner plot
    ## create corner plots
    SFfig = corner.corner(SFsamples, labels=["$M^*$", "$\phi^*$", "alpha"])#  SF
    #
    # initialize result arrays
    SFresult_means = np.array([0]*ndim)
    SFrestult_sigmas = np.array([0]*ndim)
    SFresult_means = SFsamples.mean(axis=0)
    SFresult_sigmas = SFsamples.std(axis=0)
    #display results
    # SF
    print('SF result:')
    print('Means: ',SFresult_means)
    print('Sigmas: ',SFresult_sigmas)
    #
    SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = SFsamples.mean(axis=0)
    #
    #
    #
    bin_entry1 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'SF'+delim+str(np.round(SFM_star_mcmc,decimals=3))+delim+str(np.round(SFresult_sigmas[0],decimals=3))+delim+str(np.round(SFphi_mcmc,decimals=3))+delim+str(np.round(SFresult_sigmas[1],decimals=3))+delim+str(np.round(SFalpha_mcmc,decimals=3))+delim+str(np.round(SFresult_sigmas[2],decimals=3))
    writer = '%s'%str(bin_entry1)+'\n'
    f.write(writer)
    #
    #
    t1 = time.time()       # stop timer
    total_time3 = t1-t0     # calculate time elapsed
    hours = total_time3 / 3600
    print('Corner plot created, results displayed.')
    print('Time elapsed: ', total_time3, ' seconds, or: ',hours, ' hours.') # display time to run
    #
    total_run_time = total_time1 + total_time2 + total_time3
    hours = total_run_time / 3600
    print('Total run time for SF: ', total_run_time, ' seconds, or: ',hours, ' hours.') # display time to run
#
## Q LOOP
if Q_flag ==1:
    # initialize walkers
    #ndim, nwalkers, nsteps = 5, 10000, 200
    print('Starting Q loop for ',nwalkers,' walkers and ',nsteps,' steps...')
#
    t0 = time.time()       # start timer to run emcee.EnsembleSampler()
    ## setup initial positions of walkers in a Guassian ball around the least-squares position    
    Qpos = Qresult['x'] + 1e-4*np.random.randn(nwalkers, ndim)
    # setup the emcee sampler
    Qsampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(Q_midbins_mcmc, Q_smf_mcmc,Q_error_mcmc))
    Qsampler.run_mcmc(Qpos, nsteps)  #run sampler from initial position "pos" for "nsteps" steps
    #
    t1 = time.time()       # stop timer
    total_time1 = t1-t0     # calculate time elapsed
    hours = total_time1 / 3600
    print('Q loop done.')
    print('Time elapsed for Q: ', total_time1, ' seconds, or: ',hours, ' hours.') # display time to run emcee Sampler   
    #
    ### delete burn-in (by specifying which cells to keep from *sampler.chain) and reshape for plotting corner plots
    print('Starting chain flattening...')
    t0 = time.time()       # start timer to flatten chain
    Qsamples = Qsampler.chain[:, :, :].reshape((-1, ndim))
    t1 = time.time()       # stop timer
    total_time2 = t1-t0     # calculate time elapsed
    hours = total_time2 / 3600
    print('Chain flatenning done.')
    print('Time elapsed: ', total_time2, ' seconds, or: ',hours, ' hours.') # display time to run
    #
    print('Creating corner plot & displaying results...')
    t0 = time.time()       # start timer to create corner plot
    ## create corner plots
    Qfig = corner.corner(Qsamples, labels=["$M^*$", "$\phi^*_1$", "$\alpha_1$", "$\phi^*_2$", "$\alpha_2$"])#  Q
    #
    # initialize result arrays
    Qresult_means = np.array([0]*ndim)
    Qrestult_sigmas = np.array([0]*ndim)
    Qresult_means = Qsamples.mean(axis=0)
    Qresult_sigmas = Qsamples.std(axis=0)
    #display results
    # SF
    print('Q result:')
    print('Means: ',Qresult_means)
    print('Sigmas: ',Qresult_sigmas)
    #
    QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = Qsamples.mean(axis=0)
    #
    ## print OUTPUT to file
    bin_entry1 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'Q'+delim+str(np.round(QM_star_mcmc,decimals=3))+delim+str(np.round(Qresult_sigmas[0],decimals=3))+delim+str(np.round(Qphi_mcmc,decimals=3))+delim+str(np.round(Qresult_sigmas[1],decimals=3))+delim+str(np.round(Qalpha_mcmc,decimals=3))+delim+str(np.round(Qresult_sigmas[2],decimals=3))
    writer = '%s'%str(bin_entry1)+'\n'
    f.write(writer)
    #
    #
    t1 = time.time()       # stop timer
    total_time3 = t1-t0     # calculate time elapsed
    hours = total_time3 / 3600
    print('Corner plot created, results displayed.')
    print('Time elapsed: ', total_time3, ' seconds, or: ',hours, ' hours.') # display time to run
    #
    total_run_time = total_time1 + total_time2 + total_time3
    hours = total_run_time / 3600
    print('Total run time for Q: ', total_run_time, ' seconds, or: ',hours, ' hours.') # display time to run
#
## TOTAL POPULATION LOOP
if T_flag ==1:
    # initialize walkers
    #ndim, nwalkers, nsteps = 3, 1000, 2000
    print('Starting Total loop for ',nwalkers,' walkers and ',nsteps,' steps...')
#
    t0 = time.time()       # start timer to run emcee.EnsembleSampler()
    ## setup initial positions of walkers in a Guassian ball around the least-squares position    
    Tpos = Tresult['x'] + 1e-4*np.random.randn(nwalkers, ndim)
    # setup the emcee sampler
    Tsampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(SF_midbins, total_smf, total_error))
    Tsampler.run_mcmc(Tpos, nsteps)  #run sampler from initial position "pos" for "nsteps" steps
    #
    t1 = time.time()       # stop timer
    total_time1 = t1-t0     # calculate time elapsed
    hours = total_time1 / 3600
    print('Total loop done.')
    print('Time elapsed for Total: ', total_time1, ' seconds, or: ',hours, ' hours.') # display time to run emcee Sampler   
    #
    ### delete burn-in (by specifying which cells to keep from *sampler.chain) and reshape for plotting corner plots
    print('Starting chain flattening...')
    t0 = time.time()       # start timer to flatten chain
    Tsamples = Tsampler.chain[:, :, :].reshape((-1, ndim))
    t1 = time.time()       # stop timer
    total_time2 = t1-t0     # calculate time elapsed
    hours = total_time2 / 3600
    print('Chain flatenning done.')
    print('Time elapsed: ', total_time2, ' seconds, or: ',hours, ' hours.') # display time to run
    #
    print('Creating corner plot & displaying results...')
    t0 = time.time()       # start timer to create corner plot
    ## create corner plots
    Tfig = corner.corner(Tsamples, labels=["$M^*$", "$\phi^*$", "alpha"])#  SF
    #
    # initialize result arrays
    Tresult_means = np.array([0]*ndim)
    Trestult_sigmas = np.array([0]*ndim)
    Tresult_means = Tsamples.mean(axis=0)
    Tresult_sigmas = Tsamples.std(axis=0)
    #display results
    # SF
    print('Total result:')
    print('Means: ',Tresult_means)
    print('Sigmas: ',Tresult_sigmas)
    #
    TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = Tsamples.mean(axis=0)
    #
    #
    ## print OUTPUT to file
    bin_entry1 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'Total'+delim+str(np.round(TM_star_mcmc,decimals=3))+delim+str(np.round(Tresult_sigmas[0],decimals=3))+delim+str(np.round(Tphi_mcmc,decimals=3))+delim+str(np.round(Tresult_sigmas[1],decimals=3))+delim+str(np.round(Talpha_mcmc,decimals=3))+delim+str(np.round(Tresult_sigmas[2],decimals=3))
    writer = '%s'%str(bin_entry1)+'\n'
    f.write(writer)
    #
    #
    t1 = time.time()       # stop timer
    total_time3 = t1-t0     # calculate time elapsed
    hours = total_time3 / 3600
    print('Corner plot created, results displayed.')
    print('Time elapsed: ', total_time3, ' seconds, or: ',hours, ' hours.') # display time to run
    #
    total_run_time = total_time1 + total_time2 + total_time3
    hours = total_run_time / 3600
    print('Total run time: ', total_run_time, ' seconds, or: ',hours, ' hours.') # display time to run
#  
#
#
#
#
## close the file       
f.close()
#
#
print('\n\n"emcee_chi2*.py"  terminated successfully.\n')
#
#
#                        
###################     PROGRAM END     ###################