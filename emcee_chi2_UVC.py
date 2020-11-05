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
## FLAGS: 0 = off (i.e. skip), 1 = on (i.e. execute code)
#
#
######
###### MAY NEED TO EDIT ######
######
#
SF_flag = 1          # star-forming population
Q_flag = 1           # Q pop
T_flag = 0           # total pop
#
field_flag = 1              # 0 = off, fit CLUSTER SMFs; 1 = on, fit FIELD SMFs
#
#
#### MISCELANEOUS Section: if fitting the field population, rename the arrays so you don't have to change the rest of the coded
if field_flag == 1:
    # SF_midbins = SF_field_midbins
    # SF_smf = SF_field_smf
    # SF_error = SF_field_error
    # Q_smf = Q_field_smf
    # Q_error = Q_field_error
    # total_smf = total_field_smf
    # total_error = total_field_error
    # SF_midbins = SF_field_midbins
    SF_smf = SF_field_smf_UVC
    SF_error = SF_field_smf_relerr_UVC * SF_field_smf_UVC
    Q_smf = Q_field_smf_UVC
    Q_error = Q_field_smf_relerr_UVC * Q_field_smf_UVC
    #total_smf = SF_field_smf_UVC + Q_field_smf_UVC
    #total_error = total_field_error
    #
    SF_midbins = np.delete(SF_midbins,[0,1,2,3])
    SF_smf = np.delete(SF_smf,[0,1,2,3])
    SF_error = np.delete(SF_error,[0,1,2,3])
    Q_smf = np.delete(Q_smf,[0,1,2,3])
    Q_error = np.delete(Q_error,[0,1,2,3])
    #total_smf = np.delete(total_smf,[0,1,2,3])
    #total_error = np.delete(total_error,[0,1,2,3])
#
######
###### MAY NEED TO EDIT ######
######
#
# initialize walkers
ndim, nwalkers, max_nsteps = 3, 500, 10000    # (# of parameters), (#walkers), (#steps/walker)
#
#
labels = ["M*","phi","alpha"]
## FUNCTIONS
#
#
# define likelihood function for SINGLE schechter function; based on Eq.(2) of VDB 2013 et al
def lnlike(theta, midbins, smf, smf_error):
    M_star, phi, alpha = theta
    midbins = np.array(midbins)
    smf = np.array(smf)
    smf_error = np.array(smf_error)
    model = np.array( np.log(10) * phi * (10**((midbins-M_star)*(1+alpha))) * np.exp(-10**(midbins-M_star)) )
    return -0.5*np.sum( ( (smf - model) / smf_error)**2 )
#
## define likelihood function for DOUBLE schechter function; based on Eq.(11) of Weigel et al 2016
#def lnlike2(theta, midbins, smf, smf_error):
#    M_star, phi1, alpha1, phi2, alpha2 = theta
#    model = np.log(10)*np.exp(-10**(midbins-M_star))*((phi1*(10**((midbins-M_star)*(1+alpha1))))+(phi2*(10**((midbins-M_star)*(1+alpha2)))))
#    return -0.5*np.sum((smf - model)**2 / smf_error**2)
#
#
# define "prior" function - SINGLE schechter
def lnprior(theta):
    M_star, phi, alpha = theta
    if (9.5 < M_star < 14) and (-50 < phi < 50) and (-5 < alpha < 2):
        return 0.0
    else:
        return -np.inf
#
## define "prior" function - DOUBLE schechter
#def lnprior2(theta):
#    M_star, phi1, alpha1, phi2, alpha2 = theta
#    if 9.5 < M_star < 15 and -10 < phi1 < 30 and -2 < alpha1 < -0.5 and -10 < phi2 < 30 and -2 < alpha2 < -0.5:
#        return 0.0
#    else:
#        return -np.inf
#
# combine with lnlike to write a function for the probability function - SINGLE schchter:
def lnprob(theta, midbins, smf, smf_error):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    else:
        return lp + lnlike(theta, midbins, smf, smf_error)
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
## SECTION (1): setup lists
#
### SETUP appropriate lists and convert them to arrays for use in MCMC code:
## this first bit replaces 'NaN' entries in SF_error (errors to the SMF) with 1e5. This only applies to pixels for which the SMF count = 0, hence the error (given by sqrt(smf) is NaN. It is done for the purpose of computing chi_squared later on
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
#array_list = [SF_smf_mcmc,SF_error_mcmc,SF_midbins_mcmc,Q_smf_mcmc,Q_error_mcmc,Q_midbins_mcmc]
#for x in array_list:
#    x = np.array(x)
SF_smf_mcmc = np.array(SF_smf_mcmc)
SF_error_mcmc = np.array(SF_error_mcmc)
SF_midbins_mcmc = np.array(SF_midbins_mcmc)
Q_smf_mcmc = np.array(Q_smf_mcmc)
Q_error_mcmc = np.array(Q_error_mcmc)
Q_midbins_mcmc = np.array(Q_midbins_mcmc)
#
## Re-shape arrays
SF_midbins_mcmc = SF_midbins_mcmc.reshape(len(SF_midbins_mcmc))
Q_midbins_mcmc = Q_midbins_mcmc.reshape(len(Q_midbins_mcmc))
#
#
#
### MCMC CODE:
#
## PT 1: intialize guess randomly (user-defined), then with maximum-likelihood (ML) method
#
#
## SECTION (2): inialize MCMC guess w/ least-squares estimate
#
# define guesses for initial guesses - SINGLE schechter
M_star_guess = 10.0
phi_guess = 0.005
alpha_guess = -0.6
#
theta_guess = [M_star_guess,phi_guess,alpha_guess]
# define guesses for initial guesses - DOUBLE schechter
#M_star_guess = 10.5
#phi_guess1 = 1
#alpha_guess1 = -1
#phi_guess2 = 1
#alpha_guess2 = -1
#
#
#
## open a file to print to
#
#
## check if directories to store outputs exist. if not, create them
output_dir = '/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/diagnostic_outputs/smf_fits/binning_method_%i'%membership_correction_binning_flag
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
if field_flag == 0:
    f = open('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/diagnostic_outputs/smf_fits/binning_method_%i'%membership_correction_binning_flag+'/z_spec_%.3f'%z_cutoff[0]+'_z_phot_%.3f'%z_cutoff[1]+'_walkers_%i'%nwalkers+'_smf_fits_UVC.txt','w+')
elif field_flag == 1:
    f = open('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/diagnostic_outputs/smf_fits/binning_method_%i'%membership_correction_binning_flag+'/z_spec_%.3f'%z_cutoff[0]+'_z_phot_%.3f'%z_cutoff[1]+'_walkers_%i'%nwalkers+'_field_smf_fits_UVC.txt','w+')
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
if SF_flag == 1:
    # optimize SF
    SF_nll = lambda *args: -lnlike(*args)       # single schechter fit - SF
    #
    ## fits to a single curve
    SFresult = op.minimize(SF_nll,[M_star_guess, phi_guess, alpha_guess], args=(SF_midbins_mcmc, SF_smf_mcmc, SF_error_mcmc), method='Nelder-Mead')
    SFM_ml, SFphi_ml, SFalpha_ml = SFresult['x']
    #
    print("SFResult - max. likelihood:")
    print(SFresult['x'])
    #
    # full fits
    # SF_model_ml = np.log(10)*SFphi_ml*(10**((x-SFM_ml)*(1+SFalpha_ml)))*np.exp(-10**(x-SFM_ml))
#
if Q_flag ==1:
    # optimize Q
    Q_nll = lambda *args: -lnlike(*args)       # single schechter fit - Q
    #
    ## fits to a single curve
    Qresult = op.minimize(Q_nll,[M_star_guess, phi_guess, alpha_guess], args=(Q_midbins_mcmc, Q_smf_mcmc, Q_error_mcmc), method='Nelder-Mead')
    QM_ml, Qphi_ml, Qalpha_ml = Qresult['x']
    #
    print("\nQResult - max. likelihood:")
    print(Qresult['x'])
    #
    # full fits
    # Q_model_ml = np.log(10)*Qphi_ml*(10**((x-QM_ml)*(1+Qalpha_ml)))*np.exp(-10**(x-QM_ml))
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
    print("\nTResult - max. likelihood:")
    print(Tresult['x'])
    #
    # full fits
    # T_model_ml = np.log(10)*Tphi_ml*(10**((x-TM_ml)*(1+Talpha_ml)))*np.exp(-10**(x-TM_ml))
#
#
## generate points to plot Schechter fit
# x = np.linspace(SF_midbins[0],SF_midbins[len(SF_midbins)-1],num=1000)#
#
##
#
#
###################
###################
#
## MCMC & Uncertainty estimation
#
## STAR-FORMING LOOP
if SF_flag ==1:
    #
    print('\nStarting SF loop for ',nwalkers,' walkers and ',max_nsteps,' steps...')
#
    t0 = time.time()       # start timer to run emcee.EnsembleSampler()
    ## setup initial positions of walkers in a Guassian ball around the least-squares position
    SFpos = SFresult['x'] + 1e-4*np.random.randn(nwalkers, ndim)
    #
    # Set up the backend
    # Don't forget to clear it in case the file already exists
    filename = '/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/diagnostic_outputs/smf_fits/binning_method_%i'%membership_correction_binning_flag+'/sf_smf_UVC.h5'
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)
    #
    # setup the emcee sampler
    SFsampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(SF_midbins_mcmc, SF_smf_mcmc,SF_error_mcmc))
    #
    #
    # We'll track how the average autocorrelation time estimate changes
    index = 0
    autocorr = np.empty(max_nsteps)

    # This will be useful to testing convergence
    old_tau = np.inf

    # Now we'll sample for up to max_n steps
    for sample in SFsampler.sample(SFpos, iterations=max_nsteps, progress=True):
        # Only check convergence every 100 steps
        if SFsampler.iteration % 100:
            continue

        # Compute the autocorrelation time so far
        # Using tol=0 means that we'll always get an estimate even
        # if it isn't trustworthy
        tau = SFsampler.get_autocorr_time(tol=0)
        autocorr[index] = np.mean(tau)
        index += 1
        #
        # Check convergence
        converged = np.all( (tau * 100) < SFsampler.iteration)
        converged &= np.all( (np.abs(old_tau - tau) / tau) < 0.01)
        if converged:
            print('# of iterations: %i'%SFsampler.iteration)
            break
        old_tau = tau
    #
    #
    #
    #SFsampler.run_mcmc(SFpos, max_nsteps, progress = True)  #run sampler from initial position "pos" for "max_nsteps" steps
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
    SFfig.suptitle('Star-forming')
    if field_flag == 0:
        plt.savefig('/Users/gsarrouh/Research/NSERC_2017_HFF/Plots/SMF/corner_plots/z_spec_%.3f'%z_cutoff[0]+'_z_phot_%.3f'%z_cutoff[1]+'_walkers_%i'%nwalkers+'_steps_%i'%max_nsteps+'_SF_corner_UVC.png')
    elif field_flag == 1:
        plt.savefig('/Users/gsarrouh/Research/NSERC_2017_HFF/Plots/SMF/corner_plots/z_spec_%.3f'%z_cutoff[0]+'_z_phot_%.3f'%z_cutoff[1]+'_walkers_%i'%nwalkers+'_steps_%i'%max_nsteps+'_SF_field_corner_UVC.png')
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
    ## display result w/ 1-sigma uncertainty (16th, 50th, 86th percentiles)
    from IPython.display import display, Math
    #
    for i in range(ndim):
        mcmc = np.percentile(SFsamples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        txt = "\mathrm{{{3}}} = {0:.7f}_{{-{1:.7f}}}^{{{2:.7f}}}"
        txt = txt.format(mcmc[1], q[0], q[1], labels[i])
        display(Math(txt))
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
    #ndim, nwalkers, max_nsteps = 5, 10000, 200
    print('\nStarting Q loop for ',nwalkers,' walkers and ',max_nsteps,' steps...')
#
    t0 = time.time()       # start timer to run emcee.EnsembleSampler()
    ## setup initial positions of walkers in a Guassian ball around the least-squares position
    # Qpos = theta_guess + 1e-4*np.random.randn(nwalkers, ndim)
    Qpos = Qresult['x'] + 1e-4*np.random.randn(nwalkers, ndim)
    #
    # Set up the backend
    # Don't forget to clear it in case the file already exists
    filename = '/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/diagnostic_outputs/smf_fits/binning_method_%i'%membership_correction_binning_flag+'/q_smf_UVC.h5'
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)
    #
    # setup the emcee sampler
    Qsampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(Q_midbins_mcmc, Q_smf_mcmc,Q_error_mcmc))
    #
    #
    # We'll track how the average autocorrelation time estimate changes
    index = 0
    autocorr = np.empty(max_nsteps)

    # This will be useful to testing convergence
    old_tau = np.inf

    # Now we'll sample for up to max_n steps
    for sample in Qsampler.sample(Qpos, iterations=max_nsteps, progress=True):
        # Only check convergence every 100 steps
        if Qsampler.iteration % 100:
            continue

        # Compute the autocorrelation time so far
        # Using tol=0 means that we'll always get an estimate even
        # if it isn't trustworthy
        tau = Qsampler.get_autocorr_time(tol=0)
        autocorr[index] = np.mean(tau)
        index += 1
        #
        # Check convergence
        converged = np.all( (tau * 100) < Qsampler.iteration)
        converged &= np.all( (np.abs(old_tau - tau) / tau) < 0.01)
        if converged:
            print('# of iterations: %i'%Qsampler.iteration)
            break
        old_tau = tau
    #
    #
    # Qsampler.run_mcmc(Qpos, max_nsteps, progress = True)  #run sampler from initial position "pos" for "max_nsteps" steps
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
    Qfig = corner.corner(Qsamples, labels=["$M^*$", "$\phi$", "alpha"])#, "$\phi^*_2$", "$\alpha_2$"])#  Q
    Qfig.suptitle('Quiescent')
    #
    if field_flag == 0:
        plt.savefig('/Users/gsarrouh/Research/NSERC_2017_HFF/Plots/SMF/corner_plots/z_spec_%.3f'%z_cutoff[0]+'_z_phot_%.3f'%z_cutoff[1]+'_walkers_%i'%nwalkers+'_steps_%i'%max_nsteps+'_Q_corner_UVC.png')
    elif field_flag == 1:
        plt.savefig('/Users/gsarrouh/Research/NSERC_2017_HFF/Plots/SMF/corner_plots/z_spec_%.3f'%z_cutoff[0]+'_z_phot_%.3f'%z_cutoff[1]+'_walkers_%i'%nwalkers+'_steps_%i'%max_nsteps+'_Q_field_corner_UVC.png')
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
    #
    ## display result w/ 1-sigma uncertainty (16th, 50th, 86th percentiles)
    from IPython.display import display, Math
    #
    for i in range(ndim):
        mcmc = np.percentile(Qsamples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        txt = "\mathrm{{{3}}} = {0:.7f}_{{-{1:.7f}}}^{{{2:.7f}}}"
        txt = txt.format(mcmc[1], q[0], q[1], labels[i])
        display(Math(txt))
    #
    #
    ## investigate movement of walkers
    plt.figure()
    plt.plot(Tsampler.get_chain()[:, 0, 0], "k", lw=0.5)
    #plt.xlim(0, 5000)
    #plt.ylim(-5.5, 5.5)
    plt.title("move: StretchMove", fontsize=14)
    plt.xlabel("step number")
    plt.ylabel("x");
    #
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
    #ndim, nwalkers, max_nsteps = 3, 1000, 2000
    print('\nStarting TOTAL loop for ',nwalkers,' walkers and ',max_nsteps,' steps...')
#
    t0 = time.time()       # start timer to run emcee.EnsembleSampler()
    ## setup initial positions of walkers in a Guassian ball around the least-squares position
    # Tpos = theta_guess + 1e-4*np.random.randn(nwalkers, ndim)
    Tpos = Tresult['x'] + 1e-4*np.random.randn(nwalkers, ndim)
    #
    # Set up the backend
    # Don't forget to clear it in case the file already exists
    if field_flag == 0:
        filename = '/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/diagnostic_outputs/smf_fits/binning_method_%i'%membership_correction_binning_flag+'/total_smf_UVC.h5'
    elif field_flag == 1:
        filename = '/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/diagnostic_outputs/smf_fits/binning_method_%i'%membership_correction_binning_flag+'/total_field_smf_UVC.h5'
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)
    #
    # setup the emcee sampler
    Tsampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(Q_midbins_mcmc, total_smf, total_error), backend=backend)
    #
    # We'll track how the average autocorrelation time estimate changes
    index = 0
    autocorr = np.empty(max_nsteps)

    # This will be useful to testing convergence
    old_tau = np.inf

    # Now we'll sample for up to max_n steps
    for sample in Tsampler.sample(Tpos, iterations=max_nsteps, progress=True):
        # Only check convergence every 100 steps
        if Tsampler.iteration % 100:
            continue

        # Compute the autocorrelation time so far
        # Using tol=0 means that we'll always get an estimate even
        # if it isn't trustworthy
        tau = Tsampler.get_autocorr_time(tol=0)
        autocorr[index] = np.mean(tau)
        index += 1
        #
        # Check convergence
        converged = np.all( (tau * 100) < Tsampler.iteration)
        converged &= np.all( (np.abs(old_tau - tau) / tau) < 0.01)
        if converged:
            print('# of iterations: %i'%Tsampler.iteration)
            break
        old_tau = tau
    #
    # Tsampler.run_mcmc(Tpos, max_nsteps, progress = True)  #run sampler from initial position "pos" for "max_nsteps" steps
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
    Tfig = corner.corner(Tsamples, labels=["$M^*$", "$\phi^*$", "alpha"])#  total
    Tfig.suptitle('Total')
    #
    if field_flag == 0:
        plt.savefig('/Users/gsarrouh/Research/NSERC_2017_HFF/Plots/SMF/corner_plots/z_spec_%.3f'%z_cutoff[0]+'_z_phot_%.3f'%z_cutoff[1]+'_walkers_%i'%nwalkers+'_steps_%i'%max_nsteps+'_total_corner_UVC.png')
    elif field_flag == 1:
        plt.savefig('/Users/gsarrouh/Research/NSERC_2017_HFF/Plots/SMF/corner_plots/z_spec_%.3f'%z_cutoff[0]+'_z_phot_%.3f'%z_cutoff[1]+'_walkers_%i'%nwalkers+'_steps_%i'%max_nsteps+'_total_field_corner_UVC.png')
    #
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
    ## display result w/ 1-sigma uncertainty (16th, 50th, 86th percentiles)
    from IPython.display import display, Math
    #
    for i in range(ndim):
        mcmc = np.percentile(Tsamples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        txt = "\mathrm{{{3}}} = {0:.7f}_{{-{1:.7f}}}^{{{2:.7f}}}"
        txt = txt.format(mcmc[1], q[0], q[1], labels[i])
        display(Math(txt))
    #
    #
    ## investigate movement of walkers
    plt.figure()
    plt.plot(Tsampler.get_chain()[:, 0, 0], "k", lw=0.5)
    #plt.xlim(0, 5000)
    #plt.ylim(-5.5, 5.5)
    plt.title("move: StretchMove", fontsize=14)
    plt.xlabel("step number")
    plt.ylabel("x");
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
