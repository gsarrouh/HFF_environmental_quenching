# the following script runs the "emcee" MCMC sampler, using a chi-squared metric as the cost function: chi^2 = -0.5 * ( (data - model)^2 / error^2 )
#
## (vi) fit a SCHECHTER function to scatter plot
#
import pandas as pd
import emcee
import time
import scipy.optimize as op
import scipy.integrate as integrate
import corner
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
### MCMC CODE:
#
#
######
###### MAY NEED TO EDIT ######
######
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
######
######
#
#
#
#
######
# define guesses for initial guesses - SINGLE schechter
M_star_guess = 10
phi_guess = 1
alpha_guess = -1
#
# define guesses for initial guesses - DOUBLE schechter
M_star_guess = 10
phi_guess1 = 1
alpha_guess1 = -1
phi_guess2 = 1
alpha_guess2 = -1
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
    Q_nll = lambda *args: -lnlike2(*args)       # single schechter fit - Q
    #
    ## fits to a single curve
    Qresult = op.minimize(Q_nll,[M_star_guess, phi_guess1, alpha_guess1, phi_guess2, alpha_guess2], args=(Q_midbins_mcmc, Q_smf_mcmc, Q_error_mcmc), method='Nelder-Mead')
    QM_ml, Qphi1_ml, Qalpha1_ml, Qphi2_ml, Qalpha2_ml = Qresult['x']
    #
    # full fits
    Q_model_ml = np.log(10)*np.exp(-10**(x-QM_ml))*((Qphi1_ml*(10**((x-QM_ml)*(1+Qalpha1_ml))))+(Qphi2_ml*(10**((x-QM_ml)*(1+Qalpha2_ml)))))
#
#
if T_flag ==1:
    # optimize Total
    T_nll = lambda *args: -lnlike(*args)       # single schechter fit - total
    #
    ## fits to a single curve
    Tresult = op.minimize(T_nll,[M_star_guess, phi_guess, alpha_guess], args=(SF_midbins, total_smf, total_error), method='Nelder-Mead')
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
# combine with lnlike to write a function for the probability function - DOUBLE schechter:
def lnprob2(theta, midbins, smf, error):
    lp = lnprior2(theta)
    if not np.isfinite(lp):
        return -np.inf
    else:
        return lp + lnlike2(theta, midbins, smf,error)
#
#
## STAR-FORMING LOOP
if SF_flag ==1:
    # initialize walkers
    ndim, nwalkers, nsteps = 3, 100, 500000
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
    ndim, nwalkers, nsteps = 5, 10000, 200
    print('Starting Q loop for ',nwalkers,' walkers and ',nsteps,' steps...')
#
    t0 = time.time()       # start timer to run emcee.EnsembleSampler()
    ## setup initial positions of walkers in a Guassian ball around the least-squares position    
    Qpos = Qresult['x'] + 1e-4*np.random.randn(nwalkers, ndim)
    # setup the emcee sampler
    Qsampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob2, args=(Q_midbins_mcmc, Q_smf_mcmc,Q_error_mcmc))
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
    QM_star_mcmc, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = Qsamples.mean(axis=0)
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
    ndim, nwalkers, nsteps = 3, 100, 200000
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
#################
#################  END  #################
#################
#
#
#
#
### visualize the paths of the walkers
###
#fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
#samples = SFsampler.chain             # returns a nwalkers x nsteps x ndim array, showing the position of each walker at each step
#labels = ["M_star", "phi", "alpha"]
#for i in range(ndim):
#    ax = axes[i]
#    ax.plot(samples[:, :, i], "k", alpha=0.3)
#    ax.set_xlim(0, nsteps)
#    ax.set_ylabel(labels[i])
#    ax.yaxis.set_label_coords(-0.1, 0.5)
#
#axes[-1].set_xlabel("step number");
#
#################
## get values for all parameters:
#####SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(SFsamples, [68, 95, 99.7],axis=0)))
#
# print result
#from IPython.display import display, Math
#



# 







#
#for i in range(ndim):
#    mcmc = np.percentile(SFsamples[:, i], [68, 95, 99.7])
#    q = np.diff(mcmc)
#    txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
#    txt = txt.format(mcmc[1], q[0], q[1], labels[i])
#    display(Math(txt))
#
#
## build a model for plotting
#SF_model_mcmc = np.log(10)*SFphi_mcmc*(10**((x-SFM_star_mcmc)*(1+SFalpha_mcmc)))*np.exp(-10**(x-SFM_star_mcmc))
#Q_model_mcmc = np.log(10)*Qphi_mcmc*(10**((x-QM_star_mcmc)*(1+Qalpha_mcmc)))*np.exp(-10**(x-QM_star_mcmc))
#T_model_mcmc = np.log(10)*Tphi_mcmc*(10**((x-TM_star_mcmc)*(1+Talpha_mcmc)))*np.exp(-10**(x-TM_star_mcmc))
#
# create arrays for plotting purposes to display values down to a y_min = 1
#SF_model_ml_plot = []
#SF_model_mcmc_plot = []
#x_plot_SF = []
#for ii in range(len(SF_model_ml)):
#    if SF_model_ml[ii] > 1:
#        SF_model_ml_plot.append(SF_model_ml[ii])
#        SF_model_mcmc_plot.append(SF_model_mcmc[ii])
#        x_plot_SF.append(x[ii])
# do same for Q pop
#Q_model_ml_plot = []
#Q_model_mcmc_plot = []
#x_plot_Q = []
#for ii in range(len(Q_model_ml)):
#    if Q_model_ml[ii] > 1:
#        Q_model_ml_plot.append(Q_model_ml[ii])
#        Q_model_mcmc_plot.append(Q_model_mcmc[ii])
#        x_plot_Q.append(x[ii])
#
# do same for total pop
#T_model_ml_plot = []
#T_model_mcmc_plot = []
#x_plot_T = []
#for ii in range(len(T_model_ml)):
#    if T_model_ml[ii] > 1:
#        T_model_ml_plot.append(T_model_ml[ii])
#        T_model_mcmc_plot.append(T_model_mcmc[ii])
#        x_plot_T.append(x[ii])
#