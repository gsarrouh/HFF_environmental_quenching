## (vi) fit a SCHECHTER function to scatter plot
#
import pandas as pd
import emcee
import scipy.optimize as op
import scipy.integrate as integrate


#

### MCMC CODE:
#
#
M_star_guess = 10
phi_guess = 20
alpha_guess = -1
#
## define likelihood function
## single schechter fit - SF
def schechter(theta,midbins):
    M_star, phi, alpha = theta
    log_M_star = np.log10(M_star)
    log_M = np.log10(midbins)
    schechter = np.exp(-10**(log_M - log_M_star))*(10**(log_M-log_M_star))**(1+alpha)
    return schechter

def lnlike(theta,midbins):
    M_star, phi, alpha = theta
    log_M_star = np.log10(M_star)
    log_M = np.log10(midbins)
    log_L = -(1/M_star)*np.sum(midbins) + (1 + alpha)*np.sum(np.log(midbins/M_star) - np.sum(integrate.quad(lambda log_M: schechter, np.log10(7.3),np.log(12.15), args = (theta,midbins))

#def lnlike(theta, midbins, smf, smf_error):
#    M_star, phi, alpha = theta
#    model = np.log(10)*phi*(10**((midbins-M_star)*(1+alpha)))*np.exp(-10**(midbins-M_star))
#    return np.sum((smf - model)**2 / smf_error**2)

#

# optimize
nll = lambda *args: lnlike(*args)       # single schechter fit - SF
#
## fits to a single curve
SFresult = op.minimize(nll,[M_star_guess, phi_guess, alpha_guess], args=(SF_midbins, SF_smf, SF_error), method='Nelder-Mead')
SFM_ml, SFphi_ml, SFalpha_ml = SFresult['x']
#
## generate points to plot Schechter fit
x = np.linspace(SF_midbins[0],SF_midbins[len(SF_midbins)-1],num=1000)#
#
# full fits
SF_model = np.log(10)*SFphi_ml*(10**((x-SFM_ml)*(1+SFalpha_ml)))*np.exp(-10**(x-SFM_ml))
#
###################
###################
#
## Uncertainty estimation
#
# define "prior" function
def lnprior(theta):
    M_star, phi, alpha = theta
    if 9.5 < M_star < 15 and -10 < phi < 30 and -2 < alpha < -0.5:
        return 0.0
    else:
        return -np.inf
#
# combine with lnlike to write a function for the probability function:
def lnprob(theta, midbins, smf):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    else:
        return lp + lnlike(theta, midbins, smf)
#
# initialize walkers
ndim, nwalkers, nsteps = 3, 100, 1000
## setup initial positions of walkers in a Guassian ball around the least-squares position    
SFpos = SFresult['x'] + 1e-4*np.random.randn(nwalkers, ndim)
#
# setup the emcee sampler
SFsampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(SF_midbins, SF_smf))
SFsampler.run_mcmc(SFpos, nsteps)  #run sampler from initial position "pos" for "nsteps" steps
#
## visualize the paths of the walkers
#
fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
samples = SFsampler.chain             # returns a nwalkers x nsteps x ndim array, showing the position of each walker at each step
labels = ["M_star", "phi", "alpha"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number");
#
## delete burn-in
SFsamples = SFsampler.chain[:, 50:, :].reshape((-1, ndim))
#################
## create corner plots
import corner
SFfig = corner.corner(SFsamples, labels=["$M^*$", "$\phi^*$", "alpha"])#  SF
#                      truths=[m_true, b_true, np.log(f_true)])
Qfig = corner.corner(Qsamples, labels=["$M^*_1$", "$\phi^*_1$", "alpha_1","$M^*_2$", "$\phi^*_2$", "alpha_2"]) # Q
#
#
## get values for all parameters:
M_star_mcmc, phi_mcmc, alpha_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(SFsamples, [66, 95, 99.7],
                                                axis=0)))
#
#
#
#
#
#