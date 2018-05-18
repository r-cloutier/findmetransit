from imports import *
from priors import *

def box_transit_model(theta, t):
    '''Return the box transmit model.'''
    fmodel = np.ones(t.size)
    P, T0, depth, duration = theta
    phase = foldAt(t, P, T0)
    phase[phase>.5] -= 1
    intransit = (phase*P <= .5*duration) & (phase*P >= -.5*duration)
    fmodel[intransit] = 1. - depth
    return fmodel


def lnlike(theta, bjd, f, ef):
    fmodel = box_transit_model(theta, bjd)
    return -.5*(np.sum((f-fmodel)**2 / ef**2 - np.log(1./ef**2)))


def lnprior(theta, theta0):
    P, T0, Z, D = theta
    P0, T00, Z0, D0 = theta0
    lps = np.zeros(4)
    lps[0] = lnuniform(P, P0-2*D0, P0+2*D0)
    lps[1] = lnuniform(T0, T00-2*D0, T00+2*D0)
    lps[2] = lnuniform(Z, 0, 1)
    lps[3] = lnuniform(D, 0, 1)
    ##print P, T0, Z, D, lps
    return lps.sum()


def lnprob(theta, theta0, bjd, f, ef):
    lp = lnprior(theta, theta0)
    if np.isfinite(lp):
        return lp + lnlike(theta, bjd, f, ef)
    else:
        return -np.inf


def run_emcee(theta, theta0, bjd, f, ef, initialize, nwalkers=100, burnin=200,
              nsteps=200, a=2):
    '''Run mcmc on an input light curve with no transit model.'''
    # initialize chains
    assert len(theta) == len(initialize)
    p0 = []
    for i in range(nwalkers):
    	p0.append(theta + initialize*np.random.randn(len(theta)))
    
    # initialize sampler
    args = (theta0, bjd, f, ef)
    sampler = emcee.EnsembleSampler(nwalkers, len(theta), lnprob, args=args, a=a)

    # run burnin
    print 'Running burnin...'
    t0 = time.time()
    p0,_,_ = sampler.run_mcmc(p0, burnin)
    print 'Burnin acceptance fraction is %.4f'%np.mean(sampler.acceptance_fraction)
    print 'Burnin took %.4f minutes\n'%((time.time()-t0)/60.)
    sampler.reset()

    # run MCMC
    print 'Running MCMC (RVs)...'
    p0,_,_ = sampler.run_mcmc(p0, nsteps)
    samples = sampler.chain.reshape((-1, len(theta)))
    print "Mean acceptance fraction: %.4f"%np.mean(sampler.acceptance_fraction)
    print 'Full MCMC took %.4f minutes'%((time.time()-t0)/60.)

    return sampler, samples
