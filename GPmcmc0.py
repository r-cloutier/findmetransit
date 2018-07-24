from imports import *
from priors import *

def get_model0(theta, t, f, ef):
    a, l, G, P = np.exp(theta[:4])
    k1 = george.kernels.ExpSquaredKernel(l)
    k2 = george.kernels.ExpSine2Kernel(G,P)
    gp = george.GP(a*(k1+k2))
    try:
        gp.compute(t, ef)
    except (ValueError, np.linalg.LinAlgError):
        return -np.inf
    mu, cov = gp.predict(f, t)
    sig = np.sqrt(np.diag(cov))
    return gp, mu, sig


def lnlike(theta, t, f, ef):
    gp,_,_ = get_model0(theta, t, f, ef)
    lnL = gp.lnlikelihood(f, quiet=True)
    return lnL


def lnprior(theta):
    lna, lnl, lnG, lnP = theta
    lps = np.zeros(4)
    lps[0] = lnuniform(lna, np.log(1e-5), np.log(1))
    lps[1] = lnuniform(lnl, np.log(1e-1), np.log(1e8))
    lps[2] = lnuniform(lnG, np.log(1e-3), np.log(1e3))
    lps[3] = lnuniform(lnP, np.log(1e-1), np.log(5e2))
    return lps.sum()


def lnprob(theta, t, f, ef):
    lp = lnprior(theta)
    if np.isfinite(lp):
        return lp + lnlike(theta, t, f, ef)
    else:
        return -np.inf


def run_emcee(theta, t, f, ef, initialize, nwalkers=100, burnin=200,
              nsteps=400, a=2):
    '''Run mcmc on an input light curve with no transit model.'''
    # initialize chains
    assert len(theta) == len(initialize)
    p0 = []
    for i in range(nwalkers):
    	p0.append(theta + initialize*np.random.randn(len(theta)))
    
    # initialize sampler
    args = (t, f, ef)
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
