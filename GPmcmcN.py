from imports import *
from priors import *
from linear_lnlike import box_transit_model


def get_modelN(theta, t, f, ef):
    a, l, G, P = np.exp(theta[:4])
    s = theta[4]
    k1 = george.kernels.ExpSquaredKernel(l)
    k2 = george.kernels.ExpSine2Kernel(G,P)
    gp = george.GP(a*(k1+k2))
    try:
        gp.compute(t, np.sqrt(ef**2 + s**2))
    except (ValueError, np.linalg.LinAlgError):
        return -np.inf
    Ntransits = (theta.size-5) / 4
    assert Ntransits >= 1
    fmodel = np.zeros(t.size)
    for i in range(Ntransits):
	theta_tmp = theta[5+4*i:9+4*i]   # P,T0,depth,duration
	fmodel += box_transit_model(theta_tmp, t) - 1
    mu, cov = gp.predict(f-fmodel+1, t)
    sig = np.sqrt(np.diag(cov))
    return gp, mu, sig


def lnlike(theta, t, f, ef):
    gp,_,_ = get_modelN(theta, t, f, ef)
    lnL = gp.lnlikelihood(f, quiet=True)
    return lnL


def lnprior(theta, theta0):
    '''theta0 = Pb, T0b, Zb, Db, Pc, T0c, ...'''
    Ntransits = (theta.size-5) / 4
    lna, lnl, lnG, lnP, s = theta
    lps = np.zeros(5+4*Ntransits)
    lps[0] = lnuniform(lna, np.log(1e-5), np.log(1))
    lps[1] = lnuniform(lnl, np.log(1), np.log(1e4))
    lps[2] = lnuniform(lnG, np.log(1e-2), np.log(1e2))
    lps[3] = lnuniform(lnP, np.log(1), np.log(5e2))
    lps[4] = lnuniform(s, 0, 1) if s > 0 else -np.inf
    for i in range(Ntransits):
	P0, T00, depth0, duration0 = theta0[4*i:4+4*i]
	P, T0, depth, duration = theta[5+4*i:9+4*i]
	lps[4*i+5] = lnuniform(P, P0-.5*dutation0, P0+.5*duration0)
        lps[4*i+6] = lnuniform(T0, T00-.5*duration0, T00+.5*duration)
        lps[4*i+7] = lnuniform(depth, 0, 2*depth0)
        lps[4*i+8] = lnuniform(duration, 0, 2*duration0)
    return lps.sum()


def lnprob(theta, theta0, t, f, ef):
    lp = lnprior(theta, theta0)
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
    args = (theta, t, f, ef)
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
