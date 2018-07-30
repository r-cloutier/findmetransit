from imports import *
from periodogram import compute_LSperiodogram
import GPmcmc0 as mcmc0
import GPmcmcN as mcmcN
import linear_lnlike as llnl
import rvs
from sensitivity_class import Sensitivity
from scipy.stats import normaltest

global SNRthresh
SNRthresh = 5

def download_one_fits(fname):
    '''Get the full fits file for a TOI.'''
    urlprefix = 'https://archive.stsci.edu/missions/tess/ete-6/tid/00/000/000/057/'
    hdu = fits.open(urlprefix+fname)
    try:
    	hdu.writeto('MAST_Data/%s'%fname)
    except IOError:
	pass
    return hdu


def get_lc(hdu):
    '''Get the light curve from a fits file.'''
    hdr = hdu[0].header
    t, f, ef = hdu[1].data['TIME'], \
	       hdu[1].data['PDCSAP_FLUX'], \
               hdu[1].data['PDCSAP_FLUX_ERR']
    g = np.isfinite(t) & np.isfinite(f) & np.isfinite(ef)
    norm = np.median(f[g])
    t, f, ef = t[g], f[g]/norm, ef[g]/norm
    return hdr, t, f, ef


def boxcar(t, f, ef, dt=.2, include_edges=False):
    '''Boxbar bin the light curve.'''
    Nbounds = int(np.floor((t.max()-t.min()) / dt))
    tbin, fbin, efbin = np.zeros(Nbounds-1), np.zeros(Nbounds-1), \
			np.zeros(Nbounds-1)
    for i in range(Nbounds-1):
  	inds = np.arange(t.size/Nbounds) + i*t.size/Nbounds
	tbin[i]  = t[inds].mean()
	fbin[i]  = np.median(f[inds])
	efbin[i] = np.mean(ef[inds]) / np.sqrt(inds.size)
    ##plt.plot(t, f, 'k.')
    ##plt.plot(tbin, fbin, 'bo')
    ##plt.show()
    if include_edges:
        tbin = np.append(np.append(t.min(), tbin), t.max())
	fbin = np.append(np.append(np.median(f[:10]), fbin), np.median(f[-10:]))
        efbin = np.append(np.append(np.median(ef[:10]), efbin), np.median(ef[-10:]))
    return tbin, fbin, efbin


def initialize_GP_hyperparameters(bjd, f, ef, Pindex=0):
    '''Guess initial values of the GP hyperparameters.'''
    # P and l from periodogram
    per,_,pwrn = compute_LSperiodogram(bjd, f, ef,
                                       plims=(1.1,bjd.max()-bjd.min()))
    ##plt.plot(per, pwrn, 'k-'), plt.xscale('log'), plt.show()
    Pgp, inds, i = 1, np.argsort(pwrn)[::-1], int(Pindex)
    # ensure Pgp is not within 5% of 1 or 2 days
    while np.isclose(Pgp,1,rtol=.05) or np.isclose(Pgp,2,rtol=.05):
	Pgp = per[inds[i]]
	i += 1
    lnl, lnG, lnP = np.log(Pgp*3), 0., np.log(Pgp)
    # s is just a fraction of the photometric uncertainty
    s = np.median(ef)*.1
    # a from binned LC
    tbin, fbin, efbin = boxcar(bjd, f, ef)
    lna = np.log(np.max(abs(fbin-fbin.mean())) * .75)
    return lna, lnl, lnG, lnP#, s


def do_optimize_0(bjd, f, ef, fname, N=10, Npnts=5e2, Nsig=3):
    '''First fit the PDC LC with GP and a no planet model using an optimization 
    routine and tested with N different initializations.'''
    # test various hyperparameter initializations and keep the one resulting
    # in the most gaussian-like residuals
    N = int(N)
    pvalues, thetaGPs_in, thetaGPs_out = np.zeros(N), np.zeros((N,4)), \
                                         np.zeros((N,4))
    for i in range(N):
        # get initial GP parameters (P, and l from periodogram)
        thetaGPs_in[i] = initialize_GP_hyperparameters(bjd, f, ef, Pindex=i)
        Prot = np.exp(thetaGPs_in[i,3])
        # bin the light curve
        if (Prot/4. > (bjd.max()-bjd.min())/1e2) or \
           (np.arange(bjd.min(), bjd.max(), Prot/4.).size > Npnts):
            dt = (bjd.max()-bjd.min())/Npnts 
        else: 
            dt = Prot/4.
	# trim outliers to avoid fitting deep transits
 	g = abs(f-np.median(f)) <= Nsig*f.std()
        tbin, fbin, efbin = boxcar(bjd[g], f[g], ef[g], dt=dt)
        gp,mu,sig,thetaGPs_out[i] = fit_GP_0(thetaGPs_in[i],
                                             tbin, fbin, efbin)
        # compute residuals and the normality test p-value
        _,pvalues[i] = normaltest(fbin-mu)
    # select the most gaussian-like residuals
    g = pvalues == pvalues.max()
    if g.sum() == 1:
        thetaGPin = thetaGPs_in[g]
        thetaGPout = thetaGPs_out[g]
    else:
        thetaGPin = thetaGPs_in[0]
        thetaGPout = thetaGPs_out[0]
    return thetaGPs_in, thetaGPs_out, thetaGPin.reshape(4), thetaGPout.reshape(4)


def fit_GP_0(thetaGP, tbin, fbin, efbin):
    '''optimize the hyperparameters of this quasi-periodic GP'''
    assert len(thetaGP) == 4
    a, l, G, Pgp = np.exp(thetaGP)
    k1 = george.kernels.ExpSquaredKernel(l)
    k2 = george.kernels.ExpSine2Kernel(G,Pgp)
    gp = george.GP(a*(k1+k2))
    bnds = ((-np.inf,np.inf),(-np.inf,np.inf),(-5,5),(-10,10))
    results = gp.optimize(tbin, fbin, efbin, **{'bounds':bnds})
    try:
        gp.compute(tbin, efbin)
    except (ValueError, np.linalg.LinAlgError):
        return np.repeat(None,4)
    mu, cov = gp.predict(fbin, tbin)
    sig = np.sqrt(np.diag(cov))
    return gp, mu, sig, results[0]
    

def OLDdo_mcmc_0(sens, bjd, f, ef, fname, Nmcmc_pnts=3e2, 
	         nwalkers=100, burnin=200, nsteps=400, a=2):
    '''First fit the PDC LC with a GP and no planet model.'''
    # get initial GP parameters (P, and l from periodogram)
    thetaGP = initialize_GP_hyperparameters(bjd, f, ef)
    sens.thetaGP = thetaGP
    initialize = np.array([.1,.1,.01,.1,thetaGP[4]*.1])
    assert 0 not in initialize
    Prot = np.exp(thetaGP[3])
    # bin the light curve
    if (Prot/4. > (bjd.max()-bjd.min())/1e2) or \
       (np.arange(bjd.min(), bjd.max(), Prot/4.).size > Nmcmc_pnts):
        dt = (bjd.max()-bjd.min())/Nmcmc_pnts 
    else: 
        dt = Prot/4.
    tbin, fbin, efbin = boxcar(bjd,f,ef,dt=dt)
    sampler, samples = mcmc0.run_emcee(thetaGP, tbin, fbin, efbin, initialize,
				       nwalkers=nwalkers, burnin=burnin,
				       nsteps=nsteps, a=a)
    if samples.shape[0] > 1e4:
    	inds = np.arange(samples.shape[0])
    	np.random.shuffle(inds)
    	results = mcmc0.get_results(samples[inds][:int(1e4)])
    else:
	results = mcmc0.get_results(samples)
    return sampler, samples, results


def do_mcmc_N(thetaGP, params, bjd, f, ef, Nmcmc_pnts=3e2,
              nwalkers=100, burnin=200, nsteps=400, a=2):
    '''Fit the LC with a GP and a full transit model.'''
    Ntransits = params.size / 4
    theta = params.reshape(Ntransits*4)
    assert thetaGP.size == 5
    initialize = np.array([.1,.1,.01,.1,thetaGP[4]*.1])
    for i in range(Ntransits):
	initialize = np.append(initialize, [.1,.1,params[i,2]*.1,params[i,3]*.1])
    assert 0 not in initialize
    Prot = np.exp(thetaGP[3])
    if (Prot/4. > (bjd.max()-bjd.min())/1e2) or (np.arange(bjd.min(), bjd.max(), Prot/4.).size > Nmcmc_pnts):
        dt = (bjd.max()-bjd.min())/Nmcmc_pnts 
    else: 
        dt = Prot/4.
    tbin, fbin, efbin = boxcar(bjd,f,ef,dt=dt)
    theta_full = np.append(thetaGP, theta)
    sampler, samples = mcmcN.run_emcee(theta_full, tbin, fbin, efbin,
                                       initialize, nwalkers=nwalkers,
                                       burnin=burnin, nsteps=nsteps, a=a)
    results = mcmcN.get_results(samples)
    return sampler, samples, results


def save_fits(arr, fname):
    hdu = fits.PrimaryHDU(arr)
    hdu.writeto(fname, clobber=True)


def _optimize_GP(thetaGP, x, res, ey):
    '''Optimize the GP parameters of a mean-subtracted time-series and return the 
    mean GP model.'''
    assert len(thetaGP) == 4
    a, l, G, Pgp = np.exp(thetaGP)
    k1 = george.kernels.ExpSquaredKernel(l)
    k2 = george.kernels.ExpSine2Kernel(G,Pgp)
    gp = george.GP(a*(k1+k2))
    try:
	results = gp.optimize(x, res, ey)
	results = results[0]
    	gp.compute(x, ey)
    	mu, cov = gp.predict(res, x)
    	sig = np.sqrt(np.diag(cov))
    except (ValueError, np.linalg.LinAlgError):
	gp, results, mu, sig = None, np.zeros(len(thetaGP)), np.zeros(x.size), np.zeros(x.size)
    return gp, results, mu, sig


def find_transits(sens, bjd, f, ef, thetaGP, hdr, fname, Npnts=5e2):
    '''Search for periodic transit-like events.'''
    # "detrend" the lc
    assert len(thetaGP) == 4
    Prot = np.exp(thetaGP[3])
    if (Prot/4. > (bjd.max()-bjd.min())/1e2) or (np.arange(bjd.min(), bjd.max(), Prot/4.).size > Npnts):
    	dt = (bjd.max()-bjd.min())/Npnts
    else: 
	dt = Prot/4.
    tbin, fbin, efbin = boxcar(bjd, f, ef, dt=dt, include_edges=True)
    #_, mubin, sigbin = mcmc0.get_model0(thetaGP, tbin, fbin, efbin)
    _, resultsGP, mubin, sigbin = _optimize_GP(thetaGP, tbin, fbin, efbin)
    fintmu, fintsig = interp1d(tbin, mubin), interp1d(tbin, sigbin)
    mu, sig = fintmu(bjd), fintsig(bjd)
    fcorr = f - mu + 1
    sens.bjd, sens.f, sens.ef = bjd, f, ef
    sens.mu, sens.sig, sens.fcorr = mu, sig, fcorr
    sens.resultsGP_detrend = resultsGP

    # do linear search first
    print 'Computing lnL over transit times and durations...\n'
    transit_times, durations, lnLs, depths = llnl.linear_search(bjd, fcorr, ef)
    sens.transit_times, sens.durations = transit_times, durations
    sens.lnLs_linearsearch, sens.depths_linearsearch = lnLs, depths
    sens.SNRs_linearsearch = (lnLs-np.median(lnLs, axis=0)) / llnl.MAD2d(lnLs)

    # get transit candidates and initial parameters guesses
    print 'Computing lnL over periods and mid-transit times...\n'
    Ps, T0s, Ds, Zs, lnLs_transit = llnl.compute_transit_lnL(bjd, fcorr, ef,
                                                             transit_times,
                                                             durations, lnLs,
                                                             depths, SNRthresh)

    print 'Finding transit-like events and making transit parameter guesses...\n'
    Rs = hdr['RADIUS']
    POIs, T0OIs, DOIs, ZOIs, lnLOIs, params, EBparams = llnl.identify_transit_candidates(Ps, T0s, Ds, Zs, 
									       		 lnLs_transit, 
											 durations.size, Rs,
									       		 bjd, fcorr, ef)
    sens.POIs, sens.T0OIs, sens.DOIs, sens.ZOIs, sens.lnLOIs = POIs, T0OIs, DOIs, ZOIs, lnLOIs

    # return the parameters (P,T0,D,Z=depth) of the confirmed transit candidates
    return params, EBparams
 

def estimate_box_transit_model(P, T0, Rs, t, f, ef):
    '''Estimate the transit depth and duration given P and T0. Return 
    the box transit model.'''
    phase = foldAt(t, P, T0)
    phase[phase>.5] -= 1
    intransit_approx = (phase*P <= 15./60/24) & (phase*P >= -15./60/24)
    depth = np.median(f[intransit_approx])
    duration = rvs.transit_width(P, Rs, Rs, 
				 rvs.m2Rearth(np.sqrt(depth)*rvs.Rsun2m(Rs)), 
				 0)
    model = llnl.box_transit_model((P,T0,depth,duration), t)
    return model


def is_good_star(hdr):
    if hdr['TESSMAG'] > 12:
	return False
    elif hdr['TEFF'] >= 4200:
	return False
    else:
	return True


def main(fname):
    # get fits file
    print 'Downloading fits file...\n'
    hdu = download_one_fits(fname)

    # get LC
    hdr, bjd, f, ef = get_lc(hdu)

    # only continue if it's a bright M-dwarf
    #if not is_good_star(hdr):  TEMP
#	raise ValueError('Not a star of interest.')

    # fit systematics with a GP
    print 'Fitting LC with GP alone...\n'
    fname_short = fname.replace('.fits','')
    sens = Sensitivity(fname_short)
    samplerGP, samplesGP, resultsGP = do_mcmc_0(sens, bjd, f, ef, fname_short)
    sens.add_GPsamples(samplesGP)
    sens.add_GPresults(resultsGP)
    ##save_fits(samplesGP, 'Results/%s/GP_samples'%fname_short)
    ##save_fits(resultsGP, 'Results/%s/GP_results'%fname_short)

    # search for transits in the corrected LC and get the transit parameters guesses
    print 'Searching for transit-like events...\n'
    params, EBparams = find_transits(sens, bjd, f, ef, resultsGP[0], hdr, fname_short)
    sens.params_guess = params
    sens.EBparams_guess = EBparams

    # run full mcmc with GP+transit model
    # how to do the GP with the full LC??
    #if Ntransits > 0:
 	#print 'Fitting LC with GP + %i transit models'%Ntransits
    	#sampler, samples, results = do_mcmc_N(resultsGP[0], params, bjd, f, ef)
	#save_fits(samples, 'Results/%s/full_samples'%fname_short)


if __name__ == '__main__':
    fname = sys.argv[1]
    main(fname)
