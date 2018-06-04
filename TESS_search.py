from imports import *
from periodogram import compute_LSperiodogram
import GPmcmc0 as mcmc0
import GPmcmcN as mcmcN
import linear_lnlike as llnl
import rvs
from sensitivity_class import Sensitivity

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


def initialize_GP_hyperparameters(bjd, f, ef):
    '''Guess initial values of the GP hyperparameters.'''
    # a from binned LC
    tbin, fbin, efbin = boxcar(bjd, f, ef)
    lna = np.log(np.max(abs(fbin-fbin.mean())) * .75)
    # P and l from periodogram
    per,_,pwrn = compute_LSperiodogram(bjd, f, ef, plims=(1.1,bjd.max()-bjd.min()))
    ##plt.plot(per, pwrn, 'k-'), plt.xscale('log'), plt.show()
    Pgp, inds, i = 1, np.argsort(pwrn)[::-1], 0
    while np.isclose(Pgp,1,rtol=.05) or np.isclose(Pgp,2,rtol=.05):
	Pgp = per[inds[i]]
	i += 1
    lnl, lnG, lnP = np.log(Pgp), 0., np.log(Pgp)
    s = np.median(ef)*.1
    return lna, lnl, lnG, lnP, s


def do_mcmc_0(sens, bjd, f, ef, fname, Nmcmc_pnts=3e2, 
	      nwalkers=100, burnin=200, nsteps=400, a=2):
    '''First fit the PDC LC with a GP and no planet model.'''
    thetaGP = initialize_GP_hyperparameters(bjd, f, ef)
    sens.add_thetaGP(thetaGP)
    ##save_fits(thetaGP, 'Results/%s/GP_theta'%fname)
    initialize = np.array([.1,.1,.01,.1,thetaGP[4]*.1])
    assert 0 not in initialize
    Prot = np.exp(thetaGP[3])
    if (Prot/4. > (bjd.max()-bjd.min())/1e2) or (np.arange(bjd.min(), bjd.max(), Prot/4.).size > Nmcmc_pnts):
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
    sampler, samples = mcmcN.run_emcee(theta_full, tbin, fbin, efbin, initialize,
                                       nwalkers=nwalkers, burnin=burnin,
                                       nsteps=nsteps, a=a)
    results = mcmcN.get_results(samples)
    return sampler, samples, results


def save_fits(arr, fname):
    hdu = fits.PrimaryHDU(arr)
    hdu.writeto(fname, clobber=True)


def find_transits(sens, bjd, f, ef, thetaGP, hdr, fname, Nmcmc_pnts=3e2):
    '''Search for periodic transit-like events.'''
    # "detrend" the lc
    Prot = np.exp(thetaGP[3])
    if (Prot/4. > (bjd.max()-bjd.min())/1e2) or (np.arange(bjd.min(), bjd.max(), Prot/4.).size > Nmcmc_pnts):
    	dt = (bjd.max()-bjd.min())/Nmcmc_pnts
    else: 
	dt = Prot/4.
    tbin, fbin, efbin = boxcar(bjd, f, ef, dt=dt, include_edges=True)
    _, mubin, sigbin = mcmc0.get_model0(thetaGP, tbin, fbin, efbin)
    fintmu, fintsig = interp1d(tbin, mubin), interp1d(tbin, sigbin)
    mu, sig = fintmu(bjd), fintsig(bjd)
    fcorr = f - mu + 1
    sens.add_raw_timeseries(bjd, f, ef)
    sens.add_timeseries(mu, sig, fcorr)
    ##save_fits(np.array([bjd, f, ef, mu, sig, fcorr]).T, 'Results/%s/time_series'%fname)

    # do linear search first
    print 'Computing lnL over transit times and durations...\n'
    transit_times, durations, lnLs, depths = llnl.linear_search(bjd, fcorr, ef)
    sens.add_linearsearch(transit_times, durations, lnLs, depths)
    ##save_fits(transit_times, 'Results/%s/transit_times'%fname)
    ##save_fits(durations, 'Results/%s/durations'%fname)
    ##save_fits(lnLs, 'Results/%s/lnLs_linearsearch'%fname)
    ##save_fits(depths, 'Results/%s/depths'%fname)

    # get transit candidates and initial parameters guesses
    print 'Computing lnL over periods and mid-transit times...\n'
    Ps, T0s, Ds, Zs, lnLs_transit = llnl.compute_transit_lnL(bjd, fcorr, ef, transit_times, durations, lnLs, depths, 
						     	     SNRthresh)
    sens.add_paramguesses(Ps, T0s, Ds, Zs, lnLs_transit)
    ##save_fits(Ps, 'Results/%s/P_params'%fname)
    ##save_fits(T0s, 'Results/%s/T0_params'%fname)
    ##save_fits(Ds, 'Results/%s/durations_params'%fname)
    ##save_fits(Zs, 'Results/%s/depths_params'%fname)
    ##save_fits(lnLs_transit, 'Results/%s/lnLs_transit_params'%fname)

    print 'Finding transit-like events and making transit parameter guesses...\n'
    POIs, T0OIs, DOIs, ZOIs, lnLOIs, params, EBparams = llnl.identify_transit_candidates(Ps, T0s, Ds, Zs, 
									       		 lnLs_transit, 
											 durations.size,
									       		 bjd, fcorr, ef)
    sens.add_OIs(POIs, T0OIs, DOIs, ZOIs, lnLOIs)
    ##save_fits(POIs, 'Results/%s/POIs'%fname)
    ##save_fits(T0OIs, 'Results/%s/T0OIs'%fname)
    ##save_fits(DOIs, 'Results/%s/DOIs'%fname)
    ##save_fits(ZOIs, 'Results/%s/ZOIs'%fname)
    ##save_fits(lnLOIs, 'Results/%s/lnLOIs'%fname)

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
    sens.add_params_guess(params)
    sens.add_EBparams_guess(params)
    ##save_fits(params, 'Results/%s/params_guess'%fname)

    # run full mcmc with GP+transit model
    # how to do the GP with the full LC??
    #if Ntransits > 0:
 	#print 'Fitting LC with GP + %i transit models'%Ntransits
    	#sampler, samples, results = do_mcmc_N(resultsGP[0], params, bjd, f, ef)
	#save_fits(samples, 'Results/%s/full_samples'%fname_short)


if __name__ == '__main__':
    fname = sys.argv[1]
    main(fname)
