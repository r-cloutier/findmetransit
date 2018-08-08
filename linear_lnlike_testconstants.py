from imports import *
from scipy.optimize import curve_fit
import vetting as vett
import rvs, mcmc1, batman
from scipy.interpolate import LinearNDInterpolator as lint

def lnlike(bjd, f, ef, fmodel):
    return -.5*(np.sum((f-fmodel)**2 / ef**2 - np.log(1./ef**2)))


def box_transit_model(theta, t):
    '''Return the box transmit model.'''
    fmodel = np.ones(t.size)
    P, T0, depth, duration = theta
    phase = foldAt(t, P, T0)
    phase[phase>.5] -= 1
    intransit = (phase*P <= .5*duration) & (phase*P >= -.5*duration)
    fmodel[intransit] = 1. - depth
    return fmodel


def box_transit_model_curve_fit(t, P, T0, depth, duration):
    '''Return the box transmit model.'''
    fmodel = np.ones(t.size)
    phase = foldAt(t, P, T0)
    phase[phase>.5] -= 1
    intransit = (phase*P <= .5*duration) & (phase*P >= -.5*duration)
    fmodel[intransit] = 1. - depth
    return fmodel


def box_transit_model_time(theta, t):
    '''Return the box transit model as a function of transit time 
    instead of P and T0.'''
    fmodel = np.ones(t.size)
    T, depth, duration = theta 
    intransit = (t >= T-.5*duration) & (t <= T+.5*duration)
    fmodel[intransit] = 1. - depth
    return fmodel


def get_depth_lnlike(theta, bjd, fcorr, ef, depthmax=.1, N=2e2):
    '''Given a transit time and duration, get the max lnL depth and its lnL.'''
    # sample depths
    depthmin, N = np.median(ef), int(N)
    depths = 10**(np.random.uniform(np.log10(depthmin), np.log10(depthmax), N))

    # compute lnLs of the depths given a transit time and duration
    T, duration = theta
    lnLs = np.zeros(N)
    for i in range(N):
    	fmodel = box_transit_model_time((T,depths[i],duration), bjd)
        lnLs[i] = lnlike(bjd, fcorr, ef, fmodel)

    # get maximum-likelihood depth
    g = lnLs == lnLs.max()
    if g.sum() == 1:
	lnL, D = float(lnLs[g]), float(depths[g])
    elif g.sum() == 0:
	lnL, D = np.nan, np.nan
    else:
	g = np.where(g)[0][0]
	lnL, D = lnLs[g], depths[g]
    return lnLs, depths, lnL, D


def linear_search(bjd, fcorr, ef):
    '''Evaluate the lnL as a function of transit duration and 
    transit time/epoch of a transit.'''
    # setup transit time and duration grids
    transit_times = np.arange(bjd.min(), bjd.max(), 30./60/24)
    durations = np.array([1.2,2.4,4.8]) / 24  # coarse grid in transit duration

    # get max lnL depth over the linear grid of transit times and durations
    NT, ND = transit_times.size, durations.size
    lnLs, depths = np.zeros((NT,ND)), np.zeros((NT,ND))
    for i in range(NT):
	print float(i) / NT
	for j in range(ND):
	    theta = transit_times[i], durations[j]
	    _,_,lnLs[i,j],depths[i,j] = get_depth_lnlike(theta, bjd, fcorr, ef)

    return transit_times, durations, lnLs, depths


def MAD1d(arr):
    assert len(arr.shape) == 1
    return np.median(abs(arr-np.median(arr)))

def MAD2d(arr):
    assert len(arr.shape) == 2
    mads = np.zeros(arr.shape[1])
    for i in range(mads.size):
	mads[i] = np.median(abs(arr[:,i] - np.median(arr[:,i])))
    return mads


def find_transit_parameters(bjd, fcorr, ef, 
			    transit_times, durations, lnLs, depths, 
			    SNRthresh):
    '''Find periodic transits and return initial guesses of P, T0, 
    duration, and depth.'''
    # find high S/N peaks in lnL as a function of the transit time
    SNRs = (lnLs-np.median(lnLs, axis=0)) / MAD2d(lnLs)
    gSNR = SNRs > SNRthresh

    # search over coarse duration grid
    NT, ND = transit_times.size, durations.size
    Ps_full, T0s_full, durations_full, depths_full, lnLs_full = np.zeros(0), \
						     		np.zeros(0), \
                                                     		np.zeros(0), \
                                                     		np.zeros(0), \
								np.zeros(0)
    for i in range(ND):

    	if gSNR[:,i].sum() == 0:  # no transit-like events
	    pass

   	else:  # some potential transit-like events
	    
	    # get unique approximate transit times of transit-like events
	    Ts_all = transit_times[gSNR[:,i]]
	    Ts_reduced = np.array([Ts_all[0]])
	    for j in range(1,Ts_all.size):
		if np.isclose(Ts_all[j]-Ts_reduced[-1], 0, atol=durations[i]*2):  # same transit
		    pass
		else:
		    Ts_reduced = np.append(Ts_reduced, Ts_all[j])

 	    # adjust to more accurate transit times
	    Ntransits = Ts_reduced.size
	    T0s = np.zeros(0)
	    for j in range(Ntransits):
		g = (bjd >= Ts_reduced[j]-2*durations[i]) & \
		    (bjd <= Ts_reduced[j]+2*durations[i])
		if g.sum() > 0:
		    fsmooth = gaussian_filter1d(fcorr[g], 5)
		    T0s = np.append(T0s, bjd[g][fsmooth == fsmooth.min()])
		    ##plt.plot(bjd, fcorr, '-', bjd[g], fcorr[g], 'o')
		    ##plt.plot(bjd[g], fsmooth, '-', lw=2)
		    ##plt.axvline(Ts_reduced[j]), plt.axvline(T0s[j], ls='--')
		    ##plt.show()
	    T0s = np.unique(T0s)

	    # search for all periodic transits (periods between transit events)
	    Ntransits = T0s.size
	    for j in range(Ntransits):
		for k in range(Ntransits):
		    if j != k:
			Ps_full = np.append(Ps_full, T0s[j]-T0s[k])
		      	T0s_full = np.append(T0s_full, T0s[j])
			durations_full = np.append(durations_full, durations[i])
		    	phase = foldAt(bjd, Ps_full[-1], T0s_full[-1])
		    	phase[phase > .5] -= 1
		    	intransit = (phase*Ps_full[-1] <= .5*durations_full[-1]) & \
				    (phase*Ps_full[-1] >= -.5*durations_full[-1])
                        depths_full = np.append(depths_full, 1-np.median(fcorr[intransit]))  # could be negative
			theta = Ps_full[-1], T0s_full[-1], depths_full[-1], durations_full[-1]
			fmodel = box_transit_model(theta, bjd)
			lnLs_full = np.append(lnLs_full, lnlike(bjd, fcorr, ef, fmodel))
			##plt.plot(phase, fcorr, '-'), plt.plot(phase[intransit], fcorr[intransit], 'o'), plt.show()

    # trim (must be greater than ~.5 and smaller than the baseline)
    g = (Ps_full >= .49) & (Ps_full <= (bjd.max()-bjd.min()))

    return Ps_full[g], T0s_full[g], durations_full[g], depths_full[g], lnLs_full[g]


def compute_P_lnL(bjd, fcorr, ef, theta, N=1e2):
    '''Compute lnL for each P marginalized over T0 then return the maxmimum 
    likelihood result.'''
    N = int(N)
    P, T00, Z, D = theta
    lnLs, thetas = np.zeros(N), np.zeros((N,4))
    for i in range(N):
	T0 = np.random.uniform(T00-.5*P,T00+.5*P)
	phase = foldAt(bjd, P, T0)
	phase[phase > .5] -= 1
	intransit = (phase*P <= .5*D) & (phase*P >= -.5*D)
	thetas[i] = P, T0, 1-np.median(fcorr[intransit]), D
	fmodel = box_transit_model(thetas[i], bjd)
	lnLs[i] = lnlike(bjd, fcorr, ef, fmodel)
    # get maximum lnL result
    g = lnLs == lnLs.max()
    return lnLs[g][0], thetas[g][0]


def compute_transit_lnL(bjd, fcorr, ef, transit_times, durations, lnLs, depths, SNRthresh):
    '''Get transit parameters and compute the lnL to identify transits.'''
    # get transit parameters
    assert lnLs.shape == (transit_times.size, durations.size)
    assert depths.shape == (transit_times.size, durations.size)
    Ps, T0s, Ds, Zs, lnLs_transit = find_transit_parameters(bjd, fcorr, ef,
							    transit_times,
                                                            durations, 
							    lnLs, depths,
                                                            SNRthresh)
    assert Ps.size == T0s.size
    assert Ps.size == Ds.size
    assert Ps.size == Zs.size
    assert Ps.size == lnLs_transit.size

    # compute lnL
    '''lnLs = np.zeros(Ps.size)
    for i in range(lnLs.size):
	print float(i)/Ps.size
 	theta = Ps[i], T0s[i], Zs[i], Ds[i]
	fmodel = box_transit_model(theta, bjd)
	##plt.plot(bjd, fcorr, 'o', bjd, fmodel, '-'), plt.show()
	#lnLs[i] = lnlike(bjd, fcorr, ef, fmodel)
	lnLs[i], thetaout = compute_P_lnL(bjd, fcorr, ef, theta)
	print thetaout
	Ps[i], T0s[i], Zs[i], Ds[i] = thetaout'''

    return Ps, T0s, Ds, Zs, lnLs_transit


def remove_multiple_on_lnLs(bjd, Ps, T0s, Ds, Zs, lnLs, dP=.1):
    '''remove multiple orbital periods but dont assume the shortest one is
    correct, instead select the one with the highest lnL.'''
    assert Ps.size == T0s.size
    assert Ps.size == Ds.size
    assert Ps.size == Zs.size
    assert Ps.size == lnLs.size
    to_remove = np.zeros(0)
    for i in range(Ps.size):
	Ntransits = int((bjd.max()-bjd.min()) / Ps[i])
	lim = Ntransits+1 if Ntransits+1 > 2 else 3
	for j in range(2,lim):
	    # check positive multiples
	    isclose = np.isclose(Ps, Ps[i]*j, atol=dP*2)
	    if np.any(isclose):
		iscloselnL = lnLs[isclose] <= lnLs[i]
		to_remove = np.append(to_remove, Ps[isclose][iscloselnL])
	    # check inverse multiples
	    isclose = np.isclose(Ps, Ps[i]/float(j), atol=dP*2)
	    if np.any(isclose):
		iscloselnL = lnLs[isclose] <= lnLs[i]
		to_remove = np.append(to_remove, Ps[isclose][iscloselnL])
    to_remove = np.unique(to_remove)
    assert to_remove.size <= Ps.size
    to_remove_inds = np.where(np.in1d(Ps, to_remove))[0]
    Ps_final = np.delete(Ps, to_remove_inds)
    T0s_final = np.delete(T0s, to_remove_inds)
    Ds_final = np.delete(Ds, to_remove_inds)
    Zs_final = np.delete(Zs, to_remove_inds)
    lnLs_final = np.delete(lnLs, to_remove_inds)
    return Ps_final, T0s_final, Ds_final, Zs_final, lnLs_final


def remove_multiples(bjd, Ps, T0s, Ds, Zs, lnLs, dP=.1):
    assert Ps.size == T0s.size
    assert Ps.size == Ds.size
    assert Ps.size == Zs.size
    assert Ps.size == lnLs.size
    to_remove = np.zeros(0)
    for i in range(Ps.size):
        if Ps[i] not in to_remove:
            Ntransits = int((bjd.max()-bjd.min()) / Ps[i])
            for j in range(2,Ntransits+1):
                isclose = np.isclose(Ps, Ps[i]*j, atol=dP*2)
                if np.any(isclose):
                    to_remove = np.append(to_remove, Ps[isclose])
    to_remove = np.unique(to_remove)
    assert to_remove.size <= Ps.size
    to_remove_inds = np.where(np.in1d(Ps, to_remove))[0]
    Ps_final = np.delete(Ps, to_remove_inds)
    T0s_final = np.delete(T0s, to_remove_inds)
    Ds_final = np.delete(Ds, to_remove_inds)
    Zs_final = np.delete(Zs, to_remove_inds)
    lnLs_final = np.delete(lnLs, to_remove_inds)
    return Ps_final, T0s_final, Ds_final, Zs_final, lnLs_final


def _get_LDcoeffs(Ms, Rs, Teff, Z=0):
    '''Interpolate Claret 2017 grid of limb darkening coefficients to a
    given star.'''
    # get LD coefficient grid (Z is always 0 for some reason)
    clarlogg, clarTeff, clarZ, clar_a, clar_b = \
                                    np.loadtxt('LDcoeffs/claret17.tsv',
                                               delimiter=';', skiprows=37).T

    # interpolate to get the stellar LD coefficients
    logg = np.log10(6.67e-11*rvs.Msun2kg(Ms)*1e2 / rvs.Rsun2m(Rs)**2)
    lint_a = lint(np.array([clarTeff,clarlogg]).T, clar_a)
    lint_b = lint(np.array([clarTeff,clarlogg]).T, clar_b)

    return float(lint_a(Teff,logg)), float(lint_b(Teff,logg))


def transit_model_func_curve_fit(u1, u2):
    def transit_model_func_in(bjd, P, T0, aRs, rpRs, inc):
        p = batman.TransitParams()
        p.t0, p.per, p.rp = 0., 1., rpRs
        p.a, p.inc, p.ecc = aRs, inc, 0.
        p.w, p.limb_dark, p.u = 90., 'quadratic', [u1,u2]
        phase = foldAt(bjd, P, T0)
        m = batman.TransitModel(p, phase)
        f = m.light_curve(p)
        return f
    return transit_model_func_in


def _fit_params(params, bjd, fcorr, ef, Ms, Rs, Teff):
    '''Get best-fit parameters.'''
    assert params.shape == (4,)
    P, T0, depth, duration = params
    u1, u2 = _get_LDcoeffs(Ms, Rs, Teff)
    aRs = rvs.AU2m(rvs.semimajoraxis(P,Ms,0)) / rvs.Rsun2m(Rs)
    rpRs = np.sqrt(depth)
    p0 = P, T0, aRs, rpRs, 90.
    incs = np.array([float(rvs.inclination(P,Ms,Rs,1)), float(rvs.inclination(P,Ms,Rs,-1))])
    bnds = ((P*.9, T0-P*1.1, aRs*.9, 0, incs.min()),
            (P*1.1, T0+P*1.1, aRs*1.1, 1, incs.max()))
    try:
        popt,_ = curve_fit(transit_model_func_curve_fit(u1,u2),
                           bjd, fcorr, p0=p0, sigma=ef,
                           absolute_sigma=True, bounds=bnds)
        P, T0, aRs, rpRs, inc = popt
        depth = rpRs**2
        b = rvs.impactparam_inc(P, Ms, Rs, inc)
        duration = P/(np.pi*aRs) * np.sqrt((1+np.sqrt(depth))**2 - b*b)
	func = transit_model_func_curve_fit(u1, u2)
	fmodel = func(bjd, P, T0, aRs, rpRs, inc)
    except RuntimeError:
	func = transit_model_func_curve_fit(u1, u2)
	fmodel = func(bjd, P, T0, aRs, rpRs, 90.)
    return P, T0, depth, duration, fmodel


def trim_planets(params, lnLOIs, Nplanetsmax=5):
    '''If there are too many planet candidates then remove some based on 
    their lnLs. Most are FPs.'''
    Nplanets = params.shape[0]
    assert params.shape == (Nplanets,4)
    assert lnLOIs.size == Nplanets
    if Nplanets > Nplanetsmax:
	tokeep = np.sort(np.argsort(lnLOIs)[::-1][:Nplanetsmax])  # retain the same ordering
	return params[tokeep], lnLOIs[tokeep]
    else:
	return params, lnLOIs


def confirm_transits(params, lnLs, bjd, fcorr, ef, Ms, Rs, Teff, dispersion_sig, depth_sig, bimodalfrac):
    '''Look at proposed transits and confirm whether or not a significant 
    dimming is seen.'''
    Ntransits = params.shape[0]
    assert lnLs.size == Ntransits
    paramsout, to_remove_inds = np.zeros((Ntransits,4)), np.zeros(0)
    transit_condition_scatterin_gtr_scatterout = np.zeros(Ntransits, dtype=bool)
    transit_condition_depth_gtr_rms = np.zeros(Ntransits, dtype=bool)
    transit_condition_no_bimodal_flux_intransit = np.zeros(Ntransits, dtype=bool)
    transit_condition_orbitalP_fits_in_WF = np.zeros(Ntransits, dtype=bool)
    print 'Confirming proposed transits...'
    for i in range(Ntransits):
	print float(i) / Ntransits
	# run mcmc to get best parameters for the proposed transit
	#initialize = np.array([params[i,3],params[i,3],.1*params[i,2],
        #                       .1*params[i,3]])
  	#sampler, samples = mcmc1.run_emcee(params[i], params[i], 
	#				   bjd, fcorr, ef, initialize, a=1.9)
	#results = mcmc1.get_results(samples)
	# get optimized parameters for this transit
	P, T0, depth, duration,_ = _fit_params(params[i], bjd, fcorr, ef, Ms, Rs, Teff)
	paramsout[i] = P, T0, depth, duration

	# get in and out of transit window
        phase = foldAt(bjd, P, T0)
        phase[phase > .5] -= 1
	Dfrac = .25   # fraction of the duration in-transit (should be <.5 to ignore ingress & egress)
        intransit = (phase*P >= -Dfrac*duration) & (phase*P <= Dfrac*duration)
	intransitfull = (phase*P >= -duration/2) & (phase*P <= duration/2)
	outtransit = (phase*P <= -(1.+Dfrac)*duration) | (phase*P >= (1.+Dfrac)*duration)
        #plt.plot(phase, fcorr, 'ko', phase[intransit], fcorr[intransit], 'bo'), plt.show()

        # check scatter in and out of the proposed transit to see if the transit is real
	#cond1 = np.median(fcorr[intransit]) <= np.median(fcorr[outtransit]) - dispersion_sig*MAD1d(fcorr[outtransit])
	cond1 = (np.median(fcorr[outtransit]) - np.median(fcorr[intransit])) / MAD1d(fcorr[outtransit]) > dispersion_sig
	transit_condition_scatterin_gtr_scatterout[i] = cond1
	# also check that the transit depth is significant relative to the noise
	depth = 1-np.median(fcorr[intransit])
	sigdepth = np.median(ef[intransit])
	cond2 = depth/sigdepth > depth_sig
	transit_condition_depth_gtr_rms[i] = cond2
	# ensure that the flux measurements intransit are not bimodal (ie. at depth and at f=1 which would indicate a 
	# bad period and hence a FP
	y, x = np.histogram(fcorr[intransitfull], bins=30)
        x = x[1:] - np.diff(x)[0]/2.
	cond3 = float(y[x<x.mean()].sum())/y.sum() > bimodalfrac
	transit_condition_no_bimodal_flux_intransit[i] = cond3
        # ensure that at least two transits will fit within the observing window otherwise its just a
        # single transit-like event
        cond4 = (T0-P >= bjd.min()) | (T0+P <= bjd.max())
        transit_condition_orbitalP_fits_in_WF[i] = cond4
	if cond1 and cond2 and cond3 and cond4:
	    pass
	else:
	    to_remove_inds = np.append(to_remove_inds, i)

    # remove false transits
    paramsout = np.delete(paramsout, to_remove_inds, 0)
    lnLsout = np.delete(lnLs, to_remove_inds)

    return paramsout, lnLsout, transit_condition_scatterin_gtr_scatterout, transit_condition_depth_gtr_rms, transit_condition_no_bimodal_flux_intransit, transit_condition_orbitalP_fits_in_WF


def identify_EBs(params, bjd, fcorr, ef, Rs, SNRthresh=3., rpmax=30):
    '''For each proposed planet in params, check if there is a clearly-defined 
    secondary eclipse as is indicative of a secondary eclipse. Could also have a 
    V-shaped "transit" but so highly inclined transiting planets.'''
    Nplanets = params.shape[0]
    notEB = np.ones(Nplanets)
    for i in range(Nplanets):
	# check for significant eclipse depths at two times in case T0 is near an edge of the WF
	P, T0, depth, duration = params[i]
	eclipse1 = (bjd >= T0+.5*P-.5*duration) & (bjd <= T0+.5*P+.5*duration)
	eclipse2 = (bjd >= T0-.5*P-.5*duration) & (bjd <= T0-.5*P+.5*duration)	
	outeclipse1 = (bjd >= T0+.5*P-2*duration) & (bjd <= T0+.5*P-duration)
        outeclipse2 = (bjd >= T0-.5*P-2*duration) & (bjd <= T0-.5*P-duration)

	rms_ineclipse1  = fcorr[eclipse1].std()
	rms_outeclipse1 = fcorr[outeclipse1].std()
        rms_ineclipse2  = fcorr[eclipse2].std()
        rms_outeclipse2 = fcorr[outeclipse2].std()
        if (rms_ineclipse1 >= SNRthresh*rms_outeclipse1) or (rms_ineclipse2 >= SNRthresh*rms_outeclipse2):
            notEB[i] = 0.

    # check that planets are not too big
    rpRss = np.sqrt(params[:,2])
    rps = rvs.m2Rearth(rvs.Rsun2m(rpRss*Rs))
    notEB[(rpRss > .5) | (rps > rpmax)] = 0

    # save planet and EB parameters
    notEB = notEB.astype(bool)
    params, EBparams = params[notEB], params[np.invert(notEB)]
    return params, EBparams 
