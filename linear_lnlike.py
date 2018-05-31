from imports import *
import rvs, mcmc1

global dispersion_sig
dispersion_sig = 2

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
    transit time or epoch of a transit.'''
    # setup transit time and duration grids
    transit_times = np.arange(bjd.min(), bjd.max(), 30./60/24)
    durations = np.array([1.2,2.4,4.8]) / 24   # coarse grid in transit duration

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
    #SNRs = (lnLs-np.median(lnLs, axis=0)) / np.std(lnLs, axis=0)
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
	    ##T0s = np.repeat(T0s,Ntransits).reshape(Ntransits,Ntransits)
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
                        depths_full = np.append(depths_full, 1-np.median(fcorr[intransit]))
			theta = Ps_full[-1], T0s_full[-1], depths_full[-1], durations_full[-1]
			fmodel = box_transit_model(theta, bjd)
			lnLs_full = np.append(lnLs_full, lnlike(bjd, fcorr, ef, fmodel))
			##plt.plot(phase, fcorr, '-'), plt.plot(phase[intransit], fcorr[intransit], 'o'), plt.show()

    # trim
    g = (Ps_full >= .49) & (Ps_full <= (bjd.max()-bjd.min())/3.)

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
							    transit_times, durations, 
							    lnLs, depths, SNRthresh)
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
    to_remove = np.zeros(0)
    for i in range(Ps.size):
	Ntransits = int((bjd.max()-bjd.min()) / Ps[i])
	for j in range(2,Ntransits+1):
	    isclose = np.isclose(Ps, Ps[i]*j, atol=dP*2)
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


def identify_transit_candidates(Ps, T0s, Ds, Zs, lnLs, Ndurations, bjd, fcorr, ef):
    '''Given the transit parameters and their lnLs, identify transit 
    candidates.'''
    assert Ps.size == T0s.size
    assert Ps.size == Ds.size
    assert Ps.size == Zs.size
    assert Ps.size == lnLs.size

    # remove common periods based on maximum likelihood 
    dP = .1
    sort = np.argsort(Ps)
    POIs, T0OIs, DOIs, ZOIs, lnLOIs = Ps[sort], T0s[sort], Ds[sort], Zs[sort], lnLs[sort]
    POIs_red, T0OIs_red, DOIs_red, ZOIs_red, lnLOIs_red = np.zeros(0), np.zeros(0), \
							  np.zeros(0), np.zeros(0), \
							  np.zeros(0)
    for i in range(POIs.size):
        isclose = np.isclose(POIs, POIs[i], atol=dP*2)
        if np.any(isclose):
            g = lnLOIs == lnLOIs[isclose].max()
            POIs_red = np.append(POIs_red, POIs[g])
            T0OIs_red = np.append(T0OIs_red, T0OIs[g])
            DOIs_red = np.append(DOIs_red, DOIs[g])
            ZOIs_red = np.append(ZOIs_red, ZOIs[g])
	    lnLOIs_red = np.append(lnLOIs_red, lnLOIs[g])
    _,unique = np.unique(POIs_red, return_index=True)
    POIs_red, T0OIs_red, DOIs_red, ZOIs_red, lnLOIs_red = POIs_red[unique], \
							  T0OIs_red[unique], \
                                                          DOIs_red[unique], \
                                                          ZOIs_red[unique], \
                                                          lnLOIs_red[unique]

    # remove multiple transits (i.e. 2P, 3P, 4P...)
    POIs_final, T0OIs_final, DOIs_final, ZOIs_final, lnLOIs_final = \
  	remove_multiple_on_lnLs(bjd, POIs_red, T0OIs_red, DOIs_red, ZOIs_red,
			 	lnLOIs_red)

    # get initial parameter guess for identified transits
    params = np.array([POIs_final, T0OIs_final, ZOIs_final, DOIs_final]).T

    # remove duplicates
    params = params[np.unique(params[:,0], return_index=True)[1]]

    # identify bona-fide transit-like events
    params = confirm_transits(params, bjd, fcorr, ef)

    # re-remove multiple transits based on refined parameters
    p,t0,d,z,_ = remove_multiples(bjd, params[:,0], params[:,1], params[:,3], 
				  params[:,2], np.zeros(params[:,0].size))
    params = np.array([p,t0,z,d]).T

    return POIs_final, T0OIs_final, DOIs_final, ZOIs_final, lnLOIs_final, params


def confirm_transits(params, bjd, fcorr, ef):
    '''Look at proposed transits and confirm whether or not a significant dimming 
    is seen.'''
    Ntransits = params.shape[0]
    paramsout, to_remove_inds = np.zeros((Ntransits,4)), np.zeros(0)
    print 'Confirming proposed transits...'
    for i in range(Ntransits):
	print float(i) / Ntransits
	# run mcmc to get best parameters for the proposed transit
	initialize = np.array([params[i,3],params[i,3],.1*params[i,2],.1*params[i,3]])
  	sampler, samples = mcmc1.run_emcee(params[i], params[i], 
					   bjd, fcorr, ef, initialize, a=1.9)
	results = mcmc1.get_results(samples)
	P, T0, depth, duration = results[0]
	paramsout[i] = P, T0, depth, duration
        phase = foldAt(bjd, P, T0)
        phase[phase > .5] -= 1
	Dfrac = .3   # fraction of the duration in-transit (should be <.5 to ignore ingress & egress)
        intransit = (phase*P >= -Dfrac*duration) & (phase*P <= Dfrac*duration)
	outtransit = (phase*P <= -(1.+Dfrac)*duration) | (phase*P >= (1.+Dfrac)*duration)
        ##plt.plot(phase, fcorr, 'ko', phase[intransit], fcorr[intransit], 'bo'), plt.show()
        # check scatter in and out of the proposed transit to see if the transit is real
	cond1 = np.median(fcorr[intransit]) <= np.median(fcorr[outtransit]) - dispersion_sig*MAD1d(fcorr[outtransit])
	sigdepth = np.mean(results[:,2])
  	cond2 = (1-np.median(fcorr[intransit]) <= depth+sigdepth) & \
	   	(1-np.median(fcorr[intransit]) >= depth-sigdepth)
	if cond1 and cond2:
	    pass
	else:
	    to_remove_inds = np.append(to_remove_inds, i)

    # remove false transits
    paramsout = np.delete(paramsout, to_remove_inds, 0)

    return paramsout