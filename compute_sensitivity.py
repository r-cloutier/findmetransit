from TESS_search import *
import matplotlib.gridspec as gridspec
import linear_lnlike as llnl
import batman, ellc
from massradius import radF2mass
from truncate_cmap import *
from sensitivity_class import *
from joint_LCmodel import *


# M dwarfs from the TESS CTL (candidate target list) only 
global medTmag, medTeff, medRs
medTmag, medTeff, medRs = 10.5, 3400, 0.33


def rescale_rms(mag):
    # get values from simulated TESS LC (rms is scatter not mean ef)
    rms0, mag0 = 0.0017390130051, 8.02999973
    flux_ratio = 10**(-.4*(mag0-mag))
    return rms0 * np.sqrt(flux_ratio)


def get_timeseries(mag, Teff, Rs, Ms, Ps, rpRss, add_systematic=True,
		   Ndays_field=27):
    # get WF
    fname = 'tess2019128220341-0000000005712108-0016-s_lc.fits'
    hdu = download_one_fits(fname)
    hdr, bjdorig,_,_ = get_lc(hdu)
    # extend to longer baselines (i.e. more than 27 days)
    dt = 2.  # min
    bjd = bjdorig + 0
    assert 27 <= Ndays_field <= 351
    N_to_extend_baseline = int(Ndays_field/27.)
    for i in range(N_to_extend_baseline):
    	bjd = np.append(bjd, bjdorig+(bjd.max()-bjdorig.min())+dt)    
    g = (bjd-bjd.min() <= Ndays_field)
    bjd = bjd[g]
    
    # get uncertainties
    rms = rescale_rms(mag)
    ef = np.repeat(rms, bjd.size)

    # add planet models
    params, fmodel = get_planet_model(bjd, Ps, rpRss, Rs)
    fcorr = fmodel + np.random.randn(bjd.size) * rms

    # add eclipsing binaries (sometimes)
    ebmodel, EBparams = get_EB_model(bjd, Rs, Mdwarf_binary_frac=0)
    fcorr += ebmodel

    # add trends if desired
    if add_systematic:
	a = 10**(np.random.uniform(-3,-1))
	P = 10**(np.random.uniform(0,2))
	l = P * np.random.uniform(3,10)
	G = 10**(np.random.uniform(-1,1))
	s = 0
	k1 = george.kernels.ExpSquaredKernel(l)
	k2 = george.kernels.ExpSine2Kernel(G,P)
	gp = george.GP(k1+k2)
	# get time binning and bin the light curve
    	Npnts_per_timescale, Npntsmin, Npntsmax = 8., 5e2, 1e3
    	timescale_to_resolve = P / Npnts_per_timescale
    	Ttot = bjd.max() - bjd.min()
    	if Ttot/timescale_to_resolve < Npntsmin:
            dt = Ttot / Npntsmin
    	elif Ttot/timescale_to_resolve > Npntsmax:
            dt = Ttot / Npntsmax
    	else:
            dt = timescale_to_resolve
	tbin, fbin, efbin = boxcar(bjd, fcorr, ef, dt=dt, include_edges=True, tfull=bjd)
	fint = interp1d(tbin, gp.sample(tbin))
	fact = fint(bjd)
	fact -= np.median(fact)
	fact *= a/fact.max()
	f = fcorr + fact
    else:
	fact = np.zeros(f.size)
	f = fcorr
    
    # add stochastic jumps in flux (CRs?)
    Njumps = int(np.floor(np.random.exponential(2.5)))
    inds = np.arange(bjd.size, dtype=int)
    np.random.shuffle(inds)
    fjumps = np.zeros(bjd.size)
    fjumps[inds[:Njumps]] = np.random.uniform(3,8,Njumps)*rms
    f += fjumps

    # update stellar info
    hdr['TESSMAG'], hdr['TEFF'], hdr['RADIUS'] = mag, Teff, Rs
    logg = np.log10(6.67e-11*rvs.Msun2kg(Ms)*1e2 / rvs.Rsun2m(Rs)**2)
    hdr['LOGG'] = logg

    return hdr, bjd, f, ef, fact, fjumps, fcorr, params, EBparams


def get_planet_model(bjd, Ps, rpRss, Rs, N=1e2):
    # setup LC parameters
    Ps, rpRss = np.ascontiguousarray(Ps), np.ascontiguousarray(rpRss)
    Nplanets = len(Ps)
    T0s = np.random.uniform(bjd.min(), bjd.max(), Nplanets)
    aRs, inc_degs = np.zeros(Nplanets), np.zeros(Nplanets)

    # get planet models
    fmodel = np.zeros(bjd.size)
    for i in range(Nplanets):
	aRs[i] = rvs.AU2m(rvs.semimajoraxis(Ps[i],Rs,0)) / rvs.Rsun2m(Rs)
	b = np.random.uniform(-1,1)
	inc_degs[i] = rvs.inclination(Ps[i],Rs,Rs,b)
	theta = Ps[i], T0s[i], rpRss[i], aRs[i], inc_degs[i], 0., 90.
	fmodel += batman_transit_model(theta, bjd) - 1

    params = np.array([Ps, T0s, rpRss, aRs, inc_degs]).T
    return params, fmodel + 1


def batman_transit_model(theta, bjd):
    P, T0, rpRs, aRs, inc, ecc, omega = theta
    params = batman.TransitParams()     
    params.per = 1.
    params.t0  = 0.
    params.rp  = rpRs
    params.a   = aRs
    params.inc = inc
    params.ecc = ecc
    params.w   = omega
    params.limb_dark = 'quadratic'
    params.u         = [-0.0023, 0.1513]
    phase = (bjd-T0) / P
    m = batman.TransitModel(params, phase)
    return m.light_curve(params)


# binary frac from 79 M2-M4.5 dwarfs from .1-1e4 AU (http://adsabs.harvard.edu/abs/1997AJ....113.2246R)
def get_EB_model(bjd, Rs, Mdwarf_binary_frac=.3):
    inc = np.rad2deg(np.arccos(np.random.uniform(-1,1)))
    # sample Kepler EB periods
    Ps, Teffs = np.loadtxt('EBs/Kepler_EBC_v3.dat', delimiter=',', usecols=(1,9)).T
    g, P = (Teffs<=4e3) & (Teffs>0), 0
    while P <= 0:
    	P = np.random.choice(Ps[g]) * np.random.normal(1,.1)
    R2 = np.random.uniform(.08, Rs)  # radius of the stellar companion
    b = rvs.impactparam_inc(P, Rs, Rs, inc, mp_Mearth=rvs.kg2Mearth(rvs.Msun2kg(R2)))
    
    if (abs(b) <= 1) and (np.random.rand() <= Mdwarf_binary_frac):   # eclipse -> compute the LC model
    	T0 = np.random.uniform(bjd.min(), bjd.max())
    	a = rvs.AU2m(rvs.semimajoraxis(P, Rs, rvs.kg2Mearth(rvs.Msun2kg(R2))))
	r1, r2 = rvs.Rsun2m(Rs)/a, rvs.Rsun2m(R2)/a
	assert r1 < 1
	assert r2 < 1
	sbratio = (R2/Rs)**3
    	EBlc = ellc.lc(bjd, r1, r2, sbratio, inc, t_zero=T0, period=P, q=R2/Rs, 
		       ld_1='quad', ldc_1=[-0.0023, 0.1513], 
		       ld_2='quad', ldc_2=[-0.0023, 0.1513]) - 1
	assert np.all(np.isfinite(EBlc))
        EBparams = r1, r2, sbratio, inc, T0, P, R2/Rs
	return EBlc, EBparams

    else:   # no eclipsing binary
	return np.zeros(bjd.size), tuple(np.zeros(7))


def create_summary_image(fname_short, planetcols=['b','g','r','k','c'],
			 SNRthresh=5, pltt=True, label=False):
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(4,7)
    # plot raw LC with transits highlighted
    sens = loadpickle('Results/%s/Sensitivity_class'%fname_short)
    bjd, f, ef, mu, sig, fcorr = sens.bjd, sens.f, sens.ef, sens.mu, sens.sig, sens.fcorr
    ax1 = plt.subplot(gs[1,:])
    ax1.errorbar(bjd, f, ef, fmt='k.', capsize=0, alpha=.3)
    ax1.plot(bjd, mu, 'r-', lw=1.5)
    params_true = sens.params_true
    P_true, T0_true, rpRs_true, aRs_true, inc_true = params_true[0]
    rp_true = rvs.m2Rearth(rvs.Rsun2m(rpRs_true * sens.Rs))
    Nplanets = params_true.shape[0]
    for i in range(Nplanets):
	P, T0 = params_true[i,:2]
	Ntransits = int(np.round((bjd.max()-bjd.min()) / P))
	transit_times = T0 + np.arange(-Ntransits*2, Ntransits*2)*P
	for j in range(transit_times.size):
    	    ax1.axvline(transit_times[j], ymin=0, ymax=.05, color='b', lw=2)
    ax1.set_xlim((bjd.min(), bjd.max()))
    ax1.set_ylabel('Normalized Flux')
    ax1.set_title('Raw LC with systematic model shown and real transits marked')
    for i in range(2):
    	ax1.axvline(T0_true+i*P_true, ls='--', lw=.9, color='b')

    # print stellar parameters
    try:
	ax1.text(1.01, 2.15, 'Star:', transform=ax1.transAxes)
  	ax1.text(1.01, 2.0, 'Tmag = %.2f'%sens.Tmag, transform=ax1.transAxes)
	ax1.text(1.01, 1.85, 'R$_s$ = %.2f R$_{\odot}$'%sens.Rs, transform=ax1.transAxes)
	ax1.text(1.01, 1.7, 'M$_s$ = %.2f M$_{\odot}$'%sens.Ms, transform=ax1.transAxes)
	ax1.text(1.01, 1.55, 'T$_{eff}$ = %i K'%sens.Teff, transform=ax1.transAxes)
    except AttributeError:
	pass

    # plot lnLs of the transit times
    ax0 = plt.subplot(gs[0,:])
    transit_times = sens.transit_times
    durations = sens.durations
    lnLs = sens.lnLs_linearsearch
    assert lnLs.shape == (transit_times.size, durations.size)
    SNRs = (lnLs-np.median(lnLs, axis=0)) / np.std(lnLs, axis=0)
    SNRmeds = (lnLs-np.median(lnLs, axis=0)) / MAD(lnLs)
    i = 1
    ax0.plot(transit_times, SNRs[:,i], 'k-', lw=1.4)
    ax0.plot(transit_times, SNRmeds[:,i], 'g-', lw=.9)
    ax0.axhline(SNRthresh, ls='--', color='k')
    ax0.set_xlim(ax1.get_xlim()), ax0.set_ylabel('SNR')
    ax0.set_title('lnL of transits models at transit times (duration = %.2f hrs)'%(durations[i]*24))
    for i in range(2):
        ax0.axvline(T0_true+i*P_true, ls='--', lw=.9, color='b')

    # plot the corrected light curve and found transits
    ax2 = plt.subplot(gs[2,:])
    ax2.errorbar(bjd, fcorr, ef, fmt='k.', capsize=0, alpha=.3)
    params = sens.params_guess
    Nfound = params.shape[0]
    for i in range(Nfound):
	P, T0 = params[0,:2]
	Ntransits = int(np.round((bjd.max()-bjd.min()) / P))
        transit_times = T0 + np.arange(-Ntransits*2, Ntransits*2)*P
	for j in range(transit_times.size):
	    ax2.axvline(transit_times[j], ymin=0, ymax=.05, color='g', lw=2)
    ax2.set_xlim((bjd.min(), bjd.max()))
    ax2.set_ylabel('Normalized Flux'), ax2.set_xlabel('Time [BJD-2,457,000]')
    ax2.set_title('Corrected LC with %i detected transits marked'%Nfound)

    # plot phase-folded true transit model
    ax3 = plt.subplot(gs[3:,:3])
    phase = foldAt(bjd, P_true, T0_true)
    phase[phase > .5] -= 1
    s = np.argsort(phase)
    ax3.errorbar(phase, fcorr, ef, fmt='k.', capsize=0, alpha=.3)
    theta = P_true, T0_true, rpRs_true, aRs_true, inc_true, 0., 90.
    fmodel = batman_transit_model(theta, bjd)
    ax3.plot(phase[s], fmodel[s], 'b-', lw=3)
    rp = rvs.m2Rearth(rvs.Rsun2m(rpRs_true*sens.Rs))
    xwidth = rvs.transit_width(P_true, sens.Ms, sens.Rs, rp, 0) / P_true
    ax3.set_xlim((-1.5*xwidth, 1.5*xwidth))
    ax3.set_ylabel('Normalized Flux'), ax3.set_xlabel('Phase (P=%.4f days)'%P_true)

    # plot phase_folded planet detection
    ax4 = plt.subplot(gs[3:,3:6])
    if Nfound > 0:
        P, T0, depth,_ = params[0]
	T0new = T0 + np.arange(-100,100)*P
	T0 = float(T0new[abs(T0new-T0_true) == np.min(abs(T0new-T0_true))]) # rescale T0 to the true T0
        rpRs = np.sqrt(depth)
    	phase = foldAt(bjd, P, T0)
    	phase[phase > .5] -= 1
    	s = np.argsort(phase)
    	ax4.errorbar(phase, fcorr, ef, fmt='k.', capsize=0, alpha=.3)
    	theta = P, T0, rpRs, aRs_true, inc_true, 0., 90.
    	fmodel = batman_transit_model(theta, bjd)
    	ax4.plot(phase[s], fmodel[s], 'g-', lw=3)
    	rp = rvs.m2Rearth(rvs.Rsun2m(rpRs*sens.Rs))
    	xwidth = rvs.transit_width(P, sens.Ms, sens.Rs, rp, 0) / P
    	ax4.set_xlabel('Phase (P=%.4f days)'%P)
    ax4.set_xlim((-1.5*xwidth, 1.5*xwidth))
    ax4.set_ylim(ax3.get_ylim())
    ax4.set_yticklabels('')

    # report parameters
    b_true = float(rvs.impactparam_inc(P_true, sens.Ms, sens.Rs, inc_true))
    true_label = 'True parameters:\nP = %.4f days\nT0 = %.3f BJD\nrp/Rs = %.3f\nrp = %.2f R$_{\oplus}$\na/Rs = %.2f\nb = %.2f'%(P_true, T0_true, rpRs_true, rp_true, aRs_true, b_true)
    ax4.text(1.1, .34, true_label, transform=ax4.transAxes)
    found_label = 'Detected parameters:\n'
    if Nfound > 0:
       	found_label += 'P = %.4f days\nT0 = %.3f BJD\nrp/Rs = %.3f'%(P, T0, rpRs)
    else:
	found_label += 'no planets found'
    ax4.text(1.1, .21, found_label, verticalalignment='top', transform=ax4.transAxes)

    fig.subplots_adjust(left=.08, bottom=.07, right=.87, top=.96, hspace=.25)
    if label:
	plt.savefig('plots/summary_%s.png'%fname_short)
    if pltt:
        plt.show()
    plt.close('all')


def MAD(lnLs):
    mads = np.zeros(lnLs.shape[1])
    for i in range(mads.size):
	l = lnLs[:,i]
	mads[i] = np.median(abs(l-np.median(l)))
    return mads


def compute_sensitivity(fname, Ps, rpRss, Tmag, Rs, Ms, Teff,
			add_systematic=True, Ndays_field=27):
    # create timeseries and save
    assert len(Ps) == len(rpRss)
    hdr, bjd, f, ef, fact, fjumps, fcorr, params, EBparams = \
                                get_timeseries(Tmag, Teff, Rs, Ms, Ps, rpRss,
                                               add_systematic, Ndays_field)
    fname_short = fname.replace('.fits','')
    sens = Sensitivity(fname_short)
    sens.bjd, sens.f, sens.ef, sens.fact, sens.fjumps = bjd, f, ef, fact, fjumps
    sens.DONE = False
    sens.pickleobject()

    # save true parameters for cross-checking
    sens.Tmag, sens.Rs, sens.Ms, sens.Teff = Tmag, Rs, Ms, Teff
    sens.params_true = params
    sens.params_true_labels = np.array(['Ps', 'T0s', 'rp_Rss', 'a_Rss', 'inc_degs'])
    sens.rps = rvs.m2Rearth(rvs.Rsun2m(sens.params_true[:,2]*sens.Rs))
    sens.EBr1, sens.EBr2, sens.EBsbratio, sens.EBinc, sens.EBT0, sens.EBP, sens.EBq = EBparams
    sens.pickleobject()

    # fit systematics with a GP
    print 'Fitting LC with GP alone...\n'
    #samplerGP, samplesGP, resultsGP = do_mcmc_0(sens, bjd, f, ef, fname_short)
    #sens.samplesGP = samplesGP
    #sens.resultsGP = resultsGP
    sens.thetaGPall, sens.resultsGPall, thetaGPin, thetaGPout = do_optimize_0(bjd, f, ef)
    sens.thetaGP, sens.resultsGP = thetaGPin, thetaGPout
    sens.pickleobject()

    # search for transits in the corrected LC and get the transit parameters guesses
    print 'Searching for transit-like events...\n'
    params, EBparams, maybeEBparams = find_transits(sens, bjd, f, ef,
                                                    thetaGPout, hdr,
                                                    fname_short)
    sens.params_guess = params
    sens.params_guess_labels = np.array(['Ps', 'T0s', 'depths [Z]', \
                                         'durations [D]'])
    sens.EBparams_guess, sens.maybeEBparams_guess = EBparams, maybeEBparams
    sens.pickleobject()

    # are planets detected?
    sens.is_detected = np.array([int(np.any(np.isclose(params[:,0], Ps[i], rtol=.02)))
                                 for i in range(Ps.size)])
    assert sens.is_detected.size == sens.params_true.shape[0]
    sens.is_FP = np.array([int(np.invert(np.any(np.isclose(sens.params_guess[i,0],Ps,rtol=.02))))
                           for i in range(sens.params_guess.shape[0])])
    assert sens.is_FP.size == sens.params_guess.shape[0]
    sens.paramsFP_guess = sens.params_guess[sens.is_FP.astype(bool)]
    sens.pickleobject()

    # do joint GP+transit model
    if np.any(sens.is_detected.astype(bool)):
    	params, transit_params, resultsGPfin, mufin, sigfin, LDcoeffs = joint_LC_fit(sens)
    	sens.params_guessfin = params
	sens.transit_params = transit_params
	sens.u1, sens.u2 = LDcoeffs
    	sens.resultsGPfin = resultsGPfin
    	sens.mufin, sens.sigfin = mufin, sigfin    	

    sens.DONE = True
    sens.pickleobject()


def get_completeness_grid(prefix='TOIsensitivity351', pltt=True):
    '''Get the completeness as a function of P & rp/Rs for a set of sensitivity
    calculations'''
    # get folders
    folders = np.array(glob.glob('Results/%s*/Sensitivity_class'%prefix))
    Nfolders = folders.size
    assert Nfolders > 0

    # get planet/stellar variables and detections
    Ps, rps, detected = np.zeros(0), np.zeros(0), np.zeros(0)
    Rss, Mss, Teffs, Tmags = np.zeros(0), np.zeros(0), \
                             np.zeros(0), np.zeros(0)
    foldersdet = np.zeros(0)
    PsFP, rpsFP = np.zeros(0), np.zeros(0)
    RsFP, MsFP, TeffFP, TmagFP = np.zeros(0), np.zeros(0), \
                                 np.zeros(0), np.zeros(0)
    foldersFP = np.zeros(0)
    for i in range(Nfolders):
 	print float(i)/Nfolders
        assert Ps.size == rps.size
        assert Ps.size == detected.size
	sens = loadpickle(folders[i])
	if sens.DONE:
            # get injected planet parameters
            foldersdet = np.append(foldersdet, folders[i])
	    Ps = np.append(Ps, sens.params_true[:,0])
            rps = np.append(rps, sens.rps)
            detected = np.append(detected, sens.is_detected)

            # get stellar parameters
            Rss = np.append(Rss, sens.Rs)
            Mss = np.append(Mss, sens.Ms)
            Teffs = np.append(Teffs, sens.Teff)
            Tmags = np.append(Tmags, sens.Tmag)

            # get false positive parameters
            NFPs = sens.paramsFP_guess.size
            for i in range(NFPs):
                foldersFP = np.append(foldersFP, folders[i])
                PsFP = np.append(PsFP, sens.paramsFP_guess[i,0])
                rp = rvs.m2Rearth(rvs.Rsun2m(np.sqrt(sens.paramsFP_guess[i,2])*sens.Rs))
                rpsFP = np.append(rpsFP, rp)
                RsFP = np.append(RsFP, sens.Rs)
                MsFP = np.append(MsFP, sens.Ms)
                TeffFP = np.append(TeffFP, sens.Teff)
                TmagFP = np.append(TmagFP, sens.Tmag)

	else:
	    pass

    assert Ps.size == rps.size
    assert Ps.size == detected.size
    sensgrid = Sensitivity_grid(prefix)
    sensgrid.foldersdet, sensgrid.Ps, sensgrid.rps, sensgrid.detected = foldersdet, Ps, \
                                                                        rps, detected
    sensgrid.Rss, sensgrid.Mss, sensgrid.Teffs, sensgrid.Tmags = Rss, Mss, Teffs, Tmags 
    sensgrid.foldersFP, sensgrid.PsFP, sensgrid.rpsFP = foldersFP, PsFP, rpsFP
    sensgrid.RsFP, sensgrid.MsFP, sensgrid.TeffFP, sensgrid.TmagFP = RsFP, MsFP, \
                                                                     TeffFP, TmagFP

    # get completeness
    Pgrid = np.logspace(np.log10(.5),np.log10(27.4),12)
    rpgrid = np.logspace(np.log10(.5),np.log10(15),9)
    Ndet, Ntrue = np.zeros((Pgrid.size-1, rpgrid.size-1)), \
		  np.zeros((Pgrid.size-1, rpgrid.size-1))
    for i in range(Pgrid.size-1):
	for j in range(rpgrid.size-1):
	    g = (Ps >= Pgrid[i]) & (Ps <= Pgrid[i+1]) & \
		(rps >= rpgrid[j]) & (rps <= rpgrid[j+1])
	    Ndet[i,j], Ntrue[i,j] = float(detected[g].sum()), detected[g].size

    # save grids
    sensgrid.Pgrid, sensgrid.Pgrid = Pgrid, rpgrid
    sensgrid.Ndet, sensgrid.Ntrue, sensgrid.sensitivity = Ndet, Ntrue, Ndet/Ntrue
    sensgrid.pickleobject()
    
    # plot grid
    if pltt:
	plt.figure()
	plt.pcolormesh(Pgrid, rpgrid, (Ndet/Ntrue).T, 
		       cmap=truncate_colormap(plt.get_cmap('bone_r'),.1,1),
		        vmin=0, vmax=1)
	plt.colorbar()
	plt.xlim((Pgrid.min(),Pgrid.max())), plt.ylim((rpgrid.min(),rpgrid.max()))
	plt.xscale('log'), plt.yscale('log')
	plt.xlabel('Period [days]'), plt.ylabel('r$_p$ [R$_{\oplus}$')
	plt.show()


def add_planets(Nplanets, Ps, rps, Ms, Ndays_field):
    '''add additional planet to be transiting and ensure dynamical stability.'''
    assert Ps.size == 1
    assert rps.size == 1

    Nplanets = int(Nplanets)
    if Nplanets > 1:
        # ensure stability is multiple planets
    	unstable, ind = True, 0
    	while unstable and ind < 100:
	    # sample additional planets
	    Ps = np.append(Ps[0], 10**np.random.uniform(np.log10(.5), np.log10(Ndays_field), Nplanets-1))
	    rps = np.append(rps[0], 10**np.random.uniform(np.log10(.5), np.log10(15), Nplanets-1))
	    sort = np.argsort(Ps)
	    Ps, rps = Ps[sort], rps[sort]

	    # check stability
      	    Ls = .23 * 3.828e26 * Ms**2.3
	    smas = rvs.AU2m(rvs.semimajoraxis(Ps, Ms, 0))
	    Fs = Ls / (4*np.pi*smas**2) / 1367.
    	    mps = radF2mass(rps, Fs)
	    assert Ps.size == Nplanets
            assert mps.size == Nplanets
	    unstable = False if np.all(rvs.is_Lagrangestable(Ps, Ms, mps, np.zeros(Nplanets)).astype(bool)) else True
	    ind += 1
	return Ps, rps

    elif Nplanets == 1:
	return Ps, rps

    else:
	raise ValueError('Need at least 1 planet. (Nplanets = %i)'%Nplanets)



def do_i_run(fname):
    if os.path.exists('%s/Sensitivity_class'%fname):
	self = loadpickle('%s/Sensitivity_class'%fname)
	return not self.DONE
    else:
	return True


if __name__ == '__main__':
    # create light curve for this planet
    fnum = int(sys.argv[1])
    Pin = float(sys.argv[2])
    Pout = float(sys.argv[3])
    assert Pin < Pout
    rpin = float(sys.argv[4])
    rpout = float(sys.argv[5])
    assert rpin < rpout
    Ndays_field = float(sys.argv[6])

    # cut in Tmag such that the sample contains ~1e4 stars from which we expect to detect 1e4*.01~1e2 planets
    #Tmag_fullsky = np.repeat(Tmag,12)  # Tmag only covers 1/12 of the sky
    #plt.hist(Tmag_fullsky, bins=1000, log=True, histtype='step', cumulative=True), plt.axhline(1e4), plt.show()
    Tmag, Teff, Rs, Ms = np.loadtxt('CTL/CTL_RC.dat', delimiter=',').T
    g = (Tmag < 11.1) & (Teff<4e3)
    Tmag, Teff, Rs, Ms = Tmag[g], Teff[g], Rs[g], Ms[g]

    # for this planet sample
    Nplanets, Nstars_per_star = 2, 20
    for i in range(Nstars_per_star):
        fname = 'TOIsensitivity27vett_mult%i_planet%.5d_iteration%.3d'%(Nplanets,fnum,i)
	if do_i_run('Results/%s'%fname):
            g = np.random.randint(0,Tmag.size)  # get random star from the candidate target list
    	    # get planet params
    	    Ps, rps = np.random.uniform(Pin,Pout,1), np.random.uniform(rpin,rpout,1)
    	    Ps, rps = add_planets(Nplanets, Ps, rps, Ms[g], Ndays_field)
            rpRs = rvs.Rearth2m(rps) / rvs.Rsun2m(Rs[g])
            compute_sensitivity(fname, Ps, rpRs, Tmag[g], Rs[g], Ms[g], Teff[g],
                            	Ndays_field=Ndays_field)
