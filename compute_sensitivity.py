from TESS_search import *
import matplotlib.gridspec as gridspec
import linear_lnlike as llnl
import batman, ellc
from truncate_cmap import *
from sensitivity_class import *


# M dwarfs from the TESS CTL (candidate target list) only 
global medTmag, medTeff, medRs
medTmag, medTeff, medRs = 10.5, 3400, 0.33


def scale_rms(mag):
    # get values from simulated TESS LC (rms is scatter not mean ef)
    rms0, mag0 = 0.0017390130051, 8.02999973
    flux_ratio = 10**(-.4*(mag0-mag))
    return rms0 * np.sqrt(flux_ratio)


def get_timeseries(mag, Teff, Rs, Ps, rpRss, add_systematic=True,
		   N_to_extend_baseline=0):
    # get WF
    fname = 'tess2019128220341-0000000005712108-0016-s_lc.fits'
    hdu = download_one_fits(fname)
    hdr, bjdorig,_,_ = get_lc(hdu)
    # extend to longer baselines (i.e. more than 27 days)
    dt = np.median(np.diff(bjdorig))
    bjd = bjdorig + 0
    for i in range(int(N_to_extend_baseline)):
    	bjd = np.append(bjd, bjdorig+(bjd.max()-bjdorig.min())+dt)    

    # get uncertainties
    rms = scale_rms(mag)
    ef = np.repeat(rms, bjd.size)

    # add planet models
    params, fmodel = get_planet_model(bjd, Ps, rpRss, Rs)
    fcorr = fmodel + np.random.randn(bjd.size) * rms

    # add eclipsing binaries (sometimes)
    ebmodel, EBparams = get_EB_model(bjd, Rs)
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
	tbin, fbin, efbin = boxcar(bjd, fcorr, ef, include_edges=True)
	fint = interp1d(tbin, gp.sample(tbin))
	fact = fint(bjd)
	fact -= np.median(fact)
	f = fcorr + fact * (a/fact.max())
    else:
	f = fcorr
    
    # update stellar info
    hdr['TESSMAG'], hdr['TEFF'], hdr['RADIUS'] = mag, Teff, Rs
    
    return hdr, bjd, f, ef, fcorr, params, EBparams


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
			add_systematic=True, N_to_extend_baseline=0):
    # create timeseries and save
    hdr, bjd, f, ef, fcorr, params, EBparams = get_timeseries(Tmag, Teff, Rs, Ps, rpRss,
					            	      add_systematic, 
						      	      N_to_extend_baseline)
    fname_short = fname.replace('.fits','')
    sens = Sensitivity(fname_short)
    sens.add_raw_timeseries(bjd, f, ef)

    # save true parameters for cross-checking
    sens.add_star((Tmag, Rs, Ms, Teff))
    sens.add_params_true(params)
    sens.add_EBparams_true(EBparams)
    ##save_fits(params, 'Results/%s/params_true'%fname_short)

    # fit systematics with a GP
    print 'Fitting LC with GP alone...\n'
    samplerGP, samplesGP, resultsGP = do_mcmc_0(sens, bjd, f, ef, fname_short)
    sens.add_samplesGP(samplesGP)
    sens.add_resultsGP(resultsGP)
    ##save_fits(samplesGP, 'Results/%s/GP_samples'%fname_short)
    ##save_fits(resultsGP, 'Results/%s/GP_results'%fname_short)

    # search for transits in the corrected LC and get the transit parameters guesses
    print 'Searching for transit-like events...\n'
    params = find_transits(sens, bjd, f, ef, resultsGP[0], hdr, fname_short)
    sens.add_params_guess(params)
    ##save_fits(params, 'Results/%s/params_guess'%fname_short)

    # is the planet detected?
    detected = True if np.any(np.isclose(params[:,0], Ps[0], atol=Ps[0]*.02)) else False
    sens.add_detected(np.array([detected]).astype(int))
    ##save_fits(np.array([detected]).astype(int), 'Results/%s/is_detected'%fname_short)



def get_completeness_grid(prefix='TOIsensitivity351', pltt=True):
    '''Get the completeness as a function of P & rp/Rs for a set of sensitivity
    calculations'''
    # get folders
    folders = np.array(glob.glob('Results/%s_*/Sensitivity_class'%prefix))
    Nfolders = folders.size
    assert Nfolders > 0

    # get variables and detections
    Ps, rps, detected = np.zeros(0), np.zeros(0), np.zeros(0)
    for i in range(Nfolders):
 	print float(i)/Nfolders
	sens = loadpickle(folders[i])
	try:
	    _ = getattr(sens, 'is_detected')
	    Ps = np.append(Ps, sens.params_true[0,0])
	    rp = rvs.m2Rearth(rvs.Rsun2m(sens.params_true[0,2]*sens.Rs))
            rps = np.append(rps, rp)
            detected = np.append(detected, sens.is_detected)
	except AttributeError:
	    pass

    # get completeness
    Pgrid = np.logspace(np.log10(.5),np.log10(9.5),11)
    rpgrid = np.logspace(np.log10(.5),np.log10(15),7)
    #Pgrid, rpgrid = np.unique(Ps), np.unique(rps)
    Ndet, Ntrue = np.zeros((Pgrid.size-1, rpgrid.size-1)), \
		  np.zeros((Pgrid.size-1, rpgrid.size-1))
    for i in range(Pgrid.size-1):
	for j in range(rpgrid.size-1):
	    g = (Ps >= Pgrid[i]) & (Ps <= Pgrid[i+1]) & \
		(rps >= rpgrid[j]) & (rps <= rpgrid[j+1])
	    Ndet[i,j], Ntrue[i,j] = float(detected[g].sum()), detected[g].size
    completeness = Ndet / Ntrue

    # save fits
    save_fits(Pgrid, 'Results/%s_Pgrid'%prefix)
    save_fits(rpgrid, 'Results/%s_rpgrid'%prefix)
    save_fits(Ndet, 'Results/%s_Ndet'%prefix)
    save_fits(Ntrue, 'Results/%s_Ntrue'%prefix)

    # plot grid
    if pltt:
	plt.figure()
	plt.pcolormesh(Pgrid, rpgrid, completeness.T, 
		       cmap=truncate_colormap(plt.get_cmap('bone_r'),.1,1),
		        vmin=0, vmax=1)
	plt.colorbar()
	plt.xlim((Pgrid.min(),Pgrid.max())), plt.ylim((rpgrid.min(),rpgrid.max()))
	plt.xscale('log'), plt.yscale('log')
	plt.xlabel('Period [days]'), plt.ylabel('rp')
	plt.show()



if __name__ == '__main__':
    fnum = int(sys.argv[1])
    P = float(sys.argv[2])
    rp = float(sys.argv[3])
    Tmag, Teff, Rs, Ms = np.loadtxt('CTL/CTL_RC.dat', delimiter=',').T
    # cut in Tmag such that the sample contains ~1e4 stars from which we expect to detect 1e4*.01~1e2 planets
    #Tmag_fullsky = np.repeat(Tmag,12)  # Tmag only covers 1/12 of the sky
    #plt.hist(Tmag_fullsky, bins=1000, log=True, histtype='step', cumulative=True), plt.axhline(1e4), plt.show()
    g = (Tmag < 11.1) & (Teff<4e3)
    Tmag, Teff, Rs, Ms = Tmag[g], Teff[g], Rs[g], Ms[g]
    Nplanets_per_star = 20
    for i in range(Nplanets_per_star):
	fname = 'TOIsensitivity351_%.5d_%.3d'%(fnum,i)
	g = np.random.randint(0,Tmag.size)
	rpRs = rvs.Rearth2m(rp) / rvs.Rsun2m(Rs[g])
	compute_sensitivity(fname, [P], [rpRs], Tmag[g], Rs[g], Ms[g], Teff[g],
			    N_to_extend_baseline=np.ceil(P/27.)-1)
