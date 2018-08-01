from sensitivity_class import *
import rvs, batman
from scipy.interpolate import LinearNDInterpolator as lint


def get_LDcoeffs(Teff, Ms, Rs, Z=0):
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


def transit_model_func(bjd, P, T0, aRs, rpRs, inc, **kwargs):
    p = batman.TransitParams()
    p.t0, p.per, p.rp = 0., 1., rpRs
    p.a, p.inc, p.ecc = aRs, inc, 0.
    u1, u2 = kwargs['u1'], kwargs['u2']
    p.w, p.limb_dark, p.u = 90., 'quadratic', [u1,u2]
    phase = foldAt(bjd, P, T0)
    m = batman.TransitModel(p, phase)
    f = m.light_curve(p)
    return f


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



def joint_LC_fit(sens, Nsig=3, medkernel=49):
    '''Iteratively optimize the planetary and GP model parameters'''
    # get fixed LD coefficients for this star
    u1, u2 = get_LDcoeffs(sens.Teff, sens.Ms, sens.Rs)
    lds = {'u1':u1, 'u2':u2}    

    # optimize planets in detrended LC
    Nplanets = sens.params_guess.shape[0]
    transit_model = np.zeros(sens.tbin.size)
    for i in range(Nplanets):

        # make initial transit parameter guess
        P, T0, depth = sens.params_guess[i,:3]
        aRs = rvs.AU2m(rvs.semimajoraxis(P,sens.Ms,0)) / rvs.Rsun2m(sens.Rs)
        rpRs = np.sqrt(depth)
        p0 = P, T0, aRs, rpRs, 90.
        bnds = ((P*.98, T0-.5*P, aRs*.9, 0,
                 float(rvs.inclination(P,sens.Ms,sens.Rs,1))),
                (P*1.02, T0+.5*P, aRs*1.1, .5,
                 float(rvs.inclination(P,sens.Ms,sens.Rs,-1))))

        # optimize transit model parameters
        popt,_ = curve_fit(transit_model_func_curve_fit(u1, u2),
                           sens.bjd, sens.fcorr, p0=p0, sigma=sens.ef,
                           absolute_sigma=True, bounds=bnds)
        transit_model += transit_model_func(sens.tbin, *popt, **lds) - 1
    transit_model += 1
    
    # optimize GP model in the presence of the transit model
    a, l, G, Pgp = np.exp(sens.resultsGP)
    k1 = george.kernels.ExpSquaredKernel(l)
    k2 = george.kernels.ExpSine2Kernel(G,Pgp)
    gp = george.GP(a*(k1+k2))
    # trim outliers and median filter to avoid fitting deep transits
    resbin = sens.fbin - transit_model
    g = abs(resbin-np.median(resbin)) <= Nsig*resbin.std()
    tbin, resbin, efbin = sens.tbin[g], medfilt(resbin[g],medkernel), \
                          sens.efbin[g]
    try:
        results = gp.optimize(tbin, resbin, efbin)
        resultsGP = results[0]
        gp.compute(tbin, efbin)
        mubin, covbin = gp.predict(resbin, sens.tbin)
        sigbin = np.sqrt(np.diag(covbin))
        fintmu, fintsig = interp1d(sens.tbin, mubin), \
                          interp1d(sens.tbin, sigbin)
        mu, sig = fintmu(sens.bjd), fintsig(sens.bjd)
    except (ValueError, np.linalg.LinAlgError):
        resultsGP, mu, sig = np.zeros(len(thetaGP)), np.zeros(sens.bjd.size), \
                             np.zeros(sens.bjd.size)

    # optimize planets again on the newly detrended LC
    fcorr2 = sens.f - mu + 1 if mu.sum() > 0 else sens.f - mu
    transit_params = np.zeros((Nplanets, 5))
    for i in range(Nplanets):

        # make initial transit parameter guess
        P, T0, depth = sens.params_guess[i,:3]
        aRs = rvs.AU2m(rvs.semimajoraxis(P,sens.Ms,0)) / rvs.Rsun2m(sens.Rs)
        rpRs = np.sqrt(depth)
        p0 = P, T0, aRs, rpRs, 90.
        bnds = ((P*.98, T0-.5*P, aRs*.9, 0,
                 float(rvs.inclination(P,sens.Ms,sens.Rs,1))),
                (P*1.02, T0+.5*P, aRs*1.1, .5,
                 float(rvs.inclination(P,sens.Ms,sens.Rs,-1))))

        # optimize transit model parameters
        transit_params[i],_ = curve_fit(transit_model_func_curve_fit(u1, u2),
                                        sens.bjd, fcorr2, p0=p0,
                                        sigma=sens.ef, absolute_sigma=True,
                                        bounds=bnds)
    # get final parameters
    Ps, T0s = transit_params[:,:2].T
    depths = transit_params[:,3]**2
    rps = rvs.m2Rearth(rvs.Rsun2m(transit_params[:,3]*sens.Rs))
    bs = rvs.impactparam_inc(Ps, sens.Ms, sens.Rs, transit_params[:,4])
    durations = rvs.transit_width(Ps, sens.Ms, sens.Rs, rps, bs)
    paramsout = np.array([Ps, T0s, depths, durations]).T

    return paramsout, transit_params, resultsGP, mu, sig
