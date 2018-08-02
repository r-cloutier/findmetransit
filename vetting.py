from imports import *
import rvs, batman
from scipy.interpolate import LinearNDInterpolator as lint


def identify_EBs(params, bjd, fcorr, ef, Ms, Rs, Teff,
                 SNRthresh=3., rpmax=30, detthresh=5):
    '''For each proposed planet in params, run through a variety of checks to 
    vetting the planetary candidates and identify EB false positives.'''
    Nplanets = params.shape[0]
    paramsout, isEB, maybeEB = np.zeros_like(params), np.zeros(Nplanets), \
                               np.zeros(Nplanets)
    EBconditions = np.zeros((Nplanets, 5)).astype(bool)
    EBcondition_labels = np.array(['rpRs > 0.5', 'rp>%.1f'%rpmax,
                                   'secondary eclipse detected',
                                   'transit is V-shaped',
                                   'transit duration is too long for a planet'])
    for i in range(Nplanets):

        # get best fit parameters
        paramsout[i] = _fit_params(params[i], bjd, fcorr, ef, Ms, Rs, Teff)
        
        # ensure the planet is not too big
        rpRs = np.sqrt(params[i,2])
        isEB[i] = 1 if rpRs > .5 else isEB[i]
        EBconditions[i,0] = True if rpRs > .5 else False
        rp = rvs.m2Rearth(rvs.Rsun2m(rpRs*Rs))
        isEB[i] = 1 if rp > rpmax else isEB[i]
        EBconditions[i,1] = True if rp > rpmax else False

        # check for secondary eclipse
        eclipse = _is_eclipse(params[i], bjd, fcorr, ef, detthresh)
        isEB[i] = 1 if eclipse else isEB[i]
        EBconditions[i,2] = True if eclipse else False
        
        # check for ellipsoidal variations
        #ellipsoidal = _is_ellipsoidal()
        #isEB[i] = 1 if ellipsoidal else isEB[i]
        # how can I do this without knowing the parameters of the binary?

        # flag V-shaped transits (does not implies an EB)
        Vshaped, duration = _is_Vshaped(params[i], bjd, fcorr, ef, Ms, Rs, Teff)
        maybeEB[i] = 1 if Vshaped else maybeEB[i]
        EBconditions[i,3] = True if Vshaped else False

        # is duration reasonable? # Gunther+2016
        EBduration = _is_EB_duration(duration, paramsout[i,0], Ms, Rs)
        isEB[i] = 1 if EBduration else isEB[i]
        EBconditions[i,4] = True if EBduration else False

    # save planet and EB parameters
    maybeEB[isEB==1] = 1
    isEB, maybeEB = isEB.astype(bool), maybeEB.astype(bool)
    params, EBparams, maybeEBparams = params[np.invert(isEB)], params[isEB], \
                                      params[maybeEB]
    return params, EBparams, maybeEBparams, EBconditions, EBcondition_labels



def _fit_params(params, bjd, fcorr, ef, Ms, Rs, Teff):
    '''Get best-fit parameters.'''
    assert params.shape == (4,)
    P, T0, depth, duration = params
    u1, u2 = _get_LDcoeffs(Ms, Rs, Teff)
    aRs = rvs.AU2m(rvs.semimajoraxis(P,Ms,0)) / rvs.Rsun2m(Rs)
    rpRs = np.sqrt(depth)
    p0 = aRs, rpRs, 90.
    bnds = ((aRs*.9, 0, float(rvs.inclination(P,Ms,Rs,1.1))),
            (aRs*1.1, 1, float(rvs.inclination(P,Ms,Rs,-1.1))))
    try:
        popt,_ = curve_fit(transit_model_func_curve_fit(P,T0,u1,u2),
                           bjd, fcorr, p0=p0, sigma=ef,
                           absolute_sigma=True, bounds=bnds)
        aRs, rpRs, inc = popt
        depth = rpRs**2
        b = rvs.impactparam_inc(P, Ms, Rs, inc)
        duration = P/(np.pi*aRs) * np.sqrt((1+np.sqrt(depth))**2 - b*b)
        return P, T0, depth, duration
    except RuntimeError:
        return params


def _box_transit_model(theta, bjd):
    '''Return the box transmit model.'''
    fmodel = np.ones(bjd.size)
    assert theta.shape == (4,)
    P, T0, depth, duration = theta
    phase = foldAt(t, P, T0)
    phase[phase>.5] -= 1
    intransit = (phase*P <= .5*duration) & (phase*P >= -.5*duration)
    fmodel[intransit] = 1. - depth
    return fmodel

    
def lnlike(theta, bjd, fcorr, ef):
    #P, T0, depth, duration = theta
    fmodel = _box_transit_model(theta, bjd)
    inv_sigma2 = 1. / (ef**2 + fmodel**2)
    return -0.5*(np.sum((fcorr - fmodel)**2 * inv_sigma2 - np.log(inv_sigma2)))


def _is_eclipse(params, bjd, fcorr, ef, DT):
    '''check for secondary eclipses in the light curve and return True if the 
    signal is detmerined to be from an EB'''
    assert params.shape == (4,)
    P, T0, depth, duration = params

    # get times in transit and occultation
    phase = foldAt(bjd, P, T0)
    phase[phase > .5] -= 1
    intransit = (phase >= -.5*duration/P) & (phase <= .5*duration/P)
    inoccultation = (phase <= .5*(duration/P - 1)) | \
                    (phase >= -.5*(duration/P - 1))

    # define EB criteria
    depth_transit, var_transit = depth, np.std(fcorr[intransit])**2
    depth_occultation = 1-np.median(fcorr[inoccultation])
    var_occultation = np.std(fcorr[inoccultation])**2
    return (depth_occultation / np.sqrt(var_occultation) > DT) & \
        ((depth_transit - depth_occultation) / np.sqrt(var_transit + \
                                                       var_occultation) > DT)


def _is_ellipsoidal():
    '''check for ellipsoidal variations that are indicative of a close binary. 
    see sect 8.1 Sullivan+2015'''
    return None


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


def transit_model_func_curve_fit(P, T0, u1, u2):
    def transit_model_func_in(bjd, aRs, rpRs, inc):
        p = batman.TransitParams()
        p.t0, p.per, p.rp = 0., 1., rpRs
        p.a, p.inc, p.ecc = aRs, inc, 0.
        p.w, p.limb_dark, p.u = 90., 'quadratic', [u1,u2]
        phase = foldAt(bjd, P, T0)
        m = batman.TransitModel(p, phase)
        f = m.light_curve(p)
        return f
    return transit_model_func_in
    

def _is_Vshaped(params, bjd, fcorr, ef, Ms, Rs, Teff):
    '''check if transit is V-shaped although this does not imply an EB 
    because inclined transiting planets can look like this as well.'''
    # fit transit model
    assert params.shape == (4,)
    P, T0, depth = params[:3]
    u1, u2 = _get_LDcoeffs(Ms, Rs, Teff)
    aRs = rvs.AU2m(rvs.semimajoraxis(P,Ms,0)) / rvs.Rsun2m(Rs)
    rpRs = np.sqrt(depth)
    p0 = aRs, rpRs, 90.
    bnds = ((aRs*.9, 0, float(rvs.inclination(P,Ms,Rs,1.1))),
            (aRs*1.1, 1, float(rvs.inclination(P,Ms,Rs,-1.1))))
    try:
        popt,_ = curve_fit(transit_model_func_curve_fit(P,T0,u1,u2),
                           bjd, fcorr, p0=p0, sigma=ef,
                           absolute_sigma=True, bounds=bnds)
        aRs, rpRs, inc = popt
    except RuntimeError:
        return False, 0.

    # get ingress, egress times and duration
    transit_func = transit_model_func_curve_fit(P,T0,0,0)
    fmodel = transit_func(bjd, *popt)
    phase = foldAt(bjd, P, T0)
    phase[phase>.5] -= 1
    depth = fmodel.min()
    T1, T4 = phase[(phase<=0) & (fmodel==1)].max()*P, \
             phase[(phase>=0) & (fmodel==1)].min()*P
    in_ingress = (phase<=0) & np.isclose(fmodel,depth,rtol=1e-4)
    T2 = phase[in_ingress].min()*P if in_ingress.sum() > 0 else (T4-T1)*.1
    in_egress = (phase>=0) & np.isclose(fmodel,depth,rtol=1e-4)
    T3 = phase[in_egress].max()*P if in_egress.sum() > 0 else (T4-T1)*.1
    Tingress, Tegress, duration = T2-T1, T4-T3, T4-T1
    Tedge = Tingress + Tegress
    
    # V-shaped if T2==T3 or if ingress+egress time is the same as the duration 
    if T2 == T3:
        return True, duration
    elif (Tedge >= duration*.99) & (Tedge <= duration*1.01):
        return True, duration
    else:
        return False, duration


def _is_EB_duration(duration, P, Ms, Rs, rpmax=30):
    '''If the duration is longer than the maximum allowed planet duration 
    for a given period and star, then flag as an EB.'''
    durationmax = rvs.transit_width(P, Ms, Rs, rpmax, 0)
    return duration > durationmax
