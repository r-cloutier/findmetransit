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


def transit_model_func(bjd, P, T0, aRs, rpRs, inc, u1, u2):
    p = batman.TransitParams()
    p.t0, p.per, p.rp = 0., 1., rpRs
    p.a, p.inc, p.ecc = aRs, inc, 0.
    p.w, p.limb_dark, p.u = 90., 'quadratic', [u1,u2]
    phase = foldAt(bjd, P, T0)
    m = batman.TransitModel(p, phase)
    f = m.light_curve(p)
    return f
    

def joint_LC_fit(sens):#params, thetaGP, bjd, f, fcorr, ef):
    '''Iteratively optimize the planetary and GP model parameters'''
    # optimize planets in detrended LC
    Nplanets = sens.params_guess.shape[0]
    transit_model = np.zeros(bjd.size)
    for i in range(Nplanets):
        
	curve_fit()

    return paramsout, resultsGP, mu, sig
