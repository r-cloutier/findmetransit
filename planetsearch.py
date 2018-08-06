from compute_sensitivity import *
from TESSLC_class import *
from TESS_search import *
import linear_lnlike as llnl
import batman, ellc
from massradius import radF2mass
from truncate_cmap import *
from joint_LCmodel import *


def read_data(fits_name):
    '''Get one light curve and stellar parameters.'''
    bjd, f, ef = None, None, None # TEMP
    Rs, Teff, Tmag = hdu.header['RADIUS'], hdu.headr['TEFF'], \
                     hdu.header['TESSMAG']
    Ms = None
    return bjd, f, ef, Ms, Rs, Teff, Tmag


def is_star_of_interest(Ms, Rs, Teff, Tmag):
    logg = np.log10(6.67e-10*rvs.Msun2kg(Ms)*1e2 / rvs.Rsun2m(Rs)**2)
    loggmin, Teffmax, Tmagmax = 3, 4e3, 12
    print 'stellar params: ', logg, Teff, Tmag
    if (logg >= loggmin) & (Teff <= Teffmax) & (Tmag <= Tmagmax):
        return True
    else:
        return False


def planet_search(fits_name):
    '''Run a planet search using the pipeline defined in compute_sensitivity 
    to search for planets.'''
    # get data and only run the star if it is of interest
    bjd, f, ef, Ms, Rs, Teff, Tmag = read_data(fits_name)
    assert is_star_of_interest(Ms, Rs, Teff, Tmag)
    tesslc = TESSLC(fits_name)
    tesslc.bjd, tesslc.f, tesslc.ef = bjd, f, ef
    tesslc.Ms, tesslc.Rs, tesslc.Teff, tesslc.Tmag = Ms, Rs, Teff, Tmag
    tesslc.DONE = False
    tesslc.pickleobject()

    # fit initial GP
    thetaGPall, resultsGPall, thetaGPin, thetaGPout = do_optimize_0(bjd, f, ef)
    tesslc.thetaGPall, tesslc.resultsGPall = thetaGPall, resultsGPall
    tesslc.thetaGPin, tesslc.thetaGPout = thetaGPin, thetaGPout
    tesslc.pickleobject()
 
    # search for transits in the corrected LC and get the transit parameters
    # guesses
    print 'Searching for transit-like events...\n'
    params, EBparams, maybeEBparams = find_transits(tesslc, bjd, f, ef,
                                                    thetaGPout, hdr,
                                                    fname_short)
    tesslc.params_guess = params
    tesslc.params_guess_labels = np.array(['Ps', 'T0s', 'depths [Z]', \
                                           'durations [D]'])
    tesslc.EBparams_guess, tesslc.maybeEBparams_guess = EBparams, maybeEBparams
    tesslc.pickleobject()

    # are planets detected?
    tesslc.is_detected = np.array([int(np.any(np.isclose(params[:,0], Ps[i],
                                                         rtol=.02)))
                                   for i in range(Ps.size)])
    assert tesslc.is_detected.size == tesslc.params_true.shape[0]
    tesslc.is_FP = np.array([int(np.invert(np.any(np.isclose(tesslc.params_guess[i,0],Ps,rtol=.02))))
                             for i in range(tesslc.params_guess.shape[0])])
    assert tesslc.is_FP.size == tesslc.params_guess.shape[0]
    tesslc.paramsFP_guess = tesslc.params_guess[tesslc.is_FP.astype(bool)]
    tesslc.pickleobject()

    # do joint GP+transit model
    if np.any(tesslc.is_detected.astype(bool)):
        params, transit_params, resultsGPfin, mufin, sigfin, LDcoeffs = \
                                                            joint_LC_fit(tesslc)
        tesslc.params_guessfin = params
        tesslc.transit_params = transit_params
        tesslc.u1, tesslc.u2 = LDcoeffs
        tesslc.resultsGPfin = resultsGPfin
        tesslc.mufin, tesslc.sigfin = mufin, sigfin

    tesslc.DONE = True
    tesslc.pickleobject()
    

if __name__ == '__main__':
    fits_name = sys.argv[1]
    planet_search(fits_name)
