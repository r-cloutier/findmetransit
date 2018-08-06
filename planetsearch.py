from compute_sensitivity import *
from TESSLC_class import *


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
    tesslc.pickleobject()

    # 

    

if __name__ == '__main__':
    fits_name = sys.argv[1]
    planet_search(fits_name)
