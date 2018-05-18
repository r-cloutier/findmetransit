from imports import *

def compute_LSperiodogram(t, y, ey=np.zeros(0), plims=(.5,2e2)):
    '''Compute the LS periodogram of an input timeseries.'''
    if ey.size == 0:
        ey = np.ones(t.size)
    pmin, pmax = plims
    freqs = 1./np.logspace(np.log10(pmin), np.log10(pmax), 1e3)
    power = LombScargle(t, y, ey).power(freqs)
    powernorm = power / power.std()
    return 1./freqs, power, powernorm
