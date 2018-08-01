import numpy as np
import rvs

def rad2mass1(rad):
    '''Compute the mass of a planet from its radius in Earth radii assuming
    a step-wise density of Earth-like for rp <= 1.6 Rearth and Neptune-like
    otherwise.'''
    rhoN_rhoE = 1.638 / 5.514
    rhoJ_rhoE = 1.326 / 5.514
    if rad <= 1.6:
        return rad**3
    elif rad > 1.6 and rad < 234:
        return rhoN_rhoE * rad**3
    else:
        return rhoJ_rhoE * rad**3

def rad2mass2(rp, maxE=1.5, maxN=4, mean_only=True):
    '''Convert the input radii in Earth radii to masses using the
    analytical model from Weiss & Marcy 2014.'''
    rhoJ = 1.326  # cgs
    rmsE, rmsN = 2.7, 4.7      # scatter in Earth masses
    if rp <= maxE:
        mp = .44*rp**3 + .614*rp**4
        if not mean_only:
            mptmp = mp + np.random.randn() * rmsE
            while mptmp < 0:
                mptmp = mp + np.random.randn() * rmsE
            mp = mpstmp
    elif rp <= maxN:
        mp = 2.69 * rp**(.93)
        if not mean_only:
            mptmp = mp + np.random.randn() * rmsN
            while mptmp < 0:
                mptmp = mp + np.random.randn() * rmsN
            mp = mpstmp
    else:
        mp = rvs.kg2Mearth(4*np.pi/3 * rhoJ * 1e-3 * (rvs.Rearth2m(rp)*1e2)**3)
    return mp


def radF2mass(rps, Fs=336.5):
    '''Use relation from Weiss+Marcy 2014 and the relation for giant planets which includes 
    flux from Weiss+2013. rps in Rearth radii, Fs in Earth insolation'''
    # isolate different regimes
    rps, Fs = np.ascontiguousarray(rps), np.ascontiguousarray(Fs)
    assert rps.size == Fs.size
    Fs = Fs*1367*1e7*1e-4   # erg/s/cm^2
    rocky = rps < 1.5
    small = (1.5 <= rps) & (rps < 4)
    neptune = (4 <= rps) & (rps < 13.668)
    giant = rps >= 13.668

    # compute mean mass in each regime
    mps = np.zeros(rps.size)
    mps[rocky]   = .44*rps[rocky]**3 + .614*rps[rocky]**4
    mps[small]   = 2.69 * rps[small]**(.93)
    mps[neptune] = (rps[neptune]*Fs[neptune]**.03 / 1.78)**(1/.53)
    mps[giant]   = np.random.uniform(150,2e3)

    return mps


def MR_planets():
    '''Get the masses and radii from some exoplanets.'''
    names = np.array(['Venus','Earth','Uranus','Neptune','GJ 1132b','GJ 1214b','TRAPPIST-1d',
		      'TRAPPIST-1e','TRAPPIST-1f','TRAPPIST-1g','55 Cnc e','CoRoT-7b',
		      'HD 97658b','HIP 116454b','LHS 1140b','Kepler 10b','Kepler 10c',
		      'Kepler 36b','Kepler 78b','Kepler 93b'])
    mp = np.array([.814,1,14.54,17.15,1.85,6.55,.41,.62,.68,1.34,8.38,4.45,7.55,11.8,6.65,4.61,7.4,4.45,1.86,4.02])
    emp = np.array([0,0,0,0,.28,.98,.27,.58,.18,.88,.39,.99,.8,1.33,1.82,1.33,1.2,.3,.33,.68])
    rp = np.array([.949,1,4.007,3.88,1.16,2.678,.772,.918,1.045,1.127,2.08,1.61,2.247,2.59,1.43,1.48,2.35,1.49,1.17,1.48])
    erp = np.array([0,0,0,0,.11,.13,.03,.039,.038,.041,.16,.1,.097,.179,.1,.04,.07,.035,.13,.019])
    return names, mp, emp, rp, erp
