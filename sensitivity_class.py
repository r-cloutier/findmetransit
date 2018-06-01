from imports import *

def loadpickle(fname):
    fObj = open(fname, 'rb')
    self = pickle.load(fObj)
    fObj.close()
    return self


class Sensitivity:

    def __init__(self, fname):
	self.fname_full = 'Results/%s'%fname
        self.fname = fname
    	try:
            os.mkdir(self.fname_full)
        except OSError:
            pass
	self.pickleobject()


    def add_star(self, arr):
	self.Tmag, self.Rs, self.Ms, self.Teff = arr
	self.pickleobject()

    def add_params_true(self, arr):
	self.params_true = arr
        self.pickleobject()

    def add_EBparams_true(self, arr):
  	self.EBr1, self.EBr2, self.EBsbratio, self.EBinc, self.EBT0, self.EBP, self.EBq = arr
	self.pickleobject()

    def add_samplesGP(self, arr):
	self.samplesGP = arr
        self.pickleobject()

    def add_resultsGP(self, arr):
	self.resultsGP = arr
        self.pickleobject()

    def add_params_guess(self, arr):
	self.params_guess = arr
        self.pickleobject()

    def add_detected(self, arr):
	self.is_detected = arr
        self.pickleobject()

    def add_thetaGP(self, arr):
	self.thetaGP = arr
        self.pickleobject()

    def add_raw_timeseries(self, bjd, f, ef):
    	self.bjd = bjd
	self.f = f
	self.ef = ef
	self.pickleobject()

    def add_timeseries(self, mu, sig, fcorr):
	self.mu = mu
	self.sig = sig
	self.fcorr = fcorr
        self.pickleobject()

    def add_linearsearch(self, transit_times, durations, lnLs, depths):
	self.transit_times = transit_times
	self.durations = durations
	self.lnLs_linearsearch = lnLs
	self.depths_linearsearch = depths
        self.pickleobject()

    def add_paramguesses(self, Ps, T0s, Ds, Zs, lnLs_transit):
	self.Ps = Ps
	self.T0s = T0s
	self.Ds = Ds
	self.Zs = Zs
	self.lnLs_transit = lnLs_transit
        self.pickleobject()

    def add_OIs(self, POIs, T0OIs, DOIs, ZOIs, lnLOIs):
	self.POIs = POIs
	self.T0OIs = T0OIs
	self.DOIs = DOIs
	self.ZOIs = ZOIs
	self.lnLOIs = lnLOIs
        self.pickleobject()


    def pickleobject(self):
        fObj = open('%s/Sensitivity_class'%self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close()
