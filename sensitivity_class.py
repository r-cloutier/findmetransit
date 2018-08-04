from imports import *
from truncate_cmap import *


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

    def add_EBparams_guess(self, arr):
	self.EBparams_guess = arr
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



class Sensitivity_grid:

    def __init__(self, prefix):
        self.prefix = prefix
	self.fname_full = 'Results/%s'%prefix
    	self.pickleobject()

        
    def plot_sensitivity(self, Pgrid=True, pltt=True, label=False):
 	if Pgrid:
	    xarr, zarr = self.Pgrid, self.sensitivityP
	    xlabel = 'Period [days]'
	else:
            xarr, zarr = self.Fgrid, self.sensitivityF
	    xlabel = 'Insolation [S$_{\oplus}$]'

        fig = plt.figure(figsize=(5.5,4.7))
        ax = fig.add_subplot(111)
        cax = ax.pcolormesh(xarr, self.rpgrid, zarr.T,
                            cmap=truncate_colormap(plt.get_cmap('rainbow'),0,1),
                            vmin=0, vmax=1)
        cbar_axes = fig.add_axes([.08,.1,.84,.04])
        cbar = fig.colorbar(cax, cax=cbar_axes, orientation='horizontal')
        #cbar.ax.set_xticklabels(cticklabels)
        cbar.set_label('Detection Sensitivity', labelpad=.1)
	if Pgrid:
            ax.axvline(self.Pgrid.max()/2, ls='--', c='k')
            ax.axvline(self.Pgrid.max()/3, ls='--', c='k')
            ax.set_xlim((self.Pgrid.min(),self.Pgrid.max()))
	else:
	    ax.set_xlim((self.Fgrid.max(),self.Fgrid.min()))
        ax.set_ylim((self.rpgrid.min(),self.rpgrid.max()))
        ax.set_xscale('log'), ax.set_yscale('log')
        ax.set_xlabel(xlabel), ax.set_ylabel('r$_p$ [R$_{\oplus}$]')

	fig.subplots_adjust(bottom=.26, left=.14, right=.96, top=.97)
	if label:
	    try:
		os.mkdir('plots')
	    except OSError:
		pass
	    plt.savefig('plots/sensgrid_%s.png'%self.prefix)
        if pltt:
            plt.show()
        plt.close('all')

        
    def pickleobject(self):
        fObj = open('%s_SensitivityGrid'%self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close()
