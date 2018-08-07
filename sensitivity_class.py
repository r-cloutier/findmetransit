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
	    xlabel, pnglabel = 'Period [days]', 'P'
	else:
            xarr, zarr = self.Fgrid, self.sensitivityF
	    xlabel, pnglabel = 'Insolation [S$_{\oplus}$]', 'F'

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
	    plt.savefig('plots/sensgrid%s_%s.png'%(pnglabel, self.prefix))
        if pltt:
            plt.show()
        plt.close('all')


    def plot_FP_correction(self, Pgrid=True, pltt=True, label=False):
        if Pgrid:
	    xarr, xgrid = self.Ps, self.Pgrid
	    xFP, xlabel, pnglabel = self.PsFP, 'Period [days]', 'P'
        else:
            xarr, xgrid = self.Fs, self.Fgrid
	    xFP, xlabel, pnglabel = self.FsFP, 'Insolation [S$_{\oplus}$]', 'F'

        # plot planet detections
        fig = plt.figure(figsize=(12,4))
        ax1 = fig.add_subplot(131)
        N_det, x, y = np.histogram2d(xarr[self.detected==1], self.rps[self.detected==1],
                                     bins=(xgrid, self.rpgrid))
        cax1 = ax1.pcolormesh(x, y, N_det.T, cmap=truncate_colormap(plt.get_cmap('rainbow'),0,1))
        cbar_axes1 = fig.add_axes([.07,.1,.27,.04])
        cbar1 = fig.colorbar(cax1, cax=cbar_axes1, orientation='horizontal')
        cbar1.set_label('Number of detected (simulated) planets', labelpad=.1)
        if Pgrid:
            ax1.axvline(self.Pgrid.max()/2, ls='--', c='k')
            ax1.axvline(self.Pgrid.max()/3, ls='--', c='k')
            ax1.set_xlim((self.Pgrid.min(),self.Pgrid.max()))
        else:
	    ax1.set_xlim((self.Fgrid.max(),self.Fgrid.min()))
        ax1.set_ylim((self.rpgrid.min(),self.rpgrid.max()))
        ax1.set_xscale('log'), ax1.set_yscale('log')
        ax1.set_xlabel(xlabel), ax1.set_ylabel('r$_p$ [R$_{\oplus}$]')

        # plot planet false positives
        ax2 = fig.add_subplot(132)
        N_FP, x, y = np.histogram2d(xFP, self.rpsFP, bins=(xgrid, self.rpgrid))
        cax2 = ax2.pcolormesh(x, y, np.log10(N_FP).T, cmap=truncate_colormap(plt.get_cmap('rainbow'),0,1))
        cbar_axes2 = fig.add_axes([.39,.1,.27,.04])
        cbar2 = fig.colorbar(cax2, cax=cbar_axes2, orientation='horizontal')
        cbar2.set_label('log$_{10}$ Number of false positives', labelpad=.1)
        if Pgrid:
            ax2.axvline(self.Pgrid.max()/2, ls='--', c='k')
            ax2.axvline(self.Pgrid.max()/3, ls='--', c='k')
            ax2.set_xlim((self.Pgrid.min(),self.Pgrid.max()))
        else:
	    ax2.set_xlim((self.Fgrid.max(),self.Fgrid.min()))
        ax2.set_ylim((self.rpgrid.min(),self.rpgrid.max()))
        ax2.set_xscale('log'), ax2.set_yscale('log'), ax2.set_xlabel(xlabel)

        # plot yield correction (i.e. reduce the yield because of the fraction of FPs): yield(P,rp) *= 1 - N_FP(P,rp) / (N_detected(P,rp) + N_FP(P,rp))
        ax3 = fig.add_subplot(133)
        yield_correction = 1 - N_FP / (N_det + N_FP)
        yield_correction[N_det+N_FP==0] = 1.
        cax3 = ax3.pcolormesh(xgrid, self.rpgrid, yield_correction.T,
                              cmap=truncate_colormap(plt.get_cmap('rainbow'),0,1))
        cbar_axes3 = fig.add_axes([.71,.1,.27,.04])
        cbar3 = fig.colorbar(cax3, cax=cbar_axes3, orientation='horizontal')
        cbar3.set_label('Multiplicative yield correction', labelpad=.1)
        if Pgrid:
            ax3.axvline(self.Pgrid.max()/2, ls='--', c='k')
            ax3.axvline(self.Pgrid.max()/3, ls='--', c='k')
            ax3.set_xlim((self.Pgrid.min(),self.Pgrid.max()))
        else:
	    ax3.set_xlim((self.Fgrid.max(),self.Fgrid.min()))
        ax3.set_ylim((self.rpgrid.min(),self.rpgrid.max()))
        ax3.set_xscale('log'), ax3.set_yscale('log'), ax3.set_xlabel(xlabel)
    
        fig.subplots_adjust(bottom=.29, left=.06, right=.99, top=.97, wspace=.12)
        if label:
	    try:
	        os.mkdir('plots')
	    except OSError:
	        pass
            plt.savefig('plots/FPgrid%s_%s.png'%(pnglabel, self.prefix))
        if pltt:
            plt.show()
        plt.close('all')
    
        
    def pickleobject(self):
        fObj = open('%s_SensitivityGrid'%self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close()
