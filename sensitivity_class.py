from imports import *
from truncate_cmap import *
import rvs

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


    def plot_grids(self, Pgrid=True, pltt=True, label=False):
        '''Plot the various grids of interest including the number of planet 
        detections, FPs, sensitivity etc.'''
        if Pgrid:
	    xarr, xgrid, xsens = self.Ps, self.Pgrid, self.sensitivityP
	    xFP, xlabel, pnglabel = self.PsFP, 'Period [days]', 'P'
            vert1, vert2 = self.Pgrid.max()/2, self.Pgrid.max()/3
            xlim = self.Pgrid.min(), self.Pgrid.max()
        else:
            xarr, xgrid, xsens = self.Fs, self.Fgrid, self.sensitivityF
	    xFP, xlabel, pnglabel = self.FsFP, 'Insolation [S$_{\oplus}$]', 'F'
            vert1, vert2 = .9, .22
            xlim = self.Fgrid.max(), self.Fgrid.min()
        
        fig = plt.figure(figsize=(12,9))
        ax1 = fig.add_subplot(231)
        N_det, x, y = np.histogram2d(xarr[self.detected==1], self.rps[self.detected==1],
                                     bins=(xgrid, self.rpgrid))
        cax1 = ax1.pcolormesh(x, y, N_det.T, cmap=truncate_colormap(plt.get_cmap('rainbow'),0,1))
        cbar_axes1 = fig.add_axes([.06,.96,.27,.03])
        cbar1 = fig.colorbar(cax1, cax=cbar_axes1, orientation='horizontal')
        cbar1.ax.tick_params(labelsize=9)
        cbar1.set_label('Number of detected (injected) planets', fontsize=10, labelpad=.1)
        ax1.axvline(vert1, ls='--', c='k')
        ax1.axvline(vert2, ls='--', c='k')
        ax1.set_xlim(xlim), ax1.set_ylim((self.rpgrid.min(),self.rpgrid.max()))
        ax1.set_xscale('log'), ax1.set_yscale('log')
        ax1.set_ylabel('r$_p$ [R$_{\oplus}$]', fontsize=12)
        
        ax2 = fig.add_subplot(232)
        N_FP, x, y = np.histogram2d(xFP, self.rpsFP, bins=(xgrid, self.rpgrid))
        cax2 = ax2.pcolormesh(x, y, N_FP.T, cmap=truncate_colormap(plt.get_cmap('rainbow'),0,1))
        cbar_axes2 = fig.add_axes([.38,.96,.27,.03])
        cbar2 = fig.colorbar(cax2, cax=cbar_axes2, orientation='horizontal')
        cbar2.ax.tick_params(labelsize=9)
        cbar2.set_label('Number of false positives', fontsize=10, labelpad=.1)
        ax2.axvline(vert1, ls='--', c='k')
        ax2.axvline(vert2, ls='--', c='k')
        ax2.set_xlim(xlim), ax2.set_ylim((self.rpgrid.min(),self.rpgrid.max()))
        ax2.set_xscale('log'), ax2.set_yscale('log')
        
        ax3 = fig.add_subplot(233)
        yield_correction = 1 - N_FP / (N_det + N_FP)
        yield_correction[N_det+N_FP==0] = 1.
        cax3 = ax3.pcolormesh(xgrid, self.rpgrid, yield_correction.T,
                              cmap=truncate_colormap(plt.get_cmap('rainbow'),0,1))
        cbar_axes3 = fig.add_axes([.69,.96,.27,.03])
        cbar3 = fig.colorbar(cax3, cax=cbar_axes3, orientation='horizontal')
        cbar3.ax.tick_params(labelsize=9)
        cbar3.set_label('Multiplicative yield correction', fontsize=10, labelpad=.1)
        ax3.axvline(vert1, ls='--', c='k')
        ax3.axvline(vert2, ls='--', c='k')
        ax3.set_xlim(xlim), ax3.set_ylim((self.rpgrid.min(),self.rpgrid.max()))
        ax3.set_xscale('log'), ax3.set_yscale('log')
        
        ax4 = fig.add_subplot(234)
        N_det, x, y = np.histogram2d(xarr[self.detected==1], self.rps[self.detected==1],
                                     bins=(xgrid, self.rpgrid))
        cax4 = ax4.pcolormesh(x, y, (N_det*yield_correction).T,
                              cmap=truncate_colormap(plt.get_cmap('rainbow'),0,1))
        cbar_axes4 = fig.add_axes([.06,.1,.27,.03])
        cbar4 = fig.colorbar(cax4, cax=cbar_axes4, orientation='horizontal')
        cbar4.ax.tick_params(labelsize=9)
        cbar4.set_label('**Corrected # of detected (injected) planets', fontsize=10, labelpad=.1)
        ax4.axvline(vert1, ls='--', c='k')
        ax4.axvline(vert2, ls='--', c='k')
        ax4.set_xlim(xlim), ax4.set_ylim((self.rpgrid.min(),self.rpgrid.max()))
        ax4.set_xscale('log'), ax4.set_yscale('log')
        ax4.set_xlabel(xlabel, fontsize=12)
        ax4.set_ylabel('r$_p$ [R$_{\oplus}$]', fontsize=12)
        
        ax5 = fig.add_subplot(235)
        cax5 = ax5.pcolormesh(xgrid, self.rpgrid, xsens.T,
                              cmap=truncate_colormap(plt.get_cmap('rainbow'),0,1),
                              vmin=0, vmax=1)
        cbar_axes5 = fig.add_axes([.38,.1,.27,.03])
        cbar5 = fig.colorbar(cax5, cax=cbar_axes5, orientation='horizontal')
        cbar5.ax.tick_params(labelsize=9)
        cbar5.set_label('Detection sensitivity', fontsize=10, labelpad=.1)
        ax5.axvline(vert1, ls='--', c='k')
        ax5.axvline(vert2, ls='--', c='k')
        ax5.set_xlim(xlim), ax5.set_ylim((self.rpgrid.min(),self.rpgrid.max()))
        ax5.set_xscale('log'), ax5.set_yscale('log')
        ax5.set_xlabel(xlabel, fontsize=12)
        
        ax6 = fig.add_subplot(236)
        transitprob = get_median_transitprob(self, Pgrid=Pgrid)
        cax6 = ax6.pcolormesh(xgrid, self.rpgrid, (xsens*transitprob).T,
                              cmap=truncate_colormap(plt.get_cmap('rainbow'),0,1),
                              vmin=0)
        cbar_axes6 = fig.add_axes([.69,.1,.27,.03])
        cbar6 = fig.colorbar(cax6, cax=cbar_axes6, orientation='horizontal')
        cbar6.ax.tick_params(labelsize=9)
        cbar6.set_label('Detection completeness:\nsensitivity $\cdot$ med(transit prob)', fontsize=10, labelpad=1)
        ax6.axvline(vert1, ls='--', c='k')
        ax6.axvline(vert2, ls='--', c='k')
        ax6.set_xlim(xlim), ax6.set_ylim((self.rpgrid.min(),self.rpgrid.max()))
        ax6.set_xscale('log'), ax6.set_yscale('log')
        ax6.set_xlabel(xlabel, fontsize=12)

        fig.subplots_adjust(bottom=.22, left=.06, right=.96, top=.9, hspace=.2, wspace=.18)
        if label:
	    try:
	        os.mkdir('plots')
	    except OSError:
	        pass
            plt.savefig('plots/grids%s_%s.png'%(pnglabel, self.prefix))
        if pltt:
            plt.show()
        plt.close('all')
    
        
    def pickleobject(self):
        fObj = open('%s_SensitivityGrid'%self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close()


        

def _get_median_transitprob(self, Pgrid=True):
    if Pgrid:
        Pred = 10**(np.log10(self.Pgrid[1:])-np.diff(np.log10(self.Pgrid))[0]/2)
    else:
        xred = 10**(np.log10(self.Fgrid[1:])-np.diff(np.log10(self.Fgrid))[0]/2)
        Lmed = 3.828e26 * np.median(self.Rss)**2 * (np.median(self.Teffs)/5777)**4
        amed = rvs.m2AU(np.sqrt(Lmed/(4*np.pi*1367*xred)))
        Pred = rvs.period_sma(amed, np.median(self.Mss), 0)
                           
    rpred = 10**(np.log10(self.rpgrid[1:])-np.diff(np.log10(self.rpgrid))[0]/2)
    transitprob = np.zeros_like(self.sensitivityP)
    assert transitprob.shape == (Pred.size, rpred.size)
    for i in range(Pred.size):
        for j in range(rpred.size):
            transitprob[i,j] = (rvs.Rsun2m(np.median(self.Rss)) + rvs.Rearth2m(rpred[j])) / rvs.AU2m(rvs.semimajoraxis(Pred[i],np.median(self.Mss),0))
    return transitprob
