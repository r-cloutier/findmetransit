from sensitivity_class import *
from linear_lnlike_testconstants import confirm_transits

global cols
cols = ['b','g','r']

def report_parameters(self):
    star = 'Stellar parameters:\nTmag     = %.3f\nMs       = %.3f\nRs       = %.3f\nTeff     = %i\nsig_phot = %.1f ppm'%(self.Tmag, self.Ms, self.Rs, self.Teff, self.ef[0]*1e6)
    Nplanets = self.params_true.shape[0]
    Plabels = ('%.3f, '*Nplanets)[:-2]%tuple([p for p in self.params_true[:,0]])
    rplabels = ('%.2f, '*Nplanets)[:-2]%tuple([r for r in self.rps])
    rpRslabels = ('%.5f, '*Nplanets)[:-2]%tuple([r for r in self.params_true[:,2]])
    rpRs2labels = ('%.6f, '*Nplanets)[:-2]%tuple([r**2*1e6 for r in self.params_true[:,2]])
    detlabels = ('%i, '*Nplanets)[:-2]%tuple([d for d in self.is_detected])
    planets = '\n\n%i injected planets:\nP  [days] = %s\nrp        = %s\nrp/Rs     = %s\n(rp/Rs)^2 = %s ppm\ndetected  = %s'%(Nplanets, Plabels, rplabels, rpRslabels, rpRs2labels,detlabels)
    NFPs = self.paramsFP_guess.shape[0]
    Plabels = ('%.3f, '*NFPs)[:-2]%tuple([p for p in self.paramsFP_guess[:,0]])
    FPs = '\n\n%i false positives:\nP [days] = %s'%(NFPs, Plabels)
    print star+planets+FPs
    

def planet_title(self):
    '''Define the planet title'''
    Nplanets = self.params_true.shape[0]
    title = ''
    for i in range(Nplanets):
        title += 'planet %i: %.3f days, %.2f Rearth\n'%(i+1, self.params_true[i,0], self.rps[i])
    return title[:-1]

    
def plot_raw_LC(self):
    '''Plot the raw light curve without any models removed.'''
    report_parameters(self)
    plt.errorbar(self.bjd, self.f, self.ef, fmt='k.', capsize=0,
                 elinewidth=.8)
    plt.plot(self.bjd, self.mu, 'g-')
    plt.xlabel('BJD'), plt.ylabel('Normalized flux')
    plt.xlim((self.bjd.min(),self.bjd.max()))
    plt.title(planet_title(self), fontsize=12)
    for i in range(self.params_true.shape[0]):
        P, T0 = self.params_true[i,:2]
        for k in range(-100,100):
            plt.axvline(T0+k*P, ymin=0, ymax=.05, c=cols[i], lw=2)
            plt.axvline(T0+k*P, ymin=.95, ymax=1, c=cols[i], lw=2)
    plt.show()
    plt.close('all')


def plot_corrected_LC(self):
    '''Plot the detrended light curve.'''
    report_parameters(self)
    plt.errorbar(self.bjd, self.fcorr, self.ef, fmt='k.', capsize=0,
                 elinewidth=.8)
    plt.xlabel('BJD'), plt.ylabel('Normalized flux')
    plt.xlim((self.bjd.min(),self.bjd.max()))
    plt.title(planet_title(self), fontsize=12)
    for i in range(self.params_true.shape[0]):
        P, T0 = self.params_true[i,:2]
        for k in range(-100,100):
            plt.axvline(T0+k*P, ymin=0, ymax=.05, c=cols[i], lw=2)
            plt.axvline(T0+k*P, ymin=.95, ymax=1, c=cols[i], lw=2)
    plt.show()
    plt.close('all')
    

def plot_transit_search(self, SNRthresh=5, log=False):
    '''Plot the SNR of the tophat model versus time each duration.'''
    report_parameters(self)
    Ndurations = self.SNRs_linearsearch.shape[1]
    Nplanets = self.params_true.shape[0]
    fig = plt.figure(figsize=(7.5,4+1.3*Ndurations))
    for i in range(Ndurations):
        ax = fig.add_subplot(Ndurations,1,i+1)
        yarr = np.log10(self.SNRs_linearsearch[:,i]) if log else self.SNRs_linearsearch[:,i]
        ax.plot(self.transit_times, yarr, 'k-')
        ax.axhline(SNRthresh, ls='--', c='k', lw=.9)
        ylabel = 'log$_{10}$(Transit S/N)' if log else 'Transit S/N'
        ax.set_ylabel(ylabel), ax.set_xlim((self.bjd.min(), self.bjd.max()))
        ax.set_title('duration = %.1f hours'%(self.durations[i]*24),
                     fontsize=12)
        if i < Ndurations-1:
            ax.set_xticklabels('')
        
        # plot planet markers
        for j in range(Nplanets):
            P, T0 = self.params_true[j,:2]
            for k in range(-100,100):
                plt.axvline(T0+k*P, ymin=0, ymax=.05, c=cols[j], lw=2)
                plt.axvline(T0+k*P, ymin=.95, ymax=1, c=cols[j], lw=2)
                
                
    ax.set_xlabel('Mid-transit time [BJD]')
    fig.subplots_adjust(bottom=.1, left=.1, right=.96, top=.96, hspace=.3)
    plt.show()
    plt.close('all')


def plot_phased_LC(self, P, T0, xlim=(.05,.05)):
    '''consider any planet ephemeris, partiucularly those of FPs.'''
    phase = foldAt(self.bjd, P, T0)
    phase[phase>.5] -= 1
    plt.figure(figsize=(12,4))
    plt.subplot(121)
    plt.errorbar(phase, self.fcorr, self.ef, fmt='k.', capsize=0,
		 elinewidth=.9)
    plt.xlabel('Phase [%.3f days]'%P), plt.xlim((-.5,.5))
    plt.ylabel('Normalized flux')
    plt.subplot(122)
    plt.errorbar(phase, self.fcorr, self.ef, fmt='k.', capsize=0,
                 elinewidth=.9)
    plt.xlabel('Phase [%.3f days]'%P), plt.xlim((-.05,.05))
    plt.show()
    plt.close('all')



def report_failed_planet_candidates(self):
    '''Report the potential planet candidates and if they failed certain 
    transit criteria.'''
    Nplanet_guess = self.params_guess_priorto_confirm.shape[0]
    for i in range(Nplanet_guess):
        print 'Period = %.3f days'%self.params_guess_priorto_confirm[i,0]
        print 'is median flux in-transit significantly deeper than out of transit?\n%s'%self.transit_condition_scatterin_gtr_scatterout[i] 
        print 'is transit depth S/N greater than the required threshold?\n%s\n'%self.transit_condition_depth_gtr_rms[i] 


def check_confirm_transits(self, dispersion_sig, depth_sig, bimodalfrac):
    '''Run confirm transits but with any arbitrary threshhold values.'''
    lnLOIs = self.lnLOIs[np.in1d(self.POIs, self.params_guess_priorto_confirm[:,0])]
    paramsout, lnLsout, transit_condition_scatterin_gtr_scatterout, transit_condition_depth_gtr_rms, transit_condition_no_bimodal_flux_intransit, transit_condition_orbitalP_fits_in_WF = confirm_transits(self.params_guess_priorto_confirm, lnLOIs, self.bjd, self.fcorr, self.ef, self.Ms, self.Rs, self.Teff, dispersion_sig, depth_sig, bimodalfrac)
    return self.params_guess_priorto_confirm, paramsout, lnLsout, transit_condition_scatterin_gtr_scatterout, transit_condition_depth_gtr_rms, transit_condition_no_bimodal_flux_intransit, transit_condition_orbitalP_fits_in_WF
