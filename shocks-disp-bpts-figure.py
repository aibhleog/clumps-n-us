'''

Velocity dispersion, velc

'''

from fitting_ifu_spectra import *
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors as co
import matplotlib.patheffects as PathEffects
import warnings
warnings.filterwarnings("ignore")
saveit = True


# # one galaxy for now
# target = 'SunburstArc-P2P3'


# reading in latest redshifts & setting target order to descending z
redshifts = pd.read_csv('plots-data/all_redshifts.global.txt',sep=r'\s+')
redshifts = redshifts.loc[redshifts.target!='SGAS2111'].copy()
target_order = redshifts.sort_values(by='z',ascending=False).target.values.copy().tolist()
# dropping SGAS2111

# making 2 for Sunburst
target_order[-1] += '-P1'
target_order.append(target_order[-1][:-1]+'2P3')


# looking at velocity offset from systemic and dispersion
def fwhmkms2sig(lam): return 1/2.35/3e05*((1+z)*lam) # units matching input ref wavelength
def sig2vdisp(lam): return 3e05/((1+z)*lam) # units matching input ref wavelength

colors = get_colors('cmr.guppy',6,ymax=1)





# ----------------------------------
# comparison with SDSS MaNGA eLINER
# ----------------------------------
def get_colors(cmap,size,ymax=0.8):
    cmap = plt.get_cmap(cmap)
    colors = [cmap(j) for j in np.linspace(0,ymax,size)]
    return colors

from scipy.stats import gaussian_kde
# --- extensions wanted:  EMLINE_GFLUX, EMLINE_GFLUX_IVAR
# --- slices wanted: C14[Hb], C16[O3], C20[O1], C23[Ha], C24[N2], C27+C26[S2] // 1e-17 cgs
# --- extensions wanted:  EMLINE_GSIGMA, EMLINE_GSIGMA_IVAR
# --- slices wanted: C14[Hb], C16[O3], C20[O1], C23[Ha], C24[N2], C27+C26[S2] // km/s
catpath = '/Users/tahutch1/data/catalogs/manga/'
manga = fits.open(catpath+'manga-7975-3704-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz')

# making BPT axes
manga_n2 = manga['EMLINE_GFLUX'].data[24].copy() # NII6585
manga_s2 = manga['EMLINE_GFLUX'].data[25].copy()+manga['EMLINE_GFLUX'].data[26].copy() # S2tot
manga_o1 = manga['EMLINE_GFLUX'].data[20].copy() # OI6300

manga_n2ha = np.log10(manga_n2 / manga['EMLINE_GFLUX'].data[23].copy()) # N2 / Ha
manga_s2ha = np.log10(manga_s2 / manga['EMLINE_GFLUX'].data[23].copy()) # S2tot / Ha
manga_o1ha = np.log10(manga_o1 / manga['EMLINE_GFLUX'].data[23].copy()) # O1 / Ha

manga_o3hb = np.log10(manga['EMLINE_GFLUX'].data[16].copy() /\
                      manga['EMLINE_GFLUX'].data[14].copy()) # O3 / Hb

# velocity dispersion
manga_y = manga['EMLINE_GSIGMA'].data[23].flatten() # Ha sigma, flattened to 1D
# ---------------------------------------










veloffdisp = {}

# THE FIGURE


plt.figure(figsize=(9.5,10)) # part 1
gs0 = gridspec.GridSpec(4,1,height_ratios=np.ones(4),hspace=0.02) # part 1

# plt.figure(figsize=(9.5,10*(3/4))) # part 2
# gs0 = gridspec.GridSpec(3,1,height_ratios=np.ones(3),hspace=0.02) # part 2

broad_flag = False
targets = target_order[:4] # part 1
# targets = target_order[4:] # part 2

end = len(targets) - 1
for i,target in enumerate(targets):

    # checking if broad comp
    if target == 'SGAS1050' or target[:3] == 'Sun':
        broad_flag = True
        broad = fits.open(f'plots-data/{target}-broadnarrow-disp-bpt-shocks.fits')
    
    # making sub_axes
    gs = gridspec.GridSpecFromSubplotSpec(2,4,subplot_spec=gs0[i], 
                                          width_ratios=[1,1,1,0.3],wspace=0,
                                          height_ratios=[0.3,1],hspace=0)
    c = target_order.index(target)
    # reading in data
    if 'Sun' in target: grating = 'g235h'; c = 5
    else: grating = 'g235m'
    galaxy, oldz, bigpath, grating = get_galaxy_info(target,grating)
    path = f'plots-data/cubespec-linefits/{target}/'
    
    z = redshifts.loc[redshifts.target==target.split('-')[0],'z'].values[0]
    zerr = redshifts.loc[redshifts.target==target.split('-')[0],'zerr'].values[0]

    # --- CODE PULLED FROM shocks-diagnostics.py
    # dictionary to house all the line  maps
    fluxes = {}
    
    # reading in maps for the main lines fit in each function
    # (i.e., if there were multiple lines fit, it's just the first for now)
    for key in ['o3','o1','ha','s2']:
        line = key
        if key == 'o3': line = 'o3hb'
        if (target == 'SGAS1050' or target == 'SunburstArc-P1') and key == 'o1': 
            line = 'o1s3'
            
        l = fits.getdata(path+f'{target}-{line}-params-stdpix.fits')[1]
        lerr = fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=2)[1]
        norm = fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=3)
    
        l,lerr = l*norm, lerr*norm
    
        # # doing a S/N cut first
        # l[l/lerr<2] = np.nan
    
        if key == 's2': key = 's2a'
        fluxes[key] = [l.copy(),lerr.copy()]
    
    
        if key == 'ha':
            header = fits.getheader(path+f'{target}-{line}-params-stdpix.fits')
            zs = [fits.getdata(path+f'{target}-{line}-params-stdpix.fits')[0],
                  fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=2)[0]]
            sig = [fits.getdata(path+f'{target}-{line}-params-stdpix.fits')[2],
                   fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=2)[2]]
    
            zs[0][zs[0]==np.float64(0.0)] = np.nan
            sig[0][sig[0]==np.float64(0.0)] = np.nan



    
    # now adding the more complicated line fluxes
    # ---------------------
    # --- HBETA
    key,line,indx = 'hb','o3hb',3 # O3Hb ratio
    l = fits.getdata(path+f'{target}-{line}-params-stdpix.fits')[indx]
    lerr = fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=2)[indx]
    norm = fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=3)
    new_l = fluxes['o3'][0].copy() / l.copy() # O3 / O3Hb = Hb
    new_lerr = new_l * np.sqrt( (fluxes['o3'][1].copy()/fluxes['o3'][0].copy())**2 +\
                                (lerr.copy()/l.copy())**2 )
    fluxes[key] = [new_l.copy(),new_lerr.copy()]
    
    # --- NII
    key,line,indx = 'n2','ha',3 # NIIHa ratio
    l = fits.getdata(path+f'{target}-{line}-params-stdpix.fits')[indx]
    lerr = fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=2)[indx]
    norm = fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=3)
    new_l = fluxes['ha'][0].copy() * l.copy() # Ha * NIIHa = NIIb
    new_lerr = new_l * np.sqrt( (fluxes['ha'][1].copy()/fluxes['ha'][0].copy())**2 +\
                                (lerr.copy()/l.copy())**2 )
    fluxes[key] = [new_l.copy(),new_lerr.copy()]
    
    # --- S2 b
    key,line,indx = 's2b','s2',2 # S2 ratio
    l = fits.getdata(path+f'{target}-{line}-params-stdpix.fits')[indx]
    lerr = fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=2)[indx]
    norm = fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=3)
    new_l = fluxes['s2a'][0].copy() / l.copy() # S2a / S2 ratio = S2b
    new_lerr = new_l * np.sqrt( (fluxes['s2a'][1].copy()/fluxes['s2a'][0].copy())**2 +\
                                (lerr.copy()/l.copy())**2 )
    fluxes[key] = [new_l.copy(),new_lerr.copy()]
    
    # --- COMBINED S2
    fluxes['s2'] = [fluxes['s2a'][0].copy() + fluxes['s2b'][0].copy(),
                    np.sqrt(fluxes['s2a'][1].copy()**2 + fluxes['s2b'][1].copy()**2)]
    
    # running through and doing S/N > 2
    for key in fluxes:
        snr = fluxes[key][0].copy() / fluxes[key][1].copy()
        fluxes[key][0][snr<2] = np.nan
    
    # ---------------------------------------------------
    # --- END OF CODE PULLED FROM shocks-diagnostics.py


    # PREPPING LOG BPT RATIOS
    prelog_n2ha = fluxes['n2'][0].copy() / fluxes['ha'][0].copy()
    prelog_n2haerr = prelog_n2ha * np.sqrt((fluxes['n2'][1].copy()/fluxes['n2'][0].copy())**2 + (fluxes['ha'][1].copy()/fluxes['ha'][0].copy())**2)
    n2ha = np.log10(prelog_n2ha.copy())
    n2haerr = (prelog_n2haerr/prelog_n2ha) * 0.434
    # for SGAS1050, removing severe outliers
    n2ha[n2ha<-3] = np.nan
    n2ha[n2ha>1] = np.nan
    n2ha[n2haerr==0] = np.nan

    prelog_s2ha = fluxes['s2'][0].copy() / fluxes['ha'][0].copy()
    prelog_s2haerr = prelog_s2ha * np.sqrt((fluxes['s2'][1].copy()/fluxes['s2'][0].copy())**2 + (fluxes['ha'][1].copy()/fluxes['ha'][0].copy())**2)
    s2ha = np.log10(prelog_s2ha.copy())
    s2haerr = (prelog_s2haerr/prelog_s2ha) * 0.434
    # for SGAS1050, removing severe outliers
    s2ha[s2ha<-3] = np.nan
    s2ha[s2ha>1] = np.nan
    s2ha[s2haerr==0] = np.nan

    prelog_o1ha = fluxes['o1'][0].copy() / fluxes['ha'][0].copy()
    prelog_o1haerr = prelog_o1ha * np.sqrt((fluxes['o1'][1].copy()/fluxes['o1'][0].copy())**2 + (fluxes['ha'][1].copy()/fluxes['ha'][0].copy())**2)
    o1ha = np.log10(prelog_o1ha.copy())
    o1haerr = (prelog_o1haerr/prelog_o1ha) * 0.434
    # for SGAS1050, removing severe outliers
    o1ha[o1ha<-3] = np.nan
    o1ha[o1ha>1] = np.nan
    o1ha[o1haerr==0] = np.nan


    # global bins for x axes
    bins = get_global_bins([n2ha[np.isfinite(n2ha)].copy(),
                            s2ha[np.isfinite(s2ha)].copy(),
                            o1ha[np.isfinite(o1ha)].copy()],bins=40)
    

    # velocity dispersion
    velocity_dispersion = sig[0] * sig2vdisp(.6564) # km/s
    velocity_dispersionerr = sig[1] * sig2vdisp(.6564) # km/s
    velocity_dispersion = SNRCUT(velocity_dispersion,velocity_dispersionerr,5)
    velocity_dispersion = NOZERO(velocity_dispersion.copy())


    # additional scatter kwargs
    ec = dict({'edgecolor':'k','linewidth':0.25})
    # levels for manga kde
    levels = [0.003,0.008,0.015,0.1] # sigma vs N2Ha
    level_colors = get_colors('cmr.freeze',len(levels)+2)[2:]
    
    
    
    # --- N2HA
    ax = plt.subplot(gs[4])
    ax.scatter(n2ha.flatten(),velocity_dispersion.flatten(),
               s=10,alpha=0.7,color=colors[c],**ec)
    ax.errorbar(n2ha.flatten(),velocity_dispersion.flatten(),
                yerr=velocity_dispersionerr.flatten(),xerr=n2haerr.flatten(),
                color='none',ms=14,alpha=0.3,
                ecolor='k',fmt='.',mew=2,zorder=0,lw=1)
    

    # median 
    ax.axvline(np.nanmedian(n2ha),ls='--',color=colors[c],lw=0.75,alpha=0.5,zorder=0)
    ax.axhline(np.nanmedian(velocity_dispersion),ls='--',color=colors[c],lw=0.75,alpha=0.5,zorder=0)
    
    ax.set_ylabel(r'$\sigma$ [km/s]')
    if i == end: ax.set_xlabel(r'log$_{10}$ [NII] / H$\alpha$')
    else: ax.set_xticklabels([])

    ax.set_xlim(-2.5,0.1)
    ax.set_ylim(20,290)
    xlims = np.asarray(ax.get_xlim())
    ylims = np.asarray(ax.get_ylim())

    ax.text(0.05,0.83,target,transform=ax.transAxes,fontsize=11)

    # IF BROAD COMPONENT
    if broad_flag:
        ax.scatter(broad[4].data[0].flatten(),broad[2].data[0].flatten(),zorder=0,
                   s=10,alpha=0.3,color='none',edgecolor=colors[c],lw=1)
        ax.errorbar(broad[4].data[0].flatten(),broad[2].data[0].flatten(),
                yerr=broad[2].data[1].flatten(),xerr=broad[4].data[1].flatten(),
                color='none',ms=14,alpha=0.1,
                ecolor='k',fmt='.',mew=2,zorder=-1,lw=0.5)
    

    # ----------------------------------
    # ----------------------------------
    x,y = manga_n2ha.copy().flatten(), manga_y.copy().flatten()
    # remvoing nans & infs
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask].copy()
    y = y[mask].copy()
    
    # mkaing the kde & plotting
    nbins = 50
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    k = gaussian_kde(np.vstack([x,y]))
    zval = k(np.vstack([xi.flatten(), yi.flatten()]))
    
    cf = ax.contourf(xi,yi,zval.reshape(xi.shape),zorder=0,alpha=0.75,
            levels=levels,colors=level_colors)
    # ----------------------------------
    # ----------------------------------


    # adding histogram of values
    ax_hist0 = plt.subplot(gs[0]); ax_hist0.axis('off')
    __ = ax_hist0.hist(n2ha[np.isfinite(n2ha)].flatten(),bins=bins,
                      color=colors[c],alpha=0.5)
    __ = ax_hist0.hist(n2ha[np.isfinite(n2ha)].flatten(),bins=bins,
                      histtype='step',color='k',lw=0.5,alpha=0.5)
    hist_ymax = ax_hist0.get_ylim()[1]
    if hist_ymax < np.nanmax(s2ha[np.isfinite(s2ha)]): 
        hist_ymax = np.nanmax(s2ha[np.isfinite(s2ha)])
    ax_hist0.set_ylim(0,hist_ymax)
    ax_hist0.set_xlim(xlims)

    # median 
    ax_hist0.axvline(np.nanmedian(n2ha),ls='--',color=colors[c],lw=0.75,alpha=0.5,zorder=0)

    # manga hist
    ax_hist05 = ax_hist0.twinx(); ax_hist05.axis('off')
    __ = ax_hist05.hist(x,bins=20,color=level_colors[0],alpha=0.35)
    ax_hist05.set_xlim(xlims)

    
    # IF BROAD COMPONENT
    if broad_flag:
        broad_n2ha = broad[4].data[0].copy()
        __ = ax_hist0.hist((broad_n2ha[np.isfinite(broad_n2ha)]).flatten(),
                           bins=bins,histtype='step',
                           color=colors[c],lw=1.5,alpha=0.75)
    

    
    # --- S2Ha    
    ax = plt.subplot(gs[5])
    ax.scatter(s2ha.flatten(),velocity_dispersion.flatten(),
               s=10,alpha=0.7,color=colors[c],**ec)
    ax.errorbar(s2ha.flatten(),velocity_dispersion.flatten(),
                yerr=velocity_dispersionerr.flatten(),xerr=s2haerr.flatten(),
                color='none',ms=14,alpha=0.3,
                ecolor='k',fmt='.',mew=2,zorder=0,lw=1)

    # marking the shockly likely region
    ax.axvspan(-0.5,1,alpha=0.4,color='gray',zorder=0)

    # median 
    ax.axvline(np.nanmedian(s2ha),ls='--',color=colors[c],lw=0.75,alpha=0.5,zorder=0)
    ax.axhline(np.nanmedian(velocity_dispersion),ls='--',color=colors[c],lw=0.75,alpha=0.5,zorder=0)
    
    if i == end: ax.set_xlabel(r'log$_{10}$ [SII] / H$\alpha$')
    else: ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    

    # ----------------------------------
    # ----------------------------------
    x,y = manga_s2ha.copy().flatten(), manga_y.copy().flatten()
    # remvoing nans & infs
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask].copy()
    y = y[mask].copy()
    
    # mkaing the kde & plotting
    nbins = 50
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    k = gaussian_kde(np.vstack([x,y]))
    zval = k(np.vstack([xi.flatten(), yi.flatten()]))
    
    cf = ax.contourf(xi,yi,zval.reshape(xi.shape),zorder=0,alpha=0.75,
            levels=levels,colors=level_colors)
    # ----------------------------------
    # ----------------------------------


    # IF BROAD COMPONENT
    if broad_flag:
        ax.scatter(broad[6].data[0].flatten(),broad[2].data[0].flatten(),zorder=0,
                   s=10,alpha=0.3,color='none',edgecolor=colors[c],lw=1)
        ax.errorbar(broad[6].data[0].flatten(),broad[2].data[0].flatten(),
                yerr=broad[2].data[1].flatten(),xerr=broad[6].data[1].flatten(),
                color='none',ms=14,alpha=0.1,
                ecolor='k',fmt='.',mew=2,zorder=-1,lw=0.5)
    

    # adding histogram of values
    ax_hist1 = plt.subplot(gs[1]); ax_hist1.axis('off')
    __ = ax_hist1.hist(s2ha[np.isfinite(s2ha)].flatten(),bins=bins,
                      color=colors[c],alpha=0.5)
    __ = ax_hist1.hist(s2ha[np.isfinite(s2ha)].flatten(),bins=bins,
                      histtype='step',color='k',lw=0.5,alpha=0.5)
    ax_hist1.set_ylim(0,hist_ymax)
    ax_hist1.set_xlim(xlims)

    # median 
    ax_hist1.axvline(np.nanmedian(s2ha),ls='--',color=colors[c],lw=0.75,alpha=0.5,zorder=0)
    
    # manga hist
    ax_hist15 = ax_hist1.twinx(); ax_hist15.axis('off')
    __ = ax_hist15.hist(x,bins=20,color=level_colors[0],alpha=0.35)
    ax_hist15.set_xlim(xlims)


    # IF BROAD COMPONENT
    if broad_flag:
        broad_s2ha = broad[6].data[0].copy()
        __ = ax_hist15.hist((broad_s2ha[np.isfinite(broad_s2ha)]).flatten(),
                           bins=bins,histtype='step',
                           color=colors[c],lw=1.5,alpha=0.75)

    
    
    
    # --- O1 / Ha
    ax = plt.subplot(gs[6])
    ax.scatter(o1ha.flatten(),velocity_dispersion.flatten(),
               s=10,alpha=0.7,color=colors[c],**ec)
    ax.errorbar(o1ha.flatten(),velocity_dispersion.flatten(),
                yerr=velocity_dispersionerr.flatten(),xerr=o1haerr.flatten(),
                color='none',ms=14,alpha=0.3,
                ecolor='k',fmt='.',mew=2,zorder=0,lw=1)
    
    # marking the shockly likely region
    ax.axvspan(-1.2,1,alpha=0.4,color='gray',zorder=0)

    # median 
    ax.axvline(np.nanmedian(o1ha),ls='--',color=colors[c],lw=0.75,alpha=0.5,zorder=0)
    ax.axhline(np.nanmedian(velocity_dispersion),ls='--',color=colors[c],lw=0.75,alpha=0.5,zorder=0)
    
    if i == end: ax.set_xlabel(r'log$_{10}$ [OI] / H$\alpha$')
    else: ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    if i == 0:
        txt = ax.text(-0.45,50,'MaNGA\neLIER',ha='center',color='#BBF2FC',#'w',
                fontweight='bold',fontsize=6)
        txt.set_path_effects([PathEffects.withStroke(linewidth=1.5,
                                                     foreground=level_colors[0])])

    # ----------------------------------
    # ----------------------------------
    x,y = manga_o1ha.copy().flatten(), manga_y.copy().flatten()
    # remvoing nans & infs
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask].copy()
    y = y[mask].copy()
    
    # mkaing the kde & plotting
    nbins = 50
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    k = gaussian_kde(np.vstack([x,y]))
    zval = k(np.vstack([xi.flatten(), yi.flatten()]))
    
    cf = ax.contourf(xi,yi,zval.reshape(xi.shape),zorder=0,alpha=0.75,
            levels=levels,colors=level_colors)
    # ----------------------------------
    # ----------------------------------

    
    # adding histogram of values
    ax_hist2 = plt.subplot(gs[2]); ax_hist2.axis('off')
    __ = ax_hist2.hist(o1ha[np.isfinite(o1ha)].flatten(),bins=bins,
                      color=colors[c],alpha=0.5)
    __ = ax_hist2.hist(o1ha[np.isfinite(o1ha)].flatten(),bins=bins,
                      histtype='step',color='k',lw=0.5,alpha=0.5)
    ax_hist2.set_ylim(0,hist_ymax)
    ax_hist2.set_xlim(xlims)
    # median 
    ax_hist2.axvline(np.nanmedian(o1ha),ls='--',color=colors[c],lw=0.75,alpha=0.5,zorder=0)

    # manga hist
    ax_hist25 = ax_hist2.twinx(); ax_hist25.axis('off')
    __ = ax_hist25.hist(x,bins=30,color=level_colors[0],alpha=0.35)
    ax_hist25.set_xlim(xlims)
    


    # --- velocity dispersion histogram
    ax = plt.subplot(gs[7]); ax.axis('off')
    ax.set_ylim(ylims)
    # ax.set_xscale('log')
    # median 
    ax.axhline(np.nanmedian(velocity_dispersion),ls='--',color=colors[c],lw=0.75,alpha=0.5,zorder=0)

    # manga hist
    ax5 = ax.twiny(); ax5.axis('off')
    __ = ax5.hist(y,bins=30,color=level_colors[0],alpha=0.35,
                  zorder=0,orientation='horizontal',label='MaNGA')
    ax5.set_ylim(ylims)
    # ax5.set_xscale('log')

    __ = ax.hist(velocity_dispersion.flatten(),bins=30,orientation='horizontal',
                 alpha=0.5,color=colors[c],range=(ylims),label='simple gaussian')
    __ = ax.hist(velocity_dispersion.flatten(),bins=30,orientation='horizontal',
                      histtype='step',color='k',lw=0.5,alpha=0.5,range=(ylims))

    # IF BROAD COMPONENT
    if broad_flag:
        if target == 'SGAS1050': bbins = 40
        else: bbins = 30
        broad_disp = broad[2].data[0].copy()
        __ = ax.hist((broad_disp[np.isfinite(broad_disp)]).flatten(),
                           bins=bbins,histtype='step',orientation='horizontal',
                           color=colors[c],lw=1.5,alpha=0.75,label='additional broad\ncomponent')


    # legend
    handles,labels = ax.get_legend_handles_labels()
    mhandles,mlabels = ax5.get_legend_handles_labels()
    ax.legend(np.concatenate([handles,mhandles]),np.concatenate([labels,mlabels]),
             handletextpad=0.5,handlelength=1,fontsize=10,loc=2,
             bbox_to_anchor=(0,1.2))
                             

        


    broad_flag = False


plt.savefig('plots-data/shocks-family-vdispBPT-part1.pdf')
# plt.savefig('plots-data/shocks-family-vdispBPT-part2.pdf')
plt.show()
plt.close('all')






# # getting ha+n2 broad velocity fits
# ha = fits.getdata(path+f'{target}-hadouble-velocity-params-stdpix.fits')
# haerr = fits.getdata(path+f'{target}-hadouble-velocity-params-stdpix.fits',ext=2)
# hanorm = fits.getdata(path+f'{target}-hadouble-velocity-params-stdpix.fits',ext=3)

# zs = [ha[0].copy(), haerr[0].copy()] # narrow comp
# zs_broad = [ha[1].copy(), haerr[1].copy()]
# sig = [ha[3].copy(), haerr[3].copy()]

# ha_narrow = [ha[2].copy(),haerr[2].copy()]
# ha_broad = [ha[2].copy()*ha[4].copy(),
#             (ha[2].copy()*ha[4].copy()) *\
#             np.sqrt( (haerr[2].copy()/ha[2].copy())**2 +\
#                      (haerr[4].copy()/ha[4].copy())**2 )]






'''

# clump & global measurements
integrated = pd.read_csv('plots-data/globally-integrated/ALL_fluxes_PLUS-TEMDENEBV.csv')
integ_columns = np.array(integrated.columns.tolist())
measurement = integrated.label.values
restwave = integrated.restwave.values


# -----------------------------
# pulling clump & global info!
# -----------------------------
integ_indx = [target in s for s in integ_columns] # bool
integ_targ = integrated[integ_columns[integ_indx]].copy()
global_targ = integ_targ[[f'{target}_STACK_flux',
                          f'{target}_STACK_fluxerr']].copy() # first two columns
clumps_targ = integ_targ[integ_columns[integ_indx][2:]].copy() # the rest
# ----------------------------------------------------------






'''










