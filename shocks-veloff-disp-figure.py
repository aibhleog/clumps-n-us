'''

Velocity dispersion, velc

'''

from fitting_ifu_spectra import *
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors as co
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




# -- to do north up east left
PA_SGAS1527 = 138.5 + 25 # deg
PA_SGAS1050 = 138.5 - 340 # deg
PA_CosmicEye = 138.5 - 241 - 15 # deg
PA_SGAS1429 = 138.5 + 35 # deg
# inside imshow & contour, add this kwarg: transform=transform+ax.transData


lims = {'SGAS1050':[(-105,-8),(-33,33)],
        'CosmicEye':[(-25,25),(-50,-5)],
        'SGAS1429':[(-45,0),(-37,3)],
        'SGAS1527':[(-50,0),(-32,9)],
        'SGAS1110':[None,(5,72)]}

xshift = {'CosmicEye':-0.1,'SGAS1429':-0.1,'SGAS1110':0,'SunburstArc-P1':-0.075,
          'SGAS1050':0,'SGAS1527':-0.025,'SunburstArc-P2P3':-0.025}


veloffdisp = {}

# THE FIGURE

plt.figure(figsize=(10,14))
# gs = gridspec.GridSpec(7,3,width_ratios=[1,1,1],#hspace=0,
gs = gridspec.GridSpec(7,4,width_ratios=[1,1,1,1],hspace=0,wspace=0,
                       height_ratios=[1,1,1,1,1.5,1,1.5])


for i,target in enumerate(target_order):

    transform_kwarg = dict()
    transform = plt.matplotlib.transforms.Affine2D()
    if target in ['CosmicEye','SGAS1050','SGAS1527','SGAS1429']: 
        transform.rotate_deg(eval(f'PA_{target}'))

    

    # reading in data
    if 'Sun' in target: grating = 'g235h'
    else: grating = 'g235m'
    galaxy, oldz, bigpath, grating = get_galaxy_info(target,grating)
    path = f'plots-data/cubespec-linefits/{target}/'
    
    z = redshifts.loc[redshifts.target==target.split('-')[0],'z'].values[0]
    zerr = redshifts.loc[redshifts.target==target.split('-')[0],'zerr'].values[0]
    
    fluxes = {}
    
    # --- Ha+NII
    han2 = fits.getdata(path+f'{target}-ha-params-stdpix.fits')
    han2err = fits.getdata(path+f'{target}-ha-params-stdpix.fits',ext=2)
    han2norm = fits.getdata(path+f'{target}-ha-params-stdpix.fits',ext=3)
    fluxes['ha'] = [han2[1].copy()*han2norm, han2err[1].copy()*han2norm]
    
    # --- O1
    key,line,indx = 'o1','o1',1 # O1
    if (target == 'SGAS1050' or target == 'SunburstArc-P1') and key == 'o1': 
            line = 'o1s3'
    l = fits.getdata(path+f'{target}-{line}-params-stdpix.fits')[indx]
    lerr = fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=2)[indx]
    norm = fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=3)
    fluxes[key] = [SNRCUT(l.copy(),lerr.copy(),2.5),lerr.copy()]
    
    
    # --- S2a
    key,line,indx = 's2a','s2',1 # S2 ratio
    l = fits.getdata(path+f'{target}-{line}-params-stdpix.fits')[indx]
    lerr = fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=2)[indx]
    norm = fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=3)
    fluxes[key] = [SNRCUT(l.copy(),lerr.copy(),2.5),lerr.copy()]
    
    # --- S2b
    key,line,indx = 's2b','s2',2 # S2 ratio
    l = fits.getdata(path+f'{target}-{line}-params-stdpix.fits')[indx]
    lerr = fits.getdata(path+f'{target}-{line}-params-stdpix.fits',ext=2)[indx]
    new_l = fluxes['s2a'][0].copy() / l.copy() # S2a / S2 ratio = S2b
    new_lerr = new_l * np.sqrt( (fluxes['s2a'][1].copy()/fluxes['s2a'][0].copy())**2 +\
                                (lerr.copy()/l.copy())**2 )
    fluxes[key] = [new_l.copy(),new_lerr.copy()]
    
    # --- COMBINED S2
    fluxes['s2'] = [fluxes['s2a'][0].copy() + fluxes['s2b'][0].copy(),
                    np.sqrt(fluxes['s2a'][1].copy()**2 + fluxes['s2b'][1].copy()**2)]
    
    
    # setting up other values
    zs = [han2[0].copy(), han2err[0].copy()] # narrow comp
    sig = [han2[2].copy(), han2err[2].copy()]
    
          
    
    # masking out pixels with S/N<3 from Ha
    zs[0][(fluxes['ha'][0]/fluxes['ha'][1])<5] = np.nan
    sig[0][(fluxes['ha'][0]/fluxes['ha'][1])<5] = np.nan
    
    velocity_offset = veloff(z,NOZERO(zs[0])) # km/s
    velocity_lower = CAPIT(CAPIT(NOZERO(veloff(z,zs[0]-zs[1])),1e3),-1e3,sign='<')
    velocity_upper = CAPIT(CAPIT(NOZERO(veloff(z,zs[0]+zs[1])),1e3),-1e3,sign='<')
    velocity_offseterr = abs(velocity_upper-velocity_lower)
    
    velocity_dispersion = sig[0] * sig2vdisp(.6564) # km/s
    velocity_dispersionerr = sig[1] * sig2vdisp(.6564) # km/s
    
    # doing a S/N cut for velocity dispersion, specifically
    velocity_dispersion = SNRCUT(velocity_dispersion,velocity_dispersionerr,3)
    velocity_dispersion = NOZERO(velocity_dispersion.copy())
    
    
    # --- VELOCITY OFFSET
    ax = plt.subplot(gs[0+i*4]); ax.axis('off')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    im = ax.imshow(velocity_offset,cmap='cmr.fusion_r',clim=(-200,200),transform=transform+ax.transData)
    # ax.set_ylabel(target)
    ax.text(-0.15+xshift[target],0.5,target,va='center',rotation=90,transform=ax.transAxes)

    # adding the transform if it's the right object
    # if target == 'CosmicEye': transform_kwarg = dict({'transform':transform+ax.transData})
        

    if i == 0: 
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="10%", pad=0.05)
        
        cbar = plt.colorbar(im,location='top',cax=cax)
        ticks = cbar.get_ticks() #np.array([-150,0,150])
        cbar.set_ticks(ticks) # just making sure they stay what they are
        
        labels = [str(round(t)) for t in ticks]
        for j,lab in enumerate(labels): 
            if ticks[j] < 0: labels[j] = '$-$' + str(abs(round(ticks[j])))
            # if j == 0: labels[j] = r'$\leq$' + labels[j]
            # if j == len(labels)-1: labels[j] = r'$\geq$' + lab
            
        cbar.set_ticklabels(labels,fontsize=13)

    if i == 6:
        ax.set_title('velocity offset [km/s]',y=-0.3,fontsize=14)
        
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="7%", pad=0.05)
        
        cbar = plt.colorbar(im,location='bottom',cax=cax)
        ticks = cbar.get_ticks() #np.array([-150,0,150])
        cbar.set_ticks(ticks) # just making sure they stay what they are
        
        labels = [str(round(t)) for t in ticks]
        for j,lab in enumerate(labels): 
            if ticks[j] < 0: labels[j] = '$-$' + str(abs(round(ticks[j])))
            # if j == 0: labels[j] = r'$\leq$' + labels[j]
            # if j == len(labels)-1: labels[j] = r'$\geq$' + lab
            
        cbar.set_ticklabels(labels,fontsize=13)
    
    # adding contours of Ha
    x0,y0 = np.arange(0,zs[0].shape[1]),np.arange(0,zs[0].shape[0])
    g = np.meshgrid(x0,y0)
    colors = ['k']
    contours = ax.contour(g[0],g[1],fluxes['ha'][0],
                          [np.nanmedian(fluxes['ha'][0])+np.nanstd(fluxes['ha'][0])*3],
                          origin='upper',alpha=0.8,
                          linewidths=1.8,colors=colors,
                          transform=transform+ax.transData)

    if target in lims.keys():
        ax.set_xlim(lims[target][0])
        ax.set_ylim(lims[target][1])


    
    
    # setting up the pride colormap to clip off both dark ends
    clipped_pride = co.ListedColormap([plt.get_cmap('cmr.pride')(j) 
                                           for j in np.linspace(0.1,0.9,1000)])
    

    # --- VELOCITY DISPERSION
    ax = plt.subplot(gs[1+i*4]); ax.axis('off')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    im = ax.imshow(velocity_dispersion,cmap=clipped_pride,clim=(75,175),transform=transform+ax.transData)

    # adding the transform if it's the right object
    # if target == 'CosmicEye': transform_kwarg = dict({'transform':transform+ax.transData})
        

    if i == 0: 
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="10%", pad=0.05)
        
        cbar = plt.colorbar(im,location='top',cax=cax)
        ticks = np.arange(100,175,50)
        cbar.set_ticks(ticks) # just making sure they stay what they are
        
        labels = [str(round(t)) for t in ticks]
        cbar.set_ticklabels(labels,fontsize=13)
        
    if i == 6: 
        ax.set_title(r'$\sigma_{H\alpha}$ [km/s]',y=-0.3,fontsize=14)
        
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="7%", pad=0.05)
        
        cbar = plt.colorbar(im,location='bottom',cax=cax)
        ticks = np.arange(100,175,50)
        cbar.set_ticks(ticks) # just making sure they stay what they are
        
        labels = [str(round(t)) for t in ticks]
        cbar.set_ticklabels(labels,fontsize=13)
    
    # adding contours of Ha
    x0,y0 = np.arange(0,sig[0].shape[1]),np.arange(0,sig[0].shape[0])
    g = np.meshgrid(x0,y0)
    colors = ['k']
    contours = ax.contour(g[0],g[1],fluxes['ha'][0],
                          [np.nanmedian(fluxes['ha'][0])+np.nanstd(fluxes['ha'][0])*3],
                          origin='upper',alpha=0.8,
                          linewidths=1.8,colors=colors,
                          transform=transform+ax.transData)
    
    if target in lims.keys():
        ax.set_xlim(lims[target][0])
        ax.set_ylim(lims[target][1])


    
    # O1 / Ha
    prelog_o1ha = fluxes['o1'][0].copy() / fluxes['ha'][0].copy()
    prelog_o1haerr = prelog_o1ha * np.sqrt((fluxes['o1'][1].copy()/fluxes['o1'][0].copy())**2 + (fluxes['ha'][1].copy()/fluxes['ha'][0].copy())**2)
    o1ha = np.log10(prelog_o1ha.copy())
    o1haerr = (prelog_o1haerr/prelog_o1ha) * 0.434
    
    # S2 / Ha
    prelog_s2ha = fluxes['s2'][0].copy() / fluxes['ha'][0].copy()
    prelog_s2haerr = prelog_s2ha * np.sqrt((fluxes['s2'][1].copy()/fluxes['s2'][0].copy())**2 + (fluxes['ha'][1].copy()/fluxes['ha'][0].copy())**2)
    s2ha = np.log10(prelog_s2ha.copy())
    s2haerr = (prelog_s2haerr/prelog_s2ha) * 0.434
    
    
    # ratio = 's2ha'
    # ratio = 'o1ha'
    rationames = [r'log$_{10}$ [SII] / H$\alpha$',r'log$_{10}$ [OI] / H$\alpha$']
    for v,ratio in enumerate(['s2ha','o1ha']):

        # --- FLUX RATIO
        ax = plt.subplot(gs[2+v+i*4]); ax.axis('off')
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        im = ax.imshow(eval(ratio),clim=(-3,0.7),transform=transform+ax.transData)
        if i == 6: 
            ax.set_title(rationames[v],y=-0.3,fontsize=14)
            
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("bottom", size="7%", pad=0.05)
            cbar = plt.colorbar(im,location='bottom',cax=cax)
            ticks = cbar.get_ticks()
            labels = [str(round(t)) for t in ticks]
            for j,l in enumerate(labels):
                if ticks[j] < 0: labels[j] = '$-$' + str(abs(round(ticks[j])))
            cbar.set_ticklabels(labels,fontsize=13)
    
        if i == 0: 
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("top", size="10%", pad=0.05)
            cbar = plt.colorbar(im,location='top',cax=cax)
            ticks = cbar.get_ticks()
            labels = [str(round(t)) for t in ticks]
            for j,l in enumerate(labels):
                if ticks[j] < 0: labels[j] = '$-$' + str(abs(round(ticks[j])))
            cbar.set_ticklabels(labels,fontsize=13)

        # adding the transform if it's the right object
        # if target == 'CosmicEye': transform_kwarg = dict({'transform':transform+ax.transData})
        
        
        # adding contours of Ha
        x0,y0 = np.arange(0,sig[0].shape[1]),np.arange(0,sig[0].shape[0])
        g = np.meshgrid(x0,y0)
        colors = ['k']
        contours = ax.contour(g[0],g[1],fluxes['ha'][0],
                              [np.nanmedian(fluxes['ha'][0])+np.nanstd(fluxes['ha'][0])*3],
                              origin='upper',alpha=0.8,
                              linewidths=1.8,colors=colors,
                              transform=transform+ax.transData)
    
        # adding a footprint
        footprint = SNRCUT(fluxes['ha'][0].copy(),fluxes['ha'][1].copy())/fluxes['ha'][1].copy()
        footprint = NOZERO(np.isfinite(footprint).astype(float))
        # ax.imshow(footprint,zorder=0,clim=(0,1),cmap='Grays',alpha=0.4,
        #           transform=transform+ax.transData)


        if target in lims.keys():
            ax.set_xlim(lims[target][0])
            ax.set_ylim(lims[target][1])


    # SNR cut for vdisp
    velocity_dispersion = SNRCUT(velocity_dispersion.copy(),velocity_dispersionerr,3)
    

    # writing to prep for other figure
    veloffdisp[target] = {'veloff':velocity_offset.copy(),
                          'veldisp':velocity_dispersion.copy()}




plt.savefig('plots-data/shocks-family-velocity.pdf')
plt.show()
plt.close('all')






# # getting ha+n2 broad velocity fits
# han2 = fits.getdata(path+f'{target}-han2double-velocity-params-stdpix.fits')
# han2err = fits.getdata(path+f'{target}-han2double-velocity-params-stdpix.fits',ext=2)
# han2norm = fits.getdata(path+f'{target}-han2double-velocity-params-stdpix.fits',ext=3)

# zs = [han2[0].copy(), han2err[0].copy()] # narrow comp
# zs_broad = [han2[1].copy(), han2err[1].copy()]
# sig = [han2[3].copy(), han2err[3].copy()]

# ha_narrow = [han2[2].copy(),han2err[2].copy()]
# ha_broad = [han2[2].copy()*han2[4].copy(),
#             (han2[2].copy()*han2[4].copy()) *\
#             np.sqrt( (han2err[2].copy()/han2[2].copy())**2 +\
#                      (han2err[4].copy()/han2[4].copy())**2 )]

















