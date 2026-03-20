'''

STEP2:  integrating the spectra using the pixel masks from Step 1


'''

from fitting_ifu_spectra import *
import warnings
warnings.filterwarnings("ignore")
# importing the jwst_templates code
sys.path.append('../')
from jwst_templates import continuum


target = 'SGAS1110'
saveit = True 
grating = 'prism'


# galaxy info
galaxy, z, path, grating = get_galaxy_info(target,grating)
if grating == 'prism': extra = '-prism'
else: extra = ''


# making helper function for integration that's a tad more flexible
def integrate_spec(s3d,mask):
    # broadcasting mask to IFU cube
    longmask = np.broadcast_to(mask, s3d.shape)
    # integrating to 1D
    array = np.nansum(s3d * longmask, axis=(1,2))
    return array.copy()

def convert_MJy_sr_to_fnu(dataframe,pixar_sr):
    dataframe[['fnu','fnuerr']] *= pixar_sr * 1e6 * 1e-23 # MJy/sr -> erg/s/cm2/Hz
    return dataframe.copy()


def get_no_line(wl,z_sys):
    '''
    using the flag_near_lines function from the jwst_templates.continuum module
    to mask out emission lines so no_line can flag only "continuum" regions
    '''
    spec = pd.DataFrame({'wave':wl,'flam':np.zeros(len(wl)),'flamerr':np.zeros(len(wl))})
    spec['no_line'] = spec.index.values.copy() # easier to mess with this way
    LL = continuum.load_default_linelist(v2mask=1200)  # for R~1000
    
    continuum.flag_near_lines(spec,z_sys,LL,colv2mask='v2mask')
    spec['no_line'] = False
    spec.loc[spec.linemask == True,'no_line'] = True
    return spec.no_line.values.copy()






# reading in cube g235m
# filename = galaxy['grating']['g235m']['clipped']
filename = path+galaxy['grating']['g235m']['filename'][:-5]+'-FSbkgd.fits'
data, header = fits.getdata(filename,header=True)
error = fits.getdata(filename,ext=2)
wav1 = np.arange(header['CRVAL3'],  # wavelength in micron
            header['CRVAL3']+(header['CDELT3']*(len(data)-1)), # need for g235m
            header['CDELT3'])
pixar_sr1 = header['PIXAR_SR'] # pixel area, in sr

# reading in cube prism
# filename2 = galaxy['grating']['prism']['clipped']
filename2 = path+galaxy['grating']['prism']['filename'][:-5]+'-FSbkgd.fits'
data2, header2 = fits.getdata(filename2,header=True)
error2 = fits.getdata(filename2,ext=2)
wav2 = np.arange(header2['CRVAL3'],  # wavelength in micron
            header2['CRVAL3']+(header2['CDELT3']*len(data2)), 
            header2['CDELT3'])
pixar_sr2 = header2['PIXAR_SR'] # pixel area, in sr





# CLUMPS DEFINED BY GOURAV & BRIAN
clumps = pd.read_csv('plots-data/clumps/region-files/'+\
                     # 'sgas2111_clumps_nircam_rgb-shiftedNIRSpec-pix.txt',sep=',',
                     f'{target.lower()}_clumps_nircam_rgb-shiftedNIRSpec-shiftingindividual-pix.txt',
                     sep=',',names=['x','y','w'])
clumps[['x','y']] -= 1 # subtract 1 from x,y for DS9 --> python
clumps['h'] = clumps.w.values.copy()
clumps['a'] = 0 # degree
clumps[['w','h']] *= 2 # radius --> diameter b/c Ellipse not Circle
clumps['ID'] = [int(f+1) for f in clumps.index.values]

if target == 'SGAS2111':
    # dropping other source clump
    clumps.drop(index=37,inplace=True)

# dropping any with negative X,Y or over max X,Y values (not in IFS FOV)
clumps = clumps.query('x > 0 and y > 0'+\
                     f'and x < {data.shape[2]}'+\
                     f'and y < {data.shape[1]}').copy()




cmap = plt.get_cmap('cmr.guppy')
colors = [cmap(j) for j in np.linspace(0.05,0.95,len(clumps))]


# RUNNING THROUGH MASKS MAKING SPECTRA

spect,ylims = [],[]
plt.figure(figsize=(12,35))
gs = gridspec.GridSpec(len(clumps),1,height_ratios=np.ones(len(clumps)),hspace=0)

for i,j in enumerate(clumps.index.values.copy()):
    ID = int(clumps.loc[j,'ID'])

    # clump mask g235m
    maskfile = f'plots-data/clumps/{target}-clump{ID}-complete-mask.fits'
    mask = fits.getdata(maskfile)
    num_pix = np.nansum(mask)
    print(f'for clump {j}, g235m, total pixels {round(num_pix,3)}',end='; ')

    # clump mask prism
    maskfile2 = f'plots-data/clumps/{target}-clump{ID}-complete-mask-prism.fits'
    mask2 = fits.getdata(maskfile2)
    num_pix2 = np.nansum(mask2)
    print(f'prism, total pixels {round(num_pix2,3)}')

    

    # LETS FUCKING GO
    ax = plt.subplot(gs[i])


    # first grating
    flux = integrate_spec(data,mask)
    ferr = np.sqrt(integrate_spec(error**2,mask))
    spec1 = pd.DataFrame({'wave':wav1,'fnu':flux,'fnuerr':ferr})
    spec1 = convert_MJy_sr_to_fnu(spec1.copy(),pixar_sr1)
    
    # second grating
    flux = integrate_spec(data2,mask2)
    ferr = np.sqrt(integrate_spec(error2**2,mask2))
    spec2 = pd.DataFrame({'wave':wav2,'fnu':flux,'fnuerr':ferr})
    spec2 = convert_MJy_sr_to_fnu(spec2.copy(),pixar_sr2)


    final_spec = spec1.copy()
    final_spec2 = spec2.copy()

    
    ax.step(final_spec2.wave,final_spec2.fnu,where='mid',alpha=0.8,
            label=f'Clump {ID}',lw=1,color=colors[i])
    # spect.append(final_spec2.copy())

    if saveit == True:
        final_spec.to_csv(f'plots-data/clumps/{target}-clump{ID}-spectrum-g235m.txt',
                          sep='\t',index=False)
        final_spec2.to_csv(f'plots-data/clumps/{target}-clump{ID}-spectrum-prism.txt',
                          sep='\t',index=False)

    med_y = np.nanmedian(final_spec2.loc[:100,'fnu'])
    ax.set_ylim(med_y-med_y*0.8,)
    ylims.append(ax.get_ylim())
    
    ax.axvline(.3727*(1+z),zorder=0,color='k',alpha=0.4,ls=':',lw=1)
    ax.axvline(.4863*(1+z),zorder=0,color='k',alpha=0.4,ls=':',lw=1)
    ax.axvline(.4960*(1+z),zorder=0,color='k',alpha=0.4,ls=':',lw=1)
    ax.axvline(.5008*(1+z),zorder=0,color='k',alpha=0.4,ls=':',lw=1)
    ax.axvline(.6564*(1+z),zorder=0,color='k',alpha=0.4,ls=':',lw=1)

    
    leg = ax.legend(fontsize=11,loc=2,handlelength=0.5,handletextpad=1.1)
    leg.legend_handles[0].set_linewidth(7)
    # ax.set_title(method)
    ax.set_yscale('log')
    # ax.set_ylim(0,)
    if i == len(clumps)-1: ax.set_xlabel('observed wavelength [microns]')
    else: ax.set_xticklabels([])
    if i == int(len(clumps)/2): ax.set_ylabel('flux density [erg/s/cm$^2$/Hz]')


plt.tight_layout()
plt.savefig(f'plots-data/clumps/{target}-clumps-spectra.pdf')
plt.show()
plt.close('all')




# # saving new clumps file
# clumps = clumps[[ 'ID','median_mu','avg_mu','x', 'y', 'w', 'h', 'a']].copy()
# # clumps.to_csv('clumps/clumps_f115w_wcsmatch_pixelcoords-magmedian.txt',
# clumps.to_csv('clumps/clump_regions_f460m_v1-pix-magmedian.txt',
#               sep='\t',index=False)












sys.exit(0)


'''

Older code, keeping for future reference.  You should ignore :)

'''



if method == 'divavgmu-sublocalbcg':
    
    plt.figure(figsize=(12,35))
    gs = gridspec.GridSpec(len(clumps),1,height_ratios=np.ones(len(clumps)),hspace=0)
    
    for jj in clumps.index.values.copy():
        ax = plt.subplot(gs[jj])
        # renaming one clump to match photometry numbering
        j = jj
    
        # clump mask
        maskfile = f'plots-data/clumps/waz-clump{j}-complete-mask.fits'
        mask = fits.getdata(maskfile)
        num_pix = np.nansum(mask)
        print(f'for clump {j}',end=', ')
    
        # identifying which BCG region corresponds to the slice
        bcg_slices_filler = bcg_slices.copy()
        bcg_slices_filler[mask==0] = np.nan
        regions_overlapped = np.array(list(set(bcg_slices_filler[np.isfinite(bcg_slices_filler)].flatten())))
        print(f'bcg regions overlapped with: {regions_overlapped}')
        # reading in those region(s)
        bcg_slices_spec = []
        for r in regions_overlapped:
            filler = pd.read_csv(f'plots-data/measuring-bcg/region-masks/region{int(r)}-middle-spectrum.txt',sep='\t')
            bcg_slices_spec.append(filler)
    
        # LOCAL BCG TO SUBTRACT
        bcg_final = bcg_pix[['wave']].copy()
        if len(bcg_slices_spec) > 1:
            bcg_final['fnu_smooth'] = np.nanmean([
                bcg_slices_spec[0].fnu_smooth.values.copy(),
                bcg_slices_spec[1].fnu_smooth.values.copy()],axis=0)
            bcg_final['fnuerr_smooth'] = np.sqrt(
                bcg_slices_spec[0].fnuerr_smooth.values.copy()**2 +
                bcg_slices_spec[1].fnuerr_smooth.values.copy()**2)

            ax.step(bcg_final.wave,bcg_slices_spec[0].fnu_smooth,
                    where='mid',alpha=0.8,lw=1,label=f'local BCG region slice {int(regions_overlapped[0])}')
            ax.step(bcg_final.wave,bcg_slices_spec[1].fnu_smooth,
                    where='mid',alpha=0.8,lw=1,label=f'local BCG region slice {int(regions_overlapped[1])}')
            
        else:
            bcg_final['fnu_smooth'] = bcg_slices_spec[0].fnu_smooth.values.copy()
            bcg_final['fnuerr_smooth'] = bcg_slices_spec[0].fnuerr_smooth.values.copy()

            ax.step(bcg_final.wave,bcg_slices_spec[0].fnu_smooth,
                    where='mid',alpha=0.8,lw=1,label=f'local BCG region slice {int(regions_overlapped[0])}')
    
            
        ax.step(bcg_final.wave,bcg_final.fnu_smooth,where='mid',
                alpha=0.8,color='gray',label=f'local BCG for Clump {jj}')
        
        # ax.set_ylim(ylims[jj])
    
        leg = ax.legend(fontsize=11,loc=1,handlelength=0.5,
                        frameon=True,handletextpad=1.1)
        for l in leg.legend_handles: l.set_linewidth(7)
        ax.set_yscale('log')

        if jj == len(clumps)-1: ax.set_xlabel('observed wavelength [microns]')
        else: ax.set_xticklabels([])
        if jj == int(len(clumps)/2): ax.set_ylabel('flux density [erg/s/cm$^2$/Hz]')

        if saveit == True:
            bcg_final.to_csv(f'plots-data/clumps/waz-clump{j}-spectrum-{method}-localBCG.txt',
                          sep='\t',index=False)

    
    plt.tight_layout()
    plt.savefig(f'plots-data/clumps/waz-all-clumps-spectra-{method}-localBCG.pdf')
    plt.show()
    plt.close('all')

    
