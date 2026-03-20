'''

plotting up the 5 clumps pedram will include in his lensing paper


'''
import reproject
import warnings
import matplotlib.patheffects as PathEffects
warnings.filterwarnings("ignore")
from fitting_ifu_spectra import *
galaxy, z, datapath, grating = get_galaxy_info('SGAS1110',grat='g235m')

# reading in data
contmap,header = fits.getdata(galaxy['grating']['prism']['cont_map'],header=True)

path = 'plots-data/clumps/'
clumps = [15,17,23,24,25]
lens_clump_names = [26.2,25.2,22.2,23.2,24.2]

spec_g235m = pd.DataFrame({})

for i in range(len(clumps)):
    filler_g235m = pd.read_csv(path+f'SGAS1110-clump{clumps[i]}-spectrum-g235m.txt',sep=r'\t')

    if i == 0:
        spec_g235m['wave'] = filler_g235m.wave.values.copy()

    spec_g235m[f'clump{clumps[i]}'] = filler_g235m.fnu.values.copy()



# looking at spectra
cmap = plt.get_cmap('cmr.guppy')
colors = [cmap(j) for j in np.linspace(0.05,0.6,len(clumps))] 
colors = colors[::-1]


plt.figure(figsize=(13,4.5))
gs = gridspec.GridSpec(1,2,width_ratios=[3,1],wspace=0.02)


# THE 1D G235M SPECTRA
ax = plt.subplot(gs[0])

for i in range(len(clumps)):
    ax.step(spec_g235m.wave,spec_g235m[f'clump{clumps[i]}'],
            where='mid',alpha=0.65,
            zorder=5,color=colors[i])


names = [r'H$\beta$','','[OIII]','HeI',r'H$\alpha$','[SII]']
for i,l in enumerate([.4864,.4960,.5008,.5877,.6564,.6717]):
    ax.axvline(l*(1+z),color='k',
               ls=':',zorder=-10,lw=0.85)
    ax.text(l*(1+z)+0.0045,1e-29,names[i],rotation=90)

ax.set_ylabel(r'$f_\nu$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]')
ax.set_xlabel('observed wavelength [microns]')
ax.set_yscale('log')
ax.set_ylim(5e-32,)
ax.set_xlim(spec_g235m.wave.min(),2.38)

ax.text(1.85,2e-29,'z$_{spec}$ = '+f'{z}',
        fontfamily='serif',fontsize=17)




# THE 2D IMAGE WITH REGIONS OVERLAID

# reprojecting contmap into a north up frame
# using an older reduction for its WCS
wcs = WCS(header)
__,headerNorthUp = fits.getdata(datapath+'SGAS1110/L3/old_pmaps/Level3_SGAS1110_NOBG_OUTLIER_XY0p03_prism-clear_s3d.fits',header=True)
wcsNorthUp = WCS(headerNorthUp).dropaxis(-1)
contmapNorthUp, footprint = reproject.reproject_interp((contmap,header),wcsNorthUp)

ax = plt.subplot(gs[1],projection=wcsNorthUp)
im = ax.imshow(contmapNorthUp,norm='log',
               cmap='viridis',clim=(8,500))
ax.axis('off')

# adding regions on top
from astropy.coordinates import Angle, SkyCoord
from regions import CircleSkyRegion,RectangleSkyRegion

# CLUMPS DEFINED BY GOURAV & BRIAN
tclumps = pd.read_csv('plots-data/clumps/region-files/'+\
                     f'SGAS1110_clumps_nircam_rgb-shiftedNIRSpec-shiftingindividual-pix.txt',
                     sep=',',names=['x','y','w'])
tclumps[['x','y']] -= 1 # subtract 1 from x,y for DS9 --> python
tclumps = tclumps.loc[[14,16,22,23,24]].copy()
tclumps['r'] = tclumps.w * 0.03 # arcsec/pixel

for i,j in enumerate(tclumps.index.values):
    x,y = tclumps.loc[j,['x','y']].values.copy()
    coord = wcs.pixel_to_world_values([[x,y]])[0]
    center = SkyCoord(ra=coord[0],dec=coord[1],unit='deg')
    radius = Angle(tclumps.loc[j,'r'],'arcsec')
    sky_region = CircleSkyRegion(center, radius)
    pixel_region = sky_region.to_pixel(wcsNorthUp)
    pixel_region.plot(ax=ax,color='w',lw=2.75) # outline
    pixel_region.plot(ax=ax,color=colors[i],lw=2) # color

    # adding region ID
    xshift,yshift = -7,0
    if 1 < i < 4: yshift = -1.5 + (i-3)*2.5 # 24.2 and 23.2 are close together
    print(yshift)
    txt = ax.text(pixel_region.center.x+xshift,
            pixel_region.center.y+yshift,
            lens_clump_names[i],ha='right',
            fontsize=14,color=colors[i],
            fontweight='medium')
    txt.set_path_effects([PathEffects.withStroke(linewidth=0.8, foreground='w')])




# adding scale bar
bar = SkyCoord(ra=167.58397,dec=64.997177,unit='deg')
bar_sky = RectangleSkyRegion(center=bar,
                width=0.5*u.arcsec,
                height=0.02*u.arcsec)#,
                # angle=5 * u.deg)
pixel_bar = bar_sky.to_pixel(wcsNorthUp)
pixel_bar.plot(ax=ax,color='w')

ax.text(0.105,0.24,'0.5"',transform=ax.transAxes,
        color='w',fontsize=14,fontweight='medium')


ax.set_xlim(25.5,145.5)
ax.set_ylim(50.5,220.5)


plt.savefig('plots-data/pedram-SGAS1110-zspec-confirm.pdf')
plt.show()
plt.close('all')










