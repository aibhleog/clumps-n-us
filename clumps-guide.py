'''

Making the guide figure for Gourav


'''

from fitting_ifu_spectra import *
from matplotlib.patches import Ellipse,Rectangle
import matplotlib.patheffects as PathEffects

target = 'SGAS2111'
grating = 'prism'


def add_ellipse(x,y,w,h,a,alph=0.3,c='r',rad=0):
    ellipseN = Ellipse((x,y),
                      width = w,
                      height = h,
                      angle = a,
                      alpha = alph,
                      facecolor = c)
    return ellipseN



# galaxy info
galaxy, z, fullpath, grating = get_galaxy_info(target,grating)
if grating == 'prism': extra = '-prism'
elif extra == None: extra = ''

# continuum-only map
cont_filename = f'plots-data/{target}-continuum-only-{grating}-map.fits'
cont_map_org = fits.getdata(cont_filename)




# CLUMPS DEFINED BY GOURAV & BRIAN
clumps = pd.read_csv('plots-data/clumps/region-files/'+\
                     # 'sgas2111_clumps_nircam_rgb-shiftedNIRSpec-pix.txt',sep=',',
                     f'{target.lower()}_clumps_nircam_rgb-shiftedNIRSpec-shiftingindividual{extra}-pix.txt',
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
                     f'and x < {cont_map_org.shape[1]}'+\
                     f'and y < {cont_map_org.shape[0]}').copy()








# plotting up just the 2D with regions
# --------------------------------------
cmap = plt.get_cmap('cmr.guppy')
colors = [cmap(j) for j in np.linspace(0.05,0.95,len(clumps))] 

plt.figure(figsize=(8,6))
ax = plt.gca()

clims = (-1e2,3e2)

ax.imshow(cont_map_org,origin='lower',clim=clims)
ax.axis('off')

ax.set_xlim(ax.get_xlim())
ax.set_ylim(ax.get_ylim())

# adding clumps to mask
ellipses = []
for i,j in enumerate(clumps.index.values):
    x,y,w,h,a,ID = clumps.loc[j].values # subtract 1 from x,y for DS9 --> python
    ellipse = add_ellipse(x,y,w=w+1,h=h+1,a=a,alph=1,c=colors[i])
    ellipses.append(ellipse)

    xx,yy = clumps.loc[j,'x'],clumps.loc[j,'y']
    if j == 45: xx,yy = clumps.loc[j,'x']+2,clumps.loc[j,'y']-3
    elif j == 32: xx,yy = clumps.loc[j,'x']-2,clumps.loc[j,'y']+3
    txt = ax.text(xx,yy,int(ID),color=colors[i])
    txt.set_path_effects([PathEffects.withStroke(linewidth=1.5, foreground='w')])

    
# ax.legend(handletextpad=0.35,fontsize=14,labelspacing=0.3)

# # legend order
# handles,labels = ax.get_legend_handles_labels()
# ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
#           fontsize=11,handlelength=0.9,handletextpad=0.35,labelspacing=0.3,
#           loc=2,bbox_to_anchor=(0,0.87))
# ax.legend
ax.set_title(target)

for ellipse in ellipses:
    ax.add_patch(ellipse) # make sure patch is added to plot last


plt.tight_layout()
plt.savefig(f'plots-data/clumps/{target}-clumps-guide.pdf')
plt.show()
plt.close('all')



sys.exit(0)


# setting up code for the few clumps highlighted in the survey paper

# I think originally identified in NIRSpec, there were 8 clump
# regions that Gourav had given me to make clump spectra
# now remaking some figures and will be giving him updated spectra 
clumps_wanted = np.array([34,35,32,26,27,42,43,45])
clumps_wanted -= 1 # from IDs-->index


# ----- code in terminal after running ha-balmerbreak.py
clims = (1e-20,2e-18)
cmap = plt.get_cmap('cmr.guppy')
colors = [cmap(j) for j in np.linspace(0.05,0.95,len(clumps_wanted))] 

plt.figure(figsize=(8,6))
ax = plt.subplot(111,projection=wcsNorthUp); ax.axis('off')

im = ax.imshow(cont,norm='log',cmap='inferno',clim=clims)
ax.set_title('Stellar Continuum',fontsize=14,y=1.03)
ax.set_xlim(xlims)
ax.set_ylim(ylims)



# adding northup
N = SkyCoord(ra='21h11m18.5s',dec='-01d14m30.5s')
E = SkyCoord(ra='21h11m18.484s',dec='-01d14m30.25s')
N_sky = RectangleSkyRegion(center=N,
                width=0.5*u.arcsec,
                height=0.02*u.arcsec)
E_sky = RectangleSkyRegion(center=E,
                width=0.5*u.arcsec,
                height=0.02*u.arcsec,angle=90*u.deg)
pixel_N = N_sky.to_pixel(wcsNorthUp)
pixel_E = E_sky.to_pixel(wcsNorthUp)
pixel_N.plot(ax=ax,color='k')
pixel_E.plot(ax=ax,color='k')
ax.text(0.72,0.115,'E',ha='right',transform=ax.transAxes,
        fontsize=12,fontweight='semibold')
ax.text(0.82,0.24,'N',ha='right',transform=ax.transAxes,
        fontsize=12,fontweight='semibold')

# adding scale bar
bar = SkyCoord(ra='21h11m18.78s',dec='-01d14m26.75s')
bar_sky = RectangleSkyRegion(center=bar,
                width=0.5*u.arcsec,
                height=0.02*u.arcsec)
pixel_bar = bar_sky.to_pixel(wcsNorthUp)
pixel_bar.plot(ax=ax,color='k')
ax.text(0.17,0.84,'0.5"',ha='center',transform=ax.transAxes,
        fontsize=12,fontweight='roman')



# ----- from spec-for-pedram ----------------
# reading in headerinfo
__,header = fits.getdata(galaxy['grating']['prism']['clipped'],header=True)
wcs = WCS(header).dropaxis(-1) # IFU aligned

# adding regions on top
from astropy.coordinates import Angle, SkyCoord
from regions import CircleSkyRegion,RectangleSkyRegion

# CLUMPS DEFINED BY GOURAV & BRIAN
tclumps = pd.read_csv('plots-data/clumps/region-files/'+\
                     f'{target.lower()}_clumps_nircam_rgb-shiftedNIRSpec-shiftingindividual{extra}-pix.txt',
                     sep=',',names=['x','y','w'])
tclumps[['x','y']] -= 1 # subtract 1 from x,y for DS9 --> python
tclumps = tclumps.loc[clumps_wanted].copy()
tclumps['r'] = tclumps.w * 0.03 # arcsec/pixel

for i,j in enumerate(clumps_wanted):
    x,y = tclumps.loc[j,['x','y']].values.copy()
    coord = wcs.pixel_to_world_values([[x,y]])[0]
    center = SkyCoord(ra=coord[0],dec=coord[1],unit='deg')
    radius = Angle(tclumps.loc[j,'r'],'arcsec')
    sky_region = CircleSkyRegion(center, radius)
    pixel_region = sky_region.to_pixel(wcsNorthUp)
    pixel_region.center.y += 1.5 # shifting all north a bit
    pixel_region.center.x += 0.5 # shifting all east a bit
    pixel_region.plot(ax=ax,color='w',lw=1.75) # outline
    # pixel_region.plot(ax=ax,color='w',lw=2.75) # outline
    # pixel_region.plot(ax=ax,color=colors[i],lw=2) # color


    xshift,yshift = -4.5,0
    if j == 33 or j == 26: yshift = 3
    if j == 42 or j == 25: xshift,yshift = -2,5
    txt = ax.text(pixel_region.center.x+xshift,
            pixel_region.center.y+yshift,
            # i,ha='right',
            clumps_wanted[i]+1,ha='right',
            fontsize=16,color='k',#colors[i],
            fontweight='medium')
    txt.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='w')])
    # txt.set_path_effects([PathEffects.withStroke(linewidth=0.8, foreground='w')])
# --------------------------------------------


plt.savefig('plots-data/SGAS2111-survey-clumps.pdf')
plt.savefig('plots-data/SGAS2111-survey-clumps.png',transparent=True)
plt.show()
plt.close('all')
