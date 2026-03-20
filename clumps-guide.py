'''

Making the guide figure for Gourav


'''

from fitting_ifu_spectra import *
from matplotlib.patches import Ellipse,Rectangle
import matplotlib.patheffects as PathEffects

target = 'SGAS1110'
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
















