'''

fancy new clump spectra 

STEP1:  making the masks including partial coverage
        this takes a WHILE to run.


'''

from fitting_ifu_spectra import *
from matplotlib.patches import Ellipse,Rectangle
# importing the jwst_templates code
sys.path.append('../')
from jwst_templates import spec as jspec


# galaxy
target = 'SGAS1110'
saveit = True 
grating = 'g235m'


# galaxy info
galaxy, z, path, grating = get_galaxy_info(target,grating)
if grating == 'prism': extra = '-prism'
else: extra = ''


# just wrote this into a function to make life easier
def get_points(mypatch,size,radius=0.5):
    '''
    INPUTS:
    >> mypatch ---- the matplotlib patch shape
    >> size ------- integer, the size in pixels of one 
                    side of the "galaxy" image
                    
    OUTPUTS:
    >> points ----- the list of valid x,y coordinates
                    that overlap with the patch
    '''
    # create a list of possible coordinates
    x,y = np.arange(0,size),np.arange(0,size)

    g = np.meshgrid(x,y)
    coords = list(zip(*(c.flat for c in g)))

    # create the list of valid coordinates (from patch)
    points = np.vstack([p for p in coords if mypatch.contains_point(p, radius=radius)])
    return np.array(points)



def region_overlap(points1,points2):
    '''
    assumes you want to know how many points2 are in points1, with the final array being a binary sort of points2
    https://stackoverflow.com/questions/55434338/check-if-array-is-part-of-a-bigger-array
    '''
    isit = np.zeros(len(points2)).astype(int)
    for i,coord in enumerate(points2):
        if (points1 == coord).all(axis=1).any():
            isit[i] = 1
    return isit
    
    

def rebin(a,newshape):
    '''
    Rebin an array to a new shape.
    '''
    assert len(a.shape) == len(newshape)
    
    slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
    coordinates = np.mgrid[slices]
    indices = coordinates.astype('i')   #choose the biggest smaller integer index
    return a[tuple(indices)]


def mapping_pixel(x,y,scale=1):
    '''
    Takes a pixel coordinate from the original map and creates a 
    Rectangle patch that has nXn points (based upon the scale #).
    
    '''
    new_x,new_y = x*scale-0.5,y*scale-0.5
    return Rectangle((new_x,new_y), 1*scale, 1*scale, alpha=0.3,facecolor='r')


def get_pixel_coverage_mask(x,y):
    # using pixel coordinates to make a pixel
    rectangle = mapping_pixel(x,y,scale) # in scaled pixel coords, 1 pixel
    # counting points covered
    pix_points = get_points(rectangle,len(rebin_galaxy),radius=0.5/scale)
    overlap_index = region_overlap(check_rebin_points,pix_points)
    total_covered = len(pix_points[overlap_index==1])
    coverage = total_covered / scale**2
    return coverage



# reading in cube
# filename = galaxy['grating'][grating]['clipped']
filename = path+galaxy['grating'][grating]['filename'][:-5]+'-FSbkgd.fits'
data, header = fits.getdata(filename,header=True)

# defining things
sli = galaxy['grating'][grating]['slice']
data_slice = data[sli]



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
                     f'and x < {data.shape[2]}'+\
                     f'and y < {data.shape[1]}').copy()


def add_ellipse(data,x,y,w,h,a,alph=0.3,c='r',rad=0):
    ellipseN = Ellipse((x,y),
                      width = w,
                      height = h,
                      angle = a,
                      alpha = alph,
                      facecolor = c)

    pointsN = get_points(ellipseN,np.nanmax(data.shape),radius=rad)
    # cont_map[pointsN[:,1],pointsN[:,0]] = np.nan # cutting out pixels
    return ellipseN,pointsN



scale = 5 # makes a 5X5 grid for each pixel

# rebinning to more points for the refined part!
rebin_data_slice = rebin(data_slice.copy(), np.asarray(data_slice.shape)*scale)


# regions
ellipses = []

# all coverages
all_coverage = np.zeros_like(data_slice)

# running through clumps:
for jj in clumps.index.values:
    print(f'{clumps.index.values.tolist().index(jj)}/{len(clumps)}')
    x,y,w,h,a,ID = clumps.loc[jj].values 

    # the same aperture but scaled to the new grid
    # but with ~regular radius
    __,check_rebin_points = add_ellipse(rebin_data_slice,
                                (x+0.5)*scale-(1/scale),(y+0.5)*scale-(1/scale),
                                w=w*scale,h=h*scale,a=a,rad=1)
    
    # original aperture, larger radius
    ellipse,check_points = add_ellipse(data_slice,x,y,w=w,h=h,a=a,rad=1)
    ellipses.append(ellipse)

    # empty array
    coverage = np.zeros_like(data_slice)
    
    # running through all of the coordinates in the ORIGINAL grid
    for i,coord in enumerate(check_points):
        print(f'At coordinate {i}/{len(check_points)}',end=',')
        
        x,y = coord # in original pixel coords
        rectangle = mapping_pixel(x,y,scale) # in scaled pixel coords, 1 pixel
        
        # counting points covered
        pix_points = get_points(rectangle,np.nanmax(rebin_data_slice.shape),radius=0.5/scale)
        overlap_index = region_overlap(check_rebin_points,pix_points)
        total_covered = len(pix_points[overlap_index==1])
        coverage[y,x] = total_covered / scale**2 # fraction coverage
        all_coverage[y,x] = total_covered / scale**2 # fraction coverage

        print(f'coverage: {coverage[y,x]}')
    
    print(f'total coverage: {np.nansum(coverage.flatten())}')
        
    if saveit == True:
        j = jj
        hdu = fits.PrimaryHDU(coverage.copy())
        hdu.writeto(f'plots-data/clumps/{target}-clump{int(ID)}-complete-mask{extra}.fits',overwrite=True)
        print('saved mask')
    
    print(end='\n\n')
    

coverage[coverage==0] = np.nan # we don't care about the non-aperture pixels



plt.figure(figsize=(8,6))
plt.axis('off')

im = plt.imshow(all_coverage,origin='lower',cmap='Greens')
plt.colorbar(im)

for i,indx in enumerate(clumps.index.values):
    plt.text(clumps.loc[indx,'x'],clumps.loc[indx,'y'],indx)

    # add patch
    plt.gca().add_patch(ellipses[i])

plt.tight_layout()
# plt.savefig(f'plots-data/clumps/clumps-pixel-coverage{extra}.pdf')
plt.show()
plt.close('all')















