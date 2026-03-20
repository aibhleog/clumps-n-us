'''
This module take an input pandas dataframe spectrum and does things with it.
My brain works in pandas dataframes, sorry. :P

The functions within this module each have summaries within them of the
inputs required (and optional) and then what the function outputs.


Last updated:  June 2024


'''

__author__ = 'Taylor Hutchison'
__email__ = 'astro.hutchison@gmail.com'

import warnings
import numpy as np

from astropy import units as u
from astropy.nddata import StdDevUncertainty
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.modeling import models
from astropy.stats import sigma_clip
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cmasher
import pandas as pd
import astropy.io.fits as fits
from astropy.wcs import WCS
import sys,json


def veloff(z1,z2):
    delv = 2.998e5 * (z2 - z1) / (1 + z1)
    return delv # in km/s


def get_galaxy_info(target,grat):
    # reading in file that has all of the galaxy values
    with open('plots-data/galaxies.txt') as f:
        data = f.read()

    # reconstructing dictionary
    galaxies = json.loads(data)
    path = galaxies['path']

    galaxy = galaxies[target] # getting galaxy info
    z = galaxy['z']
    grating = grat
        
    return galaxy,z,path,grating

        

def get_coords(cube_shape):
    '''
    INPUTS:
    >> cube_shape --- a tuple with 3 numbers, representing a cube with 3 dim.
                      can be calculated doing cube.shape
                      
    OUTPUTS:
    >> coordinates -- a list of coordinates from the cube dimensions   
    
    '''
    # making list of spaxel coordinates based on data shape
    x0,y0 = np.arange(0,cube_shape[2]),np.arange(0,cube_shape[1])
    g = np.meshgrid(x0,y0)
    coords = list(zip(*(c.flat for c in g)))
    return coords
        
        
        
def get_mask(galaxy,array_2d=False,layers=False,subpix=False,extra=None):
    '''
    INPUTS:
    >> galaxy -------- the name of the galaxy mask I want
    
    OUTPUTS:
    >> coordinates --- a list of coordinates from the mask,
                       specifying where the galaxy light is
    '''
    
    if extra == None: extra = ''
    # for the sub-pixel drizzling
    if subpix == True:
        extra += '-0.05pix'
    
    # just the full galaxy mask
    if layers == False:
        try:
            galaxy_mask = fits.getdata(f'plots-data/{galaxy}-mask{extra}.fits')

            # makings a list of coordinates
            coordinates = list(zip(*np.where(galaxy_mask == 1)))

            if array_2d == False: return coordinates
            else: return galaxy_mask

        except:
            print("\nWrong file and/or file doesn't exist yet.",end='\n\n')
            sys.exit(0) # kills script

            
    # want the mask slices instead
    # note that the first slice is the full galaxy mask
    else:
        try:
            filename = f'plots-data/{galaxy}-mask-layers{extra}.fits'
            filename2 = f'plots-data/{galaxy}-mask-layers{extra}.txt'
            
            galaxy_mask = fits.getdata(filename)
            mask_layers_info = np.loadtxt(filename2,delimiter='\t')
            
            mask_layers = []
                
            with fits.open(filename) as hdul:
                for i in range(len(hdul)):
                    # map layer
                    galaxy_mask = hdul[i].data

                    # makings a list of coordinates
                    coordinates = list(zip(*np.where(galaxy_mask == 1)))

                    if array_2d == False: mask_layers.append(coordinates)
                    else: mask_layers.append(galaxy_mask)
                
            return mask_layers, mask_layers_info
                

        except:
            print("\nWrong layers file and/or file doesn't exist yet.",end='\n\n')
            sys.exit(0) # kills script


def get_colors(cmap,size):
    cmap = plt.get_cmap(cmap)
    colors = [cmap(j) for j in np.linspace(0,0.8,size)]
    return colors
        
        
def gaussian(xaxis, mean, A1, sig):
    '''
    Model to used in fitting single line.
    
    '''
    g1 = A1*np.exp(-np.power(xaxis - mean, 2.) /( 2 * np.power(sig, 2.)))
    return g1
        
        
        
def gaussian_doublet(xaxis, mean, A1, A2, sig, sep):
    '''
    Model to be used in fitting blended doublet.
    
    '''
    g1 = A1*np.exp(-np.power(xaxis - mean, 2.) /( 2 * np.power(sig, 2.)))
    
    mean2 = mean + sep
    g2 = A2*np.exp(-np.power(xaxis - mean2, 2.) /( 2 * np.power(sig, 2.)))
    
    return g1+g2
        
        
        
        
        
def spec_wave_range(spec,wave_range,index=False):
    '''
    INPUTS:
    >> spec  --------  a pandas dataframe set up with columns
                       "wave", "flam", "ferr"
    >> wave_range  --  a list of 2 numbers, min & max wavelengths
                       describing the wavelength range of the line
    >> index --------  a flag that if True only returns index vals
                       
    OUTPUTS:
    >> results ------  a "zoomed in" dataframe
    '''
    
    # zooms in to specified wavelength range
    wave_query = f'{wave_range[0]}<wave<{wave_range[1]}'
    spec = spec.query(wave_query).copy()
    
    # if I only want the index values, nothing else
    if index == True:
        spec = spec.index.values
    
    return spec


def get_spec(x,y,dd,h,ee=None,verbose=False,cgs=False):
    '''
    Makes spectrum from one spaxel or many spaxels
    
    INPUTS:
    >> x,y  -----------  either 1 value each or 1 array of values each
    >> dd,h  ----------  data cube and header info
    >> ee (opt) -------  optional error cube (default None)
    >> verbose (opt) --  if True, will return pix coords used 
    >> cgs (opt) ------  if True, will return spectrum in cgs units
                       
    OUTPUTS:
    >> spec ---------  spectrum dataframe in [microns], [My/sr] (or cgs)
    '''
    d = dd.copy()
    
    # checking if 1 pixel or many
    if type(x) != int and type(x) != np.int64:
        mask = np.zeros_like(d[0])

        if len(x) == len(y):
            mask[y,x] = 1
        else:
            s = []
            for xx in x: 
                for yy in y: 
                    s.append([xx,yy])
            s = np.asarray(s)
            mask[s[:,1],s[:,0]] = 1
            x,y = s[:,0],s[:,1]
            
        longmask = np.broadcast_to(mask, d.shape)
        d[longmask<1] = np.nan
        dat = np.nansum(d,axis=(1,2))

        if type(ee) != type(None): e = ee.copy(); err = np.nansum(e,axis=(1,2))
        else: err = np.zeros_like(dat)
            
    else:
        dat = [float(f) for f in d[:,int(y),int(x)].copy()]
        
        if type(ee) != type(None): 
            e = ee.copy()
            err = [float(f) for f in e[:,int(y),int(x)].copy()]
        else: err = np.zeros_like(dat)

    # rest of prep
    wave = np.arange(h['CRVAL3'], 
                     h['CRVAL3']+(h['CDELT3']*len(d)), 
                     h['CDELT3'])
    
    if len(wave)>len(dat): wave = wave[:-1].copy()
    
    spec = pd.DataFrame({'wave':wave,'flam':dat,'flamerr':err})

    if cgs == True:
        spec = convert_MJy_cgs(spec.copy())

    
    if verbose == True: return spec.copy(), [x,y]
    else: return spec.copy()


def convert_MJy_sr_to_MJy(spec):
    '''
    The spectra from the reduced data are in MJy/sr.
    Converting to MJy so they can be converted to cgs units.
    
    INPUTS:
    >> spec  --------  a pandas dataframe set up with columns
                       "wave", "flam", "ferr"
                       
    OUTPUTS:
    >> spec ---------  the same pandas dataframe but in MJy    
    '''
    # taking nominal pixel area from FITS header for data
    pix_area = 2.35040007004737E-13 # in steradians
    
    # converting spectrum flam 
    spec['flam'] *= pix_area # MJy/sr --> MJy

    # converting spectrum error
    spec['flamerr'] *= pix_area # MJy/sr --> MJy
    
    return spec.copy()
    


def convert_MJy_cgs(spec):
    '''
    INPUTS:
    >> spec  --------  a pandas dataframe set up with columns
                       "wave", "flam", "ferr"; assumes wave is 
                       in microns and flam, ferr are in MJy
    OUTPUTS:
    >> spec ---------  the same pandas dataframe but in cgs
    '''
    # converting from MJy/sr to MJy
    spec = convert_MJy_sr_to_MJy(spec.copy()) 
    
    # converting spectrum flam to cgs units
    spec['flam'] *= 1e6 # MJy --> Jy
    spec['flam'] *= 1e-23 # Jy --> erg/s/cm/Hz (fnu)
    spec['flam'] *= 2.998e18 / (spec.wave.values*1e4)**2 # fnu --> flam
    
    # converting spectrum error to cgs units
    spec['flamerr'] *= 1e6 # MJy --> Jy
    spec['flamerr'] *= 1e-23 # Jy --> erg/s/cm/Hz (fnu)
    spec['flamerr'] *= 2.998e18 / (spec.wave.values*1e4)**2 # fnu --> flam
    
    return spec.copy()


def moving_average(a, n=3):
    '''
    INPUTS:
    >> a  ----------  an array of values
    >> n (opt) -----  window size for moving average (for the
                      initial continuum fit)
                       
    OUTPUTS:
    >> avg ---------  a moving average array the same length as "a"
    
    
    The original code for the moving average can be found here:
https://stackoverflow.com/questions/14313510/how-to-calculate-moving-average-using-numpy/54628145
    
    '''
    ret = np.nancumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    
    avg = ret[n-1:]/n
    avg = np.concatenate((np.ones(6)*avg[0], # to make same len as spec
                          avg))
    
    return avg




# CODE FOR JULISSA, PULLED FROM THE JWST-TEMPLATES/jwst_templates GITHUB REPO
# ----------------------------------------


def integrate1D_mask(cubefile, maskfile, operation=np.nansum):
    '''
    INPUTS:
    >> cubefile  ----  the file name+location of the IFU cube FITS file
    >> maskfile  ----  the file name+location of the mask FITS file
    >> operation  ---  the numpy operation to perform (default nansum)
    
    OUTPUTS:
    >> df_cgs  ------  the pandas dataframe spectrum with columns
                       wave, flam, and flamerr (converted into cgs units)
    
    ---------
    Simplified from JWST-TEMPLATES/jwst_templates respository
    https://github.com/JWST-Templates/jwst_templates/blob/main/spec.py
    '''
    # reading in IFU cube and header info
    hdu = fits.open(cubefile)
    header = hdu[1].header
    dat = hdu[1].data  # in MJy/sr
    err = hdu[2].data  # error, same units as dat
    
    # using header info
    pix_area = header['PIXAR_SR'] # pixel area, in sr
    wl = np.arange(header['CRVAL3'],  # wavelength in micron
            header['CRVAL3']+(header['CDELT3']*len(dat)), 
            header['CDELT3'])
    if len(wl) > len(dat):
        wl = wl[:-1] # deal with weird arange glitch adding an extra wl value


    # reading in the mask of the galaxy (0=not galaxy, 1=galaxy)
    mask = fits.getdata(maskfile) 
    if mask.shape != dat[0].shape:
        mask = mask.T
    
    
    # broadcasting mask to IFU cube
    longmask = np.broadcast_to(mask, dat.shape)
    
    # using chosen operation on flux density and associated error
    flux = operation(dat * longmask, axis=(1,2))
    error = np.sqrt(operation(err**2 * longmask, axis=(1,2)))
    
    # making dataframe of final spectrum
    df = pd.DataFrame({'wave':wl, 'flam':flux, 'flamerr':error})
    df_cgs = convert_jwst_to_cgs(df, pix_area)
        
    return df_cgs


def convert_jwst_to_cgs(spec,pix_area):
    '''
    INPUTS:
    >> spec  --------  a pandas dataframe set up with columns
                       "wave", "flam", "ferr"; assumes wave is 
                       in microns and flam, ferr are in MJy
    OUTPUTS:
    >> spec ---------  the same pandas dataframe but in cgs
    '''    
    # converting flam and error to MJy
    spec['flam'] *= pix_area # MJy/sr --> MJy
    spec['flamerr'] *= pix_area # MJy/sr --> MJy
    
    # converting spectrum flam to cgs units
    spec['flam'] *= 1e6 # MJy --> Jy
    spec['flam'] *= 1e-23 # Jy --> erg/s/cm/Hz (fnu)
    spec['flam'] *= 2.998e18 / (spec.wave.values*1e4)**2 # fnu --> flam
    
    # converting spectrum error to cgs units
    spec['flamerr'] *= 1e6 # MJy --> Jy
    spec['flamerr'] *= 1e-23 # Jy --> erg/s/cm/Hz (fnu)
    spec['flamerr'] *= 2.998e18 / (spec.wave.values*1e4)**2 # fnu --> flam
    
    return spec.copy()

