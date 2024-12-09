import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Ellipse, Circle
from astropy.modeling.functional_models import Ellipse2D
from scipy.optimize import curve_fit
import glob
import pdb
from astropy.visualization import simple_norm, ZScaleInterval
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization.stretch import LinearStretch, LogStretch
import warnings
warnings.filterwarnings('ignore')
plt.style.use(["science", "default"])

alma_images = sorted(glob.glob("/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/all_ALMA/*.fits"))
rms_dict = {'CDFS-18': 2.54e-05,
 'CDFS-19': 2.37e-05,
 'ES1_z_0.88': 1.79e-05,
 'ES1_z_0.99': 1.78e-05,
 'ES1_z_0.99b': 1.79e-05,
 'ES1_z_1.04': 1.69e-05,
 'ES1_z_1.38': 1.71e-05,
 'ES1_z_1.40': 1.65e-05,
 'ES1_z_1.60': 1.68e-05,
 'ES1_z_1.65': 1.72e-05,
 'ES1_z_1.70': 1.75e-05,
 'ES1-12': 2.15e-05,
 'ES1-18': 2.23e-05,
 'ES1-25': 2.2e-05,
 'ES1-26': 2.19e-05,
 'ES1-34': 2.6e-05,
 'ES1-35': 2.18e-05,
 'XMM_113': 1.9e-05,
 'XMM_z_0.81': 1.92e-05,
 'XMM_z_0.9': 1.9e-05,
 'XMM_z_1.0': 1.85e-05,
 'XMM-11': 2.3e-05,
 'XMM-19': 2.6e-05,
 'XMM-27': 2.2e-05,
 'XMM-29': 2.4e-05,
 'XMM-30': 2.21e-05}

def image_processing(image):
    hdu = fits.open(image)[0]
    header = hdu.header
    data = hdu.data[0,0,:,:]
    wcs = WCS(header, naxis = 2)
    data[np.where(np.isnan(data))] = 0
    return header, data, wcs

def plotting_img_psf(alma, psf, name):
    h1,d1,w1 = image_processing(alma)
    h2,d2,w2 = image_processing(psf)
    fig = plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(121, projection=w1)
    ax2 = fig.add_subplot(122, projection=w2, sharey = ax1)

    fwhm_maj1 = h1['BMAJ']*u.deg
    fwhm_min1 = h1['BMIN']*u.deg
    pa1 = h1['BPA'] + 90
    fwhm_maj2 = h2['BMAJ']*u.deg
    fwhm_min2 = h2['BMIN']*u.deg
    pa2 = h2['BPA'] + 90

    ra = h1['CRVAL1']*u.deg
    dec = h1['CRVAL2']*u.deg

    center = SkyCoord(ra, dec, frame = 'fk5')
    extent1 = SkyCoord(ra - fwhm_maj1, dec + fwhm_min1, frame='fk5')
    extent2 = SkyCoord(ra - fwhm_maj2, dec + fwhm_min2, frame='fk5')

    centerpix = w1.world_to_pixel(center)
    extpix1 = w1.world_to_pixel(extent1)
    extpix2 = w2.world_to_pixel(extent2)

    xlim = (centerpix[0]-100,centerpix[0]+100)
    ylim = (centerpix[1]-100,centerpix[1]+100)

    arad1, brad1 = int(extpix1[0] - centerpix[0]), int(extpix1[1] - centerpix[1])
    arad2, brad2 = int(extpix2[0] - centerpix[0]), int(extpix2[1] - centerpix[1])

    beam1 = Ellipse((xlim[0]+arad1*0.75, ylim[0]+brad1*0.75),width=arad1, height=brad1, angle=pa1, facecolor='yellow', edgecolor='white', zorder=11)
    beam2 = Ellipse((xlim[0]+arad2*0.75, ylim[0]+brad2*0.75),width=arad2, height=brad2, angle=pa2, facecolor='yellow', edgecolor='white', zorder=11)

    im1 = ax1.imshow(d1, origin = 'lower', interpolation='None',cmap = 'Greys_r',norm = simple_norm(d1, stretch = 'log', min_cut = 0, max_cut = 0.001))
    ax1.contour(d1, levels=np.array([2,3,4,5,6,7,8,9,10])*rms_dict[name], colors='red', linewidth = 1, zorder = 1, alpha = 0.5)
    im2 = ax2.imshow(d2, origin = 'lower', interpolation='None',cmap = 'Greys_r',norm = simple_norm(d2, stretch = 'log', min_cut = 0, max_cut = 10))
    ax1.add_patch(beam1)
    ax2.add_patch(beam2)
    ax1.set_ylabel("Dec")
    ax1.set_xlabel("RA")
    ax2.set_xlabel("RA")
    ax2.coords[1].set_ticks_visible(False)
    ax2.coords[1].set_ticklabel_visible(False)
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)
    fig.suptitle(name)
    fig.savefig("/Users/arames52/Desktop/Primary Project/postage_stamps/bcgs_psf/" + name + ".jpg", dpi = 300)

for alma in alma_images:
    name = alma.split("/")[-1].split("_ALMA.fits")[0]
    psf = "/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/ALMA_psfs/"+ name + ".fits"
    plotting_img_psf(alma, psf, name)
