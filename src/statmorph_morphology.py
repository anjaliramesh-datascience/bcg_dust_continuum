import glob
import statmorph
import photutils
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import fits
import scipy.ndimage as ndi
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from matplotlib.patches import Rectangle, Ellipse, Circle
from astropy.visualization import simple_norm, ZScaleInterval
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization.stretch import LinearStretch, LogStretch
from statmorph.utils.image_diagnostics import make_figure

rms_dict = {'CDFS-18': 2.54e-05,'CDFS-19': 2.37e-05,'ES1_z_0.88': 1.79e-05,'ES1_z_0.99': 1.78e-05,'ES1_z_0.99b': 1.79e-05,'ES1_z_1.04': 1.69e-05,'ES1_z_1.38': 1.71e-05,'ES1_z_1.40': 1.65e-05,'ES1_z_1.60': 1.68e-05,'ES1_z_1.65': 1.72e-05,'ES1_z_1.70': 1.75e-05,'ES1-12': 2.15e-05,'ES1-18': 2.23e-05,'ES1-25': 2.2e-05,'ES1-26': 2.19e-05,'ES1-34': 2.6e-05,'ES1-35': 2.18e-05,'XMM_113': 1.9e-05,'XMM_z_0.81': 1.92e-05,'XMM_z_0.9': 1.9e-05,'XMM_z_1.0': 1.85e-05,'XMM-11': 2.3e-05,'XMM-19': 2.6e-05,'XMM-27': 2.2e-05,'XMM-29': 2.4e-05,'XMM-30': 2.21e-05}
size_dict = {'ES1-34': (150, 150),'CDFS-18': (150, 150),'ES1-18': (100, 100),'ES1_z_0.88': (75, 75),'ES1_z_0.99': (100, 100),
'ES1_z_0.99b': (60, 60),'ES1_z_1.04': (150, 150),'ES1_z_1.38': (75, 75),'ES1_z_1.40': (100, 100),'ES1_z_1.60': (150, 150),
'ES1_z_1.65': (100, 100),'ES1_z_1.70': (150, 150),'XMM-11': (100, 100),'XMM-29': (150, 150),'XMM_113': (100, 100),
'XMM_z_0.81': (75, 75),'XMM_z_0.9': (100, 100),'XMM_z_1.0': (130, 130)}

alma_images = sorted(glob.glob("/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/ALMA/*.fits"))

sersic_n = {}
r20 = {}
flag_sersic = {}
morph_nopsf = {}
morph_psf = {}

def image_processing(path, name):
    """
    input : path to fits image
    returns : 8"x8" cutout of image, clean beam constructed from major and minor FWHM
    """
    hdu = fits.open(path)[0]
    header = hdu.header
    data = hdu.data[0,0,:,:]
    wcs = WCS(header, naxis = 2)
    data[np.where(np.isnan(data))] = 0

    position = SkyCoord(header['CRVAL1']*u.deg, header['CRVAL2']*u.deg, frame = 'fk5')
    bmaj = header['BMAJ']*u.deg.to(u.arcsec)/0.045
    bmin = header['BMIN']*u.deg.to(u.arcsec)/0.045
    pa = (header['BPA'] + 90)*u.deg.to(u.rad)

    size = size_dict[name]

    image = Cutout2D(data, position = position, size = size, wcs = wcs).data
    weightmap = np.full(image.shape, rms_dict[name])
    psf = Gaussian2DKernel(x_stddev = gaussian_fwhm_to_sigma*bmin, y_stddev = gaussian_fwhm_to_sigma*bmaj, theta = pa, x_size = image.shape[0], y_size = image.shape[1]).array

    return image, weightmap, psf

def compute_sersic(image, weightmap, psf):

    threshold = photutils.detect_threshold(image, 1.5, background = None)
    npixels = 5
    segm = photutils.detect_sources(image, threshold, npixels)
    label = np.argmax(segm.areas) + 1
    segmap = segm.data == label
    segmap_float = ndi.uniform_filter(np.float64(segmap), size=20)
    segmap = segmap > 0.5
    source_morphs = statmorph.source_morphology(image, segmap, weightmap = weightmap, psf = psf)
    morph = source_morphs[0]
    sersic_n[name] = morph.sersic_n
    flag_sersic[name] = morph.flag_sersic
    r20[name] = morph.r20

    return morph


def morphology_without_psf(image, weightmap):
    """
    Running statmorph.
    1) Detecting source with photutils
    2) Creating segmentation map
    3) Calling statmorph - input : cutout image, segmap and gain = 1
    """
    threshold = photutils.detect_threshold(image, 1.5, background = 0)
    npixels = 5
    segm = photutils.detect_sources(image, threshold, npixels)
    label = np.argmax(segm.areas) + 1
    segmap = segm.data == label
    segmap_float = ndi.uniform_filter(np.float64(segmap), size=20)
    segmap = segmap_float > 0.5
    source_morphs = statmorph.source_morphology(image, segmap, weightmap = weightmap)
    morph = source_morphs[0]

    return morph

def plotting_sersic_profile(morph, psf, name):

    # With PSF
    theta_vec = np.linspace(0.0, 2.0*np.pi, 200)
    image = np.float64(morph._cutout_stamp_maskzeroed)  # skimage wants double
    ny, nx = image.shape
    xc, yc = morph._xc_stamp, morph._yc_stamp  # centroid
    y, x = np.mgrid[0:ny, 0:nx]
    xcs, ycs = morph._sersic_model.x_0.value, morph._sersic_model.y_0.value
    sersic_model = morph._sersic_model(x, y)
    residual = image - sersic_model
    if morph.sky_sigma > 0:
        sersic_model += np.random.normal(scale=morph.sky_sigma, size=(ny, nx))
    sersic_model = np.reshape(sersic_model, (ny,nx))
    fig,ax = plt.subplots(1,5, figsize = (10,5))
    ax[0].imshow(image, origin = 'lower', cmap = 'Greys_r', vmin = 0, vmax = 0.2)
    ax[0].set_title(name, size = 10)
    ax[1].imshow(psf, origin = 'lower', cmap = 'Greys_r', vmin = 0, vmax = 0.001)
    ax[1].set_title("PSF", size = 10)
    ax[3].imshow(sersic_model, origin = 'lower', cmap = 'Greys_r', vmin = 0, vmax = 0.2)
    ax[3].set_title("Sersic Model + PSF", size = 10)
    ax[4].imshow(residual, origin = 'lower', cmap = 'Greys_r', vmin = 0, vmax = 0.2)
    ax[4].set_title("Residual", size = 10)
    ax[3].plot(xcs, ycs, 'bo', markersize=1, label='Sérsic Center')
    R = float(nx**2 + ny**2)
    theta = morph.sersic_theta
    x0, x1 = xcs - R*np.cos(theta), xcs + R*np.cos(theta)
    y0, y1 = ycs - R*np.sin(theta), ycs + R*np.sin(theta)
    ax[3].plot([x0, x1], [y0, y1], 'r--', lw=1, label='Major Axis (Sérsic)')
    # Half-radius ellipse
    a = morph.sersic_rhalf
    b = a * (1.0 - morph.sersic_ellip)
    xprime, yprime = a*np.cos(theta_vec), b*np.sin(theta_vec)
    x = xcs + (xprime*np.cos(theta) - yprime*np.sin(theta))
    y = ycs + (xprime*np.sin(theta) + yprime*np.cos(theta))
    ax[3].plot(x, y, 'r', label='Half-Light Ellipse (Sérsic)')
    ax[3].set_xlim(-0.5, nx-0.5)
    ax[3].set_ylim(-0.5, ny-0.5)

    # Without PSF
    morph1 = morph_nopsf[name]
    image = np.float64(morph1._cutout_stamp_maskzeroed)  # skimage wants double
    ny, nx = image.shape
    xc, yc = morph1._xc_stamp, morph1._yc_stamp  # centroid
    y, x = np.mgrid[0:ny, 0:nx]
    xcs, ycs = morph1._sersic_model.x_0.value, morph1._sersic_model.y_0.value
    sersic_model = morph1._sersic_model(x, y)
    residual = image - sersic_model
    if morph1.sky_sigma > 0:
        sersic_model += np.random.normal(scale=morph1.sky_sigma, size=(ny, nx))
    sersic_model = np.reshape(sersic_model, (ny,nx))

    ax[2].imshow(sersic_model, origin = 'lower', cmap = 'Greys_r', vmin = 0, vmax = 0.2)
    ax[2].set_title("Sersic Model", size = 10)
    ax[2].plot(xcs, ycs, 'bo', markersize=1, label='Sérsic Center')
    R = float(nx**2 + ny**2)
    theta = morph1.sersic_theta
    x0, x1 = xcs - R*np.cos(theta), xcs + R*np.cos(theta)
    y0, y1 = ycs - R*np.sin(theta), ycs + R*np.sin(theta)
    ax[2].plot([x0, x1], [y0, y1], 'r--', lw=1, label='Major Axis (Sérsic)')
    # Half-radius ellipse
    a = morph1.sersic_rhalf
    b = a * (1.0 - morph1.sersic_ellip)
    xprime, yprime = a*np.cos(theta_vec), b*np.sin(theta_vec)
    x = xcs + (xprime*np.cos(theta) - yprime*np.sin(theta))
    y = ycs + (xprime*np.sin(theta) + yprime*np.cos(theta))
    ax[2].plot(x, y, 'r', label='Half-Light Ellipse (Sérsic)')
    ax[2].set_xlim(-0.5, nx-0.5)
    ax[2].set_ylim(-0.5, ny-0.5)

    ax[3].legend(fontsize=5)

    for axs in ax:
        axs.set_xticks([])
        axs.set_yticks([])
    return fig
