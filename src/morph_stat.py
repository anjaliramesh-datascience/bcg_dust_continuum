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


def morphology(path, size):
    name = path.split("/")[-1].split("_ALMA.fits")[0]
    hdu = fits.open(path)[0]
    header = hdu.header
    data = hdu.data[0,0,:,:]
    wcs = WCS(header, naxis = 2)
    data[np.where(np.isnan(data))] = 0
    bmaj = header['BMAJ']*u.deg.to(u.arcsec)/0.045
    bmin = header['BMIN']*u.deg.to(u.arcsec)/0.045
    pa = (header['BPA'] + 90)*u.deg.to(u.rad)
    position = SkyCoord(header['CRVAL1']*u.deg, header['CRVAL2']*u.deg, frame = 'fk5')
    image = Cutout2D(data, position = position, size =size, wcs = wcs).data
    rms = rms_dict[name]
    weight_map = np.full(image.shape, rms)
    psf = Gaussian2DKernel(x_stddev = gaussian_fwhm_to_sigma*bmin, y_stddev = gaussian_fwhm_to_sigma*bmaj, theta = pa, x_size = image.shape[0], y_size = image.shape[1]).array
    threshold = photutils.detect_threshold(image, 1.5, background = None)
    npixels = 5
    segm = photutils.detect_sources(image, threshold, npixels)
    label = np.argmax(segm.areas) + 1
    segmap = segm.data == label
    segmap_float = ndi.uniform_filter(np.float64(segmap), size=20)
    segmap = segmap > 0.5
    source_morphs = statmorph.source_morphology(image, segmap, weightmap = weight_map, psf = psf)
    morph = source_morphs[0]
    fig = make_figure(morph)
    return morph
