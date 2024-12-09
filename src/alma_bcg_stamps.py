import numpy as np
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

BCGs = ['CDFS-18','CDFS-19','ES1-12','ES1-18','ES1-25','ES1-26','ES1-34','ES1-35',
'ES1_z_0.88','ES1_z_0.99','ES1_z_0.99b','ES1_z_1.04','ES1_z_1.38','ES1_z_1.40','ES1_z_1.60',
'ES1_z_1.65','ES1_z_1.70','XMM-11','XMM-19','XMM-27','XMM-29','XMM-30','XMM_113',
'XMM_z_0.81','XMM_z_0.9','XMM_z_1.0']
RMS = [0.0000218,0.0000218,0.0000215,0.0000203,0.0000192,0.0000227,0.0000217,0.0000197,
0.0000161,0.0000158,0.0000164,0.000015,0.0000154,0.0000165,0.0000158658,0.0000163,0.0000174,
0.0000219,0.0000202,0.000018,0.0000212,0.00002,0.0000175,0.0000168,0.0000178,0.000018]
S_max = [0.00386,0.0000898,0.0000645,0.00011,0.000141,0.0000681,0.000312,0.0000577,
0.0000606,0.0000943,0.0000876,0.000152,0.000221,0.000192,0.0001181054,0.000983,0.000288,
0.000114,0.0000606,0.000054,0.000147,0.00013,0.000106,0.0000874,0.0000942,0.00013]

# Dictionary of RMS values for making S/N contours
rms_dict = {}
smax_dict = {}
for key,value1, value2 in zip(BCGs, RMS, S_max):
    rms_dict[key] = value1
    smax_dict[key] = value2

alma_images = sorted(glob.glob("/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/all_ALMA/*.fits"))
sparcs_images = sorted(glob.glob("/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/SpARCS_z/*.fits"))

def image_processing(image_path, name, img_type):
    hdu = fits.open(image_path)[0]
    if img_type == 'alma':
        data = hdu.data[0,0,:,:]
        header = hdu.header
        wcs = WCS(header, naxis = 2)
    elif img_type == 'z':
        data = hdu.data
        header = hdu.header
        wcs = WCS(header)
    data[np.where(np.isnan(data))] = 0
    return data,header,wcs

fig = plt.figure(figsize = (20,20))
for alma in alma_images:
    name = alma.split("/")[-1].split("_ALMA.fits")[0]
    data, header, wcs = image_processing(alma, name, 'alma')
    ax = fig.add_subplot(4,7, alma_images.index(alma)+1, projection = wcs)
    fwhm_maj = header['BMAJ']*u.deg
    fwhm_min = header['BMIN']*u.deg
    pa = header['BPA'] + 90
    ra = header['CRVAL1']*u.deg
    dec = header['CRVAL2']*u.deg

    center = SkyCoord(ra, dec, frame = 'fk5')
    extent = SkyCoord(ra - fwhm_maj, dec + fwhm_min, frame='fk5')
    centerpix = wcs.world_to_pixel(center)
    extpix = wcs.world_to_pixel(extent)

    xlim = (centerpix[0]-50,centerpix[0]+50)
    ylim = (centerpix[1]-50,centerpix[1]+50)

    arad, brad = int(extpix[0] - centerpix[0]), int(extpix[1] - centerpix[1])
    beam = Ellipse((xlim[0]+arad*0.75, ylim[0]+brad*0.75),width=arad, height=brad, angle=pa, facecolor='yellow', edgecolor='white', zorder=11)
    levels = np.linspace(rms_dict[name]*2.5, smax_dict[name], 10)
    ax.imshow(data, origin = 'lower', cmap = 'Greys_r', interpolation = 'None', vmin = 0, vmax = 0.0001)
    ax.contour(data, levels = levels, colors='green', linewidths = 1, zorder = 1, alpha = 1)
    ax.add_patch(beam)
    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_ticks_visible(False)
    lon.set_ticklabel_visible(False)
    lat.set_ticks_visible(False)
    lat.set_ticklabel_visible(False)
    lon.set_axislabel('')
    lat.set_axislabel('')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

# plt.tight_layout()
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.1,
                    hspace=0.1)
fig.savefig("alma_greyscale_contours.jpg")


fig = plt.figure(figsize = (20,20))
for alma, sparcs in zip(alma_images, sparcs_images):
    name = alma.split("/")[-1].split("_ALMA.fits")[0]
    alma_data, alma_header, alma_wcs = image_processing(alma, name, 'alma')
    z_data, z_header, z_wcs = image_processing(sparcs, name, 'z')
    ax = fig.add_subplot(4,7, alma_images.index(alma)+1, projection = alma_wcs)
    fwhm_maj = alma_header['BMAJ']*u.deg
    fwhm_min = alma_header['BMIN']*u.deg
    pa = alma_header['BPA'] + 90
    ra = alma_header['CRVAL1']*u.deg
    dec = alma_header['CRVAL2']*u.deg

    center = SkyCoord(ra, dec, frame = 'fk5')
    extent = SkyCoord(ra - fwhm_maj, dec + fwhm_min, frame='fk5')
    centerpix = alma_wcs.world_to_pixel(center)
    extpix = alma_wcs.world_to_pixel(extent)

    xlim = (centerpix[0]-50,centerpix[0]+50)
    ylim = (centerpix[1]-50,centerpix[1]+50)

    arad, brad = int(extpix[0] - centerpix[0]), int(extpix[1] - centerpix[1])
    beam = Ellipse((xlim[0]+arad*0.75, ylim[0]+brad*0.75),width=arad, height=brad, angle=pa, facecolor='yellow', edgecolor='white', zorder=11)
    levels = np.linspace(rms_dict[name]*2.5, smax_dict[name], 10)

    if name.startswith("XMM"):
        vmax = 80
    else:
        vmax = 300
    ax.imshow(z_data, origin = 'lower', cmap = 'Greys_r', interpolation = 'None', vmin = 0, vmax = vmax)
    trans = ax.get_transform(z_wcs)
    ax.contour(alma_data, levels = levels, colors='red', linewidths = 1, zorder = 1, alpha = 1)

    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_ticks_visible(False)
    lon.set_ticklabel_visible(False)
    lat.set_ticks_visible(False)
    lat.set_ticklabel_visible(False)
    lon.set_axislabel('')
    lat.set_axislabel('')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
