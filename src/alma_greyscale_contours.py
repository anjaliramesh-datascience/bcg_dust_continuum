import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Ellipse, Circle, Arrow
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

redshift = [0.813, 0.69484, 0.74, 0.8068, 1.08807, 0.85192, 1.7, 1.7, 0.56395, 1.19372, 0.9188, 1.04,
            1.38, 1.4, 1.6, 1.65, 1.7, 0.79, 1.04789, 1.38, 1.45, 1.45, 1.6, 0.7827207, 0.84957, 0.5352355]
rms_dict = {}
smax_dict = {}
z_dict = {}
for key,value1, value2, value3 in zip(BCGs, RMS, S_max, redshift):
    rms_dict[key] = value1
    smax_dict[key] = value2
    z_dict[key] = value3

alma_images = sorted(glob.glob("ALMA_fits_files/*.fits"))
contour_steps = {'CDFS-18':25, 'CDFS-19':1, 'ES1-12':1,'ES1-26':1, 'ES1-35':1, 'XMM-19':1, 'XMM-27':1,
                'ES1-18':1, 'ES1-25':1, 'ES1_z_0.88':1, 'ES1_z_0.99':1, 'ES1_z_0.99b':1,
                'ES1_z_1.60':1, 'XMM-11':1, 'XMM-29':1, 'XMM-30':1,'XMM_113':1, 'XMM_z_0.81':1,
                'XMM_z_0.9':1, 'XMM_z_1.0':1, 'ES1-34':2, 'ES1_z_1.04':2, 'ES1_z_1.38':2, 'ES1_z_1.40':2,
                'ES1_z_1.70':4, 'ES1_z_1.65':10}
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

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

fig = plt.figure(figsize = (18,10))
for alma in alma_images:
    name = alma.split(".fits")[0].split("ALMA_fits_files/")[-1]
    data, header, wcs = image_processing(alma, name, 'alma')
    ax = fig.add_subplot(4,7, alma_images.index(alma)+1, projection = wcs)

    bmaj = header['BMAJ']*u.deg.to(u.arcsec)/0.045
    bmin = header['BMIN']*u.deg.to(u.arcsec)/0.045
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
    beam = Ellipse((xlim[0]+arad, ylim[0]+brad),width=bmaj, height=bmin, angle=pa, facecolor='yellow', edgecolor='white', zorder=11)

    levels = np.arange(rms_dict[name]*2.5, smax_dict[name], contour_steps[name]*rms_dict[name], dtype=None)
    steps = contour_steps[name]

    length = 10*(1/(0.045*(cosmo.kpc_proper_per_arcmin(z_dict[name])/60))).value
    kpc_line = Arrow(x = xlim[0] + 5, y = ylim[1]-10, dx = length, dy = 0.05, width = 0.5, color = 'yellow')
    line_text = "10 kpc"
    ax.text(0.1, 0.98, line_text, horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, color = 'white', alpha = 1)
    text = str(round(steps,2)) + "$\sigma$"
    ax.imshow(data, origin = 'lower', cmap = 'Greys_r', interpolation = 'None', vmin = 0, vmax = 0.0001)
    ax.contour(data, levels = levels, colors='green', linewidths = 1, zorder = 1, alpha = 1)
    ax.add_patch(beam)
    ax.add_patch(kpc_line)
    ax.text(1, 0.01, text,
            horizontalalignment='right', verticalalignment='bottom',
            transform=ax.transAxes, color = 'white', alpha = 1)
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
                    wspace=0.03,
                    hspace=0.05)
fig.savefig("alma_postage_stamps.png", dpi = 400, bbox_inches = 'tight')
