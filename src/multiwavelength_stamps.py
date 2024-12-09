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
import pandas as pd
from astropy.table import Table
from astropy.visualization import simple_norm, ZScaleInterval
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization.stretch import LinearStretch, LogStretch
import warnings
warnings.filterwarnings('ignore')
plt.style.use(["science", "stylesheet.txt"])


# Storing the image paths as list
alma_images = sorted(glob.glob("/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/all_ALMA/*.fits"))
sparcs_images = sorted(glob.glob("/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/SpARCS_z/*.fits"))
mips_images = sorted(glob.glob("/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/IRAC/mips_24_micron/*.fits"))
irac_images = sorted(glob.glob("/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/IRAC/irac_3.6_micron/*.fits"))

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
contour_steps = {'CDFS-18':25, 'CDFS-19':1, 'ES1-12':1,'ES1-26':1, 'ES1-35':1, 'XMM-19':1, 'XMM-27':1,
                'ES1-18':1, 'ES1-25':1, 'ES1_z_0.88':1, 'ES1_z_0.99':1, 'ES1_z_0.99b':1,
                'ES1_z_1.60':1, 'XMM-11':1, 'XMM-29':1, 'XMM-30':1,'XMM_113':1, 'XMM_z_0.81':1,
                'XMM_z_0.9':1, 'XMM_z_1.0':1, 'ES1-34':2, 'ES1_z_1.04':2, 'ES1_z_1.38':2, 'ES1_z_1.40':2,
                'ES1_z_1.70':4, 'ES1_z_1.65':10}
redshift = [0.813, 0.69484, 0.74, 0.8068, 1.08807, 0.85192, 1.7, 1.7, 0.56395, 1.19372, 0.9188, 1.04,
            1.38, 1.4, 1.6, 1.65, 1.7, 0.79, 1.04789, 1.38, 1.45, 1.45, 1.6, 0.7827207, 0.84957, 0.5352355]
rms_dict = {}
smax_dict = {}
z_dict = {}
for key,value1, value2, value3 in zip(BCGs, RMS, S_max, redshift):
    rms_dict[key] = value1
    smax_dict[key] = value2
    z_dict[key] = value3

irac_sources_path = "/Users/arames52/Desktop/CASA Imaging/irac_sources/"
mips_sources_path = "/Users/arames52/Desktop/CASA Imaging/mips_sources/"
hermes_path = "/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/Herschel/catalog/"

def image_processing(image, img_type = None):
    hdu = fits.open(image)[0]
    if img_type == 'alma':
        header = hdu.header
        data = hdu.data[0,0,:,:]
        wcs = WCS(header, naxis = 2)
    else:
        header = hdu.header
        data = hdu.data
        wcs = WCS(header)
    data[np.where(np.isnan(data))] = 0
    #         data /= np.amax(data)
    return header, data, wcs

def plotting_stamps(z_band, irac, mips, alma, name):
    header_z, data_z, wcs_z = image_processing(z_band)
    header_irac, data_irac, wcs_irac = image_processing(irac)
    header_mips, data_mips, wcs_mips = image_processing(mips)
    header_alma, data_alma, wcs_alma = image_processing(alma, 'alma')
    fig = plt.figure(dpi = 300)
    ax1 = fig.add_subplot(141, projection=wcs_alma)
    ax2 = fig.add_subplot(142, projection=wcs_alma,sharey=ax1)
    ax3 = fig.add_subplot(143, projection=wcs_alma,sharey=ax1)
    ax4 = fig.add_subplot(144, projection=wcs_alma,sharey=ax1)
    z_trans = ax1.get_transform(wcs_z)
    irac_trans = ax2.get_transform(wcs_irac)
    mips_trans = ax3.get_transform(wcs_mips)
#    bmaj = header_alma['BMAJ']*u.deg.to(u.arcsec)/0.045
#    bmin = header_alma['BMIN']*u.deg.to(u.arcsec)/0.045
    fwhm_maj = header_alma['BMAJ']*u.deg
    fwhm_min = header_alma['BMIN']*u.deg
    pa = header_alma['BPA'] + 90
    ra = header_alma['CRVAL1']
    dec = header_alma['CRVAL2']
    center = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5')
    extent = SkyCoord(ra*u.deg - fwhm_maj, dec*u.deg + fwhm_min, frame='fk5')
    centerpix = wcs_alma.world_to_pixel(center)
    extpix = wcs_alma.world_to_pixel(extent)
    xlim = (centerpix[0]-150,centerpix[0]+150)
    ylim = (centerpix[1]-150,centerpix[1]+150)
    arad, brad = (extpix[0] - centerpix[0]), (extpix[1] - centerpix[1])
    beam1 = Ellipse((xlim[0]+arad*0.75, ylim[0]+brad*0.75),width=bmaj, height=bmin, angle=pa, facecolor='yellow', edgecolor='white', zorder=11, lw = 0.2)
    beam2 = Ellipse((xlim[0]+arad*0.75, ylim[0]+brad*0.75),width=bmaj, height=bmin, angle=pa, facecolor='yellow', edgecolor='white', zorder=11, lw = 0.2)
    beam3 = Ellipse((xlim[0]+arad*0.75, ylim[0]+brad*0.75),width=bmaj, height=bmin, angle=pa, facecolor='yellow', edgecolor='white', zorder=11, lw = 0.2)
    beam4 = Ellipse((xlim[0]+arad*0.75, ylim[0]+brad*0.75),width=bmaj, height=bmin, angle=pa, facecolor='yellow', edgecolor='white', zorder=11, lw = 0.2)
    if name.startswith("XMM"):
        vmax = 80
    else:
        vmax = 300
    levels = np.arange(rms_dict[name]*3, smax_dict[name], contour_steps[name]*rms_dict[name], dtype=None)

    # irac sources
    if irac_sources_path + name + ".csv" not in sorted(glob.glob(irac_sources_path + "*.csv")):
        irac_sources_coord = []
    else:
        irac_df = pd.read_csv(irac_sources_path + name + ".csv")
        irac_sources_coord = list(zip(irac_df.ra, irac_df.dec))
        
        
    c = SkyCoord(ra*u.deg, dec*u.deg)
    
    if name.startswith("ES1"):
        file = "L6-ELAIS-S1-SWIRE_xID24_DR3.fits"
    elif name.startswith("XMM"):
        file = "L6-XMM-LSS-SWIRE_xID24_DR3.fits"
    elif name.startswith("CDFS"):
        file = "L5-CDFS-SWIRE_xID24_DR3.fits"
    _f = fits.open(hermes_path + file)
    _hermes = Table(_f[1].data).to_pandas()
    _h_catalog = SkyCoord(ra=np.array(_hermes['RA'])*u.degree, dec=np.array(_hermes['Dec'])*u.degree, frame = 'icrs')
    idx = c.separation(_h_catalog) < 20*u.arcsec
    c_matches = _h_catalog[idx]

    steps = contour_steps[name]
    im1 = ax1.imshow(data_z, transform = z_trans, origin = 'lower', interpolation='None',cmap = 'viridis',vmin = 0, vmax = vmax)
    ax1.contour(data_alma, levels=levels, colors='red', linewidths = 0.5, zorder = 1)
    im2 = ax2.imshow(data_irac, transform = irac_trans, origin = 'lower', interpolation='None',cmap = 'viridis',vmin = 0, vmax = 1)
    
    ax2.contour(data_alma, levels=levels, colors='red', linewidths = 0.5, zorder = 1)
    for coord in irac_sources_coord:
        c = SkyCoord(coord[0]*u.deg, coord[1]*u.deg, frame = 'fk5')
        c_pix = wcs_irac.world_to_pixel(c)
        ax2.scatter(c_pix[0], c_pix[1], marker = '+', color = 'blue', zorder = 5, transform = irac_trans, s = 10, lw= 0.8, label = 'irac sources')
    im3 = ax3.imshow(data_mips, transform = mips_trans, origin = 'lower', interpolation='None',cmap = 'viridis',vmin = 0, vmax = 0.5)
    
    ax3.contour(data_alma, levels=levels, colors='red', linewidths = 0.5, zorder = 1)
    for coord in c_matches:
        c_pix = wcs_mips.world_to_pixel(coord)
        ax3.scatter(c_pix[0], c_pix[1], marker = '+', color = 'blue', zorder = 5, transform = mips_trans, s = 10, lw= 0.8, label = 'mips sources')
    im4 = ax4.imshow(data_alma, origin = 'lower', interpolation='None',cmap = 'viridis',vmin = 0, vmax = 0.0001)
    ax4.contour(data_alma, levels=levels, colors = 'red', linewidths = 0.5, zorder = 1, alpha = 1)
    text = str(round(steps,2)) + "$\sigma$"
    ax4.text(1, 0.01, text,
            horizontalalignment='right', verticalalignment='bottom',
            transform=ax4.transAxes, color = 'white', alpha = 1)

    ax1.add_patch(beam1)
    ax2.add_patch(beam2)
    ax3.add_patch(beam3)
    ax4.add_patch(beam4)
    ax1.set_ylabel("Dec")
    ax1.set_xlabel("RA")
    ax2.set_xlabel("RA")
    ax3.set_xlabel("RA")
    ax4.set_xlabel("RA")

    ax1.coords[1].set_ticks_visible(False)
    ax1.coords[1].set_ticklabel_visible(False)
    ax2.coords[1].set_ticks_visible(False)
    ax2.coords[1].set_ticklabel_visible(False)
    ax3.coords[1].set_ticks_visible(False)
    ax3.coords[1].set_ticklabel_visible(False)
    ax4.coords[1].set_ticks_visible(False)
    ax4.coords[1].set_ticklabel_visible(False)
    ax1.coords[0].set_ticks_visible(False)
    ax1.coords[0].set_ticklabel_visible(False)
    ax2.coords[0].set_ticks_visible(False)
    ax2.coords[0].set_ticklabel_visible(False)
    ax3.coords[0].set_ticks_visible(False)
    ax3.coords[0].set_ticklabel_visible(False)
    ax4.coords[0].set_ticks_visible(False)
    ax4.coords[0].set_ticklabel_visible(False)
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)
    ax3.set_xlim(xlim)
    ax3.set_ylim(ylim)
    ax4.set_xlim(xlim)
    ax4.set_ylim(ylim)
    ax1.set_title("SpARCS z-band", size = 10)
    ax2.set_title("IRAC 3.6 $\mu$m", size = 10)
    ax3.set_title("MIPS 24 $\mu$m", size = 10)
    ax4.set_title("ALMA 1.2 mm", size = 10)
    # fig.suptitle(name, y = 0.70)
    return fig



for z, irac, mips, alma in zip(sparcs_images, irac_images, mips_images, alma_images):
    name = alma.split("/")[-1].split(".fits")[0]
    print(name)
    fig = plotting_stamps(z, irac, mips, alma, name)
    fig.savefig("/Users/arames52/Desktop/postage_stamps/" + name + "_stamps.jpg",bbox_inches ='tight', dpi = 300)
