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
import warnings
warnings.filterwarnings('ignore')
plt.style.use(["default"])
#rms of alma bcg images
rms_dict = {'CDFS-18': 2.54e-05,'CDFS-19': 2.37e-05,'ES1_z_0.88': 1.79e-05,'ES1_z_0.99': 1.78e-05,'ES1_z_0.99b': 1.79e-05,'ES1_z_1.04': 1.69e-05,'ES1_z_1.38': 1.71e-05,'ES1_z_1.40': 1.65e-05,'ES1_z_1.60': 1.68e-05,'ES1_z_1.65': 1.72e-05,'ES1_z_1.70': 1.75e-05,'ES1-12': 2.15e-05,'ES1-18': 2.23e-05,'ES1-25': 2.2e-05,'ES1-26': 2.19e-05,'ES1-34': 2.6e-05,'ES1-35': 2.18e-05,'XMM_113': 1.9e-05,'XMM_z_0.81': 1.92e-05,'XMM_z_0.9': 1.9e-05,'XMM_z_1.0': 1.85e-05,'XMM-11': 2.3e-05,'XMM-19': 2.6e-05,'XMM-27': 2.2e-05,'XMM-29': 2.4e-05,'XMM-30': 2.21e-05}

 # sparcs z band coordinates
z_coords = {'CDFS-18': {'ra': 53.31966, 'dec': -26.89633},'CDFS-19': {'ra': 53.48372, 'dec': -27.25943},'XMM-11': {'ra': 36.0999, 'dec': -2.96742},'XMM-19': {'ra': 35.7742, 'dec': -4.22653},'XMM-27': {'ra': 36.56509, 'dec': -4.94118},'XMM-29': {'ra': 36.87717, 'dec': -4.53427},'XMM-30': {'ra': 34.68846, 'dec': -5.71579},'ES1-12': {'ra': 9.76266, 'dec': -44.36449},'ES1-18': {'ra': 7.03082, 'dec': -43.34483},'ES1-25': {'ra': 6.97814, 'dec': -43.2798},'ES1-26': {'ra': 8.92828, 'dec': -43.04841},'ES1-34': {'ra': 9.21062, 'dec': -42.11781},'ES1-35': {'ra': 9.34235, 'dec': -43.50759},'XMM_z_0.81': {'ra': 34.87042, 'dec': -4.11672},'XMM_z_0.9': {'ra': 36.1217, 'dec': -4.17031},'XMM_z_1.0': {'ra': 35.6521, 'dec': -3.84174},'XMM_113': {'ra': 36.10969, 'dec': -3.39178},'ES1_z_0.88': {'ra': 8.81961, 'dec': -42.68544},'ES1_z_0.99': {'ra': 9.41931, 'dec': -43.6527},'ES1_z_0.99b': {'ra': 8.68045, 'dec': -42.45776},'ES1_z_1.04': {'ra': 8.8564, 'dec': -42.23985},'ES1_z_1.38': {'ra': 8.64012, 'dec': -42.49283},'ES1_z_1.40': {'ra': 8.29814, 'dec': -43.49844},'ES1_z_1.60': {'ra': 9.49156, 'dec': -42.98388},'ES1_z_1.65': {'ra': 7.96871, 'dec': -44.43894},'ES1_z_1.70': {'ra': 7.34152, 'dec': -43.44523}}

# Storing the image paths as list
alma_images = sorted(glob.glob("/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/all_ALMA/*.fits"))
sparcs_images = sorted(glob.glob("/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/SpARCS_z/*.fits"))
mips_images = sorted(glob.glob("/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/IRAC/mips_24_micron/*.fits"))
irac_images = sorted(glob.glob("/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/IRAC/irac_3.6_micron/*.fits"))

# Function for creating the postage stamps
def plotting_stamps(z_band, irac, mips, alma, name):
    """
    z_band : path to z-band image
    irac: path to irac 3.6 micron image
    mips : path to 24 micron mips image
    alma : path to alma image
    """
    # opening all the fits images and processing
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

    bmaj = header_alma['BMAJ']*u.deg.to(u.arcsec)/0.045
    bmin = header_alma['BMIN']*u.deg.to(u.arcsec)/0.045
    fwhm_maj = header_alma['BMAJ']*u.deg
    fwhm_min = header_alma['BMIN']*u.deg
    pa = header_alma['BPA'] + 90
    ra = header_alma['CRVAL1']*u.deg
    dec = header_alma['CRVAL2']*u.deg

    center = SkyCoord(ra, dec, frame = 'fk5')
    extent = SkyCoord(ra - fwhm_maj, dec + fwhm_min, frame='fk5')
    centerpix = wcs_alma.world_to_pixel(center)
    extpix = wcs_alma.world_to_pixel(extent)

    xlim = (centerpix[0]-100,centerpix[0]+100)
    ylim = (centerpix[1]-100,centerpix[1]+100)

    arad, brad = (extpix[0] - centerpix[0]), (extpix[1] - centerpix[1])
    beam1 = Ellipse((xlim[0]+arad*0.75, ylim[0]+brad*0.75),width=bmaj, height=bmin, angle=pa, facecolor='yellow', edgecolor='white', zorder=11, lw = 0.2)
    beam2 = Ellipse((xlim[0]+arad*0.75, ylim[0]+brad*0.75),width=bmaj, height=bmin, angle=pa, facecolor='yellow', edgecolor='white', zorder=11, lw = 0.2)
    beam3 = Ellipse((xlim[0]+arad*0.75, ylim[0]+brad*0.75),width=bmaj, height=bmin, angle=pa, facecolor='yellow', edgecolor='white', zorder=11, lw = 0.2)
    beam4 = Ellipse((xlim[0]+arad*0.75, ylim[0]+brad*0.75),width=bmaj, height=bmin, angle=pa, facecolor='yellow', edgecolor='white', zorder=11, lw = 0.2)
    if name.startswith("XMM"):
        vmax = 80
    else:
        vmax = 300
    z_coord = SkyCoord(ra = z_coords[name]['ra']*u.degree, dec = z_coords[name]['dec']*u.degree, frame = 'icrs')
    z_pix = wcs_alma.world_to_pixel(z_coord)
    im1 = ax1.imshow(data_z, transform = z_trans, origin = 'lower', interpolation='None',cmap = 'Greys_r',vmin = 0, vmax = vmax)
    ax1.contour(data_alma, levels=np.array([2,3,4,5,6,7,8,9,10])*rms_dict[name], colors='red', linewidths = 0.5, zorder = 1)
    # ax1.scatter(centerpix[0], centerpix[1], marker = '+', color = 'blue', linewidth = 1, label = "3.6 micron")
    # ax1.scatter(z_pix[0], z_pix[1], marker = '+', color = 'green', zorder = 4, label = "Z band")

    im2 = ax2.imshow(data_irac, transform = irac_trans, origin = 'lower', interpolation='None',cmap = 'Greys_r',vmin = 0, vmax = 1)
    ax2.contour(data_alma, levels=np.array([2,3,4,5,6,7,8,9,10])*rms_dict[name], colors='red', linewidths = 0.5, zorder = 1)
    # ax2.scatter(centerpix[0], centerpix[1], marker = '+', color = 'blue', linewidth = 1)

    im3 = ax3.imshow(data_mips, transform = mips_trans, origin = 'lower', interpolation='None',cmap = 'Greys_r',vmin = 0, vmax = 0.3)
    ax3.contour(data_alma, levels=np.array([2,3,4,5,6,7,8,9,10])*rms_dict[name], colors='red', linewidths = 0.5, zorder = 1)
    # ax3.scatter(centerpix[0], centerpix[1], marker = '+', color = 'blue', linewidth = 1)

    im4 = ax4.imshow(data_alma, origin = 'lower', interpolation='None',cmap = 'Greys_r',vmin = 0, vmax = 0.0001)
    ax4.contour(data_alma, levels=np.array([2,3,4,5,6,7,8,9,10])*rms_dict[name], colors = 'red', linewidths = 0.5, zorder = 1, alpha = 1)
    # ax4.scatter(centerpix[0], centerpix[1], marker = '+', color = 'blue', linewidth = 1)

    # legend = ax1.legend()
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
    # ax1.coords[0].set_axislabel('')
    # ax2.coords[0].set_axislabel('')
    # ax3.coords[0].set_axislabel('')
    # ax4.coords[0].set_axislabel('')
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
    fig.suptitle(name, y = 0.70)
    return fig


for z, irac, mips, alma in zip(sparcs_images, irac_images, mips_images, alma_images):
    name = alma.split("/")[-1].split("_ALMA.fits")[0]
    fig = plotting_stamps(z, irac, mips, alma, name)
    fig.savefig("/Users/arames52/Desktop/Primary Project/postage_stamps/" + name + "_alma.jpg",bbox_inches ='tight', dpi = 300)
    break
