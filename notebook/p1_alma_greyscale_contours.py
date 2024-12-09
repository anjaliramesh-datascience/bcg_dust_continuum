import pandas as pd
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
from astropy.cosmology import FlatLambdaCDM
import os

cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
path = "/Volumes/T7/BCG_Dust_Continuum/"
os.chdir(path)

alma_image_directory = 'data/tclean_output/natural_imaging/'
bcg_info = pd.read_csv("data/bcg_info/bcg_basic_data.csv", header=0)

file_ids = bcg_info['file_id']
rms = bcg_info['RMS']
sn = bcg_info['S/N']
step = bcg_info['contour_steps']
z = bcg_info['redshift']
spec_z = bcg_info['spec_z']

def get_data(file_id):
    
    image_path = alma_image_directory + file_id + "_natural.fits"
    hdu = fits.open(image_path)[0]
    header = hdu.header
    data = hdu.data[0,0,:,:]
    wcs = WCS(header, naxis = 2)
    data = np.nan_to_num(data)
    
    return data,header,wcs

def get_beam(header):

    wcs = WCS(header, naxis=2)
    bmaj = header['BMAJ']*u.deg.to(u.arcsec)/0.045
    bmin = header['BMIN']*u.deg.to(u.arcsec)/0.045
    fwhm_maj = header['BMAJ']*u.deg
    fwhm_min = header['BMIN']*u.deg
    pa = header['BPA'] + 90
    ra = header['CRVAL1']*u.deg
    dec = header['CRVAL2']*u.deg
    center = SkyCoord(ra, dec, frame='fk5')
    extent = SkyCoord(ra - fwhm_maj, dec + fwhm_min, frame='fk5')
    centerpix = wcs.world_to_pixel(center)
    extpix = wcs.world_to_pixel(extent)

    xlim = (centerpix[0]-50, centerpix[0]+50)
    ylim = (centerpix[1]-50, centerpix[1]+50)

    arad, brad = int(extpix[0] - centerpix[0]), int(extpix[1] - centerpix[1])
    beam = Ellipse((xlim[0]+arad, ylim[0]+brad),width=bmaj, height=bmin, angle=pa, facecolor='yellow', edgecolor='white', zorder=11)

    return xlim, ylim, beam

def get_kpc_line(xlim, ylim, z):

    length = 10*(1/(0.045*(cosmo.kpc_proper_per_arcmin(z)/60))).value
    kpc_line = Arrow(x = xlim[0] + 5, y = ylim[1]-10, dx = length, dy = 0.05, width = 0.5, color = 'yellow')

    return kpc_line

def get_contour_color(specz):

    if specz == 0:
        contour_color = 'red'
    else:
        contour_color = 'green'

    return contour_color

def get_contour_levels(rms, sn, step):

    return np.arange(rms*2.5, sn*rms, step*rms, dtype=None)

def make_postage_stamps():

    # Create a figure with a grid of subplots
    fig, axs = plt.subplots(4,7, figsize=(18, 10), constrained_layout=True)
    num_subplots = 26  # Adjust to create 26 subplots
    # Flatten the axs array to iterate through subplots
    axs = axs.flatten()

    for i, file_id in enumerate(file_ids):
        
        data, header, wcs = get_data(file_id)
        ax = axs[i]
        # Set the WCS for this subplot
        ax.wcs = wcs
        
        xlim, ylim, beam = get_beam(header)
        levels = get_contour_levels(rms[i], sn[i], step[i])
        kpc_line = get_kpc_line(xlim, ylim, z[i])   
        contour_color = get_contour_color(spec_z[i])
        
            
        ax.imshow(data, origin = 'lower', cmap = 'Greys_r', interpolation = 'None', vmin = 0, vmax = 0.0001)
        ax.contour(data, levels = levels, colors=contour_color, linewidths = 1, zorder = 1, alpha = 1)
        ax.add_patch(beam) 
        ax.add_patch(kpc_line)

        if spec_z[i] == 0:
            z_text = "$z_p$="+str(round(z[i],2))
        else:
            z_text = "$z_s$="+str(round(z[i],2))

        ax.text(0.98,0.98, z_text, horizontalalignment='right', verticalalignment='top',transform=ax.transAxes, color='white', alpha=1)

        if step[i] > 1:
            text = str(step[i]) + "$\sigma$"
            ax.text(1, 0.01, text, horizontalalignment='right',
                    verticalalignment='bottom', transform=ax.transAxes, color='white', alpha=1)
        
        ax.set_title(f'BCG-{i + 1}')

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([])
        ax.set_yticks([])

    axs[0].text(0.1, 0.98, "10 kpc", horizontalalignment='left', verticalalignment='top',
            transform=axs[0].transAxes, color = 'white', alpha = 1)

    axs[26].set_visible(False)
    axs[27].set_visible(False)

    fig.suptitle("Postage stamps of our sample BCGs at 230GHz", size=16)
        
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.03,
                        hspace=0.05)  

    plt.show()

    return fig

fig = make_postage_stamps()
fig.savefig("plots/alma_postage_stamps.png", dpi = 400,bbox_inches = 'tight')



