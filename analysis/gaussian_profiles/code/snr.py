import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import glob

def read_fits(image, img_type = None):
    hdu = fits.open(image)[0]
    if img_type == 'alma':
        header = hdu.header
        data = hdu.data[0,0,:,:]
        wcs = WCS(header, naxis = 2)
    else:
        header = hdu.header
        data = hdu.data
        wcs = WCS(header)
    data = np.nan_to_num(data)
#    data /= np.amax(data)
    return header, data, wcs

def snr_rms(data, cent_coord, r, i_r, o_r):
    x = cent_coord[0]
    y = cent_coord[1]
    imin = x - o_r
    imax = x + o_r + 1
    jmin = y - o_r
    jmax = y + o_r + 1
    annulus = []
    circle = []
    for i in np.arange(imin, imax).astype(int):
        for j in np.arange(jmin, jmax).astype(int):
            ij = np.array([i,j])
            dist = np.linalg.norm(ij - np.array((x,y)))
            if dist > i_r and dist <= o_r:
                annulus.append(data[i][j])
            if dist <= r:
                circle.append(data[i][j])
    annulus = np.array(annulus)
    circle = np.array(circle)
    signal = np.max(circle)
    rms = np.sqrt(np.mean(annulus**2))
    snr = signal/rms
    return snr, rms
