import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import glob

data_path = "/Users/arames52/Research/Data/Images/ALMA/"

residual_path = "/Users/arames52/Research/Data/GnS_profile/Gaussian/Residuals/"

rms_image = {}
rms_residual = {}
snr = {}

images = sorted(glob.glob(data_path + "*.fits"))
residuals = sorted(glob.glob(residual_path + "*.fits"))

def ra_dec(file):
    hdu = fits.open(file)[0]
    header = hdu.header
    ra = header['CRVAL1']
    dec = header['CRVAL2']
    return ra,dec

for file1, file2 in zip(images, residuals):
    bcg_name = file1.split("_natural.fits")[0].split(data_path)[-1]
    # bcg_name = file.split("_residual.image.fits")[0].split(residual_path)[-1]
    # residual_file = residual_path + bcg_name + "_residual.image"
    r,d = ra_dec(file1)
    ra = str(r) + "deg"
    dec = str(d) + "deg"
    r1 = "3arcsec"
    r2 = "6arcsec"
    # r3 = "2arcsec"
    region1 = "annulus[[" + ra + "," + dec + "]," + "["+ r1 +","+r2+ "]]"
    region2 = "circle[[" + ra + "," + dec + "]," + r1 + "]"

    result1 = imstat(imagename = file1, region = region1)
    signal = imstat(imagename = file1, region = region2)
    result2 = imstat(imagename = file2, region = region2)

    rms = result1['rms'][0]
    rms_residual[bcg_name] = result2['rms'][0]
    max_signal = signal['max'][0]
    snr[bcg_name] = max_signal/rms
    rms_image[bcg_name] = rms

for key in rms_image.keys():
    print(key, rms_image[key], snr[key])

with open("/Users/arames52/Research/Data/rms.pkl", "wb") as f:
    pickle.dump(rms_image, f)

with open("/Users/arames52/Research/Data/snr.pkl", "wb") as f:
    pickle.dump(snr, f)
