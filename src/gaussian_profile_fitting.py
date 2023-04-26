from astropy.io import fits
import sys
sys.path.append('/Users/arames52/bcg_dust_continuum/src/')
import bcg_parameter_file as p
from data_ingestion import bcg_regions_load as rl
import glob

region_parameters = rl(p.region_parameter_file_path)
bcgs = list(region_parameters['id'])
ras = list(region_parameters['ra'])
decs = list(region_parameters['dec'])
radius = list(region_parameters['radius'])

def region_params(i):
    x, y = str(ras[i]) + "deg", str(decs[i]) + "deg"
    r = str(round(radius[i],2)) + "arcsec"
    region = "circle[["+ x +","+y+ "]," + r+ "]"
    return region

natural_imfit_results = {}

for i in range(len(bcgs)):
    bcg_image_file = glob.glob(p.alma_images_path + bcgs[i] + "*.fits")[0]
    bcg_region = region_parameters(i)
    natural_imfit_results[bcgs[i]] = imfit(bcg_image_file, region = bcg_region)

with open(p.data_path + "natural_imfit_results.pkl", "wb") as f:
    pickle.dump(natural_imfit_results, f)


