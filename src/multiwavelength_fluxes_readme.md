# Web Scraping Multi-Wavelength BCG Fluxes

## Sources
1. DES query interface - https://datalab.noirlab.edu/query.php?name=des_dr2.main
    - Python module astro-datalab can be used to query the data
2. Spitzer SWIRE - https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd - for direct download
    - However, we can use Python's curl function to directly download from the catalog database
3. Herschel HerMES - https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd
    - urllib.request can be used to download give the catalog url (specified in params.py)
4. ALMA flux - Calculated from IMFIT

## Input files/parameters
1. Coordinate file of sources in the format name|ra|dec
    - Read as csv file
2. HerMES Tar file download link

# Code Description

## multiwavelength_catalog_download.py

 **Replace /Users/arames52/ by the directory in which you have stored the bcg_dust_continuum folder**

### download_des_catalog(ra, dec, fov)
- Given the ra, dec and fov (default value 0.005 degrees) of target object, we download the magnitudes of sources within 20 arcseconds of the bcg coordinate.
- The output is a dictionary with keys being bcg name and values being dataframe of object magnitudes
- This is stored as a pickle file

### download_swire_catalog(name, ra, dec, bcg)
- Given the catalog name, ra, dec and id of the target source, we download the flux of objects within 20 arcseconds of the target.
- The output file is stored as csv files for individual bcgs under data/catalogs/SWIRE_catalog/

### download_hermes_catalog(catalog_url)
- We download the DR3 xID24 catalog of all HerMES sources
- This is stored for each field (CDFS, ES1, XMM)


## multiwavelength_catalog_compilation.py

Reads the downloaded DES, SWIRE, HerMES and ALMA data and forms a combined multiwavelength flux catalog for all 26 BCGs


### catalog_matching(c, df)
- Takes in the bcg skycoord and catalog dataframe and returns the source and nearest neighbours
### mag_to_mJy()
- Converts DES magnitudes to mJy units
        *  mJy is the desired unit for SED fitting

# To Run the files - 

python multiwavelength_catalog_download.py
python muliwavelength_catalog_compilation.py
