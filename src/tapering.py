import glob

imsize=[864,864]
cell='0.045arcsec'
specmode='mfs'
gridder='standard'
deconvolver='hogbom'
niter=5000
usemask='user'
field='0'
pbcor=True
threshold='0.03mJy'
interactive=False
robust = 0.5
weighting = 'natural'

# Change these for tapering with different on-sky FWHM
uvtaper = '6arcsec'
output = '6"/'

path = "/Users/arames52/Research/Data/MS/"
op_path = "/Users/arames52/Research/Data/CASA_Outputs/UVtaper/"
ms_files = sorted(glob.glob(path + "*.ms"))

for file in ms_files:
    vis = file
    name = file.split(path)[-1].split("_calibrated.ms")[0]
    mask = path + "Masks/" + name + "_natural.mask"
    imagename = op_path + output + name
    tclean(vis = vis, imagename = imagename, imsize = imsize, cell = cell,
    specmode = specmode, gridder = gridder, deconvolver = deconvolver,
    weighting = weighting, niter = niter, usemask = usemask, mask = mask, field = field,
    pbcor = pbcor, robust = robust, threshold = threshold, interactive = interactive,
    uvtaper = uvtaper)
    exportfits(imagename = imagename + ".image.pbcor", fitsimage = imagename + ".fits")
