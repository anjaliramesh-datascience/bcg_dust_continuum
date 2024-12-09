vis = "/Users/arames52/Research/Data/MS/ES1_z_1.40_calibrated.ms"
output_directory = "/Users/arames52/Research/Data/"
imsize = [864,864]
field = 'ES1_z_1.4'
cell = '0.045arcsec'
gridder = 'standard'
deconvolver = 'hogbom'
specmode = 'mfs'
niter = 5000
pbcor = True
usemask = 'auto-multithresh'
threshold = "0.03mJy"
interactive = True
weighting = 'natural'
robust = 0.5
sidelobethreshold=3.0
noisethreshold=5.0
lownoisethreshold=1.5
negativethreshold=7.0
minbeamfrac=0.3
intent = 'OBSERVE_TARGET#ON_SOURCE'
imagename = output_directory + "ES1_z_1.40_natural"
tclean(vis = vis, imagename = imagename, imsize = imsize, intent = intent, cell = cell,
specmode = specmode, gridder = gridder, deconvolver = deconvolver,usemask = usemask,
weighting = weighting, niter = niter, field = field,
pbcor = pbcor, robust = robust, threshold = threshold, interactive = interactive,
sidelobethreshold=sidelobethreshold, noisethreshold=noisethreshold, lownoisethreshold=lownoisethreshold,
negativethreshold=negativethreshold,minbeamfrac=minbeamfrac)

usemask = "user"
mask = output_directory + "ES1_z_1.40_natural.mask"
weight = ["briggs", "briggs", "natural", "natural", "natural", "natural"]
taper = ["", "", "0.7arcsec", "1arcsec", "2arcsec", "3arcsec"]
r = [0, 0.5, 0.5, 0.5, 0.5, 0.5]
file_postfix = ["briggs_0", "briggs_0.5", "taper_0.7", "taper_1", "taper_2", "taper_3"]

for i in range(len(weight)):
    imagename = output_directory + "ES1_z_1.40_" + file_postfix[i]
    tclean(vis = vis, imagename = imagename, imsize = imsize, intent = intent, cell = cell,
    specmode = specmode, gridder = gridder, deconvolver = deconvolver,usemask = usemask, mask = mask,
    weighting = weight[i], niter = niter, field = field,
    pbcor = pbcor, robust = r[i], threshold = threshold, interactive = False, uvtaper = taper[i])
