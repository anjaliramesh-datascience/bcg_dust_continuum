imsize = [864,864]
field = 'ES1_z_1.4'
cell = '0.045arcsec'
gridder = 'standard'
deconvolver = 'hogbom'
specmode = 'mfs'
niter = 5000
pbcor = True
usemask = 'auto-multithresh'
threshold = "0.06mJy"
interactive = True
weighting = 'natural'
robust = 0.5
sidelobethreshold=3.0
noisethreshold=5.0
lownoisethreshold=1.5
negativethreshold=7.0
minbeamfrac=0.3
intent = 'OBSERVE_TARGET#ON_SOURCE'

tclean(vis = vis, imagename = imagename, imsize = imsize, intent = intent, cell = cell,
specmode = specmode, gridder = gridder, deconvolver = deconvolver,usemask = usemask,
weighting = weighting, niter = niter, field = field,
pbcor = pbcor, robust = robust, threshold = threshold, interactive = interactive,
sidelobethreshold=sidelobethreshold, noisethreshold=noisethreshold, lownoisethreshold=lownoisethreshold,
negativethreshold=negativethreshold,minbeamfrac=minbeamfrac)
