import glob
import matplotlib.pyplot as plt
import pickle
import numpy as np

path = "/Users/arames52/Desktop/CASA Imaging/MS_files/measurement_sets/"
# bcgs = ['CDFS-18', 'ES1-18', "ES1_z_0.88", "ES1-25"]
component_files = sorted(glob.glob(path + "*.cl"))
bcgs = ['CDFS-18','ES1-18','ES1-25','ES1_z_0.99','ES1_z_0.99b','ES1_z_1.04','ES1_z_1.38','ES1_z_1.60','ES1_z_1.65','ES1_z_1.70','XMM-113','XMM-11','XMM-29','XMM_z_0.81','XMM_z_0.9', 'XMM_z_1.0']
with open("imfit_results.pkl", "rb") as f:
    imfit_results = pickle.load(f)

imfit_flux = imfit_results['natural']

fig,ax = plt.subplots(1,1, figsize = (8,6))
for bcg in bcgs:
    cl.open(path + bcg + "_uv.cl")
    fit = cl.getcomponent(0)
    uv_flux = fit['flux']['value'][0]*1000
    uv_error = fit['flux']['error'][0]*1000
    image_flux = imfit_flux[bcg]['deconvolved']['component0']['flux']['value'][0]*1000
    image_error = imfit_flux[bcg]['deconvolved']['component0']['flux']['error'][0]*1000
    # fwhm = imfit_flux[bcg]['results']['component0']['beam']['beamarcsec']['major']['value']
    ax.scatter(uv_flux, image_flux, marker = '*', color = 'red', s = 10)
lims = [
np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
ax.plot(lims, lims, ':', color = 'black',alpha=0.75, zorder=0)

median_uv_error = np.median(uv_error)
median_image_error = np.mean(image_error)
# text1 = "{:.3f} median uv flux error \n".format(median_uv_error)
text2 = "{:.3f}mJy mean image flux error".format(median_image_error)

text = text1 + text2
ax.text(x = 0.5, y = 0.1, s = text2, transform = ax.transAxes)
    # ax.errorbar(uv_flux, image_flux, yerr = image_error, xerr = uv_error, ls = 'none', ecolor = 'black', alpha = 0.5,solid_capstyle='projecting', capsize=5, elinewidth = 1)
    # ax.scatter(uv_flux, uv_flux, marker = marker, color = 'blue', label = 'uv flux', s = 10)
    # ax.errorbar(fwhm, uv_flux, yerr = uv_error,ls = 'none', ecolor = 'black', alpha = 0.5,solid_capstyle='projecting', capsize=5, elinewidth = 1)
    # quad_err = np.sqrt(image_error**2 + uv_error**2)
    # print(quad_err)
    # fd = np.abs((uv_flux-image_flux))/quad_err
    # text = bcg + " " + str(round(fd, 2)) + "$\sigma$ (" + '{:0.2e}'.format(quad_err) + ")"
    # if bcg == 'ES1-18':
    #     xpos = fwhm-0.05
    #     ypos = image_flux-0.15
    # else:
    #     xpos = fwhm-0.05
    #     ypos = image_flux
    #
    # ax.text(xpos, ypos, s = text)

ax.set_xlabel("UV Flux [mJy]")
ax.set_ylabel("Image Flux [mJy]")
# handles, labels = plt.gca().get_legend_handles_labels()
# by_label = dict(zip(labels,handles))
# plt.legend(by_label.values(), by_label.keys())
plt.savefig("/Users/arames52/Desktop/CASA Imaging/Plots/uv_image_flux.jpg", dpi = 300)
