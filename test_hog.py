import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from PIL import Image
from skimage.feature import hog
from skimage import data, color, feature, exposure
import time

start=time.time()

filen='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2765_2017_06_02/'+'win12_02_03_5.0X_bright.tif'

#for fun...
#filen='/Users/wheatley/Documents/Photos/xmas/G0207507.jpg'

#### finest grids ####
#'win34_02_07_5.0X_bright.tif' #even without blur, sigma=3 is too high a threshold. 1 seems ok...2 is better (at least at straighter lines)

#### fine grids ####
#'win12_02_03_5.0X_bright.tif' # with blur, no edges found for either sigma. without blur, sigma=3 does good job of removing bridges
#'win12_02_02_5.0X_bright.tif' # how does it deal with obvious defects? not bad... the 'bubbles' underneath add some distortion though (quantify)

#### coarse grids ####
#'win11_02_08_5.0X_bright.tif' - seems to work well. Compare histograms with and without gaussian blur to see how necessary that is
#'win11_01_01_5.0X_bright.tif' # Gaussian blur doesn't really help. Algorithm misses a lot near the corners... 
imraw = Image.open(filen)
im=np.array(imraw)
im = exposure.equalize_hist(im)

#np.zeros((128, 128))
#im[32:-32, 32:-32] = 1

#im = ndi.rotate(im, 15, mode='constant')
#im = ndi.gaussian_filter(im, 4)


####### QUESTIONS #########
# How important is it to filter noise via Gaussian? the technique calls it essential...
# How to calculate appropriate Gaussian to use based on structure size and pixel size? eg. for the finest grids, the # of pixels in a slat is almost equal to 4 so a 4 pixel Gaussian smooths them out entirely! Maybe this is included in the Canny theory...

############################

#im += 0.2 * np.random.random(im.shape)

# Compute the Canny filter for two values of sigma
#edges1 = feature.canny(im)
edges2 = feature.canny(im, sigma=3)

imedge = edges2
rawimage=imraw

#raw image
#fd_raw, hog_image_raw = hog(rawimage, orientations=9, pixels_per_cell=(16, 16),
#                    cells_per_block=(20, 20), visualise=True)
#edge detected image
fd_edge, hog_image_edge = hog(imedge, orientations=36, pixels_per_cell=(6, 16),
                    cells_per_block=(20, 20), visualise=True)
print 'fn took %.2f seconds' % (time.time() - start)

fig2, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(8, 8),
                                    sharex=True, sharey=True)

#fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=2, ncols=2, figsize=(8, 8), sharex=True, sharey=True)
#print np.shape(fd_raw)

#ax1.axis('off')
#ax1.imshow(imraw, cmap=plt.cm.gray)
#ax1.set_title('Input image, raw')
#ax1.set_adjustable('box-forced')

# Rescale histogram for better display
#hog_image_rescaled_raw = exposure.rescale_intensity(hog_image_raw, in_range=(0, 0.02))

#ax2.axis('off')
#ax2.imshow(hog_image_rescaled_raw, cmap=plt.cm.gray)
#ax2.set_title('Histogram of Oriented Gradients')
#ax1.set_adjustable('box-forced')

ax3.axis('off')
ax3.imshow(imedge, cmap=plt.cm.gray)
ax3.set_title('Input image, edge detection')
ax3.set_adjustable('box-forced')

# Rescale histogram for better display
hog_image_rescaled_edge = exposure.rescale_intensity(hog_image_edge, in_range=(0, 0.02))

ax4.axis('off')
ax4.imshow(hog_image_rescaled_edge, cmap=plt.cm.gray)
ax4.set_title('Histogram of Oriented Gradients')
ax1.set_adjustable('box-forced')

fig2.show()

