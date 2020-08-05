"""
===================
Canny edge detector
===================

The Canny filter is a multi-stage edge detector. It uses a filter based on the
derivative of a Gaussian in order to compute the intensity of the gradients.The
Gaussian reduces the effect of noise present in the image. Then, potential
edges are thinned down to 1-pixel curves by removing non-maximum pixels of the
gradient magnitude. Finally, edge pixels are kept or removed using hysteresis
thresholding on the gradient magnitude.

The Canny has three adjustable parameters: the width of the Gaussian (the
noisier the image, the greater the width), and the low and high threshold for
the hysteresis thresholding.

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from PIL import Image
from skimage import feature, exposure


# Generate noisy image of a square
filen='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2765_2017_06_02/'+'win11_01_01_5.0X_bright.tif'

#for fun...
#filen='/Users/wheatley/Documents/Photos/xmas/G0207507.jpg'
#filen='/Users/wheatley/Documents/Photos/xmas/GOPR7982.jpg'
#filen='/Users/wheatley/Documents/Photos/IMG_20161227_055911795.jpg'


#### finest grids ####
#'win34_02_07_5.0X_bright.tif' #even without blur, sigma=3 is too high a threshold. 1 seems ok...2 is better (at least at straighter lines)

#### fine grids ####
#'win12_02_03_5.0X_bright.tif' # with blur, no edges found for either sigma. without blur, sigma=3 does good job of removing bridges
#'win12_02_02_5.0X_bright.tif' # how does it deal with obvious defects? not bad... the 'bubbles' underneath add some distortion though (quantify)

#### coarse grids ####
#'win11_02_08_5.0X_bright.tif' - seems to work well. Compare histograms with and without gaussian blur to see how necessary that is
#'win11_01_01_5.0X_bright.tif' # Gaussian blur doesn't really help. Algorithm misses a lot near the corners... 
imraw = Image.open(filen)#.convert('1').convert('L')
im=np.array(imraw)
#im = exposure.rescale_intensity(im, in_range=(0, 0.02))

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
edges1 = feature.canny(im)
edges2 = feature.canny(im, sigma=3)

# display results
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(3, 8),
                                    sharex=True, sharey=True)

ax1.imshow(im, cmap=plt.cm.gray)
ax1.axis('off')
ax1.set_title('noisy image', fontsize=14)

ax2.imshow(edges1, cmap=plt.cm.gray)
ax2.axis('off')
ax2.set_title('Canny filter, $\sigma=1$', fontsize=14)

ax3.imshow(edges2, cmap=plt.cm.gray)
ax3.axis('off')
ax3.set_title('Canny filter, $\sigma=3$', fontsize=14)

fig.tight_layout()

fig.show()

