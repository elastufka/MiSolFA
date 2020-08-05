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
from matplotlib import cm
from skimage import feature, exposure
from skimage.transform import  (hough_line, hough_line_peaks,
                               probabilistic_hough_line)



# Generate noisy image of a square
filen='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2765_2017_06_02/'+'win12_02_02_5.0X_bright.tif'

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
imraw = Image.open(filen) #.rotate(45)#.convert('1').convert('L')
im=np.array(imraw)
#im = exposure.rescale_intensity(im, in_range=(0, 0.02))

#max contrast
im = exposure.equalize_hist(im)

# Equalization
#selem = disk(30)
#img_eq = rank.equalize(img, selem=selem)

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

#h1, theta1, d1 = hough_line(edges1)
#h2, theta2, d2 = hough_line(edges2)

lines1 = probabilistic_hough_line(edges1, threshold=10, line_length=50,
                                 line_gap=2)
lines2 = probabilistic_hough_line(edges2, threshold=10, line_length=50,
                                 line_gap=2)


fig, axes = plt.subplots(1, 3, figsize=(15, 6),
                         subplot_kw={'adjustable': 'box-forced'}, sharex=True, sharey=True)
ax = axes.ravel()

ax[0].imshow(imraw, cmap=cm.gray)
ax[0].imshow(np.ma.masked_where(edges2 == 0,edges2),cmap=cm.autumn)
ax[0].set_title('Input image')
ax[0].set_axis_off()

#ax[1].imshow(edges1, cmap=cm.gray)
for line in lines1:
    p0, p1 = line
    ax[1].plot((p0[0], p1[0]), (p0[1], p1[1]), color='r')
ax[1].set_xlim((0, np.shape(imraw)[0]))
ax[1].set_ylim((np.shape(imraw)[1], 0))
ax[1].set_title('Probabilistic Hough')

#ax[2].imshow(edges2, cmap=cm.gray)
for line in lines2:
    p0, p1 = line
    ax[2].plot((p0[0], p1[0]), (p0[1], p1[1]),color='r')
ax[2].set_xlim((0, np.shape(imraw)[0]))
ax[2].set_ylim((np.shape(imraw)[1], 0))
ax[2].set_title('Probabilistic Hough')

#ax[1].imshow(edges1, cmap=cm.gray)
#for _, angle, dist in zip(*hough_line_peaks(h1, theta1, d1)):
#    y0 = (dist - 0 * np.cos(angle)) / np.sin(angle)
#    y1 = (dist - np.shape(imraw)[1] * np.cos(angle)) / np.sin(angle)
#    ax[1].plot((0, np.shape(imraw)[1]), (y0, y1), '-r')
#ax[1].set_xlim((0, np.shape(imraw)[1]))
#ax[1].set_ylim((np.shape(imraw)[0], 0))
#ax[1].set_axis_off()
#ax[1].set_title('Detected lines, 1 sigma')

#ax[2].imshow(edges2, cmap=cm.gray)
#for _, angle, dist in zip(*hough_line_peaks(h2, theta2, d2)):
#    y0 = (dist - 0 * np.cos(angle)) / np.sin(angle)
#    y1 = (dist - np.shape(imraw)[1] * np.cos(angle)) / np.sin(angle)
#    ax[2].plot((0, np.shape(imraw)[1]), (y0, y1), '-r')
#ax[2].set_xlim((0, np.shape(imraw)[1]))
#ax[2].set_ylim((np.shape(imraw)[0], 0))
#ax[2].set_axis_off()
#ax[2].set_title('Detected lines, 3 sigma')

for a in ax:
    a.set_axis_off()
    a.set_adjustable('box-forced')


fig.tight_layout()
fig.show()

#calculate angles, make histogram
def get_length(line):
    deltax=line[1][0]-line[0][0]
    deltay=line[1][1]-line[0][1]
    length=(deltax**2 + deltay**2)**0.5
    return length

def get_angle(line):
    deltax=line[1][0]-line[0][0]
    deltay=line[1][1]-line[0][1]
    theta=np.arctan(deltay/deltax)
    thetadeg=theta*180./np.pi
    return thetadeg

nom = -45. #nominal angle

theta1,theta2=[],[]
for l1,l2 in zip(lines1,lines2):
    if get_length(l1) >= 50. :
        theta1.append(get_angle(l1))
    if get_length(l2) >= 50. :        
        theta2.append(get_angle(l2))
    
#make a histogram

fig,ax=plt.subplots()

foo=ax.hist(theta1,np.linspace(-90,90,90), facecolor='b',alpha=.6)
foo=ax.hist(theta2,np.linspace(-90,90,90), facecolor='g',alpha=.6)

ax.set_xlim([-90,90])
ax.set_yscale('log')
ax.set_ylim([1,1000])
fig.show()

#plot outilers
out1,out2=[],[]
for l1,l2 in zip(lines1,lines2):
    if get_angle(l1) !=-45.0:
       #theta1.append(get_angle(l1))
       out1.append(l1)
    if get_angle(l2) !=-45.0:
       #theta2.append(get_angle(l2))
       out2.append(l2)

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(5,6), sharex=True,sharey=True)
ax1.imshow(imraw, cmap=cm.gray)
ax2.imshow(imraw, cmap=cm.gray)
ax1.imshow(np.ma.masked_where(edges2 == 0,edges2),cmap=cm.cool,alpha=.7)
ax2.imshow(np.ma.masked_where(edges2 == 0,edges2),cmap=cm.cool,alpha=.7)

for line in out1:
    p0, p1 = line
    ax1.plot((p0[0], p1[0]), (p0[1], p1[1]), color='r')
ax1.set_xlim((0, np.shape(imraw)[0]))
ax1.set_ylim((np.shape(imraw)[1], 0))
ax1.set_title('Probabilistic Hough')

for line in out2:
    p0, p1 = line
    ax2.plot((p0[0], p1[0]), (p0[1], p1[1]), color='r')
ax2.set_xlim((0, np.shape(imraw)[0]))
ax2.set_ylim((np.shape(imraw)[1], 0))
ax2.set_title('Probabilistic Hough')

ax1.set_axis_off()
ax1.set_adjustable('box-forced')
ax2.set_axis_off()
ax2.set_adjustable('box-forced')

fig.show()
