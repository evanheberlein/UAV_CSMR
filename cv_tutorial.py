import cv2 as cv
import sys

img = cv.imread(cv.samples.findFile("C:\\Users\\evanh\\Box\\Cornell\\Fall_2021\\Carp_UAV\\0036SET\\002\\IMG_0400_6.tif"))
if img is None:
    sys.exit("Could not read image!")

cv.imshow("Display window", img)
k = cv.waitKey(0)

if k == ord("s"):
    cv.imwrite("IR_400.tif", img)

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

tif_file = mpimg.imread('C:\\Users\\evanh\\Box\\Cornell\\Fall_2021\\Carp_UAV\\0036SET\\002\\IMG_0400_6.tif')
imarray = np.array(tif_file)
imgplot = plt.imshow(tif_file, cmap = 'jet')
