# # From http://tsaith.github.io/combine-images-into-a-video-with-python-3-and-opencv-3.html
# import cv2 as cv
# import argparse
# import os
# import glob
#
# # Constructing argument parser
# # ap = argparse.ArgumentParser()
# # ap.add_argument("-ext", "--extension", required=False, default='tif', help="extension name. default is 'tif'.")
# # ap.add_argument("-o", "--output", required=False, default='output.mp4', help="output video file")
# # args = vars(ap.parse_args())
#
# # Arguments to parser - need recursive filepath: use glob https://docs.python.org/3/library/glob.html
# # dir_path = glob.glob(r"C:/Users/evanh/Box/Cornell/Fall_2021/Carp_UAV/0036SET/**/*6.tif", recursive = True)
# paths = [os.path.normpath(i) for i in glob.glob(r"C:/Users/evanh/Box/Cornell/Fall_2021/Carp_UAV/0036SET/**/*6.tif",
#                                                 recursive = True)] # .normpath converts all / to \\ - still readable?
#
# # path = "C:/Users/evanh/Box/Cornell/Fall_2021/Carp_UAV/0036SET"
# # fname = []
# # for root, d_names, f_names in os.walk(path): # os.walk():https://www.bogotobogo.com/python/python_traversing_directory_tree_recursively_os_walk.php
# # 	for f in f_names:
# #         if fname.endswith('6.tif')
# # 		    fname.append(os.path.join(root, f))
# # ext = args['extension']
# output = 'carp_uav_ir.mp4'
# # images = []
# # for f in glob._listdir(dir_path): # glob:https://www.geeksforgeeks.org/how-to-use-glob-function-to-find-files-recursively-in-python/
# #     if f.endswith(ext):
# #         images.append(f)
#
# # Extract dimensions from first image
# image_path = paths[0]
# frame = cv.imread(image_path)
# cv.imshow('video', frame)
# height, width, channels = frame.shape
#
# # Set video codec
# fourcc = cv.VideoWriter_fourcc(*'mp4v')
# out = cv.VideoWriter(output, fourcc, 20.0, (width, height))
#
# for image in paths:
#
#     frame = cv.imread(image)
#
#     out.write(frame)
#
#     cv.imshow('video', frame)
#     if (cv.waitKey(1) & 0xFF) == ord('q'): # Hit q to exit
#         break
#
# # Release output
# out.release()
# cv2.destroyAllWindows()
#
# # print("The output video is {}".format(output))
#
# # Can wite .mp4 file, but get error '!_src.empty() in function 'cv::cvtColor'', tif has no default color?
# # https://docs.opencv.org/3.4/d8/d01/group__imgproc__color__conversions.html#ga397ae87e1288a81d2363b61574eb8cab cvtColor doc

# Final tutorial used:  https://stackoverflow.com/questions/34975972/how-can-i-make-a-video-from-array-of-images-in-matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation
import matplotlib.image as mpimg
import os
import glob
import numpy as np

# Initialize empty arrays/figures
frames = []
img = []
fig = plt.figure()

path = [os.path.normpath(i) for i in glob.glob(r"C:/Users/evanh/Box/Cornell/Fall_2021/Carp_UAV/0036SET/**/*6.tif",
                                                recursive = True)] # .normpath converts all / to \\ - still readable?for root, d_names, f_names in os.walk(path): # os.walk():https://www.bogotobogo.com/python/python_traversing_directory_tree_recursively_os_walk.php
path = path[0:990] # Set range here
for p in path:
    tif_file = mpimg.imread(p)
    img.append(np.array(tif_file))

for i in range(len(path)):
    frames.append([plt.imshow(img[i],
                              cmap = 'jet',
                              # cmap=cm.Greys_r,
                              animated = True)])
    plt.colorbar()
    plt.show()
    #plt.close('all') # Not working? Figs stacking up
    print(i)

ani = animation.ArtistAnimation(fig, frames, interval = 50, blit=True, repeat_delay = 1000)
ani.save('full_carp.mp4')

