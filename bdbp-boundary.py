#!/usr/bin/env python
# Extract boundary information from c. elegans image.
# Usage: bdbp-boundary.py HDF5FILE FRAMENUMBER OUTPUTDIR
# Output: A "here's the worm" mask.
# XXX: work in progress

import math
import numpy
import hdf5lflib

import cv
import cv2
import matplotlib.pyplot as plt
import matplotlib.patches

#various file processing/OS things
import os
import sys
import tables

def processFrame(i, node, outputBase, ar, cw):
    uvframe = hdf5lflib.compute_uvframe(node, ar, cw)

    plt.figure()
    imgplot = plt.imshow(uvframe, cmap=plt.cm.gray)
    plt.show()

    # Smooth twice
    uvframe = cv2.medianBlur(uvframe, 3)
    uvframe = cv2.medianBlur(uvframe, 3)

    plt.figure()
    imgplot = plt.imshow(uvframe, cmap=plt.cm.gray)
    plt.show()

    # Threshold
    background_color = uvframe.mean()
    foreground_i = uvframe > background_color
    uvframe[foreground_i] = 255.
    uvframe[numpy.invert(foreground_i)] = 0.

    plt.figure()
    imgplot = plt.imshow(uvframe, cmap=plt.cm.gray)
    plt.show()

def processFile(filename, outputDirectoryPath, frameNo):
    h5file = tables.open_file(filename, mode = "r")
    ar = h5file.get_node('/', '/autorectification')
    cw = h5file.get_node('/', '/cropwindow')
    outputBase = outputDirectoryPath + os.path.splitext(os.path.basename(filename))[0]
    processFrame(frameNo, h5file.get_node('/', '/images/' + str(frameNo)), outputBase, ar, cw)
    return True

if __name__ == '__main__':
    filename = sys.argv[1]
    frameNo = int(sys.argv[2])
    outputDirectoryPath = sys.argv[3]
    if not processFile(filename, outputDirectoryPath, frameNo):
        sys.exit(1)
