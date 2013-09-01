#!/usr/bin/env python
# Extract boundary information from c. elegans image.
# Usage: bdbp-boundary.py HDF5FILE FRAMENUMBER OUTPUTDIR
# Output: A "here's the worm" mask.
# XXX: work in progress

import math
import numpy
import numpy.ma as ma
import scipy.ndimage as ndimage
import scipy.ndimage.morphology
import hdf5lflib

import cv
import cv2
import matplotlib.pyplot as plt
import matplotlib.patches

#various file processing/OS things
import os
import sys
import tables

def print_mask(mask):
    """
    A debug print function that prints out a bool array in asciiart.
    """
    for row in mask:
        for item in row:
            sys.stdout.write('#' if item else '.')
        sys.stdout.write('\n')

def edge_dist_if_within(edgedists, coord):
    """
    Return edge distance at a given coord (possibly ma.masked) or
    ma.masked if we are accessing outside of the image.
    """
    if coord[0] < 0 or coord[1] < 0:
        return ma.masked
    try:
        return edgedists[tuple(coord)]
    except IndexError:
        return ma.masked

def computeEdgeDistances(uvframe):
    """
    Create and return a 2D matrix @edgedists as a companion to @uvframe,
    containing for each pixel a distance to the nearest edge (more precisely,
    the nearest 0-valued pixel).

    We compute @edgedists in a floodfill fashion spreading from zero-areas
    to the middle of one-areas iteratively, with distances approximated
    on the pixel grid.
    """
    # edgedists is a masked array, with only already computed values unmasked;
    # at first, uvframe == 0 already are computed (as zeros)
    edgedists = ma.array(numpy.zeros(uvframe.shape, dtype = numpy.float), mask = (uvframe > 0))
    #numpy.set_printoptions(threshold=numpy.nan)
    #print edgedists

    flood_spread = scipy.ndimage.morphology.generate_binary_structure(2, 2)
    neighbor_ofs = [[-1,-1],[-1,0],[-1,1], [0,-1],[0,0],[0,1],  [1,-1],[1,0],[1,1]]
    s2 = math.sqrt(2)
    neighbor_dist = [s2,1,s2, 1,0,1, s2,1,s2]

    while ma.getmaskarray(edgedists).any():
        # scan masked area for any elements that have unmasked neighbors
        done_mask = numpy.invert(ma.getmaskarray(edgedists))
        todo_mask = done_mask ^ scipy.ndimage.binary_dilation(done_mask, flood_spread)
        #print_mask(todo_mask)
        for i in numpy.transpose(numpy.nonzero(todo_mask)):
            neighbor_val = ma.array([
                    edge_dist_if_within(edgedists, i + ofs) + dist
                        for ofs, dist in zip(neighbor_ofs, neighbor_dist)
                ])
            newdist = ma.min(neighbor_val)
            # We assert that this update never affects value other fields
            # visited later in this iteration of floodfill
            edgedists[tuple(i)] = newdist

    return edgedists.data

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

    # Annotate with information regarding the nearest edge
    edgedists = computeEdgeDistances(uvframe)

    fig, axes = plt.subplots(ncols = 2)
    axes[0].imshow(uvframe, cmap=plt.cm.gray)
    axes[1].imshow(edgedists)
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
