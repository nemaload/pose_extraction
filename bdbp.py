#!/usr/bin/env python
# Extract boundary information from c. elegans image.
# Usage: bdbp-boundary.py HDF5FILE FRAMENUMBER OUTPUTDIR
# Output: A "here's the worm" mask.
# XXX: work in progress

import math
import random

import numpy
import numpy.ma as ma
import scipy.ndimage as ndimage
import scipy.ndimage.morphology
import hdf5lflib

import networkx as nx

import cv
import cv2
import matplotlib.pyplot as plt
import matplotlib.patches
from matplotlib.path import Path

#various file processing/OS things
import os
import sys
import tables


NUM_SAMPLES = 80


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

def display_graph(ax, graph, points):
    verts = []
    codes = []
    for i,j in graph.edges():
        verts.append([points[i][1], points[i][0]])
        codes.append(Path.MOVETO)
        verts.append([points[j][1], points[j][0]])
        codes.append(Path.LINETO)
    path = Path(verts, codes)
    patch = matplotlib.patches.PathPatch(path, facecolor='none', edgecolor='blue', lw=1)
    ax.add_patch(patch)

def display_path(ax, pathlist, points):
    verts = []
    codes = []
    for i in pathlist:
        verts.append([points[i][1], points[i][0]])
        codes.append(Path.LINETO)
    codes[0] = Path.MOVETO
    path = Path(verts, codes)
    patch = matplotlib.patches.PathPatch(path, facecolor='none', edgecolor='green', lw=1)
    ax.add_patch(patch)


def computeEdgeDistances(uvframe):
    """
    Create a 2D matrix @edgedists as a companion to @uvframe,
    containing for each pixel a distance to the nearest edge (more precisely,
    the nearest 0-valued pixel).

    We compute @edgedists in a floodfill fashion spreading from zero-areas
    to the middle of one-areas iteratively, with distances approximated
    on the pixel grid.

    We return a tuple (edgedists, edgedirs), where edgedirs contains information
    about the relative offset of the nearest edge piece.
    """
    # edgedists is a masked array, with only already computed values unmasked;
    # at first, uvframe == 0 already are computed (as zeros)
    edgedists = ma.array(numpy.zeros(uvframe.shape, dtype = numpy.float), mask = (uvframe > 0))
    edgedirs = ma.array(numpy.zeros(uvframe.shape, dtype = (numpy.float, 2)), mask = [[[j,j] for j in i] for i in uvframe > 0])
    #numpy.set_printoptions(threshold=numpy.nan)
    #print edgedists
    #print edgedirs

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
            nearestnei = ma.argmin(neighbor_val)

            # We assert that this update never affects value other fields
            # visited later in this iteration of floodfill
            edgedists[tuple(i)] = neighbor_val[nearestnei]

            nearestneicoord = i + neighbor_ofs[nearestnei]
            #print "-", nearestneicoord, edgedirs[tuple(nearestneicoord)]
            edgedirs[tuple(i)] = edgedirs[tuple(nearestneicoord)] + tuple(neighbor_ofs[nearestnei])
            #print "+", i, edgedirs[tuple(i)]

    return (edgedists.data, edgedirs.data)

def sampleRandomPoint(uvframe):
    """
    Return a coordinate tuple of a random point with non-zero value in uvframe.
    """
    while True:
        c = (random.randint(0, uvframe.shape[0]-1), random.randint(0, uvframe.shape[1]-1))
        if uvframe[c] > 0:
            return c

def pointsToBackbone(points, uvframe):
    # Generate a complete graph over these points,
    # weighted by Euclidean distances
    g = nx.Graph()
    g.add_nodes_from(range(len(points)))
    for i in range(len(points)):
        for j in range(i+1, len(points)):
            # TODO: scipy's cpair? but we will need to construct
            # a graph anyway
            g.add_edge(i, j, {'weight': math.pow(points[i][0]-points[j][0], 2) + math.pow(points[i][1]-points[j][1], 2)})

    # Reduce the complete graph to MST
    gmst = nx.minimum_spanning_tree(g)

    # Show the MST
    # f = plt.figure()
    # imgplot = plt.imshow(uvframe, cmap=plt.cm.gray)
    # display_graph(f.add_subplot(111), gmst, points)
    # plt.show()

    # Diameter of the minimum spanning tree will generate
    # a "likely pose walk" through the graph
    tip0 = max(nx.single_source_dijkstra_path_length(gmst, 0).items(), key=lambda x:x[1])[0] # funky argmax
    (tip1_lengths, tip1_paths) = nx.single_source_dijkstra(gmst, tip0)
    tip1 = max(tip1_lengths.items(), key=lambda x:x[1])[0]
    backbone = tip1_paths[tip1]

    return backbone

def gradientAscent(edgedists, edgedirs, point):
    """
    We want to move the point along the gradient from the edge of the worm
    to the center. However, simple non-guided gradient ascend will obviously
    make all the points converge in some middle point; we do not want to
    move along the A-P axis. Therefore, we instead move _from_ the nearest
    edge.
    """
    #print edgedists[tuple(point)], edgedirs[tuple(point)], max(abs(edgedirs[tuple(point)]))
    walkDir = edgedirs[tuple(point)] / max(abs(edgedirs[tuple(point)]))
    #print point, walkDir
    bestDist = edgedists[tuple(point)]
    bestPoint = point
    # From now on, point may be a non-integer; however we always return an int
    while point > [0,0] and point < edgedists.shape:
        intpoint = [round(point[0]), round(point[1])]
        if edgedists[tuple(intpoint)] < bestDist:
            break
        bestDist = edgedists[tuple(intpoint)]
        bestPoint = intpoint
        point = [point[0] - walkDir[0], point[1] - walkDir[1]]
        #print ">", bestPoint, bestDist, point, edgedists[round(point[0]), round(point[1])]
    return bestPoint

def poseExtract(uvframe, edgedists, edgedirs):
    """
    Output a sequence of coordinates of pose curve control points.
    """
    # Pick a random sample of points
    points = [sampleRandomPoint(uvframe) for i in range(NUM_SAMPLES)]

    # Generate a backbone from the points set
    backbone = pointsToBackbone(points, uvframe)
    print backbone

    # Show the backbone
    f = plt.figure()
    imgplot = plt.imshow(uvframe, cmap=plt.cm.gray)
    display_path(f.add_subplot(111), backbone, points)
    plt.show()

    # Refine points on backbone by fixed-direction gradient ascend
    # over edgedists
    for i in backbone:
        #print "---", i, points[i]
        points[i] = gradientAscent(edgedists, edgedirs, points[i])
        #print "->", points[i]

    # Show the backbone
    f = plt.figure()
    imgplot = plt.imshow(uvframe, cmap=plt.cm.gray)
    display_path(f.add_subplot(111), backbone, points)
    plt.show()

    # TODO: Redo the complete graph - MST - diameter with final graph
    # to get fine tracing

    # TODO: Extend tips by slowest-rate gradient descent
    return backbone

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
    (edgedists, edgedirs) = computeEdgeDistances(uvframe)

    fig, axes = plt.subplots(ncols = 2)
    axes[0].imshow(uvframe, cmap=plt.cm.gray)
    axes[1].imshow(edgedists)
    plt.show()

    print poseExtract(uvframe, edgedists, edgedirs)

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
