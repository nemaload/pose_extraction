#!/usr/bin/python
#
# Usage: ./raw2png.py inputFile outputFile width height

import sys

import numpy
import scipy.misc

inputFile = sys.argv[1]
outputFile = sys.argv[2]
width = int(sys.argv[3])
height = int(sys.argv[4])

f = open(inputFile, 'rb')
frame = numpy.fromfile(f, dtype = 'short').reshape((height, width))
f.close()

scipy.misc.imsave(outputFile, frame)
