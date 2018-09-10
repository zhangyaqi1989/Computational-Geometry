#!/usr/bin/env python3
"""
written by Yaqi Zhang (zhang623@wisc.edu)
University of Wisconsin-Madison

"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull


def hw2():
    if len(sys.argv) != 2:
        print("Usage: python {} <filename>".format(sys.argv[0]))
        sys.exit(1)
    points = np.loadtxt(sys.argv[1])
    hull = ConvexHull(points)
    print("convex hull area = {:0.6f}".format(hull.volume))

if __name__ == '__main__':
    hw2()
