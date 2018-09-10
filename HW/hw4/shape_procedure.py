#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
implement shape procedure for shape reconstruction
"""

import sys
sys.path.insert(0, '../lib/')
import numpy as np
import matplotlib.pyplot as plt
from triangulation_DS import Triangulation


def shape():
    """apply crust algorithm to a set of points defined in file"""
    if len(sys.argv) != 2:
        print("Usage: python {} <filename>".format(sys.argv[0]))
        sys.exit(1)
    points = np.loadtxt(sys.argv[1])
    tri = Triangulation(points)
    figure = plt.figure(figsize=(8, 8))
    plt.subplot(1, 1, 1)
    tri.shape_procedure()
    plt.show()


if __name__ == '__main__':
    shape()
