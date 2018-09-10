#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
implement crust algorithm for shape reconstruction
"""

import sys
sys.path.insert(0, '../lib/')
import numpy as np
import matplotlib.pyplot as plt
from triangulation_DS import Triangulation


def crust():
    """apply crust algorithm to a set of points defined in file"""
    if len(sys.argv) != 2:
        print("Usage: python {} <filename>".format(sys.argv[0]))
        sys.exit(1)
    points = np.loadtxt(sys.argv[1])
    tri = Triangulation(points)
    figure = plt.figure(figsize=(8, 8))
    plt.subplot(1, 1, 1)
    tri.crust()
    plt.show()


if __name__ == '__main__':
    crust()
