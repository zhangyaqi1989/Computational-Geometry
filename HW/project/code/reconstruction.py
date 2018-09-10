#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
implement marching squares algorithm to reconstruct
the level set of an implicit function
"""

# standard library
import sys
import re

# third party library
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

# local modules
sys.path.insert(0, '../../lib/')
from marching_square import compute_grids, make_function,\
            marching_squares, construct_paths, is_simple_polygon


def user_input():
    """ handle user input """
    if len(sys.argv) != 2:
        print("Usage: python {} <expression_file>".format(sys.argv[0]))
        sys.exit(1)
    with open(sys.argv[1], 'r') as in_file:
        lines = [line.strip() for line in in_file.readlines() if line.strip()]
    expression = lines[0]
    limits = [float(token) for token in re.split(r',\s*', lines[1])]
    h = float(lines[2])
    return expression, limits, h


def main():
    expression, limits, h = user_input()
    # print(expression, limits, h)
    nx, ny = compute_grids(limits, h)
    fun = make_function(expression)
    segs = marching_squares(fun, limits, nx=nx, ny=ny, plot_points=False)
    paths = construct_paths(segs)
    fig, ax = plt.subplots(figsize=(8, 8))
    patches = []
    for path in paths:
        patch = Polygon(path, True)
        patches.append(patch)
    patch_collection = PatchCollection(patches)
    ax.add_collection(patch_collection)
    if len(paths) == 0:
        print("there is no closed level set found in the given domain")
        sys.exit(1)
    if len(paths) > 1:
        print("The constructed boundary is not a simple polygon")
        plt.show() # draw the paths
        sys.exit(1)
    polygon = paths[0]
    # check polygon is a simple polygon
    if not is_simple_polygon(polygon):
        print("The constructed boundary is not a simple polygon")
        plt.show() # draw the nonsimple polygon
        sys.exit(1)
    # np.savetxt('polygon-1.txt', polygon, delimiter='\t')

    # plot the constructed boundary
    closed_polygon = np.vstack((polygon, polygon[0:1, :]))
    ax.plot(closed_polygon[:, 0], closed_polygon[:, 1], 'b-o')
    ax.axis('equal')
    plt.show()


if __name__ == "__main__":
    main()
