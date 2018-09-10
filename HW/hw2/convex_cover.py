#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
compute convex hull, bounding box, diameter and width
of a set of 2D points
"""

import sys
sys.path.insert(0, '../lib/')
import numpy as np
import matplotlib.pyplot as plt
from convex_hull import convex_hull, bounding_rectangle, compute_polygon_area,\
        compute_diameter, compute_width


def convex_cover():
    '''
    read a set of points from file and compute
    convex hull, bounding rectangle, diameter and width of the points
    and visualize
    '''
    if len(sys.argv) != 2:
        print("Usage: python {} <filename>".format(sys.argv[0]))
        sys.exit(1)
    points = np.loadtxt(sys.argv[1])
    hull = convex_hull(points, 1)
    print("convex hull area: {:0.6f}".format(compute_polygon_area(hull)))
    rectangle, area = bounding_rectangle(hull)
    print("bounding rectangle area: {:0.6f}".format(area))
    D, pair = compute_diameter(points)
    (x0, y0), (x1, y1) = pair
    print("diameter: {:0.6f}".format(D))
    width, _, _, ends = compute_width(points)
    print("width = {:0.6f}".format(width))
    visualize(points, hull, rectangle, pair, ends)


def visualize(points, hull, rectangle, D_pair, W_ends):
    # 1. scatter points
    plt.scatter(points[:, 0], points[:, 1], color='g')
    # 2. plot convec hull
    M = len(hull)
    for i in range(M):
        j = (i + 1) % M
        plt.plot([hull[i][0], hull[j][0]], [hull[i][1], hull[j][1]], 'r--')
    # 3. plot bounding rectangle
    for i in range(4):
        j = (i + 1) % 4
        plt.plot([rectangle[i][0], rectangle[j][0]], [rectangle[i][1], rectangle[j][1]], 'b--')
    # 4. plot diameter
    (x0, y0), (x1, y1) = D_pair
    plt.plot([x0, x1], [y0, y1], 'k--')
    plt.scatter([x0, x1], [y0, y1], color="r")
    # 5. plot width
    plt.plot(W_ends[:, 0], W_ends[:, 1], 'g--')
    plt.scatter(W_ends[:, 0], W_ends[:, 1], color="r")
    plt.axis("equal")
    plt.show()


if __name__ == '__main__':
    convex_cover()
