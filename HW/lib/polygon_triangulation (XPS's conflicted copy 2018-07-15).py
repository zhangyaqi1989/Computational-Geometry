#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
implements polygon triangulation (ear cut algorithm)
"""

# standard library
import sys
from collections import namedtuple
##################################

# 3rd party library
import numpy as np
import matplotlib.pyplot as plt
from llist import dllist, dllistnode

# local module
from algorithm import is_between, is_left_on, is_left, intersect


Point2D = namedtuple('Point2D', ['x', 'y'])


def diagonal(idx_a, idx_b, vllist):
    """returns true iff (idx_a, idx_b) is a proper internal diagonal"""
    return is_in_cone(idx_a, idx_b, vllist) and is_in_cone(idx_a, idx_b, vllist) \
        and diagonalie(idx_a, idx_b, vllist)


def ear_init(vllist):
    """initialize ear condition for all vertices"""
    ear_dic = {}
    n = len(vllist)
    for i in range(n):
        idx_next = (i + 1) % n
        idx_prev = i - 1
        ear_dic[vllist[i]] = diagonal(idx_prev, idx_next, vllist)
    return ear_dic


def polygon_triangulate(points):
    """ return n - 2 triangles store in the triangulation tri
        each triangle has the representation as [(x1, y1), (x2, y2), (x3, y3)]
        (x1, y1) <--> (x2, y2) is the diagonal
    """
    vllist = dllist([Point2D(x, y) for x, y in points])
    n = len(vllist)
    ear_dic = ear_init(vllist)
    tri = []
    while n > 3:
        earfound = False
        for i in range(n):
            if ear_dic[vllist[i]]:
                earfound = True
                idx_v2 = i
                idx_v3 = (idx_v2 + 1) % n
                idx_v4 = (idx_v3 + 1) % n
                idx_v1 = idx_v2 - 1
                idx_v0 = idx_v1 - 1
                p1, p2, p3 = vllist[idx_v1], vllist[idx_v2], vllist[idx_v3]
                tri.append(((p1.x, p1.y), (p3.x, p3.y), (p2.x, p2.y)))
                # update ear condition of diagonal endpoints
                ear_dic[vllist[idx_v1]] = diagonal(idx_v0, idx_v3, vllist)
                ear_dic[vllist[idx_v3]] = diagonal(idx_v1, idx_v4, vllist)
                # cut off the ear index_v2
                del vllist[idx_v2]
                n -= 1
                break
    return tri


def diagonalie(idx_a, idx_b, vllist):
    """returns true iff (idx_a, idx_b) is a proper iternal or external
       diagonal of vllist
    """
    a = vllist[idx_a]
    b = vllist[idx_b]
    n = len(vllist)
    for i in range(n):
        c = vllist[i]
        c1 = vllist[(i + 1) % n]
        # print(c, c1)
        if (not are_points_identical(a, c)) and\
           (not are_points_identical(a, c1)) and\
           (not are_points_identical(b, c)) and\
           (not are_points_identical(b, c1)) and\
           intersect(a, b, c, c1):
            return False
    return True


def is_in_cone(idx_a, idx_b, vllist):
    """returns true iff diagonal (idx_a, idx_b) is strictly internal to
       the polygon in sht neighborhood of the endpoint
    """
    n = len(vllist)
    idx_a0 = idx_a - 1
    idx_a1 = (idx_a + 1) % n
    a = vllist[idx_a]
    b = vllist[idx_b]
    a0 = vllist[idx_a0]
    a1 = vllist[idx_a1]
    # 1. if a is a convex vertex
    if is_left_on(a, a1, a0):
        return is_left(a, b, a0) and is_left(b, a, a1)

    # 2. else a is reflex
    return not (is_left_on(a, b, a1) and is_left_on(b, a, a0))


def are_points_identical(a, b, tol=1e-12):
    """returns true if point a and b are identical considering numerical error
    >>> are_points_identical(Point2D(1, 2), Point2D(1, 2))
    True
    >>> are_points_identical(Point2D(1, 2), Point2D(1 + 1e-13, 2))
    True
    >>> are_points_identical(Point2D(1, 2), Point2D(1 + 1e-11, 2))
    False
    """
    return np.isclose(a.x, b.x, rtol=0, atol=tol) and np.isclose(
        a.y, b.y, rtol=0, atol=tol)


def read_vertices(filename):
    """read vertices from file"""
    coords = np.loadtxt(filename)
    return coords


def plot_polygon(vllist):
    """plot polygon defined by vertices llist"""
    fig, ax = plt.subplots()
    points = [(point.x, point.y) for point in vllist]
    points.append(points[0])
    xs, ys = zip(*points)
    ax.plot(xs, ys, 'b-o')
    plt.show()


def plot_polygon_triangulation(vllist, tri):
    """plot polygon defined by vertices llist and its triangulation"""
    fig, ax = plt.subplots()
    points = [(point.x, point.y) for point in vllist]
    points.append(points[0])
    xs, ys = zip(*points)
    ax.plot(xs, ys, 'b-o')
    for (x1, y1), (x2, y2), (_, _) in tri:
        ax.plot([x1, x2], [y1, y2], 'g--')
    plt.show()


def main():
    """main function"""
    if len(sys.argv) != 2:
        print("Usage: >> python {} <points_file>".format(sys.argv[0]))
        sys.exit(1)
    points_file = sys.argv[1]
    coords = read_vertices(points_file)
    tri = polygon_triangulate(coords)
    plot_polygon_triangulation(dllist([Point2D(x, y) for x, y in coords]), tri)


def doctest():
    """run doctest"""
    import doctest
    doctest.testmod()


if __name__ == "__main__":
    main()
    # doctest()
