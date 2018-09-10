#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
implements arbitrary triangulation of a set of 2D points
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from convex_hull import convex_hull

def plot_polygon(points):
    xs, ys = zip(*points)
    plt.plot(xs, ys, 'b-o')
    plt.plot([xs[-1], xs[0]], [ys[-1], ys[0]], 'b-o')


def _print_is_lefts(triangle, point):
    """for debugging purpose"""
    for i in range(3):
        j = (i + 1) % 3
        print(_is_left(triangle[i], triangle[j], point))


def _is_left(A, B, point, tol=1e-16):
    x0, y0 = A
    x1, y1 = B
    x, y = point
    v = (x1 - x0) * (y - y0) - (y1 - y0) * (x - x0)
    if abs(v) < tol:
        return 0
    elif v > 0:
        return 1
    else:
        return -1


def is_ccw(points):
    """check if points are in ccw"""
    n = len(points)
    area = 0
    xs, ys = zip(*points)
    for i in range(n):
        j = (i + 1) % n
        area += xs[i] * ys[j] - xs[j] * ys[i]
    return area > 0


def is_in(tri, point):
    """check if point is in a triangle or on an edge"""
    n = len(tri)
    s = set()
    for i in range(n):
        j = (i + 1) % n
        pos = _is_left(tri[i], tri[j], point)
        if pos == 0:
            return 0, i
        s.add(pos)
    if len(s) == 1:
        return 1, None
    return -1, None


def in_triangle(triangle, point):
    """test if a point is in triangle or not"""
    count = 0
    for i in range(3):
        j = (i + 1) % 3
        isl = _is_left(triangle[i], triangle[j], point)
        if isl > 0:
            count = count + 1
        elif isl < 0:
            return -1, -1
        else:
            v = 3 - (i + j)
    if count == 3:
        return 1, -1
    elif count == 2:
        return 0, v
    else:
        return -1, -1


def triangulation(in_points):
    if type(in_points) == np.ndarray:
        points = list(zip(in_points[:, 0], in_points[:, 1]))
    elif type(in_points) == list:
        points = in_points[:]
    else:
        print("unsupported points type {:s}".format(type(in_points)))
        sys.exit(1)
    hull = convex_hull(points, alg=2)
    triangles = []
    n = len(hull)
    for i in range(n-2):
        if _is_left(hull[0], hull[i+1], hull[i+2]) == 0:
            inner_points.append(hull[i + 1])
        triangle = (hull[0], hull[i+1], hull[i+2])
        if not is_ccw(triangle):
            triangle = triangle[::-1]
        triangles.append(triangle)
    inter_points = set(points) - set(hull)
    # print(inter_points)
    # triangles = set(triangles)
    # print(triangles)
    for point in inter_points:
        ones = 0
        tmp = []
        rms = []
        for triangle in triangles:
            p, v = in_triangle(triangle, point)
            if p == 1:
                rms.append(triangle)
                for j in range(3):
                    small_t = [triangle[j], triangle[(j+1)%3], point]
                    if not is_ccw(small_t):
                        small_t = small_t[::-1]
                    tmp.append(small_t)
                break
            elif p == 0:
                ones = ones + 1
                rms.append(triangle)
                vertex = triangle[v]
                # del triangle[v]
                for j in range(3):
                    if j == v:
                        continue
                    small_t = [triangle[j], vertex, point]
                    if not is_ccw(small_t):
                        small_t = small_t[::-1]
                    tmp.append(small_t)
                if ones == 2:
                    break
        triangles.extend(tmp)
        # print(tmp)
        for triangle in rms:
            triangles.remove(triangle)
    return triangles


def triangulation_test(n, lim_min=0, lim_max=5):
    points = np.random.uniform(lim_min, lim_max, (n, 2))
    # points = np.loadtxt("debug.txt")
    points = np.loadtxt("point-files/testPoints_heart.txt", delimiter=",")
    # points = np.loadtxt("../hw3-1/points-4.txt")
    triangles = triangulation(points)
    for triangle in triangles:
        plot_polygon(triangle)
    plt.show()

if __name__ == "__main__":
    n = 20
    triangulation_test(n)
