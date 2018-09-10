#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
Define some utility functions
"""

import sys
import random
import math
from collections import defaultdict
import pdb
import numpy as np
import matplotlib.pyplot as plt
from convex_hull import graham_scan, jarvis, monotone_chain

def init_points(n, xmin=0, xmax=1, ymin=0, ymax=1):
    xs = np.random.uniform(xmin, xmax, n)
    ys = np.random.uniform(ymin, ymax, n)
    return (xs, ys)


def shuffle_polygon(points):
    xs, ys = zip(*points)
    lines = []
    n = len(xs)
    for i in range(n):
        j = (i + 1) % n
        lines.append(((xs[i], ys[i]), (xs[j], ys[j])))
    random.shuffle(lines)
    return lines

def hw2():
    n = 10
    x0, y0 = 3, 4
    r = 5
    ts = np.linspace(0, np.pi * 2 - 0.1, n)
    xs = x0 + np.cos(ts) * r
    ys = y0 + np.sin(ts) * r
    points = list(zip(list(xs), list(ys)))
    # points = [(0, 0), (1, 0), (1, 1), (0, 1)]
    n = len(points)
    lines = shuffle_polygon(points)
    (x0, y0), (x1, y1) = lines[0]
    polygon = [(x0, y0), (x1, y1)]
    used = [0] * n
    used[0] = 1
    # pdb.set_trace()
    for _ in range(n - 1):
        p0 = polygon[-1]
        for i, line in enumerate(lines):
            if used[i]:
                continue
            else:
                p1, p2 = line
                if is_point_close(p0, p1):
                    used[i] = 1
                    polygon.append(p2)
                    break
                elif is_point_close(p0, p2):
                    used[i] = 1
                    polygon.append(p1)
                    break
        else:
            print("fail")
    print(polygon)

def test_convex_hull():
    n_points = 30
    xs, ys = init_points(n_points)
    points = list(zip(xs, ys))
    plt.xlim([-0.2, 1.2])
    plt.ylim([-0.2, 1.2])
    plt.scatter(xs, ys, color='b')
    # hull = graham_scan(points)
    # hull = jarvis(points)
    hull = monotone_chain(points)
    M = len(hull)
    for i in range(M):
        j = (i + 1) % M
        plt.plot([hull[i][0], hull[j][0]], [hull[i][1], hull[j][1]], 'r--')
    plt.show()


def is_point_close(p0, p1, rel_tol=1e-6):
    x0, y0 = p0
    x1, y1 = p1
    return (math.isclose(x0, x1, rel_tol=rel_tol) and 
            math.isclose(y0, y1, rel_tol=rel_tol))

if __name__ == "__main__":
    # test_convex_hull()
    hw2()
