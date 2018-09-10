#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
define some basic algorithms for lines
"""

import numpy as np
import matplotlib.pyplot as plt
from geometry import Point, Line

def is_left(line, point, tol=1e-12):
    (x0, y0), (x1, y1) = line.get_coords()
    x, y = point.get_coords()
    v = (x1 - x0) * (y - y0) - (y1 - y0) * (x - x0)
    if abs(v) < tol:
        return 0
    elif v > 0:
        return 1
    else:
        return -1


def is_line_cross(line1, line2, tol=1e-12):
    v1 = is_left(line1, line2.p0, tol)
    v2 = is_left(line1, line2.p1, tol)
    if v1 * v2 <= 0:
        return True
    else:
        return False

def is_line_parallel(line1, line2, tol=1e-12):
    (x1, y1), (x2, y2) = line1.get_coords()
    (x3, y3), (x4, y4) = line2.get_coords()
    v = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(v) < tol:
        return True
    else:
        return False


def is_line_segments_intersect(line1, line2, tol=1e-12):
    (x1, y1), (x2, y2) = line1.get_coords()
    (x3, y3), (x4, y4) = line2.get_coords()
    a = is_left(line1, Point(x3, y3), tol) * is_left(line1, Point(x4, y4), tol)
    b = is_left(line2, Point(x1, y1), tol) * is_left(line2, Point(x2, y2), tol)
    return a == -1 and b == -1


def compute_line_intersect(line1, line2, tol=1e-12):
    (x1, y1), (x2, y2) = line1.get_coords()
    (x3, y3), (x4, y4) = line2.get_coords()
    denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(denominator) < tol:
        return None
    else:
        numerator1 = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)
        numerator2 = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)
        intersect_point = Point(numerator1/denominator, numerator2/denominator)
        return intersect_point


if __name__ == "__main__":
    line1 = Line(0, 0, 1, 1)
    point = Point(np.sqrt(2), np.sqrt(2))
    print(is_left(line1, point))
    line2 = Line(1, 0, 1, 1)
    print(is_line_parallel(line1, line2))
    print(compute_line_intersect(line1, line2))
    line3 = Line(1, 0, 0, 1)
    line4 = Line(0, 0, 0.2, 0.2)
    print(is_line_segments_intersect(line4, line3))
    print(compute_line_intersect(line3, line4))

