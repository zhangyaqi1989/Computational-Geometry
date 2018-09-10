#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
implements some helper functions used by
polygon triangulation (ear cut algorithm)
"""

# standard library
import sys
from collections import namedtuple

# third party library
import numpy as np
import matplotlib.pyplot as plt


Point2D = namedtuple('Point2D', ['x', 'y'])


def xor(x, y):
    """exclusive or: true off exactly one argument is true
    >>> xor(1, 0)
    True
    >>> xor(0, 1)
    True
    >>> xor(0, 0)
    False
    >>> xor(1, 1)
    False
    """
    return bool(x) != bool(y)


def is_left(a, b, c, tol=1e-12):
    """returns true iff c is strictly to the left of the directed line
       through a to b
    >>> is_left(Point2D(0, 0), Point2D(1, 1), Point2D(1, 0))
    False
    >>> is_left(Point2D(0, 0), Point2D(1, 1), Point2D(0, 1))
    True
    """
    return compute_area_sign(a, b, c, tol) > 0


def is_left_on(a, b, c, tol=1e-12):
    """returns true if c is strictly to the left
       or on the directed line through a to b
    >>> is_left_on(Point2D(0, 0), Point2D(1, 1), Point2D(1, 0))
    False
    >>> is_left_on(Point2D(0, 0), Point2D(1, 1), Point2D(0, 1))
    True
    >>> is_left_on(Point2D(0, 0), Point2D(1, 1), Point2D(0.5, 0.5))
    True
    """
    return compute_area_sign(a, b, c, tol) >= 0


def compute_area_sign(a, b, c, tol=1e-12):
    """return 1 if area is positive, 0 if area is zero, -1 if area is negative
    >>> compute_area_sign(Point2D(0, 0), Point2D(1, 1), Point2D(0, 1))
    1
    >>> compute_area_sign(Point2D(0, 0), Point2D(1, 1), Point2D(1, 0))
    -1
    >>> compute_area_sign(Point2D(0, 0), Point2D(1, 1), Point2D(0.5, 0.5))
    0
    """
    area2 = compute_area2(a, b, c)
    if area2 > tol:
        return 1
    elif area2 < -tol:
        return -1
    else:
        return 0


def compute_area2(a, b, c):
    """returns twice the signed area of the triangle
       determined by a, b, c. the area is positive if
       a, b, c are oriented ccw, negative is cw, and
       zero if the points are collinear
    >>> compute_area2(Point2D(0, 0), Point2D(1, 0), Point2D(1, 1))
    1
    >>> compute_area2(Point2D(0, 0), Point2D(1, 0), Point2D(2, 0))
    0
    >>> compute_area2(Point2D(1, 1), Point2D(1, 0), Point2D(0, 0))
    -1
    """
    return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)


def intersect(a, b, c, d):
    """returns true iff segments ab and cd intersect, properly or improperly
    >>> intersect(Point2D(0, 0), Point2D(1, 1), Point2D(1, 0), Point2D(0, 1))
    True
    >>> intersect(Point2D(0, 0), Point2D(1, 1), Point2D(1, 0), Point2D(2, 0))
    False
    >>> intersect(Point2D(0, 0), Point2D(1, 1), Point2D(1, 1), Point2D(2, 2))
    True
    """
    if intersect_prop(a, b, c, d):
        return True
    elif is_between(a, b, c) or is_between(a, b, d) or\
            is_between(c, d, a) or is_between(c, d, b):
                return True
    else:
        return False


def intersect_prop(a, b, c, d):
    """returns true iff ab properly intersects cd: they
       share a point interior to both segments. The
       properness of the intersection is ensured by using
       strict leftness
    >>> intersect_prop(Point2D(0, 0), Point2D(1, 1), Point2D(1, 0), Point2D(0, 1))
    True
    >>> intersect_prop(Point2D(0, 0), Point2D(1, 1), Point2D(1, 0), Point2D(2, 0))
    False
    """
    if is_collinear(a, b, c) or is_collinear(a, b, d) or is_collinear(c, d, a)\
            or is_collinear(c, d, b):
                return False
    return xor(is_left(a, b, c), is_left(a, b, d)) and xor(
            is_left(c, d, a), is_left(c, d, b))


def is_collinear(a, b, c, tol=1e-12):
    """returns true if c is strictly on
       the directed line through a to b
    >>> is_collinear(Point2D(0, 0), Point2D(1, 1), Point2D(1, 0))
    False
    >>> is_collinear(Point2D(0, 0), Point2D(1, 1), Point2D(0, 1))
    False
    >>> is_collinear(Point2D(0, 0), Point2D(1, 1), Point2D(0.5, 0.5))
    True
    """
    return compute_area_sign(a, b, c, tol) == 0


def is_between(a, b, c):
    """ returns true iff point c lies on the closed segment ab.
        first check that c is collinear with and b
    >>> is_between(Point2D(0, 0), Point2D(1, 1), Point2D(1, 0))
    False
    >>> is_between(Point2D(0, 0), Point2D(1, 0), Point2D(1, 0))
    True
    >>> is_between(Point2D(0, 0), Point2D(1, 0), Point2D(2, 0))
    False
    >>> is_between(Point2D(0, 0), Point2D(0, 1), Point2D(0, 0.5))
    True
    >>> is_between(Point2D(0, 0), Point2D(0, 1), Point2D(0, 2))
    False
    """
    if not is_collinear(a, b, c):
        return False
    if abs(a.x - b.x) > abs(a.y - b.y):
        return (a.x <= c.x and c.x <= b.x) or (a.x >= c.x and c.x >= b.x)
    else:
        return (a.y <= c.y and c.y <= b.y) or (a.y >= c.y and c.y >= b.y)


def doctest():
    """run doctest"""
    import doctest
    doctest.testmod()


if __name__ == "__main__":
    doctest()
