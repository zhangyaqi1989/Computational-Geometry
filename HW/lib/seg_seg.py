#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
implements line segment and line segment intersection
"""

# standard library
import sys

# third party library
import numpy as np
import matplotlib.pyplot as plt

# local library
from algorithm import Point2D, is_collinear, is_between, doctest


def seg_seg_intersect(a, b, c, d):
    """ return
        'e', point : the segments collinearly overlap, sharing a point
        'v', point : the endpoint of one segment is on the other segment
                     but 'e' does not hold
        '1', point : the segments intersect properly (i.e. they share a point 
                     and neither 'v' nor 'e' holds)
        '0', None or point : the segments do not intersect
    >>> seg_seg_intersect(Point2D(0, 0), Point2D(1, 1), Point2D(1, 0), Point2D(0, 1))
    ('1', Point2D(x=0.5, y=0.5))
    >>> seg_seg_intersect(Point2D(0, 0), Point2D(1, 1), Point2D(1, 0), Point2D(0.9, 0.1))
    ('0', Point2D(x=0.5000000000000001, y=0.5000000000000001))
    >>> seg_seg_intersect(Point2D(0, 0), Point2D(1, 1), Point2D(1, 0), Point2D(0.5, 0.5))
    ('v', Point2D(x=0.5, y=0.5))
    >>> seg_seg_intersect(Point2D(0, 0), Point2D(1, 1), Point2D(1, 1), Point2D(2, 2))
    ('e', Point2D(x=1, y=1))
    """
    denom = a.x * (d.y - c.y) + b.x * (c.y - d.y) + d.x * (b.y - a.y) \
                + c.x * (a.y - b.y)
    if denom == 0:
        return _parallel_intersect(a, b, c, d)
    num = a.x * (d.y - c.y) + c.x * (a.y - d.y) + d.x * (c.y - a.y)
    if num == 0.0 or num == denom:
        code = 'v'
    s = num / denom
    num = - (a.x * (c.y - b.y) + b.x * (a.y - c.y) + c.x * (b.y - a.y))
    if num == 0.0 or num == denom:
        code = 'v'
    t = num / denom
    if 0.0 < s < 1.0 and 0.0 < t < 1.0:
        code = '1'
    elif s < 0.0 or s > 1.0 or t < 0.0 or t > 1.0:
        code = '0'
    x = a.x + s * (b.x - a.x)
    y = a.y + s * (b.y - a.y)
    return code, Point2D(x=x, y=y)


def _parallel_intersect(a, b, c, d):
    """ returns 'e', Point when ab and cd collinearly overlap sharing
        an point, otherwise return '0', None
    >>> _parallel_intersect(Point2D(0, 0), Point2D(1, 1), Point2D(1, 1), Point2D(2, 1))
    ('e', Point2D(x=1, y=1))
    >>> _parallel_intersect(Point2D(0, 0), Point2D(1, 1), Point2D(1, 0), Point2D(2, 1))
    ('0', None)
    """
    if not is_collinear(a, b, c):
        return '0', None
    if is_between(a, b, c):
        point = c
    if is_between(a, b, d):
        point = d
    if is_between(c, d, a):
        point = a
    if is_between(c, d, b):
        point = b
    return 'e', point


if __name__ == "__main__":
    doctest()
