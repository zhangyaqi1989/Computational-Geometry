#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
define some basic geometry class Point, Line, Polygon
"""

from math import isclose

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        # return f'({self.x}, {self.y})'
        return "(%f, %f)" % (self.x, self.y)

    def get_coords(self):
        return self.x, self.y

    def __eq__(self, other):
        # return self.x == other.x and self.y == other.y
        return isclose(self.x, other.x) and isclose(self.y, other.y)

class Line:
    def __init__(self, x0, y0, x1, y1):
        self.p0 = Point(x0, y0)
        self.p1 = Point(x1, y1)

    def get_coords(self):
        return self.p0.get_coords(), self.p1.get_coords()

    def __str__(self):
        return (str(self.p0) + " --> " + str(self.p1))


class Polygon:

    def __init__(self, xs, ys):
        assert(len(xs) == len(ys))
        self.xs = xs[:]
        self.ys = ys[:]
        self.n_vertices = len(xs)

    def __str__(self):
        lst = []
        for x, y in zip(self.xs, self.ys):
            lst.append("(%f, %f)" % (x, y))
            lst.append(" --> ")
        return "".join(lst[:-1])

    def compute_area(self):
        area = 0
        n = self.n_vertices
        for i in range(n):
            j = (i + 1) % n
            area += xs[i] * ys[j] - xs[j] * ys[i]
        return 0.5 * abs(area)

if __name__ == "__main__":
    p = Point(1, 2)
    print(p)
    line = Line(0, 0, 1, 1)
    print(line)
    xs = [0, 1, 1]
    ys = [0, 0, 1]
    polygon = Polygon(xs, ys)
    print(polygon)
    print(polygon.compute_area())
