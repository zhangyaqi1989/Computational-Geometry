#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
Billiard game simulation
"""

import sys
from fractions import gcd
sys.path.insert(0, '../lib/')
import numpy as np
import matplotlib.pyplot as plt
from geometry import Point, Line
from geometry_algorithm import compute_line_intersect, is_line_segments_intersect

class Billiard:
    def __init__(self, m, n, px, py, qx, qy, dx, dy):
        '''init board, calculate path and intersections'''
        self.m = m
        self.n = n
        self.pq = Line(px, py, qx, qy)
        xs = [0, m, m, 0, 0]
        ys = [0, 0, n, n, 0]
        self.boundarys = []   # board boundary
        self.vertices = []    # vertices on board
        self.norms = ([np.array([0, 1]), np.array([-1, 0]), np.array([0, -1]), 
                np.array([1, 0])])
        self.points = [Point(0, 0)] # path points
        self.n_inters = 0
        self.inter_points = []
        for i in range(4):
            line = Line(xs[i], ys[i], xs[i + 1], ys[i + 1])
            point = Point(xs[i], ys[i])
            self.boundarys.append(line)
            self.vertices.append(point)
        x, y = 0, 0
        finished = False
        max_points = m + n - 2 + 1
        curr = 0 # avoid inf loop
        while not finished and curr <= max_points:
            curr += 1
            line = Line(x, y, x + dx, y + dy)
            for index, bound_line in enumerate(self.boundarys):
                point = compute_line_intersect(line, bound_line)
                if point is not None and point != self.points[-1]:
                    temp_x, temp_y = point.get_coords()
                    if 0 <= temp_x <= m and 0 <= temp_y <= n:
                        if is_line_segments_intersect(Line(x, y, point.x, point.y), self.pq):
                            inter_point = compute_line_intersect(Line(x, y, point.x, point.y), self.pq)
                            self.n_inters += 1
                            self.inter_points.append(inter_point)
                        x, y = temp_x, temp_y
                        self.points.append(point)
                        # general bounding alg
                        v_in = np.array([dx, dy])
                        norm = self.norms[index]
                        para = np.dot(v_in, norm) * norm
                        perp = v_in - para
                        v_out = perp - para
                        dx, dy = v_out
                        if point in self.vertices:
                            finished = True
                        break


    def _plot_line(self, point1, point2):
        '''plot line (point1, point2)'''
        x1, y1 = point1.get_coords()
        x2, y2 = point2.get_coords()
        plt.plot([x1, x2], [y1, y2], 'b-')


    def show(self):
        ''' create animation '''
        plt.axis("equal")
        # 1. plot bound
        plt.plot([0, self.m, self.m, 0, 0], [0, 0, self.n, self.n, 0])
        plt.draw()
        plt.pause(0.2)
        # 2. plot path
        n_points = len(self.points)
        for i in range(n_points - 1):
            self._plot_line(self.points[i], self.points[i + 1])
            plt.draw()
            plt.pause(0.001)
        # 3. plot PQ
        (px, py), (qx, qy) = self.pq.get_coords()
        plt.plot([px, qx], [py, qy], 'g--')
        # 4. plot intersect points
        for point in self.inter_points:
            plt.scatter(point.x, point.y, color='r')
            plt.draw()
        for point in self.points:
            print(point)
        print(self.n_inters)
        # show some info if m and n are not relatively prime
        if gcd(self.m, self.n) != 1:
            print("m and n are not relatively prime")
        plt.axis('equal')
        plt.show()


def hw1():
    '''
    create a billiard table and simulate bounding ball
    output bounding points and intersection number
    and visualize
    '''
    if len(sys.argv) != 7:
        print("Usage: >> python {} <m> <n> <px> <py> <qx> <qy>".format(sys.argv[0]))
        sys.exit(1)
    else:
        m, n, px, py, qx, qy = (float(arg) for arg in sys.argv[1:])
    # assume 45 degree shooting angle
    dx = 1
    dy = 1
    billiard = Billiard(m, n, px, py, qx, qy, dx, dy)
    billiard.show()


if __name__ == "__main__":
    hw1()
