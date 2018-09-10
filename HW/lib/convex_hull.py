#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
implement different convex hull algorithms,
diameter, width of a set of 2D points
"""

import sys
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from polygon import rotating_calipers


def compute_polygon_area(points):
    ''' shoelace algorithm to compute polygon area 
        points are ordered
    '''
    if type(points) == list:
        xs, ys = zip(*points)
        n = len(xs)
    elif type(points) == np.ndarray:
        xs = points[:, 0]
        ys = points[:, 1]
        n = points.shape[0]
    else:
        print("unsupported points type {:s}".format(type(points)))
        sys.exit(1)
    area = 0
    for i in range(n):
        j = (i + 1) % n
        area += xs[i] * ys[j] - xs[j] * ys[i]
    return 0.5 * np.abs(area)


def convex_hull(points, alg=1):
    """ compute convex hull using specified alg """
    points_t = type(points)
    if points_t != np.ndarray and points_t != list:
        print("unsupported points type {:s}".format(points_t))
        sys.exit(1)
    if points_t == np.ndarray:
        input_points = list(zip(points[:, 0], points[:, 1]))
    else:
        input_points = points[:]
    fun_dict = {1:graham_scan, 2:jarvis, 3:monotone_chain}
    try:
        fun = fun_dict[alg]
    except (KeyError):
        print("unsupported algorithem type {:d}".format(alg))
        print("alg = 1: graham scan")
        print("alg = 2: jarvis march")
        print("alg = 3: monotone chain")
        sys.exit(1)
    return fun(input_points)


def _is_left(A, B, point, tol=1e-12):
    """ test  if point in on the left of AB """
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


def jarvis(points):
    """ jarvis march CH alg """
    n_points = len(points)
    left_index = 0
    for i in range(1, n_points):
        if points[i][0] < points[left_index][0]:
            left_index = i
    point_on_hull = points[left_index]
    i = 0
    hull = []
    while True:
        hull.append(point_on_hull)
        end_point = points[0]
        for j in range(1, n_points):
            if (end_point == point_on_hull) or _is_left(hull[i], end_point, points[j]) == 1:
                end_point = points[j]
        i = i + 1
        point_on_hull = end_point
        if end_point == hull[0]:
            break
    return hull


def _ccw(P1, P2, P3):
    """ test if P1-P2-P3 is ccw """
    return (P2[0] - P1[0]) * (P3[1] - P1[1]) - (P2[1] - P1[1]) * (P3[0] - P1[0])


def _cosine(origin, point):
    """ used for sorting in graham scan """
    vector = np.array(point) - np.array(origin)
    return vector[0] / np.linalg.norm(vector)


def graham_scan(points):
    """ graham scan CH alg """
    n_points = len(points)
    point = min(points, key=lambda x : (x[1], x[0]))
    points.remove(point)
    points.sort(key=lambda x : _cosine(point, x), reverse=True)
    points = [points[-1], point] + points
    M = 1
    for i in range(2, n_points + 1):
        while _ccw(points[M - 1], points[M], points[i]) <= 0:
            if M > 1:
                M -= 1
                continue
            elif i == N:
                break
            else:
                i += 1
        M += 1
        x, y = points[i]
        points[i] = points[M]
        points[M] = x, y
    hull = points[1 : M + 1]
    return hull


def monotone_chain(points):
    """ monotone chain CH alg """
    points = sorted(set(points))
    if len(points) <= 1:
        return points
    lower = []
    for point in points:
        while len(lower) >= 2 and _ccw(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)

    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and _ccw(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def compute_line_intersect(A, B, C, D, tol=1e-12):
    """ compute intersect points of AB and CD """
    x1, y1 = A
    x2, y2 = B
    x3, y3 = C
    x4, y4 = D
    denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(denominator) < tol:
        return None
    else:
        numerator1 = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)
        numerator2 = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)
        intersect_point = np.array([numerator1/denominator, numerator2/denominator])
        return intersect_point


def axes_aligned_bounding_rectangle(points):
    """ compute axes axes aligned bounding rectangle """
    if type(points) != np.ndarray:
        points = np.array(points)
    xmin, ymin = points.min(axis=0)
    xmax, ymax = points.max(axis=0)
    xs = [xmin, xmax, xmax, xmin]
    ys = [ymin, ymin, ymax, ymax]
    box = list(zip(xs, ys))
    return box


def bounding_rectangle(hull):
    """ compute bounding rectangle """
    hull = np.array(hull)
    M = hull.shape[0]
    area = float('Inf')
    width = float('Inf')
    for i in range(M):
        j = (i + 1) % M
        A = hull[i, :]
        B = hull[j, :]
        height = -float('Inf')
        left = float('Inf')
        right = -float('Inf')
        unit = B - A
        unit = unit / norm(unit)
        for idx in range(M):
            C = hull[idx, :]
            vector = C - A
            para_len = np.dot(vector, unit)
            if para_len > right:
                right = para_len
                best_right = idx
            if para_len < left:
                left = para_len
                best_left = idx
            para = para_len * unit
            perp = vector - para
            if norm(perp) > height:
                height = norm(perp)
                best_height = idx
        curr_area = (right - left) * height
        # print("Height = {}".format(height))
        # print("curr_area = " + str(curr_area))
        if curr_area < area:
            area = curr_area
            record = [best_left, best_right, best_height, i, area]
    best_left, best_right, best_height, i, area = record
    x0, y0 = hull[i, :]
    x1, y1 = hull[(i + 1) % M, :]
    # print(width)
    # plt.plot([x0, x1], [y0, y1], 'r-')
    # plt.scatter(hull[best_height, 0], hull[best_height, 1], color='blue')
    # plt.scatter(hull[best_left, 0], hull[best_left, 1], color='black')
    # plt.scatter(hull[best_right, 0], hull[best_right, 1], color='red')
    A = hull[i, :]
    B = hull[(i + 1) % M, :]
    LP = hull[best_left, :]
    RP = hull[best_right, :]
    dx, dy = B - A
    TP = hull[best_height, :]
    LL = compute_line_intersect(A, B, LP, LP + np.array([-dy, dx]))
    LR = compute_line_intersect(A, B, RP, RP + np.array([-dy, dx]))
    UR = compute_line_intersect(TP, TP + np.array([dx, dy]), RP, RP + np.array([-dy, dx]))
    UL = compute_line_intersect(TP, TP + np.array([dx, dy]), LP, LP + np.array([-dy, dx]))
    # plt.scatter(LL[0], LL[1], color='yellow')
    # plt.scatter(LR[0], LR[1], color='yellow')
    # plt.scatter(UR[0], UR[1], color='yellow')
    # plt.scatter(UL[0], UL[1], color='yellow')
    # x1, y1 = LL
    # x2, y2 = LR
    # x3, y3 = UR
    # x4, y4 = UL
    # plt.plot([x1, x2, x3, x4, x1], [y1, y2, y3, y4, y1], 'b--')
    return np.array([LL, LR, UR, UL]), area


def compute_diameter(points):
    """ compute diameter of a set of points """
    pairs = rotating_calipers(points[:])
    D, pair = max([(((p[0] - q[0])**2 + (p[1] - q[1])**2), (p, q)) for \
            p, q in pairs])
    return (np.sqrt(D), pair)


def compute_width(points):
    """ compute width of a set of points """
    hull = convex_hull(points, 1)
    hull = np.array(hull)
    M = hull.shape[0]
    best_width = float('Inf')
    for i in range(M):
        j = (i + 1) % M
        A = hull[i, :]
        B = hull[j, :]
        curr_best_width = -float('Inf')
        unit = B - A
        unit = unit / norm(unit)
        for idx in range(M):
            C = hull[idx, :]
            vector = C - A
            para_len = np.dot(vector, unit)
            para = para_len * unit
            perp = vector - para
            if norm(perp) > curr_best_width:
                curr_best_width = norm(perp)
                best_width_idx = idx
        if best_width > curr_best_width:
            best_width = curr_best_width
            record = [best_width, i, best_width_idx]
    best_width, i, best_width_idx = record
    x0, y0 = hull[i, :]
    x1, y1 = hull[(i + 1) % M, :]
    A = hull[i, :]
    B = hull[(i + 1) % M, :]
    dx, dy = B - A
    TP = hull[best_width_idx, :]
    mid_point = compute_line_intersect(A, B, TP, TP + np.array([-dy, dx]))
    return best_width, i, best_width_idx, np.array([TP, mid_point])


def main(n):
    """ generate n 2D points randomly 
        compute convex hull, bounding box, axes aligned bounding box,
        diameter and width.
    """
    # np.random.seed(1)
    points = np.random.uniform(0, 5, (20, 2)) # test ndarray
    xs = points[:,0]
    ys = points[:,1]
    points = list(zip(xs, ys)) # test points
    hull = convex_hull(points, 3)
    print("convex hull area = {:0.6f}".format(compute_polygon_area(hull)))
    plt.scatter(xs, ys, color='g')
    M = len(hull)
    rectangle, area = bounding_rectangle(hull)
    print("bounding rectangle area = {:0.6f}".format(area))
    box = axes_aligned_bounding_rectangle(points)
    print("bounding rectangle area (axes aligned) = {:0.6f}".format(compute_polygon_area(box)))
    for i in range(M):
        j = (i + 1) % M
        plt.plot([hull[i][0], hull[j][0]], [hull[i][1], hull[j][1]], 'r--')
    for i in range(4):
        j = (i + 1) % 4
        plt.plot([rectangle[i][0], rectangle[j][0]], [rectangle[i][1], rectangle[j][1]], 'b--')
        plt.plot([box[i][0], box[j][0]], [box[i][1], box[j][1]], 'g--')
    pairs = rotating_calipers(points[:])
    D, pair = compute_diameter(points)
    (x0, y0), (x1, y1) = pair
    plt.plot([x0, x1], [y0, y1], 'k--')
    plt.scatter([x0, x1], [y0, y1], color="red")
    print("diameter = {:0.6f}".format(D))
    best_width, i, best_width_idx, ends = compute_width(points)
    print("width = {:0.6f}".format(best_width))
    plt.plot(ends[:, 0], ends[:, 1], 'g--')
    plt.scatter(ends[:, 0], ends[:, 1], color="r")
    plt.axis('equal')
    plt.show()


def test_CH():
    points = np.loadtxt("../hw3-1/points-4.txt")
    xs = points[:,0]
    ys = points[:,1]
    plt.scatter(xs, ys, color='blue')
    points = list(zip(xs, ys)) # test points
    hull = convex_hull(points, alg=2)
    M = len(hull)
    for i in range(M):
        j = (i + 1) % M
        plt.plot([hull[i][0], hull[j][0]], [hull[i][1], hull[j][1]], 'r--o')
    plt.show()


if __name__ == "__main__":
    n = 30
    main(n)
    # test_CH()
