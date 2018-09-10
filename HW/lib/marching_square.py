#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
implements the whole shape processing pipeline
includes boundary reconstruction (marching squares),
simple polygon validation, point membership check,
minimum distance, polygon triangulation and
medial axis
"""

# standard library
import sys
from collections import defaultdict
import math

# 3rd party library
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import path
from skimage.morphology import skeletonize
from skimage.util import invert
from scipy.spatial import Delaunay
from sympy import lambdify
from sympy.parsing.sympy_parser import parse_expr
import seaborn as sns

from scipy.spatial import Voronoi, voronoi_plot_2d

# local library
from triangulation_DS import Triangulation
from polygon_triangulation import polygon_triangulate
from seg_seg import seg_seg_intersect, Point2D
from polygon_delaunay_triangulation import plot_triangulation, polygon_dt
from polygon_medial_axis import compute_polygon_medial_axis, plot_polygon_medial_axis

# sns.set()


def f4(x, y):
    return -(x**2 + y**2)**2 + 2*1.05**2*(x**2 - y**2) + (1 - 1.05**4)


def f1(x, y):
    return 5*x**3 - 17.3*y**2 + np.sin(x*y)

def f3(x, y):
    return x**2/4 + y**2 - 9

def f2(x, y):
    """circle example"""
    return x**2 + y**2 - 4


def compute_root(fun, xs, ys, idx):
    """use linear interpolation to compute root on square"""
    i, j = idx, (idx + 1) % 4
    x1, x2 = xs[i], xs[j]
    y1, y2 = ys[i], ys[j]
    z1, z2 = fun(x1, y1), fun(x2, y2)
    assert(z1 * z2 <= 0)
    rate = (0 - z1) / (z2 - z1)
    x0 = rate * (x2 - x1) + x1
    y0 = rate * (y2 - y1) + y1
    return (x0, y0)


def marching_squares(fun, limits, nx=40, ny=40, plot_points=True):
    """implements marching squares algorithm
       to compute zero contour of function fun
       limits = [xmin, xmax, ymin, ymax]
       nx and ny defines the number of sample points on x-axis and y-axis
       plot_points is for debug purpose
    """
    xmin, xmax, ymin, ymax = limits
    xs = np.linspace(xmin, xmax, nx)
    ys = np.linspace(ymin, ymax, ny)
    zs = np.zeros((nx, ny))
    simple_cases = ['1110', '1101', '1011', '0111', '0001', '0010', '0100', \
            '1000', '1100', '1001', '0011', '0110']
    case_dic = {}
    segs = []
    for case in simple_cases:
        lst = []
        for i in range(4):
            j = (i + 1) % 4
            if case[i] != case[j]:
                lst.append(i)
        case_dic[case] = tuple(lst)
    for i in range(nx):
        for j in range(ny):
            zs[i][j] = fun(xs[i], ys[j])
    if plot_points:
        points = [(xs[i], ys[j], 'r') if zs[i][j] > 0 else (xs[i], ys[j], 'b')\
                for i in range(nx) for j in range(ny)]
        red_points = [point for point in points if point[-1] == 'r']
        blue_points = [point for point in points if point[-1] == 'b']
        points_xs, points_ys, colors = zip(*red_points)
        plt.plot(points_xs, points_ys, 'ro')
        points_xs, points_ys, colors = zip(*blue_points)
        plt.plot(points_xs, points_ys, 'bo')
    for i in range(nx - 1):
        for j in range(ny - 1):
            idx_xs = [i, i+1, i+1, i]
            idx_ys = [j+1, j+1, j, j]
            temp_xs = [xs[idx] for idx in idx_xs]
            temp_ys = [ys[idx] for idx in idx_ys]
            temp_zs = ['1' if zs[x_i][y_i] > 0 else '0' for x_i, y_i in zip(idx_xs, idx_ys)]
            temp_zs = ''.join(temp_zs)
            if temp_zs in case_dic:
                segs.append([compute_root(fun, temp_xs, temp_ys, idx) for idx in case_dic[temp_zs]])
            elif temp_zs in ['1010', '0101']:
                cx, cy = 0.5 * (xs[i] + xs[i + 1]), 0.5 * (ys[i] + ys[i + 1])
                cz = fun(cx, cy)
                if cz > 0:
                    if temp_zs == '1010':
                        segs.append([compute_root(fun, temp_xs, temp_ys, idx) for idx in (0, 1)])
                        segs.append([compute_root(fun, temp_xs, temp_ys, idx) for idx in (2, 3)])
                    else:
                        segs.append([compute_root(fun, temp_xs, temp_ys, idx) for idx in (0, 3)])
                        segs.append([compute_root(fun, temp_xs, temp_ys, idx) for idx in (1, 2)])
                else:
                    if temp_zs == '1010':
                        segs.append([compute_root(fun, temp_xs, temp_ys, idx) for idx in (0, 3)])
                        segs.append([compute_root(fun, temp_xs, temp_ys, idx) for idx in (1, 2)])
                    else:
                        segs.append([compute_root(fun, temp_xs, temp_ys, idx) for idx in (0, 1)])
                        segs.append([compute_root(fun, temp_xs, temp_ys, idx) for idx in (2, 3)])
    return segs


def is_point_close(p0, p1, rel_tol=1e-12):
    """check if two points are close"""
    x0, y0 = p0
    x1, y1 = p1
    return (math.isclose(x0, x1, rel_tol=rel_tol) and 
            math.isclose(y0, y1, rel_tol=rel_tol))


def point_membership_check(polygon, point):
    """check the point is in/out/on the polygon which
       is defined by a list of ordered vertices
    """
    count = 0
    x, y = point
    n = len(polygon)
    assert(n >= 3)
    for i in range(n):
        j = (i + 1) % 3
        if is_point_close(point, polygon[i]) or is_point_close(point, polygon[j]):
            return "on"
        x1, y1 = polygon[i]
        x2, y2 = polygon[j]
        if y >= min(y1, y2) and y <= max(y1, y2):
            if not np.isclose(y1, y2, atol=1e-12):
                x_intersect = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            else:
                if min(x1, x2) < x < max(x1, x2):
                    return 'edge'
                else:
                    continue
            if np.isclose(x, x_intersect, atol=1e-12):
                return 'edge'
            if x < x_intersect and  min(y1, y2) < y:
                count += 1
        elif np.isclose(y, y1, atol=1e-12) and np.isclose(y, y2, atol=1e-12):
            if min(x1, x2) < x < max(x1, x2):
                return 'edge'
    if count % 2 == 0:
        return 'out'
    else:
        return 'in'


def construct_paths(segs):
    """constructed paths from segs"""
    n = len(segs)
    assert(n > 0)
    paths = []
    path = []
    used = [0] * n
    for _ in range(n):
        if not path: # path is empty
            for i, seg in enumerate(segs):
                if used[i]:
                    continue
                else:
                    used[i] = 1
                    P1, P2 = seg
                    path.append(P1)
                    path.append(P2)
                    break
        else:
            P0 = path[-1]
            for i, seg in enumerate(segs):
                if used[i]:
                    continue
                else:
                    P1, P2 = seg
                    if is_point_close(P0, P1):
                        used[i] = 1
                        path.append(P2)
                        break
                    elif is_point_close(P0, P2):
                        used[i] = 1
                        path.append(P1)
                        break
            else:
                print("Polygon construction failed")
        if is_point_close(path[0], path[-1]):
            path = np.array(path)
            path = path[:-1]
            area = 0
            xs, ys = path[:, 0], path[:, 1]
            m = path.shape[0]
            for i in range(m):
                j = (i + 1) % m
                area += xs[i] * ys[j] - xs[j] * ys[i]
            if area < 0:
                path = np.flipud(path)
            paths.append(path)
            path = []
    return paths


def compute_point_to_line_segment(P1, P2, P0):
    """compute point to line distance
       compute distance between P0 and line segment P1P2
    """
    u = (P2 - P1) / LA.norm(P2 - P1, 2) # unit vector of P1P2
    v = P0 - P1
    proj = np.dot(u, v)
    if proj <= 0:
        return LA.norm(P0 - P1, 2)
    elif proj >= 1:
        return LA.norm(P0 - P2, 2)
    else:
        return LA.norm(v - proj*u, 2)


def compute_point_to_polygon(polygon, P):
    """compute the minimum distance from point to polygon
       polygon is defined by ordered points on its boundary
    """
    min_distance = float('inf')
    n = polygon.shape[0]
    for j in range(n):
        k = (j + 1) % n
        P1 = polygon[j]
        P2 = polygon[k]
        distance = compute_point_to_line_segment(P1, P2, P)
        if distance < min_distance:
            min_distance = distance
    return min_distance


def compute_grids(limits, h):
    """ compute nx and ny based on limits and h
    """
    xmin, xmax, ymin, ymax = limits
    nx = np.ceil((xmax - xmin)/(h/np.sqrt(2)))
    ny = np.ceil((ymax - ymin/(h/np.sqrt(2))))
    return int(nx), int(ny)


def is_simple_polygon(polygon):
    """ check if whether a polygon is simple polygon or not """
    n = polygon.shape[0]
    for i in range(n):
        j = (i + 1) % n
        A = polygon[i, :]
        PA = Point2D(A[0], A[1])
        B = polygon[j, :]
        PB = Point2D(B[0], B[1])
        start = (j + 1) % n
        for p in range(n - 3):
            end = (start + 1) % n
            C = polygon[start, :]
            D = polygon[end, :]
            PC = Point2D(C[0], C[1])
            PD = Point2D(D[0], D[1])
            code, _ = seg_seg_intersect(PA, PB, PC, PD)
            if code != '0':
                # for debug
                '''
                fig, ax = plt.subplots()
                ax.plot([PA.x, PB.x], [PA.y, PB.y])
                ax.plot([PC.x, PD.x], [PC.y, PD.y])
                plt.show()
                print(PA, PB, PC, PD)
                '''
                return False
            start = (start + 1) % n
    return True


def project_pipeline(expression, limits, h, points_query):
    ## 0. set up figure axes
    xmin, xmax, ymin, ymax = limits
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    axes = axes.ravel() # flatten axes
    ax1, ax2, ax3, ax4 = axes
    for ax in axes:
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])

    ## 1. use marching square to reconstruct boundary
    ## 1.1 plot the boundary
    nx, ny = compute_grids(limits, h)
    fun = make_function(expression)
    segs = marching_squares(fun, limits, nx=nx, ny=ny, plot_points=False)
    for seg in segs:
        (x1, y1), (x2, y2) = seg
        ax1.plot([x1, x2], [y1, y2], 'b-o')
    ax1.axis('equal')
    ax1.set_title('Construction')
    paths = construct_paths(segs)

    ## 1.2 plot the interior
    patches = []
    for temp_path in paths:
        patch = Polygon(temp_path, True)
        patches.append(patch)
    patch_collections = PatchCollection(patches)
    ax1.add_collection(patch_collections)

    ## 2. validate whether the polygon is a simple polygon or not
    if len(paths) == 0:
        print("there is no closed level set found in the given domain")
        sys.exit(1)
    if len(paths) > 1:
        print("The constructed boundary is not a simple polygon")
        plt.show() # draw the paths
        sys.exit(1)
    polygon = paths[0]
    # check polygon is a simple polygon
    if not is_simple_polygon(polygon):
        print("The constructed boundary is not a simple polygon")
        plt.show() # draw the nonsimple polygon
        sys.exit(1)

    # plot the constructed boundary
    closed_polygon = np.vstack((polygon, polygon[0:1, :]))
    ax2.plot(closed_polygon[:, 0], closed_polygon[:, 1], 'b-o')
    ax2.set_title('Polygon Construction')
    ax2.axis('equal')
    points = polygon
    n = points.shape[0]

    ## 3. point membership check and minimum distance
    for i in range(points_query.shape[0]):
        P = points_query[i]
        min_distance = compute_point_to_polygon(points, P)
        print('{}: {}, {:0.6f}'.format(P, point_membership_check(points, P), min_distance))

    ## 4. Delaunay triangulation
    # 4.1 this is an arbitary polygon triangulation
    '''
    tri = polygon_triangulate(points)
    for (x1, y1), (x2, y2), (x3, y3) in tri:
        ax3.plot([x1, x2], [y1, y2], 'g--')
    temp_points = np.vstack((points, points[0:1, :]))
    ax3.plot(temp_points[:, 0], temp_points[:, 1], 'b-o')
    '''
    # 4.2 this is delaunay triangulation of a set of points
    '''
    tri = Delaunay(points)
    ax3.triplot(points[:, 0], points[:, 1], tri.simplices.copy())
    ax3.plot(points[:, 0], points[:, 1], 'o')
    '''
    # 4.3 this is polygon delaunay triangulation
    # write out the polygon
    # np.savetxt('expression-boundary.txt', points, delimiter='\t')
    tri_dt = polygon_dt(points)
    plot_triangulation(tri_dt, ax3)
    ax3.axis('equal')
    ax3.set_title('Delaunay Triangulation')

    ## 5. approximate medial axis using voronoi diagram
    # this is our voronoi
    '''
    tri = Triangulation(points)
    tri.compute_voronoi()
    tri.plot_medial_axis()
    '''
    # scipy voronoi
    # vor = Voronoi(points)
    # voronoi_plot_2d(vor, ax=ax4)
    medial_axis = compute_polygon_medial_axis(points, h=0.5)
    plot_polygon_medial_axis(points, medial_axis, ax=ax4)
    ax4.axis('equal')
    ax4.set_title('Medial Axis Approximation')
    # ax4.axis("off")
    '''
    original, skeleton = polygon_skeleton(points, limits)
    # ax4.imshow(original, cmap=plt.cm.gray)
    ax4.imshow(invert(skeleton), cmap=plt.cm.gray)
    ax4.set_xticks([])
    ax4.set_yticks([])
    '''
    plt.show()


def user_input():
    """handle user input and read files"""
    if len(sys.argv) != 3:
        print("Usage: >> python {} <expression_file> <points_file>".format(sys.argv[0]))
        sys.exit(1)
    expression_file = sys.argv[1]
    points_file = sys.argv[2]
    points_query = np.loadtxt(points_file)
    with open(expression_file, 'r') as in_file:
        lines = [line for line in in_file.readlines() if len(line.strip()) != 0]
    assert(len(lines) >= 3)
    expression = lines[0].strip()
    limits = [float(token) for token in lines[1].split(', ')]
    h = float(lines[2])
    return expression, limits, h, points_query


def make_function(expression):
    """create a function object f(x, y) from an expression"""
    fun = lambdify([parse_expr('x', evaluate=False), parse_expr('y', evaluate=False)], \
            parse_expr(expression, evaluate=False))
    return fun


def test():
    expression, limits, h, points_query = user_input()
    project_pipeline(expression, limits, h, points_query)


def polygon_skeleton(verts, limits):
    """return the original and skeleton image of the polygon
       polygon is defined by the vertices
    """
    xmin, xmax, ymin, ymax = limits
    polygon = path.Path(verts)
    xs = np.arange(xmin, xmax, 0.05)
    ys = np.arange(ymin, ymax, 0.05)
    points = [(x, y) for x in xs for y in ys]
    zs = polygon.contains_points(points)
    zs = zs.reshape(len(xs), len(ys))
    zs = zs.astype(int)
    original = zs.T
    skeleton = skeletonize(original)
    return original, skeleton


def test_marching_squares():
    # limits = [-2.5, 2.5, -2.5, 2.5]
    limits = [-7, 7, -4, 4]
    xmin, xmax, ymin, ymax = limits
    # limits = [-2, 2, -2, 2]
    nx, ny = 8, 8
    fun = f3
    segs = marching_squares(fun, limits, nx=nx, ny=ny, plot_points=False)
    paths = construct_paths(segs)
    assert(len(paths) == 1)
    verts = paths[0]
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(121)
    xs, ys = zip(*verts)
    ax.plot(xs, ys, "b-o")
    ax.plot([xs[0], xs[-1]], [ys[0], ys[-1]], "b-o")
    ax.axis("equal")

    # sample pixels
    '''
    polygon = path.Path(verts)
    xs = np.arange(xmin, xmax, 0.01)
    ys = np.arange(ymin, ymax, 0.01)
    points = [(x, y) for x in xs for y in ys]
    zs = polygon.contains_points(points)
    zs = zs.reshape(len(xs), len(ys))
    zs = zs.astype(int)
    original = zs.T
    skeleton = skeletonize(image)
    '''
    original, skeleton = polygon_skeleton(verts, limits)
    ax = fig.add_subplot(122)
    ax.imshow(original, cmap=plt.cm.gray)
    ax.imshow(skeleton, cmap=plt.cm.gray)
    ax.axis("equal")
    plt.show()


if __name__ == "__main__":
    test()
    # test_marching_squares()
