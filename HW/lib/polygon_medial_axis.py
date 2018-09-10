#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
use voronoi diagram to approximate polygon medial axis
"""

# standard library
import sys
import math

# third party library
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from matplotlib import path

# local module
from triangulation_DS import compute_circumcircle
from polygon_delaunay_triangulation import plot_polygon, polygon_dt,\
        plot_triangulation


def sample_points_on_boundary(vs, h):
    """ given a polygon defined by ccw vs
        return a sampled boundary
    """
    n_vs = vs.shape[0]
    seg_lengths = [LA.norm(vs[i, :] - vs[(i + 1) % n_vs, :], 2) for i in range(n_vs)]
    new_vs = []
    for i, seg_length in enumerate(seg_lengths):
        dx, dy = vs[(i + 1) % n_vs, :] - vs[i, :]
        n_splits = int(seg_length // h + 1)
        dx, dy = dx / n_splits, dy / n_splits
        x, y = vs[i, :]
        for j in range(n_splits):
            new_vs.append((x, y))
            x, y = x + dx, y + dy
    return np.array(new_vs)


def compute_polygon_medial_axis(vertices, h=0.5):
    """ compute the approximated medial axis of a polygon defined 
        by ccw vertices
    """
    new_vertices = sample_points_on_boundary(vertices, h)
    tri_dt = polygon_dt(new_vertices)
    # plot_triangulation(tri_dt) # pass
    # we have vertices and triangles in tri_dt
    # we need a dictionary that maps edge to face
    visited_edges = set()
    edge_to_face_dic = {}
    vs = tri_dt['vertices']
    for idx_tri, triangle in enumerate(tri_dt['triangles']):
        for i in range(3):
            j = (i + 1) % 3
            idx_a, idx_b = triangle[i], triangle[j]
            if idx_b < idx_a:
                idx_a, idx_b = idx_b, idx_a
            if not (idx_a, idx_b) in visited_edges:
                visited_edges.add((idx_a, idx_b))
                edge_to_face_dic[(idx_a, idx_b)] = [idx_tri]
            else:
                edge_to_face_dic[(idx_a, idx_b)].append(idx_tri)
    # check edge_to_face_dic is correct
    '''
    for (idx_a, idx_b), tri_lst in edge_to_face_dic.items():
        xa, ya = vs[idx_a, :]
        xb, yb = vs[idx_b, :]
        if len(tri_lst) == 1: # boundary edge
            style = 'b-o'
        elif len(tri_lst) == 2: # internal edge
            style = 'r-o'
        else:
            assert(False)
        plt.plot([xa, xb], [ya, yb], style)
    '''
    # compute all the circumcircles
    voronoi_vertices = []
    for triangle in tri_dt['triangles']:
        temp_vs = vs[triangle, :]
        A, B, C = temp_vs[0, :], temp_vs[1, :], temp_vs[2, :]
        P, _ = compute_circumcircle(A, B, C)
        voronoi_vertices.append(P)
    polygon = path.Path(vs) # check
    ins = polygon.contains_points(voronoi_vertices)
    medial_axis = []
    for (idx_a, idx_b), tri_lst in edge_to_face_dic.items():
        if len(tri_lst) == 2: # internal edge
            idx_1, idx_2 = tri_lst
            if ins[idx_1] and ins[idx_2]: # voronoi edge is in polygon
                x1, y1 = voronoi_vertices[idx_1]
                x2, y2 = voronoi_vertices[idx_2]
                medial_axis.append(((x1, y1), (x2, y2)))
        else:
            assert(len(tri_lst) == 1)
    return medial_axis


def plot_polygon_medial_axis(vertices, medial_axis, ax=None):
    """ plot the polygon defined by vertices and its medial axis"""
    if not ax:
        fig, ax = plt.subplots(figsize=(8, 8))
    for (x1, y1), (x2, y2) in medial_axis:
        ax.plot([x1, x2], [y1, y2], 'r--')
    plot_polygon(vertices, ax=ax)


def main(vertices):
    """ show compute polygon medial axis example """
    # new_vs = sample_points_on_boundary(vertices, h=0.2)
    # plot_polygon(new_vs)
    medial_axis = compute_polygon_medial_axis(vertices)
    plot_polygon_medial_axis(vertices, medial_axis)
    plt.axis('equal')
    plt.show()


if __name__ == "__main__":
    vertices = np.array([(0, 0), (4, 0), (4, 1), (0, 1)])
    main(vertices)
