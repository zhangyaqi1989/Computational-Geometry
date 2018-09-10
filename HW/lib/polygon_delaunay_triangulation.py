#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
use triangle library to do constrained delaunay
triangulation of a polygon
"""


# standard library
import sys

# third party library
import numpy as np
import matplotlib.pyplot as plt
import triangle


def make_face(vertices):
    """ make face based on polygon ccw vertices """
    n_vertices = len(vertices)
    face = {}
    face['vertices'] = vertices
    segs = [(i, (i + 1) % n_vertices) for i in range(n_vertices)]
    face['segments'] = np.array(segs)
    return face


def plot_polygon(vertices, ax=None):
    """ plot polygon """
    if not ax:
        fig, ax = plt.subplots(figsize=(8, 8))
    assert(vertices.shape[0] >= 3)
    vs = np.vstack((vertices, vertices[0:1, :]))
    ax.plot(vs[:, 0], vs[:, 1], 'b-o', linewidth=2)


def plot_triangulation(tri, ax=None):
    """ plot the triangulation """
    if not ax:
        fig, ax = plt.subplots(figsize=(8, 8))
    vs = tri['vertices']
    for triangle in tri['triangles']:
        plot_polygon(vs[triangle, :], ax=ax)
    ax.axis('equal')


def polygon_dt(vertices):
    """ create delaunay triangulation of a polygon defined by ccw vertices"""
    face = make_face(vertices)
    tri = triangle.triangulate(face, 'p')
    keys_keep = ['vertices', 'triangles']
    keys_delete = []
    for key, value in tri.items():
        if key not in keys_keep:
            keys_delete.append(key)
    for key in keys_delete:
        del tri[key]
    return tri


def main():
    points = [(0, 0), (2, 0), (2, 2), (1, 2), (1, 1), (0, 1)]
    vertices = np.array(points)
    # face = make_face(vertices)
    # tri = triangle.triangulate(face, 'p')
    tri = polygon_dt(vertices)
    fig, ax = plt.subplots(figsize=(8, 8))
    plot_triangulation(tri, ax)
    plt.show()


if __name__ == "__main__":
    main()
