#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
Define a Triangulation Data Structure which can do
naive triangulation, Delaunay triangulation and voronoi diagram
and crust algorithm and shape procedure
"""

import sys
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay, Voronoi, voronoi_plot_2d
from convex_hull import convex_hull
from triangulation import triangulation, is_in


def _is_point_in_shape(points, point):
    """check if point is in the polygon bounded by points"""
    n = len(points)
    for i in range(n):
        j = (i + 1) % n
        A, B = points[i], points[j]
        if _is_left(A, B, point) < 0:
            return False
    return True


def _print_dict(d):
    """print dictinary"""
    keys = list(d.keys())
    keys.sort()
    for key in keys:
        print("{} --> {}".format(key, d[key]))


def _point_in_hull(hull, point):
    """check if a point is in hull or not"""
    n = len(hull)
    for i in range(n):
        j = (i + 1) % n
        if _is_left(hull[i], hull[j], point) < 0:
            return False
    return True


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


def compute_circumcircle(A, B, C):
    """compute circumcircle of an triangle ABC"""
    (x1, y1), (x2, y2), (x3, y3) = A, B, C
    A = np.array([[x3 - x1, y3 - y1],[x3 - x2, y3 - y2]])
    Y = np.array([(x3**2 + y3**2 - x1**2 - y1**2), (x3**2 + y3**2 - x2**2 - y2**2)])
    if np.linalg.det(A) == 0:
        return False
    Ainv = np.linalg.inv(A)
    X = 0.5 * np.dot(Ainv,Y)
    x, y = X[0], X[1]
    r = np.sqrt((x - x1)**2 + (y - y1)**2)
    return (x, y), r


def _plot_circumcircle(A, B, C, D):
    """compute and plot circumcircle"""
    (x1, y1), (x2, y2), (x3, y3) = A, B, C
    (x, y), r = compute_circumcircle(A, B, C)
    n = 100
    ts = np.linspace(0, 2*np.pi, n)
    xs = x + r * np.cos(ts)
    ys = y + r * np.sin(ts)
    plt.plot(xs, ys, 'b-')
    plt.scatter([x1, x2, x3, x], [y1, y2, y3, y], color='r')
    x4, y4 = D
    plt.scatter(x4, y4, color='g')
    plt.axis('equal')


def _is_left(A, B, point, tol=1e-12):
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


def _ccw(A, B, C):
    if _is_left(A, B, C) == -1:
        return (A, C, B)
    else:
        return (A, B, C)


def _is_incircumcircle(A, B, C, D):
    A, B, C = _ccw(A, B, C)
    (x1, y1), (x2, y2), (x3, y3), (x4, y4) = A, B, C, D
    M = np.array([[x1 - x4, y1 - y4, (x1 - x4)**2 + (y1 - y4)**2], \
            [x2 - x4, y2 - y4, (x2 - x4)**2 + (y2 - y4)**2], \
            [x3 - x4, y3 - y4, (x3 - x4)**2 + (y3 - y4)**2]])
    return np.linalg.det(M) > 0


class Triangulation:
    """ a triangulation data struture which supports
        delaunay triangulation, voronoi diagram, crust algorithm
        shape procedure of a set of point
    """

    def __init__(self, points):
        self.vertices = np.array(points)
        self.simplices = None # same as numpy
        self.vertex_faces_dict = defaultdict(set)
        self.vertex_edges_dict = defaultdict(set)
        self.edge_faces_dict = defaultdict(set)
        self.v_dict = {}
        n = self.vertices.shape[0]
        for i in range(n):
            row = self.vertices[i, :]
            self.v_dict[row[0], row[1]] = i
        triangles = triangulation(self.vertices) # naive triangulation
        n_tris = len(triangles)
        self.simplices = np.zeros((n_tris, 3), dtype=np.int32)
        self.simplices_DT = None
        self.tri_set = set()
        for i, triangle in enumerate(triangles):
            P0, P1, P2 = triangle
            idx0 = self.v_dict[P0]
            idx1 = self.v_dict[P1]
            idx2 = self.v_dict[P2]
            self.simplices[i, 0] = idx0
            self.simplices[i, 1] = idx1
            self.simplices[i, 2] = idx2
            tri = (idx0, idx1, idx2)
            self.edge_faces_dict[self.order_edge(idx0, idx1)].add(tri)
            self.edge_faces_dict[self.order_edge(idx1, idx2)].add(tri)
            self.edge_faces_dict[self.order_edge(idx2, idx0)].add(tri)
            self.vertex_faces_dict[idx0].add(tri)
            self.vertex_faces_dict[idx1].add(tri)
            self.vertex_faces_dict[idx2].add(tri)
            self.tri_set.add(tri)
        for edge in self.edge_faces_dict.keys():
            idx0, idx1 = edge
            self.vertex_edges_dict[idx0].add(edge)
            self.vertex_edges_dict[idx1].add(edge)
        # _print_dict(self.vertex_faces_dict)
        # _print_dict(self.vertex_edges_dict)


    def get_nedges(self):
        return len(self.edge_faces_dict)


    def get_nvertices(self):
        return self.vertices.shape[0]


    def get_nfaces(self, face='naive'):
        if face.lower() == 'naive':
            return self.simplices.shape[0]
        elif face.lower() == 'delaunay':
            return self.simplices_DT.shape[0]
        else:
            print("Unsupported face type {}".format(face))
            sys.exit(1)


    def plot_naive_triangulation(self):
        plt.scatter(self.vertices[:, 0], self.vertices[:, 1], color='red')
        plt.triplot(self.vertices[:, 0], self.vertices[:, 1], self.simplices.copy(), color='blue')


    def plot_delaunay(self):
        if self.simplices_DT is None:
            self.create_delaunay()
        # ref_tri = Delaunay(self.vertices) # check
        # self.check_delaunay()
        plt.scatter(self.vertices[:, 0], self.vertices[:, 1], color='red')
        # plt.triplot(self.vertices[:, 0], self.vertices[:, 1], ref_tri.simplices.copy(), color='green')
        plt.triplot(self.vertices[:, 0], self.vertices[:, 1], self.simplices_DT.copy(), color='blue')


    def check_delaunay(self):
        for key, value in self.edge_faces_dict.items():
            print("{} --> {}".format(key, value))
        edges = list(self.edge_faces_dict.keys())
        for edge in edges:
            tris = self.edge_faces_dict[edge]
            if len(tris) == 1:
                continue
            elif len(tris) == 2:
                tri1, tri2 = tris
                idx1, idx2, idx3 = tri1
                A, B, C = self.vertices[idx1, :], self.vertices[idx2, :], self.vertices[idx3, :]
                idxp = list(set(tri2) - set(edge))[0]
                P = self.vertices[idxp, :]
                is_in = _is_incircumcircle(A, B, C, P)
                if is_in:
                    print("bad")
            else:
                # print("edge is " + edge)
                print(edge, tris)


    def get_vertex(self, index):
        return self.vertices[index, 0], self.vertices[index, 1]


    def get_edge(self, edge_indexes):
        return (self.get_vertex(edge_indexes[0]), self.get_vertex(edge_indexes[1]))


    def get_triangle(self, tri_indexes):
        return (self.get_vertex(tri_indexes[0]), self.get_vertex(tri_indexes[1]), self.get_vertex(tri_indexes[2]))


    def print_edges(self):
        """print edges"""
        for edge in self.edge_faces_dict.keys():
            print("{}, {}".format(self.v_str(edge[0]), self.v_str(edge[1])))


    def v_str(self, index):
        """change vertex to string"""
        return "({0:0.4f}, {1:0.4f})".format(self.vertices[index, 0], self.vertices[index, 1])


    def print_faces(self, face="naive"):
        """print faces"""
        if face.lower() == "naive":
            simplicies = self.simplices
        elif face.lower() == "delaunay":
            simplicies = self.simplices_DT
        else:
            print("Unsupported face type {}".format(face))
            sys.exit(1)
        nfaces = simplicies.shape[0]
        v = self.vertices
        for i in range(nfaces):
            face = simplicies[i, :]
            print("{}, {}, {}".format(self.v_str(face[0]), \
                    self.v_str(face[1]), self.v_str(face[2])))


    def create_delaunay(self):
        """create delaunay triangulation using edge flipping"""
        edges = list(self.edge_faces_dict.keys())
        edge_enstack = {}
        for edge in edges:
            edge_enstack[edge] = 1
        count = 0
        while len(edges) != 0:
            edge = edges.pop()
            edge_enstack[edge] = 0
            if edge not in self.edge_faces_dict:
                continue
            tris = self.edge_faces_dict[edge]
            if len(tris) == 1:
                continue
            elif len(tris) == 2:
                tri1, tri2 = tris
                idx1, idx2, idx3 = tri1
                A, B, C = self.vertices[idx1, :], self.vertices[idx2, :], self.vertices[idx3, :]
                idxp = list(set(tri2) - set(edge))[0]
                P = self.vertices[idxp, :]
                is_in = _is_incircumcircle(A, B, C, P)
                if is_in:
                    idxq = list(set(tri1) - set(edge))[0]
                    self.flip_and_update(tris, edge, idxp, idxq)
                    # edges.append(self.order_edge(idxp, idxq))
                    # assert(not edge in self.edge_faces_dict)
                    idxr, idxs = edge
                    cands = [(idxp, idxr), (idxq, idxr), (idxp, idxs), (idxq, idxs)]
                    cands = [self.order_edge(idxa, idxb) for idxa, idxb in cands]
                    for cand in cands:
                        if cand in edge_enstack:
                            if edge_enstack[cand] == 0:
                                edge_enstack[cand] = 1
                                edges.append(cand)
                        else:
                            edge_enstack[cand] = 1
                            edges.append(cand)
                    # edges.append(self.order_edge(idxp, idxr))
                    # edges.append(self.order_edge(idxq, idxr))
                    # edges.append(self.order_edge(idxp, idxs))
                    # edges.append(self.order_edge(idxq, idxs))
            else: # for debugging
                pass
                # print("hello")
                # print(edge, tris)
        self.simplices_DT = np.array(list(self.tri_set))
        # update vertex_edges_dict and vertex_faces_dict
        self.vertex_edges_dict.clear()
        self.vertex_faces_dict.clear()
        for edge in self.edge_faces_dict.keys():
            idx0, idx1 = edge
            self.vertex_edges_dict[idx0].add(edge)
            self.vertex_edges_dict[idx1].add(edge)
        for tri in self.tri_set:
            idx0, idx1, idx2 = tri
            self.vertex_faces_dict[idx0].add(tri)
            self.vertex_faces_dict[idx1].add(tri)
            self.vertex_faces_dict[idx2].add(tri)
        # _print_dict(self.vertex_faces_dict)
        # _print_dict(self.vertex_edges_dict)


    def plot_medial_axis(self):
        """plot approximate medial axis using voronoi diagram"""
        xmin, xmax = round(self.vertices[:, 0].min()), round(self.vertices[:, 0].max())
        ymin, ymax = round(self.vertices[:, 1].min()), round(self.vertices[:, 1].max())
        margin = min(ymax - ymin, xmax - xmin) / 10
        xmin, xmax = xmin - margin, xmax + margin
        ymin, ymax = ymin - margin, ymax + margin
        plt.xlim([xmin, xmax])
        plt.ylim([ymin, ymax])
        # print("[%f, %f] X [%f, %f]" % (xmin, xmax, ymin, ymax))
        # 1. plot the vertices
        plt.scatter(self.vertices[:, 0], self.vertices[:, 1], color="blue")
        # 2. plot delaunay diagram
        # plt.triplot(self.vertices[:, 0], self.vertices[:, 1], self.simplices_DT.copy(), 'b--')
        # 3. plot voronoi vertices
        # plt.scatter(self.voronoi_vertices[:, 0], self.voronoi_vertices[:, 1], color="orange")
        vs = self.voronoi_vertices # alias
        # 4. plot voronoi edges
        for idxa, idxb in self.voronoi_edges:
            A, B = vs[idxa], vs[idxb]
            if _is_point_in_shape(self.vertices, A) and _is_point_in_shape(self.vertices, B):
                plt.plot([vs[idxa, 0], vs[idxb, 0]], [vs[idxa, 1], vs[idxb, 1]], color='red')
            # else:
            #    plt.plot([vs[idxa, 0], vs[idxb, 0]], [vs[idxa, 1], vs[idxb, 1]], color='blue')
        # margin = min(ymax - ymin, xmax - xmin) / 10
        # 5. plot voronoi rays
        '''
        for (xc, yc), (dx, dy) in self.voronoi_rays:
            if xmin <= xc <= xmax and ymin <= yc <= ymax:
                if abs(dx) > abs(dy):
                    xfar = xmax + margin if dx > 0 else xmin - margin
                    yfar = yc + (xfar - xc) / dx * dy
                else:
                    yfar = ymax + margin if dy > 0 else ymin - margin
                    xfar = xc + (yfar - yc) / dy * dx
                plt.plot([xc, xfar], [yc, yfar], "b--")
        '''

    def plot_voronoi(self):
        """plot voronoi diagram"""
        xmin, xmax = round(self.vertices[:, 0].min()), round(self.vertices[:, 0].max())
        ymin, ymax = round(self.vertices[:, 1].min()), round(self.vertices[:, 1].max())
        xmin, xmax = xmin - 1, xmax + 1
        ymin, ymax = ymin - 1, ymax + 1
        plt.xlim([xmin, xmax])
        plt.ylim([ymin, ymax])
        # print("[%f, %f] X [%f, %f]" % (xmin, xmax, ymin, ymax))
        # 1. plot the vertices
        plt.scatter(self.vertices[:, 0], self.vertices[:, 1], color="blue")
        # 2. plot delaunay diagram
        # plt.triplot(self.vertices[:, 0], self.vertices[:, 1], self.simplices_DT.copy(), 'b--')
        # 3. plot voronoi vertices
        plt.scatter(self.voronoi_vertices[:, 0], self.voronoi_vertices[:, 1], color="orange")
        vs = self.voronoi_vertices # alias
        # 4. plot voronoi edges
        for idxa, idxb in self.voronoi_edges:
            plt.plot([vs[idxa, 0], vs[idxb, 0]], [vs[idxa, 1], vs[idxb, 1]], color='red')
        margin = min(ymax - ymin, xmax - xmin) / 10
        # 5. plot voronoi rays
        for (xc, yc), (dx, dy) in self.voronoi_rays:
            if xmin <= xc <= xmax and ymin <= yc <= ymax:
                if abs(dx) > abs(dy):
                    xfar = xmax + margin if dx > 0 else xmin - margin
                    yfar = yc + (xfar - xc) / dx * dy
                else:
                    yfar = ymax + margin if dy > 0 else ymin - margin
                    xfar = xc + (yfar - yc) / dy * dx
                plt.plot([xc, xfar], [yc, yfar], "r--")


    def compute_voronoi_vertices(self):
        """compute voronoi vertices"""
        if not hasattr(self, 'voronoi_vertices'):
            if self.simplices_DT is None:
                self.create_delaunay()
            nfaces = self.simplices_DT.shape[0]
            self.voronoi_vertices = np.zeros((nfaces, 2))
            self.faces_dict = {}
            for i in range(nfaces):
                idx_a, idx_b, idx_c = self.simplices_DT[i, :]
                self.faces_dict[idx_a, idx_b, idx_c] = i
                A, B, C = self.vertices[idx_a, :], self.vertices[idx_b, :], self.vertices[idx_c, :]
                P, _ = compute_circumcircle(A, B, C)
                self.voronoi_vertices[i, :] = np.array(P)


    def _point_distance(self, P1, P2):
        """compute point distance"""
        x1, y1 = P1
        x2, y2 = P2
        return np.sqrt((x1 - x2)**2 + (y1 - y2)**2)


    def shape_procedure(self):
        """implement shape procedure"""
        xmin, xmax = round(self.vertices[:, 0].min()), round(self.vertices[:, 0].max())
        ymin, ymax = round(self.vertices[:, 1].min()), round(self.vertices[:, 1].max())
        xmin, xmax = xmin - 1, xmax + 1
        ymin, ymax = ymin - 1, ymax + 1
        plt.xlim([xmin, xmax])
        plt.ylim([ymin, ymax])
        self.compute_voronoi_vertices()
        vs = self.voronoi_vertices # alias
        hull = convex_hull(self.vertices, 1)
        self.voronoi_edges = []
        for edge, face in self.edge_faces_dict.items():
            face = list(face)
            p, q = self.get_vertex(edge[0]), self.get_vertex(edge[1])
            if len(face) == 2:
                face_idxa, face_idxb = self.faces_dict[face[0]], self.faces_dict[face[1]]
                a, b = vs[face_idxa, :], vs[face_idxb, :]
                m = 0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])
                if _is_left(p, q, a) * _is_left(p, q, b) == -1:
                    if 0.5 * self._point_distance(a, b) > self._point_distance(p, m):
                        plt.plot([p[0], q[0]], [p[1], q[1]], 'b-')
                    else:
                        plt.plot([a[0], b[0]], [a[1], b[1]], 'g-')
                else:
                    plt.plot([a[0], b[0]], [a[1], b[1]], 'g-')

            else:
                face_idx = self.faces_dict[face[0]]
                i = hull.index(p)
                j = (i + 1) % len(hull)
                if hull[j] != q:
                    p, q = q, p
                a = self.voronoi_vertices[face_idx, :]
                if _is_left(p, q, a) == 1:
                    plt.plot([p[0], q[0]], [p[1], q[1]], 'b-')


    def compute_voronoi(self):
        """compute voronoi diagram based on delaunay triangulation"""
        self.compute_voronoi_vertices()
        vs = self.voronoi_vertices # alias
        # compute voronoi edge
        hull = convex_hull(self.vertices, 1)
        self.voronoi_edges = []
        self.voronoi_rays = []
        for edge, face in self.edge_faces_dict.items():
            face = list(face)
            if len(face) == 2:
                face_idxa, face_idxb = self.faces_dict[face[0]], self.faces_dict[face[1]]
                self.voronoi_edges.append((face_idxa, face_idxb))
            else: # boundary face
                face_idx = self.faces_dict[face[0]]
                xc, yc = vs[face_idx, 0], vs[face_idx, 1]
                (x1, y1), (x2, y2) = self.get_edge(edge)
                i = hull.index((x1, y1))
                j = (i + 1) % len(hull)
                if hull[j] != (x2, y2):
                    (x1, y1), (x2, y2) = (x2, y2), (x1, y1)
                dx, dy = x2 - x1, y2 - y1
                h, v = -dy, dx
                intersect_point = compute_line_intersect((x1, y1), (x2, y2), (xc, yc), (xc + h, yc + v))
                isl =  _is_left((x1, y1), (x2, y2), (xc, yc))
                if isl == 1:
                    h = intersect_point[0] - xc
                    v = intersect_point[1] - yc
                elif isl == 0: # check this
                    h = -dy
                    v = dx
                else:
                    h = xc - intersect_point[0]
                    v = yc - intersect_point[1]
                L = np.sqrt(h**2 + v**2)
                cos = h / L
                sin = v / L
                self.voronoi_rays.append(((xc, yc), (cos, sin)))


    def edges_equal(self, edge1, edge2, tol=1e-12):
        (x1, y1), (x2, y2) = edge1
        (x3, y3), (x4, y4) = edge2
        return (abs(x1 - x3) <= tol and abs(y1 - y3) <= tol and abs(x2 - x4) <= tol and abs(y2 - y4) <= tol) \
                or (abs(x1 - x4) <= tol and abs(y1 - y4) <= tol and abs(x2 - x3) <= tol and abs(y2 - y3) <= tol)


    def crust(self):
        """implement crust algorithm"""
        self.compute_voronoi_vertices()
        plt.scatter(self.vertices[:, 0], self.vertices[:, 1], color="red")
        vertices = np.vstack((self.vertices, self.voronoi_vertices))
        if self.simplices_DT is None: # not yet delaunay
            self.create_delaunay()
        lst = []
        for i in range(vertices.shape[0]):
            x, y = vertices[i, :]
            flag = False
            for xx, yy in lst:
                if abs(xx - x) <= 1e-12 and abs(yy - y) <= 1e-12:
                    flag = True
                    break
            if not flag:
                lst.append([x, y])
        vertices = np.array(lst)
        temp_tri = Triangulation(vertices.copy())
        temp_tri.create_delaunay()
        # temp_tri.plot_delaunay()
        old_edges = set()
        for edge_indexes in self.edge_faces_dict.keys():
            edge = self.get_edge(edge_indexes)
            old_edges.add(edge)
        new_edges = set()
        for edge_indexes in temp_tri.edge_faces_dict.keys():
            edge = temp_tri.get_edge(edge_indexes)
            new_edges.add(edge)
        # edges = old_edges.intersection(new_edges)
        # self.crust_edges = list(edges)
        edges = []
        for old_edge in old_edges:
            for new_edge in new_edges:
                if self.edges_equal(old_edge, new_edge):
                    edges.append(old_edge)
                    break
        self.crust_edges = edges
        for (x1, y1), (x2, y2) in self.crust_edges:
            plt.plot([x1, x2], [y1, y2], 'b-')


    def ccw(self, idx1, idx2, idx3):
        A, B, C = self.vertices[idx1, :], self.vertices[idx2, :], self.vertices[idx3, :]
        if _is_left(A, B, C) == -1:
            return (idx1, idx3, idx2)
        else:
            return (idx1, idx2, idx3)


    def flip_and_update(self, tris, edge, idxp, idxq):
        """flip an edge and update"""
        # print(edge)
        idxr, idxs = edge
        self.edge_faces_dict.pop(edge, None)
        for tri in tris:
            if tri in self.tri_set:
                self.tri_set.remove(tri)
            for i in range(3):
                j = (i + 1) % 3
                idxi, idxj = tri[i], tri[j]
                edge_temp = self.order_edge(idxi, idxj)
                if edge_temp in self.edge_faces_dict and tri in self.edge_faces_dict[edge_temp]:
                    self.edge_faces_dict[edge_temp].remove(tri)
                    if len(self.edge_faces_dict[edge_temp]) == 0:
                        self.edge_faces_dict.pop(edge_temp, None)
        new_tri1 = self.ccw(idxr, idxp, idxq)
        new_tri2 = self.ccw(idxs, idxp, idxq)
        for tri in [new_tri1, new_tri2]:
            self.tri_set.add(tri)
            for i in range(3):
                j = (i + 1) % 3
                idxi, idxj = tri[i], tri[j]
                edge_temp = self.order_edge(idxi, idxj)
                self.edge_faces_dict[edge_temp].add(tri)


    def order_edge(self, idx0, idx1):
        return (idx0, idx1) if idx0 < idx1 else (idx1, idx0)


    def plot(self):
        plt.plot(self.vertices[:, 0], self.vertices[:, 1], 'o')


def test_DS():
    n = 40
    # np.random.seed(10)
    points = np.random.uniform(0, 5, (n, 2))
    tri = Triangulation(points)
    figure = plt.figure(figsize=(24, 8))
    plt.subplot(1, 3, 1)
    tri.plot_naive_triangulation()
    plt.subplot(1, 3, 2)
    tri.plot_delaunay()
    plt.subplot(1, 3, 3)
    tri.compute_voronoi()
    tri.plot_voronoi()
    # plt.subplot(1, 4, 4)
    vor = Voronoi(tri.vertices)
    voronoi_plot_2d(vor)
    # tri.crust()
    # tri.print_edges()
    # tri.print_faces()
    # tri.plot()
    plt.show()


def test_voronoi():
    n = 10
    # np.random.seed(10)
    points = np.random.uniform(0, 5, (n, 2))
    tri = Triangulation(points)
    figure = plt.figure(figsize=(8, 8))
    plt.subplot(1, 1, 1)
    tri.compute_voronoi()
    tri.plot_voronoi()
    vor = Voronoi(tri.vertices)
    voronoi_plot_2d(vor)
    # tri.crust()
    # tri.print_edges()
    # tri.print_faces()
    # tri.plot()
    plt.show()



def test_crust_and_shape():
    # filename = "testPoints_heart.txt"
    # filename = "test-file-1.txt"
    filename = "point-files/testPoints_circle.txt"
    # filename = "testPoints_random20.txt"
    # filename = "testPoints_astroid.txt"
    points = np.loadtxt(filename, delimiter=',')
    # filename = "test-file-1.txt"
    # filename = "test-file-2.txt"
    # points = np.loadtxt(filename)
    # n = 20
    # points = np.random.uniform(0, 5, (n, 2))
    tri = Triangulation(points)
    # tri.plot_naive_triangulation()
    # tri.create_delaunay()
    # tri.plot_delaunay()
    tri.crust()
    plt.axis('equal')
    plt.show()


def test_circumcircle():
    # A, B, C = np.array([0, 0]), np.array([1, 0]), np.array([1, 1])
    # D = np.array([0, 1.5])
    points = np.random.uniform(0, 5, (4, 2))
    A, B, C = points[0, :], points[1, :], points[2, :]
    D = points[3, :]
    _plot_circumcircle(A, B, C, D)
    print(_is_incircumcircle(A, B, C, D))
    plt.show()

if __name__ == "__main__":
    # test_DS()
    # test_circumcircle()
    # test_crust_and_shape()
    test_voronoi()
