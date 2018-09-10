#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
implements some algorithms about polygon
"""

import numpy as np
import matplotlib.pyplot as plt

def _ccw(P1, P2, P3):
    return (P2[0] - P1[0]) * (P3[1] - P1[1]) - (P2[1] - P1[1]) * (P3[0] - P1[0])


def _monotone_chain(points):
    if type(points) == np.ndarray:
        xs = list(points[:, 0])
        ys = list(points[:, 1])
        points = list(zip(xs, ys))
    upper = []
    lower = []
    points = sorted(set(points))
    for point in points:
        while len(upper) > 1 and _ccw(upper[-2], upper[-1], point) >= 0:
            upper.pop()
        while len(lower) > 1 and _ccw(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        upper.append(point)
        lower.append(point)
    return upper, lower


def rotating_calipers(points):
    upper, lower = _monotone_chain(points)
    i = 0
    j = len(lower) - 1
    pairs = []
    while i < len(upper) - 1 or j > 0:
        pairs.append((upper[i], lower[j]))
        if i == len(upper) - 1:
            j -= 1
        elif j == 0:
            i += 1

        elif (upper[i+1][1] - upper[i][1]) * (lower[j][0] - lower[j-1][0]) > \
                (lower[j][1] - lower[j-1][1]) * (upper[i+1][0] - upper[i][0]):
            i += 1
        else:
            j -= 1
    return pairs

if __name__ == "__main__":
    xmin, xmax = 0, 5
    ymin, ymax = 0, 5
    # np.random.seed(1)
    n = 20
    xs = np.random.uniform(xmin, xmax, n)
    ys = np.random.uniform(ymin, ymax, n)
    points = list(zip(list(xs), list(ys)))
    naive_pairs = [(points[i], points[j]) for i in range(n) for j in range(n) if j > i]
    # print(naive_pairs)
    naive_distances = [(x0 - x1)**2 + (y0 - y1)**2 for ((x0, y0), (x1, y1)) in naive_pairs]
    print("max distance is {:0.6f}".format(max(naive_distances)))
    upper, lower = _monotone_chain(points[:])
    # print(len(upper))
    # print(len(lower))
    pairs = rotating_calipers(points[:])
    D, pair = max([(((p[0] - q[0])**2 + (p[1] - q[1])**2), (p, q)) for \
            p, q in pairs])
    for point in points:
        x, y = point
        plt.scatter(x, y, color='g')
    '''
    for point in upper:
        x, y = point
        plt.scatter(x, y, color='b')
    for point in lower:
        x, y = point
        plt.scatter(x, y, color='r')
    for pair in pairs:
        (x0, y0), (x1, y1) = pair
        plt.plot([x0, x1], [y0, y1], 'b--')
    '''
    (x0, y0), (x1, y1) = pair
    dx, dy = x1 - x0, y1 - y0
    dx, dy = -dy, dx
    plt.plot([x0 - dx, x0, x0 + dx], [y0 - dy, y0, y0 + dy], '--')
    plt.plot([x1 - dx, x1, x1 + dx], [y1 - dy, y1, y1 + dy], '--')
    plt.plot([x0, x1], [y0, y1], '-')
    print('Diameter is {0:f}'.format(D))
    distances = [(x0 - x1)**2 + (y0 - y1)**2 for ((x0, y0), (x1, y1)) in pairs]
    # print(distances)
    plt.axis('equal')
    hull = lower[:]
    xs, ys = zip(*hull)
    plt.plot(xs, ys, 'b--')
    hull = upper[:]
    xs, ys = zip(*hull)
    plt.plot(xs, ys, 'b--')
    plt.show()
