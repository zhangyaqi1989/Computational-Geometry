#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
create input points file of different closed curves
"""

import numpy as np
import matplotlib.pyplot as plt

def heart1(ts):
    """the first heart on
       http://mathworld.wolfram.com/HeartCurve.html
    """
    rs = 1 - np.sin(ts)
    xs = rs * np.cos(ts)
    ys = rs * np.sin(ts)
    return (xs, ys)


def heart2(ts):
    """the 6th heart on
       http://mathworld.wolfram.com/HeartCurve.html
    """
    xs = 16 * (np.sin(ts) ** 3)
    ys = 13 * np.cos(ts) - 5 * np.cos(2 * ts) - 2 * np.cos(3 * ts) - np.cos(4 * ts)
    return (xs, ys)


def bifolium(ts, a=1):
    """http://mathworld.wolfram.com/Bifolium.html"""
    rs = 4 * a * np.sin(ts) ** 2 * np.cos(ts)
    xs = rs * np.cos(ts)
    ys = rs * np.sin(ts)
    return (xs, ys)


def cayley_sextic(ts, a=1):
    """http://mathworld.wolfram.com/CayleysSextic.html"""
    rs = 4 * a * np.cos(ts / 3) ** 3
    xs = rs * np.cos(ts)
    ys = rs * np.sin(ts)
    return (xs, ys)


def cochleoid(ts, a=1):
    """http://mathworld.wolfram.com/Cochleoid.html"""
    rs = a * np.sin(ts) / ts
    xs = rs * np.cos(ts)
    ys = rs * np.sin(ts)
    return (xs, ys)


def write_coords_to_file(filename, xs, ys):
    """write coordinates to file"""
    n = min(len(xs), len(ys))
    with open(filename, 'w') as out_file:
        for i in range(n):
            out_file.write("{:f}\t{:f}\n".format(xs[i], ys[i]))
    print("Write {:d} points to {:s}".format(n, filename))


def generate_coords(ts, func):
    return func(ts)


if __name__ == "__main__":
    ts = np.arange(0, np.pi * 2, 0.1)
    # xs, ys = generate_coords(ts, heart2)
    # xs, ys = generate_coords(ts, bifolium)
    # xs, ys = generate_coords(ts, cayley_sextic)
    xs, ys = generate_coords(ts, cochleoid)
    # write_coords_to_file("bifolium.txt", xs, ys)
    # write_coords_to_file("cayley_sextic.txt", xs, ys)
    write_coords_to_file("cochleoid.txt", xs, ys)
    """
    plt.scatter(xs, ys)
    plt.axis('equal')
    plt.show()
    """
