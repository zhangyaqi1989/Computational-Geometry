#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
use skimage skeleton to approximate medial axis
"""

import numpy as np
from skimage.morphology import skeletonize
from skimage import data
import matplotlib.pyplot as plt
from skimage.util import invert


def f3(x, y):
    return x**2/4 + y**2 - 9

def f2(x, y, a=1):
    return (x**2 - a**2)*(x - a)**2 + (y**2 - a**2)**2


if __name__ == "__main__":
    # xmin, xmax = -7, 7
    # ymin, ymax = -4, 4
    xmin, xmax = -2, 2
    ymin, ymax = -2, 2
    nx, ny = 500, 500
    xs = np.linspace(xmin, xmax, nx)
    ys = np.linspace(ymin, ymax, ny)
    X, Y = np.meshgrid(xs, ys)
    # f_vec = np.vectorize(f3)
    f_vec = np.vectorize(f2)
    Z = f_vec(X, Y)
    Z[Z <= 0] = 0
    Z[Z > 0] = 1
    image = invert(Z)
    # image = invert(data.horse())

    # perform skeletonization
    skeleton = skeletonize(image)

    # display results
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4),
                             sharex=True, sharey=True)

    ax = axes.ravel()

    ax[0].imshow(image, cmap=plt.cm.gray)
    ax[0].axis('off')
    ax[0].set_title('original', fontsize=20)

    ax[1].imshow(skeleton, cmap=plt.cm.gray)
    ax[1].axis('off')
    ax[1].set_title('skeleton', fontsize=20)

    fig.tight_layout()
    plt.show()
