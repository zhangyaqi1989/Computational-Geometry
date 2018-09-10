#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
the whole pipeline of shape processing
"""

# standard library
import sys
sys.path.insert(0, '../../lib/')

# third party library
import numpy as np
import matplotlib.pyplot as plt

# local modules
from marching_square import user_input, project_pipeline

def main():
    expression, limits, h, points_query = user_input()
    project_pipeline(expression, limits, h, points_query)

if __name__ == "__main__":
    main()
