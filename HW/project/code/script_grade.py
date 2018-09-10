#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
a grade script for the project
"""

# standard library
import sys
import subprocess

# third party library
# import sh
# import numpy as np
# import matplotlib.pyplot as plt

code_commands_dic = {'m1' : ('medial_axis.py', 'poly18.txt'),
                     'm2' : ('medial_axis.py', 'poly_snake.txt'),
                     't1' : ('triangulation.py', 'poly18.txt'),
                     't2' : ('triangulation.py', 'poly_snake.txt'),
                     'p1' : ('pmc.py', 'polygon_1.txt', 'points_1.txt'),
                     'p2' : ('pmc.py', 'polygon_2.txt', 'points_2.txt'),
                     'r1' : ('reconstruction.py', 'expression_1.txt'),
                     'r2' : ('reconstruction.py', 'expression_3.txt'),
                     's1' : ('shape_processing.py', 'expression_1.txt', 'points_1.txt'),
                     's2' : ('shape_processing.py', 'expression_2.txt', 'points_1.txt'),
                    }


def main():
    assert(len(sys.argv) == 2)
    code = sys.argv[1]
    if code in code_commands_dic:
        cmd = ['python']
        cmd.extend(code_commands_dic[code])
        output = subprocess.check_output(cmd, timeout=10)
        if output:
            print(output.decode(), end='')
        # output = sh.python(*code_commands_dic[code])
        # contents = output.stdout.decode()
        # print(contents, end='')


if __name__ == "__main__":
    main()
