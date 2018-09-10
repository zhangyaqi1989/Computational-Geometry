#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
# format all txt files to
# x1 y1
# ...
# xn yn
# delimiter = "\t"
##################################

import os
import re

def format_inputfile(filename):
    """format inputfile and save to new file"""
    lines = []
    with open(filename, 'r') as in_file:
        for line in in_file.readlines():
            if not line.strip(): # empty
                continue
            else:
                tokens = re.split(",\s* |\s*", line.strip())
                lines.append("\t".join(tokens))
    name, ext = filename.rsplit('.', 1)
    outfile_name = name + "-new.txt"
    with open(outfile_name, 'w') as out_file:
        for line in lines:
            out_file.write(line + "\n")

def main():
    """format all the txt files in the current directory"""
    files = [item for item in os.listdir('.') if item.endswith('.txt')]
    for item in files:
        format_inputfile(item)


if __name__ == "__main__":
    main()
