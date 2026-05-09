#!/usr/bin/env python3
import sys

infile_name = sys.argv[1]
outfile_name = sys.argv[2]


with open(infile_name, 'r') as infile, open(outfile_name, 'w') as outfile:
    for line in infile:
        print('chr' + line, end='', file=outfile)
