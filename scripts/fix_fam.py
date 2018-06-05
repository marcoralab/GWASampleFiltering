#!/usr/bin/env python3

import sys

with open(sys.argv[1], "r") as input:
    with open(sys.argv[2], "w") as output:
        for line in input:
            row = line.split()
            split2 = row[1].split("_")
            row[0] = split2[0]
            row[1] = "_".join(split2[1:])
            print("\t".join(row), file=output)
