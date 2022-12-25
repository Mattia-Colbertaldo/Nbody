#!/usr/bin/env python3
import numpy as np
from sys import argv, exit
import os

if __name__ == "__main__":
    if len(argv) < 2:
        exit("Usage is check.py <file_1> <file_2>")
    input_file1 = argv[1]
    input_file2 = argv[2]

    parts1 = np.loadtxt(input_file1)
    parts2 = np.loadtxt(input_file2)

    num_parts = int(parts1[0][0])

    for r in range(1, num_parts):
        for c in range(0,3):
            if (parts1[r][c] != parts2[r][c]):
                exit("Not equal...  (ノಠ ∩ಠ)ノ彡( \o°o)\    [line " + str(r) + "]")
    print("Equal! (╯°□°）╯︵ ┻━┻ ")