import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from itertools import groupby
from sys import argv, exit

plt.style.use('_mpl-gallery')

# Make data
np.random.seed(19680801)
n = 100
rng = np.random.default_rng()
xs = rng.uniform(23, 32, n)
ys = rng.uniform(0, 100, n)
zs = rng.uniform(-50, -25, n)
###

# Parse Command Line Arguments
if len(argv) < 3:
    exit("Usage is render.py <input file> <output file> [cutoff]")
input_file = argv[1]
output_file = argv[2]

# Read data
with open(input_file, 'r') as f:
    lines = []
    for line in f:
        lines.append([float(i) for i in line.rstrip().split(" ")])
    print(lines)

    # Plot
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.scatter(xs, ys, zs)

    ax.set(xticklabels=[],
        yticklabels=[],
        zticklabels=[])

    plt.show()