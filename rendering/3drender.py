import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
from sys import argv, exit
import pandas as pd

plt.style.use('_mpl-gallery')

# Parse Command Line Arguments
if len(argv) < 2:
    exit("Usage is render.py <input file> <num parts>")
input_file = argv[1]
num_parts = int(argv[2])

savefreq = 1
nsteps = 1000
nsaves = int(nsteps / savefreq)

density = 0.0005
size = np.sqrt(density * num_parts);

# Read data
with open(input_file, 'r') as f:
    parts = np.loadtxt(input_file)

print(len(parts[:,0]))
print(len(parts[:,1]))
print(len(parts[:,2]))

t = np.array([np.ones(num_parts)*i for i in range(nsaves+1)]).flatten()
print(t)
print(len(t))
df = pd.DataFrame({"time": t ,"x" : parts[:,0], "y" : parts[:,1], "z" : parts[:,2]})

def update_graph(num):
    data=df[df['time']==num]
    graph._offsets3d = (data.x, data.y, data.z)
    title.set_text(str(num_parts) + ' particles, time={}ms'.format(num*10))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_zlim3d(0,size)
title = ax.set_title('NBody Simulation')

data=df[df['time']==0]
graph = ax.scatter(data.x, data.y, data.z)

ani = matplotlib.animation.FuncAnimation(fig, update_graph, nsaves+1, 
                               interval=20, blit=False)

# Save as mp4. This requires mplayer or ffmpeg to be installed
# anim.save('lorentz_attractor.mp4', fps=15, extra_args=['-vcodec', 'libx264'])


plt.show()