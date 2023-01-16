import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
from sys import argv, exit
import pandas as pd

  
plt.style.use('_mpl-gallery')

# Parse Command Line Arguments
if len(argv) < 1:
    exit("Usage is render.py <input file> <num parts>")
input_file = argv[1]

parts = np.loadtxt(input_file)

# Getting info
num_parts = int(parts[0][0])
size = parts[0][1]
nsteps = parts[0][2]
savefreq = 1
nsaves = int(nsteps / savefreq)

t = np.array([np.ones(num_parts)*i for i in range(nsaves+1)]).flatten()
df = pd.DataFrame({"time": t ,"x" : parts[1:,0], "y" : parts[1:,1], "z" : parts[1:,2]})

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

anim = matplotlib.animation.FuncAnimation(fig, update_graph, nsaves+1, 
                            interval=10, blit=False)

writervideo = matplotlib.animation.FFMpegWriter(fps=60)
anim.save('prova.mp4', writer=writervideo)
plt.close()

# Save as mp4. This requires mplayer or ffmpeg to be installed (NOT WORKING)
# ani.save('lorentz_attractor.mp4', fps=60, extra_args=['-vcodec', 'libx264'])

