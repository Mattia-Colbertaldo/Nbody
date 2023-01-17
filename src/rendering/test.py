import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import cnames
from matplotlib import animation
from sys import argv, exit

plt.style.use('_mpl-gallery')

# Parse Command Line Arguments
if len(argv) < 2:
    exit("Usage is render.py <input file> <num parts>")
input_file = argv[1]
num_parts = int(argv[2])

# Read data
with open(input_file, 'r') as f:
    parts = np.loadtxt(input_file)
    print(parts[:])

    # Set up figure & 3D axis for animation
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1], projection='3d')
    ax.axis('off')

    # choose a different color for each trajectory
    colors = plt.cm.jet(np.linspace(0, 1, num_parts))

    # set up lines and points
    lines = sum([ax.plot([], [], [], '-', c=c)
                for c in colors], [])
    pts = sum([ax.plot([], [], [], 'o', c=c)
            for c in colors], [])

    # prepare the axes limits
    ax.set_xlim((-25, 25))
    ax.set_ylim((-35, 35))
    ax.set_zlim((5, 55))

    # set point-of-view: specified by (altitude degrees, azimuth degrees)
    ax.view_init(30, 0)

    # initialization function: plot the background of each frame
    def init():
        for line, pt in zip(lines, pts):
            line.set_data([], [])
            line.set_3d_properties([])

            pt.set_data([], [])
            pt.set_3d_properties([])
        return lines + pts

    # animation function.  This will be called sequentially with the frame number
    def animate(i):
        print(i)
        # we'll step two time-steps per frame.  This leads to nice results.
        i = (2 * i) % parts.shape[1]

        for line, pt, part in zip(lines, pts, parts):
            print(part[:i])
            x, y, z = part[:i].T
            line.set_data(x, y)
            line.set_3d_properties(z)

            pt.set_data(x[-1:], y[-1:])
            pt.set_3d_properties(z[-1:])

        ax.view_init(30, 0.3 * i)
        fig.canvas.draw()
        return lines + pts

    # instantiate the animator.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=1000, interval=30, blit=True)

    # Save as mp4. This requires mplayer or ffmpeg to be installed
    # anim.save('lorentz_attractor.mp4', fps=15, extra_args=['-vcodec', 'libx264'])
    plt.show()