import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import pylab
import numpy as np
from numpy import ma
from mpl_toolkits.axes_grid1 import make_axes_locatable

x1, x2, y1, y2 = 40.0, 40.40, 1.30, 1.70
width = x2-x1
normal_depth = 0.00798
output_step = 2.0
end_time = 12+output_step
Lx = 42.0
output_series = np.arange(0, end_time, step=2, dtype=int)
profile_series = np.arange(0, end_time, step=2, dtype=int)
for i in output_series:
    filename = "out-h-"+str(i)
    picname = "g1-"+str(i)+".png"
    m = np.loadtxt(fname=filename, delimiter=" ")
    m1 = m[:, 0]
    m2 = m[:, 1]
    m3 = m[:, 2]
    fig, ax = plt.subplots()
    ax.set_aspect(1)
    # change the figure size
    plt.rcParams['figure.figsize'] = [40, 21]
    # set the font globally
    plt.rcParams.update({'font.family': 'Tex Gyre Pagella','font.size':45 })
    # plt.xlim([20.0, 40.0])
    # plt.ylim([0.0, 3.0])
    triang = tri.Triangulation(m1, m2)
    tcf = ax.tricontourf(triang, m3, levels=20)
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="50%", pad=0.2)
    #fig.colorbar(tcf, cax=cax)
    #fig.colorbar(tcf)
    plt.colorbar(tcf,fraction=0.01, pad=0.02)
    # ax.tricontour(triang, m3, colors = 'k')
    #ax.tricontour(triang, m3)
    rectangle = plt.Rectangle((x1, y1), width, width, facecolor="k", fill=True)
    ax.add_patch(rectangle)
    plt.savefig(picname, bbox_inches='tight')

for i in profile_series:
    filename = "profile-"+str(i);
    picname = "profile1-"+str(i)+".png";
    m = np.loadtxt(fname=filename, delimiter=" ");
    m1 = m[:,0];
    m2 = m[:,1];
    m11 = m1[(m2<10.0) & ((m1>x1)|(m1<x2))];
    m22 = m2[(m2<10.0) & ((m1>x1)|(m1<x2))];
    fig, ax = plt.subplots()
    plt.xlim(0.0, Lx)
    plt.ylim(0.0, 1.15*np.max(m22[0:-2]))
    _ = ax.plot(m11, m22, 'b-o', markersize=15)
    plt.savefig(picname)

# zoomed-in contour for the block
for i in output_series:
    filename = "out-h-"+str(i)
    picname = "g2-"+str(i)+".png"
    m = np.loadtxt(fname=filename, delimiter=" ")
    m1 = m[:, 0]
    m2 = m[:, 1]
    m3 = m[:, 2]
    fig, ax = plt.subplots()
    # change the figure size
    plt.rcParams['figure.figsize'] = [15, 12]
    # set the font globally
    plt.rcParams.update({'font.family': 'Tex Gyre Pagella','font.size':17 })
    plt.xlim([39.0, 42.0])
    plt.ylim([0.0, 3.0])
    triang = tri.Triangulation(m1, m2)
    tcf = ax.tricontourf(triang, m3, levels=20)
    fig.colorbar(tcf)
    # ax.tricontour(triang, m3, colors = 'k')
    #ax.tricontour(triang, m3)
    rectangle = plt.Rectangle((x1, y1), width, width, facecolor="k", fill=True)
    ax.add_patch(rectangle)
    plt.savefig(picname)
