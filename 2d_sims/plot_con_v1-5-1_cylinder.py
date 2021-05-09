import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import pylab
import numpy as np
from numpy import ma
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True

x_center, y_center, cylinder_radius = 20.20, 1.50, 0.20
normal_depth = 0.00798
output_step = 2.0
end_time = 36+output_step
Lx = 22.0
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
    ax.set_yticks(np.arange(1.50,(4.50+1),1))
    ax.set_aspect(1)
    # change the figure size
    plt.rcParams['figure.figsize'] = [40, 30]
    # set the font globally
    plt.rcParams.update({'font.family': 'Tex Gyre Pagella','font.size':42 })
    ax.set_xlabel(r'{$x$ (\rm{m})}')
    ax.set_ylabel(r'{$y$ (\rm{m})}')
    plt.xlim([0.0, Lx])
    plt.ylim([0.0, 6.0])
    triang = tri.Triangulation(m1, m2)
    levels = np.arange(start=0.0, stop=0.077, step=0.007)
    #tcf = ax.tricontourf(triang, m3, levels=levels, vmin=0.0, vmax=0.077, cmap='bwr')
    tcf = ax.tricontourf(triang, m3, levels=20, vmin=0.0, cmap='bwr')
    #tcf.set_clim(0.0, 0.70)
    divider = make_axes_locatable(ax)
    cax = divider.new_vertical(size='22.5%', pad=1.0)
    fig.add_axes(cax)
    #cmap = mpl.cm.bwr
    #bounds = np.arange(start=0.0, stop=0.070, step=0.007)
    #norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')
    fig.colorbar(tcf, cax=cax, ax=ax, orientation='horizontal')
    #fig.colorbar(tcf, cax=cax)
    #fig.colorbar(tcf)
    #plt.colorbar(tcf,fraction=0.01, pad=0.02, orientation='horizontal')
    # ax.tricontour(triang, m3, colors = 'k')
    #ax.tricontour(triang, m3)
    cylinder = plt.Circle((x_center, y_center), cylinder_radius, facecolor="k", fill=True)
    ax.add_patch(cylinder)
    plt.savefig(picname, bbox_inches='tight')

#for i in profile_series:
    #filename = "profile-"+str(i);
    #picname = "profile1-"+str(i)+".png";
    #m = np.loadtxt(fname=filename, delimiter=" ");
    #m1 = m[:,0];
    #m2 = m[:,1];
    #m11 = m1[(m2<10.0) & ((m1>x1)|(m1<x2))];
    #m22 = m2[(m2<10.0) & ((m1>x1)|(m1<x2))];
    #fig, ax = plt.subplots()
    #plt.xlim(0.0, Lx)
    #plt.ylim(0.0, 1.15*np.max(m22[0:-2]))
    #_ = ax.plot(m11, m22, 'b-o', markersize=15)
    #plt.savefig(picname)

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
    plt.rcParams.update({'font.family': 'Tex Gyre Pagella','font.size':25})
    plt.rc('axes', labelsize=32) 
    ax.set_xlabel(r'{$x$ (\rm{m})}')
    ax.set_ylabel(r'{$y$ (\rm{m})}')
    plt.xlim([4.0, 7.0])
    plt.ylim([1.50, 4.50])
    triang = tri.Triangulation(m1, m2)
    tcf = ax.tricontourf(triang, m3, levels=20, cmap='bwr')
    fig.colorbar(tcf)
    # ax.tricontour(triang, m3, colors = 'k')
    #ax.tricontour(triang, m3)
    cylinder = plt.Circle((x_center, y_center), cylinder_radius, facecolor="k", fill=True)
    ax.add_patch(cylinder)
    plt.savefig(picname)
