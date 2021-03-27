import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
output_step = 2.0;
end_time = 40+output_step;
output_series = np.arange(0, end_time, step=2, dtype=np.int);
for i in output_series:
    filename = "out-h-"+str(i);
    picname = "g1-"+str(i)+".png";
    m = np.loadtxt(fname=filename, delimiter=" ");
    m1 = m[:,0];
    m2 = m[:,1];
    m3 = m[:,2];
    fig, ax = plt.subplots()
    # plt.xlim([20.0, 40.0])
    # plt.ylim([0.0, 3.0])
    triang = tri.Triangulation(m1, m2)
    tcf = ax.tricontourf(triang, m3)
    fig.colorbar(tcf)
    # ax.tricontour(triang, m3, colors = 'k')
    ax.tricontour(triang, m3)
    plt.savefig(picname)
