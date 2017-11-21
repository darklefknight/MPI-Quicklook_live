import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import numpy as np
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

high = 0.001
low = 0.003

x_max = .005

fig, ax12 = plt.subplots()

norm = Normalize(vmin=0,vmax=x_max)
cmap=cm.get_cmap("rainbow")


x_ticks = np.linspace(0,x_max,5)



ax12.plot((-1, low), (1, 1), lw=5, color="black")
ax12.plot((-1, high), (2, 2), lw=5, color="black")


ax12.scatter(low, 1, s=300, color=cmap(norm(low)), lw=3, zorder=10, edgecolors="black")
ax12.scatter(high, 2, s=300, color=cmap(norm(high)), lw=3, zorder=10, edgecolors="black")

ax12.set_xticklabels('')
ax12.set_xticks(x_ticks)
x_offset = x_max/8
ax12.set_xticks([x+x_offset for x in x_ticks], minor=True)
ax12.set_xticklabels("")


ax12.tick_params(axis="x", width=0, color="white",which="both")
ax12.tick_params(axis="y", width=0, color="white")

ax12.set_yticks([1, 2]
                )
ax12.set_yticklabels([])
ax12.grid(axis="x",ls="dashed")

ax12.set_xlim(0, x_max)
ax12.set_ylim(0, 3)

ax12.text(0.02, 0.38, 'Near Surface',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax12.transAxes,
        color='black', fontsize=15)

ax12.text(0.02, 0.72, 'Higher Atmosphere',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax12.transAxes,
        color='black', fontsize=15)

ax12.text(0.07, 0.05, 'Clear       Light Dust       Dusty    Very Dusty',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax12.transAxes,
        color='black', fontsize=15)
