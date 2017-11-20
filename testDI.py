import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import numpy as np
import matplotlib.patches as mpatches

high = 0.001
low = 0.004

norm = Normalize(vmin=0,vmax=0.005)
cmap=cm.get_cmap("rainbow")

y_labels = ["Near Surface","Higher Atmosphere"]
x_ticks = np.linspace(0,0.005,4)
x_labels = ["Clear","Light Dust","Dusty","Very Dusty"]

fig, ax = plt.subplots()

ax.patches.Fac
ax.plot((-1,low),(1,1),lw=5,color="black")
ax.plot((-1,high),(2,2),lw=5,color="black")


ax.scatter(low,1,s=300,color=cmap(norm(low)),lw=3,zorder=10,edgecolors="black")
ax.scatter(high,2,s=300,color=cmap(norm(high)),lw=3,zorder=10,edgecolors="black")

ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels)

ax.tick_params(axis="y",width=0,color="white")

ax.set_yticks([1,2])
ax.set_yticklabels(y_labels)
ax.set_xlim(0,0.005)
ax.set_ylim(0,3)
