import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# The data
x = np.linspace(0, 10, 1000)
y = np.sin(2 * np.pi * x)

# The colormap
cmap = cm.jet

# Create figure and axes
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(1, 1, 1)

c = np.linspace(0, 10, 1000)
ax.scatter(x, y, c=c, cmap=cmap)