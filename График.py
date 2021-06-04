import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

fig = plt.figure()
ax = Axes3D(fig)
X = np.arange(-2,2,0.1)
Y = np.arange(-2,2,0.1)
X,Y = np.meshgrid(X,Y)

def f(x,y):
    # return (1-x) ** 2 + 100 * (y - x ** 2) ** 2
    return 2 * (x ** 2) + 2 * y ** 2 - 2 * y * x - 4 * x - 6 * y

surf = ax.plot_surface(X,Y,f(X,Y), rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()