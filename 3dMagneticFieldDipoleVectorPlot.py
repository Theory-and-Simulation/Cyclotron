
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

#Kind of gives a 3d representation of dipole field using vectors

fig = plt.figure()
ax = fig.gca(projection='3d')

x, y, z = np.meshgrid(np.arange(-100, 100, 20),
                      np.arange(-100, 100, 20),
                      np.arange(-100, 100, 50))
#Meshgrid gives puts x, y, and z into a grid format
#Arange gives the range of values for x, y, and z, the first value being the minimum and the second value being the maximum, the third value is step size
#The minumum and maximum values are straightforward
#The stepsize effectively tells where the vector arrows are allowed to be plotted using this format... a vector arrow will show up for every .4 value in the z axis

r=np.linalg.norm([x,y,z])
u =3*x*z
v = 3*y*z
w = (2*(z**2)-(x**2)-(y**2))
print(r)
#Using the actual equation with r**-5 will effectively send u, v, and w to 0, thus killing the vectors in the plot, since r is large compared to x, y, and z, at least for the current range

ax.quiver(x, y, z, u, v, w, length=.00125)
#The first 3 arguments (in this case x, y, and z), give the bases
#The next 3 arguments (in this case u, v, and w), give the vector components (for x, y, and z respectively)
#Length is the length of the graphical vector arrow drawn in the plot

plt.show()