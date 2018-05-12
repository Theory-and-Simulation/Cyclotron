
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

#Kind of gives a 3d representation of dipole field using vectors

fig = plt.figure()
ax = fig.gca(projection='3d')

x, y, z = np.meshgrid(np.arange(-10, 10, 2),
                      np.arange(-10, 10, 2),
                      np.arange(-10, 10, 5))
#Meshgrid gives puts x, y, and z into a grid format
#Arange gives the range of values for x, y, and z, the first value being the minimum and the second value being the maximum, the third value is step size
#The minumum and maximum values are straightforward
#The stepsize effectively tells where the vector arrows are allowed to be plotted using this format... a vector arrow will show up for every .4 value in the z axis





#for x in range(-10, 11):
#	if x%10==0:
#		print(x)
#	for y in range(-10, 11):
#		if y%10==0:
#			print(y)
#	for z in range(-10,11):
#		if z%10==0:
#			print(z)

r=np.linalg.norm([x,y,z])
u =((np.sqrt(x**2+y**2+z**2))**-5)*3*x*z
v = ((np.sqrt(x**2+y**2+z**2))**-5)*3*y*z
w = ((np.sqrt(x**2+y**2+z**2))**-5)*(2*(z**2)-(x**2)-(y**2))
#Using the actual equation with r**-5 will effectively send u, v, and w to 0, thus killing the vectors in the plot, since r is large compared to x, y, and z, at least for the current range

ax.quiver(x, y, z, u, v, w, length=400)
#The first 3 arguments (in this case x, y, and z), give the bases
#The next 3 arguments (in this case u, v, and w), give the vector components (for x, y, and z respectively)
#Length is the length of the graphical vector arrow drawn in the plot

#Distance vector r is only being checked once, for loop is probably needed
#Actually wait, just express r in terms of x,y,z in the equations for u,v,w

#The plot looks much uglier, but the relative size of the vectors look correct, so this should give a generally alright view of the field magnitude

plt.show()