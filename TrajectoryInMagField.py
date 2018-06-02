import numpy as np
import math as math 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from decimal import Decimal
from mayavi import mlab

speed_of_light = 3.0E08
E_field = [5000000.0, 0.0, 0.0] 
B0=3.07*(10**-5)
Re=6.3781*(10**6)
INITPOSIT=[2.0, 3.0, 1.0]
class Particle: 
	def __init__(self, position, velocity, mass, charge): 
		self.position = position
		self.velocity = velocity
		self.mass = mass
		self.charge = charge
	
	def radius(self,b): #input magnetic field as argument in array
		b = np.linalg.norm(b)
		speed = np.linalg.norm(self.velocity)
		radius = speed * self.mass / (self.charge * b)
		return radius
		
	def period(self):
		field = np.linalg.norm(B_field)
		period = 2.0 * math.pi * self.mass /( field * self.charge )
		return period
		


B_field = [((np.sqrt(INITPOSIT[0]**2+INITPOSIT[1]**2+INITPOSIT[2]**2))**-5)*B0*3*INITPOSIT[0]*INITPOSIT[2], -((np.sqrt(INITPOSIT[0]**2+INITPOSIT[1]**2+INITPOSIT[2]**2))**-5)*B0*3*INITPOSIT[1]*INITPOSIT[2], -((np.sqrt(INITPOSIT[0]**2+INITPOSIT[1]**2+INITPOSIT[2]**2))**-5)*B0*(2*(INITPOSIT[2]**2)-(INITPOSIT[0]**2)-(INITPOSIT[1]**2))]

# choose a proton as the particle, describe in 3 Dimension 
proton = Particle([5.0, 0.0, 1.00], [-2.0, 6.0, -1.0], 1.67E-27, +1.60E-19)



###################################################################################
#	Given position and velocity, the function returns acceleration as an array in 3D
#	The magnetic field is applied where x > 0 or  x < - gap_size
#	The electric field is applied where x < 0 and x > - gap_size
###################################################################################

jumps =0
jumps_max = 46
i = 0

def acceleration( charge_to_mass_ratio, r, v ):
	global i, jumps # if you use the global statement, the variable will become available 
	a=charge_to_mass_ratio*np.cross(v, B_field)

	return a

###############################################################################
#  
#    Solving Equations of Motions Using RK4
#
###############################################################################


# initializes array that will hold the velocity of particle 
# each time it gains velocity (crosses gap)
velo = np.zeros((jumps_max + 1,))
velo[0] = np.linalg.norm(proton.velocity)


def rk4(particle, iterations, desired_value):
	#initialize arrays to be appended by rk4 results
	RK4_pos = []
	RK4_vel = []

	
	n = 100
	h = particle.period() / n
	charge_to_mass_ratio = particle.charge / particle.mass
	
	ini_pos = np.array(particle.position) 
	ini_vel = np.array(particle.velocity) 
	
	for i in range(iterations + 10):
		i += 1
		p_i = ini_pos
		v_i = ini_vel
	
		k1 = h * (1.60E-19/1.67E-27)*np.cross(v_i,[((np.sqrt(p_i[0]**2+p_i[1]**2+p_i[2]**2))**-5)*B0*3*p_i[0]*p_i[2], -((np.sqrt(p_i[0]**2+p_i[1]**2+p_i[2]**2))**-5)*B0*3*p_i[1]*p_i[2], -((np.sqrt(p_i[0]**2+p_i[1]**2+p_i[2]**2))**-5)*B0*(2*(p_i[2]**2)-(p_i[0]**2)-(p_i[1]**2))])
		l1 = h * v_i
		
		
		k2 = h * (1.60E-19/1.67E-27)*np.cross(v_i+k1/2.0,[((np.sqrt((np.array(p_i)+np.array(l1)/2.0)[0]**2+(p_i+l1/2.0)[1]**2+(np.array(p_i)+np.array(l1)/2.0)[2]**2))**-5)*B0*3*(np.array(p_i)+np.array(l1)/2.0)[0]*(np.array(p_i)+np.array(l1)/2.0)[2], -((np.sqrt((np.array(p_i)+np.array(l1)/2.0)[0]**2+(np.array(p_i)+np.array(l1)/2.0)[1]**2+(np.array(p_i)+np.array(l1)/2.0)[2]**2))**-5)*B0*3*(np.array(p_i)+np.array(l1)/2.0)[1]*(np.array(p_i)+np.array(l1)/2.0)[2], -((np.sqrt((np.array(p_i)+np.array(l1)/2.0)[0]**2+(np.array(p_i)+np.array(l1)/2.0)[1]**2+(np.array(p_i)+np.array(l1)/2.0)[2]**2))**-5)*B0*(2*((np.array(p_i)+np.array(l1)/2.0)[2]**2)-((np.array(p_i)+np.array(l1)/2.0)[0]**2)-((np.array(p_i)+np.array(l1)/2.0)[1]**2))])
		l2 = h * (v_i + (k1 * 0.5))
		

		k3 = h * (1.60E-19/1.67E-27)*np.cross(v_i+k2/2.0,[((np.sqrt((np.array(p_i)+np.array(l2)/2.0)[0]**2+(np.array(p_i)+np.array(l2)/2.0)[1]**2+(np.array(p_i)+np.array(l2)/2.0)[2]**2))**-5)*B0*3*(np.array(p_i)+np.array(l2)/2.0)[0]*(np.array(p_i)+np.array(l2)/2.0)[2], -((np.sqrt((np.array(p_i)+np.array(l2)/2.0)[0]**2+(np.array(p_i)+np.array(l2)/2.0)[1]**2+(np.array(p_i)+np.array(l2)/2.0)[2]**2))**-5)*B0*3*(np.array(p_i)+np.array(l2)/2.0)[1]*(np.array(p_i)+np.array(l2)/2.0)[2], -((np.sqrt((np.array(p_i)+np.array(l2)/2.0)[0]**2+(np.array(p_i)+np.array(l2)/2.0)[1]**2+(np.array(p_i)+np.array(l2)/2.0)[2]**2))**-5)*B0*(2*((np.array(p_i)+np.array(l2)/2.0)[2]**2)-((np.array(p_i)+np.array(l2)/2.0)[0]**2)-((np.array(p_i)+np.array(l2)/2.0)[1]**2))])
		l3 = h * (v_i + (k2 * 0.5))
	
		k4 = h * (1.60E-19/1.67E-27)*np.cross(v_i+k3,[((np.sqrt((np.array(p_i)+np.array(l3))[0]**2+(np.array(p_i)+np.array(l3))[1]**2+(np.array(p_i)+np.array(l3))[2]**2))**-5)*B0*3*(np.array(p_i)+np.array(l3))[0]*(np.array(p_i)+np.array(l3))[2], -((np.sqrt((np.array(p_i)+np.array(l3))[0]**2+(np.array(p_i)+np.array(l3))[1]**2+(np.array(p_i)+np.array(l3))[2]**2))**-5)*B0*3*(np.array(p_i)+np.array(l3))[1]*(np.array(p_i)+np.array(l3))[2], -((np.sqrt((np.array(p_i)+np.array(l3))[0]**2+(np.array(p_i)+np.array(l3))[1]**2+(np.array(p_i)+np.array(l3))[2]**2))**-5)*B0*(2*((np.array(p_i)+np.array(l3))[2]**2)-((np.array(p_i)+np.array(l3))[0]**2)-((np.array(p_i)+np.array(l3))[1]**2))])
		l4 = h * (v_i + k3)
		

		ini_vel = ini_vel + (k1 + 2.0 * (k2 + k3) + k4) / 6.0
		ini_pos = ini_pos + (l1 + 2.0 * (l2 + l3) + l4) / 6.0
		vel_mag = np.linalg.norm(ini_vel)
		RK4_vel.append(vel_mag)
		RK4_pos.append(ini_pos)
			
	#Controls which array function returns, resets jump count
	if desired_value == 'velocity': 
		value = RK4_vel
	elif desired_value == 'position': 
		value = RK4_pos
	
	return value



#######################################################################################
#  
#    Since this program graphs many figures, it is more convenient and efficient
#		to create a plot function. 
#		This function can accomodate 2D and 3D plots, note the first argument
#		must be a 3 element list, whose entries are labels for x,y,z coordinates, respectively.
#
#########################################################################################

def plot(coord_labels, title, dimension, x, y, z = 0):
	
	fig = plt.figure()
	
	if dimension == '2d':
		ax = fig.add_subplot(1,1,1)
		ax.plot(x, y )
	elif dimension == '3d':
		ax = fig.add_subplot(111,projection = dimension)
		ax.plot3D(x, z, y)
		ax.set_zlabel(coord_labels[2])
		
	ax.set_xlabel(coord_labels[0])
	ax.set_ylabel(coord_labels[1])
	ax.set_title(title)
	
		


###############################################################################
#  
#    Plotting Trajectory of Particle
#
###############################################################################


def position_plot(particle):

	x = []
	y = []
	z = []
	
	x.append(particle.position[0])
	y.append(particle.position[1])
	z.append(particle.position[2])
	
	results = rk4(proton , 30000 ,'position')	
	
	#extracts x,y,z coordinates
	for i in results:
		x.append(i[0])
		y.append(i[1])
		z.append(i[2])
		
	coordlabel = ["X-axis (m)","Y-axis (m)","Z-axis (m)"]
	plot(coordlabel, "Particle Trajectory in the Earth's Magnetic Field", '3d', x,y,z)

	
#Note:axes are switched up

###############################################################################
#  
#    velocity	vs. time plot
#
###############################################################################


def velocity_plot(particle):
	n = 4
	h = particle.period() / n
	
	v = rk4(proton , 500000 ,'velocity')	
	t = np.linspace(0, len(v) * h, len(v))
	
	
	coordlabel = ["Time (s)","Velocity (m/s)","Z-axis (m)"]
	plot(coordlabel, "Beam Velocity as a Function of Time", '2d', t,v)
	

	#prints final velocity
	final_velo = velo[len(velo) - 1]
	print final_velo

	#prints final time
	final_t = t[len(velo) - 1]
	print final_t





	


jumps = 0
position_plot(proton)

plt.show()


