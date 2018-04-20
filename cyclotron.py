###################################################################################
#	Cyclotron Accelerator:
#		Particle is subjected to constant magnetic field within the dee's that serves
#		to keep particle in semi-circular motion. Electric field is applied within
#		gaps which accelerates particle, increasing its energy.
#	Lorentz Force Equation:
#		Motion is described by F = q [ E + v x B]
###################################################################################


import numpy as np
import math as math 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from decimal import Decimal


speed_of_light = 3.0E08
E_field = [5000000.0, 0.0, 0.0] 
B_field = [0.0, 0.0, -1.5]
bmag = np.linalg.norm(B_field)

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
		



# choose a proton as the particle, describe in 3 Dimension 
proton = Particle([0.00, 0.0, 0.0], [0.05*speed_of_light, 0.0, 0.0], 1.67E-27, +1.60E-19)

# size of gap between deed
gap_size = .5 * proton.radius(B_field)


###################################################################################
#	Given position and velocity, the function returns acceleration as an array in 3D
#	The magnetic field is applied where x > 0 or  x < - gap_size
#	The electric field is applied where x < 0 and x > - gap_size
###################################################################################

jumps =0
jumps_max = 46
i = 0

def acceleration( charge_to_mass_ratio, r, v ):
	global i, jumps                # if you use the global statement, the variable will become available 
	if jumps >= jumps_max:			# "outside" the scope of the function, effectively becoming a global variable
		# No acceleration
		a = 0.0
	elif r[0] >= 0 or r[0] <= -gap_size:
		a = np.cross(v, B_field) 
		a = a * charge_to_mass_ratio
		if i:  #Executes every time count is nonzero, occurs after particle enters dee from gap
			jumps += 1
			velo[jumps] = np.linalg.norm(v) #saves new particle speed after gap cross
			i = 0
	else: 	
		a = np.array(E_field)	
		a = a * charge_to_mass_ratio
		if r[1] > 0:
			a = -a
		i += 1
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
	
	n = 400
	h = particle.period() / n
	charge_to_mass_ratio = particle.charge / particle.mass
	
	ini_pos = np.array(particle.position) 
	ini_vel = np.array(particle.velocity) 
	
	for i in range(iterations + 10):
		i += 1
		p_i = ini_pos
		v_i = ini_vel
	
		k1 = h * acceleration( charge_to_mass_ratio, p_i, v_i )
		l1 = h * v_i
		
		
		k2 = h * acceleration( charge_to_mass_ratio, ini_pos + (l1 * 0.5), ini_vel + (k1 * 0.5) )
		l2 = h * (ini_vel + (k1 * 0.5))
		

		k3 = h * acceleration( charge_to_mass_ratio, ini_pos + (l2 * 0.5), ini_vel + (k2 * 0.5) )
		l3 = h * (ini_vel + (k2 * 0.5))
	
		k4 = h * acceleration( charge_to_mass_ratio, ini_pos + l3, ini_vel + k3 )
		l4 = h * (ini_vel + k3)
		

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
	
	results = rk4(proton , 10000 ,'position')	
	
	#extracts x,y,z coordinates
	for i in results:
		x.append(i[0])
		y.append(i[1])
		z.append(i[2])
		
	coordlabel = ["X-axis (m)","Y-axis (m)","Z-axis (m)"]
	plot(coordlabel, "Particle Trajectory in Cyclotron Accelerator", '3d', x,y,z)

	

###############################################################################
#  
#    velocity	vs. time plot
#
###############################################################################


def velocity_plot(particle):
	n = 400
	h = particle.period() / n
	
	v = rk4(proton , 10000 ,'velocity')	
	t = np.linspace(0, len(v) * h, len(v))
	
	
	coordlabel = ["Time (s)","Velocity (m/s)","Z-axis (m)"]
	plot(coordlabel, "Beam Velocity as a Function of Time", '2d', t,v)
	

	#prints final velocity
	final_velo = velo[len(velo) - 1]
	print final_velo

	#prints final time
	final_t = t[len(velo) - 1]
	print final_t

###############################################################################
#  
#    velocity vs radius
#
###############################################################################


def vel_rad_plot(particle):


	radius = np.zeros((len(velo),))

	for item in range(1, len(velo)+1):
		radius[item-1] = velo[item-1] * particle.mass / (particle.charge * bmag)
		
	
	coordlabel = ["Radius (m)","Velocity (m/s)","Z-axis (m)"]
	plot(coordlabel, "Beam Velocity at Various Extraction Radius", '2d', radius,velo)
	


velocity_plot(proton) 
jumps = 0
position_plot(proton)
vel_rad_plot(proton)
plt.show()



