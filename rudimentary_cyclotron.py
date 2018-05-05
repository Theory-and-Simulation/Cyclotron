import numpy as np
import math as math
import matplotlib.pyplot as plt
from decimal import Decimal

final = 4
iteration = 10
c = 3.0e8
gap = 1

class Particle:
	
	#def __init__(self, ms, chrg, elec, mag, pos, vel):
	def __init__(self, ms, chrg, elec, mag):
		self.pMass = ms
		self.pCharge = chrg
		self.bField = np.array(mag)
		self.eField = np.array(elec)

	"""accelertaion function:
	@parameter: rPos is the position of the particle
	@parameter: vPos is the velocity of the particle
	returns: the acceleration at the given coordinate given the properties of the particle
	"""
	def acceleration(self, rPos, vPos):
		global iteration
		q_by_m = (self.pCharge/self.pMass)
		if iteration == final:
			#if we have reached the final iteration, then we have turned off all the acceleration
			return 0

		elif 0>=rPos[0] or rPos[0]>=gap:
			#if the particle is in magnetic field, and we simply give the acceleration due to magnetic field
			return q_by_m*np.array(np.cross(vPos, self.bField))

		else:
			#otherwise the particle has the coordinates between 0 and 1, there is electric field
			if rPos[1]<0:
				#if the particle is below the origin, then the acceleration is positive
				return q_by_m*self.eField
			else:
				#if the particle is below the origin, then the acceleration is negative
				return -q_by_m*self.eField/self.pMass

e = [5000000.0,0.,0.]
b = [0.,0.,-1.5]
mass = 1.67e-27
charge = 1.6e-19

proton = Particle(mass, charge, e, b)

# setting the upper and lower limits of the partition a <= t <= b

a = 0.0
b = 10*2*math.pi * proton.pMass / (np.linalg.norm(b) * proton.pCharge)
 
 
# Setting the number of partitions and calculating the width of each partion
n = 4000
h = (b-a)/n

#
t = np.arange( a, b+h, h )
r = np.zeros((n+1,3))
v = np.zeros((n+1,3))

t[0]= 0.0
r[0]= [0.0, 0.0, 0.0]
v[0]= [0.05*c, 0.0, 0.0]


"""*
function derivR is the derivative dr/dt as a function of time, r and v
@parameter tVal is the time
@parameter rVal is the position vector
@parameter vVal is the velocity vector
@returns the value of derivative at the given point
"""
def deriR(tVal, rVal, vVal):
	return vVal

"""*
function derivV is the derivative dv/dt as a function of time, r and v
@parameter tVal is the time
@parameter rVal is the position vector
@parameter vVal is the velocity vector
@returns the value of derivative at the given point
"""
def deriV(t, r, v):
	return proton.acceleration( r, v)

"""*
function rungeKutta approtimates r on the interval [t(k-1), t(k)]
@parameter tVal is the time point at t(k-1)
@parameter rVal is the position vector at t(k-1)
@parameter vVal is the velocity vector at t(k-1)
@parameter h = t(k) - t(k-1)
@returns an arrar with two vectors. first vector is approtimation of 
"""
def rungeKutta(tVal, rVal, vVal, h):
	#calculates k1, k2, k3, k4, l1, l2, l3 and l4

	k1 = h*deriR(tVal, rVal, vVal)
	l1 = h*deriV(tVal, rVal, vVal)

	k2 = h*deriR(tVal + h/2, rVal + k1/2, vVal + l1/2)
	l2 = h*deriV(tVal + h/2, rVal + k1/2, vVal + l1/2)

	k3 = h*deriR(tVal + h/2, rVal + k2/2, vVal + l2/2)
	l3 = h*deriV(tVal + h/2, rVal + k2/2, vVal + l2/2)

	k4 = h*deriR(tVal + h, rVal + k3, vVal + l3)
	l4 = h*deriV(tVal + h, rVal + k3, vVal + l3)
	
	return [(k1+ 2*(k2 + k3) + k4)/6, (l1+ 2*(l2 + l3) + l4)/6]


"""
for loop runs to update values of r and v
over the different partitions of T = [a, b]
with step h using the rungeKutta method
"""
for i in range(1, n+1):
	increase = rungeKutta(t[i-1], r[i-1], v[i-1], h)
	r[i] = r[i-1] + increase[0]
	v[i] = v[i-1] + increase[1]


"""
function func returns the value of the actual function at a given value t
"""
# def func( t ):
#     return np.exp(t)
x = np.zeros((n+1,))
y = np.zeros((n+1,))

for i in range (0, n+1):
	x[i] = r[i][0]
	y[i] = r[i][1]




plt.plot(x, y, label='RK4')
plt.xlabel('x') 
plt.ylabel('y')
plt.legend(loc=2)
plt.grid()


plt.show()
