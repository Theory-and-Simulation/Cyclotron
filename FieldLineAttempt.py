# -*- coding: utf-8 -*-
import numpy as np
from scipy import special


#Fair warning, I had trouble getting mlab and mayavi components to install correctly, so I used 3rd party software called Enthought Canpoy that included the whole thing


#Here are the relevant equations for the magnitude of the magnetic field of a dipole in cartesian coordinates

#Bx =((np.sqrt(x**2+y**2+z**2))**-5)*3*x*z
#By = ((np.sqrt(x**2+y**2+z**2))**-5)*3*y*z
#Bz = ((np.sqrt(x**2+y**2+z**2))**-5)*(2*(z**2)-(x**2)-(y**2))

#Magnitude=np.sqrt((Bx**2)+(By**2)+(Bz**2))

#All this is according to this paper, https://arxiv.org/pdf/1112.3487.pdf
#which could've sworn was in our repository already (I might have missed it)
#An additional source corroborates this https://ccmc.gsfc.nasa.gov/RoR_WWW/presentations/Dipole.pdf
#This additional source expresses the terms differently, specifically using M instead of B0 and Re, but
#should be equivalent

#Now adjusting for the Earth's dipole specifically
#Bx =-((np.sqrt(x**2+y**2+z**2))**-5)*B0*(Re**3)*3*x*z
#By = -((np.sqrt(x**2+y**2+z**2))**-5)*B0*(Re**3)*3*y*z
#Bz = -((np.sqrt(x**2+y**2+z**2))**-5)*B0*(Re**3)*(2*(z**2)-(x**2)-(y**2))

#Where B0 is the measured field strength at the magnetic equator, B0=3.07*(10**-5)T
#And where Re is the radius of the Earth, Re=6.3781*(10**6)m
#With these adjustments made, the expression for magnitude is still the same: np.sqrt((Bx**2)+(By**2)+(Bz**2))






#Below is what I've managed so far with creating a graphical model, I'll upload pictures of the model as well


#### Calculate the field ####################################################

x, y, z = [e.astype(np.float32) for e in
            np.ogrid[-3.7:3.7:10j, -3.7:3.7:10j, -3.7:3.7:10j]]

B0=3.07*(10**-5)
Re=6.3781*(10**6)


# These are the equations for the x,y,z components of the magnetic field
#Important bit: We'll take Re=1 for now, which will thus set our scale of plot to multiples of Re
#I might have made an error in this simplification, I'll try and double check this later
#Upon double checking I believe the problem is with how x, y, and z are defined... I will expand more on this at the end
Bx =-((np.sqrt(x**2+y**2+z**2))**-5)*B0*3*x*z
By = -((np.sqrt(x**2+y**2+z**2))**-5)*B0*3*y*z
Bz = -((np.sqrt(x**2+y**2+z**2))**-5)*B0*(2*(z**2)-(x**2)-(y**2))
#Norming these components to get the magnetic field magnitude
np.sqrt((Bx**2)+(By**2)+(Bz**2))




# Free memory early
del x, y


del z







#### Visualize the field ####################################################
from mayavi import mlab
fig = mlab.figure(1, size=(400, 400), bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
#I believe the values for the above are unimportant

field = mlab.pipeline.vector_field(Bx, By, Bz)


# Unfortunately, the above call makes a copy of the arrays, so we delete
# this copy to free memory.
del Bx, By, Bz

# Create a sphere
r = 1 #Radius is 1, or 1 Re, as noted before
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0:pi:51j, 0:2 * pi:51j]

x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 400))
mlab.clf()

#Placing and centering the Earth
mlab.mesh(x+5.5, y+5.5, z+5.5, colormap='ocean')

#The entire above section starting from the "Create a sphere" comment is making a representation of the earth


magnitude = mlab.pipeline.extract_vector_norm(field)



field_lines = mlab.pipeline.streamline(magnitude, seedtype='line',
                                        integration_direction='both',
                                        seed_resolution=9,
                                        colormap='rainbow',
                                        #vmin=0, vmax=1
                                        )
#Integration direction must be kept to both, setting it to "forward" or "backward" will only give half the field
#I commented out vmin and vmax, as it rescaled the field line strength to a max of 1 and a min of 0                                        
#Increasing Seed resolution gives larger numbers of field lines                                       
                                        
             
#Plotting axes and using an additional outline to help visualize the axes                                         
mlab.axes(field, extent=(0,10, 0,10, 0,10))
mlab.outline(field, extent=(0,10, 0,10, 0,10))

#Various tweaks here
field_lines.stream_tracer.maximum_propagation = 100.
field_lines.seed.widget.point1 = [0, 0, 0]
field_lines.seed.widget.point2 = [0, 0, 0]
field_lines.seed.widget.resolution = 100
field_lines.seed.widget.enabled = False
#Deleting the first line of the above messes up the plot, the rest have to do with a fancy widget that's been set to not display for convenience's sake

#Note: When in the mayavi viewer, click on the upper left button that says "View the Mayavi pipeline"; head to colors and legends for the streamline tab, and click on show legend
#This will give legend for the colorcoding of the vector field strength
#Can also adjust how axes look
#Adjust number of labels, font size
#Also checking the streamline section and going to the StreamTracer tab gives integration options which can smooth the field lines 
#I've set the initial integration step to 0.1, integration step unit to 1, and total steps to 20,000 in the pictures



#The magnitudes displayed for the color map apparently depend on the ogrid set up for x, y, and z at the very beginning of this. I've tweaked the values such that they're kind of
#reasonable, but obviously I would like to do plot this properly. The actual equations themselves are correct, yet for some reason x, y, z = [e.astype(np.float32) for e in
#            np.ogrid[-3.7:3.7:10j, -3.7:3.7:10j, -3.7:3.7:10j]] this section is messing things up






mlab.show()