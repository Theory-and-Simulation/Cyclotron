# -*- coding: utf-8 -*-
import numpy as np
from scipy import special

#I've made a number of improvements on the code so far. The most obvious is that the field lines are 3d now, and more closely resemble what you would see out of a
#representation of a dipole. Other fixes include changing the way x, y, and z were defined. mgrid is now used instead of ogrid. While I'm not entirely sure how
#this solved the problem of the vector magnitudes being incorrect, it does anyways. With our current settings, we can see that the field at the north pole (assuming
#the geographical north pole and the magnetic north pole are the same) is equal to 6.14E-5 T, which is exactly what we would expect given our equations. This makes me
#much more confident in using this graphical representation as a model. The field lines cut off a bit unevenly, but overall it seems fine. And most importantly the field
#magnitudes make sense




#Fair warning, I had trouble getting mlab and mayavi components to install correctly, so I used 3rd party software called Enthought Canopy that included the whole thing


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

a = -4
b = 4
ds = 1

x, y, z = np.mgrid[a:b:ds, a:b:ds, a:b:ds];
#This is working better than ogrid... a and b seem to not change anything in the magnitude, only the position of our field in the plot
#ds=1 is great
#ds=.1 multiplies the field mag by 1000
#ds=10 puts everything to something along the lines of 6.68E-8, which from our largest magnitude at ds=1, seems to indicate that multiplying ds by 10 divided the
#field mag by 1000
#Looking at our equation for the magnetic field, I believe this makes some amount of sense, at least for the field at the north pole
#The magnetic field at the north pole is equal to 6.14E-5. but with ds=0.1, it is set to 6.14E-2... 
#If one looks at the magnetic field Bz, (which is the only component at the north pole), we see that the our singular distance component, z, is to the order of -3
#This I believe, is the connection between ds and the changing of the vector magnitudes... thus we ought to be able to safely set ds=1, as that is the "true" setting for our
#values... everything is scaled to 1
#ds also changes the position of our field in the plot, but I'm not too worried about that since we're setting it equal to 1 anyways

B0=3.07*(10**-5)
Re=6.3781*(10**6)


# These are the equations for the x,y,z components of the magnetic field
#Important bit: We'll take Re=1 for now, which will thus set our scale of plot to multiples of Re
Bx =-((np.sqrt(x**2+y**2+z**2))**-5)*B0*3*x*z
By = -((np.sqrt(x**2+y**2+z**2))**-5)*B0*3*y*z
Bz = -((np.sqrt(x**2+y**2+z**2))**-5)*B0*(2*(z**2)-(x**2)-(y**2))
#Norming these components to get the magnetic field magnitude
np.sqrt((Bx**2)+(By**2)+(Bz**2))




# Free memory early
del x, y, z







#### Visualize the field ####################################################
from mayavi import mlab
fig = mlab.figure(1, size=(400, 400), bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
#The values for the are some optional settings for the figure, like the initial size of the window, and the color

field = mlab.pipeline.vector_field(Bx, By, Bz)


# Unfortunately, the above call makes a copy of the arrays, so we delete
# this copy to free memory.
del Bx, By, Bz

# Create a sphere
r = 1 #Radius is 1, or 1 Re, as noted before
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0:pi:180j, 0:2 * pi:180j]

x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 400))
mlab.clf()

#Placing and centering the Earth
mlab.mesh(x+5, y+5, z+5, colormap='ocean')

#The entire above section starting from the "Create a sphere" comment is making a representation of the earth


magnitude = mlab.pipeline.extract_vector_norm(field)


field_lines = mlab.pipeline.streamline(magnitude, seedtype='sphere',
                                        seed_scale=1,
                                        integration_direction='both',
                                        seed_resolution=10,
                                        colormap='rainbow',
                                        #vmin=0, vmax=1
                                        )
#Seedtype is important, we set it to sphere so that we can see the lines flowing outward and back in through the sphere                                       
#Integration direction must be kept to both, setting it to "forward" or "backward" will only give half the field
#I commented out vmin and vmax, as it rescaled the field line strength to a max of 1 and a min of 0                                        
#Increasing Seed resolution gives larger numbers of field lines, but it's better to override them imo using the streamline->seed tab in the in viewer tool                              
                                        
field_lines.seed.widget.center = [5, 5, 5]    #This places our seed widget, else it would be in some other (nonoptimal spot, though we can move it within the viewer)   
field_lines.seed.widget.radius = 1 #Initial size of seed widget


#Adjusting various integration options here, so that our lines are relatively smooth          
field_lines.stream_tracer.initial_integration_step = 0.1     
field_lines.stream_tracer.integration_step_unit = 1  
field_lines.stream_tracer.maximum_number_of_steps = 25000    

      
                  
#Plotting axes and using an additional outline to help visualize the axes                                        
mlab.axes(field, extent=(0,10, 0,10, 0,10), nb_labels=11)
mlab.outline(field, extent=(0,10, 0,10, 0,10))

#I believe propagation here affects how far the lines can extend, while I'm not sure about the previous statement, I am sure that you don't want to set this to 0 or something like that
field_lines.stream_tracer.maximum_propagation = 100.

#Note: When in the mayavi viewer, click on the upper left button that says "View the Mayavi pipeline"; head to colors and legends for the streamline tab, and click on show legend
#This will give a legend for the colorcoding of the vector field strength
#Can also adjust how axes look
#Adjust number of labels, font size, etc.
#Also checking the streamline section and going to the StreamTracer tab gives integration options which can smooth the field lines
#The size of the seed widget can also be set, which will extend or contract the field lines...

#So the viewer is kinda screwy, you have to finangle it a little...

###Settings adjusted for the pictures###
#The radius of the seed sphere = 1.88801633371
#Phi and Theta resolution has been set to 8 and 16 respectively






mlab.show()