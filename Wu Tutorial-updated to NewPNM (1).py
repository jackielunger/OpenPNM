
# coding: utf-8

## Using OpenPNM to recreate Wu's "Determination of oxygen effective diffusivity in porous gas diffusion layer using a three-dimensional pore network model"

# We are going to use OpenPNM to regenerate data found by Wu et al in the paper "Determination of oxygen effective diffusivity in porous gas diffusion layer using a three-dimensional pore network model".  This will serve as a introduction to OpenPNM, and instill an ease with the structure and use of OpenPNM for future simulations.  

### Getting started

# It is assumed that the reader has already downloaded or cloned the most recent version of OpenPNM available through github at https://github.com/PMEAL/OpenPNM.  To make importing OpenPNM as simple as possible, we can add OpenPNM to our pythonpath.  This can be done using the following line of code typed into the command line (replacing /path/to/OpenPNM with the proper path where OpenPNM is saved):  

# In[ ]:

export PYTHONPATH=$PYTHONPATH:/path/to/OpenPNM


# Now that we have added OpenPNM to the python path, we can import it.  We will also need matplotlib.pyplot if we wish to generate graphs of our data in IPython Notebook, and numpy if we wish to perform math operations on our data later (such as log()).

# In[1]:

cd OpenPNM


# In[2]:

import OpenPNM
import matplotlib.pyplot as plt
import numpy as np


### Generating the pore network and adding a geometry object

# The first step will be to create a pore network we wish to run simulations on.  To regenerate Wu's data, we need to create a pore network that is nxnx2n, where n can be 8, 10, 12, 14, 16, 18, or 20. Wu also specifies that the lattice parameter of the network should be 25e-6.  OpenPNM makes creating a pore network easy.  First, we select the desired network topology (in this case, cubic), then call the generate() method.  Note that divisions[ax,by,cz] supply the number of pores in each of the x, y, and z directions of the lattice.  

# In[3]:

n = 8
Lc = 25e-6

pn = OpenPNM.Network.Cubic(name = 'Wu')

#code to run if we want to set add_boundaries to be False
pn.generate(divisions = [n,n,2*n], add_boundaries= False, lattice_spacing = [Lc])

#code to run if we want to set add_boundaries to be True
#pn.generate(divisions = [n,n,2*n], add_boundaries= True, lattice_spacing = [Lc])


# OpenPNM makes it very easy to visualize the network we have generated through the "Visualization" methods.  We can create vtk files to be viewed using ParaView (downloadable at http://www.paraview.org/download/.  It is suggested that version 3.98 is downloaded instead of 4.1).  If we were able to visualize our pore network model it would appear like this:

# In[4]:

from IPython.display import display
from IPython.display import Image

i = Image(url = 'http://i.imgur.com/ILg7ZdJ.png')
display(i)


# Next, we must create a geometry object so that each pore and throat can be given properties.  For pores, this means diameter and volume, while for throats this means diameter, length, and volume.  We ensure that the geometry will be set for all locations of the network by sending all pore locations as a parameter.
# 
# Afer setting up the geometry object we are ready to add methods to our geometry object for calculating pore and throat dimensional values.  After we add all the methods, we must remember to "regenerate" so that these values are calculated and stored.  Note that if we don't want the logger to show us its progress, we can add 'loglevel = 30' to the parameters of each function.
# 
# The order in which we add each method is important, as the 'sphere' model for calculating pore volumes requires the diameters to be already known.  Similarly, throat lengths and diameters must be set before the throat volumes can be calculated.

# In[5]:

#code to run if boundaries was set to false
geo = OpenPNM.Geometry.GenericGeometry(name = 'wu_geometry', network = pn)
geo.set_locations(pores = pn.pores('all'), throats = 'all')

#code to run if boundaries was set to True
#pn.generate(divisions = [n,n,2*n], add_boundaries= True, lattice_spacing = [Lc], loglevel = 30)
#geo = OpenPNM.Geometry.GenericGeometry(name = 'wu_geometry', network = pn)
#geo.set_locations(pores = pn.pores('internal'), throats = 'all')
#boun = pn.add_geometry(subclass='Boundary',name='boun')
#boun.set_locations(pores=pn.pores('boundary'))

low = .5e-6
high = 9.5e-6

geo.add_method(prop='pore_diameter',model='constant', value = 24e-6)
geo.add_method(prop='pore_volume',model='sphere')
geo.add_method(prop='throat_diameter', model='uniform_distribution', low = low, high = high)
geo.add_method(prop='throat_length',model='straight')
geo.add_method(prop='throat_volume',model='cylinder')

pn.regenerate_geometries()


# Now we can use methods in Tools to return information about our pores and throats.  Note that the printed throat diameters are between 1e-6 and 1.9e-5 (twice our chosen minimum and maximum radii), so we know that the uniform_distribution method is working correctly.  We used the 'straight' model for calculating throat lengths, which simply calculated this length based on the pore diameter and lattice parameter.  Because the lattice parameter is 25e-6 and our pore diameters are 24e-6, our throat lengths should all be 1e-6 which we can check by printing these values.

# In[6]:

throat_diameters = pn.get_throat_data(prop = 'diameter') #if you do not specificy locations, all locations are returned
throat_lengths = pn.get_throat_data(prop = 'length')
print(throat_diameters)
print(throat_lengths)


### Adding fluid objects and methods

# Next, we have to set up our fluids.  We could set up air and water as generic fluids and add methods to each, or we could use the air and water fluids that are already exisiting.  We will use the already existing fluids to make our lives easier.  Again, we will use the regenerate method to make sure that the values for the fluids are calculated and set using the methods chosen (which in this case are preset).  

# In[7]:

air = OpenPNM.Fluids.Air(network = pn, name = 'air')
water = OpenPNM.Fluids.Water(network = pn, name = 'water')
pn.regenerate_fluids()


### Adding Physics objects and methods

# The next step will be to set up physics objects and add the proper methods.  However, Wu does not use the simple bulk_diffusion model we already have to calculate diffusive conductance.  The method we currently have calculates the diffusive conductance accross a conduit (half of one pore, the connecting throat, and half of the next pore) instead of just accross a throat.  We are assuming that Wu calculates the diffusive conductance simply accross a throat.    
# 
# Before we add methods to our physics objects, we should write a method for calculating diffusive conductance that follows Wu's model.  This will not be encorporated into OpenPNM, but can still be used to calculate our diffusive_conductance.  We will call this bulk_diffusion_wu, and it will appear as follows:

# In[3]:

import scipy as sp

def bulk_diffusion_wu(physics,
                      network,
                      fluid,
                      geometry,
                      propname,
                      diffusivity = 'diffusivity',
                      molar_density = 'molar_density',
                      throat_diameter = 'diameter',
                      throat_length = 'length',
                      pore_diameter = 'diameter',
                      **params):
    r"""
        Calculate the diffusive conductance of throats in network (instead of a
        conduit) based on the areas
        
        Parameters
        ----------
        network : OpenPNM Network Object
        
        fluid : OpenPNM Fluid Object
        The fluid of interest
        
        Notes
        -----
        This function requires that all the necessary fluid properties already be
        calculated.
        
        """
    #ct = fluid.get_data(prop='molar_density',throats='all',mode='interpolate')
    #Interpolate pore values to throats
    DABt = fluid.get_data(prop='diffusivity',throats='all',mode='interpolate')
    #Find g for full throat
    tdia = network.get_throat_data(prop=throat_diameter)
    tlen = network.get_throat_data(prop=throat_length)
    gt = (sp.pi*DABt*tdia**2)/(tlen*4)
    g = gt[geometry.throats()]
    #check for occupancy of connected pores
    try: 
        fluid['pore.occupancy']
        connected_pores = network.find_connected_pores(fluid.throats())
        s = []
        for item in connected_pores:
            s_item = False
            for pore in item:
                if fluid['pore.occupancy'][pore]:
                    s_item = True
            s.append(s_item)
        s=sp.array(s)
        g = g * s + g*(-s)/1e3
    except: pass    
    fluid.set_data(prop=propname,throats=geometry.throats(),data=g)


# The code within the try statement is important for after invasion percolation is run and some pores and throats are occupied by invading fluid.  If both neighboring pores of a throat are filled with water, then the conductance in the throat must be set to a much smaller value.  
# 
# Now we can set up our physics objects and add the methods with the desired models.  Then we need to "regenerate" to make sure these values are calculated.  Air's diffusive conductance is set a little differently, because we have written our own method that does not exist in OpenPNM.  Note that calling regenerate_physics() will not re-calculate air diffusive conductance.

# In[9]:

phys_water = OpenPNM.Physics.GenericPhysics(network=pn,fluid=water, geometry = geo, name='standard_water_physics')
phys_air = OpenPNM.Physics.GenericPhysics(network=pn,fluid=air, geometry = geo, name='standard_air_physics')

phys_water.add_method(prop='capillary_pressure', model='washburn') #accounts for cylindrical throats
phys_water.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_water.add_method(prop='diffusive_conductance', model='bulk_diffusion', shape = 'circular')
phys_air.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')

bulk_diffusion_wu(physics = phys_air, network = pn, fluid = air, geometry = geo, propname = 'diffusive_conductance')
pn.regenerate_physics()


### Running the invasion percolation algorithm

# Now we are ready to run invasion percolation.  We need to set up our parameters, namely inlets, outlets, and end_condition.  The end_condition parameter, 'breakthrough', specifies that invasion percolation will stop as soon as 1 of the outlet pores is filled.  If we specified 'total' instead, invasion percolation would proceed until all outlet pores are full.
# 
# We will save a vtk file with the saved data so we can view what we've done.  To view this, you can open the file inside ParaView.  

# In[10]:

inlets = pn.get_pore_indices(labels = ['bottom']) #can put in brackets so the whole bottom of the lattice is considered 1 inlet
outlets = pn.get_pore_indices(labels = ['top'])
end_condition = 'breakthrough'

OP_1 = OpenPNM.Algorithms.InvasionPercolation(network = pn, name = 'OP_1',loglevel=30)
OP_1.run(invading_fluid = water, defending_fluid = air, inlets = inlets, outlets = outlets, end_condition = end_condition)
OP_1.update()
#should be uncommented if want to run fickian diffusion on an air-filled lattice
#OP_1.update(IPseq = 0)  

vis = OpenPNM.Visualization.VTK()
vis.write(filename = 'test.vtp', network=pn,fluids=[air,water])


# If all has gone well, we should be able to watch the invasion percolation take place by opening our vtk file in ParaView.  After adding a threshold, we can watch an animation of the invasion that will appear as follows:

# In[2]:

from IPython.display import YouTubeVideo
display(YouTubeVideo('0iSuypRaT7A'))


### Running the Fickian Diffusion algorithm

# Next, we can set up our Fickian Diffusion algorithm.  Before we run it, we also need to set up our top and bottom boundary.  To continue following Wu's paper, we want to make the top boundary a plane 1/4 from the top, and the bottom boundary 1/4 from the bottom.  This way the Fickian Algorithm is only calculated on the center half.
# 
# diffusive conductance needs to be recalculated between when invasion percolation is run and effective diffusivity is calculated.  OpenPNM has been written to take care of this without our intervention.  However, because we are using a method for calculating effective diffusivity that is not part of OpenPNM, we must remember to recalculate diffusive conductance ourselves.

# In[12]:

bulk_diffusion_wu(physics = phys_air, network = pn, fluid = air, geometry = geo, propname = 'diffusive_conductance')
Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loglevel = 30, loggername = 'Fickian', name = 'fickian_alg', network = pn)

A = pn._Nx**2

z_dimension = int(pn.domain_size(dimension = 'height')/Lc) #number of pores in the z direction
quarter_layer = z_dimension/4 #estimates which layer marks 1/4 up the lattice
pore_number = int(quarter_layer*A) #gives the first pore in the layer 1/4 up the lattice

bottom_boundary = list(range(pore_number, pore_number + A))
top_boundary = list(range(pn.num_pores() - pore_number, pn.num_pores() - pore_number +A))

pn.set_pore_info(label='bound1', locations = bottom_boundary)
pn.set_pore_info(label='bound2', locations = top_boundary)


# To check that our boundaries have been set properly, we can view the network in ParaView.  The images of the boundaries should appear as below:

# In[13]:

#top boundary
i = Image(url = 'http://i.imgur.com/RknaFjl.png')
display(i)
#bottom boundary
g = Image(url = 'http://i.imgur.com/oxF407s.png')
display(g)


### Calculating effective diffusivity

# The next step will be to use the effective_diffusivity() method to calculate the effective diffusivity of the network after invasion is at completion.  Note that we do not need to run Fickian diffusion before we calculate effective diffusivity, because calling the effective diffusivity method runs the diffusion itself.  We must remember to input our boundaries as parameters, so that effective_diffusivity is only calculated accross the section we have chosen instead of over the whole network.

# In[14]:

effective_diffusivity = Fickian_alg.effective_diffusivity(fluid = air, boundary_1 = bottom_boundary, boundary_2 = top_boundary)
bulk_diffusivity = air.get_pore_data(prop = 'diffusivity')
normal_diffusivity = effective_diffusivity/bulk_diffusivity


# In[15]:

print(normal_diffusivity)


### Calculating saturation

# In Wu's first graph, he compares saturation and effective_diffusivity/bulk_diffusivity.  The only value that we are missing is saturation, which we can acquire with the following loop.  We only wish to find the saturation in the section that we used to calculate effective diffusion.  Therefore, we only check pores between the last pore in the bottom boundary and the first pore in the top boundary.  Note that the volume of throats is neglected.

# In[16]:

final_pores = water.get_pore_data('occupancy')*1
pore_volumes = pn.get_pore_data(prop = 'volume')

sum_volume = 0
filled_volume = 0

for i in range(bottom_boundary[-1], top_boundary[0]):
    sum_volume += pore_volumes[i]
    if final_pores[i] != 0:
        filled_volume += pore_volumes[i]

saturation = filled_volume/sum_volume 
print(saturation)


### Writing a method for running simulations

# If we want to run many simulations in a row, it makes more sense to put all the code generated so far into one function.  We can call this function run_simulation(...), and give it n, low, high, and end_condition as paramaters (remember that n helps with the size of the network, low and high give the min and max throat radii values, and end_condition determines when we want invasion percolation to terminate).

# In[5]:

def run_simulation(n, low, high, run_dry = False):
    
    Lc = 25e-6
    n = n
    low = low
    high = high
    
    #making a network that is [n, n, 2n] pores in dimension
    pn = OpenPNM.Network.Cubic(name = 'Wu', loglevel = 30)
    pn.generate(divisions = [n,n,2*n], add_boundaries= False, lattice_spacing = [Lc], loglevel = 30)
    locations = pn.get_pore_indices()
    geo = OpenPNM.Geometry.GenericGeometry(name = 'wu_geometry', network = pn)
    geo.set_locations(pores = pn.pores('all'), throats = 'all')
        
    #all methods must be added to geo (note 'even_distribution' for throat diameter)
    geo.add_method(prop='pore_diameter',model='constant', value = 24e-6)
    geo.add_method(prop='pore_volume',model='sphere')
    geo.add_method(prop='throat_diameter', model='uniform_distribution', low = low, high = high)
    geo.add_method(prop='throat_length',model='straight')
    geo.add_method(prop='throat_volume',model='cylinder')
    pn.regenerate_geometries()
    
    #Setting up fluids
    air = OpenPNM.Fluids.Air(network = pn, name = 'air')
    water = OpenPNM.Fluids.Water(network = pn, name = 'water')
    pn.regenerate_fluids()
    
    phys_water = OpenPNM.Physics.GenericPhysics(network=pn,fluid=water, geometry = geo, name='standard_water_physics')
    phys_air = OpenPNM.Physics.GenericPhysics(network=pn,fluid=air, geometry = geo, name='standard_air_physics')
    phys_water.add_method(prop='capillary_pressure', model='washburn')
    phys_water.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
    phys_water.add_method(prop='diffusive_conductance', model='bulk_diffusion', shape = 'circular')
    phys_air.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
    #phys_air.add_method(prop='diffusive_conductance',model='bulk_diffusion')
    bulk_diffusion_wu(physics = phys_air, network = pn, fluid = air, geometry = geo, propname = 'diffusive_conductance')
    phys_water.capillary_pressure()
    phys_water.hydraulic_conductance()
    phys_water.diffusive_conductance()
    phys_air.hydraulic_conductance()
    
    inlets = [pn.get_pore_indices(labels = ['bottom'])]
    outlets = pn.get_pore_indices(labels = ['top'])
   
    inlets = pn.get_pore_indices(labels = ['bottom']) #in brackets so the whole bottom of the lattice is considered 1 pore
    outlets = pn.get_pore_indices(labels = ['top'])
    end_condition = 'breakthrough'

    OP_1 = OpenPNM.Algorithms.InvasionPercolation(network = pn, name = 'OP_1',loglevel=30)
    OP_1.run(invading_fluid = water, defending_fluid = air, inlets = inlets, outlets = outlets, end_condition = end_condition)
    OP_1.update()
    
    if(run_dry == True):
        OP_1.update(IPseq = 0)
    Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian', name = 'fickian_alg', network = pn)
    
    #must be recalculated after invasion percolation runs
    bulk_diffusion_wu(physics = phys_air, network = pn, fluid = air, geometry = geo, propname = 'diffusive_conductance')
    
    #set labels for top boundary
    #set labels for bottom boundary
    A = pn._Nx**2

    z_dimension = int(pn.domain_size(dimension = 'height')/Lc) #number of pores in the z direction
    quarter_layer = z_dimension/4 #estimates which layer marks 1/4 up the lattice
    pore_number = int(quarter_layer*A) #gives the first pore in the layer 1/4 up the lattice

    bottom_boundary = list(range(pore_number, pore_number + A))
    top_boundary = list(range(pn.num_pores() - pore_number, pn.num_pores() - pore_number +A))
    
    #effective_diffusivity has been rewritten to use the correct length.  Length being used is now half! perfect!
    effective_diffusivity = Fickian_alg.effective_diffusivity(fluid = air, boundary_1 = bottom_boundary, boundary_2 = top_boundary)
    bulk_diffusivity = air.get_pore_data(prop = 'diffusivity')
    normal_diffusivity = effective_diffusivity/bulk_diffusivity
    
    #calculation of saturation
    final_pores = water.get_pore_data('occupancy')*1
    pore_volumes = pn.get_pore_data(prop = 'volume')

    sum_volume = 0
    filled_volume = 0

    for i in range(bottom_boundary[-1], top_boundary[0]):
        sum_volume += pore_volumes[i]
        if final_pores[i] != 0:
            filled_volume += pore_volumes[i]
        
    saturation = filled_volume/sum_volume 
    result = [saturation, normal_diffusivity[0]]
    return result


### Writing a method for saving simulation data

# It is helpful to save the results generated for each simulation inside a data file so that we can see the progress as the simulation runs.  We can write a saving_simulation_data() function for this.  At the same time, we can pass a variable to define the number of times we want to run a particular kind of simulation, and have saving_simulation_data() save data for multiple runs.  In addition to being written to a file, the data should also be saved into x_values and y_values arrays, to make it easy to graph later.

# In[6]:

def saving_simulation_data(n, low, high, number_times_run, x_values, y_values, run_dry = False):
    x_values = x_values
    y_values = y_values
    
    for x in range(number_times_run):
        data = run_simulation(n, low, high, run_dry = run_dry)
        x_values.append(data[0])
        y_values.append(data[1])
        
        if(x == 0):
            f = open('simulationResults', 'w')
        else:
            f = open('simulationResults', 'a')
        f.write('\nsimulation: ')
        f.write(str(x))
        f.write('\nsaturation: ')
        f.write(str(data[0]))
        f.write('\nnormalized diffusivity: ')
        f.write(str(data[1]))
        f.close()


### Generating a graph of saturation versus normalized diffusivity

# Generating Wu's first graph should now be trivial.  First, we must gather the data by running simulations with n = 14, low = .5e-6, high = 9.5e-6.  By running simulations with different end_conditions, we make sure that saturation varies enough to show the relationship between saturation and normalized diffusivity.  Wu suggests using end_condition values of 0, breakthrough (1 pore out of total outlets), .25, .5, and 1.  By running more than 1 simulation of each type, we make sure to generate enough data to be sure of the trend.  
# 
# Note that the more simulations that are run, the longer this process will take.  To adequately regenerate Wu's data, we would want to run each of these different end_condition simulations upwards of 20 times.  For the sake of saving time, we can run each only 2 times to generate a graph that shows the trend.  With number_times_run for each of the lines below at 2, the simulation will still take a few minutes.  The progress of the simulations can be watched inside the file that the data is being saved to.

# In[7]:

x_values1 = []
y_values1 = []

saving_simulation_data(n = 14, low = .5e-6, high = 9.5e-6, number_times_run = 20, x_values = x_values1, y_values = y_values1)


# the data now has x values saved to x_values and y values saved to y_values.  We've already imported matplotlib.pyplot, so only a few lines of code are needed to generate a graph of the data.  
# 
# One thing to note is that when using IPython notebook, installed with anaconda, using matplotlib will

# In[8]:

plt.plot(x_values1, y_values1, 'ro')
plt.title('saturation versus normalized diffusivity')
plt.xlabel('saturation')
plt.ylabel('normalized diffusivity')
plt.show()


# The graph generated thus far should look something like this:

# In[3]:

i = Image(url = 'http://i.imgur.com/CbuVrur.png')
display(i)


### Generating a graph of N (network size parameter) versus f(epsilon)

# Wu's second graph shows how N affects f(epsilon).  f(epsilon) is the same as normalized diffusivity when saturation is 0, so we can run many tests on dry networks. After running many trials, we will have x_values and y_values saved.  we want to use the values that N takes as our x_values, so we must also generate n_values to make graphing possible.
# 
# Again, processing this information will take longer if x is higher.  We can run simulations with x = 2 that will run much faster, but generate much less data.  Running with x = 5 will take just a few minutes at most.

# In[15]:

x_values = []
y_values = []
n_values = []

x = 5
saving_simulation_data(n = 8, low = .5e-6, high = 9.5e-6, run_dry = True, number_times_run = x, x_values = x_values, y_values = y_values)
for i in range(x):
    n_values.append(8)
saving_simulation_data(n = 10, low = .5e-6, high = 9.5e-6, run_dry = True, number_times_run = x, x_values = x_values, y_values = y_values)
for i in range(x):
    n_values.append(10)
saving_simulation_data(n = 12, low = .5e-6, high = 9.5e-6, run_dry = True, number_times_run = x, x_values = x_values, y_values = y_values)
for i in range(x):
    n_values.append(12)
saving_simulation_data(n = 14, low = .5e-6, high = 9.5e-6, run_dry = True, number_times_run = x, x_values = x_values, y_values = y_values)
for i in range(x):
    n_values.append(14)
saving_simulation_data(n = 16, low = .5e-6, high = 9.5e-6, run_dry = True, number_times_run = x, x_values = x_values, y_values = y_values)
for i in range(x):
    n_values.append(16)
saving_simulation_data(n = 20, low = .5e-6, high = 9.5e-6, run_dry = True, number_times_run = x, x_values = x_values, y_values = y_values)
for i in range(x):
    n_values.append(20)


# Now we plot n_values versus y_values, so that we get a graph of N versus f(epsilon).  we can use plt.axis(x_min, x_max, y_min, y_max) to control the axis on the graph, to make it more clear that f(epsilon) should not change with N.  The graph should also make clear that the standard deviation of f(epsilon) decreases as N increases.  This makes sense, as increasing N also increases the number of throats.
# 
# The plot axis have been adjusted to make it more obvious that F(epsilon) barely changes with N.  

# In[16]:

plt.plot(n_values, y_values, 'ro')
plt.title('N versus F(epsilon)')
plt.xlabel('N')
plt.ylabel('F(epsilon)')
plt.axis([6, 22, .06, .1])
plt.show()


# This graph should match the one below:

# In[19]:

i = Image(url = 'http://i.imgur.com/X0oKL3j.png')
display(i)


### Generating a graph of saturation versus g(s)

# Wu's third graph plots saturation versus g(s).  g(s)f(epsilon) = normalized_diffusivity, so g(s) = normalized_diffusivity/f(epsilon).  Because we are not varying the set up of our network, f(epsilon) will be constant.  Luckily, our second graph calculates this value many times for us.  We should use the average of this value in our calculation of g(s).  Lastly, we can graph g(s) using our x_values1 from our first graph, and our g(s) values calculated from our y_values1 and average f(epsilon).    

# In[28]:

#find average value for f(epsilon)
sum = 0
count = 0
for x in y_values:
    sum += x
    count += 1
average_f = sum/count

#prints graph for g(s) 

g_values = list(range(len(y_values1)))

for x in range(len(y_values1)):
    g_values[x] = y_values1[x]/average_f
plt.plot(x_values1, g_values, 'ro')
plt.title('saturation versus g(s)')
plt.xlabel('saturation')
plt.ylabel('g(s)')
plt.axis([0, .7, 0, 1])
plt.show()


# The graph generated should mimic the one below:

# In[29]:

i = Image(url = 'http://i.imgur.com/YMTmKQ1.png')
display(i)


### Generating a graph of saturation versus alpha

# Finally, we want to print a graph of saturation versus alpha.  This gives us an idea of the equation for g(s) (g(s) = (1-s)^alpha).  The following prints an alpha value for every data point we have gathered.  Wu et al find alpha by using a best fit curve, but this simplification is enough to show that we have alpha values close to those Wu generated.  

# In[31]:

alpha = []
for x in range(len(g_values)):
    alpha_value = np.log(g_values[x])/np.log(1 - x_values1[x])
    alpha.append(alpha_value)
plt.plot(x_values1, alpha, 'ro')
plt.title('alpha versus saturation')
plt.xlabel('saturation')
plt.ylabel('alpha')
#plt.axis([0, .5, 0, 3])
plt.show()


# The graph generated should look something like this:

# In[32]:

i = Image(url = 'http://i.imgur.com/6HHbcBP.png')
display(i)


### Discrepencies between our data and Wu's data

# Our values aren't quite the same as those generated by Wu et al.  there are several reasons for this-
# 
# 1) the effective diffusivity calculations done by OpenPNM assume a binary system.  Our code finds diffusion of oxygen through nitrogen, while Wu et al find diffusion for oxygen alone.  Being able to completely copy Wu's equations would require much editing to OpenPNM, and it has been decided that maintaining only a binary system way of calculating effective diffusivity is preferable.
# 
# 2) There are some assumptions that had to be made that we cannot be sure matched perfectly with Wu's assumptions.  One example of this is temperature.  We have assumed the oxygen Wu was running through the network was at 273 degrees Kelvin, but we cannot be certain that this is the value he was using.
# 
# 3) There are some slight calculation differences Wu uses because his throats are circular instead of square in cross-section and his pores are spherical instead of cubic.  This mostly effects fluid dynamics calculations including diffusive conductance, hydraulic conductance, and capillary pressure.
