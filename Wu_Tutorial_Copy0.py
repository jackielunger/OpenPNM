# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>


import OpenPNM
import matplotlib.pyplot as plt
import numpy as np

# <codecell>

n = 8
Lc = 25e-6

pn = OpenPNM.Network.Cubic(name = 'Wu', loglevel = 30)
pn.generate(divisions = [n,n,2*n], add_boundaries= True, lattice_spacing = [Lc], loglevel = 30)

# <codecell>

locations = pn.get_pore_indices()
geo = OpenPNM.Geometry.GenericGeometry(name = 'wu_geometry', network = pn)
geo.set_locations(pores = pn.pores('all'), throats = 'all')

# <codecell>

low = 4e-6
high = 6e-6

geo.add_method(prop='pore_diameter',model='constant', value = 24e-6)
geo.add_method(prop='pore_volume',model='sphere')
geo.add_method(prop='throat_diameter', model='uniform_distribution', low = low, high = high)
geo.add_method(prop='throat_length',model='straight')
geo.add_method(prop='throat_volume',model='cylinder')

#geo.regenerate()
pn.regenerate_geometries()

# <codecell>

#pn['pore.volume'][pn['pore.boundary']]=1e-20
#pn['throat.length'][pn['throat.boundary']]=1e-20

# <codecell>

throat_diameters = pn.get_throat_data(prop = 'diameter') #if you do not specificy locations, all locations are returned
throat_lengths = pn.get_throat_data(prop = 'length')
pore_diameters = pn.get_pore_data(prop = 'diameter')
pore_volumes = pn.get_pore_data(prop = 'volume')
throat_volumes = pn.get_throat_data(prop = 'volume')
print(pore_diameters)
print(throat_diameters)
print(throat_lengths)
print(throat_volumes)
print(pore_volumes)

# <codecell>

air = OpenPNM.Fluids.Air(network = pn, name = 'air')
water = OpenPNM.Fluids.Water(network = pn, name = 'water')
pn.regenerate_fluids()

# <codecell>

phys_water = OpenPNM.Physics.GenericPhysics(network=pn,fluid=water, geometry = geo, name='standard_water_physics')
phys_air = OpenPNM.Physics.GenericPhysics(network=pn,fluid=air, geometry = geo, name='standard_air_physics')

phys_water.add_method(prop='capillary_pressure', model='purcell', r_toroid=1e-5)
phys_water.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_water.add_method(prop='diffusive_conductance', model='bulk_diffusion', shape = 'circular')
phys_air.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_air.add_method(prop='diffusive_conductance',model='bulk_diffusion', shape = 'circular')

# <codecell>

pn.regenerate_physics()

# <codecell>

inlets = pn.get_pore_indices(labels = ['bottom']) #in brackets so the whole bottom of the lattice is considered 1 pore
outlets = pn.get_pore_indices(labels = ['top'])
end_condition = 'breakthrough'

OP_1 = OpenPNM.Algorithms.InvasionPercolation(network = pn, name = 'OP_1',loglevel=30)
OP_1.run(invading_fluid = water, defending_fluid = air, inlets = inlets, outlets = outlets, end_condition = end_condition)
OP_1.update()

# <codecell>

#OpenPNM.Visualization.VTK.write(network=pn,filename = 'wu.vtp', fluids = [air, water])

# <codecell>

Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian', name = 'fickian_alg', network = pn)

z_dimension = int(pn.domain_size(dimension = 'height')/Lc) #number of pores in the z direction
quarter_layer = z_dimension/4 #estimates which layer marks 1/4 up the lattice
pore_number = quarter_layer*len(outlets) #gives the first pore in the layer 1/4 up the lattice

bottom_boundary = range(pore_number, pore_number + len(outlets))
top_boundary = range(pn.num_pores() - pore_number, pn.num_pores() - pore_number +len(outlets))

# <codecell>

pn.set_pore_info(label = 'bound1', locations = bottom_boundary)
pn.set_pore_info(label = 'bound2', locations = top_boundary)
print(pn.labels())

# <codecell>

pn.set_pore_info(label='bound1', locations = bottom_boundary)
pn.set_pore_info(label='bound2', locations = top_boundary)

BC1_pores = pn.get_pore_indices(labels= ['bound2']) #top boundary
BC2_pores = pn.get_pore_indices(labels= ['bound1']) #bottom boundary

Fickian_alg.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .6, pores = BC1_pores)
Fickian_alg.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .2, pores = BC2_pores)

# <codecell>

Fickian_alg.run(active_fluid=air)
Fickian_alg.update()

# <codecell>

#OpenPNM.Visualization.Vtp.write(filename='test.vtp',fluids=[air,water],network=pn)

# <codecell>

effective_diffusivity_original = Fickian_alg.effective_diffusivity(fluid = air) 
effective_diffusivity_original_X = Fickian_alg.effective_diffusivity(fluid = air, direction = 'x') 
effective_diffusivity_original_Y = Fickian_alg.effective_diffusivity(fluid = air, direction = 'y') 
effective_diffusivity_original_Z = Fickian_alg.effective_diffusivity(fluid = air, direction = 'x') 

# <codecell>

front_boundary = pn.get_pore_indices(labels = 'front')
back_boundary = pn.get_pore_indices(labels = 'back')
effective_diffusivity_X = Fickian_alg.effective_diffusivity(fluid = air, boundary_1 = front_boundary, boundary_2 = back_boundary)

left_boundary = pn.get_pore_indices(labels = 'left')
right_boundary = pn.get_pore_indices(labels = 'right')
effective_diffusivity_Y = Fickian_alg.effective_diffusivity(fluid = air, boundary_1 = left_boundary, boundary_2 = right_boundary)

bottom_boundary = pn.get_pore_indices(labels = 'bottom')
top_boundary = pn.get_pore_indices(labels = 'top')
effective_diffusivity_Z = Fickian_alg.effective_diffusivity(fluid = air, boundary_1 = bottom_boundary, boundary_2 = top_boundary)

# <codecell>

#print effective_diffusivity_original
#print effective_diffusivity_original_X
#print effective_diffusivity_original_Y
#print effective_diffusivity_original_Z
#
#print effective_diffusivity_X
#print effective_diffusivity_Y
#print effective_diffusivity_Z

# <codecell>

bottom_boundary = range(pore_number, pore_number + len(outlets))
top_boundary = range(pn.num_pores() - pore_number, pn.num_pores() - pore_number +len(outlets))

effective_diffusivity_Z_mod = Fickian_alg.effective_diffusivity(fluid = air, boundary_1 = bottom_boundary, boundary_2 = top_boundary)

# <codecell>

bulk_diffusivity = air.get_pore_data(prop = 'diffusivity')
normal_diffusivity = effective_diffusivity_Z/bulk_diffusivity
#print normal_diffusivity

normal_diffusivity_mod = effective_diffusivity_Z_mod/bulk_diffusivity
#print normal_diffusivity_mod

# <codecell>

np.sum(water['throat.occupancy'][pn.find_neighbor_throats(pn.get_pore_indices(labels='bottom'))])

# <codecell>

#print min(pn['pore.coords'][pn.find_neighbor_pores(pn.get_pore_indices(labels='bottom'))][2])
#print min(pn['pore.coords'][pn.get_pore_indices(labels='internal')][2])

# <codecell>

#print right_boundary

# <codecell>

left_boundary = pn.get_pore_indices(labels = 'right')
#print left_boundary

# <codecell>

pore_1 = left_boundary[0]
pore_2 = right_boundary[0]
pore_1_coords = pn.get_pore_data(prop = 'coords', locations = pore_1)
pore_2_coords = pn.get_pore_data(prop = 'coords', locations = pore_2)
#print pore_1_coords
#print pore_2_coords
for i in range(3):
            if pore_1_coords[i]!= pore_2_coords[i]:
#                print "i: ", i
                break

# <codecell>

#print left_boundary
#print [item + pn._Nx for item in left_boundary]

# <codecell>

#print effective_diffusivity_original
#print Fickian_alg.effective_diffusivity(fluid = air, direction = 'x') 
#print effective_diffusivity_X
#print effective_diffusivity_Y
#print effective_diffusivity_Z

# <codecell>

#print pn._Nx

# <codecell>

Fickian_alg.get_pore_info(label = "Dirichlet", return_indices = True)

# <codecell>

#print effective_diffusivity
#print bulk_diffusivity
#print normal_diffusivity

# <codecell>

pn.domain_size(dimension = 'top')

# <codecell>

pn._Lx

# <codecell>

(pn._Lc*8)**2

# <codecell>

#print pn._Lc

# <codecell>

#print pn._Lx/pn._Lc

# <codecell>

