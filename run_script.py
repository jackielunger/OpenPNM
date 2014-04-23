import OpenPNM
import scipy as sp

#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(name='cubic_1',loglevel=20).generate(divisions=[15, 15, 15], lattice_spacing=[0.0001],add_boundaries=True)

#==============================================================================
'''Build Geometry'''
#==============================================================================
geom = OpenPNM.Geometry.Stick_and_Ball(network=pn, name='geom')
geom.regenerate()

#==============================================================================
'''Build Fluids'''
#==============================================================================
air = pn.add_fluid('air')
N2 = pn.add_fluid('N2')
N2.apply_conditions(MW = 0.028)
O2 = pn.add_fluid('O2')
O2.apply_conditions(MW = 0.032)
air.components = {}
air.components['N2'] = N2
air.components['O2'] = O2
air.composition = {}
air.composition['O2'] = 0.21
air.composition['N2'] = 0.79

air.add_method(prop='molar_mass',model='mixture')
air.molar_mass()



#==============================================================================
'''Build Physics Objects'''
#==============================================================================
phys_water = OpenPNM.Physics.GenericPhysics(network=pn, fluid=water,geometry=geom,name='phys_water')
phys_water.add_method(prop='capillary_pressure', model='washburn')
phys_water.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_water.add_method(prop='diffusive_conductance', prop_name='gdAB', model='bulk_diffusion', diffusivity='DAB')
phys_water.regenerate()

phys_air = OpenPNM.Physics.GenericPhysics(network=pn, fluid=air,geometry=geom, name='phys_air')
phys_air.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_air.add_method(prop='diffusive_conductance', model='bulk_diffusion')
phys_air.regenerate()

#==============================================================================
'''Begin Simulations'''
#==============================================================================
'''Perform a Drainage Experiment (OrdinaryPercolation)'''
#------------------------------------------------------------------------------
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(loglevel=20,loggername='OP',name='OP_1',network=pn)
a = pn.get_pore_indices(labels=['bottom','boundary'],mode='intersection')
OP_1.setup(invading_fluid='water',defending_fluid='air',inlets=a,npts=20)
OP_1.run()
#OP_1.plot_drainage_curve()

#------------------------------------------------------------------------------
'''Perform Fickian Diffusion'''
#------------------------------------------------------------------------------
Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loglevel=20, loggername='Fickian', name='Fickian_alg',network=pn)
# Assign Dirichlet boundary conditions to top and bottom surface pores
BC1_pores = pn.get_pore_indices(labels=['top','boundary'],mode='intersection')
Fickian_alg.set_pore_info(label='Dirichlet', locations=BC1_pores)
Fickian_alg.set_pore_data(prop='BCval', data=0.6, locations=BC1_pores)
BC2_pores = pn.get_pore_indices(labels=['bottom','boundary'],mode='intersection')
Fickian_alg.set_pore_info(label='Dirichlet', locations=BC2_pores)
Fickian_alg.set_pore_data(prop='BCval', data=0.2, locations=BC2_pores)

# Updating data based on the result of Percolation Algorithms
OP_1.update(Pc=3000)
# Run simulation
Fickian_alg.run(active_fluid=air)
Fickian_alg.update()

#------------------------------------------------------------------------------
#Export to VTK
OpenPNM.Visualization.VTK().write(net=pn, fluids=[air,water])
