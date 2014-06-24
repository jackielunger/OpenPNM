r"""
===============================================================================
Submodule -- diffusive_conductance
===============================================================================

"""

import scipy as sp

def constant(physics,
             network,
             geometry,
             fluid,
             propname,
             value,
             **params):
    r"""
    Assigns specified constant value
    """
    fluid.set_data(prop=propname,throats=geometry.throats(),data=value)

def na(physics,
       network,
       geometry,
       fluid,
       propname,
       **params):
    r"""
    """
    value = -1
    fluid.set_data(prop=propname,throats=geometry.throats(),data=value)

def bulk_diffusion(physics,
                   network,
                   fluid,
                   geometry,
                   propname,
                   shape = 'square',
                   diffusivity = 'diffusivity',
                   molar_density = 'molar_density',
                   throat_diameter = 'diameter',
                   throat_length = 'length',
                   pore_diameter = 'diameter',
                   **params):
    r"""
    Calculate the diffusive conductance of conduits in network, where a 
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas

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
    #Interpolate pore fluid property values to throats
    ct = fluid.get_data(prop='molar_density',throats='all',mode='interpolate')
    DABt = fluid.get_data(prop='diffusivity',throats='all',mode='interpolate')
    #Get Nt-by-2 list of pores connected to each throat
    Ts= network.get_throat_indices()
    Ps = network.find_connected_pores(Ts,flatten=0)
    #Find g for half of pore 1
    pdia = network.get_data(prop=pore_diameter,pores='all')
    gp1 = ct*DABt*pdia[Ps[:,0]]**2/(0.5*pdia[Ps[:,0]])
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = ct*DABt*pdia[Ps[:,1]]**2/(0.5*pdia[Ps[:,1]])
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tdia = network['throat.diameter']
    tlen = network['throat.length']
    if (shape == 'circular'):
        gt = sp.pi*ct*DABt*tdia**2/(tlen*4)
    elif (shape == 'square'):
        gt = ct*DABt*tdia**2/tlen
    else:
        print('invalid shape chosen.  Either circular or square')
        return
    g = (1/gt + 1/gp1 + 1/gp2)**(-1) 
    g = g[geometry.throats()]
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
        g = g * s + g*(not s)/1.0e3
    except: pass    
    fluid['throat.'+propname] = g

