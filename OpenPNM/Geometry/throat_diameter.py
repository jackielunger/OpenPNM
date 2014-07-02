r"""
===============================================================================
Submodule -- throat_diameter
===============================================================================

"""
import scipy as sp
import scipy.stats as spst
import random as rd

def constant(geometry,
             network,
             propname,
             value,
             **params):
    r"""
    Assigns specified constant value
    """
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

def uniform_distribution(geometry,
                         network,
                         propname,
                         low = 4,
                         high = 6,
                         seed = None,
                         **params):
    
    r"""
        Assigns a uniform distribution of values
        """
    sp.random.seed(seed)
    Np = network.num_throats(geometry.name)
    value = low + sp.random.rand(Np)*(high-low) 
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

def smaller_connecting_pore(geometry,
                            network,
                            propname,
                            **params):

    r"""
    sets diameter to that of the smaller neighbor pore
    """
    
    connected_pores = network.find_connected_pores(geometry.throats())
    value = []
    for pair in connected_pores:
        value.append(min(network.get_data(pores =pair[0], prop = propname), network.get_data(pores =pair[1], prop = propname)))
    network.set_data(prop = propname, throats = geometry.throats(), data = value)   

def cylinder(geometry,
             network,
             propname,
             seed='seed',
             **params):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    """
    prob_fn = getattr(spst,params['name'])
    P = prob_fn(params['shape'],loc=params['loc'],scale=params['scale'])
    value=P.ppf(network.get_data(prop=seed,throats=geometry.throats()))
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

def cuboid(geometry,
           network,
           propname,
           **params):
    r"""
    Calculate throat diameter from seeds for a cuboidal throat
    """
    print('cuboid: nothing yet')
    
def voronoi(geometry,
            network,
            propname,
            **params):
    r"""
    Calculate throat diameter from analysis of Voronoi facets
    """
    print('voronoi: nothing yet')