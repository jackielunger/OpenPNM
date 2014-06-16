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
                         **params):
    
    r"""
        Assigns a uniform distribution of values
        """
    throats = network.get_throat_indices()
    for throat in throats: #loops through all throats and sets the throat a random radius (between high and low)
        lamd = rd.random()
        radius = low + lamd*(high - low)
        diameter = 2*radius
        network.set_throat_data(prop = propname, data = diameter, locations = throat)


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