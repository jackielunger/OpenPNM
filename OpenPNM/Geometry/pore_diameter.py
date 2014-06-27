r"""
===============================================================================
Submodule -- pore_diameter
===============================================================================

"""
import scipy as sp
import scipy.stats as spst


def constant(geometry,
             network,
             propname,
             value,
             **params):
    r"""
    Assign specified constant value
    """
    network.set_data(prop=propname,pores=geometry.pores(),data=value)
    
def weibull_cumulative(geometry,
                       network,
                       propname,
                       xmax,
                       bmin,
                       lmbda,
                       k,
                       seed = None,
                       **params):
    r"""
    Calculate pore diameter using a weibull cumulative distribution
    """
    sp.random.seed(seed)
    Np = network.num_pores(geometry.name)
    value = lmbda*(-sp.log(1-sp.random.rand(Np)*xmax))**(-1/k) + bmin 
    network.set_data(prop=propname,pores=geometry.pores(),data=value)

def sphere(geometry,
           network,
           propname,
           seed='seed',
           **params):
    r"""
    Calculate pore diameter from seed values for a spherical pore body
    """
    prob_fn = getattr(spst,params['name'])
    P = prob_fn(params['shape'],loc=params['loc'],scale=params['scale'])
    value = P.ppf(network.get_data(prop=seed,pores=geometry.pores()))
    network.set_data(prop=propname,pores=geometry.pores(),data=value)

    