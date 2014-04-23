
"""
module molar_mass
===============================================================================

"""
import scipy as sp

def constant(fluid,network,propname,value,**params):
    r"""
    Assigns specified constant value
    """
    network.set_pore_data(phase=fluid,prop=propname,data=value)

def na(fluid,network,propname,**params):
    r"""
    Assigns nonsensical, but numerical value of -1.  
    This ensurse stability of other methods 
    but introduces the possibility of being misused.
    """
    value = -1
    network.set_pore_data(phase=fluid,prop=propname,data=value)

def mixture(fluid,network,propname='MW',**params):
    r"""
    Calculates the average molecular weight of a mixture using mole fraction weighting
    """
    MW = []
    mole_frac = []
    for item in fluid.components.keys():
        MW.append(fluid.components[item].get_pore_data(prop=propname))
        mole_frac.append(fluid.composition[item])
    MW = sp.array(MW,ndmin=1)
    mole_frac = sp.array(mole_frac,ndmin=1)
    print(MW,mole_frac)
    #Ensure mole fraction sum to 1
    fsum = sp.sum(mole_frac)
    if fsum != 1: print(fluid._logger.warning('mole fractions do not add to 1, so performing normalization'))
    value = sp.sum(MW*mole_frac)/fsum
    fluid.set_pore_data(prop=propname,data=value)