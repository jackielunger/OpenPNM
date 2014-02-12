
"""
module SurfaceTension
===============================================================================

"""

import scipy as sp

def constant(fluid,network,propname,value,**params):
    r"""
    Assigns specified constant value
    """
    network.set_pore_data(phase=fluid,prop=propname,data=value)

def na(fluid,network,propname,**params):
    value = -1
    network.set_pore_data(phase=fluid,prop=propname,data=value)

def empirical(fluid,network,propname,a=[0],**params):
    r"""
    Uses a polynomial fit of empirical data to calculate property
    """
    T = network.get_pore_data(phase=fluid,prop='temperature')
    value = sp.zeros_like(T)
    for i in range(0,sp.size(a)):
        value = value + a[i]*(T**i)
    network.set_pore_data(phase=fluid,prop=propname,data=value)

def Eotvos(fluid,network,propname,k=2.1e-7,**params):
    r"""
    """
    Tc = fluid.get_pore_data(prop='Tc')
    T = network.get_pore_data(phase=fluid,prop='temperature')
    Vm = 1/network.get_pore_data(phase=fluid,prop='molar_density')
    value = k*(Tc-T)/(Vm**(2/3))
    network.set_pore_data(phase=fluid,prop=propname,data=value)

def GuggenheimKatayama(fluid,network,propname,K2=1,n=1.222,**params):
    r"""
    """
    T = network.get_pore_data(phase=fluid,prop='temperature')
    Pc = fluid.get_pore_data(prop='Pc')
    Tc = fluid.get_pore_data(prop='Tc')
    sigma_o = K2*Tc**(1/3)*Pc**(2/3)
    value = sigma_o*(1-T/Tc)**n
    network.set_pore_data(phase=fluid,prop=propname,data=value)

def BrockBird_scaling(fluid,network,propname,sigma_o=0.072,To=298.0,**params):
    r"""
    Uses Brock_Bird model to adjust surface tension from it's value at a given reference temperature to temperature of interest

    Parameters
    ----------
    fluid : OpenPNM Fluid Object

    sigma_o : float
        Surface tension at reference temperature (N/m)

    To : float
        Temperature at reference conditions (K)
    """
    Tc = fluid.get_pore_data(prop='Tc')
    Ti = network.get_pore_data(phase=fluid,prop='temperature')
    Tro = To/Tc
    Tri = Ti/Tc
    value = sigma_o*(1-Tri)**(11/9)/(1-Tro)**(11/9)
    network.set_pore_data(phase=fluid,prop=propname,data=value)



