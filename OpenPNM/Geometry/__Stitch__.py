# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 20:28:29 2013

@author: debarungupta
"""

import OpenPNM
import scipy as sp

from __Cubic__ import Cubic
from __GenericGeometry import GenericGeometry

class Stitch(GenericGeometry, Cubic): # Multiple inheritance. The functions will be searched for in Generic before Cubic. 

    def __init__(self, **kwargs):

        super(Stitch,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")

        #Instantiate pore network object
        self._net=OpenPNM.Network.GenericNetwork()
    
    @staticmethod
    def StitchNetworks(self, net1, net2):
        r"""
        Stitch two networks together. User will translate or scale beforehand. 

        Parameters
        ----------
        net1 : OpenPNM Network Object
            The network that is stiched to

        net2 : OpenPNM Network Object
            The network that is stitched

        - This should instead create a new network altogether. The new network will take pore and throat attributes from net1 and net2. 
        (Want to enstantiate new network)

        - The number of pores in total will stay the same, but the number of throats will increase. We must add these new throats by doing
        a delaunay tesselation and removing appropriate connections.

        - We have to fully define the new throats that we have generated. The connections, type, and numbering are given in the stitch throats script, 
        but we must add in seeds, diameters, lengths and volumes. Each are dependent on the last. However, the seeds are also dependent on pore seeds.
        We must add a new list of seeds to the throats somehow ?? 

        """
        self._logger.debug("Stitch Networks : Not Implemented Yet")
        
    @staticmethod
    def StitchBoundaries(self, net, **params):
        r"""
        - Enstantiate 6 NEW cubic networks. (How??)
        - Generate 6 new networks with these generated properties, as called by generate in the Generic class .
            - Generate _generate_setup  (Cubic)             [ Gives Nx, Ny, Nz, Lx, Ly, Lz, Lc]
            - Generate NEW pores through translation.       [ Gives Coords, Numbering, Type in the pore_properties dictionary ]
                - Using translate_coordinates, we will move the pores by Lc in the appropriate direction.
            - Generate pore seeds       (GenericGeometry)   [ Gives Seeds in pore_properties dictionary]
            - Generate pore_diameters   (GenericGeometry)   [ Gives diameter in pore_properties dictionary]
            - Generate pore_volumes     (GenericGeometry)   [ Gives volume in pore_properties dictionary]

            :: (Is using multiple inheritance good practice??)

        - At this point, we will have a full list of pore_properties and a empty throat_properties dictionary. We need to add throats in a similar way, but no boundary.
            - generate NEW throats      (Cubic)             [ Gives connections, type, and numbering in throat_properties]
            - generate throat Seeds     (GenericGeometry)   [ Gives Seeds in throat_properties dictionary]
            - generate throat diameter  (GenericGeometry)   [ Gives diameter in throat_properties dictionary]
            - generate throat length    (GenericGeometry)   [ Gives length in throat_properties dictionary]
            - generate throat volume    (GenericGeometry)   [ Gives volume in throat_properties dictionary]

        - We now have a full network (one out of the 6) fully defined. We should now called StitchNetworks to connect pore and throat properties.

        """
        #Instantiate pore network object. PUT THIS IN A LOOP IN THE FUTURE.
        self._net_bound_1 = OpenPNM.Network.GenericNetwork()
        self._net_bound_2 = OpenPNM.Network.GenericNetwork()
        self._net_bound_3 = OpenPNM.Network.GenericNetwork()
        self._net_bound_4 = OpenPNM.Network.GenericNetwork()
        self._net_bound_5 = OpenPNM.Network.GenericNetwork()
        self._net_bound_6 = OpenPNM.Network.GenericNetwork()




    def _stitch_pores(self, net1, net2):


        self._net.pore_properties['coords']      = sp.concatenate((net1.pore_properties['coords'],net2.pore_properties['coords']),axis = 0)
        net2.pore_properties['numbering']        = len(net1.pore_properties['numbering']) + net2.pore_properties['numbering']
        self._net.pore_properties['numbering']   = sp.concatenate((net1.pore_properties['numbering'],net2.pore_properties['numbering']),axis=0)
        self._net.pore_properties['type']        = sp.concatenate((net1.pore_properties['type'],net2.pore_properties['type']),axis = 0)
        self._net.pore_properties['seed']        = sp.concatenate((net1.pore_properties['seed'],net2.pore_properties['seed']),axis = 0)
        self._net.pore_properties['diameter']    = sp.concatenate((net1.pore_properties['diameter'],net2.pore_properties['diameter']),axis = 0)
        self._net.pore_properties['volume']      = sp.concatenate((net1.pore_properties['volume'],net2.pore_properties['volume']),axis = 0)


        
    def _stitch_throats(self, net1, net2):
        
        
        net2.throat_properties['numbering']      = len(net1.throat_properties['numbering']) + net2.throat_properties['numbering']
        self._net.throat_properties['numbering'] = sp.concatenate((net1.throat_properties['numbering'],net2.throat_properties['numbering']),axis=0)
        self._net.throat_properties['seed']      = sp.concatenate((net1.throat_properties['seed'],net2.throat_properties['seed']),axis=0)
        self._net.throat_properties['diameter']  = sp.concatenate((net1.throat_properties['diameter'],net2.throat_properties['diameter']),axis=0)
        self._net.throat_properties['volume']    = sp.concatenate((net1.throat_properties['volume'],net2.throat_properties['volume']),axis=0)
        self._net.throat_properties['length']    = sp.concatenate((net1.throat_properties['length'],net2.throat_properties['length']),axis=0)
        #self._net.throat_properties['type']
        #self._net.throat_properties['connections']
        
        #Type and connections not concatenated, they must be generated using a separate script. This is where the delaunay script comes in. 
        #Check that the types and connections are correct with VTK. 
        
        
    @staticmethod
    def translate_coordinates(self,net,displacement=[0,0,0]):
        r"""
        Translate pore network coordinates by specified amount

        Parameters
        ----------
        net : OpenPNM Network Object
            The network to which translation should be applied

        displacement : array_like
            A vector containing the amount to translate in each dimension. [0,0,0] yeilds no translation.

        """
        net.pore_properties['coords'] = net.pore_properties['coords'] + displacement

    @staticmethod
    def scale_coordinates(self,net,scale=[1,1,1]):
        r"""
        Scale pore network coordinates by specified amount

        Parameters
        ----------
        net : OpenPNM Network Object
            The network to which translation should be applied

        scale : array_like
            A vector containing the amount to scale in each dimension.  [1,1,1] yeilds no scaling.

        """
        net.pore_properties['coords'] = net.pore_properties['coords']*scale