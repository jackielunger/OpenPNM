# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 20:28:29 2013

@author: debarungupta
"""

import OpenPNM
import scipy as sp
import numpy as np
import itertools as itr
import scipy.spatial as sptl

from __Cubic__ import Cubic
#from __GenericGeometry__ import GenericGeometry

class Stitch( Cubic): # Multiple inheritance. The functions will be searched for in Generic before Cubic. 
    r"""

    The point of creating an extra class here is so that we don't have to overwrite self properties within the Cubic Class. 
    We will be able to fill in attributes of a network for each boundary to an instance of each boundary network.

    """

    def __init__(self, **kwargs):

        super(Stitch,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")
    
    
    def StitchNetworks(self, net1, net2, **params):
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
        #############################################################################
        # JUST FOR TEST PURPOSES. TRANSLATE COORDINATES MUST MOVE net2 BEFORE STITCH IS CALLED IN THE SCRIPT.
        self._generate_setup(**params)
        z_trans = net2.pore_properties['coords'][:,2].max() + self._Lc/2
        net2 = super(Cubic,self).translate_coordinates(net2, displacement = [0,0,z_trans])
        #############################################################################       


        self._net = OpenPNM.Network.GenericNetwork()        
        self._stitch_pores(net1, net2)      
        self._stitch_throats(**params)  # These parameters are of the first network. They will be used to generate the new pores
                                        # All other properties that exist from net 2 will replace the newly stitched set of throats and properties. 
        large = len(self._net.throat_properties['type'])
        small = len(net2.throat_properties['type'])
        
        self._net.throat_properties['type']         = sp.zeros(large)
        self._net.throat_properties['volume']       = sp.concatenate((self._net.throat_properties['volume'][0:large-small], net2.throat_properties['volume']))
        self._net.throat_properties['diameter']     = sp.concatenate((self._net.throat_properties['diameter'][0:large-small], net2.throat_properties['diameter']))
        self._net.throat_properties['numbering']    = sp.concatenate((self._net.throat_properties['numbering'][0:large-small], net2.throat_properties['numbering']))
        self._net.throat_properties['length']       = sp.concatenate((self._net.throat_properties['length'][0:large-small], net2.throat_properties['length']))
        self._net.throat_properties['seed']         = sp.concatenate((self._net.throat_properties['seed'][0:large-small], net2.throat_properties['seed']))
                                                
        self._logger.debug("Stitch Networks : Not Implemented Yet")
        return self._net


    def StitchBoundaries(self, net_original, **params):
        
        r"""

        - Instantiate 6 NEW cubic networks. (How??)
        - Generate 6 new networks with these generated properties, as called by generate in the Generic class .
            - Generate _generate_setup  (Cubic)             [ Gives Nx, Ny, Nz, Lx, Ly, Lz, Lc]
            - Generate NEW pores through translation.       [ Gives Coords, Numbering, Type in the pore_properties dictionary]
                - Using translate_coordinates, we will move 
                the pores by Lc in the appropriate direction.
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

        btype = params['btype']
        btype_ind = np.where(np.array(btype)==0)
        
        self._net = OpenPNM.Network.GenericNetwork()
        self._net_original = net_original
        coords = self._net_original.pore_properties['coords']
        
        for i in range(0,len(coords)):
            neighbor_coords = self._net_original.get_neighbor_pores(i)
            
            if len(neighbor_coords) != 6:
                for j in range(0,len(btype_ind[0])): # Maximum of 3 iterations.
                
                    periodic = btype_ind[0][j]
                    neighbor_ind = coords[neighbor_coords,periodic]
                    current_ind = coords[i][periodic]
                    border = np.where(neighbor_ind != current_ind)[0]
                    
                    if len(border) != 2:
                        new_coords_temp = coords[i].copy()
                        new_coords_temp[periodic] = coords[i][periodic]-(neighbor_ind[border]-current_ind)
                        coords = np.vstack((coords,new_coords_temp))
        
        
        self._net.pore_properties['coords'] = coords
        self._net.pore_properties['numbering'] = np.arange(len(coords))
        self._net.pore_properties['type']= np.zeros(len(coords),dtype=np.int8)
        self._generate_pore_seeds()                         
        self._generate_pore_diameters(params['psd_info'])   
        self._calc_pore_volumes() 
        
        self._stitch_throats( **params )
        self._set_pore_types(btype)        # The next 3 calls are specific to boundary generation, so thats why they're outside of stitch _throats. 
        self._set_throat_types()            # They set the correct types of pore and throat boundaries, and remove the external to external connections using refine. 
        self._refine_boundary_throats()
        return self._net

            
    def _set_pore_types(self, b_type):
        pore_type = np.zeros(len(self._net.pore_properties['type']))
        coords    = self._net.pore_properties['coords']
        t_count   = 0

        for i in range(2,-1,-1):
            if not b_type[i]:
                bound_1 = coords[:,i].min()
                bound_2 = coords[:,i].max()
                bound_ind_1 = np.where(coords[:,i] == bound_1)
                bound_ind_2 = np.where(coords[:,i] == bound_2)
                pore_type[bound_ind_1] = 1 + t_count
                pore_type[bound_ind_2] = 6 - t_count
            t_count = t_count + 1
        
        self._net.pore_properties['type'] = pore_type

    def _set_throat_types(self):
        r"""
        """
        
        for i in range(0,len(self._net.throat_properties['type'])):
            temp1 = self._net.pore_properties['type'][self._net.throat_properties['connections'][i,0]]
            temp2 = self._net.pore_properties['type'][self._net.throat_properties['connections'][i,1]]
            if min(temp1,temp2) > 0:
                self._net.throat_properties['type'][i] = min(temp1,temp2)

    def _refine_boundary_throats(self):
        r"""
        """
        
        mask = np.where(self._net.throat_properties['type'] == 0)
        self._net.throat_properties['volume']       = self._net.throat_properties['volume'][mask]
        self._net.throat_properties['diameter']     = self._net.throat_properties['diameter'][mask]
        self._net.throat_properties['numbering']    = self._net.throat_properties['numbering'][mask]
        self._net.throat_properties['connections']  = self._net.throat_properties['connections'][mask]
        self._net.throat_properties['length']       = self._net.throat_properties['length'][mask]
        self._net.throat_properties['seed']         = self._net.throat_properties['seed'][mask]
        self._net.throat_properties['type']         = self._net.throat_properties['type'][mask]
        
    
    def _generate_new_throats(self):
        r"""
        Stitch two networks together OR adds the boundary throats to an existing network
    
        Parameters
        ----------
        net : OpenPNM Network Object
            The network that is stiched, whos throats are being added.
    
        """
        
        pts = self._net.pore_properties['coords']
        tri = sptl.Delaunay(pts)
        I = []
        J = []
        V = []
    
        dist_comb = list(itr.combinations(range(4),2)) # [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]. This compares each colony with its neighbor with no repeats.
    
        for i in range(len(tri.simplices)):
            for j in range(len(dist_comb)):
                point_1 = tri.simplices[i,dist_comb[j][0]]
                point_2 = tri.simplices[i,dist_comb[j][1]]
                coords_1 = tri.points[point_1]
                coords_2 = tri.points[point_2]
                dist =  np.sqrt((coords_2[0] - coords_1[0]) ** 2 + (coords_2[1] - coords_1[1]) ** 2 + (coords_2[2] - coords_1[2]) ** 2)
                V.append(dist)
                I.append(point_1)
                J.append(point_2)
    
                    # Takes all the IJV values and puts it into a sparse matrix.
        spar_connections = sp.sparse.coo_matrix((np.array(V),(np.array(I),np.array(J))))
        ind = np.where(spar_connections.data < min(spar_connections.data) + 0.001)
        prelim_connections = np.array((spar_connections.row[ind],spar_connections.col[ind]))
        connections = np.zeros((len(prelim_connections[0]),2),np.int)
    
        for i in range(len(prelim_connections[0])):
            connections[i,0] = prelim_connections[0][i]
            connections[i,1] = prelim_connections[1][i]
    
        b = np.ascontiguousarray(connections).view(np.dtype((np.void, connections.dtype.itemsize * connections.shape[1])))
        _, idx = np.unique(b, return_index=True) # TAKEN FROM STACK OVERFLOW FOR THE TIME BEING.
        connections = connections[idx]
    
        self._net.throat_properties['connections'] =  connections
        self._net.throat_properties['numbering'] = np.arange(0,len(connections[:,0]))
        self._net.throat_properties['type'] = np.zeros(len(connections[:,0]),dtype=np.int8)
        
        
    def _stitch_pores(self, net1, net2):

        self._net.pore_properties['coords']      = sp.concatenate((net1.pore_properties['coords'],net2.pore_properties['coords']),axis = 0)
        net2.pore_properties['numbering']        = len(net1.pore_properties['numbering']) + net2.pore_properties['numbering']
        self._net.pore_properties['numbering']   = sp.concatenate((net1.pore_properties['numbering'],net2.pore_properties['numbering']),axis=0)
        self._net.pore_properties['type']        = sp.concatenate((net1.pore_properties['type'],net2.pore_properties['type']),axis = 0)
        self._net.pore_properties['seed']        = sp.concatenate((net1.pore_properties['seed'],net2.pore_properties['seed']),axis = 0)
        self._net.pore_properties['diameter']    = sp.concatenate((net1.pore_properties['diameter'],net2.pore_properties['diameter']),axis = 0)
        self._net.pore_properties['volume']      = sp.concatenate((net1.pore_properties['volume'],net2.pore_properties['volume']),axis = 0)
    
    def _stitch_throats(self, **params):
        
        self._generate_new_throats()    # This should give us throat connections, numbering, and type for the full matrix with Delaunay
        self._generate_throat_seeds()
        self._generate_throat_diameters(params['tsd_info'])
        self._calc_throat_lengths()
        self._calc_throat_volumes()
    
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
        return net
        
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

if __name__ == '__main__':
    test=Stitch(loggername='TestStitch')