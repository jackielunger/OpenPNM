# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 20:29:06 2013

@author: debarungupta
"""

import OpenPNM
import scipy as sp
import numpy as np
import scipy.stats as spst
import scipy.spatial as sptl
import itertools as itr
import math

from __GenericGeometry__ import GenericGeometry

class CubicBoundaries(GenericGeometry):
     def __init__(self, **kwargs):

        super(StitchCubicBoundaries,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")

        #Instantiate pore network object
        self._net=OpenPNM.Network.GenericNetwork()
        
    def generate(self, network, **params):
        self._logger.debug("Execute RunBoundaryStitch")
        