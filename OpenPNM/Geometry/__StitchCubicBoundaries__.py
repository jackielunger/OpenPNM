# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 20:29:06 2013

@author: debarungupta
"""

import OpenPNM
from __GenericGeometry__ import GenericGeometry

class CubicBoundaries(GenericGeometry):
     def __init__(self, **kwargs):
        super(CubicBoundaries,self).__init__(**kwargs)
        self._net=OpenPNM.Network.GenericNetwork()
     
     def _generate(self):
        self._logger.debug("Execute RunBoundaryStitch")    
        

if __name__ == '__main__':
    test=CubicBoundaries(loggername='TestCubicBoundaries')