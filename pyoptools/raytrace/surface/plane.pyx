#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# cython: profile=True

#------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amézquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amézquita Orozco
# Description:     Plane surface definitión module
# Symbols Defined: Plane
#------------------------------------------------------------------------------

'''Module that defines a plane surface class
'''


from numpy import array, dot, cross, inf, float64, zeros, asarray
#from enthought.traits.api import Tuple,Float
#from enthought.traits.ui.view import Group,Item

#from ray_trace.surface import Surface
from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray

from pyoptools.misc.definitions  import *
from pyoptools.misc.cmisc.cmisc cimport norm_vect
cimport numpy as np

import cython

cdef class Plane(Surface):
    '''Class to define a plane surface.

    Description:

    Plane is a class to define rectangular plane surfaces.
    The surface shape is given by the shape attribute

    Example:

        >>> ps=Plane(shape=Rectangular(size=(25,15)))
    '''

    def __init__(self,*args, **kwargs):
        Surface.__init__(self,*args, **kwargs)

    cpdef topo(self, x, y):
        return None

    #~ def __reduce__(self):
        #~
        #~ args=(self.reflectivity, self.shape)
        #~ return(type(self),args,self.__getstate__())
 #~


    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef _intersection(self,Ray A):
        """Returns the intersection point between a ray and an the plane in global 3D coordinates

        """
        # cdef np.ndarray w,retval
        cdef double epsilon=1e-6
        # get surface data
        cdef np.ndarray[np.float64_t, ndim=1] planePoint=self.shape.coord[0]
        planeNormal=self.normal(A)

        # Punto que pertenece al rayo "Origen" del rayo
        cdef np.ndarray[np.float64_t, ndim=1] rayPoint=A.pos
        #Vector paralelo a la linea
        cdef np.ndarray[np.float64_t, ndim=1] rayDirection=A.dir

        #if dot(N_,A.dir) ==0 : return inf_vect
        cdef double ndotu=dot(planeNormal,rayDirection)
        if abs(ndotu) < epsilon: return inf_vect # ray parallel or within the plane
        cdef np.ndarray w=rayPoint-planePoint
        cdef double si=dot(planeNormal,w)
        si = -si/ndotu
        cdef np.ndarray retval=w+si*rayDirection+planePoint
        return retval

    cpdef np.ndarray normal(self, ri):
        """Method that returns the normal to the surface in global coordinates
        """
        cdef np.ndarray[np.float64_t, ndim=1] v1=self.shape.coord[0]
        cdef np.ndarray[np.float64_t, ndim=1] v2=self.shape.coord[1]
        cdef np.ndarray[np.float64_t, ndim=1] v3=self.shape.coord[2]

        # Take two in-plane vectors
        p1 = v2-v1
        p2 = v3-v1
        # Get the vector perpendicular to the plane p1 and p2 span by cross-product
        N_ = cross(p1,p2)
        # normalize
        N_ = norm_vect(N_)

        N_=array(N_).astype(float64)
        return (N_)

    def _repr_(self):
        '''Return an string with the representation of the optical plane
        '''

        return "Plane(shape="+str(self.shape)+",reflectivity="+ \
                      str(self.reflectivity)+")"
