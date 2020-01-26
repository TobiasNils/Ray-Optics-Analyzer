#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# cython: profile=True

#------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amezquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Tobias Nils Ackermann
# Description:     Rectangle definition module
# Symbols Defined: Polygon
#------------------------------------------------------------------------------
#

"""Module that defines the Polygon class
"""


from pyoptools.raytrace.shape.shape cimport Shape
from matplotlib.path import Path

from numpy import arange, meshgrid, where, dot, array, arccos, sin, cos, sqrt
from numpy.linalg import norm
# cimport numpy as np


cdef class Polygon(Shape):
    ''' Class defining a polygon shape.

    coord -> tuple containing the 3 (x,y) coordinates of the corners of
             the triangle
    samples -> number of divitions per side used to sample the triangle
    '''
    # cdef public tuple coord
    # cdef public int samples
    # cdef public tuple uv_poly
    # cdef public tuple uv

    def __init__(self,coord=(array([0.,0.,0.]),array([0.,100.,0.]),array([100.,0.,0.])),samples=10,*args, **kwargs):
        Shape.__init__(self,*args, **kwargs)
        self.coord=coord
        self.samples=samples
        # try:
        #     self.uv_poly = ()
        # except Exception as e:
        #     print('Polygon.__init__: self.uv_poly = ()', e)
        # try:
        #     self.uv = ()
        # except Exception as e:
        #     print('Polygon.__init__: self.uv = ()', e)

        #Register picklable attributes
        self.addkey("coord")
        self.addkey("samples")
        # self.addkey("uv_poly")
        # self.addkey("uv")

    def __reduce__(self):

        state=None
        args=(self.coord, self.samples)
        return(type(self),args)


    cpdef hit(self, p, u,v, uv_poly):
        """Method  that returns True if a p=(x,y,z) point is inside the triangle,
        if not it returns False.
        taken from http://www.blackpawn.com/texts/pointinpoly/default.html
        """
        local_origin = self.coord[0]
        p = p-local_origin
        norm_p = norm(p)
        p_=p/norm_p
        theta = arccos(dot(p_, u))
        u_fac = norm_p*cos(theta)
        v_fac = norm_p*sin(theta)

        retval = [u_fac,v_fac]

        path = Path(list(uv_poly))
        return path.contains_point(retval, radius=10e-10)


    cpdef bint fhit(self,double px,double py,double pz):
        """This method returns TRUE if an p=(x,y,z)point is inside the surface
        aperture if not it must return FALSE.
        This is implemented for a point, in cython, to make it fast
        """
        # cdef double dot00,dot01,dot02,dot11,dot12,invDenom,u,v
        # This one needs to be optimized
        # print(self.coord[0])
        P=array((px,py,pz))
        # print(P)
        P1_uv = dot(self.uv[0],P)
        # print(P1_uv)
        P2_uv = dot(self.uv[1],P)
        # print(P2_uv)
        P_uv = [P1_uv, P2_uv]
        path = Path(self.uv_poly)
        return path.contains_point(P_uv, radius=10e-6)


    #~ cpdef polylist(self, topo): #Falta organizar el polilist
        #~ """Method that returns a tuple (point_list, poly_list) for a triangular mesh.
        #~
        #~ Attributes:
        #~ ===========
        #~
        #~ topo    Z=topo(x,y) is the function that gives the surface topography
        #~
        #~ The point list is a list of tuples (X,Y,Z) containing the coordinates of
        #~ the points used to build the surface mesh.
        #~ The poly_list is a list of tuples (n1,n2,n3,n3) containing the indices
        #~ of the points in the polylist used to build each polygon that will be
        #~ used to visualize the mesh.
        #~ """
        #~
        #~ cdef int i,j
        #~
        #~ A=array(self.coord[0])
        #~ B=array(self.coord[1])
        #~ C=array(self.coord[2])
        #~
        #~ #Get the mesh points
        #~ points=[]
        #~ for i in range(self.samples+1):
            #~ P0= A+i*(B-A)/self.samples
            #~ P1= A+i*(C-A)/self.samples
            #~ for j in range(i+1):
                #~ if i!=0:
                    #~ P=P0+(P1-P0)*j/i
                #~ else:
                    #~ P=P0
                #~ Z=topo(P[0],P[1])
                #~ points.append((P[0],P[1],Z))
                #~
        #~ from matplotlib.delaunay import delaunay
        #~
        #~ #Need to find a beter way to do this not using delaunay# or maybe to generate all using triangulations????
        #~
        #~ x=[p[0] for p in points]
        #~ y=[p[1] for p in points]
        #~ cs,e,trip,trin=delaunay(x,y)
        #~ return points, trip
    #
    # cpdef pointlist(self):
    #
    #     cdef int i,j
    #
    #     A=array(self.coord[0])
    #     B=array(self.coord[1])
    #     C=array(self.coord[2])
    #
    #     #Get the mesh points
    #     X=[]
    #     Y=[]
    #     for i in range(self.samples+1):
    #         P0= A+i*(B-A)/self.samples
    #         P1= A+i*(C-A)/self.samples
    #         for j in range(i+1):
    #             if i!=0:
    #                 P=P0+(P1-P0)*j/i
    #             else:
    #                 P=P0
    #             X.append(P[0])
    #             Y.append(P[1])
    #     return X,Y
    #
    # cpdef limits(self):
    #     """
    #     Returns the minimum limits for the aperture
    #     """
    #     dx, dy=self.size
    #     return -dx/2, dx/2, -dy/2, dy/2
