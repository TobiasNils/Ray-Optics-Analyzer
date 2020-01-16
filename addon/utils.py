# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>
import numpy as np
import bpy
from pyoptools.raytrace.component import Component
from pyoptools.raytrace.system import System
import pyoptools.raytrace.surface as surfaces
import pyoptools.raytrace.shape as shapes
import pyoptools.raytrace.ray as rays
import pyoptools.misc.cmisc as cmisc


def get_raypaths(system):
    propagated = system._p_rays
    paths = [[tuple(r.pos)]+ray2list(r) for r in propagated]
    return paths


def ray2list(ray):
    rays=[]

    P1 = ray.pos
    if len(ray.childs) > 0:
        P2 = ray.childs[0].pos
    else:
        P2 = P1 + 10. * ray.dir

    if ray.intensity != 0:

        #line=[np.array(P1),np.array(P2)]
        rays.append(tuple(P2))

    for i in ray.childs:
        rays.extend(ray2list(i))
    return rays


def evaluate_geometry():
    obj=bpy.data.objects['Cube']
    mw = obj.matrix_world

    mesh = obj.data
    mesh.calc_loop_triangles()
    loop_triangles = mesh.loop_triangles

    S_list=[]
    for tri in loop_triangles:
        tri_center = np.array(mw @ tri.center.to_3d())
        # get global vertice coordinate-vectors from mesh by their index
        tri_vertices = [mw @ mesh.vertices[index].co for index in tri.vertices]
        # convert to arrays
        vertices = [np.array(vertice.to_3d()) for vertice in tri_vertices]

        S = surfaces.Plane(shape=shapes.Triangular((vertices[0],
                                                 vertices[1],
                                                 vertices[2])))
        S_list.append(S)
    Cube = Component(surflist=S_list, material=1.5)
    S=System(complist=[Cube], n=1.4)

    return S


def lightsource(aperture, dir, halfangle = 20., type=None):
    "angles have to be entered in degrees"
    phi = np.random.uniform(0., 2*np.pi)
    angle = halfangle/360*(2*np.pi)
    if not type:
        teta = np.random.uniform(0., angle)
    try:
        if len(aperture)


def trace_rays(system, list_of_rays):
    for ray in list_of_rays:
        system.ray_add(ray)
    system.propagate()
