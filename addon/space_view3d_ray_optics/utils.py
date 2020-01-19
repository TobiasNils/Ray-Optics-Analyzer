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

import bpy
import pyoptools.raytrace.ray as rays
import pyoptools.misc.cmisc as cmisc
from pyoptools.misc.pmisc import wavelength2RGB
import chaospy

from . import draw

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def console_namespace():
    import console_python
    get_consoles = console_python.get_console
    consoles = getattr(get_consoles, "consoles", None)
    if consoles:
        for console, stdout, stderr in get_consoles.consoles.values():
            return console.locals
    return {}



class VarStates:

    @staticmethod
    def store_states():
        # Store the display states, called upon unregister the Add-on
        # This is useful when you press F8 to reload the Addons.
        # Then this function preserves the display states of the
        # console variables.
        state_props = bpy.context.window_manager.RayOpticsStatePropList
        variables = evaluate_geometry()
        for key in variables.keys():
            if key and key not in state_props:
                prop = state_props.add()
                prop.name = key
                # prop.value = variables[key]
                # prop.state = [True, False]

    @staticmethod
    def get_index(key):
        index = bpy.context.window_manager.RayOpticsStatePropList.find(key)
        return index

    @staticmethod
    def delete(key):
        state_props = bpy.context.window_manager.RayOpticsStatePropList
        index = state_props.find(key)
        if index != -1:
            state_props.remove(index)

    @staticmethod
    def toggle_display_state(key):
        state_props = bpy.context.window_manager.RayOpticsStatePropList
        if key in state_props:
            state_props[key].state[0] = not state_props[key].state[0]
        else:
            print("Odd: Can not find key %s in RayOpticsStateProps" % (key))

    @staticmethod
    def toggle_lock_state(key):
        state_props = bpy.context.window_manager.RayOpticsStatePropList
        if key in state_props:
            state_props[key].state[1] = not state_props[key].state[1]
        else:
            print("Odd: Can not find key %s in RayOpticsStateProps" % (key))


def evaluate_geometry():
    import numpy as np
    from pyoptools.raytrace.component import Component
    from pyoptools.raytrace.system import System
    import pyoptools.raytrace.surface as surfaces
    import pyoptools.raytrace.shape as shapes

    # locals = console_namespace()

    lightsources = []
    components = []
    apertures=[]
    for obj_key in bpy.data.objects.keys():
        obj=bpy.data.objects[obj_key]
        if not 'source' in obj.data.name:
            if 'n' in obj.keys():
                # for this to be true, n must be defined as custom property in Blender "Object Properties" <- NOT Mesh Properties
                mw = obj.matrix_world
                if 'reflectivity' in obj.keys():
                    reflectivity = obj['reflectivity']
                else:
                    reflectivity = 0
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

                    S = surfaces.Plane(reflectivity=reflectivity, shape=shapes.Triangular((vertices[0],
                                                             vertices[1],
                                                             vertices[2])))
                    S_list.append(S)
                C = Component(surflist=S_list, material=np.float(obj['n']))
                components.append(C)
            else:
                print(obj.name, 'not included in ray-trace; no refractive index defined.')
        else:
            mw = obj.matrix_world
            mesh = obj.data
            mesh.calc_loop_triangles()
            loop_triangles = mesh.loop_triangles

            for tri in loop_triangles:
                tri_vertices = [mw @ mesh.vertices[index].co for index in tri.vertices]
                # convert to arrays
                vertices = [np.array(vertice.to_3d()) for vertice in tri_vertices]
                apertures.append(vertices)
    try:
        Sys=System(complist=components, n=bpy.data.worlds['World']['n'])
    except KeyError:
        Sys=System(complist=components, n=1.0)

    # pass objects to python console namespace for further processing

    return Sys, apertures

def lightsource(aperture, dist=None):
    """aperture is a defined by a mesh or part of a mesh from Blender.
       Half-angle has to be entered in degrees"""
    import chaospy, numpy as np
    settings = bpy.context.window_manager.RayOpticsProp
    halfangle = settings.halfangle
    n_rays = settings.n_rays
    invert_direction = settings.invert_direction

    list_of_rays = []
    # aperture = [np.array(point) for point in aperture]
    areas = []
    for vertice_group in aperture:
        A,B,C = vertice_group
        v0=C-A
        v1=B-A
        v1xv0=np.cross(v1,v0)
        area_i = .5*sum([v1xv0[i]**2 for i in range(len(v1xv0))] )
        areas.append(area_i)

    for i,vertice_group in enumerate(aperture):
        # adapt ray number to relative area of triangle to total source area
        ni_rays = int(n_rays*areas[i]/sum(areas))

        A,B,C = vertice_group
        v0=C-A
        v1=B-A
        N = np.cross(v0,v1)
        # normalize N
        N_ = N/np.sqrt(sum([N[i]**2 for i in range(len(N))]))
        if invert_direction:
            N_ = -N_

        # normalize an in-plane vector
        v0_ = v0/np.sqrt(sum([v0[i]**2 for i in range(len(v0))]))

        phi_dist = np.around(chaospy.Uniform(0., 2*np.pi).sample(n_rays), 10)
        if not dist:
            # type specifies the distribution to use for random polar angle
            # convert half-angle to radians
            dist = chaospy.Uniform
            args=(0., halfangle/360*(2*np.pi))
        elif dist==chaospy.Normal:
            # let half-angle coincide with 3sigma
            args=(0., halfangle/360*(2*np.pi)/3)
        theta_dist = np.around(dist(*args).sample(n_rays), 10)

        for i in range(ni_rays):
            # create random anker point on plane
            while True:
                # create random point on plane as ray origin
                P = A+np.random.uniform(0,1)*v0 + np.random.uniform(0,1)*v1
                # check if it is in the triangle
                v2=P-A
                dot00=np.dot(v0,v0)
                dot01=np.dot(v0,v1)
                dot02=np.dot(v0,v2)
                dot11=np.dot(v1,v1)
                dot12=np.dot(v1,v2)

                invDenom=1./(dot00 * dot11 - dot01 * dot01)

                u = (dot11 * dot02 - dot01 * dot12) * invDenom
                v = (dot00 * dot12 - dot01 * dot02) * invDenom

                # repeat until in triangle
                if (u > 0) and (v > 0) and (u + v < 1):
                    break
            # create random azimuth angle
            phi = phi_dist[i]
            theta = theta_dist[i]
            # calculate polar angle transform
            v  = np.cos(theta)*N_ + np.sin(theta)*v0_
            # rotate around azimuth angle phi, using Rodrigues' rotation formula
            v = v*np.cos(phi) + np.cross(N_, v)*np.sin(phi) + N_*np.dot(N_,v)*(1-np.cos(phi))
            R = rays.Ray(pos=P, dir=v, intensity=1.0*np.cos(theta))
            list_of_rays.append(R)

    return list_of_rays



def trace_rays(system):
    # for ray in list_of_rays:
    #     system.ray_add(ray)
    print('... propagating rays ...')
    system.propagate()
    print('Ray tracing finished.')


def get_propagated_rays():
    vars = bpy.context.window_manager.RayOpticsVars['items']
    propagated = vars['OpticalSystem']._p_rays
    rays=[]
    for ray in propagated:
        ray = AttrDict({'wavelength':ray.wavelength,
                        'path':[AttrDict({'vertex':ray.pos,
                                            'intensity':ray.intensity}),
                                *ray2list(ray)]})
        rays.append(ray)
    return rays

def ray2list(ray):
    ray_segments=[]

    P1 = ray.pos
    if len(ray.childs) > 0:
        P2 = ray.childs[0].pos
    else:
        P2 = P1 + 10. * ray.dir

    if ray.intensity >= 0.1:
        segment = AttrDict({'vertex':tuple(P2), 'intensity':ray.intensity})
        ray_segments.append(segment)

    for i in ray.childs:
        ray_segments.extend(ray2list(i))

    return ray_segments
