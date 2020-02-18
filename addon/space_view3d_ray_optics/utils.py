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
import chaospy, random

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
        Sys, sources = evaluate_geometry()
        variables = {'OpticalSystem':Sys, 'LightSources':sources}
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

    locals = console_namespace()

    lightsources = []
    apertures=[]
    components = {}
    for obj_key in bpy.data.objects.keys():
        obj=bpy.data.objects[obj_key]
        if obj.visible_get():
            if not 'source' in obj.data.name:
                # for this to be true, n must be defined as custom property in Blender "Object Properties" <- NOT Mesh Properties
                mw = obj.matrix_world
                if 'reflectivity' in obj.keys():
                    reflectivity = obj['reflectivity']
                else:
                    reflectivity = 0
                mesh = obj.data
                # mesh.calc_loop_triangles()
                # loop_triangles = mesh.loop_triangles
                S_list={}
                for poly in mesh.polygons:
                    poly_center = np.array(mw @ poly.center.to_3d())
                    # get global vertice coordinate-vectors from mesh by their index
                    poly_vertices = [mw @ mesh.vertices[index].co for index in poly.vertices]
                    # convert to arrays
                    vertices = [np.array(vertice.to_3d()) for vertice in poly_vertices]

                    S = surfaces.Plane(reflectivity=reflectivity, shape=shapes.Polygon(tuple(vertices)))
                    S_list[poly.index] = S
                try:
                    C = Component(surflist=S_list, material=np.float(obj['n']))
                except KeyError:
                    C = Component(surflist=S_list,
                                material = bpy.data.worlds['World']['n'])
                components[mesh.name]=C
            else:
                mw = obj.matrix_world
                mesh = obj.data
                # mesh.calc_loop_triangles()
                # loop_triangles = mesh.loop_triangles

                for poly in mesh.polygons:
                    poly_vertices = [mw @ mesh.vertices[index].co for index in poly.vertices]
                    # convert to arrays
                    vertices = [np.array(vertice.to_3d()) for vertice in poly_vertices]
                    apertures.append(vertices)
        else:
            print(obj.name, 'invisible -> not included in ray-trace.')

    try:
        Sys=System(complist=components, n=bpy.data.worlds['World']['n'],
                    scene=bpy.context.scene)
    except KeyError:
        Sys=System(complist=components, n=1.0, scene=bpy.context.scene)
    # pass objects to python console namespace for further processing
    locals['Sys'] = Sys # <- not working

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
        A = vertice_group[0]
        # split polygon in triangles with A as reference point
        remaining = len(vertice_group[1:])
        hit=False
        for i in range(remaining-1):

            B=np.array(vertice_group[i+1])
            C=np.array(vertice_group[i+2])

            v0=C-A
            v1=B-A

            v1xv0=np.cross(v1,v0)
            area_i = .5*sum([v1xv0[i]**2 for i in range(len(v1xv0))] )
            areas.append(area_i)

    for i,vertice_group in enumerate(aperture):
        # adapt ray number to relative area of triangle to total source area
        ni_rays = int(n_rays*areas[i]/sum(areas))

        A = vertice_group[0]
        # split polygon in triangles with A as reference point
        remaining = len(vertice_group[1:])
        hit=False
        for j in range(remaining-1):

            B=np.array(vertice_group[j+1])
            C=np.array(vertice_group[j+2])

            v0=C-A
            v1=B-A
            N = np.cross(v0,v1)
            # normalize N
            N_ = N/np.sqrt(sum([N[i]**2 for i in range(len(N))]))
            if invert_direction:
                N_ = -N_

            # normalize an in-plane vector
            v0_ = v0/np.sqrt(sum([v0[i]**2 for i in range(len(v0))]))

            phi_dist = np.around(chaospy.Uniform(0., 2*np.pi).sample(ni_rays), 10)
            if not dist:
                # type specifies the distribution to use for random polar angle
                # convert half-angle to radians
                dist = chaospy.Uniform
                args=(0., halfangle/360*(2*np.pi))
            elif dist==chaospy.Normal:
                # let half-angle coincide with 3sigma
                args=(0., halfangle/360*(2*np.pi)/3)
            theta_dist = np.around(dist(*args).sample(ni_rays), 10)

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
                R = rays.Ray(pos=P+1e-7*v, dir=v, intensity=1.0*np.cos(theta))
                list_of_rays.append(R)

    return list_of_rays

from copy import deepcopy, copy
from multiprocess import Pool

# def doWork(system):
#     while len(system._np_rays)>0:
#         ri = system._np_rays.pop(0)
#         system.propagate_ray(ri)
#         system._p_rays.append(ri)
#     return system


def trace_rays(system):
    import time
    import numpy as np
    def propagate_ray(ri):
        # create a hitlist the function can return in order to avoid shared memory objects which break multiprocessing
        surf_hits=[]
        # for every hit, the tuple of (optsurf.id, pi, ri) will be appended
        # optsurf.id yields a list of nested keys: e.g. ['Sphere', 0]
        if ri.n== None:
            ri.n=system.n

        hit = scene.ray_cast(view, ri.pos, ri.dir)
        if not hit[0]:
            return surf_hits
        # create dict to pass information down the tree (avoid recalculating)
        hit = {'bool':hit[0],
                'location':np.array(hit[1].to_3d()),
                'normal':np.array(hit[2].to_3d()),
                'surface':hit[3],
                'object':hit[4].data.name,
                }

        # Check if ray is inside of an object (if yes, adjust n)
        # get angle between surface normal (facing outside of object) and ray
        angle = np.arccos(np.dot(ri.dir,hit['normal']))
        if angle<=np.pi/2:
            # ray is inside
            ri.n = system.complist[hit['object']].n(ri.wavelength)
        else:
            ri.n = system.n
        # add hit to list of hit optical surfaces for multiprocessing
        surf_hits.append(([hit['object'], hit['surface']],
                        (hit['location'], ri)))

        # check medium after interface
        potential_next_hit = scene.ray_cast(view,
                                            hit['location']+ri.dir*1e-6, ri.dir)
        if potential_next_hit[0]:
            angle = np.arccos(np.dot(ri.dir,
                                    np.array(potential_next_hit[2].to_3d())))
            if angle<=np.pi/2:
                # ray is inside another object
                n_medium = system.complist[
                            potential_next_hit[4].data.name].n(ri.wavelength)
            else:
                n_medium = system.n
        else:
            n_medium = system.n

        # Check surface effects on ray on component basis
        C = system.complist[hit['object']]
        ri_n = C.propagate(ri, n_medium, hit)
        for i in ri_n:
            # put the rays in the childs list
            ri.add_child(i)
        # Propagate childs
        for i in ri.get_final_rays():
            if (i!=ri):
                if i.intensity>0:
                    surf_hits_i = propagate_ray(i)
                    surf_hits=surf_hits+surf_hits_i
            else:
                raise Exception("Error, a ray can not be parent and child at the same time")

        return surf_hits

    def doWork(ri):
        surf_hits = propagate_ray(ri)
        return ri, surf_hits


    settings = bpy.context.window_manager.RayOpticsProp
    scene = bpy.context.scene
    view = scene.view_layers['View Layer']

    system.update_ids()
    #mark the start time
    startTime = time.time()
    print('... propagating rays ...')

    result=[]
    while len(system._np_rays)>0:
        ri = system._np_rays.pop(0)
        entry = doWork(ri)
        result.append(entry)
    # with Pool(settings.processes) as pool:
    #     result = pool.map(doWork, system._np_rays)
        # every entry in results contains the propagated ray ri and a corresponding list of surface hits. The latter needs to be transferred to the respective optical surfaces
    for entry in result:
        # add the propagated ray to the corresponding list
        system._p_rays.append(entry[0])
        # cycle through list of surface hits and update corresponding surfaces
        for hit in entry[1]:
            optsurf = system.get_surface(hit[0])
            #print(optsurf)
            optsurf._hit_list.append(hit[1])
            #print(hit[1])
    del result
    # system.propagate(settings.processes)

    print('Ray tracing finished.')
    #mark the end time
    endTime = time.time()
    #calculate the total time it took to complete the work
    workTime =  endTime - startTime
    #print results
    print ("The job took " + str(workTime) + " seconds to complete")

def get_propagated_rays(number):
    vars = bpy.context.window_manager.RayOpticsVars['items']

    rays=[]
    total_rays = vars['OpticalSystem']._p_rays
    if len(total_rays)<=number:
        sample = total_rays
    else:
        sample = random.sample(total_rays, number)
    for ray in sample:
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
