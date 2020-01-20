import numpy as np
import bpy
# TODO: add "World" cube for limiting the simulation geometry
# TODO: scale script to include all mesh objects within "World"

# objects and respective meshes and faces which should be included in the simulation can be accessed through meshes:
import blf
import bgl

import gpu
from gpu_extras.batch import batch_for_shader
from mathutils import Vector
import time
import chaospy

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.system import System
import pyoptools.raytrace.surface as surfaces
import pyoptools.raytrace.shape as shapes
import pyoptools.raytrace.ray as rays
import pyoptools.misc.cmisc as cmisc
from pyoptools.misc.pmisc import wavelength2RGB


SpaceView3D = bpy.types.SpaceView3D
callback_handle = []

if not bpy.app.background:
    single_color_shader = gpu.shader.from_builtin('3D_UNIFORM_COLOR')
    flat_color_shader = gpu.shader.from_builtin('3D_FLAT_COLOR')
else:
  single_color_shader = None
  flat_color_shader = None


def draw_ray(coords, colors, type='LINE_STRIP'):
    batch = batch_for_shader(flat_color_shader, type, {"pos": coords, 'color':colors})
    def draw():

        bgl.glEnable(bgl.GL_BLEND)
        bgl.glEnable(bgl.GL_DEPTH_TEST)
        batch.draw(flat_color_shader)

    bpy.types.SpaceView3D.draw_handler_add(draw, (), 'WINDOW', 'POST_VIEW')


def ray2list(ray):
    ray_segments=[]

    P1 = ray.pos
    if len(ray.childs) > 0:
        P2 = ray.childs[0].pos
    else:
        P2 = P1 + 10. * ray.dir

    if ray.intensity != 0:
        segment = AttrDict({'vertex':tuple(P2), 'intensity':ray.intensity})
        ray_segments.append(segment)

    for i in ray.childs:
        ray_segments.extend(ray2list(i))

    return ray_segments

def evaluate_geometry():

    lightsources = []
    components = []
    apertures=[]
    for obj_key in bpy.data.objects.keys():
        obj=bpy.data.objects[obj_key]
        if not 'source' in obj.data.name:
            if 'n' in obj.keys():
                # for this to be true, n must be define as custom property in Blender "Object Properties" <- NOT Mesh Properties
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

                    S = surfaces.Plane(reflectivity=0, shape=shapes.Triangular((vertices[0],
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
    print('Geometries evaluated.')
    return Sys, apertures

def trace_rays(system, list_of_rays):
    for ray in list_of_rays:
        system.ray_add(ray)
    print('... propagating rays ...')
    system.propagate()
    print('Ray tracing finished.')

def get_propagated_rays(system):
    propagated = system._p_rays
    rays=[]
    for ray in propagated:
        ray = AttrDict({'wavelength':ray.wavelength,
                        'path':[AttrDict({'vertex':ray.pos,
                                            'intensity':ray.intensity}),
                                *ray2list(ray)]})
        # paths = [[tuple(r.pos)]+ray2list(r) for r in propagated]
        rays.append(ray)
    return rays

def render_rays(S, rays):
    """Convencience function combining ray-trace and drawing of the result"""
    trace_rays(S, rays)
    for ray in get_propagated_rays(S):
        color = wavelength2RGB(ray.wavelength)
        path, colors=[], []
        for segment in ray.path:
            colors.append(tuple([*color, segment.intensity]))
            path.append(segment.vertex)
        draw_ray(path, colors)
    print(len(path), len(colors))
    bgl.glDisable(bgl.GL_BLEND)
    # callback_enable()

# def live_render(Ri):
#     T = threading.Thread(target=render_rays, args=[[Ri]])
#     T.start()
#     return T


def lightsource(aperture, inspect=False, invert_direction=False, n_rays=10, halfangle = 0., dist=None):
    """aperture is a defined by a mesh or part of a mesh from Blender.
       Half-angle has to be entered in degrees"""

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
            if inspect:
                draw_ray([tuple(P), tuple(P+10*v)])

    # print(list_of_rays)
    return list_of_rays
# R=rays.Ray(pos=(0,0,-2),dir=(0,.3,1))
# R2=rays.Ray(pos=(0,0,-2),dir=(0.3,.3,1))
# R3=rays.Ray(pos=(0,0,-2),dir=(0.3,-.3,1))
#
# Ri = [R, R2, R3]

if __name__=='__main__':

    # aperture = [(-1,0,0), (1,0,0),(0,0,1)]
    Sys, aperture = evaluate_geometry()
    # set refractive index of "world"
    Sys.n = 1.0
    # define light-source
    np_rays = lightsource(aperture, invert_direction=True, n_rays=100, halfangle=45., dist=chaospy.Normal)

    render_rays(Sys, np_rays)

    # callback_enable()


    # rays = [tuple(R.pos)]+ray2list(R)
    # draw_ray(rays)
    #
    # rays2 = [tuple(R2.pos)]+ray2list(R2)
    # draw_ray(rays2)



    # draw_ray([tuple(R.pos)]+ray2list(R),'LINES')
    # draw_ray(ray2list2(R),'LINES')

    # draw_ray(rays, 'LINES')
        #=======================================================================

        # Pi = LinePlaneCollision(N, P, R.dir, R.pos)
        # Pi = S.intersection(R)
        #=======================================================================



        #=======================================================================
        #
        # now check if the point of intersection lies in the surface segment
        # if S.shape.fhit(Pi[0],Pi[1],Pi[2]):
        #     print ("intersection at", Pi)
        #     draw_ray([Pi])

    bpy.context.view_layer.update()
    # bmesh.free()
