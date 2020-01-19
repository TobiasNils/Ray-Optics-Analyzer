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
import blf
import gpu
import bgl
from gpu_extras.batch import batch_for_shader
from pyoptools.misc.pmisc import wavelength2RGB
from . import utils
# from mathutils import Vector

SpaceView3D = bpy.types.SpaceView3D
callback_handle = None

if not bpy.app.background:
    single_color_shader = gpu.shader.from_builtin('3D_UNIFORM_COLOR')
    flat_color_shader = gpu.shader.from_builtin('3D_FLAT_COLOR')
else:
  single_color_shader = None
  flat_color_shader = None


# COLOR_POINT = (1.0, 0.0, 1.0, 1)
# COLOR_LINE = (0.5, 0.5, 1, 1)
# COLOR_BOUNDING_BOX = (1.0, 1.0, 1.0, 1.0)

def tag_redraw_areas():
    context = bpy.context

    # Py cant access notifers
    for window in context.window_manager.windows:
        for area in window.screen.areas:
            if area.type in ['VIEW_3D', 'PROPERTIES']:
                area.tag_redraw()


def callback_enable():
    # if callback_handle:
    #     return

    # handle_pixel = SpaceView3D.draw_handler_add(draw_callback_px, (), 'WINDOW', 'POST_PIXEL')
    handle_view = SpaceView3D.draw_handler_add(draw_callback_view, (), 'WINDOW', 'POST_VIEW')
    callback_handle = handle_view

    tag_redraw_areas()


def callback_disable():
    # if not callback_handle:
    #     return

    handle_view = callback_handle
    # SpaceView3D.draw_handler_remove(handle_pixel, 'WINDOW')
    SpaceView3D.draw_handler_remove(handle_view, 'WINDOW')
    callback_handle = None

    tag_redraw_areas()


def draw_ray(coords, colors, type='LINE_STRIP',shader='flat'):
    if shader=='single':
        shader = single_color_shader
        batch = batch_for_shader(shader, type, {"pos": coords})
        shader.bind()
        shader.uniform_float("color", colors)
    else:
        shader = flat_color_shader
        batch = batch_for_shader(shader, type, {"pos": coords, 'color':colors})
    # def draw():

    bgl.glEnable(bgl.GL_BLEND)
    bgl.glEnable(bgl.GL_DEPTH_TEST)
    batch.draw(shader)

    # bpy.types.SpaceView3D.draw_handler_add(draw, (), 'WINDOW', 'POST_VIEW')



def render_rays():
    """Convencience function combining ray-trace and drawing of the result"""
    flat_color_shader = gpu.shader.from_builtin('3D_FLAT_COLOR')
    settings = bpy.context.window_manager.RayOpticsProp
    if not settings.ray_hide:
        return

    for ray in utils.get_propagated_rays():
        color = wavelength2RGB(ray.wavelength)
        path, colors=[], []
        for segment in ray.path:
            colors.append(tuple([*color, segment.intensity-.4]))
            path.append(segment.vertex)
        draw_ray(path, colors)
    # print(len(path), len(colors))
    bgl.glDisable(bgl.GL_BLEND)

def render_source():
    single_color_shader = gpu.shader.from_builtin('3D_UNIFORM_COLOR')
    vars = bpy.context.window_manager.RayOpticsVars['items']
    settings = bpy.context.window_manager.RayOpticsProp
    if settings.show_sources and len(vars['OpticalSystem']._p_rays)==0:
        for ray in vars['OpticalSystem']._np_rays:
            color = (*wavelength2RGB(ray.wavelength), .5)
            draw_ray([tuple(ray.pos), tuple(ray.pos+10*ray.dir)],color,
                    type='LINES',shader='single')



def draw_callback_view():
    settings = bpy.context.window_manager.RayOpticsProp

    render_source()
    render_rays()
