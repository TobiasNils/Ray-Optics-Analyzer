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

bl_info = {
    "name": "Ray-Optics",
    "author": "Tobias Nils Ackermann",
    "version": (0, 1,),
    "blender": (2, 80, 1),
    "location": "Properties: Scene > Ray-Optics Control Panel: Menu",
    "description": "Setup Ray-tracing parameters and control render in the 3D view",
    "wiki_url": "",
    "support": "TESTING",
    "category": "3D View",
}


if "bpy" in locals():
    import importlib
    importlib.reload(utils)
    importlib.reload(draw)
else:
    from . import utils
    from . import draw

import bpy
from bpy.types import (
    Operator,
    Panel,
    PropertyGroup,
    UIList,
)
from bpy.props import (
    StringProperty,
    BoolProperty,
    BoolVectorProperty,
    FloatProperty,
    IntProperty,
    PointerProperty,
    CollectionProperty,
    EnumProperty,
)


class PanelConsoleVars(Panel):
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = 'scene'
    bl_label = "RayOptics Console"
    bl_idname = "RayOptics_PT_panel_console_vars"
    bl_category = "Ray Optics"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        state_props = bpy.context.window_manager.RayOpticsStatePropList

        if len(state_props) == 0:
            box = layout.box()
            col = box.column(align=True)
            col.label(text="No vars to display")
        else:
            layout.template_list(
                'RayOpticsVarList',
                'RayOpticsStatePropList',
                bpy.context.window_manager,
                'RayOpticsStatePropList',
                bpy.context.window_manager.RayOpticsProp,
                'index',
                rows=10
            )
        col = layout.column()
        col.operator("rayoptics.geometry_to_console")
        col.prop(bpy.context.window_manager.RayOpticsProp, "n_rays")
        col.prop(bpy.context.window_manager.RayOpticsProp, "halfangle")
        col.operator("rayoptics.init_light_sources")
        col.prop(bpy.context.window_manager.RayOpticsProp, "invert_direction")
        col.prop(bpy.context.window_manager.RayOpticsProp, "show_sources")
        col.prop(bpy.context.window_manager.RayOpticsProp, "processes")
        col.operator("rayoptics.trace_rays")
        col.prop(bpy.context.window_manager.RayOpticsProp, "render_rays")
        col.prop(bpy.context.window_manager.RayOpticsProp, "ray_hide")

#
# class DeleteVar(Operator):
#     bl_idname = "rayoptics.delete_var"
#     bl_label = "Delete Var"
#     bl_description = "Remove the variable from the Console"
#     bl_options = {'REGISTER'}
#
#     key: StringProperty(name="Key")
#
#     def execute(self, context):
#         locals = utils.console_namespace()
#         utils.VarStates.delete(self.key)
#         del locals[self.key]
#         draw.tag_redraw_areas()
#         return {'FINISHED'}

#
# class ToggleDisplay(Operator):
#     bl_idname = "rayoptics.toggle_display"
#     bl_label = "Hide/Unhide"
#     bl_description = "Change the display state of the var"
#     bl_options = {'REGISTER'}
#
#     key: StringProperty(name="Key")
#
#     def execute(self, context):
#         utils.VarStates.toggle_display_state(self.key)
#         draw.tag_redraw_areas()
#         return {'FINISHED'}


# class ToggleLock(Operator):
#     bl_idname = "rayoptics.toggle_lock"
#     bl_label = "Lock/Unlock"
#     bl_description = "Lock the var from being deleted"
#     bl_options = {'REGISTER'}
#
#     key: StringProperty(name="Key")
#
#     def execute(self, context):
#         utils.VarStates.toggle_lock_state(self.key)
#         draw.tag_redraw_areas()
#         return {'FINISHED'}
#


class Geometry2Console(Operator):
    bl_idname = "rayoptics.geometry_to_console"
    bl_label = "Inspect Geometry"
    bl_description = "Inspect Geometry and pass Optical System information to the Console"
    bl_options = {'REGISTER'}

    def execute(self, context):
        Sys, apertures = utils.evaluate_geometry()
        bpy.context.window_manager.RayOpticsVars['items']['OpticalSystem'] = Sys
        bpy.context.window_manager.RayOpticsVars['items']['apertures'] = apertures
        # locals['apertures'] = apertures
        # draw.tag_redraw_areas()
        return {'FINISHED'}

class InitLightSources(Operator):
    bl_idname = "rayoptics.init_light_sources"
    bl_label = "Initialize lightsources"
    bl_description = 'Create lightsources from any mesh with "source" in its name'
    bl_options = {'REGISTER'}

    def execute(self, context):
        import chaospy
        vars = bpy.context.window_manager.RayOpticsVars['items']
        _np_rays = utils.lightsource(vars['apertures'],
                                    dist=chaospy.Normal)
        vars['OpticalSystem'].reset()
        vars['OpticalSystem']._np_rays = _np_rays
        draw.callback_enable()
        draw.tag_redraw_areas()
        return {'FINISHED'}


class Raytrace(Operator):
    bl_idname = "rayoptics.trace_rays"
    bl_label = "Trace rays"
    bl_description = "Trace rays through Geometry"
    bl_options = {'REGISTER'}

    def execute(self, context):
        vars = bpy.context.window_manager.RayOpticsVars['items']
        utils.trace_rays(vars['OpticalSystem'])
        # locals['apertures'] = apertures
        draw.tag_redraw_areas()
        return {'FINISHED'}
# def menu_func_cleanup(self, context):
#     self.layout.operator("rayoptics.cleanup_console", text="Clear Math Vis")

#
# def console_hook():
#     # utils.VarStates.store_states()
#     draw.tag_redraw_areas()
#     context = bpy.context
#     for window in context.window_manager.windows:
#         window.screen.areas.update()
#
#
# def call_console_hook(self, context):
#     console_hook()


class RayOpticsStateProp(PropertyGroup):
    ktype: StringProperty()
    state: BoolVectorProperty(default=(False, False), size=2)


class RayOpticsVarList(UIList):
    bl_idname = "RayOptics_UL_RayOpticsVarList"

    def draw_item(self,
                  context,
                  layout,
                  data,
                  item,
                  icon,
                  active_data,
                  active_propname
                  ):

        col = layout.column()
        key = item.name
        ktype = item.ktype
        is_visible = item.state[0]
        is_locked = item.state[1]

        row = col.row(align=True)
        row.label(text='%s - %s' % (key, ktype))

        icon = 'RESTRICT_VIEW_OFF' if is_visible else 'RESTRICT_VIEW_ON'
        prop = row.operator("rayoptics.toggle_display", text='', icon=icon, emboss=False)
        prop.key = key

        icon = 'LOCKED' if is_locked else 'UNLOCKED'
        prop = row.operator("rayoptics.toggle_lock", text='', icon=icon, emboss=False)
        prop.key = key

        if is_locked:
            row.label(text='', icon='BLANK1')
        else:
            prop = row.operator("rayoptics.delete_var", text='', icon='X', emboss=False)
            prop.key = key


class RayOptics(PropertyGroup):

    index: IntProperty(
        name="index"
    )
    ray_hide: BoolProperty(
        name="Show rays",
        default=True,
        description="Show the rays propagating the system",
        # update=call_console_hook
    )
    n_rays: IntProperty(
        name="rays to trace",
        min = 0, max=100000,default=1000,
        description="Total number of rays to propagate through the system",
        # update=call_console_hook
    )
    processes: IntProperty(
        name="processes to spawn",
        min = 1, max=256,default=8,
        description="Total number of processes used in ray-propagation calculation",
        # update=call_console_hook
    )
    render_rays: IntProperty(
        name="rendered rays",
        min = 1, max=100000,default=250,
        description="Number of rays to render in 3DView",
        # update=call_console_hook
    )
    halfangle: FloatProperty(
        name="halfangle",
        min=0, max=90.,default=1,
        description=""
    )
    show_sources: BoolProperty(
        name="Inspect Sources",
        default=True,
        description="Draw tentative source representations",
        # update=call_console_hook
    )
    invert_direction: BoolProperty(
        name="Invert source direction",
        default=False,
        description="",
        # update=call_console_hook
    )

classes = (
    PanelConsoleVars,
    # DeleteVar,
    # ToggleDisplay,
    # ToggleLock,
    # ToggleMatrixBBoxDisplay,
    Geometry2Console,
    InitLightSources,
    Raytrace,
    RayOpticsStateProp,
    RayOpticsVarList,
    RayOptics,
)


def register():
    from bpy.utils import register_class

    draw.callback_enable()
    # import console_python
    # console_python.execute.hooks.append((console_hook, ()))
    for cls in classes:
        bpy.utils.register_class(cls)
    bpy.types.WindowManager.RayOpticsVars=EnumProperty(items={})[1]

    bpy.types.WindowManager.RayOpticsProp = PointerProperty(type=RayOptics)
    bpy.types.WindowManager.RayOpticsStatePropList = CollectionProperty(type=RayOpticsStateProp)
    # bpy.types.CONSOLE_MT_console.prepend(menu_func_cleanup)



def unregister():
    from bpy.utils import unregister_class

    draw.callback_disable()

    # import console_python
    # console_python.execute.hooks.remove((console_hook, ()))
    # bpy.types.CONSOLE_MT_console.remove(menu_func_cleanup)
    del bpy.types.WindowManager.RayOpticsProp
    del bpy.types.WindowManager.RayOpticsStatePropList

    for cls in classes:
        unregister_class(cls)
