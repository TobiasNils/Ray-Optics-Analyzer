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


from numpy import sqrt, arccos, pi, cos, sin, array, dot
def uv_convert(point, refPoint):
    p = point-refPoint
    # calc length of relative coord
    norm_p = sqrt(sum([p[i]**2 for i in range(3)]))
    # get angle theta between p and u
    theta = arccos(dot(p/norm_p, u))
    # Problem: dot product is symmetric around pi. -> compare p with v to distinguish angles < pi from those > pi.
    phi = arccos(dot(p/norm_p, v))
    if phi <= pi/2:
        # 0 <= theta  <= pi
        pass
    else:
        # pi < theta < 2pi
        theta = theta + pi

    # project p on u and v using trigonometric relations
    u_fac = norm_p*cos(theta)
    v_fac = norm_p*sin(theta)


    # the plane point can now be expressed as
    # p = refPoint + u_fac*u + v_fac*v
    # local coordinates are therefore
    retval = array([u_fac, v_fac])
    return retval


from code import InteractiveConsole

if __name__=='__main__':

    import bpy
    locals = console_namespace()


    #for obj in bpy.data.meshes:

        #obj.calc_loop_triangles()
        #locals[obj.name] = obj

    locals['variables'] = bpy.context.window_manager.RayOpticsVars['items']

    Sys = bpy.context.window_manager.RayOpticsVars['items']['OpticalSystem']
    locals['Sys'] = bpy.context.window_manager.RayOpticsVars['items']['OpticalSystem']


    #triangles = [tri for tri in obj.loop_triangles]
    #locals['tri'] = triangles[1]

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    
    from scipy.ndimage.filters import gaussian_filter
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.patches import Polygon
    from scipy.spatial import ConvexHull

    unit_name = bpy.context.scene.unit_settings.length_unit
    if unit_name=='METERS':
        unit={'symbol':'m', 'scale':1.}
    elif unit_name=='CENTIMETERS':
        unit = {'symbol':'cm', 'scale':1e-2}
    elif unit_name=='MILLIMETERS':
        unit={'symbol':'mm', 'scale':1e-3}
    elif unit_name=='MICROMETERS':
        unit={'symbol':'µm', 'scale':1e-6}

    unit = AttrDict(unit)

    mesh = bpy.context.active_object.data
    m = bpy.context.active_object.matrix_world

    bpy.ops.object.mode_set(mode='OBJECT')
    for poly in mesh.polygons:
        if poly.select:
            print('selected', poly.index)

            poly_center = np.array(m @ poly.center.to_3d())

            # now get corresponding optical surface from System
            optsurf = Sys[bpy.context.active_object.data.name][poly.index]
            # use the function inplane_vectors to get u,v and uv-coords from the optical surface
            u,v,uv_poly = optsurf.inplane_vectors(poly_center)

            outline= ConvexHull(uv_poly)
            outline = Polygon([uv_poly[i] for i in outline.vertices], closed=True, edgecolor='w', facecolor='none', alpha=.3)

            print(len(uv_poly))
            uv_poly = [[uv[0] for uv in uv_poly],
                        [uv[1] for uv in uv_poly]]


            #transform hit coords to uv-system
            uv_hits = [uv_convert(hit[0], poly_center) for hit in optsurf.hit_list]

            hitsU = [hit[0] for hit in uv_hits]
            hitsV = [hit[1] for hit in uv_hits]
            # add polygon vertices with weights=0
            hitsU = hitsU+uv_poly[0]
            hitsV = hitsV + uv_poly[1]


            # todo: adjust bins according to dimensions
            ui, vi = uv_poly[0],uv_poly[1]
            extent = [min(ui),max(ui),min(vi),max(vi)]
            udim = abs(max(ui)-min(ui))
            vdim = abs(max(vi)-min(vi))
            v_bins = 1e3
            bin_area = (vdim*unit.scale/v_bins)**2
            print('area', poly.area)
            # define weights  for 1 W = 1000 mW total lightsource power distributed evenly on _p_rays and per bin_area
            uv_weights = [hit[1].intensity*1/len(Sys._p_rays)/bin_area for hit in optsurf.hit_list]
            # add "hits" with zero weight on corners of polygon to extend the histogram area
            uv_weights = uv_weights+ [0. for i in range(len(uv_poly[0]))]

            heatmap, xedges, yedges = np.histogram2d(hitsU, hitsV, bins=(udim*v_bins/vdim,v_bins), weights=uv_weights)
            heatmap = gaussian_filter(heatmap, sigma=8)

            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]]



            fig, ax = plt.subplots()

            map = ax.imshow(heatmap.T, extent=extent,origin='lower', cmap=cm.jet)
            #ax.plot(uv_poly[0],uv_poly[1], 'wo', markersize=5)
            fig.colorbar(map, label='$W/m^2$')
            ax.add_patch(outline)

            incident_flux = sum([hit[1].intensity*1/len(Sys._p_rays) for hit in optsurf.hit_list])

            ax.set_xlabel('u-dimension [%s]\n\n Total Incident Flux %f W, %i·10$^3$ Incident Rays'%(unit.symbol, incident_flux, len(Sys._p_rays)/1000))
            ax.set_ylabel('v-dimension [%s]'%unit.symbol)
            ax.grid(True, color='w', linestyle='--', alpha=.1)
            ax.set_title('Irradiance Map for incident Flux')
            ax.autoscale(tight=False)
            ax.set(aspect=1)
            fig.tight_layout()

            plt.show()
            #print('uv-coords of polygon', uv_poly)
            #print('n hits', len(optsurf.hit_list))

            #x = [vec[0] for vec in optsurf.shape.coord]
            #y = [vec[1] for vec in optsurf.shape.coord]
            #z = [vec[2] for vec in optsurf.shape.coord]


            #fig = plt.figure()
            #ax = Axes3D(fig)
            #ax.add_collection3d(Poly3DCollection([zip(x, y, z)]))

    bpy.ops.object.mode_set(mode='EDIT')
