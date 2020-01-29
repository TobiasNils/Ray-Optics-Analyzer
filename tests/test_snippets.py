rays=[1,2,3,4,5,6,7,8,9,10,11,12,13,14]
#
# new_rays = [rays.pop(0) for i in range(len(rays)-5)]
#
import time
def doWork(ri):
    time.sleep(1)
    # print(ri)
    return ri

from multiprocess import Pool, cpu_count
p = Pool(14)
# import numpy as np
#mark the start time
startTime = time.time()
result = p.map(doWork, rays)
#mark the end time
endTime = time.time()
#calculate the total time it took to complete the work
workTime =  endTime - startTime
#print results
print ("The job took " + str(workTime) + " seconds to complete")
result

liste =[]
test =[2]
liste.append(*test)
liste


import copy
S = {'test':[]}
class SystemImage():
    def __init__(self, System):
        self.S = copy.deepcopy(System)
    def add(self, arg):
        self.S['test'].append(arg)

S2 = SystemImage(S)
S2.add(5)

S2.S['test']
S['test']

import random
import chaospy
import numpy as np
import matplotlib.pyplot as plt
random.sample(rays, 3)

normal  = chaospy.Normal(.5, np.sqrt(.01))
normal
p = plt.hist(np.around(normal.sample(10000),10),100, density=True)


from matplotlib.path import Path


def normal():
    """Method that returns the normal to the surface in global coordinates
    """
    v1 = polygon[0]
    v2 = polygon[1]
    v3 = polygon[2]
    # Take two in-plane vectors
    p1 = v2-v1
    p2 = v3-v1
    # Get the vector perpendicular to the plane p1 and p2 span by cross-product
    N_ = np.cross(p1,p2)
    # normalize
    N_ = N_/np.linalg.norm(N_)

    N_=np.array(N_).astype(np.float64)
    return (N_)

def inplane_vectors():
    """Method that returns the normal to the surface in global coordinates
    """
    v1=polygon[0]
    v2=polygon[1]
    # Take an in-plane vectors and the normal
    u = v2-v1
    u = u/np.linalg.norm(u)
    N_ = normal()
    # Get the already normalized vector perpendicular by cross-product
    v = np.cross(u,N_)

    # create empty list to append vertices in uv-coordinates
    uv_poly = []
    # add v1 as local (0,0) and skip it in subsequent loop
    uv_poly.append(np.array([0.,0.]))
    # cdef np.ndarray p,retval
    # cdef double norm_p,theta,u_fac,v_fac
    for p in polygon[1:]:
        norm_p = np.sqrt(sum([p[i]**2 for i in range(3)]))
        theta = np.arccos(np.dot(p/norm_p, u))
        u_fac = norm_p*np.cos(theta)
        v_fac = norm_p*np.sin(theta)

        # the plane point can now be expressed as
        # p = v1 + u_fac*u + v_fac*v
        # local coordinates are therefore
        retval = np.array([u_fac, v_fac])
        uv_poly.append(retval)
    assert len(uv_poly)==len(polygon)
    u=np.array(u).astype(np.float64)
    v=np.array(v).astype(np.float64)
    return (u, v, tuple(uv_poly))

def hit(p, u,v, uv_poly):
    """Method  that returns True if a p=(x,y,z) point is inside the triangle,
    if not it returns False.
    taken from http://www.blackpawn.com/texts/pointinpoly/default.html
    """
    local_origin = polygon[0]
    p = p-local_origin
    norm_p = np.linalg.norm(p)
    p_=p/norm_p
    theta = np.arccos(np.dot(p_, u))
    u_fac = norm_p*np.cos(theta)
    v_fac = norm_p*np.sin(theta)

    retval = [u_fac,v_fac]

    path = Path(list(uv_poly))
    return path.contains_point(retval, radius=10e-10)

def _intersection(pos, dir):
    """Returns the intersection point between a ray and an the plane in global 3D coordinates

    """
    # cdef np.ndarray w,retval
    epsilon=1e-6
    # get surface data
    planePos=polygon[0]
    planeNormal=normal()

    rayDirection = np.array(dir)/np.linalg.norm(np.array(dir))
    rayPos = np.array(pos)/np.linalg.norm(np.array(pos))

    #if dot(N_,A.dir) ==0 : return inf_vect
    ndotu=np.dot(planeNormal,rayDirection)
    if np.abs(ndotu) < epsilon: return np.array([np.inf for i in range(3)]) # ray parallel or within the plane
    w=rayPos-planePoint
    si=np.dot(planeNormal,w)
    si = -si/ndotu
    retval=w+si*rayDirection+planePoint

    # the following should be true for all 3 coordinates:
    # retval[0] = rayPos[0] + fac*rayDirection[0]
    for i in range(3):
        fac = (retval[i]-rayPos[i])/rayDirection[i]
        print(fac)
    # if fac negative, the intersection lies behind the ray -> return inf
    if fac<0: return np.array([np.inf for i in range(3)])
    else: return retval



polygon = [np.array([0.,0.,0.]),np.array([0.,5.,10.]),np.array([10.,5.,10.])]
plt.plot(polygon, 'bv')

u, v, uv_poly = inplane_vectors()

u
v
uv_poly
plt.scatter([vec[0] for vec in uv_poly],[vec[1] for vec in uv_poly] )

p = _intersection([10.,4.99999999999,0.], [0.,0.,1.])
p

hit(p, u, v, uv_poly)
