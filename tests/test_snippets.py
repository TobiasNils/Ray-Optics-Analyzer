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

polygon = [np.array([0.,0.]),np.array([100.,0.]),np.array([0.,100.])]
path = Path(polygon)

p =np.array([102.,0])
path.contains_point(p, radius=10e-6)
path.get_extents().xmin




for i, vertex in enumerate(polygon):
    vec = polygon[i-1]-vertex

from pyoptools.raytrace.surface import Surface

import pyoptools.raytrace.surface as surfaces
import pyoptools.raytrace.shape as shapes
surfaces.Plane(shape=shapes.Polygon(polygon))
Path(((0,0),(0,100),(100,0)))
type(shapes.Polygon)

import numpy as np
t = np.array


class test():
    def __init__(self):
        self.status = 'initalised.'
        print(self.status)

t = test()
import copy



t2 = copy.copy(t)

t2.status = 'test'
t2.status

t.status
t.status = 'test2'

t2.t


a = ([1,2], [2,3])
a[0]

from multiprocess import Pool




p = Pool(2)
help(p.map)
