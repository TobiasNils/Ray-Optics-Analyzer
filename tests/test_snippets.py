import numpy as np
from scipy.spatial import ConvexHull
rays=[1,2,3,4,5,6,7,8,20,9,10,11,12,13,14,14,15,15]
2*1e-3

import matplotlib.pyplot as plt
help(np.histogram2d)

# new_rays = [rays.pop(0) for i in range(len(rays)-5)]
class Obj():
    values = []
    def f(self, ri):
        time.sleep(1)
        self.values.append(ri)
        # print(ri)
        return ri, ri

import time
O = Obj()

from multiprocess import Pool, cpu_count

p = Pool(14)

result = p.map(O.f, rays)
result
O.values = [result[1] for r in result]
O.values

list(result)

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
random.sample(rays, 3)


from multiprocess import Pool, cpu_count
# cpus = cpu_count()
with Pool(processes) as pool:

    propagating_rays = []
    for i in range(len(self._np_rays)):
        ri=self._np_rays.pop(0)
        propagating_rays.append(ri)
    result = pool.map(self.doWork, propagating_rays)

self._p_rays = result
del propagating_rays
# ray_sets)
# while len(self._np_rays)>0:
#     ri=self._np_rays.pop(0)
#     self.propagate_ray(ri)
#     self._p_rays.append(ri)

def doWork(self, ri):
    # for ri in ray_set:
    self.propagate_ray(ri)
    return ri


def test():

    for i in range(1):
        if i == 0:
            return i
            break

t = test()
t = [1,2,3]
t




import chaospy
help(chaospy.WrappedCauchy)
help(chaospy.Normal)

from scipy.stats import skewnorm
import numpy as np
s=np.around(chaospy.Normal(0.5, .025).sample(10000),2)


plt.hist(s, bins=100)
plt.xlim((-5,5))
plt.show()


from scipy.stats import skewnorm, cauchy
import matplotlib.pyplot as plt

import random
#Plot histogram to check skewness

dist=chaospy.Cauchy

ni_rays=10000
mu = .5
linewidth=.05
skew=2

sigma=linewidth/2
args=(mu, sigma)

spectral_vals = np.around(dist(*args).sample(ni_rays), 3)
# cauchy.rvs(loc=mu, scale=linewidth, size=1000)
random.choices(spectral_vals,k=10)

def truncated_cauchy_rvs(loc=0, scale=1, a=0, b=1, skew=.02, size=None):
    """
    Generate random samples from a truncated Cauchy distribution.

    `loc` and `scale` are the location and scale parameters of the distribution.
    `a` and `b` define the interval [a, b] to which the distribution is to be
    limited.

    With the default values of the parameters, the samples are generated
    from the standard Cauchy distribution limited to the interval [-1, 1].
    """
    ua = np.arctan((a - loc)/scale)/np.pi + 0.5
    ub = np.arctan((b - loc)/scale)/np.pi + 0.5
    U = np.random.uniform(ua, ub, size=size)
    rvs =  loc + scale * np.tan(np.pi*(U - 0.5))
    return rvs

x = truncated_cauchy_rvs(loc=0.5, scale=linewidth/2, size=ni_rays)
ignore = plt.hist(x, bins=100, density=True, color='green', alpha=0.25)
y = truncated_cauchy_rvs(loc=0.65, scale=linewidth/2, size=ni_rays)
ignore = plt.hist(y, bins=100, density=True, color='red', alpha=0.25)
plt.grid(alpha=0.5)


random.choice(x)
# del rays
# new_rays = np.linspace(1,100, 1000)
# new_rays
# processes = 4
# ray_number = int(len(rays)/processes)
# ray_number
#
# from threading import Thread
# threads=[]
# for i in range(len(new_rays)):
#     t = Thread(target=doWork, kwargs={'ri':new_rays[i]})
#     threads.append(t)
#     t.start()
# for t in threads:
#     t.join()
# # t = Thread(target=doWork, kwargs={'rays':new_rays[ray_number*(processes-1):]})
# # t.start()
#
#
import pathos

import random
import time
import sys
from multiprocessing import Pool
from threading import Thread
random.seed()
class Dummy():
    def __init__(self, size):
        self.size=size
    def genList(self, N):
        randomList = []

        #initialize random list with values between 0 and 100
        for i in range(N):
            randomList.append(random.randint(0,10))

        return randomList
    #return the sum of all elements in the list
    #This is the same as "return sum(inList)" but in long form for readability and emphasis
    def sumList(self, inList):
        finalSum = 0

        #iterate over all values in the list, and calculate the cummulative sum
        for value in inList:
            finalSum = finalSum + value*self.size
        return finalSum
    def doWork(self, N):
        #create a random list of N integers
        myList = self.genList(N)
        finalSum = self.sumList(myList)
        return finalSum

if __name__ == '__main__':
    if len(sys.argv) == 2 and sys.argv[1].isdigit():
        N = int(sys.argv[1])
        #mark the start time
        startTime = time.time()

        d = Dummy(.5)
        #create a process Pool with 4 processes
        cpus = 8
        # pool = Pool(processes=cpus)
        threads = []
        for i in range(cpus):
            t = Thread(target=d.doWork, args=[int(N/cpus)])
            threads.append(t)
            t.start()
        for t in threads:
            t.join()
        #map doWork to availble Pool processes
        # results = pool.map(d.doWork, [int(N/8),int(N/8),int(N/8),int(N/8),int(N/8), int(N/8), int(N/8), int(N/8)])

        #sum the partial results to get the final result
        # finalSum = sumList(results)

        #mark the end time
        endTime = time.time()
        #calculate the total time it took to complete the work
        workTime =  endTime - startTime

        #print results
        print ("The job took " + str(workTime) + " seconds to complete")
        # print ("The final sum was: " + str(finalSum))
    else:
        exit(-1)
