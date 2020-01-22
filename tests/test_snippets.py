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
