rays=[1,2,3,4,5,6,7,8,9,10,11,12,13,14]

new_rays = [rays.pop(0) for i in range(len(rays))]
rays
import numpy as np
new_rays = np.linspace(1,100, 1000)
new_rays
processes = 4
ray_number = int(len(rays)/processes)
ray_number

from threading import Thread
import time
def doWork(ri):

    time.sleep(1)
    print(ri)
#mark the start time
startTime = time.time()
threads=[]
for i in range(len(new_rays)):
    t = Thread(target=doWork, kwargs={'ri':new_rays[i]})
    threads.append(t)
    t.start()
for t in threads:
    t.join()
# t = Thread(target=doWork, kwargs={'rays':new_rays[ray_number*(processes-1):]})
# t.start()

#mark the end time
endTime = time.time()
#calculate the total time it took to complete the work
workTime =  endTime - startTime
#print results
print ("The job took " + str(workTime) + " seconds to complete")



import random
import time
import sys
from multiprocessing import Pool
random.seed()
def genList (size):
    randomList = []

    #initialize random list with values between 0 and 100
    for i in range(size):
        randomList.append(random.randint(0,10))

    return randomList
#return the sum of all elements in the list
#This is the same as "return sum(inList)" but in long form for readability and emphasis
def sumList(inList):
    finalSum = 0

    #iterate over all values in the list, and calculate the cummulative sum
    for value in inList:
        finalSum = finalSum + value
    return finalSum
def doWork(N):
    #create a random list of N integers
    myList = genList (N)
    finalSum = sumList(myList)
    return finalSum
if __name__ == '__main__':
    if len(sys.argv) == 2 and sys.argv[1].isdigit():
        N = int(sys.argv[1])
        #mark the start time
        startTime = time.time()

        #create a process Pool with 4 processes
        cpus = 4
        pool = Pool(processes=cpus)


        #map doWork to availble Pool processes
        results = pool.map(doWork, [i for i in range(N)])

        #sum the partial results to get the final result
        finalSum = sumList(results)

        #mark the end time
        endTime = time.time()
        #calculate the total time it took to complete the work
        workTime =  endTime - startTime

        #print results
        print ("The job took " + str(workTime) + " seconds to complete")
        print ("The final sum was: " + str(finalSum))
    else:
        exit(-1)
