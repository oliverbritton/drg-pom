# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 14:04:36 2016

@author: Oliver Britton
"""

from multiprocessing import Process
from multiprocessing import Pool
import time

def f(name):
    print 'hello', name

def g(x):
    start = time.time()
    for i in range(x):
        i**i
    end = time.time()
   # print "Time to run is %.1f seconds" % (end-start)
    
def TestProcesses(x,y):	
    start = time.time()
    count = 0
    while count < 4:
    	jobs = []
	print "Round we go!\n"
    	for i in range(x):
            p = Process(target=g, args=(y,))
            p.start()        
            jobs.append(p)
        
    	for job in jobs:
        	job.join()
		count +=1
    end = time.time()	
    print "Time taken with %i cores is %f seconds.\n" % (x, (end-start))    


if __name__ == '__main__':
    y = 10000
    for i in range(1,5):
	print i
	TestProcesses(i,y)

#    start = time.time()
#    p1 = Process(target=g, args=(y,))
#    p2 = Process(target=g, args=(y,))
#    p3 = Process(target=g, args=(y,))
#    p4 = Process(target=g, args=(y,))
#    p5 = Process(target=g, args=(y,))
#    p1.start()
#    p2.start()
#    p3.start()
#    p4.start()
#    p5.start() 
#    p1.join()
#    p2.join()
#    p3.join()
#    p4.join()
#    p5.start()
#    p5.join()
#    end = time.time()
#    print "Total time taken is %f seconds.\n" % (end-start)
    
    

        
