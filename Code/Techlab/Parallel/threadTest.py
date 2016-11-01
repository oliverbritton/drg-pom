# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 17:38:27 2016

@author: Oliver Britton
"""

import thread
import time
import random

def runOften(threadName, sleepTime):
    i=0    
    while i < 10:
        i += 1
        time.sleep(sleepTime)
        print "%s\n" % (threadName)
        
#def runRandomly(threadName, sleepTime):
#    while 1 < 2:
#        time.sleep(sleepTime)

try:
    thread.start_new_thread(runOften,("Often runs",2))
    thread.start_new_thread(runOften,("Less often runs",5))
    thread.start_new_thread(runOften,("Fast and random",random.random()))
except Exception, e:
    print str(e)
    
    
