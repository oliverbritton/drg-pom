# -*- coding: utf-8 -*-
"""
Load top level directory and any other setup information 

Created on Tue Mar 08 10:19:01 2016

@author: comra
"""

import os
import sys
from winsound import Beep
import time

def FindProjectDir():
        
    f = open(os.path.join(start,'E:\\CLPC48\\Neuron Project\\project.hostconfig'),'r')
    pattern = ': '
    # Can read in more lines here as needed
    line = f.read()
    f.close()
    components = line.split(pattern)
    assert(len(components) == 2)
    projectDir = components[1]
    return projectDir

def GetProjectDir():
    return 'E:\\CLPC48\\Neuron Project'

def SetPaths(projectDir):
    codeDir = os.path.join(projectDir,'Code')
    figsDir = os.path.join(projectDir,'Figures') 
    simDir = os.path.join(projectDir,'Simulations')
    dataDir = os.path.join(projectDir,'Data')

    sys.path.append(projectDir)
    sys.path.append(codeDir)    
    sys.path.append(figsDir)     
    sys.path.append(simDir)
    sys.path.append(dataDir)

def GetNrnChannelDir():
    return os.path.join('Code','Models','Currents','Prototypes')

def StartScript():
    projectDir = FindProjectDir()
    SetPaths(projectDir)    
    
def get_my_path():
    import fake
    path = str(fake).split()[3][1:-9]
    os.remove( os.path.join( path, 'fake.pyc' ) )
    return path

def rdy():
    for i in range(3,8):
        Beep(i*100,450)
        time.sleep(0.01)
    
""" Run Setup """
projectDir = GetProjectDir()
SetPaths(projectDir)
# __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
#print __location__
#print __file__

import inspect
#print inspect.getfile(inspect.currentframe()) # script filename (usually with path)
#print os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script di
#
#print "All systems nominal."