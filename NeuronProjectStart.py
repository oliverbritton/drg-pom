# -*- coding: utf-8 -*-
"""
Load top level directory and any other setup information 

Created on Tue Mar 08 10:19:01 2016

@author: comra
"""

import os
import sys


def GetProjectDir():
    f = open('project.hostconfig','r')
    pattern = ': '
    # Can read in more lines here as needed
    line = f.read()
    f.close()
    components = line.split(pattern)
    assert(len(components) == 2)
    projectDir = components[1]
    return projectDir

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
    return os.path.join('Code','Currents','Prototypes')

""" Run Setup """
projectDir = GetProjectDir()
SetPaths(projectDir)

print "All systems nominal."