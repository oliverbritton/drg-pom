# -*- coding: utf-8 -*-
"""
Load top level directory and any other setup information 

Created on Tue Mar 08 10:19:01 2016
"""

import os
import sys
import time
import platform
import socket

try:
    from winsound import Beep
except ModuleNotFoundError:
    def Beep(a=None,b=None):
        print("Beep: {}, {}.".format(a,b))

def FindProjectDir():
    " Function not in use as we don't use StartScript or hostconfig anymore "
    raise Error('FindProjectDir not currently in use.')
    f = open(os.path.join(start,'C:\\Users\\comra\\Dropbox\\Backups\\Neuron\\project.hostconfig'),'r')
    pattern = ': '
    # Can read in more lines here as needed
    line = f.read()
    f.close()
    components = line.split(pattern)
    assert(len(components) == 2)
    projectDir = components[1]
    return projectDir

def GetProjectDir(): 
    return get_project_path()


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

''' neuron_paths code '''

def get_pc_name():
    n1 = platform.node()
    n2 = socket.gethostname()
    if "COMPUTERNAME" in os.environ.keys():
        n3 = os.environ["COMPUTERNAME"]
    elif "HOSTNAME" in os.environ.keys():
        n3 = os.environ["HOSTNAME"]
    else:
        n3 = None

    if n1 == n2 == n3:
        return n1
    elif n1 == n2:
        return n1
    elif n1 == n3:
        return n1
    elif n2 == n3:
        return n2
    else:
        raise Exception("Computernames are not equal to each other")

def get_project_path():
    pc_name = get_pc_name()
    if pc_name == "clpc403.cs.ox.ac.uk":
        path = "/home/scratch/olibri/Dropbox/Backups/Neuron"  
    elif pc_name == "DESKTOP-DJR4JS3":
        path = "C:\\Users\\comra\\Dropbox\\Backups\\Neuron" 
    elif (pc_name == "OLIVERBRITTON") | (pc_name == "OliverBritton"):
        path = "E:\\CLPC48\\Neuron Project" # Could put Dropbox here instead
    else:
        raise Exception("Name: {} not found".format(pc_name))
    return path

def get_nb_path():
    pc_name = get_pc_name()
    if pc_name == "clpc403.cs.ox.ac.uk":
        path = "/home/scratch/olibri/Dropbox/Backups/Python/DRG"  
    elif pc_name == "DESKTOP-DJR4JS3":
        path = "C:\\Users\\comra\\Dropbox\\Backups\\Python\\DRG" 
    elif (pc_name == "OLIVERBRITTON") | (pc_name == "OliverBritton"):
        path = "C:\\Dropbox\\Backups\\Python\\DRG" 
    else:
        raise Exception("Name: {} not found".format(pc_name))
    return path
    
def setup_paths():
    path = get_project_path()
    paths = [path, os.path.join(path, "Code")]
    for path in paths:
        sys.path.append(path)
    
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

