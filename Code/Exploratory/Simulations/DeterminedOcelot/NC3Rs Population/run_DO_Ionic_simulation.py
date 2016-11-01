# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 14:21:37 2016

@author: Oliver Britton
"""

# Test population of models loop
# 2/2/2016 
# Oliver Britton

# A test of the ability to run population of model simulations in neuron through python
import os
import sys
import time
#import pdb

# Simulation 
curDirectory = os.getcwd()
os.chdir('E:\\CLPC48\\Neuron Project\\Code\\Models\\Currents\\Prototypes')
from neuron import h
h('load_file("nrngui.hoc")')
os.chdir(curDirectory)

sys.path.append('E:\CLPC48\Neuron Project\Code\Methods')
import PopulationOfModels as pom


""" Setup """
# Models
import MakeDeterminedOcelot

if __name__ == '__main__':
    configFilename = "E:\\CLPC48\\Neuron Project\\Simulations\\Input\\cfg\\DO\\DO_NC3Rs.cfg"
    pattern = ": "
    #cfg = pom.ParseConfigFile(configFilename,pattern)
    # Initialise hoc
    #h('load_file("nrngui.hoc")')
    
    # --- Main loop setup --- 

    start = time.time()
    pom.RunPopulationOfModels(configFilename,pattern)
    
    # --- Cleanup ---
    # Write log file - list of parameters, stimulus protocols, date of simulation, model name, and list of all input files and all output files
    """ TODO!
    WriteLogFile()
    """
    end = time.time()
    print "Time taken is: %f seconds." % (end-start)


