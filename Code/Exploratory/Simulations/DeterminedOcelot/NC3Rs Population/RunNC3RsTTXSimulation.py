# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 15:04:20 2016

@author: Oliver Britton
"""


# TTX run of NC3Rs population
# 15/12/2016
# Oliver Britton


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
    #configFilename = "E:\\CLPC48\\Neuron Project\\Simulations\\Input\\cfg\\DO\\DO_NC3Rs_TTX.cfg" # TTX 100% Nav 1.7 block
    configFilename = "E:\\CLPC48\\Neuron Project\\Simulations\\Input\\cfg\\DO\\DO_NC3Rs_TTX_10.cfg"
    pattern = ": "
    #cfg = pom.ParseConfigFile(configFilename,pattern)
    # Initialise hoc
    #h('load_file("nrngui.hoc")')
    
    # --- Main loop setup --- 

    start = time.time()
    #pom.RunPopulationOfModels(configFilename, pattern, parameter_modifiers={0:0.0}) # TTX - model as 100% decrease in Nav 1.7
    pom.RunPopulationOfModels(configFilename, pattern, parameter_modifiers={0:0.1}) # TTX - model as 90% decrease in Nav 1.7
    
    # --- Cleanup ---
    # Write log file - list of parameters, stimulus protocols, date of simulation, model name, and list of all input files and all output files
    """ TODO!
    WriteLogFile()
    """
    end = time.time()
    print "Time taken is: %f seconds." % (end-start)


