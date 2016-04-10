# Test population of models loop
# 2/2/2016 
# Oliver Britton

# A test of the ability to run population of model simulations in neuron through python
import os
import sys
#import datetime
#import re
#from matplotlib import pyplot
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
    configFilename = "E:\\CLPC48\\Neuron Project\\Simulations\\Input\\cfg\\secondTest.cfg"
    pattern = ": "
    #cfg = pom.ParseConfigFile(configFilename,pattern)
    # Initialise hoc
    #h('load_file("nrngui.hoc")')
    
    # --- Loop setup --- 
    
    # Open config file
    start = time.time()
    pom.RunPopulationOfModels(configFilename,pattern)
    
        
    #%% Read parameter file in
    #parameters = pom.ReadParameterFile(cfg['parameterFilename'])
    
    # --- Main solver loop --- # TO DO! Parallelise
    #for modelNum,parameterSet in enumerate(parameters):
    
        # Initialise new model
    #    model = pom.GetModel(cfg['modelName'])
        
        # Update parameters
        
        # Run simulation protocol
    #    pom.RunSimulation(model,parameterSet,cfg['modelName'],cfg['protocol'],cfg['outputDirectory'],cfg['prefix'],modelNum)
        
        # For each simulation:
        # Reinitialise model
        # Setup output vectors
        # Setup simulation parameters
        # Run simulation
        # Open output file
        # Write to output and close
    
    # --- Cleanup ---
    # Write log file - list of parameters, stimulus protocols, date of simulation, model name, and list of all input files and all output files
    """ TODO!
    WriteLogFile()
    """
    end = time.time()
    print "Time taken is: %f seconds." % (end-start)


