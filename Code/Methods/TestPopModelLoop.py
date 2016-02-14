# Test population of models loop
# 2/2/2016 
# Oliver Britton

# A test of the ability to run population of model simulations in neuron through python
import os
import sys
import datetime
import thread
import re
from matplotlib import pyplot
import time

# Simulation 
curDirectory = os.getcwd()
os.chdir('E:\\CLPC48\\Neuron Project\\Code\\Models\\Currents\\Prototypes')
from neuron import h
h('load_file("nrngui.hoc")')
os.chdir(curDirectory)

sys.path.append('E:\CLPC48\Neuron Project\Code\Methods')
import PopulationOfModels as pom

# Models
import MakeDeterminedOcelot

# Initialise hoc
h('load_file("nrngui.hoc")')

# --- Loop setup ---

# Open config file
configFilename = "E:\\CLPC48\\Neuron Project\\Simulations\\Input\\cfg\\secondTest.cfg"
configFile = pom.ReadTextFile(configFilename)


# Define output directory from config file
# Format expected is: output directory: blah/blah

pattern = ": " # pattern to split lines between data for program and identifier
numOutputDirectories = 0
numParameterFiles = 0
numModels = 0
numSimulationNames = 0
numPrefixes = 0
numProtocols = 0

#% REFACTOR - turn into function if we use it again
for line in configFile:
    if re.search('output',line):
        numOutputDirectories += 1
        outputDirectory = pom.ParseConfigLine(line,pattern)
        
    if re.search('parameter',line):
        numParameterFiles += 1
        parameterFilename = pom.ParseConfigLine(line,pattern)     
        
    if re.search('model',line):
        numModels += 1
        modelName = pom.ParseConfigLine(line,pattern)
    
    if re.search('simulation',line):
        numSimulationNames += 1
        simulationName = pom.ParseConfigLine(line,pattern)
        
    if re.search('prefix',line):
        numPrefixes += 1
        prefix = pom.ParseConfigLine(line,pattern)

    if re.search('protocol',line):
        numProtocols +=1
        protocol = pom.ParseConfigLine(line,pattern)
        
assert numOutputDirectories == 1, 'numOutputDirectories'
assert numParameterFiles == 1, 'numParameterFiles'
assert numModels == 1, 'numModels'
assert numSimulationNames == 1, 'numSimulationNames'
assert numPrefixes == 1, 'numPrefixes'
assert numProtocols == 1, 'numProtocols'

    
#%% Read parameter file in
parameters = pom.ReadParameterFile(parameterFilename)
start = time.time()

parameters1 = parameters[0:49]
parameters2 = parameters[50:99]

# --- Main solver loop --- # TO DO! Parallelise
for modelNum,parameterSet in enumerate(parameters):

    # Initialise new model
    model = pom.GetModel(modelName)
    
    # Update parameters
    
    # Run simulation protocol
    pom.RunSimulation(model,parameterSet,modelName,protocol,outputDirectory,prefix,modelNum)
    
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


    thread.start_new_thread(runOften,("Often runs",2))

    
    


