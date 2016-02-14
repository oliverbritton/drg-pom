# -*- coding: utf-8 -*-
"""
Created on Mon Feb 08 12:14:25 2016
Functions for running a population of models loop

@author: Oliver Britton
"""
import os
import datetime
from matplotlib import pyplot
curDirectory = os.getcwd()
os.chdir('E:\\CLPC48\\Neuron Project\\Code\\Models\\Currents\\Prototypes')
from neuron import h
os.chdir(curDirectory)

import MakeDeterminedOcelot

# FUNCTIONS

def ReadTextFile(filename):
    f = open(filename,"r") #opens file with name of "test.txt"
    list = []
    for line in f:   
        list.append(line.strip())
        
    f.close()
    return list

def ReadParameterFile(filename):
    
    listOfLists = []
    with open(filename) as f:
        for line in f:
            #inner_list = [elt.strip() for elt in line.split(',')]
            # in alternative, if you need to use the file content as numbers
            innerList = [float(elt.strip()) for elt in line.split(' ')]
            listOfLists.append(innerList)
                 
        f.close()
        return listOfLists
        
def ReadTraceFile(filename):
    lines = ReadTextFile(filename)
    t = []    
    v = []
    for line in lines:
        values = line.split(' ')
        t.append(float(values[0]))
        v.append(float(values[1]))
    # Return a dictionary to allow us to add in other currents and state variables with names
    # later on
    return{'t':t, 'v':v}        
        
# Parse a line of a configuration file where a pattern is found, strip out the identifier through pattern and return the data
def ParseConfigLine(line,pattern):
    components = line.split(pattern)
    assert(len(components) == 2)
    return components[1]

def WriteLogFile(outputDirectory,parameters,protocols,modelName,simulationName,configFilename,parameterSetFilename,outputFilenames):

    # Open log file
    # First check output directory ends with file separator
    outputDirectory = FormatOutputDirectory(outputDirectory)    

    f = file.open("%s %s.log" % (outputDirectory,simulationName) ,'w')

    # Write model name, stimulation protocol name and date 

    f.write("Simulation: %s on model: %s, " % (simulationName,modelName) )
    date = datetime.datetime.now() 
    f.write("performed on: %s/%s/%s at %s:%s\n" % (date.day,date.month,date.year,date.hour,date.minute) )

    # Write input file names
    f.write("Config file: %s\n" % configFilename)
    f.write("Parameter file: %s\n" % parameterSetFilename)

    # Write list of output files
    for i in outputFilenames:
        f.write('%s\n' % i)

    # Close
    f.close()
 
def TestWriteLogFile():
	# Write a log file to test directory
	# Read it back in and check it matches line by line
	# Delete log file
    return

# Ramping stimulus protocol	
def RampStimProtocol(amp,t,dur,delay):
		#Protocol description:
		# 100? ms delay
	if t < delay:
		IStim = 0
	elif t > delay+dur:
		IStim = 0
	elif (t >= delay) & (t <= delay+dur):
		IStim = amp*(t-delay)/dur
	return IStim
	
def SquareStimProtocol(amp,t,dur,delay):
	#FIX COLONS
	if t < delay:
		IStim = 0
	elif t > delay+dur:
		IStim = 0
	elif (t >= delay) & (t <= delay+dur):
		IStim = amp
	return IStim
	

	
def TestStimProtocols():
	amp = 2
	delay = 100
	dur = 10
	
	# Test SquareStimProtocol
	# before stimulation
	t = 0
	assert(SquareStimProtocol(amp,t,dur,delay) == 0)
	# during stim
	t = 105
	assert(SquareStimProtocol(amp,t,dur,delay) == 1)
	t = 110
	assert(SquareStimProtocol(amp,t,dur,delay) == 2)
	# after stim`
	t = 111
	assert(SquareStimProtocol(amp,t,dur,delay) == 0)
	# Test SquareStimProtocol
	t = 0
	assert(SquareStimProtocol(amp,t,dur,delay) == 0)
	# during stim
	t = 105
	assert(SquareStimProtocol(amp,t,dur,delay) == 2)
	t = 110
	assert(SquareStimProtocol(amp,t,dur,delay) == 2)
	# after stim`
	t = 111
	assert(SquareStimProtocol(amp,t,dur,delay) == 0)
 
def GetModel(modelName):
    
    if modelName == 'DeterminedOcelot':
        model = MakeDeterminedOcelot.MakeDeterminedOcelot()
    else:
        assert False, 'Model name not found in GetModel!'    
    
    return model

def SetModelParameters(model,parameters,modelName):

    if modelName == 'DeterminedOcelot':
        
        numParameters = 6
        assert len(parameters) == numParameters, 'number of parameters is wrong' 
        # Sodium conductances
        model.gnabar_nav17vw = model.gnabar_nav17vw*parameters[0]
        model.gnabar_nav18hw = model.gnabar_nav18hw*parameters[1]
        model.gnabar_nav19hw = model.gnabar_nav19hw*parameters[2]
        
        # Potassium conductances
        model.gkbar_kdrtf = model.gkbar_kdrtf*parameters[3]
        model.gkbar_katf = model.gkbar_katf*parameters[4]
        model.gkbar_kmtf = model.gkbar_kmtf*parameters[5]
    
    else:     
        assert False, 'Model name not found in SetModelParameters!'
    
    return

    
	
    # Run simulation protocol
def RunSimulation(model,parameters,modelName,protocol,outputDirectory,prefix,modelNum):
    
    curDirectory = os.getcwd()
    os.chdir('E:\\CLPC48\\Neuron Project\\Code\\Models\\Currents\\Prototypes')
    from neuron import h
    import neuron
    h('load_file("nrngui.hoc")')
    os.chdir(curDirectory)
    
    # For each simulation:
    # Reinitialise model
    model = GetModel(modelName)
    
    if not protocol == 'default':
        assert False, 'Unsupported protocol'

    SetModelParameters(model,parameters,modelName)    
    
    # Setup output vectors
    v = h.Vector()
    t = h.Vector()
    ina_vec = h.Vector()
    icurr_vec = h.Vector()

    # Setup simulation parameters
    stim = h.IClamp(model(0.5))
    stim.delay = 100
    stim.dur = 700
    stim.amp = 0.5 # nA (1 nA = 100 pA)

    v.record(model(0.5)._ref_v, sec=model)
    t.record(h._ref_t)
    #ina_vec.record(cell(0.5)._ref_ina)
    #icurr_vec.record(cell(0.5)._ref_ina_nav18hw, sec=model)  
    h.finitialize(-65) # Vital! And has to go after record 
    tstop = 1000.0
    
    # Run simulation
    neuron.run(tstop)
    
    # Write output
    WriteSimulationOutput(outputDirectory,prefix,modelNum,t,v)


def WriteSimulationOutput(outputDirectory,prefix,modelNum,t,v):

    assert len(v) == len(t), 't and v vector length mismatch'      
    
    # Check output directory is in the correct form
    outputDirectory = FormatOutputDirectory(outputDirectory)
    
    # Open output
    filename = outputDirectory + prefix + str(modelNum) + '.dat'
    f = open(filename, "w")

    # Write AP, two columns, end with line break 

    for i in range(len(t)):
        f.write(str(t[i]) + " " + str(v[i]) + "\n")

    f.close()

def StepProtocolRun():
    
        # Run the step protocol
    for j, amp in enumerate(stimAmps):
        
        stim = h.IClamp(cell(0.5)) 
        stim.delay = 5
        stim.dur = stimDuration
        stim.amp = amp 
        
        v_vec = h.Vector()
        t_vec = h.Vector()
            
        
        v_vec.record(cell(0.5)._ref_v, sec=cell)
        t_vec.record(h._ref_t)
        
        h.finitialize(-65) # Vital! And has to go after record   
        tstop = 1000.0
        neuron.run(tstop)
        
def FormatOutputDirectory(outputDirectory):
    if not outputDirectory.endswith(os.sep):
        outputDirectory += os.sep
    return outputDirectory
    
def PlotTrace(filename):
    
    return
    
# Load Neuron so that all the density mechanisms get imported properly
def LoadNeuron():
    curDirectory = os.getcwd()
    os.chdir('E:\\CLPC48\\Neuron Project\\Code\\Models\\Currents\\Prototypes')
    from neuron import h
    h('load_file("nrngui.hoc")')
    os.chdir(curDirectory)