# -*- coding: utf-8 -*-
"""
Created on Mon Feb 08 12:14:25 2016
Functions for running a population of models loop

@author: Oliver Britton
"""
import os
import datetime
import re
#import pdb
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#from multiprocessing import Pool
from multiprocessing import Process
from functools import partial
import NeuronProjectStart
import Methods.Biomarkers.NeuronBiomarkers as nb
import Methods.Biomarkers.DavidsonBiomarkers as db
# Load all channel models and return to the right directory
" TODO - can I refactor this to use load neuron library?"
curDirectory = os.getcwd()
projectDir = NeuronProjectStart.GetProjectDir()
nrnChannelDir = NeuronProjectStart.GetNrnChannelDir()
os.chdir(os.path.join(projectDir,nrnChannelDir))
from neuron import h
os.chdir(curDirectory)



# FUNCTIONS

def ReadTextFile(filename):

    with open(filename,'rt') as f:
        text = f.readlines()
    return text

" To do - replace with pandas "
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
        
" To do - replace with acs.readtrace "
def ReadTraceFile(filename,lineskip=0):
    filename = filename.rstrip('\n') # Chomp off new line chars
    lines = ReadTextFile(filename)
    t = []    
    v = []
    for i, line in enumerate(lines):
        if i >= lineskip:
          t.append(float(line.split()[0]))
          v.append(float(line.split()[1]))
    # Return a dictionary to allow us to add in other currents and state variables with names
    # later on
    return{'t':t, 'v':v}     
    
def read_trace(filename, skiprows=0):
    # Refactored to use numpy from pom.ReadTraceFile,
    data = np.loadtxt(filename,skiprows=skiprows)
    trace = {'t':data[:,0], 'v':data[:,1]}
    return trace
    
def PlotTraceFile(filename):
    data = ReadTraceFile(filename)
    plt.plot(data['t'],data['v'])
    return
        
# Parse a line of a configuration file where a pattern is found, strip out the identifier through pattern and return the data
def ParseConfigLine(line,pattern):
    line = line.rstrip('\n') # Remove new line
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
 
def GetModel(model_name):
    
    if model_name == 'DeterminedOcelot':
        import MakeDeterminedOcelot
        model = MakeDeterminedOcelot.MakeDeterminedOcelot()
    elif model_name == 'DeterminedOcelotIonic':
        import MakeDeterminedOcelot
        model = MakeDeterminedOcelot.MakeDeterminedOcelotIonic()
    else:
        assert False, 'Model name not found in GetModel!'    
    
    return model

def SetModelParameters(model,parameters,model_name):

    if (model_name == 'DeterminedOcelot') or (model_name == 'DeterminedOcelotIonic'):
        
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
def RunSimulation(model, parameters, modelName, protocol, outputDirectory, prefix, modelNum, calibrationBiomarkerFile, allBiomarkerFile):
    # TODO biomarkers should be returned and output of biomarkers done somewhere else
    curDirectory = os.getcwd()
    projectDir = NeuronProjectStart.GetProjectDir()
    nrnChannelDir = NeuronProjectStart.GetNrnChannelDir()
    os.chdir(os.path.join(projectDir,nrnChannelDir))
    from neuron import h
    import neuron
    h('load_file("nrngui.hoc")')
    os.chdir(curDirectory)
    
    # For each simulation:
    # Reinitialise model
    model = GetModel(modelName)
    SetModelParameters(model,parameters,modelName)    
    
    protocolData = GetSimulationProtocol(protocol)
    numSimulations = len(protocolData.index)
    

    for simulation in range(numSimulations):
    
        # Setup output vectors
        v = h.Vector()
        t = h.Vector()
        #ina_vec = h.Vector()
        #icurr_vec = h.Vector()
    
        # Setup simulation parameters - MAKE INTO FUNCTION
        # SetStimulus()
        stim = h.IClamp(model(0.5))
        stim.delay = protocolData.loc[simulation]['stim start']
        stim.dur = protocolData.loc[simulation]['stim duration']
        stim.amp = protocolData.loc[simulation]['stim amp']/1000.0 # nA (1 nA = 1000 pA) - protocols are given in pA
    
        v.record(model(0.5)._ref_v, sec=model)
        t.record(h._ref_t)
        #ina_vec.record(cell(0.5)._ref_ina)
        #icurr_vec.record(cell(0.5)._ref_ina_nav18hw, sec=model)  
        h.finitialize(-65) # Vital! And has to go after record 
        tstop = protocolData.loc[simulation]['duration']
        
        # Run simulation
        neuron.run(tstop)
        
        # Calculate biomarkers
        t_np = np.array(t) # Transform into numpy arrays to make behaviour more understandable
        v_np = np.array(v)
        biomarkers = CalculateBiomarkers(t_np,v_np,stim.amp,modelNum)  # modelNum is the index variable
        if stim.amp == 0:
            quiescentRMP = biomarkers['RMP']
        # TODO sort out whether stim amp and index should be called here
        # organise better -  maybe a master function that calls calculate biomarkers and also gets stim amp and index?

        # Write biomarkers
        nb.WriteBiomarkers(biomarkers,allBiomarkerFile)        
        
        if biomarkers['numAPs'] > 0:
            rheobaseOnThisBeat = True
        else:
            rheobaseOnThisBeat = False
        
        # If this is the last simulation or we found rheobase, then write the trace     
        # TODO - write code to allow us to check for rheobase, but still get biomarkers for increasing
        # stim amplitude to see what happens above threshold
        if (rheobaseOnThisBeat):
            # Write output
            biomarkers['StepRheobase'] = stim.amp
            biomarkers['RMP'] = quiescentRMP # Write the RMP with the quiescent value
            nb.WriteBiomarkers(biomarkers,calibrationBiomarkerFile)
            WriteSimulationOutput(outputDirectory,prefix,modelNum,t,v)
            break
        
        if (simulation == numSimulations-1): # -1 because range gives 0 -> n-1 if numSimulations = n
        # We haven't found rheobase
            nb.WriteBiomarkers(biomarkers,calibrationBiomarkerFile)
            WriteSimulationOutput(outputDirectory,prefix,modelNum,t,v)            
            
    
    return    
    
# Run parallel simulation
def RunParallelSimulation(modelNumsAndParameters,modelName,protocol,outputDirectory,prefix):
    
#    curDirectory = os.getcwd()
#    os.chdir('E:\\CLPC48\\Neuron Project\\Code\\Models\\Currents\\Prototypes')
#    from neuron import h
    import neuron
#    h('load_file("nrngui.hoc")')
#    os.chdir(curDirectory)  
    
    # Create model
    model = GetModel(modelName)
    
    if not protocol == 'default':
        assert False, 'Unsupported protocol'
    
    # Unpack model number and parameteres    
    modelNum = int(modelNumsAndParameters[0])
    parameters = modelNumsAndParameters[1:]
    
    SetModelParameters(model,parameters,modelName)    
    
    # Setup output vectors
    v = h.Vector()
    t = h.Vector()
    #ina_vec = h.Vector()
    #icurr_vec = h.Vector()
    
    # Setup simulation parameters - MAKE INTO FUNCTION
    # SetStimulus()
    stim = h.IClamp(model(0.5))
    stim.delay = 100
    stim.dur = 700
    stim.amp = 0.5 # nA (1 nA = 1000 pA)
    
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
    return

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
    return

def StepProtocolRun(model, stimAmps, stimDuration):
    
        # Run the step protocol
    for j, amp in enumerate(stimAmps):
        
        stim = h.IClamp(model(0.5)) 
        stim.delay = 5
        stim.dur = stimDuration
        stim.amp = amp 
        
        v_vec = h.Vector()
        t_vec = h.Vector()
            
        
        v_vec.record(model(0.5)._ref_v, sec=model)
        t_vec.record(h._ref_t)
        
        h.finitialize(-65) # Vital! And has to go after record   
        tstop = 1000.0
        neuron.run(tstop)
        # TODO Write to file
        
def FormatOutputDirectory(outputDirectory):
    if not outputDirectory.endswith(os.sep):
        outputDirectory += os.sep
    return outputDirectory
    
def PlotTrace(filename):
    
    return
    
# Load Neuron so that all the density mechanisms get imported properly
" Deprecated - can use the loadneuron module to load density mechanisms now."
def LoadNeuron():
    curDirectory = os.getcwd()
    projectDir = NeuronProjectStart.GetProjectDir()
    nrnChannelDir = NeuronProjectStart.GetNrnChannelDir()
    os.chdir(os.path.join(projectDir,nrnChannelDir))
    from neuron import h
    h('load_file("nrngui.hoc")')
    os.chdir(curDirectory)
    return

def ParseConfigFile(configFilename,pattern):
    
    configFile = ReadTextFile(configFilename)
    # Define output directory from config file
    # Format expected is: output directory: blah/blah
    
    #pattern = ": " # pattern to split lines between data for program and identifier
    numOutputDirectories = 0
    numParameterFiles = 0
    numModels = 0
    numSimulationNames = 0
    numPrefixes = 0
    numProtocols = 0

    " To do - could refactor this to make the dict first, then fill it in with each line. "
    " This would be neater as we could have a list of tokens to look for, and then loop over them "
    " rather than having a separate if statement for each one. "    
    for line in configFile:
        if re.search('output',line):
            numOutputDirectories += 1
            outputDirectory = ParseConfigLine(line,pattern)
            
        if re.search('parameter',line):
            numParameterFiles += 1
            parameterFilename = ParseConfigLine(line,pattern)     
            
        if re.search('model',line):
            numModels += 1
            modelName = ParseConfigLine(line,pattern)
        
        if re.search('simulation',line):
            numSimulationNames += 1
            simulationName = ParseConfigLine(line,pattern)
            
        if re.search('prefix',line):
            numPrefixes += 1
            prefix = ParseConfigLine(line,pattern)
    
        if re.search('protocol',line):
            numProtocols +=1
            protocol = ParseConfigLine(line,pattern)
            
    assert numOutputDirectories == 1, 'numOutputDirectories'
    assert numParameterFiles == 1, 'numParameterFiles'
    assert numModels == 1, 'numModels'
    assert numSimulationNames == 1, 'numSimulationNames'
    assert numPrefixes == 1, 'numPrefixes'
    assert numProtocols == 1, 'numProtocols'
    
    return {'outputDirectory':outputDirectory, 'parameterFilename':parameterFilename,
            'modelName':modelName, 'simulationName':simulationName, 'prefix':prefix,
            'protocol':protocol}
            
def RunPopulationOfModels(configFilename, pattern, parameter_modifiers={}):
    
    cfg = ParseConfigFile(configFilename,pattern)
    # Initialise hoc
   # h('load_file("nrngui.hoc")')
    """ Can I just use LoadNeuron here? TODO """
    curDirectory = os.getcwd()    
    projectDir = NeuronProjectStart.GetProjectDir()
    nrnChannelDir = NeuronProjectStart.GetNrnChannelDir()
    os.chdir(os.path.join(projectDir,nrnChannelDir))    
    from neuron import h
    import neuron
    h('load_file("nrngui.hoc")')
    os.chdir(curDirectory)    
    
    # --- Loop  --- 
    #Read parameter file in
    parameters = ReadParameterFile(cfg['parameterFilename'])
    start = time.time()
    # --- Main solver loop --- 
    
    # Set up output
    outputDirectory = cfg['outputDirectory']
    prefix = cfg['prefix']
    
    " !!!TO CHANGE - replace with dataframes "
    calibrationBiomrakerFile = open(outputDirectory + prefix + '_biomarkers.dat','w') # File for rheobase biomarkers only for calibration
    allBiomarkerFile = open(outputDirectory + prefix + '_allbiomarkers.dat','w') # File for every biomarker
    nb.WriteHeader(calibrationBiomarkerFile)
    nb.WriteHeader(allBiomarkerFile)
    " !!!To CHANGE "

    for modelNum,parameterSet in enumerate(parameters):    
        # Initialise new model
        model = GetModel(cfg['modelName'])       
        # Modify parameters (to do, include this info in CFG)       
        if parameter_modifiers:
            assert type(parameter_modifiers) == dict
            for param_index, param_val in parameter_modifiers.iteritems():
                print("Old parameter set was {}".format(parameterSet))
                parameterSet[param_index] = param_val
                print("Modified parameter set is {}".format(parameterSet))
        # Run simulation protocol
        RunSimulation(model,parameterSet,cfg['modelName'],cfg['protocol'],cfg['outputDirectory'],cfg['prefix'],modelNum,calibrationBiomarkerFile,allBiomarkerFile)
        if (modelNum % 100) == 0:
            print("Completed model: {}".format(modelNum))
    
    # Clean up
    calibrationBiomarkerFile.close()
    allBiomarkerFile.close()
    
    end = time.time()
    return (end-start)
            
            
def RunParallelPopulationOfModels(configFilename,pattern,numProcessors):
    
    cfg = ParseConfigFile(configFilename,pattern)
    # Initialise hoc
   # h('load_file("nrngui.hoc")')

    curDirectory = os.getcwd()
    """ Todo loadNeuron"""
    projectDir = NeuronProjectStart.GetProjectDir()
    nrnChannelDir = NeuronProjectStart.GetNrnChannelDir()
    os.chdir(os.path.join(projectDir,nrnChannelDir))
    from neuron import h
    import neuron
    h('load_file("nrngui.hoc")')
    os.chdir(curDirectory)    
    
    # --- Loop  --- 
    #Read parameter file in
    parameters = ReadParameterFile(cfg['parameterFilename'])
#    start = time.time()
    # --- Main solver loop --- # TO DO! Parallelise
    
    # Divide the number of models up between processors    
    
# Might be better to use a pool    
# start numProcessors worth of processes

    modelNums = range(len(parameters))
    modelNums = np.reshape(modelNums,[len(parameters),1])    
    ###pdb.set_trace()
    # Make the iterable matrix
    parametersAndModelNums = np.concatenate((modelNums,parameters),1)
    
    RunPartialParallelSimulation = partial(RunParallelSimulation,modelName = cfg['modelName'],protocol = cfg['protocol'],outputDirectory = cfg['outputDirectory'], prefix = cfg['prefix'])   
#    pdb.set_trace()
    
    for i in parametersAndModelNums:
        
#       RunPartialParallelSimulation(i)
       p1 = Process(target=RunPartialParallelSimulation, args=(i,))
       i[0] += 100
       p2 = Process(target=RunPartialParallelSimulation, args=(i,))
       i[0] += 100
       p3 = Process(target=RunPartialParallelSimulation, args=(i,))
       i[0] += 100
       p4 = Process(target=RunPartialParallelSimulation, args=(i,))
       p1.start()
       p2.start()
       p3.start()
       p4.start()
       p1.join()
       p2.join()
       p3.join()
       p4.join()
    # Set up a pool and run simulations
#    pool =  Pool(numProcessors)
#    pool.map(RunPartialParallelSimulation, parametersAndModelNums)
#    pool.close() 
#    pool.join()

 #p5 = Process(target=g, args=(y,))
    
#    for modelNum,parameterSet in enumerate(parameters):    
#        # Initialise new model
#        model = GetModel(cfg['modelName'])       
#        # Run simulation protocol
#        RunSimulation(model,parameterSet,cfg['modelName'],cfg['protocol'],cfg['outputDirectory'],cfg['prefix'],modelNum)
#        
#    end = time.time()
    return
    
def GenerateSimulationProtocol():   
    
    """ Basically we have a list of protocol functions that return the simulations to run
    this function interfaces with those functions, or whatever we want to use in the future
    and sends a dictionary back to the main run simulation function with the details it needs to 
    run all the simulations """

    # Setup simulation parameters - MAKE INTO FUNCTION
    # SetStimulus()
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

def GetSimulationProtocol(protocol):
    # Protocol list defines the available protocols
    protocolList = {}
    protocolList['default'] = 'Basic step protocol in line with Davidson et al. 2014, PAIN'
    assert protocolList.has_key(protocol), "Stimulus protocol not found"
    
    if protocol == 'default':
        # Create a data object containing all the stimulus properties
        # Simulations: 
        # 1: D 1000 ms Stim 1 start = 0 amp = 0 duration = 0
        # 2: D 1000 ms Stim 1 Start 100 Amp 50 Duration 800 
        # 3: Same as 2 but amp+50
        
#        amps = range(0,3501,50)
        amps = range(0,51,10)        
        numSimulations = len(amps)
        durations = [1000]*numSimulations
        starts = [100]*numSimulations
        starts[0] = 0
        stimDurations = [800]*numSimulations
        stimDurations[0] = 0
        
        data = {'duration':durations, 'stim start':starts, 'stim amp':amps, 'stim duration': stimDurations}
        protocolData = pd.DataFrame(data,columns = ['duration', 'stim start', 'stim amp', 'stim duration'])
        
    return protocolData
    
def CalculateBiomarkers(t,v,stimAmp,index):
    # Calculate a standard set of biomarkers and return them
    threshold = 0 # mV
    timeThreshold = 5 # ms
    dvdtThreshold = 5 # mV/ms
    traces = nb.SplitTraceIntoAPs(t,v,threshold,timeThreshold) 
    numAPs = traces['numAPs']
    
    biomarkers = {}
    for name in db.biomarkerNames:
        biomarkers[name] = 'N/A'        
        
    biomarkers['Index'] = index
    biomarkers['numAPs'] = numAPs
    biomarkers['stimAmp'] = stimAmp  
    
    RMP = nb.CalculateRMP(traces)
    biomarkers['RMP'] = RMP
    if numAPs > 0:
    
        # Need to calculate biomarkers for each AP in the trace
        APPeak = nb.CalculateAPPeak(traces)
        APRise = nb.CalculateAPRiseTime(traces,dvdtThreshold)
        APSlopeMin, APSlopeMax = nb.CalculateAPSlopeMinMax(traces)
        APWidth = nb.CalculateAPFullWidth(traces,threshold)
        AHPAmp = nb.CalculateAHPAmp(traces,dvdtThreshold)
        
        biomarkers['APPeak'] = APPeak
        biomarkers['APRise'] = APRise
        biomarkers['APSlopeMin'] = APSlopeMin
        biomarkers['APSlopeMax'] = APSlopeMax
        biomarkers['APWidth'] = APWidth
        biomarkers['AHPAmp'] = AHPAmp
        

        
        #TestRheobase([trace,trace],[50,100])
 
    return biomarkers
    
def CalibrateBiomarkers(modelBiomarkers,calibrationRanges,biomarkerNames):
    
    # iterate over models
#    for biomarker in 
    
    # iterate over each biomarker to create array of 1s and 0s
    
    # at end of each model, if all checked biomarkers are within range, add to list of indices
    
    # return list of calibrated indices
    return 
