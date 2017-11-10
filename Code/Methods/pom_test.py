# -*- coding: utf-8 -*-
"""
Created on Mon Feb 08 12:14:25 2016
Functions for running a population of models loop

@author: Oliver Britton
"""
import os
import datetime
import time
import re
#import pdb
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from IPython.display import clear_output

#from multiprocessing import Pool
from multiprocessing import Process
from functools import partial
import NeuronProjectStart
import Methods.Biomarkers.NeuronBiomarkers as nb
import Methods.Biomarkers.DavidsonBiomarkers as db
import Methods.simulation_helpers as sh
from neuron import h
import neuron

# Load ion channel models
import Methods.Simulations.loadneuron as ln
ln.load_neuron_mechanisms(verbose=False)

# FUNCTIONS

def load_file(filename):

    with open(filename,'rt') as f:
        text = f.readlines()
    return text

def load_parameters(filename):
    parameters = pd.read_csv(filename, sep=',', header=None)
    return parameters
    
" To do - replace with pandas "
def ReadParameterFile(filename):
    print "ReadParameterFile deprecated - use load_parameters instead"
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
    assert False, "Deprecated - use read_trace instead of ReadTraceFile"
    
def read_trace(filename, skiprows=0):
    # Refactored to use numpy from pom.ReadTraceFile,
    data = np.loadtxt(filename,skiprows=skiprows)
    trace = {'t':data[:,0], 'v':data[:,1]}
    return trace
    
def PlotTraceFile(filename):
    assert False, "Deprecated - use plot_trace"
    
def plot_trace(filename):
    data = read_trace(filename)
    plt.plot(data['t'],data['v'])
        
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

        na_accu 
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
    
    configFile = load_file(configFilename)
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
            numProtocols += 1
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

class PopulationOfModels(object):

    def __init__(self, 
        name,
        sim_protocols,
        model_details,
        parameter_filename=None, 
        parameter_set_details=None,
        biomarker_filename=None,
    ):
        """
        Steps for running a population of models simulation:
        1. Read in parameters or generate and save a new set
        2. Set up simulation protocols (initialise from inputs)
        3. For each parameter set:
            4. Build a model
            5. Run the model using simulation protocols
            6. Calculate biomarkers from each simulation
            7. Collect all biomarkers from that model and add to results
        8. Collate and possibly save all biomarkers and any required simulation data
        9. (Optional) Analyse biomarkers and present summary data and/or visualisations if asked for
        """
        " To dos: "
        " 1. Support non-conductance parameters through updating sh.build_model "
        " 2. Parallelised simulations. "
        " 3. Support multi-compartment models and simulations (possibly in a separate class?). "
                
        # Set population name and mechanism and parameter names
        self.name = name
        self.model_details = model_details
        self.setup_mechanisms() # Get mechanism names, parameter names, parameter designations and parameter ranges

        # Load parameter values
        if (parameter_filename is not None) & (parameter_set_details is not None): # Check only one option for parameters is set
            raise ValueError("Both parameter filename and parameter details have been provided, choose one or the other.")
            
        self.setup_parameters(parameter_filename=parameter_filename, parameter_set_details=parameter_set_details)
        self.model_description = self.get_model_description()
        
        # Initialise results data structures
        self.setup_results()
        self.traces = {}
       
        # Setup calibration
        self.calibration_complete = False # Has calibration been performed?
        self.calibration_applied = False # Has calibration been applied on to results?
        self.calibration = None # Dataframe showing of each parameter set and each biomarker showing whether each parameter set passeed calibration for each biomarker
        self.calibrated_indices = None # Indices of parameter sets that have passed calibration to for all tested biomarkers
            
        # -- Setup simulation protocols --
        # Stimulus protocols
        self.setup_simulation_protocols(sim_protocols=sim_protocols)
        self.current_set = None
        self.celsius = 32. # In line with Davidson et al., 2014, PAIN
        
        # Drug blocks and other stuff
        " To do! "
        self.blocks = {}
        # Blocks format - {'Simulation name': {'GNa':0.5, 'GKdr':0.1}}
        # To do - add in blocks to build_sim_protocols function and here when necessary
            
        
    " --- Setup functions --- "
    def setup_results(self):
        """
        Setup the main results dataframe to store parameters, biomarkers and model metadata. Load parameters and model information into the dataframe as part of the setup process.
        """
        self.results = self.parameters.copy(deep=True)
        param_names = self.results.columns
        arrays = [['Parameters']*len(param_names), param_names]
        columns = pd.MultiIndex.from_arrays(arrays, names=['',''])
        self.results.columns = columns
        self.results.index.name = 'Model'

    
    def setup_mechanisms(self):
        " Load mechanism details from model details "
        self.mechanisms = self.model_details['mechanisms']
        self.mechanism_names = self.mechanisms.keys() # gives the list of mechanisms to insert
        self.parameter_names = [] # Public-facing name of each parameter (e.g. GNav17, GKdr)
        self.parameter_designations = {} # Map of parameter names to designations of each parameter in neuron mod file
        for mechanism in self.mechanisms:
            self.parameter_names.extend(self.mechanisms[mechanism].keys())
            self.parameter_designations.update(self.mechanisms[mechanism])
            
        assert len(self.mechanism_names) == len(set(self.mechanism_names)), "Mechanism names are not all unique."
        assert len(self.parameter_names) == len(set(self.parameter_names)), "Parameter names are not all unique."
        assert len(self.parameter_designations) == len(set(self.parameter_designations)), "Parameter designations are not all unique."
        
    
    def setup_parameters(self, parameter_filename=None, parameter_set_details=None):
        """ Load or generate a parameter set for simulations. """    
        # If parameter_filename is provided, load parameters from there "
        if parameter_filename: 
            parameters, header = self.load_parameter_set(parameter_filename, load_comment=True, comment='#')    
        # If there is no parameter filename and parameter set details have been provided, generate a new parameter set using the parameters, and optionally save it. "
        elif parameter_set_details:
            num_models = parameter_set_details['num_models']
            parameter_data = parameter_set_details['parameter_data']
            save = parameter_set_details['save']
            output_filename = parameter_set_details['output_filename']
            if isinstance(parameter_data, list):
                minimum = parameter_set_details['minimum']
                maximum = parameter_set_details['maximum']
            else:
                minimum = None
                maximum = None
            parameters, header = sh.build_parameter_set(num_models=num_models, parameter_data=parameter_data, minimum=minimum, maximum=maximum, filename=output_filename, save=save)
        else:
            raise ValueError('No filename to load or details provided to set up parameters.')
        # Finally, set parameters
        self.parameters = parameters
        self.parameter_details = header
    
    def setup_simulation_protocols(self, sim_protocols):
        
        if 'biomarker_protocol' in sim_protocols.keys():
            self.biomarker_protocol = sim_protocols['biomarker_protocol']
        else:
            self.biomarker_protocol = None
        
        self.delay = sim_protocols['delay']
        self.amp = sim_protocols['amp']
        self.dur = sim_protocols['dur']
        self.interval = sim_protocols['interval']
        self.num_stims = sim_protocols['num_stims']
        self.stim_func = sim_protocols['stim_func']
        self.t_stop = sim_protocols['t_stop']
        self.v_init = sim_protocols['v_init']
        self.currents_to_record = sim_protocols['currents_to_record']
    
    def setup_current_recording(self):
        """ Temp function to remind me to do current recording using exec somehow """
        if self.currents_to_record:
            self.currents = {}
            for current in self.currents_to_record:
                currents[current] = h.Vector()
                # Use exec to record the right current variable (e.g. ik_kdrtf)
                exec('currents[current].record(self.active_model(0.5)._ref_{}, sec=active_model'.format(current))
    
    def load_calibration_ranges(self, calibration_ranges=None, calibration_filename=None):
        """ Load calibration ranges from a dataframe or a file """
        if calibration_filename:
            calibration_ranges = pd.read_csv(calibration_filename)
            self.calibration_ranges = calibration_ranges
        elif calibration_ranges:
            assert type(calibration_ranges) == pd.DataFrame, "Supplied calibration ranges are not a pandas DataFrame"
            self.calibration_ranges = calibration_ranges
        else:
            raise ValueError("No calibration range or filename supplied.")
    
    def get_model_description(self):
        """
        Construct a string describing the model attached to a parameter set.
        Model description format:
        'Population of models with currents: x, y and z, varying parameters i, j, k, l,m,n,o. Population class constructed on -date-.'
        """
        mechanisms = ', '.join(self.mechanism_names)
        parameter_details = self.parameter_details
        date = time.strftime("%d/%m/%Y")

        model_description = 'Population of models with mechanisms: {0}. Parameters: {1}. Population class initialized on {2}.'.format(mechanisms, parameter_details, date)
        return model_description
            
    " --- Simulation functions --- "
    
    def run_simulations(self, simulation_name='Standard', sampling_freq_hz=20000, plot=False, save=False):
        """
        Runs simulations on all models in results, whether that is an initial sample or a calibrated population.
        
        """
        # TODO: Add parallelisation here - see parallelisation in https://github.com/mirams/PyHillFit/blob/master/python/PyHillFit.py
        
        # Check simulation isn't already here
        if simulation_name in self.results.columns.levels[0]:
            raise ValueError('simulation_name is already present in results.')
        else:
            biomarker_results = pd.DataFrame(index=self.results.index)
        parameters = self.results['Parameters']

        start_time = time.time()
        num_sims = len(self.results.index)
        if plot:
            plt.figure(figsize=(15,15))
        if save:
            self.traces[simulation_name] = {}
            

        for i, model_idx in enumerate(self.results.index):   
            now = time.time()
            output_frequency = 100. # Number of outputs to make - 100 equals every 1%
            reporting_period = np.ceil(len(self.results.index)/output_frequency)
            if i > 0:
                if (i == len(self.results.index)-1) | (i%reporting_period == 0):
                    clear_output()
                    print("Simulation {} of {} started. Time taken = {:.1f}s. Estimated remaining time = {:.1f}s.".format(i+1, num_sims, now-start_time, (num_sims-i)*(now-start_time)/i))
            else:
                print("Simulation set of {} simulations begun.".format(num_sims))
                
            # --- Setup model ----
            self.active_parameters = parameters.loc[model_idx]
            mechanisms = {self.parameter_designations[param]:val for param, val in self.active_parameters.iteritems()} # iteritems over pd.Series to get key value pairs and generate a dict from them
            # To do - allow different cellular parameters and ionic concs to be input here or in bespoke_simulation, as well as variation in parameters
            
            # Build model using dict with full parameter names (e.g. gbar_nav17vw not nav17vw)
            self.active_model = sh.build_model(mechanisms=mechanisms, mechanism_names=self.mechanism_names, conductances=None, mechanism_is_full_parameter_name=True)

            
            # Set temperature
            h.celsius = self.celsius # In line with Davidson et al., 2014, PAIN
            # Set ionic conditions (update reversal potentials in line with ionic conditions)
            # To do - check which ions exist and only set conditions for them
            oldstyle = h.ion_style("k_ion", 1, 2, 1, 1, 0,sec=self.active_model)
            oldstyle = h.ion_style("na_ion", 1, 2, 1, 1, 0,sec=self.active_model)
            
            # --- Run simulations ---
            
            # Rheobase simulation
            rheobase = nb.CalculateRheobase(self.active_model, amp_step=0.1, amp_max=5, make_plot=False,)
            
            # Only continue if we've found a rheobase
            if ~np.isnan(rheobase):
                # All other biomarkers simulation at rheobase
                # Set stimuli
                stims = []
                for stim_idx in range(self.num_stims):
                    stims.append(self.stim_func(0.5, sec=self.active_model))
                    stims[-1].dur = self.dur
                    stims[-1].amp = rheobase
                    stims[-1].delay = self.delay + stim_idx*(self.dur + self.interval)
                
                " To do - roll these into the pom class (sh.set_vt and sh.record_currents) "
                #v,t = sh.set_vt(cell=self.active_model)
                v = h.Vector()
                v.record(self.active_model(0.5)._ref_v, sec=self.active_model)
                t = h.Vector()
                t.record(h._ref_t)
                currents = sh.record_currents(cell=self.active_model, current_set=None)
                
                h.finitialize(self.v_init) # Vital! And has to go after record
                
                # --- Run simulation ---
                neuron.run(self.t_stop)
                v,t = np.array(v), np.array(t)
                if save:
                    self.traces[simulation_name][model_idx] = {'t':t, 'v':v}
                
                if plot:
                    #plt.subplot(10,10,i+1)
                    plt.plot(t,v)
                    plt.title("Model: {}".format(model_idx))
                    
                # Sampling
                if sampling_freq_hz == 20000:
                    " To do - use delta t between each element of t to calculate frequency, then downsample to required frequency. But at the moment we just need to match Davidson et al. (20 kHz). "
                    t = t[::2]; v = v[::2] # 20 kHz
                else:
                    raise ValueError("Sampling frequencies other than 20 kHz not supported yet.")
                    
                # --- Analyse simulation for biomarkers ---
                
                traces = nb.SplitTraceIntoAPs(t,v)
                biomarkers = nb.calculate_simple_biomarkers(traces,self.active_model)          
            else: # No rheobase found
                biomarkers = {}
                
            # RMP
            rmp_t,rmp_v = sh.simulation(amp=.0,dur=3000,delay=0,interval=0,num_stims=1,t_stop=3000.,make_plot=False,model=self.active_model)
            # Could change sampling freq here but probably not necessary for RMP
            rmp_traces = nb.SplitTraceIntoAPs(rmp_t,rmp_v)
            biomarkers['RMP'] = np.mean(nb.CalculateRMP(rmp_traces))
            # Add rheobase
            biomarkers['Rheobase'] = rheobase
            
            # Save biomarkers into results
            for biomarker in biomarkers:
                biomarker_results.loc[model_idx, biomarker] = biomarkers[biomarker]
        
        # At end of all simulations, build a multiarray with the simulation name and the biomarker results.
        self.biomarker_results = pd.DataFrame(columns=pd.MultiIndex.from_product([[simulation_name], biomarker_results.columns]), index=self.results.index)
        # Insert results
        self.biomarker_results.loc[:] = biomarker_results
        # Then concatenate it into the main results dataframe.
        self.results = pd.concat([self.results, self.biomarker_results], axis=1)
        
        print("Simulation {} complete in {} s.".format(simulation_name,time.time()-start_time))
            
    def run_bespoke_simulation(self, sim_name, sim_protocols, biomarkers_to_calculate):
        """
        Run a non-standard simulation, e.g. with drug block, or supra-threshold stimuli and save
        """
        
        # Build model, checking for optional other parameters in sim protocols
        
        # Setup stim protocols from sim_protocols
        
        # Do the simulation on the current parameter sets
        
        # Calculate the biomarkers
        
        # Save to results under "sim_name"
    
    " --- Calibration functions --- "
    
    def calibrate_population(self, biomarker_names, simulation_name, calibration_ranges='Davidson', stds=None, verbose=True):
        """ 
        Calibrate the current parameter sets, saving the data on calibration criteria passing to a dataframe which is set as the
        active calibration.
        """
        
        # First check that we have a set of results (parameters and biomarkers) that contain data for the appropriate biomarkers under the right simulation conditions
        
        assert type(simulation_name) == str # We're not supporting multiple simulation names yet
        
        """
        Calibration process:
        1. Copy the appropriate biomarkers
        2. Get the appropriate ranges
        3. Perform the calibration (see acc.calibrate for calibration function)
        4. Pick out the indices of the models that passed all calibration checks.
        5. Save the calibration results so they can be visualised. 
        """
        # Perform calibration using method defined by inputs (currently only supporting Davidson std dev based calibration
        if calibration_ranges == 'Davidson':
            if stds == None:
                stds = 1.
            ranges = db.CalibrationData(num_stds = stds)
            
            results = self.results.copy(deep=True)
            # Build dataframe to store calibration results
            calibration = pd.DataFrame(index=results.index, columns=biomarker_names)  # If we allow multiple calibrations this dataframe needs to include multiple subsections
            
            for biomarker in biomarker_names:
                minimum = ranges.loc[biomarker]['Mean'] - ranges.loc[biomarker]['Std']
                maximum = ranges.loc[biomarker]['Mean'] + ranges.loc[biomarker]['Std']
                
                # Get calibration stats (from original - self.results, as results copy is modified)
                if verbose:
                    num_in_range = ((self.results[simulation_name, biomarker] >= minimum) & (self.results[simulation_name, biomarker] <= maximum)).sum()
                    num_over =  (self.results[simulation_name, biomarker] >= maximum).sum()
                    num_under = (self.results[simulation_name, biomarker] < minimum).sum()
                    num_nans = np.isnan(self.results[simulation_name, biomarker]).sum()
                    total = num_in_range + num_over + num_under + num_nans
                    print("Biomarker: {}. Num in range: {}. Num over max: {}. Num under min: {}. Nans: {}. Total {}.".format(biomarker, num_in_range, num_over, num_under, num_nans, total))
                    
                # Do the actual calibration
                single_biomarker_calibration = ((results[simulation_name, biomarker] >= minimum) & (results[simulation_name, biomarker] <= maximum))
                results = results[single_biomarker_calibration]
                calibration[biomarker] = single_biomarker_calibration # If we allow multiple calibrations this dataframe needs to include multiple subsections
                
        else:
            raise ValueError('Calibration_ranges value not understood.')

        # Move all results to storage and bring in the new calibrated results. 
        if self.calibration_complete == False: # First calibration - so results are fresh
            self.uncalibrated_results = self.results.copy(deep=True)
            self.calibration_complete = True 
        self.calibration = calibration # Dataframe showing of each parameter set and each biomarker showing whether each parameter set passeed calibration for each biomarker
        self.calibrated_indices = results.index # Indices of parameter sets that have passed calibration to for all tested biomarkers
        self.results = results
        self.calibration_applied = True
        return

        
    def calibrate(self, biomarkers,ranges):
        """ Simple function to perform calibration of a set of biomarkers against a set of ranges.
        biomarkers is a pd.DataFrame of biomarker values 
        ranges is a pd.Dataframe of biomarker ranges 
        output is a boolean pd.Series for whether each biomarker set was within all ranges or not 
        """
        calibration = pd.Series(True,index=biomarkers.index) 
        for biomarker in biomarkers: 
            calibration  = calibration & (biomarkers[biomarker] <= ranges[biomarker].loc['max'])
            calibration  = calibration & (biomarkers[biomarker] >= ranges[biomarker].loc['min'])                
        return calibration
        
    def revert_calibration(self, preserve_calibration=False):
        " Revert population to uncalibrated state, optionally saving the calibration and calibrated results in a separate dataframe "
        if self.calibration_applied & self.calibration_complete:
            if preserve_calibration:
                self.calibrated_results = self.results.copy(deep=True)
                self.reverted_calibration = self.calibration.copy(deep=True)            
            self.results = self.uncalibrated_results.copy(deep=True) # Recover uncalibrated results
            # Clean up
            self.calibration = None
            self.uncalibrated_results = None
            self.calibration_complete = False
            self.calibration_applied = False
            
        else:
            raise ValueError('Calibration is not complete and applied.')
            
    " --- Analysis functions --- "
    def plot_parameters(self):
        pass
        
    def plot_biomarkers(self,name):
        pass
        
    def something_with_clustering(self):
        pass
        
    def parameter_correlations(self):
        pass
    
    
    " --- Storage functions --- "
    def pickle_pom(self, filename=None):
    # Serialise pom object as a pickle 
        if filename == None:
            filename = '{}.pickle'.format(self.name)
        # Remove hoc objects
        self.active_model = None
        if self.stim_func == h.IClamp:
            self.stim_func = 'IClamp'
        elif self.stim_func == h.IRamp:
            self.stim_func = 'IRamp'
        else:
            self.stim_func = 'None'
            print("Stim func type not found.")
        # Then pickle
        with open(filename, 'wb') as f:
            pickle.dump(self,f)
            
    
    def load_parameter_set(self, filename, load_comment=False, comment='#'):
        """
        Load a parameter set CSV file and return as a dataframe. Optionally, also
        return any commented lines at the top of the file. 
        """

        parameters = pd.read_csv(filename, sep=',', comment=comment,index_col=0)  
        if load_comment:
            comments = ''
            with open(filename, 'r') as f:
                while True:
                    line = f.readline()
                    if line[0] == comment:
                        comments += ' {}'.format(line)
                    else:
                        break # Just return commented lines
            return parameters, comments
        else:
            return parameters
        
    def apply_calibration(self):
        """
        Update self.results to only include the results for the calibrated parameters. Store
        the old results in self.uncalibrated_results.
        """
        if self.calibration_complete:
            self.uncalibrated_results = self.results.copy(deep=True)
            self.results=self.results.loc[self.calibration.all(axis=1)]
            self.calibration_applied = True
        else:
            print("Calibration not performed yet.")

    def remove_calibration(self):
        """
        Remove the last applied calibration and store the calibrated results if required.
        Applying then removing a calibration should return self.results back to where it started, however calibrated_results and uncalibrated_results will now exist as separate dataframes. This may be unwanted behaviour in the future, so there may need to be a cleanup function.
        """
        if self.calibration_applied:
            self.calibrated_results = self.results.copy(deep=True)
            self.results = self.uncalibrated_results.copy(deep=True)
            self.calibration_applied = False
        else:
            print("Calibration not currently applied to results.")
        
        
def load(filename):
    # Load population of models from a pickled file
    with open(filename, 'rb') as f:
        pom = pickle.load(f)
    if pom.stim_func == 'IClamp':
        pom.stim_func = h.IClamp
    elif pom.stim_func == 'IRamp':
        pom.stim_func = h.IRamp
    return pom
    
def make_pom_name(name):
    # Make the simulation name
    date = time.strftime("%d%m%y")
    return '{}_{}'.format(name,date)
        
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
    calibrationBiomarkerFile = open(outputDirectory + prefix + '_biomarkers.dat','w') # File for rheobase biomarkers only for calibration
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
