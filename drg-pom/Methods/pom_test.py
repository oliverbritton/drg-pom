# -*- coding: utf-8 -*-
"""
Created on Mon Feb 08 12:14:25 2016
Functions for running a population of models loop

@author: Oliver Britton
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='whitegrid') # For nicer plots
import multiprocessing as mp

import collections
import datetime
import pickle
import time

from functools import partial
from IPython.display import clear_output


import NeuronProjectStart
import Methods.Biomarkers.NeuronBiomarkers as nb
import Methods.Biomarkers.DavidsonBiomarkers as db
import Methods.simulation_helpers as sh
from neuron import h
import neuron

# We seem to need to load now to allow parallel scripts to work
# Perhaps there is a way to reload if we need to?
import Methods.Simulations.loadneuron as ln
ln.load_neuron_mechanisms(verbose=True)

# ---FUNCTIONS---

def load_parameters(filename):
    parameters = pd.read_csv(filename, sep=',', header=None)
    return parameters
            
            
def read_trace(filename, skiprows=0):
    # Refactored to use numpy from pom.ReadTraceFile,
    # read_trace and load_trace are synonyms for the same function
    data = np.loadtxt(filename,skiprows=skiprows)
    trace = {'t':data[:,0], 'v':data[:,1]}
    return trace
  
load_trace = read_trace # use either name to load/read traces

   
def save_trace(trace, filename):
    """ TODO - make this able to handle current recordings - e.g. save as csv with column headers
    using pandas
    """
    """
    if type(trace) == dict:
        data = np.column_stack(trace['t'], trace['v'])
        #data = np.column_stack(data,
    elif type(trace) == np.array:
        data = trace # assume formatting is correct
    elif (type(trace) == tuple) | (type(trace) == list): # assume list in order: t,v
        data = np.column_stack([trace[i] for i in range(len(trace))])
    else:
        raise TypeError("Can't parse type of trace")
    np.savetxt(filename, data)
    """
    # Save as a pickle
    with open(filename, 'wb') as f:
        pickle.dump(trace,f)
        
        
def plot_trace(filename):
    data = read_trace(filename)
    plt.plot(data['t'],data['v'])
             
             
def load(filename):
    # Load population of models from a pickled file
    with open(filename, 'rb') as f:
        pom = pickle.load(f)
    """
    if pom.stim_func == 'IClamp':
        pom.stim_func = h.IClamp
    elif pom.stim_func == 'IRamp':
        pom.stim_func = h.IRamp
    """
    return pom
    
    
def make_pom_name(name):
    # Make the simulation name
    date = time.strftime("%d%m%y")
    return '{}_{}'.format(name,date)

def allowed_save_types():
    return ('fig', 'trace', 'both', 'all')
    
def process_save_type(save_type):
    if save_type in ['fig', 'trace']:
        return [save_type]
    elif save_type == 'both':
        return ['fig', 'trace']
    elif save_type == 'all':
        return ['fig', 'trace']
    else:
        raise ValueError('save_type {} not found'.format(save_type))
    
def get_function_args(function):
    """ Just a reminder of how to get the args of a function """
    import inspect
    #@PYTHON3 - change to getfullargspec
    args  = inspect.getargspec(function)[0]
    return args
   
"""   --- SIMULATION FUNCTIONS ---
 These are set outside of the Simulation class because of multiprocessing.
 Functions sent to other processes have to be serializable (picklable), and 
 member functions aren't. So the functions are at the top level of pom's scope, 
 which makes them picklable.
"""

def simulate_iclamp(sim_id,
                                        biomarker_names,
                                        mechanisms, 
                                        mechanism_names,
                                        ions,
                                        t_stop,
                                        outputs,
                                        amp,
                                        dur,
                                        delay,
                                        interval,
                                        num_stims,
                                        sampling_freq,
                                        stim_func,
                                        v_init,
                                        celsius,
                                        options,
                                        metadata,
                                        flags,):
    """
    Runs a full simulation in IClamp mode (finds rheobase, finds rmp, and calculates biomarkers at rheobase)
    
    Parameters
    ----------------
    sim_id (int) - identifier for the simulation
    biomarker_names (list of strings) 
    mechanisms (dict) - mechanism details to construct model channels from
    mechanism_names - (list of strings) -  ion channel names
    ions (list of strings)
    t_stop (float, ms)
    outputs (list of strings) - options to specify output including 'ik' (total k+ current biomarkers)
    amp (float, nA)
    dur (float, ms)
    delay (float, ms)  
    interval (float, ms)       
    num_stims (int)       
    sampling_freq (float, Hz)
    stim_func (hoc stim function - e.g. h.IClamp, or string, e.g. 'h.IClamp')
    v_init (mV)
    celsius (oC)
    options (dict)
    metadata(dict)
    flags(dict) - flags for special options such as passing extra data to the simulation
    
    Returns
    -----------
    """
    # General setup

    #with open('test.txt','w') as f:
    #    f.write(str(mechanisms))
    import neuron
    from neuron import h
    sim_type = 'iclamp'
    
    # Need to use strings as can't serialize hoc objects
    # for parallelization
    if stim_func == 'h.IClamp':
        stim_func = h.IClamp
    elif stim_func == 'h.IRamp':
        stim_func = h.IRamp
    else:
        raise ValueError('stim_func: {} not found'.format(stim_func))

    # Initialise results @TODO - shared code between simulations - combine into function
    results = {}
    results['sim_id'] = sim_id
    results_names = biomarker_names 
    for name in results_names:
        results[name] = np.nan # @TODO - do this better to support biomarkers with variable lengths
    results['Firing pattern'] = np.nan
   
    """ Setup model """    

    # Build model with separate mechanisms and mechanism parameter names
    cell = sh.build_model(mechanisms=mechanisms, mechanism_names=mechanism_names, conductances=None, mechanism_is_full_parameter_name=True)
    h.celsius = celsius

    # @TODO - put setting ionic conditions into a function to share with VClamp
    # Set ionic conditions
    ion_mechanisms = sh.get_dynamic_ion_mechanisms()   
    for ion in ion_mechanisms:
        if ion in ions:
            cell.insert(ion_mechanisms[ion])
 
    
    # --- Run simulations ---
    # Assemble simulation kwargs
    sim_kwargs = {
        'ions':ions,
        't_stop':t_stop,
        'dur':dur,
        'delay':delay,
        'interval':interval,
        'num_stims':num_stims,
        'stim_func':stim_func,
        } 

    
    # Rheobase simulation
    rheobase = nb.calculate_rheobase(cell, amp_step=0.1, amp_max=5.0, 
        make_plot=False, sim_kwargs=sim_kwargs)

    # Courser grained rheobase if the first one fails. Amp max of 50 corresponds
    # to maximum ramp stimulus amplitude used in the Nav 1.8 study.
    if rheobase is nb.RHEO_FAIL:
        rheobase = nb.calculate_rheobase(cell, amp_step=1.0, amp_max=50.0,
            make_plot=False, sim_kwargs=sim_kwargs)
   
    # Only continue if we've found a rheobase
    if rheobase is not nb.RHEO_FAIL:
        rheobase_found = True
    else:
        rheobase_found = False # Or set all other biomarkers to np.nan?
    
    # TODO - just call sh.simulation(**sim_kwargs), adding extra args to those above as needed
    if rheobase_found or (amp is not None) or ("rheobase_sim_for_stim_amp" in flags):    
        stims = []
        for stim_idx in range(num_stims):
            stims.append(stim_func(0.5, sec=cell))
            stims[-1].dur = dur
            stims[-1].delay = delay + stim_idx*(dur + interval)
            # Set amp
            if amp is None:
                # Normal condition: rheobase from this simulation is amp
                if "rheobase_sim_for_stim_amp" not in flags:
                    stims[-1].amp = rheobase # If amp is None use rheobase as amp
                # Alternatively if flag is set: use the rheobase from another simulation as amp
                else:
                    stims[-1].amp = flags["rheobase_sim_for_stim_amp"]
            else:
                # Predetermined amplitude
                stims[-1].amp = amp
 
        v,t = sh.set_vt(cell=cell)
        vectors = sh.record_currents(cell, outputs)

        h.finitialize(v_init) # Vital! And has to go after record
        neuron.run(t_stop)

        # Clean up outputs after simulation
        v,t = np.array(v), np.array(t)

        # Recast any recorded currents so they're not mutable
        # by future simulations
        vectors = sh.recast_recorded_currents(vectors)

        # Sampling for v, currents, concs
        if sampling_freq == 20000: # Hz
        # TODO - use delta t between each element of t to calculate frequency,
        # then downsample to required frequency.
        # But at the moment we just need to match Davidson et al. (20 kHz).
        # Neuron records at 40 kHz so just take every second element
            t = t[::2]; v = v[::2] # 20 kHz

        # Currents
            for cur_name, cur in vectors.items():
                for cur_component in cur:
                    vectors[cur_name][cur_component] = vectors[cur_name][cur_component][::2]

            # To do - concs
        else:
            raise ValueError("Sampling frequencies other than 20 kHz not supported yet.")
        
        # --- Analyse simulation for biomarkers ---
        traces = nb.split_trace_into_aps(t,v)
        # This is the borked bit - TODO: fix this
        biomarkers = nb.calculate_simple_biomarkers(traces, cell)          
       
        try:
            # Check flags for alternate biomarker calculations
            for flag in flags:
                if flag == "ramp_threshold_sim_for_width":
                    _threshold = flags[flag]
                    biomarkers['APFullWidth'] = nb.average_biomarker_values(
                            nb.calculate_ap_width(traces, alpha=0.0, threshold=_threshold, method='voltage'),
                            how_to_handle_nans="return")
                    biomarkers['APHalfWidth'] = nb.average_biomarker_values(
                            nb.calculate_ap_width(traces, alpha=0.5, threshold=_threshold, method='voltage'),
                            how_to_handle_nans="return")
        except Exception as e:
            with open("{}_{}.txt".format(sim_id,stim_func),'w') as f:
                f.write(str(e))
                biomarkers['APFullWidth'] = "FAILURE"
                biomarkers['APHalfWidth'] = "FAILURE"

        #print("biomarkers.keys: {}".format(biomarkers.keys())) #debug
        #print("results.keys: {}".format(results.keys())) # debug

        # Firing pattern
        
        if num_stims == 1:
            stim_start = delay
            stim_end = delay + dur
            try:
                biomarkers['Firing pattern'] = nb.determine_firing_pattern(traces, stim_start, stim_end)
            except Exception as e:
                with open("{}_{}".format(sim_id,str(stim_func)),'w') as f:
                    f.write(str(e))
        else:
            biomarkers['Firing pattern'] = "Couldn't determine as num stims = {}, not 1".format(num_stims)
        
        for result in biomarkers:
            results[result] = biomarkers[result]
            
        # Build trace for output
        trace = {'t':t, 'v':v}
        for vector in vectors:
            trace[vector] = vectors[vector]
            
    elif rheobase_found == False:
        trace = None # pass an empty trace
        # We have already set the biomarkers to np.nan by default
    else: 
        raise ValueError('rheobase_found is not valid')

    # RMP
    rmp_out = sh.simulation(amp=0.0,dur=3000,delay=0,interval=0,num_stims=1,t_stop=3000.,make_plot=False,model=cell)
    
    # Could change sampling freq here but probably not necessary for RMP
    rmp_traces = nb.SplitTraceIntoAPs(rmp_out['t'],rmp_out['v'])
    
    results['RMP'] = np.mean(nb.CalculateRMP(rmp_traces))

    # Add rheobase
    results['Rheobase'] = rheobase

    # Output @TODO - wrap into function as code is mostly shared with simulate_vclamp
    plot = options['plot']
    save = options['save']
    
    if save == True:
        save_type = options['save_type']
        save_types = process_save_type(save_type)
        if 'fig' in save_types:
            # Save results so we can aggregate multiple simulations and plot
            results['trace'] = trace
        if 'trace' in save_types:
            filename = '{}.pickle'.format(sim_id)
            for name in ['simulation_name', 'population_name']: 
            # Last element in list will be the front-most part of the filename
                if name in metadata.keys():
                    filename = '{}_{}'.format(metadata[name], filename)
            save_trace(trace, filename)

    if plot: #@TODO 
        pass
    dump_exception(results)
    return results


def simulate_vclamp(sim_id,
                                    biomarker_names,
                                    mechanisms,
                                    mechanism_names,
                                    ions,
                                    t_stop,
                                    outputs,
                                    hold,
                                    steps,
                                    delta,
                                    durations,
                                    celsius,
                                    options,
                                    metadata,
                                    flags,):
    """
    Simulation of standard voltage clamp protocol.
    Example implementing following protocol for K+ recordings from Anabios:
    mech = {'kdrtf': 1.0, 'nav18hw': 2.0,}
    mech_names =['kdrtf', 'nav18hw']
    simulate_vclamp(mechanisms=mech, mechanism_names=@TODO, hold=-70.0, steps=[-80,80], delta=10.0, durations=[1100,500,500], outputs=['v', 'ik'])

    1. Hold at -70 mV for 1100 ms
    2. Pulse from -80 to 80 mV with delta=10mV
    3. Hold at -70 mV for 500 ms
    
    Parameters
    ----------------
    sim_id (int) - identifier for the simulation
    biomarker_names (list of strings) 
    mechanisms - dict of mechanism details to construct model channels from
    mechanism_names - names of each mechanism (e.g. ion channel names)
    ions (list of strings) - ions needed in model (e.g. ['Na', 'Ca', 'K'])
    t_stop (float, ms)
    outputs (list of strings) - options to specify output including 'ik' (total k+ current biomarkers)
    hold (float, mV) - holding potential
    delta (float, mV or None) - size of steps or None to specify steps manually
    steps (list of floats, mV or list) start and end of voltage clamp if delta is specified or list of every step
    durations (3 element list of floats, ms) - start of each stage of voltage clamp 
    celsius (oC)
    options (dict)
    metadata(dict)
    flags(dict) - flags for special options such as passing extra data to the simulation
    
    Returns
    -----------
    results - dict, contents - depends on "outputs"
    results['sim_id'] - identifier for the simulation
    results['ik_steps_absmax'] - array of absolute peak k+ currents at each step
    @TODOs
    Biomarkers of curve fitting
    """
    # General setup

    import neuron
    from neuron import h
    
    sim_type = 'vclamp'
    
    # Output options
    plot = options['plot']
    save = options['save']
    if save == True: 
        save_type = options['save_type']
    
    # Setup voltage clamp steps
    if delta == None:
        assert len(steps) > 0 # Use steps as given 
    else:
        steps = np.arange(steps[0], steps[1]+(0.1*delta),delta)
    
    # Initialise results
    results = {}
    results['sim_id'] = sim_id
    results_names =  biomarker_names 
    for name in results_names:
        results[name] = np.nan # @TODO - do this better to support biomarkers with variable lengths
    
    
    
    # Initialise storage for biomarkers in vclamp loop
    if 'ik' in outputs: 
        ik_steps_absmax = np.zeros(len(steps))
        ik_steps_absmax[:] = np.nan  
     # Simulation loop through steps
    for i, step in enumerate(steps):
        cell = sh.build_model(mechanisms=mechanisms, 
                                            mechanism_names=mechanism_names, 
                                            conductances=None, 
                                            mechanism_is_full_parameter_name=True)
        
        h.celsius = celsius

        # @TODO - put setting ionic conditions into a function to share with IClamp
        # Add ionic concentration mechanisms
        if 'K' in ions:
            oldstyle = h.ion_style("k_ion", 1, 2, 1, 1, 0,sec=cell)
        if 'Na' in ions:
            oldstyle = h.ion_style("na_ion", 1, 2, 1, 1, 0,sec=cell)
        if 'Ca' in ions:
            oldstyle = h.ion_style("ca_ion", 1, 2, 1, 1, 0,sec=cell)

        v,t = sh.set_vt(cell=cell)
        vectors = sh.record_currents(cell, outputs)

        clamp = h.VClamp(0.5, sec=cell)
        clamp.amp[0] = hold
        clamp.dur[0] = durations[0]

        clamp.amp[1] = step
        clamp.dur[1] = durations[1]

        clamp.amp[2] = hold
        clamp.dur[2] = durations[2]

        """
        if len(durations) > 3:
            t_stop = durations[3]
        else: 
            t_stop  = durations[0] + durations[1] + durations[2]
        """
        h.finitialize(hold) # Vital! And has to go after record
        neuron.run(t_stop)
         
        # ---Compute and store outputs---
        
        # Save output for pickling if required
        if i == 0: 
            trace = {'t':t}
        trace['v_{}'.format(step)] = v 
        for vector in vectors:
            trace['{}_{}'.format(vector,step)] = vectors[vector]
        
        # K+ biomarkers
        # @TODO - for now just return a vector of the peak current from each step
        if 'ik' in outputs:
            ik_steps_absmax[i] = nb.absmax(vectors['ik'])
        
    # @TODO! IMPORTANT. make optional output handling better and inside a function not manual and repeated 3 times
    if 'ik' in outputs: 
        results['ik_steps_absmax'] = ik_steps_absmax
            # @TODO Add in ih (na component)
            #ik = np.array(ik) + np.array(ih_na)
            
    # Output - code shared with iclamp (TODO - roll into one function)
    plot = options['plot']
    save = options['save']
    if save == True: 
        save_type = options['save_type']
        save_types = process_save_type(save_type)
        if 'fig' in save_types:
            results['trace'] = trace
        if 'trace' in save_types:
            filename = '{}.pickle'.format(sim_id)
            for name in ['simulation_name', 'population_name']: 
            # Last element in list will be the front-most part of the filename
                if name in metadata.keys():
                    filename = '{}_{}'.format(metadata[name], filename)

            save_trace(trace, filename)

    
    if plot: #@TODO 
        pass        
    
    return results
    
# ---CLASSES----
    
class PopulationOfModels(object):
    
    def __init__(self, 
        name,
        model_details,
        simulation_protocols=None,
        parameter_filename=None, 
        parameter_set_details=None,
    ):
        """
        Main class for running populations of models in NEURON
        
        Parameters
        ----------------
        name: str - name to identify population
        
        model_details: dict, format is:
        {'mechanisms':{mechanism_details}}
        where mechanism_details=
        {'neuron_name_of_current':{'name of each parameter of that current to vary':'name of that parameter in neuron'}}
        E.g.: 
        model_details = {'mechanisms':{}}
        model_details['mechanisms']['nav17vw_named'] = {'GNav17':'gbar_nav17vw_named'}
 
        sim_protocols: dict, for format see Simulation class for each simulation function
        
        NOTE: Only one of parameter_filename and parameter_set_details can be anything other than None, to uniquely define the parameter set of the population.
        
        parameter_filename: str, default None - use this to load an existing parameter set
        
        parameter_set_details: dict, format is:
            'num_models':int
            'parameter_data': list of parameter names
            e.g. ['GNav17', 'GNav18']
            OR
            dict of parameter names and sampling ranges
            e.g. {'GNav17':[0.0,0.4], 'GNav18':[0.,4.0]}
            'minimum': None (non-uniform parameter scaling) or float (uniform parameter scaling)- if None have to provide dict for parameter data
            'maximum': See minimum.
            'save': bool - whether to save parameter set independent of population
            'output_filename': string - where to save parameter set
        
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
        
        Example: see examples @TODO: make examples and example dir in source
        
        To dos: 
        1. Support non-conductance parameters through updating sh.build_model "
        2. Parallelised simulations. "
        3. Support multi-compartment models and simulations (possibly in a separate class?). 
        """
                
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
        self.simulations = {} # storage for simulation classes
        self.traces = {}
       
        # Setup calibration
        self.calibration_complete = False # Has calibration been performed?
        self.calibration_applied = False # Has calibration been applied on to results?
        self.calibration = None # Dataframe showing of each parameter set and each biomarker showing whether each parameter set passeed calibration for each biomarker
        self.calibrated_indices = None # Indices of parameter sets that have passed calibration to for all tested biomarkers
            
        # -- Setup simulation protocols --
        # Stimulus protocols
        self.simulation_protocols = simulation_protocols 
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
        """ Load mechanism details from model details dictionary """

        self.mechanisms = self.model_details['mechanisms']
        self.mechanism_names = list(self.mechanisms.keys()) # gives the list of mechanisms to insert, # List needed for python3  - can't copy dict keys
        
        # Test that each parameter variable has a corresponding mechanism 
        for mechanism, parameters in self.mechanisms.items():
            for _, param_var_name in parameters.items():
                # NEURON format for parameters is to put the mechanism name on the end of the parameter, so we check the last len(mechanism) characters of param name
                if mechanism != param_var_name[-len(mechanism):]: 
                    raise NameError('Parameter name {} is not consistent with associated mechanism name {}.'.format(param_var_name, mechanism))
        
        self.parameter_names = [] # Public-facing name of each parameter (e.g. GNav17, GKdr)
        self.parameter_designations = {} # Map of parameter names to designations of each parameter in neuron mod file
        for mechanism in self.mechanisms:
            self.parameter_names.extend(list(self.mechanisms[mechanism].keys())) 
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
    
    def setup_simulation(self, name, simulation_type, protocols, options=None, rerun=False):
        """
        Initialises a simulation 
        
        Parameters
        ----------------
        name: str, default 'sim'
        simulation_type: str, default 'IClamp', options: 'IClamp', 'VClamp', 'Test' (case insensitive)
        protocols: dict, default None, contents dependent on 'type'.
        options
        
        """
        # Check for name collision
        if (name in self.simulations.keys()):
            name_collision = True
            if rerun == False:
                raise ValueError('Simulation name {} is already present.'.format(name))              
        else:
            name_collision = False
            
        # Assemble all the details about the simulation to give to the object
        protocols['simulation_type'] = simulation_type
        self.simulations[name] = Simulation(name, protocols, options=options, population=self)
        # Finally, set simulation rerun and collision properties
        self.simulations[name].rerun = rerun
        self.simulations[name].name_collision = name_collision
    
    def run_simulation( self, 
                        name, 
                        simulation_type,
                        protocols=None,
                        cores=1,
                        plot=False, save=False, save_type='fig', 
                        benchmark=True,
                        rerun=False,
                        ):
        """
        Runs simulations on all models in results, whether that is an initial sample or a calibrated population.
        
        Example (in iPython):
        test_pop_$date$.py:
        # read in args
        
        # create population and simulation parameters
        pop.run_simulation(simulation_parameters, sim_name='Test')
        
        clean_notebook.py
        %run test_pop.py cores=4
        
        Parameters
        ----------------
        simulation_protocols: dict of protocol for simulation - see Simulation class for formats
        
        Returns
        -----------
        Nothing, concats simulation results to self.results, unless rerun is set to true,
        in which case it overwrites results.
        
        """
        
        # Set up the simulation
        
        if protocols == None:
            protocols = self.simulation_protocols
        assert(protocols != None), "Simulation protocol not set."
        
        parameters = self.results['Parameters']
        
        if rerun == False:
            self.setup_simulation(name, simulation_type, protocols, rerun=rerun) 
        elif rerun == True:
            self.simulations[name].rerun = True
            self.simulations[name].name_collision = True
            self.simulations[name].reset_simulation()
        else: 
            raise ValueError('rerun not bool')
            
        sim = self.simulations[name] # simulation object
        sim.pom_simulation(
            simulation_type=simulation_type,
            cores=cores,
            plot=plot,
            save=save, 
            save_type=save_type,
            benchmark=benchmark,
            rerun=rerun,)
        print("Sim results:\n {}".format(sim.results)) # debug
        # Add simulation results to the main pom dataframe
        sim_results = sim.results
        formatted_results = pd.DataFrame(
            columns=pd.MultiIndex.from_product([[name], sim_results.columns]), 
            index=self.results.index)
        formatted_results[:] = sim_results
        # Check for typing and set type if all entries are the same type
        """
        for result in formatted_results:
            types = [type(i) for i in formatted_results[result]
            if all([i == types[0] for i in types]):
                formatted_results[result] = formatted_results[result].astype(types[0])
        """
            
        # If no name collision add new columns, otherwise replace the old columns 
        if self.simulations[name].name_collision == False:
            self.results = pd.concat([self.results, formatted_results], axis=1)
        elif self.simulations[name].name_collision == True:
            #print("columns in self.results: {}\n columns in formatted_results: {}".format(self.res))
            _df = formatted_results[name].copy()
            self.results[name] = _df
        else:
            raise ValueError('name_collision not set')
            
        """@TODO  - if sim.results is really big and we run lots of simulations,
            we might want to delete it. We can use:
            pop.simulations['test'].results.info() to debug memory and
            pop.simulations['test'].results.memory_usage().sum() to measure it.
        """

        
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
                minimum = ranges.loc[biomarker]['Min']
                maximum = ranges.loc[biomarker]['Max']
               #minimum = ranges.loc[biomarker]['Mean'] - ranges.loc[biomarker]['Std']
               #maximum = ranges.loc[biomarker]['Mean'] + ranges.loc[biomarker]['Std']
                
                # Get calibration stats (from original - self.results, as results copy is modified)
                if verbose:
                    num_in_range = ((self.results[simulation_name, biomarker] >= minimum) & (self.results[simulation_name, biomarker] <= maximum)).sum()
                    num_over =  (self.results[simulation_name, biomarker] >= maximum).sum()
                    num_under = (self.results[simulation_name, biomarker] < minimum).sum()
                    num_nans = pd.isnull(self.results[simulation_name, biomarker]).sum() # np.isnan breaks on dtype=object so use pd.isnull instead
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
        """
        if self.stim_func == h.IClamp:
            self.stim_func = 'IClamp'
        elif self.stim_func == h.IRamp:
            self.stim_func = 'IRamp'
        else:
            self.stim_func = 'None'
            print("Stim func type not found.")
        """
        # Then pickle
        with open(filename, 'wb') as f:
            pickle.dump(self,f)
    
    save_pom = pickle_pom # Alias

    
    def save_parameter_set(self, filename, save_comment=True, comment='#'):
        """
        Save parameter set and details to a csv file
        """
        with open(filename, 'w') as f:
            comment_char = comment
            f.write('{} {}\n'.format(comment_char, self.parameter_details))
            # Round to 6 d.p., good for up to approx a million models
            # with scaling factors of order 1, and is pandas default.
            parameter_sets = self.parameters.round(6) 
            parameter_sets.to_csv(f)

    
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
                        comments += line.lstrip(' ' + comment) # Remove leading whitespace and comment char
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


    " --- Utility functions --- "
    def replace_mechanism(self, old_mech_name, new_mech_name):

        # Replace old (e.g.) wt mechanism in population with new (e.g. mutant) mechanism
        for i, mech_name in enumerate(self.mechanism_names):
            if mech_name == old_mech_name:
                self.mechanism_names[i] = new_mech_name
        old_params = self.mechanisms.pop(old_mech_name, None) # Delete old entry and get value
        new_params = {}
        for name, param in old_params.items():
            new_params[name] = param.replace(old_mech_name, new_mech_name)        
            self.parameter_designations[name] = param.replace(old_mech_name, new_mech_name)
        self.mechanisms[new_mech_name] =  new_params
        # Model details will auto update
        self.model_description = self.model_description.replace(old_mech_name, new_mech_name)
        
    def rename(self, name):
        self.name = name

    def get_simulation_names(self):
        return sorted(list(self.simulations.keys()))
    
            
class Simulation(object):
    """
    Main class for running all kinds of simulations
    
    Assumptions:
    One instance of the Simulation class will run one simulation set and not be reused.
    When run_simucalculate
    one instance of each biomarker for each parameter set.
    For each parameter set, the base model structure will be the same.
    
    Parameters
    -----------------
    sim_name: str, default 'simulation'
    population: PopulationOfModels, default None
    
    
    Examples
    -----------
    Use for a single simulation:
    @To do
    Use as part of a POM simulation:
    @To do
    To generate the format for protocols:
    sim = Simulation(name='test',protocols=None, population=pop)
    sim.empty_simulation_protocol(sim_name)
    
    """
    def __init__(self, name, protocols, options=None, population=None):
        self.name = name
        self.protocols = protocols
        self.traces = {} # To store traces for plotting or analysis
        self.plotting_pointer = 0
        self.num_subplots = 100
        
        # storage to define whether simulation is a rerun in the population or not
        self.rerun = None
        self.name_collision = None
       
        # Population-specific setup
        if population != None:
            self.population = population # Set reference to population
            self.simulation_ran = False # Check if simulation has been run - only run once
        else:
            # @TODO: Add in support for single cell simulations
            raise ValueError('Population must be supplied to run simulation, for single cell simulations use simulation_helpers module.')
            
        # Options - output, saving, logging etc.
        if options == None:
            self.options = options
        else:
            assert False, "options not enabled in constructor" # @TODO do this better
    

    def reset_simulation(self):
        """ 
        Reset simulation for a rerun.
        """
        self.traces = {}
        self.plotting_pointer = 0
        
        
    def pom_simulation(self, simulation_type, cores=1, plot=False, save=False,
            save_type='fig', benchmark=True, rerun=False):
        """
        Run simulations 
        
        Parameters
        -----------------
        simulation_type: string
        cores: int, default 1
        plot: bool, default False
        save: bool, default False
        save_type: str, default 'fig', options aredefined by allowed_save_types()
        benchmark: bool, default True
        
        Returns
        -----------
        Dataframe of biomarker results indexed by their index in the population
        associated with this simulation class.
        """
        if (self.simulation_ran == True) & (self.rerun == False):
            raise ValueError("This simulation has already been ran, cannot run same simulation again"
            "with the same simulator with name: {} without specifying simulation is a rerun.".format(self.name))
            
        if self.population == None:
            raise ValueError('No population specified - cannot run population of models simulation!')
        
        
        # Construct options and metadata and add to protocol
        assert(save_type in allowed_save_types()), "Save type not recognised."
        options = {'plot':plot, 'save':save, 'save_type':save_type}
        self.options = options # save options so we can use them to process output
        self.protocols['options'] = options
        # Metadata is used for saving and plotting
        metadata = {'population_name':self.population.name, 'simulation_name':self.name}
        self.protocols['metadata'] = metadata

        # Simulation setup
        self.set_simulation_function(simulation_type)
        num_sims = len(self.population.results.index)
        self.results = pd.DataFrame(
                index=self.population.results.index, 
                columns=self.get_biomarker_names(simulation_type, self.protocols['outputs']), 
                )
                                                       
        # ---- Main parallel simulation loop ---- 
        self.model_indices = self.population.results.index.tolist()
        self.model_indices.sort() # No matter the order in the dataframe we use a consistent order across simulations, unless sort changes.
        self.model_indices = tuple(self.model_indices) # Make order immutable

        assert(len(self.model_indices) == len(set(self.model_indices))), "Duplicates in model indices." 
        print("Simulation set of {} simulations begun.".format(num_sims))
        start  = time.time()
        pool = mp.Pool(processes=cores)
        
        #num_per_batch = 4
        #num_runs = int(np.ceil(float(len(self.model_indices)/num_per_batch))
        #for run in range(num_runs):

        self.count = 0
        for i, model_idx in enumerate(self.model_indices):   
            mechanisms = self.build_parameters(model_idx)
            # @TODO - allow different cellular parameters and ionic concs to be input, as well as variation in parameters
            
            sim_kwargs = self.build_simulation_kwargs(
                    sim_id = model_idx,
                    simulation_type = simulation_type,
                    simulation_parameters = self.protocols,
                    mechanisms = mechanisms)
            pool.apply_async(self.simulation_function, kwds=sim_kwargs, callback=self.log_result)

            # @TODO - not sure how to estimate time to completion when we're using a Pool.
            """
            if benchmark:
                output_frequency = 1 # Number of outputs to make - 100 equals every 1%
                now = time.time()
                reporting_period = np.ceil(len(self.model_indices)/output_frequency)
                if i > 0:
                    if (i == len(self.results.index)-1) | (i%reporting_period == 0):
                        #clear_output()
                        pass
                    print("Simulation {} of {} started. Time taken = {:.1f}s. Estimated remaining time = {:.1f}s.".format(i+1, num_sims, now-start, (num_sims-i)*(now-start)/i))
            """
         # Close and unblock pool, lock simulation
        pool.close()
        pool.join()

        
        
        if benchmark: print("Simulations finished in {} s".format(time.time()-start))
        
        # Clear out any remaining figs
        if (save == True) & ('fig' in process_save_type(save_type)) & (len(self.traces) > 0):
            plotted = True
            timeout_counter = 0
            while(plotted):
                filename =  '{0}_{1}_traces_{2}.png'.format(self.population.name, self.name, self.count)
                plotted = self.trace_plot(self.num_subplots, filename, force_plot=True)
                if plotted: self.count += 1
                timeout_counter += 1
                if timeout_counter > 10*num_sims/self.num_subplots: raise RuntimeError('Plotting timed out at end of simulation.') # Timeout if we wait for 10x the max number of plots
        
        # End and lock simulation
        self.simulation_ran = True


    def set_simulation_function(self, simulation_type):   
        """
        Function to set simulation 
        """
        simulation_type = simulation_type.lower()
        if simulation_type == 'iclamp':
            self.simulation_function = simulate_iclamp
        elif simulation_type == 'vclamp':
            self.simulation_function = simulate_vclamp
        else: 
            raise ValueError('simulation type: {} is not found.'.format(simulation_type))


    def build_parameters(self, model_idx):
        """
        Construct parameter data for a model including parameter changes if required
        @TODO! - Think about how we want to use mechanisms and mechanism names - is there a better format?
        @TODO - Understand how to handle multiprocessing errors so we can return an error message if the
        parameter names and mechanism names don't line up.
        @TODO - We have we have population.parameter_designations, population.mechanism_names and
        population.results['Parameters']....should we consolidate at least the parameters and mechanism
        names into a single data structure?

        Format of parameter_scaling is a dict with format:
        {parameter:scaling_factor}
        where parameter is the public parameter name shown in pop.results
        and where scaling factor is a number that the base parameter value will be multiplied by
        """

        # Make working copy of base parameter set. Important to copy otherwise we may modify population.results
        active_parameters = self.population.results['Parameters'].loc[model_idx].copy()

        # If we have parameters to scale, find the NEURON parameter name and apply scaling
        # TODO - should we use a better name than 'parameter_scaling'? Maybe allow several names like:
        # parameter_scaling, drug_block, drug block, parameter change, etc.? (But if there's more than
        # one of these names in protocols we throw an error). 
        if 'parameter_scaling' in self.protocols:
            for param, scaling_factor in self.protocols['parameter_scaling'].items():
                # Check parameter is present and apply scaling factor
                assert param in self.population.parameter_designations, (
                "Parameter {} not present in parameter designations".format(param))
                active_parameters.loc[param] *= scaling_factor

        # Finally, build mechanisms from parameters using parameter_designations to translate from public parameter
        # names to NEURON names
        mechanisms = {self.population.parameter_designations[param]:val for param, val in active_parameters.items()} 
        return mechanisms
              
    def get_simulation_protocol_names(self):
        """
        Gets list of parameter names for each simulation 
        """
        names = {}
        names['always'] = ['sim_id', 'options', 'metadata',] # always assign these for every simulation
        names['clamp'] = ['mechanisms', 'mechanism_names', 'biomarker_names', 'ions', 't_stop', 'outputs', 'celsius'] # Common to IClamp and VClamp
        names['iclamp'] = ['amp', 'dur', 'delay', 'interval', 'num_stims', 'stim_func', 'v_init', 'sampling_freq']
        names['vclamp'] = ['hold', 'steps', 'delta', 'durations']
        names['optional'] = ['flags', 'parameter_scaling'] # use these if provided but don't require them
        return names
        

    def build_simulation_kwargs(self, sim_id, simulation_type, simulation_parameters, mechanisms):
        """
        Setup dict of keyword arguments to feed to simulation function
        """
        simulation_type = simulation_type.lower() 
        assert (simulation_type in ['iclamp', 'vclamp', 'test']), 'Simulation type: {} not found'.format(simulation_type)
        
        # Special additions
        simulation_parameters['sim_id'] = sim_id
        simulation_parameters['mechanisms'] = mechanisms
        simulation_parameters['mechanism_names'] = self.population.mechanism_names
        simulation_parameters['biomarker_names'] = self.get_biomarker_names(simulation_type, simulation_parameters['outputs'])
        
        names = self.get_simulation_protocol_names()
        kwargs = {}
        
        #  Mechanism and mechanism names parameters common to both _clamps need some special case treatment
        if simulation_type in ['iclamp', 'vclamp']:
            for kwarg in names['clamp']:
                kwargs[kwarg] = simulation_parameters[kwarg]
            assert set(kwargs.keys()) == set(names['clamp']), "Parameters common to iclamp+vclamp missing."
            
        # Now do specific parameters
        for kwarg in names[simulation_type]:
            kwargs[kwarg] = simulation_parameters[kwarg]

        # Check for optional parameters and process them
        for kwarg in names['optional']:
            # Flags
            if kwarg == 'flags':
                if kwarg in simulation_parameters:
                    kwargs['flags'] = self.process_flags(simulation_parameters['flags'], sim_id)
                else:
                    kwargs['flags'] = {} # Set empty dict if no flags provided


            
        # Finally, assign the "always" names
        for kwarg in names['always']:
            kwargs[kwarg] = simulation_parameters[kwarg]
            
        assert set(kwargs.keys()) == set(get_function_args(self.simulation_function)), "Kwargs don't match function argument list"  
        return kwargs
        
        
    def empty_simulation_protocol(self, simulation_type):
        """ 
        Return dict with empty values to show required structure of simulation protocol for
        a given simulation type for user filling out.
        Therefore, we don't include parameters about mechanisms or sim_ids that are filled out internally.
        """
        simulation_type = simulation_type.lower()
        names = self.get_simulation_protocol_names()
        kwargs = {}

        if simulation_type in ['iclamp', 'vclamp']:
            #kwargs['mechanisms'] = None # These are specified by the population
            #kwargs['mechanism_names'] = None
            for kwarg in [name for name in names['clamp'] if name not in ['mechanisms', 'mechanism_names', 'biomarker_names']]:
                kwargs[kwarg] = None
            
            #assert set(kwargs.keys()) == set(names['clamp']), "Parameters common to iclamp+vclamp missing."
            
        # Now do specific parameters
        for kwarg in names[simulation_type]:
            kwargs[kwarg] = None
            
        # finally, assign the "always" names
        for kwarg in names['always']:
            if kwarg not in ['sim_id', 'metadata', 'options']: # sim_id and metadta are an internal variable and data structure, while options is built from the function call inputs, so none are directly user-specified
                kwargs[kwarg] = None
        return kwargs  
    
    
    def get_simulation_protocol_defaults(self):
        """ 
        Get default values for simulation protocols.
        This could be useful if some are missing in the input? 
        Or do we want to enforce having all the inputs?
        I think we want to force all inputs for research purposes - makes every 
        simulation parameter explicit.
        """
        assert False, "Not implemented to encourage simulation inputs to be explicitly documented in run scripts."
            
        
    def get_biomarker_names(self, sim_type, outputs):
        """ 
        Generate the list of biomarker names we'll calculate.
        Currently a bit messy as to add a new biomarker we need to add it here
        and in the simulation code. But it should cause an error if I miss one or the other.
        
        Returns
        -----------
        List of biomarker names for the simulation type given, for the output codes specified.
        """   
        biomarker_names = {}
        
        # Common iclamp and vclamp names
        biomarker_names['iclamp'] = ['APFullWidth', 'APPeak', 'APRiseTime', 'APSlopeMin', 'APSlopeMax',
                'AHPAmp', 'AHPTau', 'AHPTrough', 'ISI', 'RMP', 'Rheobase', 'Firing pattern']
        biomarker_names['vclamp'] = []
        
        # Set biomarkers that are only calculated when requested through outputs
        if 'ik' in outputs:
            biomarker_names['iclamp'].extend([]) # Nothing yet
            biomarker_names['vclamp'].extend(['ik_steps_absmax'])
        if 'ina' in outputs:
            biomarker_names['iclamp'].extend([]) # Nothing yet
            biomarker_names['vclamp'].extend([]) # Nothing yet
            
        return biomarker_names[sim_type.lower()]


    def process_flags(self, flags, sim_id):
        """
        Processes any flags given to the simulation and returns the correct information.
        """
        processed_flags = {}
        for flag, val in flags.items():
            assert flag in self.allowed_flags(), "Flag {} not in allowed_flags".format(flag)
            if flag == "ramp_threshold_sim_for_width":
                sim_name = val
                # Checks
                assert self.population != None, "Need a population"
                assert sim_name in self.population.results.columns, "sim_name {} not found".format(sim_name)
                assert sim_id in self.population.results.index, "sim_id {} not found".format(sim_id)

                threshold = self.population.results.at[sim_id, (sim_name, 'Threshold')]
                processed_flags[flag] = threshold

            elif flag == "rheobase_sim_for_stim_amp":
                sim_name = val
                # Check self.population != None, "Need a population"
                assert sim_name in self.population.results.columns, "sim_name {} not found".format(sim_name)
                assert sim_id in self.population.results.index, "sim_id {} not found".format(sim_id)

                amp = self.population.results.at[sim_id, (sim_name, 'Rheobase')]
                processed_flags[flag] = amp
                
            else:
                raise ValueError("Flag {} not supported.".format(flag))

        return processed_flags

    def allowed_flags(self):
        """
        Returns a list of allowed strings to be passed to simulations as flags
        """
        allowed_flags = []
        # Provide the name of a previously run simulation to use that simulation for threshold
        # in calculation of AP full width following Davidson et al. PAIN 2014. 
        allowed_flags.append("ramp_threshold_sim_for_width")
        # Use rheobase from a previous simulation as the stimulus amplitude
        allowed_flags.append("rheobase_sim_for_stim_amp")
        return allowed_flags
        
    " Output functions "

    def log_test(self, result):
        for biomarker in self.results.columns:
            self.results.loc[0, biomarker] = result
    

    def log_result(self, result):
        """
        Stores the result from a simulation using callback functionality. 
        """
        keys = result.keys()
        biomarker_names = [key for key in keys if key not in ['sim_id', 'trace']]
        sim_id = result['sim_id']
        for biomarker in biomarker_names:
            # @TODO - do this when we create the results df to avoid doing it all the time!
            #if isinstance(result[biomarker], collections.Sequence):
            #self.results[biomarker] = self.results[biomarker].astype(object)

               self.results.at[sim_id, biomarker] = result[biomarker] 
               # Use at not loc because loc breaks with lists (like firing patterns)
               # at can only access one value whereas loc could access multiple values
               # in this case we are always only wanting to put a result into a single cell
               # see: https://stackoverflow.com/questions/26483254/python-pandas-insert-list-into-a-cell
            
        # If we are doing figure plotting of multiple traces
        # store the trace in temporary storage
        if (self.options['save'] == True) & ('fig' in process_save_type(self.options['save_type'])):
                self.traces[sim_id] = result['trace']
                
        # Check for plotting
        filename = filename =  '{0}_{1}_traces_{2}.png'.format(self.population.name, self.name, self.count)
        plotted = self.trace_plot(self.num_subplots, filename)
        if plotted: 
                self.count += 1 
                
    
    def trace_plot(self, traces_to_plot, filename, force_plot=False, require_contiguous=True):
        """  
        If we want to save subplots of traces, decide if we have enough saved traces to plot
        First implementation was to just look for contiguous blocks - but what if we plot 1->101 
        and miss trace 0 in the first?
        Second implementation (CURRENTLY USED): demand we start with the first sim_id, then when we've done a plot
        move the pointer to sim_id+traces_to_plot+1.
        """
        # Only do anything if the right save options are set
        if (self.options['save'] == True) & ('fig' in process_save_type(self.options['save_type'])):
            subplot_dim = int(np.ceil(np.sqrt(traces_to_plot))) # Number of subplots (square)

            # Check whether there are enough traces to make a plot
            # as we work through the stored plots
            make_plot = False
            
            # If we're forcing a plot at the end of the simulation plot regardless as long as there are traces
            num_traces = len(self.traces)
            if (force_plot == True) & (num_traces > 0):
                make_plot = True
                
            # Else, check we have enough traces stored
            # If so, if we require them to be contiguous check from the first trace (by sim_id)
            if (num_traces >= traces_to_plot):
                if require_contiguous == True:
                    ids = self.traces.keys()
                    ids = sorted(ids)
                    required_format = self.model_indices[self.plotting_pointer:self.plotting_pointer+traces_to_plot]
                    # Change this if statement if we want to be able to plot blocks beyond the first
                    if all([i == j for i,j in zip(ids[0:traces_to_plot], required_format)]):  
                        make_plot = True
                        self.plotting_pointer += traces_to_plot
                else:
                    make_plot = True
            
            if make_plot == True:
                # Get the first traces_to_plot worth of traces, or if there's not enough, all of them
                # a_list[0:big number] will return to the end of the list and stop
                ids = sorted(self.traces.keys())
                ids = ids[0:traces_to_plot]
                
                outputs = self.protocols['outputs']
                # @TODO - use outputs to decide what else to return
                # For now we plot the first column as x, and each other column as y
                # We could include the outputs in the title or a legend in the first plot
                
                # Plotting
                fig = plt.figure(figsize=(40,40))
                legend_plotted = False;
                for i, sim_id in enumerate(ids):
                
                    # Trace plots
                    plt.subplot(subplot_dim,subplot_dim,i+1)
                    trace = self.traces[sim_id]
                    # Only plot if trace is dictionary or np array (trace may be None if no valid trace produced)
                    if isinstance(trace, np.ndarray):
                        for j in range(1,trace.shape[1]):
                            plt.plot(trace[:,0], trace[:,j]) # Col 0 is time                          
                    elif isinstance(trace, dict):
                        for key in [key for key in trace.keys() if key != 't']: # DRY
                            _data = trace[key]
                            # If data is a dictionary assume it's an ionic current 
                            # and try to reshape each entry in dictionary into an np.ndarray (from a hoc vector)
                            # Plot if the result has the same length as time
                            if isinstance(_data, dict):
                                for _current_type, _current in _data.items():
                                    _ndarray_current = np.array(_current)
                                    if len(trace['t']) == len(_ndarray_current):
                                        plt.plot(trace['t'],_ndarray_current)
                            # If it's not a dictionary, plot it if it's the same length as time
                            else:
                                if len(trace['t']) == len(_data):
                                    plt.plot(trace['t'], _data)

                    # Legend
                    if legend_plotted == False:
                        if isinstance(trace,dict):
                            legend = [key for key in trace.keys() if key != 't']
                        else:
                            legend = ['Vm'] # @TODO - use outputs to append to here
                        plt.legend(legend)
                        legend_plotted=True     
                    
                    # Annotate with model id and firing pattern
                    firing_pattern = self.results.at[sim_id, 'Firing pattern']
                    plt.title('{}: {}'.format(sim_id, firing_pattern))
                    #plt.title(sim_id) 
                
                # Save figure
                fig.savefig(filename, dpi=150)
                plt.close(fig)
                fig = None
                print("Saved to {}, number of traces stored = {}".format(filename, len(self.traces)))
                sys.stdout.flush()
                # Clear traces that were plotted after plotting
                for id in ids:
                    del self.traces[id]
                return True
                
            elif make_plot == False:
                return False
            else: 
                raise ValueError("make_plot is invalid value: {}".format(make_plot))
                
       
# --- Utils ---

def set_dataframe_types(df):
    """
    Set all columns in the DataFrame where all data has the same type to that type
    """
    typed_df = df.copy()
    for column in typed_df:
        types = [type(i) for i in typed_df[column]]
        if all([i == types[0] for i in types]):
            typed_df[column] = typed_df[column].astype(types[0])
    return typed_df

def dump_exception(e, filename="debug.txt"):
    with open(filename, "w") as f:
        f.write(str(e))
