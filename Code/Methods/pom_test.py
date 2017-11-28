# -*- coding: utf-8 -*-
"""
Created on Mon Feb 08 12:14:25 2016
Functions for running a population of models loop

@author: Oliver Britton
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import datetime
import time
#import pdb
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
    if type(trace) == dict:
        data = np.column_stack(trace['t'], trace['v'])
    elif type(trace) == np.array:
        data = trace # assume formatting is correct
    elif (type(trace) == tuple) | (type(trace) == list): # assume list in order: t,v
        data = np.column_stack([trace[0], trace[1]])
    else:
        raise TypeError("Can't parse type of trace")
    np.savetxt(filename, data)
        
def plot_trace(filename):
    data = read_trace(filename)
    plt.plot(data['t'],data['v'])
             
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

# ---CLASSES----
    
class PopulationOfModels(object):

    def __init__(self, 
        name,
        model_details,
        sim_protocols=None,
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
        self.simulation_protocols = sim_protocols 
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
    
    def setup_simulation(self, name, type, protocols, options=None):
        """
        Initialise a set of simulations
        
        Parameters
        ----------------
        name: str, default 'sim'
        type: str, default 'IClamp', options: 'IClamp', 'VClamp', 'Test' (case insensitive)
        protocols: dict, default None, contents dependent on 'type'.
        
        Returns
        -----------
        Nothing, creates self.simulations[name].
        
        """
        if name in self.simulations.keys():
            raise ValueError('Simulation name {} is already present.'.format(name))
        else:
            # Assemble all the details about the simulation to give to the object
            protocols['simulation_type'] = type
            self.simulations[name] = Simulation(name, protocols, options=options, population=self)
    
    def run_simulation(self, 
    name='sim', 
    simulation_type='IClamp',
    simulation_protocols=None, 
    sampling_freq_hz=20000, 
    plot=False, save=False, benchmark=True):
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
        
        """
        # TODO: Add parallelisation here - see parallelisation in https://github.com/mirams/PyHillFit/blob/master/python/PyHillFit.py
        
        # Set up the simulation
        if simulation_protocols == None:
            simulation_protocols = self.simulation_protocols
        assert(simulation_protocols != None), "Simulation protocol not set."
        self.setup_simulation(name, type, protocols) 
        sim = self.simulations[name]
        parameters = self.results['Parameters']
        
        
        
        
        # @START TRANSFERRING TO Simulation.pom_simulation() here
        start_time = time.time()
        num_sims = len(self.results.index)
        if plot:
            plt.figure(figsize=(15,15))
        if save:
            self.traces[simulation_name] = {}
        # @done up to here
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
           
            """
            Simulation code transferred to simulate_iclamp in Simulation
            """
           
           
        
        # At end of all simulations, build a multiarray with the simulation name and the biomarker results.
        self.biomarker_results = pd.DataFrame(columns=pd.MultiIndex.from_product([[simulation_name], biomarker_results.columns]), index=self.results.index)
        # Insert results
        self.biomarker_results.loc[:] = biomarker_results
        # Then concatenate it into the main results dataframe.
        self.results = pd.concat([self.results, self.biomarker_results], axis=1)
        
        print("Simulation {} complete in {} s.".format(simulation_name,time.time()-start_time))
        
        #@FINISH TRANSFERRING to Simulation.pom_simulation here
            
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
    
    """
    def __init__(self, name, protocols, options=None, population=None):
    
        self.name = name
        self.protocols = protocols
        self.traces = {} # To store traces for plotting or analysis
        
        # Population-specific setup
        if population != None:
            self.population = population # Set reference to population
            self.results = pd.DataFrame(index=self.population.results.index) 
            self.simulation_ran = False # Check if simulation has been run - only run once
        else:
            # @TODO: Add in support for single cell simulations
            raise ValueError('Population must be supplied to run simulation, for single cell simulations use simulation_helpers module.')
            
        # Options - output, saving, logging etc.
        if options != None:
            pass
    
    def test_model(self, gna, gk):
        cell = h.Section()
        cell.insert('hh')
        cell.gnabar_hh *= gna
        cell.gkbar_hh *= gk
        return cell
        
    def simulate_test(self, conductances):
        """
        Test parallel simulation
        """
        
        cell = self.test_model(conductances[0],conductances[1])
        t = h.Vector()
        v = h.Vector()
        t.record(h._ref_t)
        v.record(cell(0.5)._ref_v, sec=cell)

        stim = h.IClamp(0.5, sec=cell)
        stim.dur = 500
        stim.amp = 5.0
        stim.delay = 500
        tstop = 1500.

        h.finitialize(-65.) # must go after record
        neuron.run(tstop)
        output = (np.array(t),np.array(v))
        
        return "{},{}".format(mechanisms[0], mechanisms[1]), output
          
    def simulate_iclamp(self, 
                                        sim_id,
                                        mechanisms, 
                                        mechanism_names,
                                        amp,
                                        dur,
                                        delay,
                                        interval,
                                        num_stims,
                                        stim_func,
                                        outputs):
        """
        Runs a full simulation in IClamp mode (finds rheobase, finds rmp, and calculates biomarkers at rheobase)
        
        Parameters
        ----------------
        
        Returns
        -----------
        """
    
        # Setup conditions
        
        # Use sh.simulation to run each simulation
        """
         self.active_parameters = parameters.loc[model_idx]
            mechanisms = {self.parameter_designations[param]:val for param, val in self.active_parameters.iteritems()} # iteritems over pd.Series to get key value pairs and generate a dict from them
            # To do - allow different cellular parameters and ionic concs to be input here or in bespoke_simulation, as well as variation in parameters
            
            # Build model using dict with full parameter names (e.g. gbar_nav17vw not nav17vw)
            self.active_model = sh.build_model(mechanisms=mechanisms, mechanism_names=self.mechanism_names, conductances=None, mechanism_is_full_parameter_name=True)

            " @To do - wrap all the simulations into a member function "
            # Set temperature
            h.celsius = self.celsius # In line with Davidson et al., 2014, PAIN
            # Set ionic conditions (update reversal potentials in line with ionic conditions)
            " @To do - check which ions exist and only set conditions for them "
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
            """
        
    def simulate_vclamp(self, 
                                        sim_id,
                                        mechanisms,
                                        mechanism_names,
                                        ions,
                                        hold,
                                        steps,
                                        delta,
                                        durations,
                                        t_stop,
                                        outputs):
        """
        Simulation of standard voltage clamp protocol.
        Example implementing following protocol for K+ recordings from Anabios:
        simulate_vclamp(mechanisms=@TODO, mechanism_names=@TODO, hold=-70.0, steps=[-80,80], delta=10.0, durations=[1100,500,500], outputs=['v', 'ik'])

        1. Hold at -70 mV for 1100 ms
        2. Pulse from -80 to 80 mV with delta=10mV
        3. Hold at -70 mV for 500 ms
        
        Parameters
        ----------------
        sim_id - int/float/str - identifier for the simulation 
        model_data - dict - data to construct model from
        hold - holding potential - float (mV)
        delta - size of steps - float (mV) or None to specify steps manually
        steps - start and end of voltage clamp if delta is specified - 2 element list of floats (mV)
        durations - start of each stage of voltage clamp - list of floats/ints 
        outputs - list of strings, options to include are  'v' (voltage), 'ik' (total k+ current)....
        @TODO - include more options.
        
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
        
        # Setup voltage clamp steps
        if delta == None:
            assert len(steps) > 0 # Use steps as given 
        else:
            steps = np.arange(steps[0], steps[1]+(0.1*delta),delta)
        

        # Initialise results
        results = {}
        results['sim_id'] = sim_id
        if 'ik' in outputs:
            results['ik_steps_absmax'] = np.zeros(len(steps))
        
        # Simulation loop through steps
        for i, step in enumerate(steps):
            # Build model using dict with full parameter names (e.g. gbar_nav17vw not nav17vw)
            cell = sh.build_model(mechanisms=mechanisms, 
                                                mechanism_names=mechanism_names, 
                                                conductances=None, 
                                                mechanism_is_full_parameter_name=True)
            
            # Setup this simulation
            self.setup_simulation_conditions()
            
            # @TODO - put setting ionic conditions into a function to share with IClamp
            if 'K' in ions:
                oldstyle = h.ion_style("k_ion", 1, 2, 1, 1, 0,sec=cell)
            if 'Na' in ions:
                oldstyle = h.ion_style("na_ion", 1, 2, 1, 1, 0,sec=cell)
            if 'Ca' in ions:
                oldstyle = h.ion_style("ca_ion", 1, 2, 1, 1, 0,sec=cell)

            v,t = sh.set_vt(cell=cell)
            vectors = self.setup_output_recording(cell, outputs)
            
            clamp = h.VClamp(0.5, sec=cell)
            clamp.amp[0] = hold
            clamp.dur[0] = durations[0]

            clamp.amp[1] = step
            clamp.dur[1] = durations[1]

            clamp.amp[2] = hold
            clamp.dur[2] = durations[2]
            
            if len(durations) > 3:
                t_stop = durations[3]
            else: 
                t_stop  = durations[0] + durations[1] + durations[2]

            h.finitialize(hold) # Vital! And has to go after record
            # Run simulation
            neuron.run(t_stop)
             
            # ---Compute and store outputs---
            
            # @TODO - process save and plot options
            # if self.save: ...
            # if self.plot: ...
            
            
            # K+ biomarkers
            # @TODO - for now just return a vector of the peak current from each step
            if 'ik' in outputs: 
                results['ik_steps_absmax'][i] = nb.absmax(vectors['ik'])
                # @TODO Add in ih (na component)
                #ik = np.array(ik) + np.array(ih_na)

        return results
            
    def pom_simulation(self, cores=1, plot=False, save=False, save_type="fig", benchmark=True):
        """
        Run  simulations 
        
        Parameters
        -----------------
        parameters: pd.DataFrame
        simulation_parameters: dict
        cores: int, default 1
        plot: bool, default False
        save: bool, default False
        save_type: str, default 'fig', options are...
        benchmark: bool, default True
        
        Returns
        -----------
        
        """
        if self.simulation_ran == True:
            print("Error: this simulation has already been ran, cannot run same simulation again with the same simulator with name: {}.".format(self.name))
            return
            
        if self.population == None:
            raise ValueError('No population specified - cannot run population of models simulation!')
        
        # Startup options
        start = time.time()
        output_frequency = 100. # Number of outputs to make - 100 equals every 1%
        if plot == True:
            pass
            # @TODO
            #plt.figure(figsize=(15,15))
        # @ TODO
        if save == True:
            if save_type in ['fig', 'trace', 'both', 'all']:
                pass # we're ok
            else:
                raise ValueError('save_type not found.')
        
        # Simulation setup:
        simulation_type = self.protocols['simulation_type']
        self.set_simulation_function(simulation_type)
        num_sims = len(self.population.results.index)
        
        # ---- Main parallel simulation loop ----
        if __name__ == '__main__':
            pool = mp.Pool(processes=cores)
            model_indices = self.population.results.index
            assert(len(model_indices) == len(set(model_indices))), "Duplicates in model indices." 
            for i, model_idx in enumerate(model_indices):   
                # @TODO - this probs has parallel bugs
                if benchmark:
                    now = time.time()
                    reporting_period = np.ceil(len(model_indices)/output_frequency)
                    if i > 0:
                        if (i == len(self.results.index)-1) | (i%reporting_period == 0):
                            clear_output()
                            print("Simulation {} of {} started. Time taken = {:.1f}s. Estimated remaining time = {:.1f}s.".format(i+1, num_sims, now-start_time, (num_sims-i)*(now-start_time)/i))
                    else:
                        print("Simulation set of {} simulations begun.".format(num_sims))
                        
            # --- Setup model ----
            active_parameters = parameters.loc[model_idx]
            # iteritems over pd.Series to get key value pairs and generate a dict from them
            mechanisms = {self.population.parameter_designations[param]:val for param, val in self.active_parameters.iteritems()} 
            
            # To do - allow different cellular parameters and ionic concs to be input here or in bespoke_simulation, as well as variation in parameters
            
            sim_kwargs= self.build_simulation_kwargs(sim_id=model_idx, 
                                                                                      simulation_type=simulation_type, 
                                                                                      simulation_parameters=self.protocols, 
                                                                                      mechanisms=mechanisms)
            
            # Run simulation function in parallel
            pool.apply_async(simulation_function, kwds=sim_kwargs, callback=self.log_result)
            pool.close()
            pool.join()
            
        
        # Clean up and lock simulation
        self.simulation_ran = True
            
            
    def set_simulation_function(self, simulation_type):   
        """
        Function to set simulation 
        """
        simulation_type = simulation_type.lower()
        if simulation_type == 'iclamp':
            self.simulation_function = self.simulate_iclamp
        elif simulation_type == 'vclamp':
            self.simulation_function = self.simulate_vclamp
        elif simulation_type == 'test':
            self.simulation_function = self.simulate_test
        else: 
            raise ValueError('simulation type: {} is not found.'.format(simulation_type))

            
    def build_simulation_kwargs(self, sim_id, simulation_type, simulation_parameters, mechanisms):
        """
        Setup dict of keyword arguments to feed to simulation function
        """
        def assign_kwargs(simulation_type):
            kwargs = {}
            # First do parameters common to IClamp and VClamp and check we get them all
            if simulation_type in ['iclamp', 'vclamp']:
                kwargs['mechanisms'] = mechanisms
                kwargs['mechanism_names'] = self.population.mechanism_names
                kwargs['ions'] = simulation_parameters['ions']
                assert set(kwargs.keys()) == set(names['clamp']), "Parameters common to iclamp+vclamp missing."
                
            # Now do specific parameters
            for kwarg in names[simulation_type]:
                kwargs[kwarg] = simulations_parameters[kwarg]
                
            # finally, always assign sim_id for every simulation_type
            kwargs['sim_id'] = sim_id
            return kwargs
        
        simulation_type = simulation_type.lower() # convert to lower case
        assert (simulation_type in ['iclamp', 'vclamp', 'test']), 'Simulation type: {} not found'.format(simulation_type)
        
        names = self.get_simulation_parameter_names()
        kwargs = assign_kwargs(simulation_type)
        return kwargs
            
            
    def setup_simulation_protocols(self, sim_protocols):
        """
        OLD CODE NOT NEEDED PROBABLY?
        
        Parameters
        -----------------
        sim_protocols: dict containing fields describing all aspects of simulations protocols
        see code for necessary and optional fields. 
        @TODO: make more protocol fields optional with defaults (e.g. v_init at -65)
        """
        print('OLD CODE NOT NEEDED PROBABLY?')
        if 'biomarker_protocol' in sim_protocols.keys():
            self.biomarker_protocol = sim_protocols['biomarker_protocol']
        else:
            self.biomarker_protocol = None
            
        if 'ions' in  sim_protocols.keys():
            self.ions = sim_protocols['ions']
        
        self.delay = sim_protocols['delay']
        self.amp = sim_protocols['amp']
        self.dur = sim_protocols['dur']
        self.interval = sim_protocols['interval']
        self.num_stims = sim_protocols['num_stims']
        self.stim_func = sim_protocols['stim_func']
        self.t_stop = sim_protocols['t_stop']
        self.v_init = sim_protocols['v_init']
        self.currents_to_record = sim_protocols['currents_to_record']
        
    
    def setup_simulation_conditions(self):
        """
        Do anything that needs to be done 
        """
        h.celsius = self.population.celsius
        
    def setup_output_recording(self, cell, outputs):
        """
        Create vectors for recording outputs (except voltage and time, which are givens)
        """
        vectors = {}
        if 'ik' in outputs:
        # @TODO - write more outputs and wrap into function
            vectors['ik'] = h.Vector()
            vectors['ik'].record(cell(0.5)._ref_ik, sec=cell)
        if 'ina' in outputs:
            vectors['ina'] = h.Vector()
            vectors['ina'].record(cell(0.5)._ref_ina, sec=cell)
        if 'ik_ih' in outputs:
            assert(False), 'We don\'t support the hcn current as part of ik yet.'
            
        # Could use something like this to automate:
        """ Use exec to record the right current variable (e.g. ik_kdrtf)
        For output in outputs:
        exec('vectors[output].record(self.active_model(0.5)._ref_{}, sec=cell'.format(output)) """
        return vectors
        
    
    def log_result(self, result):
        """
        Stores the result from a simulation using callback functionality. 
        @OPTIMIZE - if we had loads of biomarkers it would make sense to vectorize the for loop?
        """
        keys = result.keys
        biomarker_names = [key for key in keys if key != 'sim_id'] # could change to not in ['sim_id', 'blah'] if we want to add more keys to remove
        sim_id = results['sim_id']
        for biomarker in biomarker_names:
            self.results.loc[sim_id, biomarker] = result[biomarker]
        
        
    def get_simulation_parameter_names(self):
        """
        Gets list of parameter names for each simulation 
        """
        names = {}
        names['clamp'] = ['mechanisms', 'mechanism_names', 'ions'] # Common to IClamp and VClamp
        names['iclamp'] = ['amp', 'dur', 'delay', 'interval', 'num_stims', 'stim_func']
        names['vclamp'] = ['hold', 'steps', 'delta', 'durations']
        names['test'] = ['conductances']
        return names
