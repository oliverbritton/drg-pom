# -*- coding: utf-8 -*-
"""
Simulation helpers 

Created on Thu Jun 29 11:17:20 2017

@author: Oliver Britton
"""
import os
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import numbers
import pyDOE
import itertools

# Load Neuron and DRG code
sys.path.append('E:\\CLPC48\\Neuron Project')
sys.path.append('E:\\CLPC48\\Neuron Project\Code')
import NeuronProjectStart

from neuron import h
import neuron

# We seem to need to load now to allow parallel scripts to work
# Perhaps there is a way to reload if we need to?
import Methods.Simulations.loadneuron as ln
ln.load_neuron_mechanisms(verbose=True)

def init_model(mechanisms=[],
    L=30.,  # uM
    diam=46., # uM
    cm=1., # uF/cm^2
    Ra=123., 
    ko=3., # mM
    ki=135.,  # mM
    nao=145., # mM
    nai=5., # mM
    ):
    """
    Initialise a cylindrical cell model, with the specified ionic mechanisms embedded.
    """
    cell = h.Section()
    for mechanism in mechanisms:
        try:
            cell.insert(mechanism)
        except ValueError:
            raise ValueError("argument: {} is not a density mechanism name.".format(mechanism))
    cell.L = L # microns
    cell.diam = diam #microns
    cell.cm = cm  # uF/cm2
    cell.Ra = Ra # ohm-cm
    
    try:
        cell.ko = ko
        cell.ki = ki
    except NameError:
        cell.insert('k_dummy')
        cell.ko = ko
        cell.ki = ki
        
    try:
        cell.nao = nao
        cell.nai = nai
    except NameError:
        cell.insert('na_dummy')
        cell.nao = nao
        cell.nai = nai
    
    return cell
    
def get_init_model_values(name='first'):
    " Returns model geometry and ionic concentration data to save as metadata "
    allowed_names = ['first']
    model_values = {}
    if name not in allowed_names:
        raise ValueError('Name of init model value set not found, allowed names are: {}.'.format(', '.join(allowed_names)))
    elif name == 'first': # Uses Choi/Waxman geometry, standard cm, and Davidson ionic concentrations.
        model_values['L'] = 30. # uM
        model_values['diam'] = 46. # uM
        model_values['cm'] = 1. # I uF/cm^2
        model_values['Ra'] = 123.
        model_values['ko'] = 3. # mM
        model_values['ki'] = 135. # mM
        model_values['nao'] = 145. # mM
        model_values['nai'] = 5. # mM
        
    return model_values
        
# TODO - make the build model inputs more sensibly named, they are confusing
def build_model(mechanisms={'kdrtf':1., 'katf':1., 'nav18hw':1.}, conductances=None, mechanism_names=None,mechanism_is_full_parameter_name=False):
    """
    Mechanism names is for if we don't want to vary the conductance of every mechanism in the model, or if we want to use different parameters. 
    """
    # Dict construction
    if type(mechanisms) == dict:
        if conductances:
            raise ValueError('conductances should not be provided with a dict')
        if mechanism_names == None: # Otherwise use provided names
            mechanism_names = mechanisms.keys()
            if mechanism_is_full_parameter_name:
                # Split on underscores and throw away the first bit as that will be parameter type (e.g. gbar)
                mechanism_names = [name.split('_',1)[-1] for name in mechanism_names]

        model = init_model(mechanisms=mechanism_names) 
        for mechanism, conductance in mechanisms.items():
            if mechanism_is_full_parameter_name: # Mechanisms contains the full parameter name, not just the mechanism name
                param_setting_str = 'model.{0} *= {1}'.format(mechanism, conductance)
            else: # Mechanism is assumed to be a suffix for a conductance
                param_setting_str ='model.gbar_{0} *= {1}'.format(mechanism, conductance)
            exec(param_setting_str)
        if conductances:
            assert False
            
    # List construction
    elif type(mechanisms) == list:
        if mechanism_is_full_parameter_name:
            mechanism_names = [name.split('_')[-1] for name in mechanisms]
            model = init_model(mechanisms=mechanism_names)
        else:
            model = init_model(mechanisms=mechanisms)
        if conductances:
            for mechanism, conductance in zip(mechanisms,conductances):
                if mechanism_is_full_parameter_name:
                    exec('model.{0} *= {1}'.format(mechanism, conductance)) 
                else:
                    exec('model.gbar_{0} *= {1}'.format(mechanism, conductance)) 
    else:
        assert False
    return model
    
def set_stims(amp, dur, delay, interval, num_stims, stim_func, cell):
    """
    Setup stimulus objects for a simulation
    Currently does not work for more than one stimulus
    """
    stims = []
    for i in range(num_stims):
        stims.append(stim_func(0.5, sec=cell))
        stims[-1].dur = dur
        stims[-1].amp = amp # Same total current as for ramp
        stims[-1].delay = delay + i*(dur + interval)
        
        return stims
        
def simulation_plot(t, v, currents=None, plot_type='default'):
        
    if plot_type == 'default':
        plt.figure(figsize=(5,5))
        plt.plot(t,v); plt.ylabel('Vm (mV)'); plt.xlabel('Time (ms)')
    # TODO: Support plotting currents
    else:
        raise ValueError("Plot type: {} not supported.".format(plot_type))
        
def set_vt(cell):
    v = h.Vector()
    v.record(cell(0.5)._ref_v, sec=cell)
    t = h.Vector()
    t.record(h._ref_t)
    return v,t
    
def record_currents(cell, current_set):
    """
    Record a list of currents from the cell model 
    
    Inputs:
    cell - NEURON cell model
    current_set - list of strings
    
    Outputs:
    Output is stored in a dict of dicts like so:
    currents = {mechanisms + ionic current names}
    currents[name] = {ionic current(s)}
    currents[name][ionic current] = h.Vector
    
    So for example, if we had a cell with mechanisms: kdr, nav18hw, and hcntf, our dict would be (V for recorded vector):
    currents = {'ina':{'ina':V}, 'ik':{'ik':V}, 'kdr':{'k':V}, 'nav18hw':{'na':V}, 'hcntf':{'na':V, 'k':V} }. Note for hcntf there are two
    recorded vectors.
    So to plot every output you would go
    for current in currents:
        for i in currents[current]:
            plt.plot(t,currents[current][i]
            plt.legend.append('{}_{}'.format(i,current)) # gives things like 'ina_nav18hw' or 'ina_ina'
    """
    
    if current_set == None:
        currents = None
    else:
        currents = {}
        
        for current in current_set:
            currents[current] = {}
            # Process whole ionic currents (e.g. Ina)
            if current in get_ionic_current_list():
                current_refs = {current:'_ref_{}'.format(current)}
            # Process ion channel mechanisms (e.g. nav 1.7)
            else:
                if (current[0:2] == 'na') | ('nav' in current.lower()):
                    i_strs = ['ina']
                elif (current[0] == 'k') | (any(s in current.lower() for s in ['kdr', 'ka', 'km', 'kv'])):
                    i_strs = ['ik'] 
                elif current[0:2] == 'ca':
                    i_strs = ['ica']
                elif 'hcn' in current.lower():
                    # Special case for hcn currents as they may have two different currents
                    # Or a non-specific current
                    if current == 'hcntf':
                        i_strs = ['ina', 'ik']
                    elif current == 'ch_HCNp':
                        i_strs = ['i']
                    elif current == 'hcnkn':
                        raise ValueError('hcnkn not supported yet') # TODO support hcnkn recording
                    else:
                        raise ValueError('hcn current: {} not supported'.format(current))
                elif 'pas' in current:
                    i_strs = ['i']
                elif current[0:5] == 'iconc':
                    continue # Ignore concentration mechanisms as they don't have a current
                else:
                    raise ValueError('current: "{}" is not supported for current recording.'.format(current))
                # Finish by building current reference(s) to record from
                current_refs = {i_str:'_ref_{}_{}'.format(i_str, current) for i_str in i_strs}
            
            # Finally, set up recording for each 
            for name, current_ref in current_refs.items():
                currents[current][name] = h.Vector()
                currents[current][name].record(getattr(cell(0.5), current_ref), sec=cell)
                
    return currents
    
def record_concs(cell, concs_to_record):
    """
    Record ionic concentrations from a model
    """
    if concs_to_record == None:
        concs = None
    else:
        concs = {}
        for conc in concs_to_record:
            assert conc in get_ionic_concentration_list(), "Concentration name '{}' not in list of allowed concentrations".format(conc)
            conc_ref = '_ref_{}'.format(conc)
            concs[conc] = h.Vector()
            concs[conc].record(getattr(cell(0.5), conc_ref), sec=cell)
    return concs
    
def generate_current_set(mechanisms_to_record, record_ionic_currents=True):
    """
    Takes a list of mechanisms and use this to generate a list of currents to record
    Also, if you give it an ionic current e.g. ica, ina, ik, it will record that current.
    """
    assert type(mechanisms_to_record) == list, 'mechanisms must be in a list'
    currents_to_record = mechanisms_to_record.copy() # avoid mutating the original
    if record_ionic_currents:
        if any([i[0:2] == 'na' for i in currents_to_record]):
            currents_to_record.append('ina')
        if any([i[0] == 'k' for i in currents_to_record]):
            currents_to_record.append('ik')
        if any([i[0:2] == 'ca' for i in currents_to_record]):
            currents_to_record.append('ica')
    return currents_to_record
        
"I think this can be deleted 26/03/18 "
"""
def build_sim_protocols(amp, dur, delay, interval, num_stims=40, stim_func=h.IClamp, t_stop=1000., v_init=-65.,):
    # Constructs a dictionary of stimulation protocols for simulate_ functions
    sim_protocols = {'amp':amp, 'dur':dur, 'delay':delay, 'interval':interval, 'num_stims':num_stims, 'stim_func':stim_func, 't_stop':t_stop, 'v_init':v_init,}
    return sim_protocols
"""
    
" -- Functions for building parameter sets -- "
    
def build_parameter_set_details(num_models, num_parameters, min, max, output_filename, parameter_names=None):
    " Parameter set details ... "
    param_set_details = {}
    param_set_details['num_models'] = num_models
    param_set_details['num_parameters'] = num_parameters
    param_set_details['min'] = min
    param_set_details['max'] = max
    param_set_details['output_filename'] = output_filename
    param_set_details['parameter_names'] = parameter_names
    
    return param_set_details

def build_parameter_set(num_models, parameter_data, minimum=None, maximum=None, filename=None,  save=False):
    """
    Construct parameter sets with names for each parameter using latin hypercube sampling and optionally save them as csvs.
    Inputs:
    num_models - int, number of models to construct parameter sets for.
    parameter_data - dict or list, either a list of parameter names (uniform scaling factors over all parameters) or a dict of form {'name1':(min1, max1), 'name2':(min2, max2)...} if variable scaling factors are required.
    minimum, maximum - ints, scaling ranges if scaling ranges are to be the same over all parameters
    output_filename - string, filename to output parameter sets to if required.
    save - bool, whether sets should be saved to output_filename
    Outputs:
    parameter_set - a dataframe of labelled parameter scaling factors with rows for each model and labelled columns for each parameter.
    """
    
    parameter_names = [name for name in parameter_data] # works for lists and dicts
    num_parameters = len(parameter_names)
    # Build baseline LHS
    parameter_array = pyDOE.lhs(num_parameters, samples=num_models)
    parameter_sets = pd.DataFrame(parameter_array,columns=parameter_names)
    
    " Transform parameter set using minimum(s) and maximum(s) "
    if isinstance(parameter_data, dict):
        # If dict provided, iterate over each parameter and scale it separately
        header = 'nonuniform scaling factors - '
        for parameter_name, range in parameter_data.items():
            # Range[0] is minimum scaling factor, range[1] is maximum
            assert range[0] < range[1], "Min is not less than max."
            parameter_sets[parameter_name] *= (range[1] - range[0])
            parameter_sets[parameter_name] += range[0]
            header += '{}: ({} {}) '.format(parameter_name, range[0], range[1]) # Don't use commas as we're saving as csv
    elif isinstance(minimum, numbers.Number) & isinstance(maximum, numbers.Number) & isinstance(parameter_data, list):
        # If min and max provided - scale all parameters by them and use parameter_data as list of names
        assert minimum < maximum, "Min is not less than max."
        parameter_sets *= (maximum-minimum)
        parameter_sets += minimum
        header = 'uniform scaling factors - '
        for parameter_name in parameter_names:
            header += '{}: ({} {}) '.format(parameter_name, minimum, maximum)
    else:
        raise TypeError('Min and max ranges can either be provided as ints with a list of parameter names, or as a tuple for each parameter in parameter data (as a dict)')
    
    " Output if required "
    if save:
        # Header with scaling factors 
        with open(filename, 'w') as f:
            comment_char = '#'
            f.write('{} {}\n'.format(comment_char,header)) # Write metadata as comment
            parameter_sets = parameter_sets.round(6) # Round to 6 d.p., good for up to approx a million models with scaling factors of order 1, and is pandas default.
            parameter_sets.to_csv(f)
            
    return parameter_sets, header

def build_empty_parameter_set_details():
    parameter_set_details = {'num_models':None, 'parameter_data':None, 'minimum':None, 'maximum':None, 'save':False, 'output_filename':None}
    return parameter_set_details
    
def build_empty_model_details():
    mechanism_vals = {'nav17tf': {'GNav17':'gna_nav17tf'}, }
    # mechanism val format is {name of mechanism in neuron : {public facing parameter name:name of parameter in neuron}}
    model_details = {'mechanisms':mechanism_vals}
    return  model_details

def construct_parameter_names(list_of_terms):
    parameter_parts = list(itertools.product(*list_of_terms))
    parameter_names = ['_'.join(parameter) for parameter in parameter_parts]    
    return parameter_names
    
""" --- Simulation functions --- """
    
def build_empty_sim_protocols():
    sim_protocols_keys = ['delay', 'amp', 'dur', 'interval', 'num_stims', 'stim_func', 't_stop', 'v_init', 'currents_to_record'] 

def simulation(amp, dur, delay, interval=0, num_stims=1, stim_func=h.IClamp, mechanisms={'kdrtf':1., 'katf':3., 'nav18hw':1.}, t_stop=1000., make_plot=False, plot_type='default', model=None, ions=['Na','K'], mechanisms_to_record=None, record_ionic_currents=True, concs_to_record=None):
    """
    Simulation function for individual IClamp simulations.
    Inputs:
    amp - IClamp amplitude (units?)
    dur - duration of stimulus (ms)
    delay - delay for stimulus start (ms)
    interval - interval between stimuli (ms)
    num_stims - number of stimuli (int)
    stim_func - stimulus function, e.g. h.IClamp, h.IRamp
    mechanisms - list of mechanisms or dict of mechanisms and conductances
    mechanisms_to_record - None or a list of NEURON mechanisms and current names to record and output
    Example of a NEURON mechanism: 'nav18vw'/
    Example of current names: 'ina', 'ica', 'ik'
    ions - which ionic concentrations should be dynamic 
    Mechanisms to record - ionic currents to record
    record_ionic_currents - whether to automatically record total ionic currents for relevent mechanisms 
    (e.g. if a potassium channel current is recorded, if this is true, the total potassium current will also be recorded)
    concs_to_record - if ionic concentrations are to be recorded, provide a list of concentrations here (e.g. cai, nai)
    
    """
    
    # Build a model if one is not supplied
    if not model:
        cell = build_model(mechanisms=mechanisms)
    else:
        cell = model
        
    # Set temperature
    h.celsius = 32. # In line with Davidson et al., 2014, PAIN
    
    # Insert ion concentration mechanisms (replaces use of ion_style, which is really for having different accumulation settings in different sections in multi-section models (see Ch8, NEURON book)
    ion_mechanisms = get_dynamic_ion_mechanisms()   
    for ion in ion_mechanisms:
        if ion in ions:
            cell.insert(ion_mechanisms[ion])
          
    
    # Old ion_style implementation
    """
    if 'K' in ions:
        #oldstyle = h.ion_style("k_ion", 1, 2, 1, 1, 0,sec=cell)
    if 'Na' in ions:
        oldstyle = h.ion_style("na_ion", 1, 2, 1, 1, 0,sec=cell)
    if 'Ca' in ions:
        oldstyle = h.ion_style("ca_ion", 1, 2, 1, 1, 0,sec=cell)
    """
    
    #stims = set_stims(amp=amp, dur=dur, delay=delay, interval=interval, num_stims=num_stims, stim_func=stim_func, cell=cell)
    stims = []
    for i in range(num_stims):
        stims.append(stim_func(0.5, sec=cell))
        stims[-1].dur = dur
        stims[-1].amp = amp # Same total current as for ramp
        stims[-1].delay = delay + i*(dur + interval)

    v,t = set_vt(cell=cell)
    
    " Recording of currents and concentations (and can extend to other state vars if necessary) "
    currents = None
    current_set = None
    if mechanisms_to_record:
        current_set = generate_current_set(mechanisms_to_record)
    currents = record_currents(cell=cell, current_set=current_set)
    
    concs = None
    if concs_to_record:
        concs = record_concs(cell=cell, concs_to_record=concs_to_record)
        
    h.finitialize(-65) # Vital! And has to go after record

    # Now setup ions
#    model.ko = 30.0

    # Run simulation
    neuron.run(t_stop)
#    print(model.ko)
    if make_plot:
        # TODO - sort out simulation_plot to work with currents, currently only default
        # (plot Vm only) is supported
        simulation_plot(t, v, currents=currents, plot_type=plot_type)
        
    output = {'t':np.array(t), 'v':np.array(v)}
    if currents:
        # Recast the currents as an np array so they're not mutable by future simulations.
        recast_currents = recast_recorded_currents(currents)
        output['I'] = recast_currents
    if concs:
        print("UNTESTED USING RECASTING BUT UNTESTED")
        #recast_concs = recast_recorded_currents(concs)
        output['concs'] = concs

    return output
        
        
def simulation_for_ab(amp, dur, delay, interval, num_stims=40, stim_func=h.IClamp, mechanisms={'kdrtf':1., 'katf':3., 'nav18hw':1.}, t_stop=1000.):
    """
    Simulation with plotting for Anabios
    """
    cell = build_model(mechanisms=mechanisms)

    # Set temperature
    h.celsius = 32. # In line with Davidson et al., 2014, PAIN
    
    # Set ionic conditions
    oldstyle = h.ion_style("k_ion", 1, 2, 1, 1, 0,sec=cell)
    oldstyle = h.ion_style("na_ion", 1, 2, 1, 1, 0,sec=cell)
    
    stims = []
    for i in range(num_stims):
        stims.append(stim_func(0.5, sec=cell))
        stims[-1].dur = dur
        stims[-1].amp = amp # Same total current as for ramp
        stims[-1].delay = delay + i*(dur + interval)

    v = h.Vector()
    t = h.Vector()

    v.record(cell(0.5)._ref_v, sec=cell)
    t.record(h._ref_t)

    inav18 = h.Vector()
    ik = h.Vector()
    ikdr = h.Vector()
    ia = h.Vector()

    inav18.record(cell(0.5)._ref_ina_nav18hw, sec=cell)
    ik.record(cell(0.5)._ref_ik, sec=cell)
    ikdr.record(cell(0.5)._ref_ik_kdrtf, sec=cell)
    ia.record(cell(0.5)._ref_ik_katf, sec=cell)


    h.finitialize(-65) # Vital! And has to go after record

    # Run simulation
    neuron.run(t_stop)
    
    plt.figure(figsize=(10,5))
    plt.subplot(3,1,1)
    plt.plot(t,ik)
    plt.plot(t,ikdr)
    plt.plot(t,ia)
    plt.plot(t,inav18)

    plt.subplot(3,1,2)
    plt.plot(t,v); plt.ylabel('Vm (mV)')


    plt.subplot(3,1,3)
    plt.plot(t,ia)
    plt.plot(t,ik)
    return
    
def simulation_vclamp():
    """ Run a voltage clamp simulation """    

""" --- Utility functions ---- """
    
def get_biomarker_list(biomarker_set='default'):
    """ 
    To do - return lists of commonly used biomarker combinations to make it easier to 
    get populations of models to calculate and calibrate off of standard biomarker combinations (e.g. the biomarkers used in Davidson et al.).
    """
    
    biomarker_list = [
        '',
        '',
        '',
        '',
        '',
        '',
    ]
    return biomarker_list
    
def get_ionic_current_list():
    """
    Names of supported ionic currents (whole ions, e.g. Na, Ca, K, not from individual mechanisms
    like Nav 1.7.
    """
    return ['ina', 'ik', 'ica']
    
def get_ionic_concentration_list():
    """
    Names of supported ionic concentrations
    """
    return ['nai', 'ki', 'cai']
    
def get_dynamic_ion_mechanisms():
    """ 
    Returns dict of mechanisms for dynamic ionic concentrations.
    Format is {ion_name:mechanism_name}
    """
    ion_mechanisms = {'K':'k_conc', 'Na':'na_conc', 'Ca':'ca_conc'}
    return ion_mechanisms

def get_default_simulation_kwargs(amp=None, model=None):
    # Doesn't include amp as this needs to calculated
    default_kwargs = {
            'dur':500.,
            'delay':1000.,
            'interval':0.,
            'num_stims':1,
            't_stop':1500.,
            'mechanisms':None,
            'make_plot':False,
            'plot_type':'default'}
    if amp is not None: default_kwargs['amp'] = amp
    if model: default_kwargs['model'] = model
    return default_kwargs
    

def scale_parameter(model, parameter_name, scaling_factor):
    # Can modify parameter inplace
    exec('model.{0} *= {1}'.format(parameter_name, scaling_factor))

def recast_recorded_currents(currents):
    """ Utility function to recast recorded currents as an np.array
    so the results can't be changed by future simulations. """
    recast_currents = {}
    for cur in currents:
        recast_currents[cur] = {}
        for cur_component in currents[cur]:
            recast_currents[cur][cur_component] = np.array(currents[cur][cur_component])    
    return recast_currents

def multiply_parameter(model, parameter, multiplier):
    """ Multiply a parameter by a given value. """
    assert type(parameter) == str
    val = getattr(model,parameter)
    setattr(model,parameter,val*multiplier)

def convert_list_to_string(list_of_lists):
    """ 
    Convert a list of lists of strings into a list of string
    This is a placeholder to use for firing patterns until I recode and rerun sims
    to use strings rather than lists of strings.
    This is because lists as data elements don't play nice with pandas.
    You can't do df['column'] == ['a','b'] as it thinks you're equality checking the entire Series
    rather than element by element as something like df['column'] == 'a, b' would do.
    """
    def convert(list_of_strs):
        if isinstance(list_of_strs, list):
            return ', '.join(list_of_strs)
        else: # E.g. for nans
            return list_of_strs

    return [convert(list_of_strs) for list_of_strs in list_of_lists]


def construct_simulation_set(sim_factors):
    """
    Construct a set of simulations using all combinations of factors
    sim_factors should be a dict with the name of each parameter and it's range of values 
    use an OrderedDict if you want a specific ordering of factors in the output
    """

    sim_parameter_combinations = list(itertools.product(*[val for _,val in sim_factors.items()]))

    sims = np.array([])
    for sim_parameter_set in sim_parameter_combinations:
        sim_details = {}
        for i, name in enumerate(sim_factors):
            sim_details[name] = sim_parameter_set[i]
        sims = np.append(sims, sim_details) 
    
    return sims


def divide_into_sections(simulations, n):
    """
    Chunk set of simulations into n sections
    Returns:
    Dictionary of sim_numbers with keys from 1 to n.
    """

    # Divide simulation_parameters into n segments
    num_sims = len(simulations)
    base_num_sims_per_section = num_sims // n
    num_sections_with_one_extra = num_sims % n

    sim_nums = {}
    num_sims_allocated = 0
    for i in range(1,n+1): # Sections are not 0-indexed
        if i < num_sections_with_one_extra:
            num_sims_in_section = base_num_sims_per_section + 1
        else:
            num_sims_in_section = base_num_sims_per_section

        sim_start = num_sims_allocated
        sim_end = num_sims_allocated + num_sims_in_section

        sim_nums[i] = np.arange(sim_start,sim_end) 
        num_sims_allocated += num_sims_in_section # Keep track of which simulations are distributed
    return sim_nums


def get_simulation_section(simulations, section, n):
    """
    Get a section of a set of simulations chunked into n sections
    Section starts from 1. 
    """
    assert (section > 0) & (section <= n), "section should start from 1 and go to n"
    sim_sections = divide_into_sections(simulations, n)
    simulation_section = simulations[sim_sections[section]] # This is a bug that prevents doing the first sectioon
    return simulation_section


'''
def build_model(mechanisms=['nav17vw', 'nav18hw', 'kdrtf'], conductances=[1.0,1.0,1.0]):
    cell = h.Section()
    for mechanism, conductance in zip(mechanisms,conductances):
        #print mechanism, conductance
        #print 'cell.gbar_{0} = cell.gbar_{0}*conductance'.format(mechanism)
        cell.insert(mechanism)
        exec('cell.gbar_{0} = cell.gbar_{0}*conductance'.format(mechanism)) 

    cell.L = 30 # microns
    cell.diam = 23*2 #microns
    cell.cm = (20.2)/1000  # uF/cm2
    cell.Ra = 123 # ohm-cm

    # Ionic concs
    cell.ko = 5.0
    cell.ki = 140.0

    cell.nao = 140.0
    cell.nai = 10.0
    return cell

def build_stim_protocol(cell,delay,dur,amp,stim_type='Ramp', pos=0.5):
    # Set up stimulus mechanism
    if stim_type == 'Ramp':
        stim = h.IRamp(pos,sec=cell)
    elif stim_type == 'Step':
        stim = h.IClamp(pos,sec=cell)
    else:
        raise NameError('Stim type not found.')
        
    # Parameterise stimulus mechanism
    stim.delay = delay
    stim.dur = dur
    stim.amp = amp
    
    return stim

def run_simulation(cell,t_stop=1000.0):
    
    v = h.Vector()
    t = h.Vector()
    v.record(cell(0.5)._ref_v, sec=cell)
    t.record(h._ref_t)

    h.finitialize(-65) # Vital! And has to go after record
    neuron.run(t_stop)
    return t,v
'''
