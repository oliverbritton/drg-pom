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
# Load ion channel models
import Methods.Simulations.loadneuron as ln
ln.load_neuron_mechanisms(verbose=False)

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
                # Split on underscores and use the last element as we assume mechanism names do not contain underscores
                mechanism_names = [name.split('_')[-1] for name in mechanism_names]
        model = init_model(mechanisms=mechanism_names)
        for mechanism, conductance in mechanisms.items():
            if mechanism_is_full_parameter_name: # Mechanisms contains the full parameter name, not just the mechanism name
                exec('model.{0} *= {1}'.format(mechanism, conductance)) 
            else: # Mechanism is assumed to be a suffix for a conductance
                exec('model.gbar_{0} *= {1}'.format(mechanism, conductance)) 
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
        
def simulation_plot(t, v, currents, plot_type='default'):
    if plot_type == 'default':
    
        # Need following currents:
        # ik, ikdr, ia, inav18
        plt.figure(figsize=(5,12))
        plt.subplot(4,1,1)
        plt.plot(t,currents['ina'])
        plt.plot(t,currents['inav17'])
        plt.plot(t,currents['inav18'])
        plt.plot(t,currents['inav19'])
        plt.title('All Navs')

        plt.subplot(4,1,2)
        plt.plot(t,currents['inav17'])
        plt.plot(t,currents['inav19'])
        plt.title('1.7 and 1.9')
        
        plt.subplot(4,1,3)
        plt.plot(t,v); plt.ylabel('Vm (mV)')
        plt.title('Vm')

        plt.subplot(4,1,4)
        plt.plot(t,currents['ik'])
        plt.plot(t,currents['ikdr'])
        plt.plot(t,currents['ia'])
        plt.plot(t,currents['im'])
        plt.title('All Kvs')
        
    if plot_type == 'simple':
        plt.figure(figsize=(5,5))
        plt.plot(t,v); plt.ylabel('Vm (mV)'); plt.xlabel('Time (ms)')
        
def set_vt(cell):
    v = h.Vector()
    v.record(cell(0.5)._ref_v, sec=cell)
    t = h.Vector()
    t.record(h._ref_t)
    return v,t
    
def record_currents(cell, current_set='default'):
    currents = {}
    if current_set == 'default':
        for i in ['ina', 'inav17', 'inav18', 'inav19', 'ik', 'ikdr', 'ia', 'im']:
            currents[i] = h.Vector()

        currents['inav17'].record(cell(0.5)._ref_ina_nav17vw, sec=cell)
        currents['inav18'].record(cell(0.5)._ref_ina_nav18hw, sec=cell)
        currents['inav19'].record(cell(0.5)._ref_ina_nav19hw, sec=cell)
        currents['ina'].record(cell(0.5)._ref_ina, sec=cell)
        currents['ik'].record(cell(0.5)._ref_ik, sec=cell)
        currents['ikdr'].record(cell(0.5)._ref_ik_kdrtf, sec=cell)
        currents['ia'].record(cell(0.5)._ref_ik_katf, sec=cell)
        currents['im'].record(cell(0.5)._ref_ik_kmtf, sec=cell)
    elif current_set == 'simple':
        for i in ['ina', 'ik']:
            currents[i] = h.Vector()
        currents['ina'].record(cell(0.5)._ref_ina, sec=cell)
        currents['ik'].record(cell(0.5)._ref_ik, sec=cell)
    elif current_set == None:
        pass
    else:
        raise ValueError('Current set not in list of accepted values.')
    return currents
    
def generate_current_set(mechanisms):
    """
    Take the mechanism list of a simulation and work out which currents to record
    based on their names.
    So if there is a current that begins with k its a k current (_ref_ik_k%NAME%)
    If it begins with na its an na current.
    If it doesn't follow a know pattern don't record it.
    """
    pass
        
def build_sim_protocols(amp, dur, delay, interval, num_stims=40, stim_func=h.IClamp, t_stop=1000., v_init=-65.,):

    sim_protocols = {'amp':amp, 'dur':dur, 'delay':delay, 'interval':interval, 'num_stims':num_stims, 'stim_func':stim_func, 't_stop':t_stop, 'v_init':v_init,}
    return sim_protocols

" -- Functions related to building parameter sets -- "
    
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
        comment_char = '#'
        with open(filename, 'w') as f:
            f.write('{} {}\n'.format(comment_char,header)) # Write header as comment
            parameter_sets = parameter_sets.round(6) # Round to 6 d.p., good for up to approx a million models with scaling factors of order 1, and is pandas default.
            parameter_sets.to_csv(f, comments=comment_char)
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

def simulation(amp, dur, delay, interval, num_stims=40, stim_func=h.IClamp, mechanisms={'kdrtf':1., 'katf':3., 'nav18hw':1.}, t_stop=1000., make_plot=True, plot_type='default', model=None, ions=['Na','K'], recording_set=None):
    
    # Build a model if one is not supplied
    if not model:
        cell = build_model(mechanisms=mechanisms)
    else:
        cell = model
        
    # Set temperature
    h.celsius = 32. # In line with Davidson et al., 2014, PAIN
    
    # Set ionic conditions
    if 'K' in ions:
        oldstyle = h.ion_style("k_ion", 1, 2, 1, 1, 0,sec=cell)
    if 'Na' in ions:
        oldstyle = h.ion_style("na_ion", 1, 2, 1, 1, 0,sec=cell)
    if 'Ca' in ions:
        oldstyle = h.ion_style("ca_ion", 1, 2, 1, 1, 0,sec=cell)
        
    #stims = set_stims(amp=amp, dur=dur, delay=delay, interval=interval, num_stims=num_stims, stim_func=stim_func, cell=cell)
    stims = []
    for i in range(num_stims):
        stims.append(stim_func(0.5, sec=cell))
        stims[-1].dur = dur
        stims[-1].amp = amp # Same total current as for ramp
        stims[-1].delay = delay + i*(dur + interval)

    v,t = set_vt(cell=cell)
    
    if recording_set == True:
        current_set = generate_current_set # @To do - make it generate from mechanism list
    elif make_plot == True:
        current_set = plot_type
    else:
        current_set = None
    currents = record_currents(cell=cell, current_set=current_set)

    h.finitialize(-65) # Vital! And has to go after record
    # Run simulation
    neuron.run(t_stop)
    
    if make_plot:
        simulation_plot(t, v, currents, plot_type=plot_type)
    
    return np.array(t), np.array(v)

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

""" --- Other functions ---- """
    
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
"""
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
"""

" --- Tests --- "

def test_build_parameter_set():
    build_parameter_set(10,5,0,2,output_filename='test.csv',parameter_names=['GNa', 'GKr', 'G3', 'G4', 'G5'],save=True)
    print("Check that test.csv opens is a formatted 10 row, 5 column grid with a header of parameter names")
    return