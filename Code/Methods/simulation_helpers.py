# -*- coding: utf-8 -*-
"""
Simulation helpers 

Created on Thu Jun 29 11:17:20 2017

@author: Oliver Britton
"""

import os
import sys
import numpy as np
from matplotlib import pyplot as plt
import pyDOE

# Load Neuron and DRG code
sys.path.append('E:\\CLPC48\\Neuron Project')
sys.path.append('E:\\CLPC48\\Neuron Project\Code')
import NeuronProjectStart

from neuron import h
import neuron

def init_model(mechanisms,
    L=30., 
    diam=46., 
    cm=20.2/1000., 
    Ra=123., 
    ko=3., 
    ki=135., 
    nao=145., 
    nai=5.,
    ):
    """
    Initialise a cylindrical cell model, with the specified ionic mechanisms embedded.
    """
    cell = h.Section()
    for mechanism in mechanisms:
        cell.insert(mechanism)
    cell.L = L # microns
    cell.diam = diam #microns
    cell.cm = cm  # uF/cm2
    cell.Ra = Ra # ohm-cm
    
    cell.ko = ko
    cell.ki = ki
    cell.nao = nao
    cell.nai = nai
    return cell

def build_model(mechanisms={'kdrtf':1., 'katf':1., 'nav18hw':1.}, conductances=None):
    
    # Dict construction
    if type(mechanisms) == dict:
        model = init_model(mechanisms=mechanisms.keys())
        for mechanism, conductance in mechanisms.iteritems():
            exec('model.gbar_{0} *= {1}'.format(mechanism, conductance)) 
        if conductances:
            assert False
            
    # List construction
    elif type(mechanisms) == list:
        model = init_model(mechanisms=mechanisms)
        if conductances:
            for mechanism, conductance in zip(mechanisms,conductances):
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
        plt.figure(figsize=(10,5))
        plt.subplot(3,1,1)
        plt.plot(t,currents['ik'])
        plt.plot(t,currents['ikdr'])
        plt.plot(t,currents['ia'])
        plt.plot(t,currents['inav18'])

        plt.subplot(3,1,2)
        plt.plot(t,v); plt.ylabel('Vm (mV)')

        plt.subplot(3,1,3)
        plt.plot(t,currents['ia'])
        plt.plot(t,currents['ik'])
        
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
        for i in ['ik', 'ikdr', 'ia', 'inav18']:
            currents[i] = h.Vector()

        currents['inav18'].record(cell(0.5)._ref_ina_nav18hw, sec=cell)
        currents['ik'].record(cell(0.5)._ref_ik, sec=cell)
        currents['ikdr'].record(cell(0.5)._ref_ik_kdrtf, sec=cell)
        currents['ia'].record(cell(0.5)._ref_ik_katf, sec=cell)
        return currents
        
def build_sim_protocols(amp, dur, delay, interval, num_stims=40, stim_func=h.IClamp, t_stop=1000., v_init=-65.,):

    sim_protocols = {'amp':amp, 'dur':dur, 'delay':delay, 'interval':interval, 'num_stims':num_stims, 'stim_func':stim_func, 't_stop':t_stop, 'v_init':v_init,}
    return sim_protocols
    
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

def build_parameter_set(num_models, num_parameters, minimum, maximum, output_filename=None, parameter_names=None, save=False):
    
    parameter_sets = pyDOE.lhs(num_parameters, num_models)
    # Transform parameter set using minimumimum and maximumimum
    # Scale size of range
    parameter_sets *= (maximum-minimum)
    # Scale to right minimumimum (and right maximumimum if we scaled range right)
    parameter_sets += minimum    
    
    if output_filename:
        if parameter_names:
            assert num_parameters == len(parameter_names), "Number of parameters != length of parameter names list."
            header = ', '.join(parameter_names)
        else:
            header = ''
        
        if save:
            np.savetxt(output_filename, parameter_sets, fmt='%.3f', header=header, comment='') # Write to 3 d.p precision 
    return parameter_sets
    

def simulation(amp, dur, delay, interval, num_stims=40, stim_func=h.IClamp, mechanisms={'kdrtf':1., 'katf':3., 'nav18hw':1.}, t_stop=1000., make_plot=True, plot_type='default', model=None):
    
    # Build a model if one is not supplied
    if not model:
        cell = build_model(mechanisms=mechanisms)
    else:
        cell = model
        
    # Set temperature
    h.celsius = 32. # In line with Davidson et al., 2014, PAIN
    
    # Set ionic conditions
    oldstyle = h.ion_style("k_ion", 1, 2, 1, 1, 0,sec=cell)
    oldstyle = h.ion_style("na_ion", 1, 2, 1, 1, 0,sec=cell)
    
    #stims = set_stims(amp=amp, dur=dur, delay=delay, interval=interval, num_stims=num_stims, stim_func=stim_func, cell=cell)
    stims = []
    for i in range(num_stims):
        stims.append(stim_func(0.5, sec=cell))
        stims[-1].dur = dur
        stims[-1].amp = amp # Same total current as for ramp
        stims[-1].delay = delay + i*(dur + interval)

    v,t = set_vt(cell=cell)
    currents = record_currents(cell=cell, current_set=plot_type)

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