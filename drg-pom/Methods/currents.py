# -*- coding: utf-8 -*-
"""
currents - library for building current models and analysing both experimental and simulated voltage clamp recordings

Created on Mon Apr 10 16:30:04 2017

@author: Oliver Britton
"""

import os
import pandas as pd
import numpy as np



" Voltage clamp - set up for simulations and helper functions " 

def voltage_clamp(t,params):
    "Defines a square voltage protocol with an arbitrary number of voltage steps"
    assert len(params['t']) == len(params['v']) # check that each time point has an associated voltage point
    # Go through each threshold and set v accordingly. Inefficient if n is large. 
    for t_thresh,v_thresh in zip(params['t'],params['v']):
        if t >= t_thresh:
            v = v_thresh
    return v

def get_voltage_clamp_vector(t,voltage_func,voltage_func_params):
    """ 
    Gets the voltage clamp vector to match a given vector of times, t, using the
    voltage clamp function and parameterisation given as input.
    """
    v = np.zeros(t.shape)
    t = np.array(t)
    for i,t_val in enumerate(t):
        v[i] = voltage_func(t_val,voltage_func_params)
    return v

def make_voltage_clamp_list(start, stop, delta, index, voltages, times):
    """
    Build a list of voltage clamp parameter dictionaries for a voltage clamp step protocol with 1 variable step
    """
    # Build list of dictionaries with different voltage steps
    voltage_clamp_list = []
    for v_delta in np.arange(start,stop+0.1,delta):
        voltages[index] = v_delta
        voltage_clamp_list.append({'v':voltages[:], 't':times, 'step':v_delta}) # [:] does a shallow copy of the list so we don't change all elements
    return voltage_clamp_list
        
def predefined_voltage_clamp_list(current,protocol):
    " Build predefined voltage clamp lists for particular protocols "
    if current == 'ttx_s':
        if protocol == 'act':
            return make_voltage_clamp_list(start=-90.,stop=10.,delta=5.,index=1,voltages=[-60.,None,-60.], times=[0,100,200])            
        # These two protocols have times beginning at negative points as 0 is the point at which experimental recording starts, but the protocol begins earlier 
        elif protocol == 'fast inact':
            return make_voltage_clamp_list(start=-120.,stop=-10.,delta=5.,index=0,voltages=[None,-10.,-120.], times=[-500,100,125])
            
        elif protocol == 'slow inact':
            return make_voltage_clamp_list(start=-90.,stop=40.,delta=10.,index=0,voltages=[None,-120.,-10.,-120.], times=[-8000,0,100,125])
    return False


" Voltage clamp analysis "

" Building and simulating current models "

" Specific current models "

" Loading experimental data "

def load_voltage_clamp_recording(filename):
    return pd.read_csv(filename,header=1).set_index('Time (ms)')

def get_experimental_voltage_clamp_recording(current='ttx_s',protocol='act'):
    " Get a particular voltage clamp recording or set of voltage clamp recordings "
    
    data = False
    if current == 'ttx_s':
        filepath = 'E:\\CLPC48\\Neuron Project\\Data\\Preliminary' 
        if protocol == 'act':
            filename = os.path.join(filepath,'EP2_2014-01-21_03_Activation.csv')
            data = load_voltage_clamp_recording(filename)
        elif protocol == 'fast inact':
            filename = os.path.join(filepath,'EP2_2014-01-21_03_FastInactivation.csv')
            data = load_voltage_clamp_recording(filename)
        elif protocol == 'slow inact':
            filename = os.path.join(filepath,'EP2_2014-01-21_03_SlowInactivation.csv')
            data = load_voltage_clamp_recording(filename)

    return data # Return whatever data we have, if current/protocol were not found return false

    
    

