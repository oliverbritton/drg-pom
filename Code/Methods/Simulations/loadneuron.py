# -*- coding: utf-8 -*-
"""
NEURON Load

Functions to load Neuron properly

Created on Fri Oct 28 12:12:50 2016 @author: Oliver Britton
"""
import os
from neuron import h

if 'mechanisms_loaded' not in locals():
    mechanisms_loaded = False

" Main load function "
def load_neuron_mechanisms(mech_type='prototype'):
    """
    Load neuron mechanisms, if not already loaded.
    """
    #global mechanisms_loaded
    if are_mechanisms_loaded() == False:
        h.nrn_load_dll(os.path.join(get_mechanism_dir(),"nrnmech.dll"))
        "Mechanisms loaded!"

def get_mechanism_dir(mech_type='prototype'):
    """
    Returns the directory nrnmech.dll is stored in. 
    Inputs:
    mech_type = 'prototype' - Defines which location to look in.  
    TODO - look relative to project directory.
    """
    " Leave room for different dirs in the future "
    if mech_type == 'prototype':
        mech_dir = 'E:\\CLPC48\\Neuron Project\\Code\\Models\\Currents\\Prototypes'
    return mech_dir
    
" Prototype for checking number and names of mechanisms "
def thingy():
    """
    Remind me how to do string refs in hoc and get names of mechanisms.
    """
    from neuron import h
    mt = h.MechanismType(0)
    sref = h.ref('')
    for i in range(int(mt.count())):
        mt.select(i)
        mt.selected(sref)
        print sref[0]
    return mt.count()
    
def are_mechanisms_loaded():
    """
    Checks if mechanisms have been loaded. Assumes a certain set of mechanisms are loaded by default
    and that if a dll has been loaded it will contain other mechanisms. 
    If an empty dll is loaded this won't think it's been loaded and could cause a crash.
    """


    " Get names of mechanisms that are always loaded "
    default_mechanisms = default_mech_list() 
    
    " Get names of currently loaded mechanims "
    current_mechanisms = h.MechanismType(0)
    str_ref = h.ref('') # hoc string reference
    current_mechanism_names = []
    for i in range(int(current_mechanisms.count())):
        current_mechanisms.select(i)
        current_mechanisms.selected(str_ref)
        current_mechanism_names.append(str_ref[0])
    
    " Get only non-default mechanisms "
    non_default_mechanisms = [mech for mech in current_mechanism_names if mech not in default_mechanisms]
        
    " If non_default_mechanisms is empty we haven't loaded (or we've loaded no new mechanisms!) "
    if not non_default_mechanisms :
        # List is empty
        loaded = False
        print("No non-default mechanisms detected, ready to load nrnmech.dll.")
    else:
        loaded = True
    return loaded    


def default_mech_list():
    """
    Returns list of default mechanism names included in NEURON.
    """
    default_mechanisms = ['morphology',
                          'capacitance',
                          'pas',
                          'extracellular',
                          'fastpas',
                          'na_ion',
                          'k_ion',
                          'hh',
                         ]
    return default_mechanisms
    
" -- Tests --- "
def test_load_mechanisms():
"""
Load multiple times, if we don't crash, we're good!
"""    
    load_neuron_mechanisms()
    load_neuron_mechanisms()
    load_neuron_mechanisms()
    

     