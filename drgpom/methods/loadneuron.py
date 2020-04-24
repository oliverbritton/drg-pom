# -*- coding: utf-8 -*-
"""
loadneuron.py
Functions to load Neuron properly.
Created on Fri Oct 28 12:12:50 2016 @author: Oliver Britton
"""
import os
import sys
from neuron import h

if 'mechanisms_loaded' not in locals():
    mechanisms_loaded = False

# Get platform
if sys.platform == 'win32':
    platform = 'windows'
elif sys.platform == 'linux':
    platform = 'linux'
else:
    raise Exception('Platform {} not supported.'.format(sys.platform))
    
    
def load_neuron_mechanisms(nrnmech_path=None, verbose=False, ):
    """
    Main loading function, loads neuron mechanisms, iff they are not already loaded.
    Can only call nrn_load_dll once, calling it a second time causes a crash, so need to check
    we have loaded mechanisms before calling.
    """

    if are_mechanisms_loaded() == False:
        if nrnmech_path == None:
            nrnmech_path = os.path.join(get_mechanism_dir(), get_mech_lib_name(platform))
        else:
            # If path is provided, check it contains library name
            if get_mech_lib_name(platform) not in nrnmech_path:
                nrnmech_path =  os.path.join(nrnmech_path, get_mech_lib_name(platform))
            else:
                pass

        print('Loading {} from {}.'.format(get_mech_lib_name(platform), nrnmech_path))
        h.nrn_load_dll(nrnmech_path) # Neuron mechanism loading function

    if verbose:
        if are_mechanisms_loaded():
            print("Mechanisms are loaded.")
        else:
            print("Mechanisms NOT loaded!")

def get_mech_lib_name(platform):
    if platform == 'windows':
        mech_lib_name = 'nrnmech.dll'
    elif platform == 'linux':
        mech_lib_name = 'libnrnmech.so'
    else:
        raise Exception('Platform {} not supported.'.format(sys.platform))
    return mech_lib_name



def get_mechanism_dir():
    """
    Returns the models directory the compiled nrnmech.dll is stored in. 
    """
    cur_path = (os.path.dirname(__file__))
    mechanisms_path = os.path.abspath(os.path.join(cur_path, os.path.pardir, 'models'))
    return mechanisms_path


def check_num_name_mechs():
    """
    Example function for checking number and names of mechanisms.
    Reminds me how to do string refs in hoc and get names of mechanisms.
    """
    from neuron import h
    mt = h.MechanismType(0)
    sref = h.ref('')
    for i in range(int(mt.count())):
        mt.select(i)
        mt.selected(sref)
        print(sref[0])
    return mt.count()
    
    
def are_mechanisms_loaded(verbose=False):
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
        if verbose:
            print("No non-default mechanisms detected, ready to load nrnmech.dll.")
    else:
        loaded = True
    return loaded    

    
def get_mechanism_list():
    """ Return list of loaded mechanism names """
    current_mechanisms = h.MechanismType(0)
    str_ref = h.ref('') # hoc string reference
    current_mechanism_names = []
    for i in range(int(current_mechanisms.count())):
        current_mechanisms.select(i)
        current_mechanisms.selected(str_ref)
        current_mechanism_names.append(str_ref[0])
    return  current_mechanism_names


def default_mech_list():
    """ Returns list of default mechanism names included in NEURON. """
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
    """ Load multiple times, if we don't crash, we're good! """    
    load_neuron_mechanisms()
    load_neuron_mechanisms()
    load_neuron_mechanisms()
    print('We didn\'t crash!')

