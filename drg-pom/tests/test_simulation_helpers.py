import pytest

import Methods.simulation_helpers as sh
from neuron import h
import neuron

def test_build_parameter_set():
    sh.build_parameter_set(num_models=10, parameter_data=['GNa', 'GKr', 'G3', 'G4', 'G5'], minimum=0, maximum=2, filename='test.csv', save=True)
    print("Check that test.csv opens is a formatted 10 row, 5 column grid with a header of parameter names")
    
def test_record_concentrations():
    """
    Test recording concentrations and integration with supporting functions:
    model building, simulations, and lists of ionic concentrations and concentration mechanisms.
    """
    cell = sh.build_model()
    for __, conc_mech in sh.get_dynamic_ion_mechanisms().items():
        cell.insert(conc_mech)
    
    # Record each concentration and make sure we get some output
    concs_to_record = sh.get_ionic_concentration_list()
    concs = sh.record_concs(cell=cell, concs_to_record=concs_to_record)
    
    # Run simulation
    v = h.Vector()
    t = h.Vector()
    v.record(cell(0.5)._ref_v, sec=cell)
    t.record(h._ref_t)
    h.finitialize(-65) # Vital! And has to go after record
    t_stop = 1500
    neuron.run(t_stop)
    
    assert type(concs) == dict
    for __,conc in concs.items():
        assert len(conc) == len(t)
    
    # Record a concentration that's not allowed and make sure we get an error
    
    with pytest.raises(AssertionError):
        concs = sh.record_concs(cell=cell, concs_to_record=['not', 'allowed', 'concentration'])

    
def test_raising_errors():
    # Generic python error is Exception
    with pytest.raises(ZeroDivisionError):
        1 / 0
    
    