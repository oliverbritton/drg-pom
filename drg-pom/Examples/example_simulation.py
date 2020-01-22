"""
Example configuration file to run a new population of models simulation
"""

# Linux
import matplotlib
matplotlib.use('Agg')

import time
import sys
import os
import multiprocessing as mp


# TODO - update with new path setup cod
# Add main path
if 'Dropbox' in os.getcwd():
    sys.path.append(os.path.join(os.getcwd().split('Dropbox')[0],'Dropbox'))
else:
    raise Exception('Not running in Dropbox directory.')

# Import POM module requirements and NEURON
import neuron_paths
neuron_paths.setup_paths()
import Methods.pom_test as pom
import Methods.simulation_helpers as sh
from neuron import h

# Setup parallelisation for the required number of cores
# If not specified default to cpu count - 1
if __name__ == '__main__':
    sys.modules['__main__'].__spec__ = None # Without this something goes wrong with repeated runnings of this script in Python3
    # Set cores to cpu count or specified number
    cores = mp.cpu_count()-1
    if len(sys.argv) >= 2:
        for arg in sys.argv[1:]: 
            if 'cores=' in arg:
                core_str = arg.rsplit('=')[-1]
                cores = int(core_str)
                print("Cores = {}".format(cores))
            
    """
    Initialise population
    -------------------
    """
    start = time.time() # For benchmarking
    sim_prefix = 'example_population'
    pop_save_filename = '{}.pickle'.format(sim_prefix)
    # Don't calibrate, load full parameter set from 2_1

    parameter_names = {'GNav17':[0.0,0.4], 'GNav18':[0.,4.0], 'GNav19':[0.,4.],
            'GKdr':[0.,4.], 'GKA':[0.,40.], 'GKM':[0.,4.], 'GH':[0.,2.], 'GKleak':[0., 0.2]}

    model_details = {'mechanisms':{}}
    model_details['mechanisms']['nav17vw_named'] = {'GNav17':'gbar_nav17vw_named'}
    model_details['mechanisms']['nav18hw_named'] = {'GNav18':'gbar_nav18hw_named'}
    model_details['mechanisms']['nav19hw'] = {'GNav19':'gbar_nav19hw'}
    model_details['mechanisms']['kdrtf'] = {'GKdr':'gbar_kdrtf'}
    model_details['mechanisms']['katf'] = {'GKA':'gbar_katf'}
    model_details['mechanisms']['kmtf'] = {'GKM':'gbar_kmtf'}
    model_details['mechanisms']['hcntf'] = {'GH': 'gbar_hcntf'}
    model_details['mechanisms']['kleak'] = {'GKleak': 'gbar_kleak'}

    parameter_set_details = {}
    parameter_set_details['num_models'] = 20000 
    parameter_set_details['parameter_data'] = parameter_names 
    parameter_set_details['minimum'] = None
    parameter_set_details['maximum'] = None
    parameter_set_details['save'] = True
    parameter_set_details['output_filename'] = 'example_project_parameters.csv'    

    pop_name = sim_prefix
    do_plot = None
    save_pop = True
    save_type = 'fig'
    benchmark = True
    outputs = [] 
    rerun = False

    pop = pom.PopulationOfModels(name=pop_name, 
                                 simulation_protocols=None, 
                                 model_details=model_details,
                                 parameter_filename=None,
                                 parameter_set_details=parameter_set_details)

    """ 
    Running simulations
    """
    
    sim_name = 'ramp'
    sim_type = 'iclamp'
    protocols = {
         'amp': None,
         'celsius': 32.0,
         'delay': 500,
         'dur': 500,
         'interval': 0,
         'ions': ['Na','K'],
         'num_stims': 1,
         'outputs': [],
         'sampling_freq': 20000,
         'stim_func': 'h.IRamp',
         't_stop': 1500.,
         'v_init': -65.,
         }

    pop.run_simulation(name=sim_name, 
                   simulation_type=sim_type, 
                   protocols=protocols,
                   cores=cores, 
                   plot=do_plot, 
                   save=save_pop, 
                   save_type=save_type, 
                   benchmark=benchmark, 
                   rerun=False)


    sim_name = 'step' 
    sim_type = 'iclamp'
    protocols = {
         'amp': None,
         'celsius': 32.0,
         'delay': 500,
         'dur': 800,
         'interval': 0,
         'ions': ['Na','K'],
         'num_stims': 1,
         'outputs': [],
         'sampling_freq': 20000,
         'stim_func': 'h.IClamp',
         't_stop': 1500.,
         'v_init': -65.,
         'flags': {'ramp_threshold_sim_for_width':'ramp'} # Use the threshold from the ramp simulation to compute ap full width
         }

    pop.run_simulation(name=sim_name, 
                   simulation_type=sim_type, 
                   protocols=protocols,
                   cores=cores, 
                   plot=do_plot, 
                   save=save_pop, 
                   save_type=save_type, 
                   benchmark=benchmark, 
                   rerun=False)

    """
    Save and summarise
    """
    print(pop.results.head())
    print("Time taken on {} cores = {}s.".format(cores,time.time()-start))
    pop.pickle_pom(pop_save_filename)
    print("Current population saved to: {}".format(pop_save_filename))
