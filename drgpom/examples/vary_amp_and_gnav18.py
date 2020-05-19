"""
Example configuration file to run a set of simulations 
"""

import matplotlib
matplotlib.use('Agg')

import time
import sys
import os
import multiprocessing as mp

import drgpom as drg

# Define number of cpu cores to use
# If not specified default to cpu count - 1
if __name__ == '__main__':
    sys.modules['__main__'].__spec__ = None # Without this repeated runnings of this script in Python 3 fail
    # Set cores to cpu count or specified number
    cores = mp.cpu_count()-1
    if len(sys.argv) >= 2:
        for arg in sys.argv[1:]: 
            if 'cores=' in arg:
                core_str = arg.rsplit('=')[-1]
                cores = int(core_str)
                print("Cores = {}".format(cores))

    start = time.time() # For benchmarking
    
    " //-- Main Program begins here --// "

    " --- Load existing population of models --- "

    pop_filename = os.path.join('example_population.pickle')
    pop = drg.load(pop_initial_filename)
    name = 'vary_amp_and_gnav18'
    save_filename = '{}.pkl'.format(name)

    """
    Set up parameters for multiple simulations
    We will provide a list of GNav 1.8 scaling factors 
    and a list of stimulus amplitudes for a step stimulus function.
    Then we will run a simulation for every possible combination of the 
    two lists.
    """
    gnav18s = [0.5, 1.0, 1.5]
    step_amplitudes = [1.0, 2.0, 3.0, 4.0, 5.0] # nA
    simulation_factors = {'GNav18':gnav18s, 'Step Amp':step_amplitudes}
    simulations = drg.construct_simulation_set(simulation_factors)

    # Simulation parameters
    save_type = 'fig' # Allowed types are 'fig', 'trace', 'both', or 'none'
    save_dir = None
    benchmark = True
    rerun = False
    outputs = [] 
 
    " --- Run Simulations --- "

    " --- Simulation 1 - Example Ramp Stimulus --- "
    sim_name = 'example_ramp'
    sim_type = 'iclamp'
    protocols = {
         'amp': None,
         'celsius': 32.0,
         'delay': 500,
         'dur': 500,
         'interval': 0,
         'ions': ['Na','K'],
         'num_stims': 1,
         'outputs': outputs,
         'sampling_freq': 20000,
         'stim_func': 'h.IRamp',
         't_stop': 1500.,
         'v_init': -65.,
         }

    pop.run_simulation(name=sim_name, 
                   simulation_type=sim_type, 
                   protocols=protocols,
                   cores=cores, 
                   save_type=save_type, 
                   save_dir=save_dir, 
                   benchmark=benchmark, 
                   rerun=rerun)


    " --- Simulation 2 - Example Step Stimulus --- "
    sim_name = 'example_step' 
    sim_type = 'iclamp'
    protocols = {
         'amp': None,
         'celsius': 32.0,
         'delay': 500,
         'dur': 800,
         'interval': 0,
         'ions': ['Na','K'],
         'num_stims': 1,
         'outputs': outputs,
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
                   save_type=save_type, 
                   save_dir=save_dir, 
                   benchmark=benchmark, 
                   rerun=False)

    " --- Save Results --- "

    print(pop.results.head())
    print("Time taken on {} cores = {}s.".format(cores,time.time()-start))
    pop.save(save_filename)
    print("Current population saved to: {}".format(save_filename))

