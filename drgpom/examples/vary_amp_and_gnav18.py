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

    pop_filename = os.path.join('data', 'example_population.pkl')
    pop = drg.load(pop_filename)
    name = 'vary_amp_and_gnav18'
    pop.name = name
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
    count = 0
    for simulation in simulations:
        count += 1
        print(f"Running simulation {count} of {len(simulations)}.")
        parameter_scaling = {key:val for key,val in simulation.items() if key in pop.parameter_designations.keys()}
        step_amp  = simulation['Step Amp']

        str_parameter_scaling = '_'.join(['{}_{}'.format(key,val) for key,val in parameter_scaling.items()])
        str_step_amp = '{}_pA'.format(int(step_amp*1000)) # Convert from nA to pA

        sim_name = name + '_step_{}_{}'.format(str_step_amp, str_parameter_scaling)
        sim_type = 'iclamp'
        protocols = {
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

        # Add simulation protocols for step
        protocols['parameter_scaling'] = parameter_scaling
        protocols['amp'] = step_amp
        print("Running simulation {}, with parameter scaling: {} and amp: {}".format(sim_name, protocols['parameter_scaling'], protocols['amp']))


        pop.run_simulation(name=sim_name, 
                       simulation_type=sim_type, 
                       protocols=protocols,
                       cores=cores, 
                       save_type=save_type,
                       save_dir=save_dir,
                       benchmark=benchmark, 
                       rerun=rerun)

    " --- Save Results --- "

    print(pop.results.head())
    print("Time taken on {} cores = {}s.".format(cores,time.time()-start))
    pop.save(save_filename)
    print("Current population saved to: {}".format(save_filename))

