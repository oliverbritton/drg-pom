"""
Example configuration file to run a simulation saving the traces of all the ionic currents for each model. Be warned this will use a lot of disk space but is necessary if you want
to analyse the ionic currents and plot currentscapes. Check the comments for this script to see where to make changes to a normal simulation to save currents.
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
    name = 'example_save_all_ionic_currents'
    pop.name = name
    save_filename = '{}.pkl'.format(name)

    # Simulation parameters
    save_dir = None
    benchmark = True
    rerun = False

    " ---// Setup for recording all ionic currents \\--- "
    save_type = 'both' # Save both trace data files and figures - we need to save the trace data files to save the ionic current data
    outputs = list(pop.mechanisms.keys()) # This tells the simulator to save all "mechanisms" as outputs (mechanisms are what NEURON calls ionic current models).
    print(f"Outputs are: {outputs}") # Prints out which currents we are saving.
    " -------------------------------------------------- "

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
         'outputs': outputs, # This line is vital - here we pass the list of ionic currents to save to the simulation.
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
         'outputs': outputs, # This line is vital - here we pass the list of ionic currents to save to the simulation.
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
