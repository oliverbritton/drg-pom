"""
Example configuration file to run a new population of models simulation
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

    " --- Define and Construct Population --- "

    name = 'example_population_fullsize'
    save_filename = '{}.pkl'.format(name)

    # Parameter set settings - set either parameter_data and num_models, or parameter_filename to
    # a non-None value, as generating population from parameter data is mututally exclusive with
    # loading population parameters from an existing parameter set file. 
    num_models = 20000
    parameter_data = {'GNav17':[0.0,0.4], 'GNav18':[0.,4.0], 'GNav19':[0.,4.],
            'GKdr':[0.,4.], 'GKA':[0.,40.], 'GKM':[0.,4.], 'GH':[0.,2.], 'GKleak':[0., 0.2]}


    parameter_filename = None

    # Defining parameters to vary for each ionic current and the names of each parameter
    model_details = {'mechanisms':{}}
    model_details['mechanisms']['nav17vw_named'] = {'GNav17':'gbar_nav17vw_named'}
    model_details['mechanisms']['nav18hw_named'] = {'GNav18':'gbar_nav18hw_named'}
    model_details['mechanisms']['nav19hw'] = {'GNav19':'gbar_nav19hw'}
    model_details['mechanisms']['kdrtf'] = {'GKdr':'gbar_kdrtf'}
    model_details['mechanisms']['katf'] = {'GKA':'gbar_katf'}
    model_details['mechanisms']['kmtf'] = {'GKM':'gbar_kmtf'}
    model_details['mechanisms']['hcntf'] = {'GH': 'gbar_hcntf'}
    model_details['mechanisms']['kleak'] = {'GKleak': 'gbar_kleak'}

    # Simulation parameters
    save_type = 'fig' # Allowed types are 'fig', 'trace', 'both', or 'none'
    save_dir = None
    benchmark = True
    rerun = False
    outputs = [] 

    if parameter_data is not None:
        parameter_set_details = {}
        parameter_set_details['num_models'] = num_models
        parameter_set_details['parameter_data'] = parameter_data 
        parameter_set_details['save'] = True
        parameter_set_details['output_filename'] = 'example_population_creation_fullsize_parameters.csv'   
    else:
        parameter_set_details = None 

    pop = drg.PopulationOfModels(name=name, 
                                 simulation_protocols=None, 
                                 model_details=model_details,
                                 parameter_filename=parameter_filename,
                                 parameter_set_details=parameter_set_details)

    " --- Run Simulations --- "

    # Ramp stimulus simulation
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
                   rerun=False)

    # Step stimulus simulation
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

    " --- Population Calibration --- "

    biomarkers_to_calibrate = {
        'ramp':['Threshold'],
        'step':['AHPTau', 
                 'APPeak',
                 'APSlopeMax',
                 'APSlopeMin',
                 'APFullWidth',
                 'RMP',
                 'Rheobase'],
    }
    calibration_ranges = 'Davidson'
    std = 1.5
    verbose_calibration = False
    for sim_name, biomarker_names in biomarkers_to_calibrate.items():
        pop.calibrate_population(
                biomarker_names=biomarker_names,
                simulation_name=sim_name,
                calibration_ranges=calibration_ranges,
                stds=std,
                verbose=verbose_calibration,
                )    


    " --- Save Population --- "

    print(pop.results.head())
    print("Time taken on {} cores = {}s.".format(cores,time.time()-start))
    pop.save(save_filename)
    print("Current population saved to: {}".format(save_filename))
