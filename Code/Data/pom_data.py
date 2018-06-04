"""
pom data 

functions to load specific populations etc.
Oliver Britton 15/02/2017
"""
import os 
import Methods.pom_test as pom

def load_population(pop_name):

    if pop_name == '2_1_2':
        sim_path = 'C:\Dropbox\Backups\Python\DRG\Simulations\DRG2 - Calcium and beyond'
        pop_filename = 'drg_2_1_2.pickle'
    else:
        raise ValueError('pop_name not found')
        
    pop = pom.load(os.path.join(sim_path, pop_filename))
    return pop