"""
Analysis functions for pom data
05/09/2018

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

"""
Main features available:

Visualisation of firing pattern percentages
---
* process_sim_name 
*
*

--
3D plot workflow
--

1. Call get_firing_pattern_percentages on a PopulationOfModels to get the firing pattern percentages for each simulation.

2. Call process_firing_pattern_data() to convert the firing pattern dataframe into a hierarchical column structure and process the simulation 
names to extract simulation parameters.

3. Use the pivot_dataframe_for_heatmap(df, x, y, z, labels) to grab the required data for a surface plot.

4. Call the 3D surface plotting function make_3d_plot to plot the data.
"""


" --- Visualisation of firing pattern percentage in a population across multiple simulations --- "
    
def get_non_rheobase_sim_names(pop):
    sims = [sim_name for sim_name in pop.results.columns.levels[0] if sim_name not in ['Parameters', 'ramp', 'step']]
    return sims
    
def how_to_fill_options():
    " Returns the available options for aggregating a quantitaitve biomarker over a population of models. "
    return ['mean', 'std', 'median']

def get_amp_str_from_sim_name(name, amp_units, delimiter):
        amp_units_with_delimiter = '{}{}'.format(delimiter, amp_units)
        end = name.find(amp_units_with_delimiter)
        start = name[0:end].rfind(delimiter)+1 # Find from the end of the reduced string
        amp_string = name[start:end]
        return amp_string

def get_amplitude(name, amp_units='pA', delimiter='_'):
    # Get the stim amplitude
    amp_string = get_amp_str_from_sim_name(name, amp_units, delimiter)
    assert float(amp_string) == int(amp_string)
    amp = int(amp_string)
    return amp
    
def get_parameter_scaling(name, parameter_name, delimiter='_'):
    start = name.find(parameter_name)+len(parameter_name)+len(delimiter)
    end = name[start:].find(delimiter)
    if end == -1:
        parameter_string = name[start:]
    else:
        parameter_string = name[start:start+end]
    parameter_scaling_factor = float(parameter_string)
    return parameter_scaling_factor
    
def fill_in_grouped_sim_data(grouped_sim_data, results, biomarker,
                             sim_types, sim_type_amp_relationship,
                             sim_name_conversions,
                             how_to_fill):
    
    # Fill in the dataframe we've created


    
    # Check inputs
    assert how_to_fill in how_to_fill_options(), "Method: {} for filling in biomarker data not supported".format(how_to_fill)
    assert any([amp_scale_factor == 1 for _,amp_scale_factor in sim_type_amp_relationship.items()]), "No base value of 1 for any of the sim types in sim_type_amp_relationship"
        
    for name in grouped_sim_data.index:
        # To do here: get simulations from the conversion dict, do the appropriate averaging for the how_to_fill method, add that data in
        # Then, get gnav18 from the simulation short name and the amplitudes from the long names using the conversion dict.
        # Add these in too.
        # And we're done
        #names = 
        
        # Get simulations
        simulations = sim_name_conversions[name]
        for simulation, sim_type in simulations.items():
            # Input biomarker data
            unaggregated_data = results.loc[:,(simulation,biomarker)]
            
            if how_to_fill == 'mean':
                data = unaggregated_data.mean()
            elif how_to_fill == 'std':
                data = unaggregated_data.std()
            elif how_to_fill == 'median':
                data = unaggregated_data.median()
            else:
                raise ValueError("How to fill method: {} not supported.".format(how_to_fill))
            
            grouped_sim_data.at[name,(biomarker, sim_type)] = data

            # Input amplitudes
            amp = get_amplitude(simulation, amp_units='pA', delimiter='_')
            grouped_sim_data.at[name, ('Amp', sim_type)] = amp
            
            # Input scaling factors
            scaling_factors = list(grouped_sim_data['Scaling factors'].columns)
            
            for scaling_factor in scaling_factors:
                # Get scaling factor 
                scaling_factor_value = get_parameter_scaling(simulation, scaling_factor, delimiter='_')
                grouped_sim_data.at[name, ('Scaling factors', scaling_factor)] = scaling_factor_value
             

" Aggregating population biomarker data per simulation "
            
def make_grouped_sim_data(pop, biomarker='APFullWidth', agg='mean'):
    " Aggregate population biomarker results per simulation to analyse at a per simulation level over the ensemble population. "
    scaled_parameter = 'GNav18'
    sim_types = ['step', 'ramp'] # First entry in list is the base
    sim_type_as_base = sim_types[0] # Use step as base as ramp has amplitude x10 of step
    sim_type_amp_relationship = {'step':1, 'ramp':10}
    assert sim_type_amp_relationship[sim_type_as_base] == 1
    how_to_fill = agg


    arrays = [[biomarker]*len(sim_types )+['Amp','Amp', 'Scaling factors'],sim_types*2 + [scaled_parameter]] # Build multiarray columns
    columns = pd.MultiIndex.from_arrays(arrays, names=['',''])

    # Get list of sim_names common to step and ramp
    # Complexity: ramp and step amps are different so we will just get the name from step
    sim_names =  list(pop.simulations.keys())
    short_sim_names = []
    sim_name_conversions = {}
    for sim_type in sim_types:
        for sim_name in sim_names:
            #  # Ignore rheobase simulations and do base sim type first - important to get oreder right in sim_name_conversions
            if (sim_name not in sim_types) & (sim_type in sim_name):
                short_sim_name, _sim_type = process_sim_name(sim_name, sim_types, sim_type_amp_relationship, amp_units='pA', delimiter='_')
                assert _sim_type == sim_type

                # Build up conversion dict from short name to full names and store sim_type
                if short_sim_name in sim_name_conversions.keys():
                    sim_name_conversions[short_sim_name][sim_name] = sim_type
                else:
                    assert sim_type == sim_type_as_base, "Sim type: {}, name:{}, short_sim_name:{}".format(sim_type, sim_name, short_sim_name)
                    sim_name_conversions[short_sim_name] = {sim_name:sim_type}

                if sim_type == sim_type_as_base: # Only add step to short_sim_names to avoid adding ramp as ramp has different amplitude
                    short_sim_names.append(short_sim_name)

    short_sim_names = sorted(short_sim_names)
    assert len(short_sim_names) == len(set(short_sim_names)), "Check for repeated simulation names failed."

    # Build up grouped sim data
    grouped_sim_data = pd.DataFrame(columns=columns, index=short_sim_names)
    fill_in_grouped_sim_data(grouped_sim_data, pop.results, biomarker,
                                 sim_types, sim_type_amp_relationship,
                                 sim_name_conversions,
                                 how_to_fill=how_to_fill)


    if biomarker is not 'Firing pattern':
        grouped_sim_data = grouped_sim_data.astype(float)
    
    return grouped_sim_data


" ---Firing pattern data processing--- "

def get_firing_pattern_percentages(pop):
    """
    Get the firing pattern percentages for each simulation and each firing pattern present in the
    population of models provided as input.
    """
    
    # Initialise firing pattern percentages dataframe
    sims = get_non_rheobase_sim_names(pop)
    all_firing_patterns = get_all_firing_patterns_in_population(pop)
    firing_pattern_percentages = pd.DataFrame(index=sims, columns=all_firing_patterns)

    # Fill in percentages
    for sim in sims:
        firing_patterns = pop.results[sim]['Firing pattern'].copy()
        percentages = calculate_percentage_of_firing_patterns(firing_patterns)

        # Save percentage for each firing pattern
        for fp, perc in percentages.items():
            firing_pattern_percentages.at[sim,fp] = perc

    # Here we could fill in all the NaNs to 0
    firing_pattern_percentages = firing_pattern_percentages.fillna(0)
    return firing_pattern_percentages

def process_firing_pattern_data(firing_pattern_percentages, sim_types=['step','ramp'],
                               sim_type_amp_relationship = {'step':1, 'ramp':10}, scaled_parameter='GNav18',
                               ):
    """
    Process simulation names to extract simulation parameters and rename to remove stimulation protocol from name
    Flow:
    1. Convert sim names by removing sim type and storing the conversion between full and shortened names.
    2. Create a formatted dataframe and fill in the firing pattern percentages from the values in the original dataframe.
    3. Extract simulation parameters from the simulation name and add to the formatted dataframe.
    
    TODO: This code shares a lot of code with the functions for aggregating biomarker data. Could refactor into one set of functions
    sharing common code.
    
    
    """
    
    # TODO: Could turn these lines creating sim_names and sim_name_conversions 
    # into a function shared with similar code for biomarkers
    sim_type_as_base = sim_types[0]
    sim_names = firing_pattern_percentages.index.tolist()
    short_sim_names = [] # Simulation names without stimulation protocol
    sim_name_conversions = {} # Conversion between full and short sim names
    
    for sim_type in sim_types:
        for sim_name in sim_names:
            
            if (sim_name not in sim_types) & (sim_type in sim_name): # Remove rheobase simulations
                short_sim_name, _sim_type = process_sim_name(sim_name, sim_types, sim_type_amp_relationship, amp_units='pA', delimiter='_')
                assert _sim_type == sim_type
                
                # Create conversion between names
                if short_sim_name in sim_name_conversions.keys():
                    sim_name_conversions[short_sim_name][sim_name] = sim_type
                else:
                    assert sim_type == sim_type_as_base, (
                        "Sim type: {}, name:{}, short_sim_name:{} is not the base sim type for the sim type amp relationship.".format(
                            sim_type, sim_name, short_sim_name))
                    sim_name_conversions[short_sim_name] = {sim_name:sim_type}

                if sim_type == sim_type_as_base: # Only add step to sim_names to avoid adding ramp as ramp has different amplitude
                    short_sim_names.append(short_sim_name)
                    
    short_sim_names = sorted(short_sim_names)
    formatted_firing_pattern_data = format_firing_pattern_percentages(
                                                            firing_pattern_percentages, 
                                                            short_sim_names,
                                                            sim_name_conversions, 
                                                            scaled_parameter, 
                                                            sim_types, 
                                                            sim_type_amp_relationship,
                                                           )
    
    return formatted_firing_pattern_data

def format_firing_pattern_percentages(firing_pattern_percentages, short_sim_names, sim_name_conversions, scaled_parameter, sim_types, sim_type_amp_relationship):
    """
    Fill in a dataframe with firing pattern percentages for each simulation
    Equivalent to fill_in_grouped_sim_data() but for firing patterns not single numeric biomarkers.
    
    Copy code from fill_in_grouped_sim_data where needed but:
    1. We don't need a how to fill option as we aggregate by percentages always.
    2. We do need to fill in for all firing patterns, not just one biomarker. 
    """
    
    " Create formatted dataframe with column multiindex "
    
    firing_pattern_names = firing_pattern_percentages.columns.tolist()
    sim_type_as_base = sim_types[0] # First sim type in list is the base
    assert sim_type_amp_relationship[sim_type_as_base] == 1

    arrays = [list(np.repeat(firing_pattern_names,len(sim_types))) + ['Amp','Amp', 'Scaling factors'], 
              sim_types*(1+len(firing_pattern_names)) + [scaled_parameter]] # Build multiarray columns
    columns = pd.MultiIndex.from_arrays(arrays, names=['',''])

    formatted_firing_pattern_data = pd.DataFrame(index=short_sim_names, columns=columns)

    " Fill in firing pattern percentages; get and fill in simulation parameters from simulation names "
    for name in formatted_firing_pattern_data.index:
        simulations = sim_name_conversions[name]
        for simulation, sim_type in simulations.items():
            
            # Fill in firing pattern percentage
            for fp_name in firing_pattern_names:
                formatted_firing_pattern_data.at[name, (fp_name, sim_type)] = firing_pattern_percentages.at[simulation,fp_name]
        
            # Fill in stimulus amplitude
            amp = get_amplitude(simulation, amp_units='pA', delimiter='_')
            formatted_firing_pattern_data.at[name, ('Amp', sim_type)] = amp

            # Fill in parameter scaling
            scaling_factors = list(formatted_firing_pattern_data['Scaling factors'].columns)
            for scaling_factor in scaling_factors:
                scaling_factor_value = get_parameter_scaling(simulation, scaling_factor, delimiter='_')
                formatted_firing_pattern_data.at[name, ('Scaling factors', scaling_factor)] = scaling_factor_value
    
    return formatted_firing_pattern_data


" --- 3D visualisation--- "

def make_3d_plot(data, x, y, z, cutoff, fillna_value, labels, angle=(20,300)):
    """
    Construct a 3d plot 
    """
    df = pd.DataFrame(columns=["X","Y","Z"])
    df['X'] = data[x].copy()
    df['Y'] = data[y].copy()
    df['Z'] = data[z].copy()

    df = df.fillna(value=fillna_value,axis=1)
    df.loc[df['Z'] >= cutoff, 'Z'] = fillna_value * 3 # Change cutoff points to be very low
    #df['Z'][df['Z'] >= cutoff] = fillna_value * 3 # old redundant code setting on slice

    # Make the plot
    fig = plt.figure(figsize=(10,10))
    ax = fig.gca(projection='3d')

    " Make the plot "
    surf = ax.plot_trisurf(df['X'], df['Y'], df['Z'], cmap=plt.cm.viridis, linewidth=0.2)
    " Color bar "

    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    ax.set_zlabel(labels[2], fontsize=20, rotation=40)

    " Rotate "
    # Rotate it
    ax.view_init(angle[0], angle[1])
    plt.show()
    return fig

def pivot_dataframe_for_heatmap(df, x, y, z, labels):
    """ 
    Pivot a dataframe and drop unwanted components so the the columns are x, the rows are y, and the element values are z. 
    Labels should be axis names ordered in a list [x,y,z].
    """
    df_labels = [x,y,z]
    df = pd.DataFrame()
    
    for df_label, heatmap_label in zip(df_labels,labels):
        df[heatmap_label] = grouped_sim_data[df_label]

    # Pivot
    df = df.pivot(index=labels[1], columns=labels[0], values=labels[2])
    df = df.iloc[::-1] # Reverse index to have origin at 0 0
    return df
    

" --- Utility --- "
    
def process_sim_name(name, sim_types, sim_type_amp_relationship, amp_units='pA', delimiter='_'):
    """
    Helper function for building an index with each row being a combined step and ramp result.
    Remove reference to sim_types and trailing underscores. Can just give one sim_type if you want to avoid problem
    with different sim_types having different amplitudes and not grouping.
    """
    for sim_type in sim_types:
        if sim_type in name:
            # Remove sim type
            new_name = name.replace(sim_type + '_','')
            
            # Convert amplitude according to sim_type - amp relationship
            # by finding text before the amp unit and converting
            
            # Find amplitude
            amp_string = get_amplitude(new_name, amp_units, delimiter)
            amp = float(amp_string)
            
            # Convert amplitude by relationship
            amp_scale_down = sim_type_amp_relationship[sim_type]
            # Warning: this may break if the scale down is something that isn't a factor of the amplitude
            new_amp = amp/amp_scale_down
            assert int(new_amp) == float(new_amp), "Something going wrong with amp: {}, check process_sim_name().".format(amp)
            new_amp_string = str(int(new_amp)) # Convert to int since float == int 
            
            # Put the new_amp_string back into the new name
            new_name = new_name.replace('{}{}{}'.format(amp_string, delimiter, amp_units),
                                        '{}{}{}'.format(new_amp_string, delimiter, amp_units))
            return new_name, sim_type 
        
    return None

def get_non_rheobase_sim_names(pop):
    sims = [sim_name for sim_name in pop.results.columns.levels[0] if sim_name not in ['Parameters', 'ramp', 'step']]
    return sims
 
def how_to_fill_options():
    "Different options for filling in a summary biomarker for a whole population."
    return ['mean', 'std', 'median']

def get_amp_str_from_sim_name(name, amp_units, delimiter):
        amp_units_with_delimiter = '{}{}'.format(delimiter, amp_units)
        end = name.find(amp_units_with_delimiter)
        start = name[0:end].rfind(delimiter)+1 # Find from the end of the reduced string
        amp_string = name[start:end]
        return amp_string

def get_amplitude(name, amp_units='pA', delimiter='_'):
    " Get the stimulus amplitude from a simulation name "
    amp_string = get_amp_str_from_sim_name(name, amp_units, delimiter)
    assert float(amp_string) == int(amp_string)
    amp = int(amp_string)
    return amp
    
def get_parameter_scaling(name, parameter_name, delimiter='_'):
    start = name.find(parameter_name)+len(parameter_name)+len(delimiter)
    end = name[start:].find(delimiter)
    if end == -1:
        parameter_string = name[start:]
    else:
        parameter_string = name[start:start+end]
    parameter_scaling_factor = float(parameter_string)
    return parameter_scaling_factor

def calculate_percentage_of_firing_patterns(firing_patterns):
    """
    Calculate percentage of each firing pattern.
    firing_patterns is a series or list of lists of firing patterns
    """
    assert type(firing_patterns) == pd.Series, "firing_patterns is not a pandas Series"
    
    firing_patterns_in_this_sim = []
    for model_fps in firing_patterns: # Iterate over models
        for fp in model_fps: # Iterate over firing patterns for one model
            firing_patterns_in_this_sim.append(fp)

    set_of_firing_patterns_in_this_sim = set(firing_patterns_in_this_sim)

    # Then iterate through and count each one
    counts = {fp:0 for fp in set_of_firing_patterns_in_this_sim}
    for model_fps in firing_patterns:
        for fp in model_fps:
            counts[fp] += 1 
   
    # Then divide each by number of models
    percentages = {fp:None for fp in set_of_firing_patterns_in_this_sim}
    num_models = len(firing_patterns)
    for fp, val in counts.items():
        percentages[fp] = 100.0*float(val)/num_models
    return percentages

def get_all_firing_patterns_in_population(pop):
    """
    Get the set of all individual firing patterns found in a population
    """
    sims = get_non_rheobase_sim_names(pop)
    all_firing_patterns = []
    for sim in sims:
        df = pop.results[sim].copy()
        sim_firing_patterns = df['Firing pattern']
        
        for model_fps in sim_firing_patterns: 
            for fp in model_fps:
                all_firing_patterns.append(fp)
    
    all_firing_patterns = set(all_firing_patterns)
    return all_firing_patterns


" --- Reshape pom results --- "

def get_biomarker_values_by_simulation(results, biomarker_names, simulation_names):
    """
    Transform pom results dataframe to turn make the simulation name a categorical column rather than a multiindex level
    
    Examples
    1 - simple:
    pop = pom.load('nav18_project_2_1_1.pickle')
    get_biomarker_values_by_simulation(pop.results, biomarker_names='APPeak', simulation_names=['ramp','step'])

    2 - multiple biomarkers, list comprehension for simulation names:
    get_biomarker_values_by_simulation(pop.results, 
                                   biomarker_names=['APPeak', 'Threshold', 'Firing pattern'], 
                                   simulation_names=[name for name in pop.results.columns.levels[0] 
                                                     if name not in ['Parameters', 'step', 'ramp']]
                                  )
    """
    data = results.copy()
    data = results.swaplevel(0,1,axis=1) # Swap so top column level is biomarker not simulation
    data = data.sort_index(level=0,axis=1) # Group under each biomarker  
    data = data.loc[:,(biomarker_names, simulation_names)] # Get only needed biomarker and simulations
    data = data.stack().reset_index().set_index('Model') # Remove simulation from column index and turn it into column
    data = data.rename(columns={'':'Simulation'}) # Name the Simulation column
    return data

