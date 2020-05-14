"""
Analysis functions for pom data
05/09/2018

"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

" Creating dataframes and aggregating population biomarker data per simulation "

def fill_in_grouped_sim_data(grouped_sim_data, results, biomarker,
                             sim_types, sim_type_amp_relationship,
                             sim_name_conversions,
                             how_to_fill, filter_biomarkers_for_outliers=False, threshold=None):
    """"
    Fill in a dataframe with aggregated biomarker data
    """    


    # Check inputs
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

            if threshold:
                # Only keep rows that are not below or at threshold
                unaggregated_data = unaggregated_data[unaggregated_data <= threshold]
                
            # filtering of biomarkers for outliers
            if filter_biomarkers_for_outliers:
                if biomarker == "APHalfWidth":
                    # Remove outliers
                    outlier_definition = get_outlier_definition(biomarker) # ms
                    outliers = unaggregated_data >= outlier_definition
                    num_outliers = outliers.sum()
                    mask = ~outliers # Invert boolean series
                    unaggregated_data = unaggregated_data[mask]
                    #if num_outliers > 0:
                        #print("Removed {} AP Half Width outliers > {}".format(num_outliers, outlier_definition))

            if how_to_fill == 'mean':
                data = unaggregated_data.mean()
            elif how_to_fill == 'std':
                data = unaggregated_data.std()
            elif how_to_fill == 'median':
                data = unaggregated_data.median()
            elif how_to_fill == 'mean_fillna':
                data = unaggregated_data.fillna(0).mean()
            elif how_to_fill == 'mean_freq_fillna':
                assert biomarker == 'ISI', "Biomarker for frequency needs to be ISI not {}".format(biomarker) 
                # Convert to frequency, then fill nans with 0s and take mean
                unaggregated_data = 1000.0/unaggregated_data
                data = unaggregated_data.fillna(0).mean() 
            elif how_to_fill == 'mean_freq_dropna':
                assert biomarker == 'ISI'
                # Convert to frequency, then DROP nans and take mean
                unaggregated_data = 1000.0/unaggregated_data
                data = unaggregated_data.dropna().mean() 
            elif isinstance(how_to_fill, int):
                # Model index
                data = unaggregated_data.loc[how_to_fill]
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

def make_grouped_sim_data(pop, biomarker='APFullWidth', agg='mean', filter_outliers=False, scaled_parameters=['GNav18'], threshold=None):
    " Aggregate population biomarker results per simulation to analyse at a per simulation level over the ensemble population. "

    sim_types = ['step', 'ramp'] # First entry in list is the base
    sim_type_as_base = sim_types[0] # Use step as base as ramp has amplitude x10 of step
    sim_type_amp_relationship = {'step':1, 'ramp':10}
    assert sim_type_amp_relationship[sim_type_as_base] == 1


    grouped_sim_data = make_empty_grouped_sim_data( pop=pop,
                                                    biomarker=biomarker,
                                                    filter_outliers=filter_outliers,
                                                    scaled_parameters=scaled_parameters,
                                                    sim_types=sim_types,
                                                    sim_type_as_base=sim_type_as_base,
                                                    sim_type_amp_relationship=sim_type_amp_relationship,
                                                    )

    sim_name_conversions = make_sim_name_conversions (pop.get_simulation_names(),
                                                       sim_types,
                                                       sim_type_amp_relationship,
                                                       sim_type_as_base
                                                      )

    fill_in_grouped_sim_data(grouped_sim_data, pop.results, biomarker,
                                 sim_types, sim_type_amp_relationship,
                                 sim_name_conversions,
                                 how_to_fill=agg,
                                 filter_biomarkers_for_outliers=filter_outliers,
                                 threshold=threshold)

    if biomarker is not 'Firing pattern':
        grouped_sim_data = grouped_sim_data.astype(float)
    
    return grouped_sim_data

def make_empty_grouped_sim_data(pop, biomarker='APFullWidth', filter_outliers=False,
                                scaled_parameters=['GNav18'], sim_types=['step','ramp'],
                                sim_type_as_base='step', sim_type_amp_relationship={'step':1, 'ramp':10},
                                ):
    " Aggregate population biomarker results per simulation to analyse at a per simulation level over the ensemble population. "
    
    arrays =[[biomarker]*len(sim_types)+['Amp']*len(sim_types) + ['Scaling factors']*len(scaled_parameters),sim_types*2 + scaled_parameters] # Build multiarray columns
    columns = pd.MultiIndex.from_arrays(arrays, names=['',''])
    sim_names =  pop.get_simulation_names()
    sim_name_conversions = make_sim_name_conversions(sim_names, sim_types, sim_type_amp_relationship, sim_type_as_base)
    short_sim_names = sorted(list(sim_name_conversions.keys()))

    grouped_sim_data = pd.DataFrame(columns=columns, index=short_sim_names)
    return grouped_sim_data

def make_sim_name_conversions(sim_names, sim_types, sim_type_amp_relationship, sim_type_as_base):
    " Make conversion dict from short sim names to full sim names with sim type "
    sim_name_conversions = {}
    # Get list of sim_names common to step and ramp
    # Complexity: ramp and step amps are different so we will just get the name from step
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

    return sim_name_conversions

" Data processing "

def process_firing_pattern_data(firing_pattern_percentages, sim_types=['step','ramp'],
                                sim_type_amp_relationship = {'step':1, 'ramp':10}, scaled_parameters=['GNav18'],
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
                                                            scaled_parameters, 
                                                            sim_types, 
                                                            sim_type_amp_relationship,
                                                           )
    return formatted_firing_pattern_data

             

def format_firing_pattern_percentages(firing_pattern_percentages, short_sim_names, sim_name_conversions, scaled_parameters, sim_types, sim_type_amp_relationship):
    """
    Fill in a dataframe with firing pattern percentages for each simulation
    Equivalent to fill_in_grouped_sim_data() but for firing patterns not single numeric biomarkers.
    
    Copy code from fill_in_grouped_sim_data where needed but:
    1. We don't need a how to fill option as we aggregate by percentages always.
    2. We do need to fill in for all firing patterns, not just one biomarker. 
    """
    
    " Create formatted dataframe with column multiindex "
    
    assert len(scaled_parameters) == 1, "Multiple scaled parameters not supported by format_firing_percentages yet"
    firing_pattern_names = firing_pattern_percentages.columns.tolist()
    sim_type_as_base = sim_types[0] # First sim type in list is the base
    assert sim_type_amp_relationship[sim_type_as_base] == 1

    arrays = [list(np.repeat(firing_pattern_names,len(sim_types))) + ['Amp','Amp', 'Scaling factors'], 
              sim_types*(1+len(firing_pattern_names)) + scaled_parameters] # Build multiarray columns
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

def define_region(pop, firing_pattern_thresholds={}, other_thresholds={}, stim_type='step', verbose=False):
    """
    Use some condition e.g. % of models with a given firing pattern to define a region of simulation space to analyse further.
    """

    import operator
    opcodes  = {'>':operator.gt, '>=':operator.ge, '<':operator.lt, '<=':operator.le}

    firing_pattern_percentages = get_firing_pattern_percentages(pop) # Excludes rheobase simulations
    firing_pattern_percentages = process_firing_pattern_data(firing_pattern_percentages)
    region = {}

    # Check each simulation against all thresholds, if it passes them all then add to defined region
    simulations = firing_pattern_percentages.index # Already excluded rheobase simulations and parameters in process_firing_pattern_data
    count = 0
    for simulation in simulations:
        accept_simulation = True

        # Firing pattern thresholds
        for firing_pattern, threshold_vals in firing_pattern_thresholds.items():
            op = opcodes[threshold_vals[0]]
            threshold = threshold_vals[1]
            fp_percentage = firing_pattern_percentages.loc[simulation,(firing_pattern, stim_type)]
            # Check if threshold is accepted
            if not op(fp_percentage, threshold):
                accept_simulation = False

        # Other thresholds (aggregated biomarkers)
        for threshold, val in other_thresholds.items():
            print('Other thresholds (aggregated biomarkers) not implemented yet.')
            # Get aggregated biomarker data and do the same operation on it as for firing pattern but with biomarker values not percentages

        if accept_simulation:
            count += 1
            amplitude = firing_pattern_percentages.loc[simulation,('Amp', stim_type)]
            scaling_factors = dict(firing_pattern_percentages.loc[simulation, 'Scaling factors'])
            region[simulation] = {'Amp':amplitude}
            for scaling_factor, val in scaling_factors.items():
                region[simulation][scaling_factor] = val

    if verbose:
        print('{} simulations out of {} accepted.'.format(count, len(simulations)))
    return region

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

def pom_consistency_in_region(pop, region, firing_pattern, stim_type, amp_stim_relationship):
    """
    Plot or find mean and std dev of the distribution of firing pattern consistency in the population, 
    within a defined region. See brown notebook 1 for a sketch of histogram and cumulative consistency plots.
    """
    
    " Get consistency for each model "
    model_consistencies = pd.DataFrame(index=pop.results.index, columns=['Consistency'])
    for model in pop.results.index:
        model_consistency = model_consistency_in_region(model, pop, region, firing_pattern, stim_type, amp_stim_relationship)
        model_consistencies.at[model, 'Consistency'] = model_consistency
        
    pom_consistency = model_consistencies
    return pom_consistency

def model_consistency_in_region(model, pop, region, firing_pattern, stim_type, amp_stim_relationship):
    """
    Calculate consistency percentage for one model for the given firing pattern, population and simulation region
    """
    
    num_simulations = 0 # Count total number of simulations or could calc directly from region and check we don't get any errors with 
              # any simulation in the region
    num_consistent_simulations = 0
    for simulation in region:
        num_simulations += 1
        sim_full_name = short_name_to_full(simulation, stim_type, amp_stim_relationship, delimiter='_', amp_units='pA')
        simulation_firing_pattern = pop.results.at[model, (sim_full_name, 'Firing pattern')]
        # Check if firing pattern is part of simulation firing pattern

        # Check for nans
        if type(simulation_firing_pattern) != list:
            if np.isnan(simulation_firing_pattern):
                continue

        if firing_pattern in simulation_firing_pattern:
            num_consistent_simulations += 1
        
    model_consistency = 100.*(num_consistent_simulations/num_simulations)
    return model_consistency

def assign_subpopulation_from_region(pop, region, criteria, verbose=False):
    """
    Compute required consistencies and assign subpopulations to a population of models based on results from
    a simulation region.
    Inputs:
    pop - a PopulationOfModels class
    region - a list of simulations
    criteria - has the format: {name:[firing pattern, opcode, value]}
    E.g. {'single < 25':['single', '>=', 25]}
    
    """
    import operator
    opcodes  = {'>':operator.gt, '>=':operator.ge, '<':operator.lt, '<=':operator.le}
    
    # Build consistencies for each criteria 
    consistencies = pd.DataFrame(index=pop.results.index)
    for criterion in criteria:
        
        firing_pattern = criteria[criterion][0]
        pom_consistency = pom_consistency_in_region(pop, 
                                            region, 
                                            firing_pattern=firing_pattern, 
                                            stim_type='step', 
                                            amp_stim_relationship={'step':1,'ramp':10})
        consistencies[firing_pattern] = pom_consistency
          
    # Find the models that fulfill all consistency criteria
    models_passing_criteria = pd.DataFrame(index=consistencies.index)
    for criterion, criterion_params in criteria.items():
        firing_pattern = criterion_params[0]
        opcode = criterion_params[1]
        val = criterion_params[2]
        op = opcodes[opcode]     
        consistency = consistencies[firing_pattern]
        models_passing_criteria[criterion] = op(consistency,val)
        
    models_passing_all_criteria = models_passing_criteria.all(axis=1)
    # Filter away models that don't pass all criteria
    subpopulation = pd.DataFrame(index=pop.results.index)
    subpopulation = subpopulation[models_passing_all_criteria]
    
    if verbose:
        print('{} models out of {} in population of models are in the subpopulation'.format(
            len(subpopulation.index), len(pop.results.index))
             )
        
    return subpopulation

def calculate_percentage_of_firing_patterns(firing_patterns):
    """
    Calculate percentage of each firing pattern.
    firing_patterns is a series or list of lists of firing patterns
    """
    assert type(firing_patterns) == pd.Series, "firing_patterns is not a pandas Series"
    firing_patterns = firing_patterns.dropna() # Clean up data

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

" Data analysis "

def pom_consistency_stats(pom_consistency):
    return {'mean':pom_consistency['Consistency'].mean(), 'std':pom_consistency['Consistency'].std()}

def pom_consistency_hist(pom_consistency):
    plt.hist(pom_consistency)

" Plotting "

def make_3d_plot(data, x, y, z, cutoff, fillna_value, labels, angle=(20,300), zticks=None, title=None):
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

    # Title
    if title is not None:
        ax.set_title(title, fontsize=18, verticalalignment='baseline')

    " Make the plot pop "
    surf = ax.plot_trisurf(df['X'], df['Y'], df['Z'], cmap=plt.cm.viridis, linewidth=0.2)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    ax.set_zlabel(labels[2], fontsize=20, rotation=40)

    if zticks != None:
        ax.set_zticks(zticks)

    # Rotate it
    ax.view_init(angle[0], angle[1])
    plt.tight_layout()

    plt.show()
    return fig, ax

def visualise_region(pop, region, stim_type, scaling_factor='GNav18'):
    " Visualise region as scatter plot with one scaling factor and amp "
    firing_pattern_percentages = get_firing_pattern_percentages(pop) #  Excludes rheobase simulations
    firing_pattern_percentages = process_firing_pattern_data(firing_pattern_percentages)
    
    for sim in firing_pattern_percentages.index:
        amp = firing_pattern_percentages.loc[sim, ('Amp', stim_type)]
        sf_val = firing_pattern_percentages.loc[sim, ('Scaling factors', scaling_factor)]
        plt.scatter(sf_val, amp, color='k')
        
    # Now plot region
    for sim, sim_vals in region.items():
        amp = sim_vals['Amp']
        sf_val = sim_vals[scaling_factor]
        plt.scatter(sf_val, amp, color='r')

def pom_consistency_cumulative_plot(pom_consistency, bins=10):
    " See https://stackoverflow.com/questions/15408371/cumulative-distribution-plots-python for code source"
    data = pom_consistency['Consistency']
    values, base = np.histogram(data, bins=bins)
    cumulative = np.cumsum(values)
    plt.plot(base[:-1], cumulative)

def plot_biomarker_boxplots(pop, biomarker, stim_amps, sim_type, 
                            sim_types=['step', 'ramp'], 
                            sim_type_amp_relationship={'step':1., 'ramp':10.},
                            base_sim_type='step',
                            save=False,
                            save_path=None,
                            fixed_scaling_factors={},
                           ):
    " Plot boxplots of biomarkers for constant stim amps varying scaling factor "
    
    sim_names = pop.get_simulation_names()
    sim_name_conversions = make_sim_name_conversions(sim_names, sim_types, sim_type_amp_relationship, base_sim_type)
    for stim_amp in stim_amps:
        
        simulations = get_simulations_with_fixed_stim_amp(pop, stim_amp, sim_type, sim_type_amp_relationship, base_sim_type, fixed_scaling_factors=fixed_scaling_factors)
 
        # Build dataframe with the simulation and biomarker value of each model
        df = pd.DataFrame()
        for sim in simulations:

            biomarker_data = pd.DataFrame(columns=[biomarker])
            biomarker_data[biomarker] = pop.results[(sim, biomarker)]
            # Add the scaling factor of the simulation
            sf = get_parameter_scaling(sim, 'GNav18', delimiter='_')
            biomarker_data['Scaling factor'] = sf 
            df = pd.concat([df, biomarker_data])

        df[biomarker] = pd.to_numeric(df[biomarker])
        fig = plt.figure(figsize=(7,8))
        sns.set_context("paper", font_scale=2.2)
        sns.boxplot(x='Scaling factor', y=biomarker, data=df, whis=np.inf)
        plt.title(f'Stim amp = {stim_amp} pA.')
        #plt.xticks([0, 5, 10, 15, 20],[0, 0.5, 1.0, 1.5, 2.0])
        #plt.ylim(0,10)
        plt.ylabel('num APs')
        plt.xlabel('GNav 1.8 scaling factor')
        
        plt.tight_layout()

        for i, sf in enumerate(sorted(set(df['Scaling factor']))):
            plt.scatter(i, df[df['Scaling factor'] == sf][biomarker].mean(), s=40, c='r', marker='x')
        if save:
            if not save_path:
                save_path = f'boxplot_{biomarker}_{stim_amp}.png'
            fig.savefig(save_path, dpi=300)

def plot_biomarker_scatterplots(pop, biomarkers, stim_amps, sim_type,
                          sim_types=['step', 'ramp'], 
                          sim_type_amp_relationship={'step':1, 'ramp':10},
                          base_sim_type='step',
                          selected_scaling_factors=None,
                          plot_together=False,
                          plot=True,
                         ):
    " Plot scatter plots of two biomarkers for constant stim amps varying scaling factor "
    assert len(biomarkers) == 2
    
    sim_names = pop.get_simulation_names()
    sim_name_conversions = make_sim_name_conversions(sim_names, sim_types, sim_type_amp_relationship, base_sim_type)
    dfs = {}
    for stim_amp in stim_amps:
        
        simulations = get_simulations_with_fixed_stim_amp(pop, stim_amp, sim_type, sim_type_amp_relationship, base_sim_type)
        
        # Build dataframe with the simulation and biomarker value of each model
        df = pd.DataFrame()
        for sim in simulations:

            biomarker_data = pd.DataFrame(columns=biomarkers)
            for biomarker in biomarkers:
                biomarker_data[biomarker] = pop.results[(sim, biomarker)]
            # Add the scaling factor of the simulation
            sf = get_parameter_scaling(sim, 'GNav18', delimiter='_')
            biomarker_data['Scaling factor'] = sf 
            df = pd.concat([df, biomarker_data])
        
        df[biomarkers] = df[biomarkers].apply(pd.to_numeric)
        dfs[stim_amp] = df
        
        if plot:
            plotted = False
            scaling_factors = sorted(set(df['Scaling factor']))
            for sf in scaling_factors:
                if (selected_scaling_factors is not None) & (sf not in selected_scaling_factors):
                    continue
                df_sf = df[df['Scaling factor'] == sf]
                if plot_together & (not plotted):
                    f, ax = plt.subplots(1,figsize=(6,6))
                    ax.set(xlabel=biomarkers[0], ylabel=biomarkers[1])
                    plotted = True
                elif not plot_together:
                    f, ax = plt.subplots(1,figsize=(6,6))
                    ax.set(title=sf, xlabel=biomarkers[0], ylabel=biomarkers[1]) 

                ax.scatter(df_sf[biomarkers[0]], df_sf[biomarkers[1]])
                print(str(df_sf[biomarkers[0]].mean()) + " " + str(df_sf[biomarkers[1]].mean()))
                #ax.set(xlim=(0,100))

                # @TODO add saving
        
    return dfs

def plot_parameters(pop):
# Plot all the parameters in a population in a grid of scatterplots

    parameters = pop.results['Parameters'].copy()
    parameter_ranges = {'GNav17':[0.0,0.4], 'GNav18':[0.,4.0], 'GNav19':[0.,4.],
                'GKdr':[0.,4.], 'GKA':[0.,40.], 'GKM':[0.,4.], 'GH':[0.,2.], 'GKleak':[0., 0.2]}

    for col in parameters:
        parameters[col] /= parameter_ranges[col][1]/2

    dark_grey = np.array([1,1,1])*(32/255)
    sns.set(style="white")
    sns.set_context(context="notebook", font_scale = 3.5,)
    g = sns.pairplot(parameters, plot_kws = {'color': dark_grey}, diag_kws = {'color': dark_grey})
    # Set scatter plot ranges
    g.set(xticks=[0,2]); g.set(yticks=[0,2])
    # Set histogram ranges
    axes = g.axes
    for i in range(len(parameters.columns)):
        axes[i,i].set_xlim(0,2)

            
" Utilities "

def get_non_rheobase_sim_names(pop):
    sims = [sim_name for sim_name in pop.results.columns.levels[0] if sim_name not in ['Parameters', 'ramp', 'step']]
    return sims
    
def get_amp_str_from_sim_name(name, amp_units, delimiter):
        amp_units_with_delimiter = '{}{}'.format(delimiter, amp_units)
        end = name.find(amp_units_with_delimiter)
        start = name[0:end].rfind(delimiter)+1 # Find from the end of the reduced string
        amp_string = name[start:end]
        return amp_string

def get_amplitude(name, amp_units='pA', delimiter='_'):
    # Get the stim amplitude
    amp_string = get_amp_str_from_sim_name(name, amp_units, delimiter)
    try:
        assert float(amp_string) == int(amp_string)
    except ValueError: # If fail to convert string to float
        return np.nan

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

def get_all_firing_patterns_in_population(pop):
    """
    Get the set of all individual firing patterns found in a population
    """
    sims = get_non_rheobase_sim_names(pop)
    all_firing_patterns = []
    for sim in sims:
        df = pop.results[sim].copy()
        sim_firing_patterns = df['Firing pattern'].dropna() # Drop nans to avoid error in trying to iterate on a nan below
        
        for model_fps in sim_firing_patterns: 
            for fp in model_fps:
                all_firing_patterns.append(fp)
    
    all_firing_patterns = set(all_firing_patterns)
    return all_firing_patterns

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

def short_name_to_full(short_name, stim_type, amp_stim_relationship={'step':1, 'ramp':10}, delimiter='_', amp_units='pA'):
    """
    Convert short name of a simulation without stim type and with amp converted by amp_stim_relationship
    to the full name of the simulation.
    E.g.
    With amp_stim_relationship = {'step':1, 'ramp':10}:
    '2_1_2_400_pA_GNav18_0.1' would become:
    For step - '2_1_2_step_400_pA_GNav18_0.1'
    For ramp - '2_1_2_ramp_4000_pA_GNav18_0.1'
    """
    
    s = short_name

    amp_str = get_amp_str_from_sim_name(s, amp_units, delimiter)
    amp_val = get_amplitude(s, amp_units, delimiter)
    
    parts = s.split(delimiter+amp_str, maxsplit=1) 
    # Need maxsplit to avoid splitting on the scaling factor(s) e.g. for amp=0 and GNav18_0.0
    # Need to add delimiter to avoid splitting on the simulation name (eg simulation_5_10 for amp=0)
    parts[0] += (delimiter + stim_type + delimiter)

    stim_type_amp_str = str(amp_val * int(amp_stim_relationship[stim_type])) # Assume that all amp_stim_relationships are expressed as ints
    parts[0] += stim_type_amp_str

    full_name = ''.join(parts)
    return full_name

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

def get_simulations_with_fixed_stim_amp(pop, stim_amp, sim_type='step', sim_type_amp_relationship={'step':1., 'ramp':10}, base_sim_type='step', sim_types=['step','ramp'], fixed_scaling_factors={}):
    " Get names of all simulations in population with a given stimulus amplitude and type "
    
    sim_names = pop.get_simulation_names()
    sim_name_conversions = make_sim_name_conversions(sim_names, sim_types, sim_type_amp_relationship, base_sim_type)
    grouped_sim_data = make_empty_grouped_sim_data(pop, biomarker='')

    # Add metadata
    scaling_factors = list(grouped_sim_data['Scaling factors'].columns)
    for name in grouped_sim_data.index:

        simulations = sim_name_conversions[name]
        for simulation, _stim_type in simulations.items():

            amp = get_amplitude(simulation, amp_units='pA', delimiter='_')
            grouped_sim_data.at[name, ('Amp', _stim_type)] = amp

            for scaling_factor in scaling_factors:
                scaling_factor_value = get_parameter_scaling(simulation, scaling_factor, delimiter='_')
                grouped_sim_data.at[name, ('Scaling factors', scaling_factor)] = scaling_factor_value


    # Filter the dataframe to find the simulations we want
    simulations = list(grouped_sim_data[grouped_sim_data[('Amp', sim_type)] == stim_amp].index)
    simulations = [short_name_to_full(sim, sim_type, amp_stim_relationship=sim_type_amp_relationship) for sim in simulations]

    # Filter simulations to only include those with particular values of a scaling factor, if required
    for fixed_sf, val in fixed_scaling_factors.items():
        simulations = [sim for sim in simulations if get_parameter_scaling(sim, fixed_sf, delimiter='_') == val] 

    return simulations

" Constants "


def get_outlier_definition(biomarker):
    """
    Centralised definitions of biomarker outliers for filtering
    """

    if biomarker == "APHalfWidth":
        outlier_definition = 100.0 # ms
    else:
        raise ValueError(f"Biomarker {biomarker} not found.")
    return outlier_definition
            
def get_default_stim_type_amp_relationship(): 
    sim_type_amp_relationship={'step':1, 'ramp':10}
    return sim_type_amp_relationship
