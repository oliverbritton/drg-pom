# neuron_biomarkers.py
# calculation of AP biomarkers from neuronal voltage traces

import sys
import numpy as np
import pandas as pd
from scipy import optimize
from matplotlib import pyplot as plt

from . import davidson_biomarkers as db
from .. import simulation_helpers as sh
from .. import analysis as an

# Biomarkers to manage and analyse neuronal simulation data and potentially experimental
# data too

RHEO_FAIL = np.nan # Code to return if rheobase calculation fails.
# np.nan == np.nan returns False so use is not instead. pom code
# relies on this value being nan to interface with pandas correctly.

def calculate_biomarkers(traces, model):
    " Calculate every biomarker and output to dict "
    " TODO: Use the rheobase to work out what simulation to run to calculate biomarkers "
    " off of (at rheobase) "
    # biomarker_names = ['APFullWidth', 'APPeak', 'APRiseTime', 'APSlopeMin', 'APSlopeMax',. 'AHPAmp', 'AHPTau', 'ISI', 'RMP', 'Rheobase']

    biomarkers = calculate_simple_biomarkers(traces, model)
    biomarkers['RMP'] =  np.mean(calculate_rmp(traces))
    # Need to do rheobase separately
    biomarkers['Rheobase'] =  calculate_rheobase(model, amp_step=0.1, amp_max=5, make_plot=False,)
    
    return biomarkers
    
def average_biomarker_values(biomarkers, how_to_handle_nans='return'):
    " Average biomarker values for multiple APs while handling biomarkers that "
    if how_to_handle_nans == 'return': # Advantage is we return nan if there are any nans - good for calibration and trouble shooting - shows up weird models easily. 
        pass
    elif how_to_handle_nans == 'remove': # Risky option. Advantage is we still get a number back in mixed cases of nan and non-nan biomarkers, which is potentially risky as it hides a problem in one or more APs.
        biomarkers = np.array(biomarkers)
        if biomarkers[~np.isnan(biomarkers)].size == 0:
            return np.nan
        else:
            biomarkers = biomarkers[~np.isnan(biomarkers)]
    else:
        raise ValueError("{} is not an accepted nan handling method.".format(how_to_handle_nans))
        
    mean_result = np.mean(biomarkers)
    return mean_result
    
def calculate_simple_biomarkers(traces, model="Not needed", how_to_handle_nans='return'):
    """ Calculate every biomarker that can be calculated from a normal simulation trace and output to dict - rheobase and RMP need to be calculated separately."""
    biomarkers = {}
    # biomarker_names = ['APFullWidth', 'APPeak', 'APRiseTime', 'APSlopeMin', 'APSlopeMax',. 'AHPAmp', 'AHPTau', 'ISI', 'RMP', 'Rheobase']
    def error_handle(filename, traces): # Error handler for finding out why biomarkers are throwing errors.
        import pickle
        print(sys.exc_info())
        print(traces['numAPs'])
        plt.figure()
        for t,v in zip(traces['t'],traces['v']):
            plt.plot(t,v)
        with open(filename, 'wb') as handle:
            pickle.dump(traces, handle)
        print("Error, traces dumped to {}.".format(filename))

    try:
        biomarkers['APFullWidth'] = average_biomarker_values(calculate_ap_full_width(traces,threshold=5.,method='gradient'), how_to_handle_nans)
    except:
        error_handle('fullwidth.pickle',traces)
    
    try:
        biomarkers['APHalfWidth'] = average_biomarker_values(calculate_ap_half_width(traces,threshold=5.,method='gradient'), how_to_handle_nans)
    except:
        error_handle('halfwidth.pickle',traces)
                   
    biomarkers['APPeak'] = average_biomarker_values(calculate_ap_peak(traces),how_to_handle_nans)
    
    try:
        biomarkers['APRiseTime'] = average_biomarker_values(calculate_ap_rise_time(traces,dvdtthreshold=5),how_to_handle_nans)
    except:
        error_handle('risetime.pickle',traces)
    ap_slope_mins, ap_slope_maxs = calculate_ap_slope_min_max(traces)
    biomarkers['APSlopeMin'] = average_biomarker_values(ap_slope_mins, how_to_handle_nans)
    biomarkers['APSlopeMax'] = average_biomarker_values(ap_slope_maxs, how_to_handle_nans)
    biomarkers['Threshold'] = average_biomarker_values(calculate_threshold(traces), how_to_handle_nans)
    
    amp, tau, trough = fit_afterhyperpolarization(traces=traces,dvdt_threshold=5, ahp_model='single_exp', full_output=False)
    """
    try:
        amp, tau = fit_afterhyperpolarization(traces=traces,dvdt_threshold=5, ahp_model='single_exp', full_output=False)
    except:
        error_handle('fitahp.pickle',traces)
        amp=0
        tau=0
    """
    biomarkers['AHPAmp'] =  amp
    biomarkers['AHPTau'] =  tau
    biomarkers['AHPTrough'] = trough
    biomarkers['ISI'] = inter_spike_interval(traces)
    biomarkers['numAPs'] = traces['numAPs']
        
    return biomarkers
    
def compute_model_biomarkers(model=None, mechanisms=None, make_plot=False, sim_kwargs=None, xlims=None):
    " Find all standard biomarkers of a model or mechanism set. "

    biomarkers = {}

    if model == None:
        model = sh.build_model(mechanisms)
    # Else use model
        
    if sim_kwargs:
        sim_kwargs['model'] = model
    else:
        sim_kwargs = sh.get_default_simulation_kwargs(model=model)

    rheobase = calculate_rheobase(model, amp_step=0.1, amp_max=5, make_plot=False, sim_kwargs = sim_kwargs)

    if (isinstance(rheobase,float) == False) & (isinstance(rheobase,int) == False):
        # Rheobase not found, don't calculate other biomarkers
        find_other_biomarkers = False
    else:
        find_other_biomarkers = True

    if sim_kwargs:
        sim_kwargs['amp'] = rheobase
        sim_kwargs['model'] = model
    else: 
        sim_kwargs = sh.get_default_simulation_kwargs(amp=rheobase, model=model)
    sim_kwargs['make_plot'] = make_plot

    if find_other_biomarkers:
        output = sh.simulation(**sim_kwargs) 
        t = output['t']; v = output['v']
        t = t[::2]; v = v[::2] # 20 kHz
        traces = split_trace_into_aps(t,v)
        biomarkers = calculate_simple_biomarkers(traces,model,how_to_handle_nans='remove')

    # RMP
    rmp_kwargs = {'amp':0.0, 'dur':3000., 'delay':0., 'interval':0., 'num_stims':1, 't_stop':3000.}
    for kwarg in sim_kwargs:
        # Write in sim_kwargs where they are not already present in rmp_kwargs
        # so that non-RMP specific kwargs are consistent between simulations
        if kwarg not in rmp_kwargs:
            rmp_kwargs[kwarg] = sim_kwargs[kwarg]

    output = sh.simulation(**rmp_kwargs)
    rmp_t = output['t']; rmp_v = output['v']
    rmp_t = rmp_t[::2]; rmp_v = rmp_v[::2] # 20 kHz
    rmp_traces = split_trace_into_aps(rmp_t,rmp_v)
    rmp = np.mean(calculate_rmp(rmp_traces))

    if (make_plot & (xlims != None)):
        plt.xlim(xlims[0], xlims[1])

    # If we calculated other biomarkers, add extras calculated in separate simulations.
    # If we didn't add to empty dictionary, will leave nans when added to master dataframe
    # which is what we want.
    biomarkers['Rheobase'] = rheobase
    biomarkers['RMP'] = rmp
    return biomarkers
    
" --- Calculation and trace manipulation functions -- "

def split_trace_into_aps(t,v,threshold=0,time_threshold=5, check_voltage_gradient=True):#
    """
    Threshold is at 0 mV which can let RF to cause spurious AP detection unless
    we perform a voltage gradient check, which defaults to True.
    
    -- Old ideas to solve the spurious AP detection problem with threshold at 0 mV --
    One idea is to do split trace and then calculate AP width using a voltage threshold of something like -25 mV.
    Then, if AP width is really long (> 100 ms?), redo the calculation with a lower threshold (0 mV?). 
    If mean AP width is then < 100 ms, we use the new split. We could write a log file to say that this happened,
    with the trace in it. 

    However, that is complex and may break if something comes up I haven't thought of.
    Instead we could reset the default threshold to 0 mV but add in a gradient check on the voltage crossing from below.
    Currently a gradient threshold of 1 mV/ms seems like it should be effective although I don't have any examples of slow
    calcium initated APs to test against. 
    """
	
    # Units for defaults
    # t, time_threshold - ms
    # v, threshold - mV

    assert len(t) == len(v), "v and t length mismatch"

    crossings = []
    time_crossings = np.array([])
    
    # Looks for crossings from below
    for i,voltage in enumerate(v[:-1]):
        if (voltage < threshold) & (v[i+1] >= threshold):
            # Check local voltage gradient if neeeded, if gradient is too small ignore the crossing
            # Time window set to 1.0 to try to counteract bug with averaging too much of the pre-upstroke. 
            if (check_voltage_gradient) & (is_voltage_gradient_too_small(i, t, v, dvdt_threshold=1.0, time_window=1.0)):
                continue # Don't add the crossing if the local voltage gradient is small and we're checking for that
            crossings.append(i)
            time_crossings = np.append(time_crossings,t[i])
                
    # For each crossing, remove all instances within the time threshold, leaving only the first crossing of the threshold 
    grouped_crossings = np.zeros(np.size(crossings),float) 
    
    for i in range(len(crossings)-1):
        if grouped_crossings[i] == 0:
            nearby_crossings = np.array( (time_crossings[i+1:] - time_crossings[i]) < time_threshold )
            # Assign 
            grouped_crossings[i+1:] += nearby_crossings
            assert all(grouped_crossings < 2), "Grouped crossing grouped more than once"
            
    firstCrossIndices = np.where(grouped_crossings == 0)
    # Need to turn crossings into a numpy array to index it with np.where
    firstCrossings = np.array(crossings)[firstCrossIndices]
    numAPs = len(firstCrossings)
    assert numAPs >= 0, "Negative number of APs!"
    
    # Assign time and voltage to traces
    times = []
    voltages = []
    # If 1 or 0 APs, return 1 trace, otherwise...
#    if (numAPs == 0) | (numAPs == 1):
#        times.append(t)
#        voltages.append(v)

    
    """ 
    There are some commented assumptions about where traces begin and end here. The core idea is that all data points in the trace have to be assigned to 1 and only 1 AP. If areas of quiescence are a problem for    particular analysis methods, they will be stripped out by other specialised functions. 
    Our goal in this function is to divide up the trace without leaving any of it out, so that we have everything for any future analysis.
    """
    # If we have multiple APs, for each AP find the minimum value
    # of Vm before the next AP
    if numAPs > 0:
        startIdx = np.zeros(numAPs,int)
        endIdx = np.zeros(numAPs,int)    
        for AP in range(numAPs):
            if AP == 0:
                startIdx[0] = 0 # Start of first AP is beginning of trace
            else:
                startIdx[AP] = endIdx[AP-1]+1 # Start of all other APs is after last AP
        
            if AP == numAPs-1:
                endIdx[AP] = len(v)-1 # End of last AP is end of trace
            else:
                # Calculate end of this trace - end is minimum voltage of this trace
                # From threshold of this AP to just before beginning of next threshold
                voltageDuringCurrentAP = v[firstCrossings[AP]:firstCrossings[AP+1]]
                
                # Get index of minimum voltage AFTER the peak
                max_idx = np.argmax(voltageDuringCurrentAP)
                minVmIdx = np.argmin(voltageDuringCurrentAP[max_idx:])
                endIdx[AP] = firstCrossings[AP] + max_idx + minVmIdx # Don't think I need to minus 1 because Python indices start at 0
                
            times.append(t[startIdx[AP]:endIdx[AP]+1])
            voltages.append(v[startIdx[AP]:endIdx[AP]+1]) # Add 1 to as Python slicing ends 1 before last index
    
        for i in range(len(startIdx)-1):
            assert endIdx[i]+1 == startIdx[i+1], "startIdx and endIdx don't match up."
    # Case for no APs - numAPs causes problems here so set indices manually
    elif numAPs == 0:
        times.append(t)
        voltages.append(v)
        startIdx = np.array([0],int)
        endIdx = np.array([len(v)-1],int)
        
    
    assert startIdx[0] == 0, "First AP doesn't start at beginning of trace."
    assert endIdx[-1] == len(v)-1, "Last AP doesn't end at end of trace."
    
    return{'t':times, 'v':voltages, 'startIndices':startIdx, 'endIndices':endIdx, 'numAPs':numAPs}   
SplitTraceIntoAPs = split_trace_into_aps # Alias
    
def voltage_gradient(t,v, method='gradient'):
    # There is a gradient function in numpy to take central differences
    if method == 'gradient':
        dvdt = np.gradient(v,t)# Central differences except at end points
    elif method == 'diff':
        dvdt = np.diff(v)/np.diff(t) # Difference between adjacent points 
    else :
        raise ValueError("Method not found.")
    return dvdt
VoltageGradient = voltage_gradient # Alias


def is_voltage_gradient_too_small(i, t, v, dvdt_threshold, time_window):
    """
    Check if the voltage gradient around the threshold crossing from below at v[i] to v[i+1]
    is too small to have a reasonable likelihood of being a real AP.
    Inputs:
    i - index at which threshold is crossed, between v[i] and v[i+1]
    t,v - time and voltage arrays
    dvdt threshold - mV/ms
    time window (either side of indices i and i+1, so effective window is double the size) - ms
    """ 
    voltage_gradient_too_small = False
    
    # Get time window around v[i] and v[i+1]
    lower_t_bound = t[i] - time_window
    if lower_t_bound < 0: lower_t_bound = 0
    upper_t_bound = t[i+1] + time_window

    # Get indices of t and v that are within the window
    t = np.array(t)
    window_indices = (t >= lower_t_bound) & (t <= upper_t_bound)
    _t = t[window_indices]

    _v = v[window_indices]
    _dvdt = np.gradient(_v,_t)

    # Check mean gradient against threshold
    if np.mean(_dvdt) < dvdt_threshold:
        voltage_gradient_too_small = True

    return voltage_gradient_too_small

    
# --- Biomarkers ---
def rmp(v):
    # RMP should be calculated from a quiescent trace (no stimulus)
    # Ignore first 90% of trace to remove artifacts
    vLen = len(v)
    startIdx = 90*vLen//100
    RMP = min(v[startIdx:])
    RMPIdx = np.argmin(v[startIdx:]) 
    return RMP, RMPIdx
RMP = rmp # Alias
    
def input_res(t, v, current_injection_time):
    # Input resistance calculated from a protocol with an equilibration phase
    # to get to RMP, followed by a sustained small current input.
    # Input res then = (v[-1] - v[RMP])/(I-0)
    # In Davidson, 50 to 100 pA current pulse was used to determine input resistance.
    
    # Divide trace at injection time:
    
    
    # Get RMP at t < injection time:
    
    # Get RMP at t > injection time:
    
# Rheobase - find the first trace with an action potential
# Assumes traces are sorted in order from smallest amplitude upwards
    """ To Do - check that simulations is a bunch of simulations, not just an array """
def rheobase(simulations,amps):
    # Check amps is sorted
    for i in range(len(amps)-1):
        assert amps[i+1] > amps[i], 'Amps in rheobase biomarker not increasing monotonically!'
        
    for simulation,amp in zip(simulations,amps):
        # Search for first trace which produces an AP        
        result = SplitTraceIntoAPs(simulation['t'],simulation['v'])
        if result['numAPs'] > 0:
            return {'rheobase':amp, 'trace':simulation}       
    # If no APs found
    return {'rheobase':np.nan, 'trace':[]}
Rheobase = rheobase # Alias
    
def ap_peak(v):

    peak = max(v)
    location = np.argmax(v)
    return [peak,location]
APPeak = ap_peak # Alias
    
def threshold(t, v, dvdt_threshold=5.):
    # Calculation of threshold voltage as described in Davidson et al., 2014 PAIN using gradient
    # Threshold is in V/s - default of 5 is what was used by Davidson et al.
    dvdt = np.gradient(v,t)    
    thresholds = []
    for i, gradient in enumerate(dvdt[0:-1]):
        if (gradient < dvdt_threshold) & (dvdt[i+1] > dvdt_threshold): # Look for crossing of threshold
            thresholds.append(v[i])
    if thresholds:
        return thresholds[0] # Only use first threshold of crossing
    else:
        return np.nan
        
def ap_rise_time(t,v,threshold=5):
    """ 
    Threshold here is a dVdt threshold in mV/ms`
    Default threshold is taken from Davidson et al. 2014, PAIN
    """
    assert threshold > 0, 'Rise time threshold is a gradient threshold, should be > 0!'
    dVdt = np.gradient(v,t)
    peak = ap_peak(v)
    peak_idx = peak[1]
    peak_time = t[peak_idx]
    
    # If dVdt is a tuple, second part is gradient
    found_thresholds = []
    for i,gradient in enumerate(dVdt[0:-1]): # Is dVdt a time vector as well?
        if gradient < threshold:
            if dVdt[i+1] > threshold:
                found_thresholds.append(i)
    
    num_threshold = len(found_thresholds)
    if num_threshold == 1:
        threshold_time = t[found_thresholds[0]]
        rise_time = peak_time - threshold_time
        if rise_time < 0:
            #rise_time = 'Rise time < 0: %.3f' % rise_time
            rise_time = np.nan
#        assert rise_time >=0, 'Rise time < 0!'
    elif num_threshold == 0:
        rise_time = np.nan
    elif num_threshold > 1:
#        assert False, 'More than 1 threshold for rise time - APs may not be clearly separated.'
        # Take the first one - later ones are probably rapid spikes e.g. on the shoulder
        threshold_time = t[found_thresholds[0]]
        rise_time = peak_time - threshold_time
    return rise_time
APRiseTime = ap_rise_time # Alias
        
def ap_slope_min_max(t,v):
    dVdt = np.gradient(v,t)
    slope_min = min(dVdt)
    slope_max = max(dVdt)
    ### Need mins and maxes
    return [slope_min,slope_max]
APSlopeMinMax = ap_slope_min_max # Alias    

def ap_width(t, v, alpha, _threshold=5., threshold_type='gradient'):
    """
    Generic ap width calculating function. Alpha determines the fraction of the voltage
    gap between threshold and ap peak that is used to set the voltage threshold.

    Specifically, if Th is the calculate threshold voltage and P is the peak voltage,
    the voltage threshold used depends on alpha so that width threshold WTh = alpha*P + (1-alpha)*Th

    So for full width alpha = 0, as we just use the bare threshold Th, for half width alpha = 0.5,
    and alpha = 1 should give a width of 0 as it goes all the way to the peak.
    
    Defaults are consistent with Davidson et al. 2014 (5 mV/ms gradient to find threshold voltage)

    _threshold named to avoid overlapping with threshold function
    """ 
    
    # Calculate AP threshold and AP peak voltages
    if threshold_type == 'gradient':
        v_threshold = threshold(t, v, _threshold)
    elif threshold_type == 'voltage':
        v_threshold = _threshold
    else:
        raise ValueError("threshold type: {} not recognised".format(threshold_type))
    if np.isnan(v_threshold):
        return np.nan
    
    v_peak = ap_peak(v)[0]
    width_v_threshold = alpha * v_peak + (1.0 - alpha) * v_threshold

    # Find crossing points
    ups, downs = find_threshold_crossings(v, width_v_threshold)

    # Check we have crossings
    if ups and downs:
        last_down = downs[-1]
        first_up = ups[0]
        width = t[last_down] - t[first_up]
        return width
    else:
        return np.nan

def ap_full_width(t,v ,_threshold=5., threshold_type='gradient'):
    """
    Calculate full width of AP by one of two methods, a voltage threshold
    or a voltage/time gradient threshold
    Defaults are consistent with Davidson et al. 2014 (5 mV/ms gradient to find threshold voltage)

    _threshold named to avoid overlapping with threshold function
    """
    
    if threshold_type == 'voltage':
        assert not np.isnan(_threshold), "Threshold {} is nan".format(_threshold)
        ups, downs = find_threshold_crossings(v, _threshold)
    elif threshold_type == 'gradient':
        # Find voltage at which we cross the dvdt threshold
        dvdt = np.gradient(v,t)
        gradient_threshold = None
        for i, _ in enumerate(dvdt[:-1]):
            if (dvdt[i] < _threshold) and (dvdt[i+1] >= _threshold):
                gradient_threshold = v[i]
                break
        # Return if we don't cross the threshold
        if gradient_threshold:
            ups, downs = find_threshold_crossings(v, gradient_threshold)
        else:
            return np.nan
    else:
        raise ValueError("threshold type: {} not recognised".format(threshold_type))

    #print(arr)
    #print(ups,downs)
    num_ups = len(ups)
    num_downs = len(downs)

    if (num_ups < 1) | (num_downs < 1):
        # Not enough crossings
        full_width = np.nan
    elif (num_ups == 1) & (num_downs == 1): 
        # One crossing of threshold each way
        full_width = t[downs[0]] - t[ups[0]]
    elif (num_ups > 1) | (num_downs > 1):
        # Too many crossings
        # Find earliest crossing from below and latest crossing from above
        # to calculate full width
        first_up = ups[0]
        last_down = downs[-1]
        full_width = t[last_down] - t[first_up]

    return full_width

APFullWidth = ap_full_width # Alias


def ap_half_width(t,v, dvdt_threshold=5.):
    """
    Definition from neuroelectro.org:
    AP duration at membrane voltage halfway between AP threshold and AP peak.
    Currently only uses gradient method for finding threshold for simplicity.
    """
    # Calculate AP threshold and AP peak voltages
    v_threshold = threshold(t,v, dvdt_threshold=dvdt_threshold,)
    v_peak = ap_peak(v)[0]
    half_width_v_threshold = (v_threshold + v_peak)/2.

    # Find crossing points
    ups, downs = find_threshold_crossings(v,half_width_v_threshold)

    # Check we have crossings
    if ups and downs:
        last_down = downs[-1]
        first_up = ups[0]
        half_width = t[last_down] - t[first_up]
        return half_width
    else:
        return np.nan

def fit_afterhyperpolarization(traces, dvdt_threshold, max_time_from_peak=50., ahp_model = 'single_exp', full_output=False):
    """ 
    Gather afterhyperpolarisation regions from a set of traces and fit them to a model 
    of a single exponential (other models can be added as needed)
    
    Outputs:
    Returns either amp,tau or if full_output is selected returns five outputs.
    """

    # Fit to single exponential model
    if ahp_model == 'single_exp':
        def model(x, a, b, c):
            return a - b * np.exp(-x/c)
    else:
        raise ValueError('Model \"{}\" not valid'.format(ahp_model))

    # Return function if we have a result we can't fit a hyperpolarisation to
    def hyperpolarisation_fit_failure(full_output):
        if full_output:
            return np.nan, np.nan, np.nan, np.nan,np.nan, np.nan
        else:
            return np.nan, np.nan, np.nan
    
    # Arrange data to contain each interval between peaks (num APs > 1) or peak to end of trace (n=1) 
    """
    This caused a bug for thr last AP since when the stimulus turns off there is a bigger AHP than the AHP itself.
    Instead, for the last trace in the sequence only, check only for 50 ms after the peak. Could also check for some
    multiple of AP full width, but 50 ms should be sufficient for all but very long abnormal APs. 
    """
    num_APs = traces['numAPs']
    if num_APs < 1:
        return hyperpolarisation_fit_failure(full_output) 
    elif num_APs == 1:
        _t = traces['t'][0]
        _v = traces['v'][0]
        max_idx = np.argmax(_v)
        # Check that the peak is not right at the end of the trace
        if max_idx == len(_v)-1:
            return hyperpolarisation_fit_failure(full_output)
        # Get the period between t(peak) and t(peak) + max_time_from_peak
        t_peak = _t[max_idx]
        t_end = t_peak + max_time_from_peak
        end_idx = np.argmin(abs(_t - t_end))

        ts = [_t[max_idx+1:end_idx]] # Single element lists from just after peak to end of max_time_from_peak period
        vs = [_v[max_idx+1:end_idx]] 
    elif num_APs > 1:
        ts = []
        vs = []
        # Traces 1 to N-1
        for i in range(num_APs-1):
            _ts = [traces['t'][idx] for idx in [i, i+1]]
            _vs = [traces['v'][idx] for idx in [i, i+1]]
            max_idxs = [np.argmax(_v) for _v in _vs]
            # Concatenate the two parts of the interval from each trace
            _t_start = _ts[0][max_idxs[0]:]
            _t_end =  _ts[1][:max_idxs[1]-1]
            _v_start =_vs[0][max_idxs[0]:]
            _v_end =  _vs[1][:max_idxs[1]-1] 
            
            _t = np.concatenate([_t_start, _t_end], axis=0)
            _v = np.concatenate([_v_start, _v_end], axis=0)
            ts.append(_t)
            vs.append(_v)

        # Trace N - final trace - use same process as for when there is only 1 AP except change the index
        _t = traces['t'][num_APs-1]
        _v = traces['v'][num_APs-1]
        max_idx = np.argmax(_v)
        # Check that the peak is not right at the end of the trace
        if max_idx == len(_v)-1:
            return hyperpolarisation_fit_failure(full_output)
        # Get the period between t(peak) and t(peak) + max_time_from_peak
        t_peak = _t[max_idx]
        t_end = t_peak + max_time_from_peak
        end_idx = np.argmin(abs(_t - t_end))

        ts.append(_t[max_idx+1:end_idx])
        vs.append(_v[max_idx+1:end_idx])
        
    # For each interval attempt to fit an AHP
    amps = []
    taus = []
    troughs = []
    if full_output: 
        output_ts = []
        output_vs = []
        popts = []
    
    for i, (t, v) in enumerate(zip(ts, vs)):
        # Start from the minimum, until dvdt exceeds the threshold given as input
        min_idx  = np.argmin(v)
        dvdt = np.gradient(v,t)
        threshold_exceeded = dvdt > dvdt_threshold
        if any(threshold_exceeded):
            cutoff_idx = np.where(threshold_exceeded)[0][0] - 1
        else: # Use the whole trace
            cutoff_idx = len(t)-1
        
        t = t[min_idx:cutoff_idx]
        v = v[min_idx:cutoff_idx]
        
        # Check v and t are all there
        if any(np.isnan(v)) | any(np.isnan(t)):
            return hyperpolarisation_fit_failure(full_output)
        
        # If the membrane potential slowly monotonically decreases after a spike, then the min_idx will be the
        # last element of the trace, and so t and v will be empty.

        # Also, if the AHP part of the trace has a very small number of elements, then the calculation might not work.
        # So we will impose a minimum length threshold for the ahp. 
        # With dt = 0.05 ms after downsampling, and a real AHP tau of order 10 ms, any trace of
        # less than length 100 (5 ms) is probably no good, # and any trace less than length 10 (0.5 ms) 
        # is almost certainly not going to contain a viable AHP. 
        # This does assume the default dt though so we should check it is not larger than about 0.1 ms,
        # which would equate to a 1 ms minimum AHP duration.
        
        length_threshold = 10
        if (len(t) <= length_threshold) | (len(v) <= length_threshold):
            return hyperpolarisation_fit_failure(full_output)    
        dt = np.gradient(t)
        assert np.mean(dt) < 0.1, "dt is large, check length_threshold"

        # We can check for this and return failure if it is the case. To do: think about whether this is the best
        # way to handle the lack of an ahp. There could also be cases where we have an AP with no AHP, 
        # and these might give odd tau and amp readings that should not be averaged. For calibration this is fine, 
        # but for mechanistic investigation it might not be so good. 
        
        # Use scipy.optimise.curvefit to fit curve - another option would be to use the more complex 
        # LMFIT library, but I have no experience with it's advantages over the basic scipy lsq fit function
        if ahp_model == 'single_exp':
            t = t - t[0] # Zero out t as model assumes this
            popt, pcov = optimize.curve_fit(model, t, v) # Do the fitting
            # following neuroelectro: https://neuroelectro.org/ephys_prop/index/
            # AHP amplitude is from threshold voltage to trough
            trough = min(v)
            # Need original AP to calculate 
            thresh = threshold(traces['t'][i], traces['v'][i], dvdt_threshold) 
            ahp_amp = trough - thresh # will be -ve
            ahp_tau = popt[2]
            ahp_trough = trough
            " Alternate approaches to calculation"
            """
            calculation_method = False
            if calculation_method == "Simple":
                ahp_amp = trough
            elif calculation_method == "ToRecovery":
                ahp_amp = None
            """
                


        else:
            raise ValueError('Model \"{}\" not valid'.format(ahp_model))

        amps.append(ahp_amp)
        taus.append(ahp_tau)
        troughs.append(ahp_trough)
        if full_output:
            output_ts.append(t)
            output_vs.append(v)
            popts.append(popt)
    
    # Return non averaged output and cutoff times and voltages if full output requested "
    if full_output == True:
        return amps, taus, troughs, output_ts, output_vs, popts
        # Otherwise just return mean amplitude and time constant of the AHP "
    else:
        amp = np.mean(amps)
        tau = np.mean(taus)
        trough = np.mean(troughs)
        return amp, tau, trough

"""
def expFunc(t, amp, slope, start):
    return amp*(1 - np.exp(-slope*t)+start)    
    
maxIdx = []
maxIdx.append(np.argmax(v)) # Get idx of max(v)
maxIdx.append(np.argmax(v2))# Get idx of max(v2)

workingTime = np.concatenate((t[maxIdx[0]:],t2[:maxIdx[1]+1]),0) ### join t[maxIdx1:] up to t2[1:maxIdx[1]]
workingVoltage = np.concatenate((v[maxIdx[0]:],v2[:maxIdx[1]+1]),0) ### join
   
# AHP amplitude
amp = min(workingVoltage)
ampIdx = np.argmin(workingVoltage)

#    dvdt = VoltageGradient(workingTime[ampIdx:], workingVoltage[ampIdx:])
#    plt.plot(workingTime[ampIdx:-1],dvdt)
#    temp = np.argwhere(dvdt > dvdt_threshold) # Temp because we only need the first element
#    takeoffIdx = temp[0][0] # TODO This will break if there's no points above dvdt_threshold
#    plt.plot(workingTime[ampIdx:ampIdx+takeoffIdx],workingVoltage[ampIdx:ampIdx+takeoffIdx])
#    plt.plot(workingTime,workingVoltage)
# AHP time constant
# TO DO!!
# Look up curve fitting    

tau = 'Time constant not implemented'
return amp, tau
"""


def inter_spike_interval(traces):
    """ Calculate average interspike interval from a divided set of traces
        Total interspike interval is the time difference between the first and last peak of a trace,
        divided by the number of intervals (number of APs - 1)
    """
    numAPs = traces['numAPs']
    if numAPs < 2:
        #print('ISI cannot be calculated with < 2 APs')
        return np.nan
    else:
        # Find the peak of the first and last trace
        voltages = traces['v']
        first_spike = np.argmax(voltages[0])
        last_spike = np.argmax(voltages[-1])

        # Get the time difference
        times = traces['t']        
        time_diff = times[-1][last_spike] - times[0][first_spike]
        assert time_diff > 0, 'time_diff for ISI < 0: {}'.format(time_diff)
        # Divide by number of intervals (numAPs - 1) to get mean ISI
        inter_spike_interval = time_diff/(numAPs-1)
        return inter_spike_interval

def absmax(i):
    """
    Returns the largest absolute value present in an array in its raw form
    (e.g. in [-2, 0, 1] it returns -2, in [-2,0,3] it returns 3.)
    """
    # Use the absolute largest value in its raw form
    if max(i) > abs(min(i)):
        return max(i) 
    elif abs(min(i)) >= max(i):
        return min(i)
    else:
        raise ValueError()
        
# ---- Calculating biomarkers over multiple traces ----

def calculate_rmp(traces):
    RMPVals = []
    for i,v in enumerate(traces['v']):
        RMPValue, RMPIdx = RMP(v)
        RMPVals.append(RMPValue)
    return RMPVals
CalculateRMP = calculate_rmp # Alias
    
def calculate_input_res():
    input_res_vals = []
    for i,v in enumerate(traces['v']):
        input_res_vals.append(input_res(v))
    return input_res_vals
CalculateInputRes = calculate_input_res # Alias
    
def calculate_ramp_ap():
    """ 
    Can't remember what this biomarker was supposed to do? 
    We just run ramp simulations and calculate biomarkers on those now.
    """
    # TODO
    return 0
CalculateRampAP = calculate_ramp_ap # Alias
    

def calculate_rheobase(cell_model, amp_step=0.1, amp_max=5., make_plot=False, sim_kwargs=None, search='simple'):
    " Run a series of simulations to calculate rheobase"
    " Rheobase is defined as the threshold current for an infinite duration pulse "
    " We'll try 2 seconds "
     
    # Fill out sim_kwargs with defaults if needed
    if sim_kwargs is None:
        sim_kwargs = {} 
    default_kwargs = {'dur':500., 'delay':1000., 'interval':0., 'num_stims':1, 't_stop':1500.,
            'mechanisms':None, 'make_plot':False, 'plot_type':'default', 'model':cell_model}
    for kwarg in default_kwargs.keys():
        if kwarg in sim_kwargs.keys():
            pass
        else:
            sim_kwargs[kwarg] = default_kwargs[kwarg]

    def rheobase_simulation(amp):
        # Returns simulation amplitude if an AP is found, otherwise returns RHEO_FAIL if no APs found
        sim_kwargs['amp'] = amp
        output = sh.simulation(**sim_kwargs)
        t = output['t']; v = output['v'];
        # Look for an AP, after throwing away the delay period, leave a 1 ms run up to catch the start
        run_up  = 1.
        delay = sim_kwargs['delay']
        stim_period_indices = (t >= (delay-run_up)) # TODO - why did I put tuple here earlier?
        t = t[stim_period_indices]
        v = v[stim_period_indices]        
        traces = split_trace_into_aps(t,v,threshold=0.,time_threshold=5.)
        if traces['numAPs'] > 0: # rheobase found
            if make_plot:
                plot_traces(traces)
            rheobase = amp
            return rheobase
        else:
            return RHEO_FAIL

    amp_min = 0.
    amps = np.arange(amp_min, amp_max, amp_step) # (nA)
    
    # Two search modes 
    # 1. simple starts from amp_min and works up until it finds an AP
    # 2. divide starts from the middle and does a binary search
    # simple should be quicker when rheobase is usually very low and very few models have no rheobase
    # divide should be quicker if rheobase is distributed any other way

    if search == 'simple':
        for amp in amps:
            rheobase = rheobase_simulation(amp)
            if rheobase is not RHEO_FAIL: # Is not is used because np.nan == np.nan reutrns False
                return rheobase
        return RHEO_FAIL

    elif search == 'divide':
        # Divide and conquer algorithm using a binary search
        idx0 = 0
        idxn = len(amps) - 1
        rheobases = np.empty(len(amps))
        rheobases[:] = None

        while idx0 <= idxn:
            midval = (idx0 + idxn)// 2
            rheobase = rheobase_simulations(amps[midval])
            rheobases[midval] = rheobase
            if rheobase is not RHEO_FAIL: # Is not is used because np.nan == np.nan reutrns False
                if midval == 0:
                    # Rheobase is minimum
                    return amps[0]
                elif rheobases[midval-1] is not RHEO_FAIL: # Is not is used because np.nan == np.nan reutrns False
                    # Found minimal amp for an AP - return rheobase
                    return amps[midval]
                else:
                    # AP found but not definitely lowest amp so lower idxn
                    idxn = midval - 1
            elif rheobase is not RHEO_FAIL: # Is not is used because np.nan == np.nan reutrns False
                if midval == (len(amps) - 1):
                    # No rheobase for highest amp
                    return RHEO_FAIL
                elif isinstance(rheobases[midval+1], float):
                    # We've found highest amp with no AP, so one up is rheobase
                    return amps[midval+1]
                else:
                    # No AP found but not definitely highest amp so raise idx0
                    idx0 = midval + 1
            else:
                raise Exception('Rheobase not accepted value.' )
        raise Exception('No rheobase found')
    elif search == 'smart':
        # Simple search but after first two searches upwards we check the max value to check for
        # no rheobase. If the first 5? searches fail we switch to binary.
        # TODO
        pass
CalculateRheobase = calculate_rheobase # Alias for compatibility

def calculate_threshold(traces, dvdt_threshold=5.):
    thresholds = []
    for t,v in zip(traces['t'], traces['v']):
        thresholds.append(threshold(t, v, dvdt_threshold=dvdt_threshold,))
    return thresholds

def calculate_ap_peak(traces):
    ap_peak_vals = []
    for _,v in zip(range(len(traces['t'])),traces['v']):
        ap_peak_vals.append(ap_peak(v)[0])
    return ap_peak_vals
CalculateAPPeak = calculate_ap_peak # Alias
    
def calculate_ap_rise_time(traces,dvdtthreshold=5.):
    ap_rise_time_vals = []
    for t,v in zip(traces['t'],traces['v']):
        ap_rise_time_vals.append(ap_rise_time(t,v,dvdtthreshold))
    return ap_rise_time_vals
CalculateAPRiseTime = calculate_ap_rise_time # Alias
    
def calculate_ap_slope_min_max(traces):
    ap_slope_min_vals = []
    ap_slope_max_vals = [] 
    for t,v in zip(traces['t'],traces['v']):
        dvdt = np.gradient(v,t)
        ap_slope_min_vals.append(min(dvdt))
        ap_slope_max_vals.append(max(dvdt))
    return ap_slope_min_vals, ap_slope_max_vals
CalculateAPSlopeMinMax = calculate_ap_slope_min_max # Alias

def calculate_ap_width(traces, alpha, threshold=0, method='voltage'):
    ap_width_vals = []
    for t,v in zip(traces['t'],traces['v']):
        ap_width_vals.append(ap_width(t,v,alpha,threshold,method))
    return ap_width_vals

def calculate_ap_half_width(traces, threshold=0, method='voltage'):
    alpha = 0.5
    ap_half_width_vals = calculate_ap_width(traces,alpha,threshold,method)
    return ap_half_width_vals

def calculate_ap_full_width(traces,threshold=0, method='voltage'):
    alpha = 0.0 # Calculate at the threshold so set alpha = 0
    ap_full_width_vals = calculate_ap_width(traces,alpha,threshold,method)
    return ap_full_width_vals
CalculateAPFullWidth = calculate_ap_full_width # Alias

def calculate_ahp_amp(traces,dvdt_threshold=5):
    ahp_amp_vals = []
    if traces['numAPs'] > 1:
        for i in range(traces['numAPs']-1):
            t = traces['t'][i]
            v = traces['v'][i]
            t2 = traces['t'][i+1]
            v2 = traces['v'][i+1]
            amp, tau, trough = fit_afterhyperpolarization(t,v,t2,v2,dvdt_threshold)
            AHPAmpVals.append(amp)
    elif traces['numAPs'] == 1:
        v = traces['v'][0]
        max_idx = np.argmax(v)
        working_voltage = v[max_idx:]### join
        amp = min(working_voltage)
        ahp_amp_vals.append(amp)
    return ahp_amp_vals
CalculateAHPAmp = calculate_ahp_amp # Alias
    
def calculate_ahp_tau():
    # TODO
    return 0
CalculateAHPTau = calculate_ahp_tau # Alias

# -- Firing Patterns --
# See Balachandar and Prescott 2018 for algorithms
# TODO: Find algorithms for phasic and burst patterns

def determine_firing_pattern(traces, stim_start, stim_end):
    """
    Define firing pattern of traces as one or more of n types:
    1. Reluctant
    2. Single
    3. Tonic
    4. Delayed
    5. Gap
    6. Phasic - multi-AP firing that ends before end of stimulus
    7. Burst firing
    8. Wide
    9. Repolarisation failure
    """

    def first_spike_delay(traces, stim_start):
        # Find delay between stim start and first spike
        first_spike_v = traces['v'][0]
        first_spike_t = traces['t'][0]
        single_spike_index = ap_peak(first_spike_v)[1]
        single_spike_time = first_spike_t[single_spike_index]
        delay = single_spike_time - stim_start
        #print("delay = {}".format(delay))
        return delay

    def first_two_spikes_isi(traces):
        # Find delay between first and second spikes
        spike_times = []
        for i in [0,1]:
            spike_idx = ap_peak(traces['v'][i])[1]
            spike_times.append(traces['t'][i][spike_idx])
        
        delay = spike_times[1] - spike_times[0]
        return delay

    def second_third_spikes_isi(traces):
        # Find delay between second and third spikes
        spike_times = []
        for i in [1,2]:
            spike_idx = ap_peak(traces['v'][i])[1]
            spike_times.append(traces['t'][i][spike_idx])
        
        delay = spike_times[1] - spike_times[0]
        return delay

    def check_delayed(traces):
        # Check if firing pattern is delayed
        delayed = False
        num_aps = traces['numAPs']
        # Delayed firing pattern criterion for 1 spike:
        # Delay from stim start to first spike is > 100 ms
        if num_aps == 1:
            if first_spike_delay(traces, stim_start) > 100.0:
                delayed = True
        # Delayed firing pattern criterion for  > 1 spike:
        # Delay between stimulus start and firing first spike is > 1.5 
        # times the ISI between spikes 1 and 2.
        elif num_aps > 1:
            if first_spike_delay(traces, stim_start) > 1.5*first_two_spikes_isi(traces):
                delayed  = True
        return delayed

    def check_gap(traces):
        gap = False
        num_aps = traces['numAPs']
        # Gap firing criteria:
        # Number of spikes > 2
        # ISI between spikes 1 and 2 > 1.5 times ISI between spikes 2 and 3
        gap = False
        if num_aps > 2:
            if first_two_spikes_isi(traces) > 1.5*second_third_spikes_isi(traces):
                gap = True
        return gap    


    def check_phasic(traces, stim_end, ratio_threshold=0.25):
        """
        Phasic - firing of multiple APs followed by a period of quiescence. 
        Cases
        1. Idea is use ratio of - Time from last spike to stimulus end:time from first to last spike
        If the ratio is above some threshold.
        2. Simply time from last peak to end of stimulus compared to a threshold.
        """
        
        phasic = False
        # Characterisation cases
        case1 = True
        case2 = False
        
        # First, check we have multiple APs
        # We will class single spikes as single spikes, not phasic.
        num_aps = traces['numAPs']
        if num_aps < 2:
            return False
        
        spike_times = []
        for i in range(num_aps):
            spike_idx = ap_peak(traces['v'][i])[1]
            spike_times.append(traces['t'][i][spike_idx])
        
        # Case 1
        if case1:
            last_spike_to_stim_end = stim_end - spike_times[-1]
            # check stimulus ended before last spike, if not can't be phasic
            if last_spike_to_stim_end > 0:
                first_to_last_spike = spike_times[-1] - spike_times[0]
                assert first_to_last_spike > 0

                ratio = last_spike_to_stim_end/first_to_last_spike
                #print("Ratio = {}".format(ratio))
                if ratio > ratio_threshold:
                    phasic = True
        
        # Case 2
        if case2:
            raw_time_threshold = 50.0
            if last_spike_to_stimulus_end > raw_time_threshold:
                phasic = True

        return phasic

    def check_bursting(traces, stim_start):
        """
        Bursting - bursts of APs separated by rest periods
        Not sure how to characterize currently.
        1. Find all AP peaks.
        2. Divide trace up into quiet periods and firing periods
        Quiet period is region where distance between two APs or last AP and 
        stimulus end is greater than some multiple of the average ISI (median?).
        """
        bursting = False
        return bursting

    def check_wide(traces, mean_width_threshold=10.0):
        """
        Abnormally wide APs - feature seen when inserting hNav 1.8 into mice (Han et al. 2015).
        """
        wide = False

        # Get width of each AP using AP half width biomarker
        # Use half width as we can compare against data in Han et al. 2015
        # Figure 8D shows half-width distributions, which motivated the choice to 
        # set default width threshold for 'wide' designation to 10 ms.
        half_widths = []
        for t, v in zip(traces['t'], traces['v']):
            half_widths.append(ap_half_width(t ,v, dvdt_threshold=5.))

        if half_widths: # Check we have widths
            if np.mean(half_widths) > mean_width_threshold:
                wide = True
        return wide

    def check_rep_fail(traces, rep_fail_threshold=0.0):
        """
        Repolarisation failure - trace does not recover to a reasonably depolarised voltage
        This can be set by user but we'll start with using 0 mV as default threshold and 
        can tune as needed.
        """
        rep_fail = False

        last_trace = traces['v'][-1]
        # Check last element of last trace against threshold
        if last_trace[-1] > rep_fail_threshold:
            rep_fail = True
        return rep_fail


    firing_pattern = []
    num_aps = traces['numAPs']
    if num_aps == 0:
        firing_pattern.append('reluctant')
    elif num_aps == 1:
        firing_pattern.append('single')
        if check_delayed(traces):
            firing_pattern.append('delayed')
        if check_wide(traces):
            firing_pattern.append('wide')
        if check_rep_fail(traces):
            firing_pattern.append('rep_fail')
    elif num_aps > 1:
        firing_pattern.append('multi')
        # Determine if tonic spiking - can't be delayed, gap, phasic or repolarisation failure
        phasic = check_phasic(traces, stim_end, ratio_threshold=0.25)
        delayed = check_delayed(traces)
        gap = check_gap(traces)
        rep_fail = check_rep_fail(traces)
        if (not delayed) and (not gap) and (not phasic) and (not rep_fail):
            firing_pattern.append('tonic')
        if phasic:
            firing_pattern.append('phasic')
        if delayed:
            firing_pattern.append('delayed')
        if gap:
            firing_pattern.append('gap')
        if rep_fail:
            firing_pattern.append('rep_fail')
        # Check wide
        if check_wide(traces, mean_width_threshold=10.0):
            firing_pattern.append('wide')


    #print(" TODO:Bursting")
    return firing_pattern


    
# ---- Plotting ----

def plot_traces(traces):
        for t,v in zip(traces['t'], traces['v']):
            plt.plot(t,v)

# ---- I/O ----

def write_header(biomarker_file):
    string = 'Index'
    for biomarker in db.biomarkerNames:      
        string += (';' + biomarker)
    string += ';' + 'stimAmp'
    string += '\n'
    biomarker_file.write(string)
    return
WriteHeader = write_header # Alias
        
def write_biomarkers(biomarkers,biomarker_file):
    # Write the values of each biomarker in csv format
    string = str(biomarkers['Index'])    
    for biomarker in db.biomarkerNames:        
        string += (';' + str(biomarkers[biomarker]))
    string += (';' + str(biomarkers['stimAmp']))
    string += '\n'
    biomarker_file.write(string)
    return
WriteBiomarkers = write_biomarkers # Alias

# ---- Frequency intensity curve biomarkers ----

class FICurves(object):
    """
    Class to hold FI curve data
    Frequencies
    Amplitudes
    Results
    Which simulations go together
    
    And allow you to extract FI curves, plot them and obtain summary statistics
    """
    
    
    def __init__(self, results, simulations):
        self.results = results.copy()
        self.simulations = simulations.copy()
        self.groups = [] # Data storage for each group
        """ 
        Get all the simulations that have a constant stimulus amplitude and aren't run to rheobase 
        and group by block parameters and stimulus type TODO: And also check all other features (e.g. 
        """
        self.group_simulations()
        self.get_FI_curves()
        
        # Process results to remove parameters
        if 'Parameters' in self.results.columns.levels[0]:
            self.results = self.results.drop('Parameters',axis=1,level=0)
            self.results.columns = self.results.columns.drop('Parameters',level=0)
        
        
        """
        If we need multiple stim unit definitions:
        1. Turn above line into for loop
        2 Create a function called get_stim_amp_designations that gives all the things like nA or pA to search for
        3. Search for any of them in the simulation name and accept those that hit one and only one of them
        """

                
    def group_simulations(self):
        """
        Group simulations together that are the same except for their stimulus amplitudes
        """
        self.groups = []
        for name, sim in self.simulations.items():
            amp, shared_params = self.get_simulation_parameters(sim.protocols) 
            # Check for a fixed amplitude (not a simulation to find rheobase)
            if amp:
                # Check whether there is an existing group with matching shared parameters (scaling factors, stim function)
                group = self.check_for_existing_group(shared_params)
                if group is not None:
                    # Add sim name and amplitude to group as key, val
                    self.groups[group]['simulations'][name] = amp
                    #print("Appending: {}".format(self.groups[group]))
                else:
                    new_group = {'simulations': {name:amp}, 'shared_params':shared_params}
                    #print("Making: {}".format(new_group))
                    self.groups.append(new_group)  
                    
                    
    def get_simulation_parameters(self, sim_protocols):
        """
        Extract needed simulation parameters
        Currently: amplitude, stimulus type and scaling factors
        """          
        amp = sim_protocols['amp']
        shared_params = {}
        shared_params['stim_type'] = sim_protocols['stim_func']
        if 'parameter_scaling' in sim_protocols:
            shared_params['parameter_scaling'] = sim_protocols['parameter_scaling']
        return amp, shared_params
    
    
    def check_for_existing_group(self, shared_params):
        """
        Check groups for a group that mathches other_params.
        Returns:
        Group index if group exists
        None if group does not exist
        """
        # Check if other_params matches all other_params in any other group
        group_idx = None
        for i, group in enumerate(self.groups):
            if shared_params == group['shared_params']:
                assert group_idx is None, "group_idx should equal None, instead: {}".format(group_idx)
                group_idx = i

        return group_idx

                                
    def get_FI_curves(self):
        """
        Use the ISIs to compute firing curves for each group, for models that have non-nan ISIs
        Firing curves are calculated for each simulation group
        """
        num_groups = len(self.groups)
        assert num_groups > 0, "Num groups: {} is not > 0".format(num_groups)
        
        # Iterate through each group
        for group in self.groups:
            """
            Get frequencies for each simulation in the group and build into a dataframe
            Dataframe format:
            rows = model indices
            columns = simulation amplitudes
            """
            idx = self.results.index
            amps = [amp for amp in group['simulations'].values()]
            fi_data = pd.DataFrame(index=idx, columns=amps)
            # Populate this group's fi curve df 
            for sim_name, amp in group['simulations'].items():
                ISIs = self.results.loc[:,(sim_name,'ISI')]
                frequencies = self.calculate_frequencies(ISIs)
                fi_data.loc[:,amp] = frequencies
            
            # Save this group's data
            group['FI'] = fi_data
            
    def plot_FI_curves(self):
        """
        Plot FI curves and maybe compute some summary statistics
        Do scatter and line plots so that if we have a single datapoint for a model it still gets plotted
        """
        num_groups = len(self.groups)
        subplot_dim = int(np.ceil(np.sqrt(num_groups))) # Number of subplots (square)
        plt.figure(figsize=(10,10))
        
        for i, group in enumerate(self.groups):
            plt.subplot(subplot_dim,subplot_dim, i+1)
            fi_data = group['FI'].copy()
            
            for idx in fi_data.index:
                data = fi_data.loc[idx,:]
                plt.plot(data.index, data)
                plt.scatter(data.index, data)
                
            plt.xlim(min(fi_data.columns)-0.1, max(fi_data.columns)+0.1)
            plt.ylim(0, None)
            
            # Make a rough title
            separator = '_'
            for i, sim_name in enumerate(group['simulations']):
                if i == 0:
                    temp_title = sim_name
                else:
                    s1 = temp_title.split(separator)
                    s2 = sim_name.split(separator)
                    title_parts = [part for part in s1 if part in s2]
                    title = separator.join(title_parts)
            plt.title(title)
            plt.xlabel('I (nA)')
            plt.ylabel('f (Hz)')

    def calculate_frequencies(self, ISIs):
        return 1000.0/ISIs # Converts ms to Hz: ISI of 500 ms = 2 Hz, ISI of 2000 ms = 0.5 Hz

        
# ---- Util ----

def get_biomarker_names(biomarker_set='all'):
    ''' 
    Biomarkers TODO from neuroelectro:
    * Input resistance
    * AP Half width
    * Membrane time constant
    * Cell capacitance (fixed by simulation)
    * Maximum firing rate
    * Sag ratio
    * Adaptation ratio
    * First spike latency
    * FI slope
    * Spike rise time
    * Spontaneous firing rate
    * There are others but think they are other names for the same concepts
    '''
    if biomarker_set == 'all':
        biomarker_names = ['Threshold', 'APFullWidth', 'APPeak', 'APRiseTime', 'APSlopeMin', 'APSlopeMax', 'AHPAmp', 'AHPTau', 'AHPTrough', 'ISI', 'RMP', 'Rheobase']
    else:
        raise ValueError('biomarker_set {} not found'.format(biomarker_set))
    return biomarker_names

def find_threshold_crossings(arr, _threshold):
    """
    Find all indices at which a threshold is crossed from above and from below
    in an array. Used for finding indices to compute ap widths and half widths.
    """
    #print("threshold = {}".format(_threshold))
    ups = []
    downs = []
    for i, _ in enumerate(arr[:-1]): # Don't iterate on last element
        # Get crossings of threshold from below
        if arr[i] < _threshold:
            if arr[i+1] >= _threshold:
                ups.append(i)
        # Get crossings of threshold from above
        if arr[i] > _threshold:
            if arr[i+1] <= _threshold:
                downs.append(i)
    return ups, downs


def add_total_width_biomarker(pop, width_biomarker='APHalfWidth', filter_width=False, verbose=False):
    """
    Add a total width biomarker to each simulation in a population's results 
    """
    
    def compute_total_width(df, width_biomarker, filter_width=False):
        """
        Compute the total width of a simulation from its results dataframe
        with optional filtering out of AP Width outliers
        """
        freq = 1000.0/df['ISI']
        numAPs = df['numAPs']
        freq[numAPs == 1] = 1 # Approximation 

        width = df[width_biomarker]
        if filter_width:
            outlier_definition = an.get_outlier_definition(width_biomarker)
            width = width[width < outlier_definition]         
        total_width = width * freq
        total_width = total_width.fillna(0)
        return total_width
    
    simulations = [col for col in pop.results.columns.levels[0] if col not in ['Parameters']]
    for col in simulations:
        if verbose:
            print(col)
        total_width = compute_total_width(df=pop.results[col], 
                                          width_biomarker=width_biomarker,
                                          filter_width=filter_width)
        pop.results.loc[:, (col, 'APTotalWidth')] = total_width

