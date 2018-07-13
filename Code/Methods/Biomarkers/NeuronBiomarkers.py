import sys

import numpy as np
from scipy import optimize
import Methods.Biomarkers.DavidsonBiomarkers as db
import Methods.simulation_helpers as sh
from matplotlib import pyplot as plt

# Biomarkers to manage and analyse neuronal simulation data and potentially experimental
# data too

def calculate_biomarkers(traces, model):
    " Calculate every biomarker and output to dict "
    " TODO: Use the rheobase to work out what simulation to run to calculate biomarkers "
    " off of (at rheobase) "
    # biomarker_names = ['APFullWidth', 'APPeak', 'APRiseTime', 'APSlopeMin', 'APSlopeMax',. 'AHPAmp', 'AHPTau', 'ISI', 'RMP', 'Rheobase']

    biomarkers = calculate_simple_biomarkers(traces, model)
    biomarkers['RMP'] =  np.mean(CalculateRMP(traces))
    # Need to do rheobase separately
    biomarkers['Rheobase'] =  CalculateRheobase(model, amp_step=0.1, amp_max=5, make_plot=False,)
    
    return biomarkers
    
def average_biomarker_values(biomarkers, how_to_handle_nans='return'):
    " Average biomarker values for multiple APs while handling biomarkers that "
    if how_to_handle_nans == 'return': # Advantage is we return nan if there are any nans - good for calibration and trouble shooting - shows up weird models easily. 
        pass
    elif how_to_handle_nans == 'remove': # Risky option. Advantage is we still get a number back in mixed cases of nan and non-nan biomarkers, which is potentially risky as it hides a problem in one or more APs.
        biomarkers = biomarkers[~np.isnan(biomarkers)]
    else:
        raise ValueError('Not an accepted value.')
        
    mean_result = np.mean(biomarkers)
    return mean_result
    
def calculate_simple_biomarkers(traces, model, how_to_handle_nans='return'):
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
        biomarkers['APFullWidth'] = average_biomarker_values(CalculateAPFullWidth(traces,threshold=0), how_to_handle_nans)
    except:
        error_handle('fullwidth.pickle',traces)
    biomarkers['APPeak'] = average_biomarker_values(CalculateAPPeak(traces),how_to_handle_nans)
    
    try:
        biomarkers['APRiseTime'] = average_biomarker_values(CalculateAPRiseTime(traces,dvdtthreshold=5),how_to_handle_nans)
    except:
        error_handle('risetime.pickle',traces)
    APSlopeMinVals, APSlopeMaxVals = CalculateAPSlopeMinMax(traces)
    biomarkers['APSlopeMin'] = average_biomarker_values(APSlopeMinVals, how_to_handle_nans)
    biomarkers['APSlopeMax'] = average_biomarker_values(APSlopeMaxVals, how_to_handle_nans)
    biomarkers['Threshold'] = average_biomarker_values(calculate_threshold(traces), how_to_handle_nans)
    
    amp, tau, trough = FitAfterHyperpolarisation(traces=traces,dvdt_threshold=5, ahp_model='single_exp', full_output=False)
    """
    try:
        amp, tau = FitAfterHyperpolarisation(traces=traces,dvdt_threshold=5, ahp_model='single_exp', full_output=False)
    except:
        error_handle('fitahp.pickle',traces)
        amp=0
        tau=0
    """
    biomarkers['AHPAmp'] =  amp
    biomarkers['AHPTau'] =  tau
    biomarkers['AHPTrough'] = trough
    biomarkers['ISI'] = InterSpikeInterval(traces)
        
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
        biomarkers = calculate_simple_biomarkers(traces,model,how_to_handle_nans='return')

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

def split_trace_into_aps(t,v,threshold=20,timeThreshold=5):#
    " Threshold is at +20 mV to avoid RF causing spurious AP detection "
	
    # Units for defaults
    # t, timeThreshold - ms
    # v, threshold - mV

    assert len(t) == len(v), "v and t length mismatch"
        
    

    crossings = []
    timeCrossings = np.array([])
    
    # Check this for loop ignores last element
    # Looks for crossings from below
    for i,voltage in enumerate(v[:-1]): # Check this works
        if voltage < threshold:
            if v[i+1] >= threshold:
                crossings.append(i)
                timeCrossings = np.append(timeCrossings,t[i])
                
    # For each crossing, remove all instances within 5 ms, leaving only the first crossing of the threshold 
    groupedCrossings = np.zeros(np.size(crossings),float) 
    
    for i in range(len(crossings)-1):
        if groupedCrossings[i] == 0:
            nearbyCrossings = np.array( (timeCrossings[i+1:] - timeCrossings[i]) < timeThreshold )
            # Assign 
            groupedCrossings[i+1:] += nearbyCrossings
            assert all(groupedCrossings < 2), "Grouped crossing grouped more than once"
            
    firstCrossIndices = np.where(groupedCrossings == 0)
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

    # If we have multiple APs, for each AP find the minimum value
    # of Vm before the next AP
    
    """ 
    There are some commented assumptions about where traces begin and end here. The core idea is that all data points in the trace have to be assigned to 1 and only 1 AP. If areas of quiescence are a problem for particular analysis methods, they will be stripped out by other specialised functions. 
    Our goal in this function is to divide up the trace without leaving any of it out, so that we have everything for any future analysis.
    """
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
        
SplitTraceIntoAPs = split_trace_into_aps # alias
    
def VoltageGradient(t,v, method='gradient'):
    # There is a gradient function in numpy to take central differences
    if method == 'gradient':
        dvdt = np.gradient(v)/np.gradient(t) # Central differences except at end points
    elif method == 'diff':
        dvdt = np.diff(v)/np.diff(t) # DIfference between adjacent points 
    else :
        raise ValueError("Method not found.")
    return dvdt
    
# --- Biomarkers ---
def RMP(v):
    # RMP should be calculated from a quiescent trace (no stimulus)
    # Ignore first 90% of trace to remove artifacts
    vLen = len(v)
    startIdx = 90*vLen//100
    RMP = min(v[startIdx:])
    RMPIdx = np.argmin(v[startIdx:]) 
    return RMP, RMPIdx
    
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
def Rheobase(simulations,amps):
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
    
def APPeak(v):

    peak = max(v)
    location = np.argmax(v)
    return [peak,location]
    
def threshold(t, v, dvdt_threshold=5., method='gradient'):
    # Calculation of threshold voltage as described in Davidson et al., 2014 PAIN
    # Threshold is in V/s - default of 5 is what was used by Davidson et al.
    dvdt = VoltageGradient(t,v, method=method)    
    thresholds = []
    for i,gradient in enumerate(dvdt[0:-1]):
        if (gradient < dvdt_threshold) & (dvdt[i+1] > dvdt_threshold): # Look for crossing of threshold
            thresholds.append(v[i])
    if thresholds:
        return thresholds[0] # Only use first threshold of crossing
    else:
        return np.nan
        
    
    # Threshold here is a dVdt threshold in V/s!
    # Default threshold is taken from Davidson et al. 2014, PAIN
def APRiseTime(t,v,threshold=5):
    assert threshold > 0, 'Rise time threshold is a gradient threshold, should be > 0!'
    dVdt = VoltageGradient(t,v)
    peak = APPeak(v)
    peakIdx = peak[1]
    peakTime = t[peakIdx]
    
    # If dVdt is a tuple, second part is gradient
    foundThresholds = []
    for i,gradient in enumerate(dVdt[0:-1]): # Is dVdt a time vector as well?
        if gradient < threshold:
            if dVdt[i+1] > threshold:
                foundThresholds.append(i)
    
    numThresholds = len(foundThresholds)
    if numThresholds == 1:
        thresholdTime = t[foundThresholds[0]]
        riseTime = peakTime - thresholdTime
        if riseTime < 0:
            #riseTime = 'Rise time < 0: %.3f' % riseTime
            riseTime = np.nan
#        assert riseTime >=0, 'Rise time < 0!'
    elif numThresholds == 0:
        riseTime = np.nan
    elif numThresholds > 1:
#        assert False, 'More than 1 threshold for rise time - APs may not be clearly separated.'
        # Take the first one - later ones are probably rapid spikes e.g. on the shoulder
        thresholdTime = t[foundThresholds[0]]
        riseTime = peakTime - thresholdTime
        
    return riseTime
        
def APSlopeMinMax(t,v):
    dVdt = VoltageGradient(t,v)
    slopeMin = min(dVdt)
    slopeMax = max(dVdt)
    ### Need mins and maxes
    return [slopeMin,slopeMax]
    
def APFullWidth(t,v,threshold=0):
    ups = []
    downs = []
    for i in range(len(v)-1):
        # Find ups (cross thresh from below)
        if v[i] < threshold:
            if v[i+1] >= threshold: # Equals here so we trigger once if we flatten off at exactly threshold
                ups.append(i)
        #Find downs (cross thresh from above)
        if v[i] > threshold:
            if v[i+1] <= threshold:
                downs.append(i)
    numUps = len(ups)
    numDowns = len(downs)

    if (numUps < 1) | (numDowns < 1):
        # Not enough crossings
        fullWidth = np.nan
    elif (numUps == 1) & (numDowns == 1): 
        # One crossing of threshold each way
        fullWidth = t[downs[0]] - t[ups[0]]
        
    elif (numUps > 1) | (numDowns > 1):
        # Too many crossings
        # Find earliest crossing from below
        # and latest crossing from above
        # to calculate full width
        earliestUp = ups[0]
        latestDown = downs[-1]
        fullWidth = t[latestDown] - t[earliestUp]
    
    return fullWidth

def fit_afterhyperpolarization(traces, dvdt_threshold, ahp_model = 'single_exp', full_output=False):
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
        ts = [_t[max_idx+1:]] # Single element lists from just after peak to end of trace
        vs = [_v[max_idx+1:]] 
    elif num_APs > 1:
        ts = []
        vs = []
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
        dv = np.gradient(v)
        dt = np.gradient(t)
        dvdt = dv/dt
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
        assert np.mean(dt) < 0.1, "dt is large, check length_threshold"
        if (len(t) <= length_threshold) | (len(v) <= length_threshold):
            return hyperpolarisation_fit_failure(full_output)
            
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

FitAfterHyperpolarisation = fit_afterhyperpolarization # Alias

'''
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
#    temp = np.argwhere(dvdt > dvdtThreshold) # Temp because we only need the first element
#    takeoffIdx = temp[0][0] # TODO This will break if there's no points above dvdtThreshold
#    plt.plot(workingTime[ampIdx:ampIdx+takeoffIdx],workingVoltage[ampIdx:ampIdx+takeoffIdx])
#    plt.plot(workingTime,workingVoltage)
# AHP time constant
# TO DO!!
# Look up curve fitting    

tau = 'Time constant not implemented'
return amp, tau
'''


def InterSpikeInterval(traces):
    # Calculate average interspike interval from a divided set of traces
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
        assert time_diff > 0, 'time_diff for ISI < 0'
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
    
def CalculateInputRes():
    input_res_vals = []
    for i,v in enumerate(traces['v']):
        input_res_vals.append(input_res(v))
    return input_res_vals
    
    
def CalculateRampAP():
    # TODO
    return 0
    

def calculate_rheobase(cell_model, amp_step=0.1, amp_max=5., make_plot=False, sim_kwargs=None, search='simple'):
    " Run a series of simulations to calculate rheobase"
    " Rheobase is defined as the threshold current for an infinite duration pulse "
    " We'll try 2 seconds "
     
    RHEO_FAIL = 'no_rheobase' # Failure code for simulations with no rheobase found

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
        stim_period_indices = [t >= (delay-run_up)]
        t = t[stim_period_indices]
        v = v[stim_period_indices]        
        traces = SplitTraceIntoAPs(t,v,threshold=0.,timeThreshold=5.)
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
            if rheobase != RHEO_FAIL:
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
            if rheobase != RHEO_FAIL:
                if midval == 0:
                    # Rheobase is minimum
                    return amps[0]
                elif rheobases[midval-1] == RHEO_FAIL:
                    # Found minimal amp for an AP - return rheobase
                    return amps[midval]
                else:
                    # AP found but not definitely lowest amp so lower idxn
                    idxn = midval - 1
            elif rheobase == RHEO_FAIL:
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
        thresholds.append(threshold(t, v, dvdt_threshold=dvdt_threshold, method='gradient'))
    return thresholds

def CalculateAPPeak(traces):
    APPeakVals = []
    for _,v in zip(range(len(traces['t'])),traces['v']):
        APPeakVals.append(APPeak(v)[0])
    return APPeakVals
    
def CalculateAPRiseTime(traces,dvdtthreshold=5.):
    APRiseTimeVals = []
    for t,v in zip(traces['t'],traces['v']):
        APRiseTimeVals.append(APRiseTime(t,v,dvdtthreshold))
    return APRiseTimeVals
    
def CalculateAPSlopeMinMax(traces):
    APSlopeMinVals = []
    APSlopeMaxVals = [] 
    for t,v in zip(traces['t'],traces['v']):
        dVdt = VoltageGradient(t,v)
        APSlopeMinVals.append(min(dVdt))
        APSlopeMaxVals.append(max(dVdt))
    return APSlopeMinVals, APSlopeMaxVals
    
def CalculateAPFullWidth(traces,threshold=0):
    
    APFullWidthVals = []
    for t,v in zip(traces['t'],traces['v']):
        APFullWidthVals.append(APFullWidth(t,v,threshold))
    return APFullWidthVals
    
def CalculateAHPAmp(traces,dvdtThreshold=5):
    AHPAmpVals = []
    if traces['numAPs'] > 1:
        for i in range(traces['numAPs']-1):
            t = traces['t'][i]
            v = traces['v'][i]
            t2 = traces['t'][i+1]
            v2 = traces['v'][i+1]
            amp, tau, trough = FitAfterHyperpolarisation(t,v,t2,v2,dvdtThreshold)
            AHPAmpVals.append(amp)
    elif traces['numAPs'] == 1:
        v = traces['v'][0]
        maxIdx = np.argmax(v)
        workingVoltage = v[maxIdx:]### join
        amp = min(workingVoltage)
        AHPAmpVals.append(amp)
    
    return AHPAmpVals
    
def CalculateAHPTau():
    # TODO
    return 0
    
# ---- Plotting ----

def plot_traces(traces):
        for t,v in zip(traces['t'], traces['v']):
            plt.plot(t,v)

# ---- I/O ----

def WriteHeader(biomarkerFile):
    string = 'Index'
    for biomarker in db.biomarkerNames:      
        string += (';' + biomarker)
    string += ';' + 'stimAmp'
    string += '\n'
    biomarkerFile.write(string)
    return
        
def WriteBiomarkers(biomarkers,biomarkerFile):
    # Write the values of each biomarker in csv format
    string = str(biomarkers['Index'])    
    for biomarker in db.biomarkerNames:        
        string += (';' + str(biomarkers[biomarker]))
    string += (';' + str(biomarkers['stimAmp']))
    string += '\n'
    biomarkerFile.write(string)
    return
        
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















