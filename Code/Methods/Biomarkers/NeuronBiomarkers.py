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
    biomarkers = {}
    # biomarker_names = ['APFullWidth', 'APPeak', 'APRiseTime', 'APSlopeMin', 'APSlopeMax',. 'AHPAmp', 'AHPTau', 'ISI', 'RMP', 'Rheobase']

    biomarkers['APFullWidth'] = np.mean(CalculateAPFullWidth(traces,threshold=0))
    biomarkers['APPeak'] = np.mean(CalculateAPPeak(traces))
    biomarkers['APRiseTime'] =  np.mean(CalculateAPRiseTime(traces,dvdtthreshold=5))
    APSlopeMinVals, APSlopeMaxVals = CalculateAPSlopeMinMax(traces)
    biomarkers['APSlopeMin'] = np.mean(APSlopeMinVals)
    biomarkers['APSlopeMax'] = np.mean(APSlopeMaxVals)
    amp, tau = FitAfterHyperpolarisation(traces=traces,dvdt_threshold=5, ahp_model='single_exp', full_output=False)
    biomarkers['AHPAmp'] =  amp
    biomarkers['AHPTau'] =  tau
    biomarkers['ISI'] = InterSpikeInterval(traces)
    biomarkers['RMP'] =  np.mean(CalculateRMP(traces))
    # Need to do rheobase separately
    biomarkers['Rheobase'] =  CalculateRheobase(model, amp_step=0.1, amp_max=5, make_plot=False,)
    
    return biomarkers

def SplitTraceIntoAPs(t,v,threshold=20,timeThreshold=5):#
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
                voltageDuringCurrentAP = v[ firstCrossings[AP]:firstCrossings[AP+1] ]
                
                # Get min voltage index
                minVmIdx = np.argmin(voltageDuringCurrentAP)
                endIdx[AP] = firstCrossings[AP] + minVmIdx # Don't think I need to minus 1 because Python indices start at 0
                
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
        
    
def VoltageGradient(t,v):
    ### Is there a diff function in scipy or numpy?
    dVdt = np.zeros(len(v)-1,float)
    for i in range(len(v)-1):
        dv = v[i+1]-v[i]
        dt = t[i+1]-t[i]
        dVdt[i] = dv/dt   
    return dVdt
        
    
# --- Biomarkers ---
def RMP(v):
    # RMP should be calculated from a quiescent trace (no stimulus)
    # Ignore first 90% of trace to remove artifacts
    vLen = len(v)
    startIdx = 90*vLen/100
    RMP = min(v[startIdx:])
    RMPIdx = np.argmin(v[startIdx:]) 
    return RMP, RMPIdx
    
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
    return {'rheobase':'N/A', 'trace':[]}
    
def APPeak(v):

    peak = max(v)
    location = np.argmax(v)
    return [peak,location]
    
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
            riseTime = 'Rise time < 0: %.3f' % riseTime
#        assert riseTime >=0, 'Rise time < 0!'
    elif numThresholds == 0:
        riseTime = 'N/A'
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
        fullWidth = 'N/A'
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

def FitAfterHyperpolarisation(traces, dvdt_threshold, ahp_model = 'single_exp', full_output=False):
    """ 
    Gather afterhyperpolarisation regions from a set of traces and fit them to a model 
    of a single exponential (other models can be added as needed)
    """
    # Define model to fit AHP to
    if ahp_model == 'single_exp':
        def model(x, a, b, c):
            return a - b * np.exp(-x/c)
    else:
        assert False
    
    # Arrange data to contain each interval between peak
    num_APs = traces['numAPs']
    if num_APs < 1:
        return False
      
    elif num_APs == 1:
        # With only one AP the interval is from the peak to the end of the trace
        _t = traces['t'][0]
        _v = traces['v'][0]
        max_idx = np.argmax(_v)
        ts = [_t[max_idx+1:]] # Single element lists from just after peak to end of trace
        vs = [_v[max_idx+1:]]
        
    elif num_APs > 1:
        # Divide data into intervals between peaks for each interval
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
            
    # For the trace from each interval, fit an AHP and store parameters
    amps = []
    taus = []
    if full_output: # Storage if interval traces are requested
        output_ts = []
        output_vs = []
        popts = []
        
    for t,v in zip(ts,vs):
        # Start from the minimum, until dvdt exceeds the threshold given as input
        min_idx  = np.argmin(v)
        dvdt = np.gradient(v)/np.gradient(t)
        threshold_exceeded = dvdt > dvdt_threshold
        if any(threshold_exceeded):
            cutoff_idx = np.where(threshold_exceeded)[0][0] - 1
        else: # Use the whole trace
            cutoff_idx = len(t)-1
        
        t = t[min_idx:cutoff_idx]
        v = v[min_idx:cutoff_idx]
        
        # use scipy.optimise.curvefit to fit curve - another option would be to use the more complex 
        # LMFIT library, but I have no experience with it's advantages over the basic scipy lsq fit function
        if ahp_model == 'single_exp':
            t = t - t[0] # Zero out t as model assumes this
            popt, pcov = optimize.curve_fit(model, t, v) # Do the fitting
        else:
            assert False, "ahp_model not found"
        
        amps.append(popt[0] - popt[1]) # calculate amp
        taus.append(popt[2]) # get tau
        
        if full_output:
            output_ts.append(t)
            output_vs.append(v)
            popts.append(popt)
    
    # Return non averaged output and cutoff times and voltages if full output requested "
    if full_output == True:
        return amps, taus, output_ts, output_vs, popts
        # Otherwise just return mean amplitude and time constant of the AHP "
    else:
        amp = np.mean(amps)
        tau = np.mean(taus)
        return amp, tau
        
    #return amp, tau
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
#    temp = np.argwhere(dvdt > dvdtThreshold) # Temp because we only need the first element
#    takeoffIdx = temp[0][0] # TODO This will break if there's no points above dvdtThreshold
#    plt.plot(workingTime[ampIdx:ampIdx+takeoffIdx],workingVoltage[ampIdx:ampIdx+takeoffIdx])
#    plt.plot(workingTime,workingVoltage)
    # AHP time constant
    # TO DO!!
    # Look up curve fitting    
    
    tau = 'Time constant not implemented'
    return amp, tau
    """
def InterSpikeInterval(traces):
    # Calculate average interspike interval from a divided set of traces
    numAPs = traces['numAPs']
    if numAPs < 2:
        print('ISI cannot be calculated with < 2 APs')
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

# ---- Calculating biomarkers over multiple traces ----

def CalculateRMP(traces):
    RMPVals = []
    for i,v in enumerate(traces['v']):
        RMPValue, RMPIdx = RMP(v)
        RMPVals.append(RMPValue)
    return RMPVals
    
def CalculateInputRes():
    # TODO
    return 0
    
def CalculateRampAP():
    # TODO
    return 0
    
def CalculateRheobase(cell_model, amp_step=0.1, amp_max=5, make_plot=False,):
    " Run a series of simulations to calculate rheobase "
    " Rheobase is defined as the threshold current for an infinite duration pulse "
    " We'll try 2 seconds "
    
    amps = np.arange(0,amp_max,amp_step) # (nA)

    # Run simulations until we reach the limit or we find an AP and return rheobase
    # Return NaN if rheobase not found

    for amp in amps:
        t,v = sh.simulation(amp=amp, dur=2000., delay=1000., interval=0, num_stims=1, mechanisms=None, t_stop=3000., make_plot=False, plot_type='default', model=cell_model)
        
        # Look for an AP (after the delay), if one is found then return amp as rheobase amplitude
        traces = SplitTraceIntoAPs(t,v,threshold=0,timeThreshold=5)
        if traces['numAPs'] > 0: # rheobase found
            if make_plot:
                plot_traces(traces)
            rheobase = amp
            return rheobase

    # No APs found - return not a number
    return np.nan
    
def CalculateThreshold(model):

    return 0

def CalculateAPPeak(traces):
    APPeakVals = []
    for i,v in zip(range(len(traces['t'])),traces['v']):
        APPeakVals.append(APPeak(v)[0])
    return APPeakVals
    
def CalculateAPRiseTime(traces,dvdtthreshold=5):
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
            amp,tau = FitAfterHyperpolarisation(t,v,t2,v2,dvdtThreshold)
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
        
















