import numpy as np
import DavidsonBiomarkers as db

# Biomarkers to manage and analyse neuronal simulation data and potentially experimental
# data too


def SplitTraceIntoAPs(t,v,threshold=0,timeThreshold=5):
	
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
    # Ignore first 1% of trace to remove artifacts
    vLen = len(v)
    startIdx = vLen/100
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

def FitAfterHyperpolarisation(t,v,t2,v2,dvdtThreshold):

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

   
# Calculate average interspike interval from a divided trace calculated
# by the SplitTraces function
def InterSpikeInterval(dividedTrace):
    # Get number of spikes
    numAPs = dividedTrace['numAPs']
    if numAPs < 2:
        return 'N/A'
    else:
        # Find the peak of the first and last trace
        startIndices = dividedTrace(['startIndices'])
        endIndices = dividedTrace(['endIndices'])
        voltages = dividedTrace['voltages']
        
        firstSpike = argmax(voltages[0])
        lastSpike = argmax(voltages[-1])
        
#        # Don't think I need these
#        firstIndex = startIndices[0] + firstSpike
#        lastIndex = startIndices[-1] + lastSpike
        
        # Get the time difference
        times = dividedTrace['times']        
        timeDiff = times[-1][lastSpike] - times[0][firstSpike]
        assert timeDiff > 0, 'timeDiff for ISI < 0!'
        # Divide by (numAPs - 1)
        interSpikeInterval = timeDiff/(numAps-1)
        return interSpikeInterval

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
    
def CalculateStepRheobase():
    # TODO
    return 0
    
def CalculateThreshold():
    #TODO
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
        
















