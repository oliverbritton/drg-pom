import os
import sys
import numpy as np
import pdb

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
    if (numAPs == 0) | (numAPs == 1):
        times.append(t)
        voltages.append(v)

    # If we have multiple APs, for each AP find the minimum value
    # of Vm before the next AP
    
    """ 
    There are some commented assumptions about where traces begin and end here. The core idea is that all data points in the trace have to be assigned to 1 and only 1 AP. If areas of quiescence are a problem for particular analysis methods, they will be stripped out by other specialised functions. 
    Our goal in this function is to divide up the trace without leaving any of it out, so that we have everything for any future analysis.
    """
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
        
    assert startIdx[0] == 0, "First AP doesn't start at beginning of trace."
    assert endIdx[-1] == len(v)-1, "Last AP doesn't end at end of trace."
        
    return{'times':times, 'voltages':voltages, 'startIndices':startIdx, 'endIndices':endIdx}  
        
    
def VoltageGradient(t,v):
    ### Is there a diff function in scipy or numpy?
    dVdt = np.zeros(len(v),float)
    for i in range(len(v)-1):
        dv = v[i+1]-v[i]
        dt = t[i+1]-t[i]
        dVdt[i] = dv/dt   
    return dVdt
        
    
# --- Biomarkers ---
def RMP(t,v):
    ### Is the a min function in scipy/numpy
    return min(v)
    
def StepRheobase(tracesFromStepTest,duration):

    return
    
def APPeak(t,v):
    ### Is there a max function
    return max(v)
    
def APRiseTime(t,v,threshold):
    dVdt = VoltageGradient(t,v)
    peak = APPeak(t,v)
    peakTime = peak[1]
    
    # If dVdt is a tuple, second part is gradient
    foundThreshold = []
    for i,gradient in enumerate(dVdt[1]): # Is dVdt a time vector as well?
        if gradient < threshold:
            if dVdt[1][i+1] > threshold:
                foundThreshold.append(i)
    
    numThresholds = len(foundThresholds)
    if numThresholds == 1:
        thresholdTime = t(foundThreshold[0])
        riseTime = peakTime - thresholdTime
        assert riseTime >=0, 'Rise time < 0!'
    elif numThresholds == 0:
        riseTime = 'N/A'
    elif numThresholds > 1:
        assert False, 'More than 1 threshold for rise time - APs may not be clearly separated.'
    return riseTime
        
def APSlopeMinMax(t,v):
    dVdt = VoltageGradient(t,v)
    slopeMin = min(dVdt)
    slopeMax = max(dVdt)
    ### Need mins and maxes
    return [slopeMin,slopeMax]
    
def APFullWidth(t,v,threshold):
    ups = []
    downs = []
    for i in range(len(v)-1):
        # Find ups (cross thresh from below)
        if v[i] < threshhold:
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
        fullWidth = t[downs[0]] - trace[ups[0]]
        
    elif (numUps > 1) | (numDowns > 1):
        # Too many crossings
        # Find earliest crossing from below
        # and latest crossing from above
        # to calculate full width
        earliestUp = ups[0]
        latestDown = downs[-1]
        fullWidth = t[latestDown] - t[earliestUp]
    
    return fullWidth

def FitAfterHyperpolarisation(t,v,t2,v2):

    maxIdx = []
    maxIdx[0] = np.argmax(v) # Get idx of max(v)
    maxIdx[1]= np.argmax(v2)# Get idx of max(v2)
    
    workingTime = np.concatenate((t[maxIdx[0]:],t2[:maxIdx[1]+1]),0) ### join t[maxIdx1:] up to t2[1:maxIdx[1]]
    workingVoltage = np.concatenate((v[maxIdx[0]:],v2[:maxIdx[1]+1]),0) ### join
       
    # AHP amplitude
    amp = min(workingVoltage)
    ampIdx = np.argmin(workingVoltage)
    
    # AHP time constant
    # TO DO!!
       
    return [amp,tau]
   
 




















