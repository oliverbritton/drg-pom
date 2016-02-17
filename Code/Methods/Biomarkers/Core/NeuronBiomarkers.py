import os
import sys
import numpy as np


def SplitTraceIntoAPs(t,v):
	
    # ///Magic numbers\\\
    threshold = 0 # mV (assumed)
    timeThreshold = 0 # ms (assumed)
    
    assert len(t) == len(v), "v and t length mismatch"
    crossings = []
    timeCrossings = []
    
    # Check this for loop ignores last element
    for i,voltage in enumerate(v(1:-1)): # Check this works
        if voltage < threshold:
            if v[i+1] > threshold:
                crossings.append(i)
               timeCrossings.append(t(i))
                
    # For each crossing, remove all instances within 5 ms, leaving only the first crossing of the threshold 
    groupedCrossings = # Make array of length crossings-1
    for i ... ### WE NEED TO CHECK THIS ACTUALLY WORKS ON THE NEW DATA BEFORE PORTING
    
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
        
def APSlopeMinMax(t,v)
    dVdt = VoltageGradient(t,v)
    slopeMin = min(dVdt)
    slopeMax = max(dVdt)
    ### Need mins and maxes
    return [slopeMin,slopeMax]
    
def APFullWidth(t,v,threshold)
    ups = []
    downs = []
    for i in range(len(v)-1):
        # Find ups (cross thresh from below)
        if v[i] < threshhold
            if v[i+1] >= threshold # Equals here so we trigger once if we flatten off at exactly threshold
                ups.append(i)
        #Find downs (cross thresh from above)
        if v[i] > threshold
            if v[i+1] <= threshold
                downs.append(i)
    numUps = len(ups)
    numDowns = len(downs)

    if (numUps < 1) | (numDowns < 1):
        # Not enough crossings
        fullWidth = 'N/A'
    elif (numUps == 1) & (numDowns == 1) 
        # One crossing of threshold each way
        fullWidth = t[downs[0]] - trace[ups[0]]
        
    elif (numUps > 1) | (numDowns > 1)
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
    
    workingTime = np.concatenate((t[maxIdx[0]:],t2[:maxIdx[1]+1),0) ### join t[maxIdx1:] up to t2[1:maxIdx[1]]
    workingVoltage = np.concatenate((v[maxIdx[0]:],v2[:maxIdx[1]+1),0) ### join
   
   # AHP amplitude
   [amp ampIdx] = min(workingVoltage)
   
   # AHP time constant
   # TO DO!!
   
   return [amp,tau]
   
 




















