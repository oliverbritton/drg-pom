# -*- coding: utf-8 -*-
"""
Tests to check NeuronBiomarkers works and I can calculate neuronal biomarkers

Created on Thu Mar 10 14:57:27 2016


@author: Oliver Britton
"""
import os
import NeuronProjectStart as nrnst
import Methods.PopulationOfModels as pom
import Methods.Biomarkers.NeuronBiomarkers as nb
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

""" Current functions tested:
'APFullWidth',
 'APPeak',
 'APRiseTime',
 'APSlopeMinMax',
 'FitAfterHyperpolarisation',
 'InterSpikeInterval',
 'RMP',
 'Rheobase',
 'SplitTraceIntoAPs',
 'VoltageGradient',
 """

""" Tests """
def TestSplitTrace(trace):
    
    traces = nb.SplitTraceIntoAPs(trace['t'],trace['v'])
    for i in range(traces['numAPs']):
        plt.plot(traces['t'][i],traces['v'][i])
        
    plt.show()
    """ TO DO Test threshold and time """
    
def TestAPFullWidth(traces,threshold):    
    
    count = 0
    for t,v in zip(traces['t'],traces['v']):
        cl = sns.color_palette()[count]
        count += 1
        APFullWidth = nb.APFullWidth(t,v,threshold)
        plt.plot(t,v,color=cl)
        i = np.argmax(v)
        tPeak = t[i]        
        plt.plot([tPeak,tPeak+APFullWidth],[threshold]*2,color=cl)
    
    
    
def TestAPPeak(traces):
    for i,t,v in zip(range(len(traces['t'])),traces['t'],traces['v']):
        APPeak = nb.APPeak(v)
        plt.plot(t,v,color=sns.color_palette()[i])
        plt.scatter(t[APPeak[1]],APPeak[0],s=20,color=sns.xkcd_rgb["bright red"])
    

def TestAPRiseTime(traces,dvdtthreshold):
    for t,v in zip(traces['t'],traces['v']):
        APRiseTime = nb.APRiseTime(t,v,dvdtthreshold)
        print APRiseTime
        i = np.argmax(v)
        plt.plot(t,v)     
        plt.plot([t[i]-APRiseTime,t[i]],2*[0],color=sns.xkcd_rgb["bright red"],lw=1) 
        

def TestFitAfterHyperpolarisation(traces,dvdtThreshold):    
    
#    FitAfterHyperpolarisation(t,v,t2,v2):
    for i in range(traces['numAPs']-1):
        t = traces['t'][i]
        v = traces['v'][i]
        t2 = traces['t'][i+1]
        v2 = traces['v'][i+1]
        amp,tau = nb.FitAfterHyperpolarisation(t,v,t2,v2,dvdtThreshold)
#
#def TestInterSpikeInterval():
#
def TestRMP(traces):
    for t,v in zip(traces['t'],traces['v']):
        RMP, RMPIdx = nb.RMP(v)
        plt.plot(t,v)
        plt.scatter(t[RMPIdx],RMP,s=20,color=sns.xkcd_rgb["bright red"],lw=1) 
        
def TestRheobase(simulations,amps):
    
    rheobase = nb.Rheobase(simulations,amps)
    
    print rheobase['rheobase']
    plt.plot(rheobase['trace']['t'],rheobase['trace']['v'])
    

     
def TestVoltageGradient(trace):
    dVdt = nb.VoltageGradient(trace['t'],trace['v'])
    plt.plot(trace['t'][0:-1],dVdt)


    
""" Load test traces """

def LoadTestTrace():
    projectDir = nrnst.GetProjectDir()
    testTraceDir = os.path.join(projectDir,'Simulations','Test','NeuronBiomarker') 
    prefix  = 'SecondTest_'
    val = 30
    trace = pom.ReadTraceFile(os.path.join(testTraceDir,prefix+str(val)+'.dat'))
    return trace
    
def LoadTestRheobaseTraces():
    projectDir = nrnst.GetProjectDir()
    testTraceDir = os.path.join(projectDir,'Simulations','Test','NeuronBiomarker') 
    prefix  = 'RheobaseTest_'
    amps = range([0,4000+1,50])    
    # Set up a numpy array to hold all of the simulations
    simulations = np.array()
    
    # Iterate
    i = 1
    for amp in amps:
        filename = os.path.join(testTraceDir,prefix,str(amp),'_',str(i),'.dat')

""" Setup """
trace = LoadTestTrace()
traces = nb.SplitTraceIntoAPs(trace['t'],trace['v'])
threshold = 0
dvdtThreshold = 5

#TestSplitTrace(trace)
#TestVoltageGradient(trace)
#TestAPFullWidth(traces,threshold)
#TestAPPeak(traces)
#TestAPRiseTime(traces,dvdtThreshold)
#TestFitAfterHyperpolarisation(traces,dvdtThreshold)
#TestInterSpikeInterval
TestRMP(traces) # TODO make it average the voltage after the minimum """
#TestRheobase([trace,trace],[50,100])
print "***"
