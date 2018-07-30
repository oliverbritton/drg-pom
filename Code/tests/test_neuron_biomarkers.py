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

""" Current biomarkers and supporting functions tested:
'APFullWidth',
 'APPeak',
 'APRiseTime',
 'APSlopeMinMax'
 'fit_afterhyperpolarization',
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
    for idx,(t,v) in enumerate(zip(traces['t'],traces['v']),1):
        cl = sns.color_palette()[count]
        count += 1
        APFullWidth = nb.APFullWidth(t,v,threshold)
        plt.plot(t,v,color=cl)
        i = np.argmax(v)
        tPeak = t[i]        
        plt.plot([tPeak,tPeak+APFullWidth],[threshold]*2,color=cl)
        print("Trace {}: Width = {} ms".format(idx, APFullWidth))
    
    
    
def TestAPPeak(traces):
    for i,t,v in zip(range(len(traces['t'])),traces['t'],traces['v']):
        APPeak = nb.APPeak(v)
        plt.plot(t,v,color=sns.color_palette()[i])
        plt.scatter(t[APPeak[1]],APPeak[0],s=20,color=sns.xkcd_rgb["bright red"])
        print("Trace {}: APPeak = {} mV, {} ms".format(i+1,APPeak[0], t[APPeak[1]]))

def TestAPRiseTime(traces,dvdtthreshold):
    for idx,(t,v) in enumerate(zip(traces['t'],traces['v']),1):
        APRiseTime = nb.APRiseTime(t,v,dvdtthreshold)
        i = np.argmax(v)
        plt.plot(t,v)     
        plt.plot([t[i]-APRiseTime,t[i]],2*[0],color=sns.xkcd_rgb["bright red"],lw=1) 
        print("Trace {}: APRiseTime = {} ms".format(idx,APRiseTime))
        

def test_fit_afterhyperpolarization(traces,dvdtThreshold):    
        amps, taus, ts, vs, popts = nb.fit_afterhyperpolarization(traces,5,full_output=True)
        for amp,tau,t,v,popt in zip(amps,taus,ts,vs,popts):
            plt.figure()
            plt.title("Amp = {} mV, tau = {} ms".format(amp,tau))
            plt.plot(t,v)
            plt.plot(t, popt[0] - popt[1] * np.exp(-t/popt[2]))

def TestInterSpikeInterval(traces):
    for t,v in zip(traces['t'],traces['v']):
        plt.plot(t,v)
    inter_spike_interval = nb.InterSpikeInterval(traces)
    print('ISI = {} ms'.format(inter_spike_interval))
    print('Frequency = {} Hz'.format(1000./inter_spike_interval))

def TestRMP(traces):
    for t,v in zip(traces['t'],traces['v']):
        RMP, RMPIdx = nb.RMP(v)
        plt.plot(t,v)
        plt.scatter(t[RMPIdx],RMP,s=20,color=sns.xkcd_rgb["bright red"],lw=1) 
        
def TestRheobase(simulations,amps):
    
    rheobase = nb.Rheobase(simulations,amps)
    
    print(rheobase['rheobase'])
    plt.plot(rheobase['trace']['t'],rheobase['trace']['v'])
     
def TestVoltageGradient(trace):
    dVdt = nb.VoltageGradient(trace['t'],trace['v'])
    plt.plot(trace['t'][0:-1],dVdt)

def test_threshold(traces, dvdt_threshold=5.0):
    thresholds = nb.calculate_threshold(traces)
    for i, (t,v) in enumerate(zip(traces['t'],traces['v'])):
        plt.plot(t,v)
        plt.plot(t,[thresholds[i]]*len(t)) 
    
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
        
def TestAllBiomarkers(traces, model):
    " Calculate every biomarker and output to dict "
    " TODO: Use the rheobase to work out what simulation to run to calculate biomarkers "
    " off of (at rheobase) "
    biomarkers = {}
    # biomarker_names = ['APFullWidth', 'APPeak', 'APRiseTime', 'APSlopeMin', 'APSlopeMax',. 'AHPAmp', 'AHPTau', 'ISI', 'RMP', 'Rheobase']
    biomarkers['Threshold'] = np.mean(nb.calculate_threshold(traces, dvdt_threshold=5.))
    biomarkers['APFullWidth'] = np.mean(nb.CalculateAPFullWidth(traces,threshold=0))
    biomarkers['APPeak'] = np.mean(nb.CalculateAPPeak(traces))
    biomarkers['APRiseTime'] =  np.mean(nb.CalculateAPRiseTime(traces,dvdtthreshold=5))
    APSlopeMinVals, APSlopeMaxVals = nb.CalculateAPSlopeMinMax(traces)
    biomarkers['APSlopeMin'] = np.mean(APSlopeMinVals)
    biomarkers['APSlopeMax'] = np.mean(APSlopeMaxVals)
    amp, tau, trough = nb.fit_afterhyperpolarization(traces=traces,dvdt_threshold=5, ahp_model='single_exp', full_output=False)
    biomarkers['AHPAmp'] =  amp
    biomarkers['AHPTau'] =  tau
    biomarkers['AHPTrough'] = trough
    biomarkers['ISI'] = nb.inter_spike_interval(traces)
    biomarkers['RMP'] =  np.mean(nb.CalculateRMP(traces))
    # Need to do rheobase separately
    biomarkers['Rheobase'] =  nb.CalculateRheobase(model, amp_step=0.1, amp_max=5, make_plot=False,)
    
    return biomarkers
    

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
#test_fit_afterhyperpolarization(traces,dvdtThreshold)
#test_inter_spike_interval
#TestRMP(traces) # TODO make it average the voltage after the minimum """
#TestRheobase([trace,trace],[50,100])
#print("***")
