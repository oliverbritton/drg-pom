# -*- coding: utf-8 -*-
"""
Tests to check NeuronBiomarkers works and I can calculate neuronal biomarkers

Created on Thu Mar 10 14:57:27 2016


@author: Oliver Britton
"""
import os
import NeuronProjectStart as nrnst
import Methods.pom_test as pom
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
 'Ipom_test,
 'RMP',
 'Rheobase',
 'SplitTraceIntoAPs',
 'VoltageGradient',
 """

""" Tests """
def test_split_trace(trace):
    t,v = trace['t'], trace['v']
    traces = nb.split_trace_into_aps(t,v)
    for i in range(traces['numAPs']):
        plt.plot(traces['t'][i],traces['v'][i])
        
    plt.show()
    """ TO DO Test threshold and time """
TestSplitTrace = test_split_trace
    
def TestAPFullWidth(traces, threshold, threshold_type='voltage'):    
    """ Test ap full width biomarker """
    count = 0
    for idx,(t,v) in enumerate(zip(traces['t'],traces['v']),1):
        ap_full_width = nb.ap_full_width(t,v,threshold,threshold_type)
        t_peak = t[np.argmax(v)]        

        cl = sns.color_palette()[count]
        plt.plot(t,v,color=cl)
        if threshold_type == 'voltage':
            plt.plot([t_peak,t_peak+ap_full_width],[threshold]*2,color=cl)
        elif threshold_type == 'gradient':
            # Find voltage at which we cross the dvdt threshold
            dvdt = np.gradient(v)/np.gradient(t)
            gradient_threshold = None
            for i, _ in enumerate(dvdt[:-1]):
                if (dvdt[i] < threshold) and (dvdt[i+1] >= threshold):
                    gradient_threshold = v[i]
                    break
            if gradient_threshold:
                plt.plot([t_peak,t_peak+ap_full_width],[gradient_threshold]*2,color=cl)
        print("Trace {}: Width = {} ms".format(idx, ap_full_width))    
        count +=1 
    
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

def load_test_trace():
    projectDir = nrnst.GetProjectDir()
    testTraceDir = os.path.join(projectDir,'Simulations','Test','NeuronBiomarker') 
    prefix  = 'SecondTest_'
    val = 30
    trace = pom.read_trace(os.path.join(testTraceDir,prefix+str(val)+'.dat'))
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
trace = load_test_trace()
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
