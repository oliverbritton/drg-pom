# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:22:36 2016

@author: Oliver Britton
"""
import os
import sys
import matplotlib.pyplot as plt
sys.path.append('E:\CLPC48\Neuron Project\Code\Methods')
import PopulationOfModels as pom
import numpy as np
import Biomarkers.Core.NeuronBiomarkers as nb


# Read in a data file
plotDir = 'E:/CLPC48/Neuron Project/Simulations/Techlab/SecondTest'
plt.figure(1)

numToPlot = 100
vmaxes = []
for i in range(numToPlot):
    

    name = os.path.join(plotDir,'SecondTest_%i.dat' % i)    
    data = pom.ReadTraceFile(name)  
    t = data['t']
    v = data['v']
    
    vmax = max(v)
    vmaxes.append(vmax)
    
#    plt.subplot(10,10,i+1)
#    plt.plot(t,[vmax]*len(v))
#    plt.plot(t,v)
    
parameters = pom.ReadParameterFile('E:\CLPC48\Neuron Project\Simulations\Input\param\DeterminedOcelot_100_0_2.param')    
parameters = np.array(parameters)

names = ['Na 1.7','Na 1.8','Na 1.9','Kdr','KA','KM']

for i in range(6):
    plt.subplot(3,2,i+1)
    plt.scatter(parameters[:,i],vmaxes)
    plt.xlabel(names[i])
    #%%
order = np.argsort(parameters[:,1]) # Nav 1.8

plt.figure(1)
for i,index in enumerate(order):
    plt.subplot(10,10,i+1)
    name = os.path.join(plotDir,'SecondTest_%i.dat' % index)    
    data = pom.ReadTraceFile(name)  
    t = data['t']
    v = data['v']
    plt.plot(t,v)
    plt.ylim([-70, 80])
    
# %
    
    

    