# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:22:36 2016

@author: Oliver Britton
"""
import os
import sys
import matplotlib.pyplot as plt
import PopulationOfModels as pom

# Read in a data file
plotDir = 'E:/CLPC48/Neuron Project/Simulations/Techlab/SecondTest'
plt.figure(1)

for i in range(100):
    

    name = os.path.join(plotDir,'SecondTest_%i.dat' % i)    
    data = pom.ReadTraceFile(name)  
    t = data['t']
    v = data['v']
    
    plt.subplot(10,10,i+1)
    plt.plot(t,v)