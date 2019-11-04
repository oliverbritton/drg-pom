# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 15:18:12 2016

@author: Oliver Britton
"""

import os
import matplotlib.pyplot as plt
import PopulationOfModels as pom
import seaborn as sns
import pandas as pd


plotDir = 'E:/CLPC48/Neuron Project/Simulations/Techlab/SecondTest'

for i in range(100):
    
    name = os.path.join(plotDir,'SecondTest_%i.dat' % i)    
    data = pom.ReadTraceFile(name)  
    t = data['t']
    v = data['v']
    
    plt.subplot(10,10,i+1)
    plt.plot(t,v)
    
# Grab biomarkers
df=pd.read_csv('E:/CLPC48/Neuron Project/Simulations/Techlab/SecondTestSecondTest__biomarkers.dat',sep=';')

# To read a column, try:
df['APWidth'].values # or equivalent