# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 08:56:08 2016

@author: Oliver Britton
"""
import NeuronProjectStart
import os
import matplotlib.pyplot as plt
import Methods.PopulationOfModels as pom
import seaborn as sns
import pandas as pd
import numpy as np
import Methods.Biomarkers.DavidsonBiomarkers as dav
import ast

def FixDF(df,biomarkerList):
    
    
    for biomarker in biomarkerList:
        col = np.array(df.loc[:,biomarker])
        for i, val in enumerate(col):
            try:
                col[i] = float(val)
            except:
                col[i] = np.nan
        df[biomarker] = col
    return df
    
def CalibrateDF(df,biomarkerList,davData):
    for biomarker in biomarkerList:
        # Minimum
        minimum = davData['Min'][biomarker]
        maximum = davData['Max'][biomarker]
        df = df[df[biomarker] >= minimum]
        df = df[df[biomarker] <= maximum]
    return df    
    

plotDir = 'E:/CLPC48/Neuron Project/Simulations/Exploratory/DeterminedOcelot/NC3Rs'
prefix = 'NC3Rs'

#for i in range(100):
#    
#    name = os.path.join(plotDir,prefix + '_%i.dat' % i)    
#    data = pom.ReadTraceFile(name)  
#    t = data['t']
#    v = data['v']
#    
#    plt.subplot(10,10,i+1)
#    plt.plot(t,v)
    
# Grab biomarkers
df=pd.read_csv(os.path.join('E:/CLPC48/Neuron Project/Simulations/Exploratory/DeterminedOcelot','NC3RsNC3Rs__biomarkers_nobrackets.dat'),sep=';')

# To read a column, try:

biomarkers = ['RMP','StepRheobase','APPeak','APRise','APSlopeMin','APSlopeMax','APWidth','AHPAmp','numAPs']
data = {}

#%%

biomarkerList = ['RMP','APWidth','APSlopeMin','APSlopeMax']

d = FixDF(df,biomarkerList)
d = CalibrateDF(d,biomarkerList,dav.davData)

#%%
#df['APWidth'].values # or equivalent
def Junk():
    for biomarker in biomarkers:
        a=[]
        for i in df[biomarker].values:
            temp = i
            if (str(temp) != 'nan' and str(temp) != "['N/A']"):
                if type(temp) == list: 
                    for j in temp:
                        if j == 'N/A': 
                            continue
                    temp = temp[0]
                if type(temp) == str:
                    temp = ast.literal_eval(temp)
                    if len(temp) > 1:
                        continue # go to next
                    if type(temp) == list:
                        if type(temp[0]) == str:
                            continue
                a.append(temp)
                    
        data[biomarker] = a
        
    #b = data['APWidth']
    #c = []
    #for i in b:
    #    if type(i[0]) == float:
    #        c.append(i)        
    
    # 
    for biomarker in biomarkers:
        plt.figure()
        sns.distplot(data[biomarker])
        plt.title(biomarker)    
        
#%%
plt.figure()
sns.set(style="whitegrid", color_codes=True)

def Plot(i,d):
    count = 0
    name = os.path.join(plotDir,prefix + '_%i.dat' % idx)    
    data = pom.ReadTraceFile(name)  
    t = data['t']
    v = data['v']
    if t[np.argmax(v)] > 90:    
        plt.plot(t,v)
        count += 1
    plt.xlim([90,160])
    return count
#models = [1046,4483]
#for i in range(2):
#    plt.subplot(1,2,i+1)
#    Plot(models[i],d)
count = 0
for i in range(0,len(d),1+int(len(d)//100)):
    idx = d.iloc[i][0]
    count += Plot(idx,d)
    plt.xlabel('Time (ms)')
    plt.ylabel('Membrane potential (mV)')
