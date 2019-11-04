# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 09:58:46 2016

@author: Oliver Britton
"""
import os
import sys
import matplotlib.pyplot as plt
sys.path.append('E:\CLPC48\Neuron Project\Code\Methods')
import PopulationOfModels as pom
import numpy as np
import time

def FigSetup():
    plt.figure(figsize=(12, 9))   
    ax = plt.subplot(111) 
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)  
    return ax

def ParseSzegedTraceFile(filename):
    trace = pom.ReadTextFile(filename)
    t = []
    v = []
    for i,line in enumerate(trace):
        if i > 4:
          t.append(float(line.split()[0]))
          v.append(float(line.split()[1]))
    return{'t':t, 'v':v}
    
def ParseSimulationFile(filename):
    trace = pom.ReadTextFile(filename)
    t = []
    v = []
    for i,line in enumerate(trace):
        if i > 2:
          t.append(float(line.split()[0]))
          v.append(float(line.split()[1]))
    return{'t':t, 'v':v}


def Tableau20():
    # These are the "Tableau 20" colors as RGB.    
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    for i in range(len(tableau20)):    
        r, g, b = tableau20[i]    
        tableau20[i] = (r / 255., g / 255., b / 255.)    
    return tableau20
    
def PlotDRGSimulations():
    plt.figure(figsize=(16, 9))          
    # Remove the plot frame lines. They are unnecessary chartjunk.    
    for i in range(100):

        ax = plt.subplot(10,10,i+1)   
        ax.spines["top"].set_visible(False)    
        ax.spines["bottom"].set_visible(True)    
        ax.spines["right"].set_visible(False)    
        ax.spines["left"].set_visible(True)   
        plt.tick_params(
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off',
            left='off',
            right='off',
            labelleft='off')
        if (i > 89):
            plt.tick_params(labelbottom='on')
            plt.xticks([0,400,800])
        plt.ylim([-80,80])
        if (i % 10) == 0:
            plt.tick_params(labelleft='on')
            plt.yticks([-60,0,60])

        filename = DRGFilename + str(i) + '.dat'
        trace = pom.ReadTraceFile(filename)
        plt.plot(trace['t'],trace['v'],linewidth=1,color=tableau20[0])
        
def PlotFourExperiments():
#    fig = plt.figure(frameon=False)
    for i,tracefile in enumerate(filenamesData):
        filename = os.path.join(directory,tracefile)
        trace = ParseSzegedTraceFile(filename)
        ax= plt.subplot(2,2,i+1)
        plt.plot(trace['t'],trace['v'],linewidth=2.5,color=tableau20[6])
        plt.ylim([-90,30])
        plt.xlim([-5,500])
#        if i == 0:
#            plt.tick_params(
#                axis='both',          # changes apply to the x-axis
#                which='both',      # both major and minor ticks are affected
#                bottom='off',      # ticks along the bottom edge are off
#                top='off',         # ticks along the top edge are off
#                labelbottom='off',
#                left='off',
#                right='off',
#                labelleft='on') # labels along the bottom edge are off
#        if i == 1:
#            plt.tick_params(
#                axis='both',          # changes apply to the x-axis
#                which='both',      # both major and minor ticks are affected
#                bottom='off',      # ticks along the bottom edge are off
#                top='off',         # ticks along the top edge are off
#                labelbottom='off',
#                left='off',
#                right='off',
#                labelleft='off') # labels along the bottom edge are off                
#        if i == 2:
        plt.tick_params(
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off',
            left='off',
            right='off',
            labelleft='off') # labels along the bottom edge are off
#        if i == 3:
#            plt.tick_params(
#                axis='both',          # changes apply to the x-axis
#                which='both',      # both major and minor ticks are affected
#                bottom='off',      # ticks along the bottom edge are off
#                top='off',         # ticks along the top edge are off
#                labelbottom='on',
#                left='off',
#                right='off',
#                labelleft='off') # labels along the bottom edge are off
                
        ax.spines["top"].set_visible(False)    
        ax.spines["bottom"].set_visible(False)       
        ax.spines["right"].set_visible(False)    
        ax.spines["left"].set_visible(False)  
#        ax.get_xaxis().tick_bottom()    
#        ax.get_yaxis().tick_left() 
    #    plt.xticks(fontsize=14) 
        if (i == 2) | (i == 3):
            plt.xlabel('Time (ms)')
        if (i == 0) | (i == 2):   
            plt.ylabel('Voltage (mV)')

def PlotExperimentsAndORd():
    plt.figure(figsize=(16, 9))          
    # Remove the plot frame lines. They are unnecessary chartjunk.    
    ax = plt.subplot(111)   
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)   

    n = [10000,568] 
    n[0] = 2000 # Make it less busy!
# All pop!
#    for i in range(n[0]):
#        filename = populationFilename + str(i+1) + '.dat'
#        trace = ParseSimulationFile(filename)
#        plt.plot(trace['t'],trace['v'],linewidth=0.5,color=tableau20[1])
        
    print "Hard bit done"
                #%%
#    for i in range(n[1]):
#        filename = calibratedFilename + str(i+1) + '.dat'
#        trace = ParseSimulationFile(filename)
#        plt.plot(trace['t'],trace['v'],linewidth=1.0,color=tableau20[0])    
#        plt.plot(trace['t'],trace['v'],linewidth=0.5,color=tableau20[1])   # For plotting with just the pale color

    filenames = pom.ReadTextFile(experimentFilenames)   

    for i,tracefile in enumerate(filenames):
        if i != 39:
            filename = os.path.join(directory,tracefile)
            trace = ParseSzegedTraceFile(filename)
            plt.plot(trace['t'],trace['v'],linewidth=1.0,color=tableau20[6]) # 6 for solid red # 7 for soft
            

        
#%% Plot individual RED traces       
    for i,tracefile in enumerate(filenames):

        a = [20]         
        if any(i == j for j in a):
            filename = os.path.join(directory,tracefile)
            trace = ParseSzegedTraceFile(filename)
            plt.plot(trace['t'],trace['v'],linewidth=1.5,color=tableau20[6]) # 6 for solid red # 7 for soft

#%% Plot ORd 

        
    trace = ParseSzegedTraceFile(filenameORdBaseline)    
    plt.plot(trace['t'],trace['v'],linewidth=1.5,color='black')
    
# Plot settings
    plt.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off',
        left='off',
        right='off',
        labelleft='off') # labels along the bottom edge are off    
    plt.ylim([-90,50])
    plt.xlim([-5,500])
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (mV)')

# Experimental traces we want:
# Min APD90: 6 : E:\Test\All\090326AK_R012_T10.dat
# Long APD90: 21 : E:\Test\All\050526AA_R013_T10.dat
# Median APD90:58 : E:\Test\All\060223AA_R006_T10.dat
# Mean APD90: 48 : E:\Test\All\100424AZ_R152_T10.dat


directory = 'E:/Test/All/'
filenamesData = ('090326AK_R012_T10.dat','060223AA_R006_T10.dat','100424AZ_R152_T10.dat','050526AA_R013_T10.dat')
filenameORdBaseline = 'E:/CLPC48/Tech Lab/Output/ORdBaseline/ORdBaselineBiphasic_1.dat'
experimentFilenames = 'E:/CLPC48/Human Project/Data/ExpFiles.txt'
populationFilename = 'E:/CLPC48\Tech Lab/Output/171013GCaL/171013GCaLControl/171013GCaLControl_'
calibratedFilename = 'E:/CLPC48/Tech Lab/Output/Mk4/Control/Mk3_State_Control_'
DRGFilename = 'E:/CLPC48/Neuron Project/Simulations/Techlab/SecondTest/SecondTest_'
tableau20 = Tableau20()

""" Main Plots"""
#PlotFourExperiments()
start = time.time()
PlotExperimentsAndORd()
#PlotDRGSimulations()
end = time.time()
print "The time is %f\n." % (end-start)  

#plt.figure(figsize=(12, 9))   
#ax = plt.subplot(111) 
#ax.spines["top"].set_visible(False)    
#ax.spines["bottom"].set_visible(False)    
#ax.spines["right"].set_visible(False)    
#ax.spines["left"].set_visible(False)   
#
#start = time.time()
#for i in range(10000):
#    filename = populationFilename + str(i+1) + '.dat'
#    trace = ParseSimulationFile(filename)
#    plt.plot(trace['t'],trace['v'],linewidth=0.5,color=tableau20[1])
#    
#end = time.time()
#print "The time is %f\n." % (end-start)  


#%%
#FigSetup()


##%%
#fig, ax = plt.subplots(2,2)
#
#count = 0
#for i in range(2):
#    for j in range(2):
#        filename = os.path.join(directory,filenamesData[count])
#        trace = ParseSzegedTraceFile(filename)
#        ax[i][j].plot(trace['t'],trace['v'])
#        ax[i][j].axis('off')
#
#        plt.show()
#        count += 1
#        
##%%
#        ax = plt.subplot(221)   

#%%


