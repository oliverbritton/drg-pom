# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 10:02:44 2016

TestDeterminedOcelot

Test a 3xNa 3xK+ model

@author: Oliver Britton


"""
import os
import sys

#import neuron

# If we load nrnmech.dll twice we will crash, so make sure our working directory is clear and load manually
os.chdir('E:\\CLPC48\\Neuron Project\\Code\\Models\\Currents\\Prototypes')
from neuron import h
#h.nrn_load_dll("E:\\CLPC48\\Neuron Project\\Code\\Models\\Currents\\Prototypes\\nrnmech.dll")
os.chdir('E:\\CLPC48\\Neuron Project\\Simulations\\')

from matplotlib import pyplot
import time


h('load_file("nrngui.hoc")')

start = time.time()

cell = h.Section()

# Insert channels

# Na channels - Vasylyev Nav 1.7, Han 1.8, Huang 1.9

cell.insert('nav17vw')
cell.insert('nav18hw')
cell.insert('nav19hw')

# K channels - Kdr, KA, KM
cell.insert('kdrtf')
cell.insert('katf')
cell.insert('kmtf')

cell.insert('naiTest')


# Modify conductances
cell.gnabar_nav17vw = cell.gnabar_nav17vw*1.0
cell.gnabar_nav18hw = cell.gnabar_nav18hw*0.9
cell.gnabar_nav19hw = cell.gnabar_nav19hw*1.0


cell.gkbar_kdrtf = cell.gkbar_kdrtf*1.0
cell.gkbar_katf = cell.gkbar_katf*1.0
cell.gkbar_kmtf = cell.gkbar_kmtf*1.0

# Set up cell properties in line with Choi Waxman

# Membrane properties from Choi Waxman
cell.L = 30 # microns
cell.diam = 23*2 #microns
cell.cm = (20.2)/1000  # uF/cm2
cell.Ra = 123 # ohm-cm

# Ionic concs
cell.ko = 5.0
cell.ki = 140.0

cell.nao = 140.0
cell.nai = 10.0


#%% 
stim = h.IClamp(cell(0.5))

stim.delay = 100
stim.dur = 700
stim.amp = 0.5 # nA (1 nA = 100 pA)

v_vec = h.Vector()
t_vec = h.Vector()
ina_vec = h.Vector()
icurr_vec = h.Vector()

v_vec.record(cell(0.5)._ref_v, sec=cell)
t_vec.record(h._ref_t)
ina_vec.record(cell(0.5)._ref_ina)
icurr_vec.record(cell(0.5)._ref_ina_nav18hw, sec=cell)

h.finitialize(-65) # Vital! And has to go after record
cell.nai = 10
tstop = 1000.0


neuron.run(tstop)

end = time.time()
print end - start


pyplot.figure(figsize=(8,4)) # Default figsize is (8,6)
pyplot.plot(t_vec, v_vec)
#pyplot.plot(t_vec,ina_vec)
#pyplot.plot(t_vec,icurr_vec)

pyplot.xlim((0,h.tstop)) 
#pyplot.ylim((-90,50))
pyplot.xlabel('time (ms)')
pyplot.ylabel('mV')
pyplot.show()

#%% Output results to file

#==============================================================================
# outputFilename = 'outputSpikes.txt'
# 
# a = t_vec
# b = v_vec
# 
# f = open(outputFilename, "w")
# for i in range(len(a)):
#     f.write(str(a[i]) + " " + str(b[i]) + "\n")
#     
# f.close()
#==============================================================================

def LoopyLoops(paramList,idx):
    
    # Set up the parameters for the current idx, using paramList
    
    # Param list
    
    activeCell = MakeDeterminedOcelot()
    
    SetDeterminedOcelotParameters(activeCell)
    
    cell.gnabar_nav17vw = cell.gnabar_nav17vw*paramList[0]
    cell.gnabar_nav18hw = cell.gnabar_nav18hw*paramList[1]
    cell.gnabar_nav19hw = cell.gnabar_nav19hw*paramList[2]


    cell.gkbar_kdrtf = cell.gkbar_kdrtf*paramList[3]
    cell.gkbar_katf = cell.gkbar_katf*paramList[4]
    cell.gkbar_kmtf = cell.gkbar_kmtf*paramList[5]
    
def MakeDeterminedOcelot():
    
    cell = h.Section()

    # Insert channels
    
    # Na channels - Vasylyev Nav 1.7, Han 1.8, Huang 1.9
    cell.insert('nav17vw')
    cell.insert('nav18hw')
    cell.insert('nav19hw')
    
    # K channels - Kdr, KA, KM
    cell.insert('kdrtf')
    cell.insert('katf')
    cell.insert('kmtf')
    
    # Modify conductances
    cell.gnabar_nav17vw = cell.gnabar_nav17vw*1.0
    cell.gnabar_nav18hw = cell.gnabar_nav18hw*0.9
    cell.gnabar_nav19hw = cell.gnabar_nav19hw*1.0
    
    
    cell.gkbar_kdrtf = cell.gkbar_kdrtf*1.0
    cell.gkbar_katf = cell.gkbar_katf*1.0
    cell.gkbar_kmtf = cell.gkbar_kmtf*1.0
    
    # Set up cell properties in line with Choi Waxman
    
    # Membrane properties from Choi Waxman
    cell.L = 30 # microns
    cell.diam = 23*2 #microns
    cell.cm = (20.2)/1000  # uF/cm2
    cell.Ra = 123 # ohm-cm
    
    # Ionic concs
    cell.ko = 5.0
    cell.ki = 140.0
    
    cell.nao = 140.0
    cell.nai = 10.0

    #print cell.gnabar_nav19hw
    #print "Cell constructed!"
    return cell
    
    
def TestMakeCell():
    # Test making a cell
    newCell = MakeDeterminedOcelot()
    IncreasedL = newCell.L + 1 
    newCell = MakeDeterminedOcelot()
    assert(newCell.L + 1 == IncreasedL)
    return
    
# REFACTOR refactor this function as a method in the model class, then can be DRY compatible on list of conductances
def ScaleDeterminedOcelotParameters(cell,parameters):
	
    numParameters = 6
    assert(len(parameters) == numParameters)
    	
    # Sodium conductances
    cell.gnabar_nav17vw = cell.gnabar_nav17vw*parameters[0]
    cell.gnabar_nav18hw = cell.gnabar_nav18hw*parameters[1]
    cell.gnabar_nav19hw = cell.gnabar_nav19hw*parameters[2]
    
    # Potassium conductances
    cell.gkbar_kdrtf = cell.gkbar_kdrtf*parameters[3]
    cell.gkbar_katf = cell.gkbar_katf*parameters[4]
    cell.gkbar_kmtf = cell.gkbar_kmtf*parameters[5]
    return
    
def TestScaleDeterminedOcelotParameters():

    cell = MakeDeterminedOcelot()
    parameterScaling = [1,10,100,1000,10000,1000000]
    originalParameters = []
    newParameters = []
    
    originalParameters.append(cell.gnabar_nav17vw)	
    originalParameters.append(cell.gnabar_nav18hw)
    originalParameters.append(cell.gnabar_nav19hw)
    originalParameters.append(cell.gkbar_kdrtf)
    originalParameters.append(cell.gkbar_katf)
    originalParameters.append(cell.gkbar_kmtf)
    
    ScaleDeterminedOcelotParameters(cell,parameterScaling)
    
    newParameters.append(cell.gnabar_nav17vw)	
    newParameters.append(cell.gnabar_nav18hw)
    newParameters.append(cell.gnabar_nav19hw)
    newParameters.append(cell.gkbar_kdrtf)
    newParameters.append(cell.gkbar_katf)
    newParameters.append(cell.gkbar_kmtf)
    
    # Test that we can recover original parameters
    for i in zip(originalParameters,parameterScaling,newParameters):
            # Test that the original parameter set times the scaling factor equals the new parameter
            assert(i[0]*i[1] == i[2])
    return
            
       
