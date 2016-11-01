# -*- coding: utf-8 -*-
"""
Create a model that can change its intracellular sodium concentration
na_accu loads in sodium and produces an equal but opposite non-specific charge
so there is no effect on voltage

Created on Fri Oct 28 10:42:03 2016

@author: Oliver Britton
"""


from neuron import h
import neuron
from matplotlib import pyplot as plt
import time

import Methods.PopulationOfModels as pom

import Methods.Simulations.loadneuron as loadneuron

if 'mechanisms_loaded' not in locals():
    loadneuron.load_neuron_mechanisms()
    mechanisms_loaded = True    


start = time.time()

cell = h.Section()

cell.insert('na_accu')
cell.gna_na_accu = 100.0
cell.insert('nadifl')
print "Yay"
#cell.gnabar_nav18cw = cell.gnabar_nav18cw*0.007
#cell.gnabar_nav17vw = cell.gnabar_nav17vw*0.0000000000001


# Set up cell properties in line with Choi Waxman

# Membrane properties from Choi Waxman
cell.L = 30 # microns
cell.diam = 23*2 #microns
cell.cm = (20.2)/1000  # uF/cm2
cell.Ra = 123 # ohm-cm

# Ionic concs
try:
    cell.ko = 5.0
    cell.ki = 140.0
except NameError:
    pass

cell.nao = 140.0
cell.nai = 10.0

#%% 
stim = h.IClamp(cell(0.5))

stim.delay = 100
stim.dur = 10
stim.amp = 0.0 # nA (1 nA = 100 pA)

v_vec = h.Vector()
t_vec = h.Vector()
na_vec = h.Vector()
ena_vec = h.Vector()
print cell(0.5).nai
v_vec.record(cell(0.5)._ref_v, sec=cell)
t_vec.record(h._ref_t)
na_vec.record(cell(0.5)._ref_nai)
ena_vec.record(cell(0.5)._ref_ena)

h.finitialize(-65) # Vital! And has to go after record
cell.nai = 10.0
tstop = 1000.0

neuron.run(tstop)

end = time.time()
print end - start
print "final cell nai is {}".format(cell(0.5).nai)

plt.figure(figsize=(8,4)) # Default figsize is (8,6)
plt.plot(t_vec, na_vec)
plt.xlim((0, h.tstop)) 
#plt.ylim((0,25))
plt.xlabel('time (ms)')
plt.ylabel('nai')
plt.show()
plt.figure(figsize=(8,4)) # Default figsize is (8,6)
plt.plot(t_vec, ena_vec)
plt.xlim((0,h.tstop)) 
plt.ylabel('ena')
plt.show()

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
