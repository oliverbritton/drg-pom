# -*- coding: utf-8 -*-
"""
Created on Thu Oct 08 09:47:18 2015

testNeuron3.py

Test making a 4-current neuronal model using the
Choi Waxman currents, and plot response to a stimulus.

@author: Oliver Britton
"""

import sys # To append path

#import neuron
from neuron import h

from matplotlib import pyplot

import time

sys.path.append('E:\CLPC48\Neuron Project\Code')
sys.path.append('E:\CLPC48\Neuron Project\Code\Models\Currents\Prototypes')

h('load_file("nrngui.hoc")')
#

start = time.time()

cell = h.Section()

cell.insert('nav17vw')
cell.insert('nav18cw')
cell.insert('kdrcw')
cell.insert('kacw')

#cell.gnabar_nav18cw = cell.gnabar_nav18cw*0.007
#cell.gnabar_nav17vw = cell.gnabar_nav17vw*0.0000000000001

cell.gkbar_kdrcw

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
stim.dur = 10
stim.amp = 0.5 # nA (1 nA = 100 pA)

v_vec = h.Vector()
t_vec = h.Vector()

v_vec.record(cell(0.5)._ref_v, sec=cell)
t_vec.record(h._ref_t)
h.finitialize(-65) # Vital! And has to go after record

tstop = 1000.0


neuron.run(tstop)

end = time.time()
print end - start


pyplot.figure(figsize=(8,4)) # Default figsize is (8,6)
pyplot.plot(t_vec, v_vec)
pyplot.xlim((0,h.tstop)) 
pyplot.ylim((-90,50))
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
