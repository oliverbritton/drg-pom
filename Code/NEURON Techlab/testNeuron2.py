# -*- coding: utf-8 -*-
"""
TestNeuron2

Created on Wed Oct 07 09:57:01 2015

@author: Oliver Britton
"""

import neuron
from neuron import h
from matplotlib import pyplot


soma = h.Section(name='soma') # don't need h 

# Insert active Hodgkin-Huxley current in the soma
soma.insert('hh')
soma.gnabar_hh = 0.12  # Sodium conductance in S/cm2
soma.gkbar_hh = 0.036  # Potassium conductance in S/cm2
soma.gl_hh = 0.0003    # Leak conductance in S/cm2
soma.el_hh = -54.3     # Reversal potential in mV

# Now  insert a full channel model


stim = h.IClamp(soma(0.5))

stim.delay = 5
stim.dur = 500
stim.amp = 10

v_vec = h.Vector()
t_vec = h.Vector()

v_vec.record(soma(0.5)._ref_v, sec=soma)
t_vec.record(h._ref_t)
h.finitialize(-65) # Vital! And has to go after record

tstop = 1000.0


neuron.run(tstop)


pyplot.figure(figsize=(8,4)) # Default figsize is (8,6)
pyplot.plot(t_vec, v_vec)
pyplot.xlim((0,h.tstop)) 
pyplot.ylim((-70,50))
pyplot.xlabel('time (ms)')
pyplot.ylabel('mV')
pyplot.show()

