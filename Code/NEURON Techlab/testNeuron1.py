# -*- coding: utf-8 -*-
"""
Test Basic Neuron commands
Created on Tue Oct 06 12:14:03 2015

@author: Oliver Britton
"""
#%%
from neuron import h, gui
import neuron
import nrn

soma = h.Section(name='soma')

h.psection()

dir(soma)

soma.insert('pas') # insert passive leak channel

print "type(soma) = ", type(soma)
print "type(soma(0.5)) =", type(soma(0.5))

mech = soma(0.5).pas
print dir(mech)

print mech.g
print soma(0.5).g_pas

#%% Alpha synapse

asyn = h.AlphaSynapse(soma(0.5))
dir(asyn)

print "asyn.e", asyn.e
print "asyn.gmax", asyn.gmax
print "asyn.onset", asyn.onset
print "asyn.tau", asyn.tau

#%%
asyn.onset = 20
asyn.gmax = 1

h.psection()

v_vec = h.Vector() # Vm vector
t_vec = h.Vector() # Time stamp vector
v_vec.record(soma(0.5)._ref_v)
t_vec.record(h._ref_t)

#%% Running a simulation

tstop = 40.0
neuron.run(tstop)

#%% Plot

from matplotlib import pyplot

pyplot.figure(figsize=(8,4))
pyplot.plot(t_vec,v_vec)
pyplot.xlabel('time (ms)')
pyplot.ylabel('mV')
pyplot.show()

