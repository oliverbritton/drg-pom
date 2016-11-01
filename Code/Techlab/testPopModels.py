# -*- coding: utf-8 -*-
"""
testNeuronPopModels
Created on Thu Dec 10 15:21:57 2015

@author: Oliver Britton
"""

#Put everything together and do a simple population of models 
import os  # File system
import sys # To append path
from neuron import h
from matplotlib import pyplot
import time

workingDir = "E:\CLPC48\Neuron Project\Simulations\Techlab"
os.chdir(workingDir)

start = time.time()

# Set up neuron
h('load_file("nrngui.hoc")')

# Load parameters from file

paramFilename = "parameters.dat"

f = open(paramFilename,'r')

parameters = []
with open(paramFilename) as f:
    for line in f:
        #inner_list = [elt.strip() for elt in line.split(',')]
        # in alternative, if you need to use the file content as numbers
        inner_list = [float(elt.strip()) for elt in line.split(' ')]
        parameters.append(inner_list)
              
f.close()

# Set up stim protocol information
''' Do stim protocol in config file eventually '''

stimDuration = 800 #ms
stepSize = 0.05 #nA (50 pA steps)
stimAmps = [x * stepSize for x in range(101)]

# For each row in parameter file:
for i, parameterSet in enumerate(parameters):
    
    # Set up model
    cell = h.Section()
    
    cell.insert('nav17cw')
    cell.insert('nav18hw')
    cell.insert('kdrcw')
    cell.insert('kacw')
    
    cell.gnabar_nav17cw = cell.gnabar_nav17cw*parameterSet[0]
    cell.gnabar_nav18cw = cell.gnabar_nav18cw*parameterSet[1]
    cell.gkbar_kdrcw = cell.gkbar_kdrcw*parameterSet[2]
    cell.gkbar_kacw = cell.gkbar_kacw*parameterSet[3]
    
    # Set up cell properties in line with Choi Waxman
    # Membrane properties from Choi Waxman
    cell.L = 30 # microns
    cell.diam = 23*2 #microns
    
    #From Leipold et al. 2015: "Current and voltage recordings were 
    #obtained in the whole-cell configuration of the patch-clamp method 
    #from isolated small DRG neurons with an electrical capacitance of 
    #less than 20 pF, to restrict the recordings mostly to C fiber neurons
    cell.cm = (20.2)/1000  # uF/cm2
    cell.Ra = 123 # ohm-cm
    
    # Ionic concs
    cell.ko = 5.0
    cell.ki = 140.0
    
    cell.nao = 140.0
    cell.nai = 10.0

    # Open output
    f = open("Model_%i.txt" %i, "w")
        
    # Run the step protocol
    for j, amp in enumerate(stimAmps):
        
        stim = h.IClamp(cell(0.5)) 
        stim.delay = 5
        stim.dur = stimDuration
        stim.amp = amp 
        
        v_vec = h.Vector()
        t_vec = h.Vector()
            
        
        v_vec.record(cell(0.5)._ref_v, sec=cell)
        t_vec.record(h._ref_t)
        
        h.finitialize(-65) # Vital! And has to go after record   
        tstop = 1000.0
        neuron.run(tstop)
        
        # Write AP to end of file, two columns, end with line break
        assert(len(v_vec) == len(t_vec))    
        for k in range(len(t_vec)):
            f.write(str(t_vec[k]) + " " + str(v_vec[k]) + "\n")
    
        #f.write("\n")
    
    f.close()

# All done! Time it.
end = time.time()
print end-start #'s to run on 1 core.

# End'