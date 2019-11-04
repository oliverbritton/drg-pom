# -*- coding: utf-8 -*-
"""
Created on Mon Feb 08 11:06:00 2016

@author: Oliver Britton
"""

from neuron import h
from Methods.Simulations.loadneuron import load_neuron_mechanisms


def MakeDeterminedOcelot():
    
    load_neuron_mechanisms()    
    
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
    
def MakeDeterminedOcelotIonic():
    
    cell = MakeDeterminedOcelot()
    
    # Add ionic trackers
    cell.insert('k_conc')
    cell.insert('na_conc')    
    
    
    return cell