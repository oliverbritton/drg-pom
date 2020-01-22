# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 09:44:17 2016

@author: Oliver Britton
"""

""" Automatic CFG file builder """

import os
import sys

import NeuronStartProject

""" SETUP """

simulationType = 'Exploratory' # Or Results, or Techlab
modelName = 'DeterminedOcelot'
simulationName = 'FirstPop'
prefix = 'FirstPop_'
protocol = 'default'

""" """


projectPath = NeuronStartProject.GetProjectDir()
simPath = os.path.join(projectPath, 'Simulations')

cfgPath = os.path.join(simPath,'Input')
paramPath = os.path.join(cfgPath,'param')

