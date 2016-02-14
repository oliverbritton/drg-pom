# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:41:32 2016

@author: Oliver Britton
"""

import ParameterSets as pom


outputFilename = 'E:\\CLPC48\\Neuron Project\\Simulations\\Input\\DeterminedOcelot_100_0_2.param'

numModels = 100
numParameters = 6
minimum = 0
maximum = 2

pom.MakeLHSParameterFile(outputFilename,numModels,numParameters,minimum,maximum)