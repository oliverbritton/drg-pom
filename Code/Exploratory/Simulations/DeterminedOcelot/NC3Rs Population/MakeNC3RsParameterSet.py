# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:41:32 2016

@author: Oliver Britton
"""

import ParameterSets as pom



outputFilename = 'E:\\CLPC48\\Neuron Project\\Simulations\\Input\\param\\DO_NC3Rs_10000_0_10.param'

numModels = 10000
numParameters = 6
minimum = 0
maximum = 10

pom.MakeLHSParameterFile(outputFilename,numModels,numParameters,minimum,maximum)
    
