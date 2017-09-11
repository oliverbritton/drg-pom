# -*- coding: utf-8 -*-
"""
Constructs population of models parameter files

Created on Fri Feb 05 13:25:57 2016

@author: Oliver Britton
"""

import numpy as np
import pyDOE

def MakeLHSParameterFile(outputFilename,numModels,numParameters,minimum,maximum, save=True):

    parameterSets = pyDOE.lhs(numParameters,numModels)
    
    # Transform parameter set using minimum and maximum
    
    # Scale size of range
    parameterSets *= (maximum-minimum)
    # Scale to right minimum (and right maximum if we scaled range right)
    parameterSets += minimum     
    
    np.savetxt(outputFilename,parameterSets,fmt='%.3f') # Write to 3 d.p precision 



def build_parameter_set(num_models, num_parameters, min, max, output_filename=None, paramater_names=None, save=True):

    parameter_sets = pyDOE.lhs(num_parameters, num_models)
    # Transform parameter set using minimum and maximum
    # Scale size of range
    parameter_sets *= (max-min)
    # Scale to right minimum (and right maximum if we scaled range right)
    parameter_sets += min    
    
    if output_filename != None:
        np.savetxt(output_filename, parameter_sets, fmt='%.3f') # Write to 3 d.p precision 
    return parameter_sets