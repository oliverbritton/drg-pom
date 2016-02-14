# -*- coding: utf-8 -*-
"""
Constructs population of models parameter files

Created on Fri Feb 05 13:25:57 2016

@author: Oliver Britton
"""

import numpy as np
import pyDOE

def MakeLHSParameterFile(outputFilename,numModels,numParameters,minimum,maximum):

    parameterSets = pyDOE.lhs(numParameters,numModels)
    
    # Transform parameter set using minimum and maximum
    
    # Scale size of range
    parameterSets *= (maximum-minimum)
    # Scale to right minimum (and right maximum if we scaled range right)
    parameterSets += minimum     
    
    
    np.savetxt(outputFilename,parameterSets,fmt='%.3f') # Write to 3 d.p precision 
    return


