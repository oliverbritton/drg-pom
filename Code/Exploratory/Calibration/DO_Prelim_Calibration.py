# -*- coding: utf-8 -*-
"""
Calibration Script for DeterminedOcelot Davidson preliminary population

Created on Sat Apr 09 17:52:09 2016

@author: Oliver Britton
"""
import os
import PopulationOfModels as pom
import Data.DavidsonBiomarkers as db

parameterPath = ''
biomarkerPath = ''
prefix = ''
calibrationRanges = db.GetCalibrationRanges()

# Use Calibration function in pom to get calibrated indices


# Use pom write function to write population parameters
# and calibration file

pom.WriteCalibration()



