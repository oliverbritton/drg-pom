# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 12:56:55 2016

@author: Oliver Britton
"""

""" Calibration ranges from Davidson et al. 2014, Pain """
import pandas as pd

def CalibrationData(num_stds=1.):
#    calibrationRanges = {}
    lowerRanges = {}
    upperRanges = {}
    completeDict = {}
    for key in DavMeans:
        lowerRanges = DavMeans[key] - num_stds*DavStds[key]
        upperRanges = DavMeans[key] + num_stds*DavStds[key]
        completeDict[key] = [DavMeans[key], DavStds[key], lowerRanges, upperRanges]       
    
    davData = pd.DataFrame(completeDict,columns)
    davData = davData.unstack().unstack()
    return davData
    
def GetCalibrationRanges():
    data = CalibrationData()
    calibrationRanges = data[['Min','Max']]
    return calibrationRanges

# Discussion: how would we structure this data?
# Need to play in Pandas

biomarkerNames = ['RMP','InputRes','RampAP','StepRheobase','Threshold','APPeak','APRise','APSlopeMin',
    'APSlopeMax','APFullWidth','AHPAmp', 'AHPTrough','AHPTau','numAPs'] # TODO automate this from the list of keys in DavMeans

columns = ['Mean','Std','Min','Max']

# Make a data frame and read in and then we can store simulation results in the same format
# If we need to add more biomarkers we just extend columns on the data frame

DavMeans = {}

DavMeans['RMP'] = -62.36
DavMeans['InputRes'] = 97.51
DavMeans['RampAP'] = 2.45
DavMeans['Rheobase'] = 1.43
DavMeans['Threshold'] = -15.73
DavMeans['APPeak'] = 64.64
DavMeans['APRiseTime'] = 528.4/1000
DavMeans['APSlopeMax'] = 326.9
DavMeans['APSlopeMin'] = -100.2
DavMeans['APFullWidth'] = 4.92
DavMeans['AHPAmp'] = -52.66
DavMeans['AHPTrough'] = DavMeans['AHPAmp']
DavMeans['AHPTau'] = 26.67

DavStds = {}

DavStds['RMP'] = 23.3
DavStds['InputRes'] = 111.4
DavStds['RampAP'] = 2.24
DavStds['Rheobase'] = 1.16
DavStds['Threshold'] = 10.0
DavStds['APPeak'] = 9.38
DavStds['APRiseTime'] = 359.6/1000
DavStds['APSlopeMax'] = 169.4
DavStds['APSlopeMin'] = 78.3
DavStds['APFullWidth'] = 3.73
DavStds['AHPAmp'] = 8.73
DavStds['AHPTrough'] = DavStds['AHPAmp']
DavStds['AHPTau'] = 22.1


biomarker_names = sorted(DavMeans.keys())
biomarkerNames = biomarker_names

davData = CalibrationData()

