# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 11:12:00 2016

@author: Oliver Britton
"""

import math
import numpy as np

a = np.array([1,2,3])

weights = np.array([89.4,89.1,88.4,88.9,88.6,88.3,88.5,88.3,88.6,87.8,87.4,87.6])
days = 23-10

easyDiffs = np.diff(weights)

diffs =[j-i for i, j in zip(weights[:-1], weights[1:])]  # or use itertools.izip in py2k

mean = np.mean(diffs)

print (87.6-82.5)/(-mean)

daysToWedding = 8+29+25

alexDiff = 103.4-102.3

alexMean = alexDiff/days

alexW = 102.3 - (alexMean*daysToWedding)