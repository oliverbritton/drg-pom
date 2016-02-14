# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 10:56:08 2016

@author: Oliver Britton
"""

import numpy as np

a = [1,2,3,4,5,6,7]

b = np.array([1,2,3,4,5,6,7],type=float)

# Slice backwards from idx 3 
a[3::-1]

# Slice forwards in increments of 2:
a[::2]