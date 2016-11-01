# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 15:19:12 2015

@author: Oliver Britton
"""

''' Test Python Loops'''

# enumerate is a good time

a = [1, 10, 100]

for i, val in enumerate(a):
    print val, i
    
# Can get index and value in a vector easily (e.g. for parameters)