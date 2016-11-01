# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 11:17:27 2015

testNeuronIO
Testbed for writing Neuron output to file

@author: Oliver Britton
"""

# Load a file of parameter sets
paramFilename = "params.dat"

f = open(paramFilename,'r')


list_of_lists = []
with open(paramFilename) as f:
    for line in f:
        #inner_list = [elt.strip() for elt in line.split(',')]
        # in alternative, if you need to use the file content as numbers
        inner_list = [float(elt.strip()) for elt in line.split(' ')]
        list_of_lists.append(inner_list)
        
        
f.close()

# Set parameter values
i = 0
GNav17_temp *= list_of_lists[i][0]
GNav18_temp *= list_oflist[i][1] 

# RUN SIMULATION

# Reset model


# How to write a couple of vectors to file
#%% This way is easy

a = [1,2,3]
b = [4,5,6]

f = open("list.txt", "w")
for i in range(len(a)):
    f.write(str(a[i]) + " " + str(b[i]) + "\n")
    
f.close()
#%% This way also works but has some 

a=[1,2,3]
b=[4,5,6]

import csv
with open('text.csv', 'w') as f:
    writer = csv.writer(f, delimiter=' ')
    writer.writerows(zip(a,b))
    

