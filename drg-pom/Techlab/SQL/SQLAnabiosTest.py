# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 16:55:22 2016

@author: Oliver Britton
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 14:22:22 2016

@author: Oliver Britton
"""

### Test of SQLite
# Use SQLite 
# Database contains tables
# Tables contain data
# All data is in 1 file

import os
import sqlite3
import time
import datetime
import random
import numpy as np

# Define connection and cursor
os.chdir("E:\CLPC48\Anabios Cardiac\Data")
# If database doesn't exist, sql will create it
conn = sqlite3.connect('Summary.db')
conn.text_factory = str
c = conn.cursor()


def readFromDB():
    # * means everything
    c.execute("SELECT [APD90(ms)] FROM FullSummary WHERE treatment = 'DMSO' AND [APD90(ms)] != '' ")
#    c.execute("SELECT * FROM stuffToPlot WHERE value > 2 AND value < 8 AND keyword='Python'")
#    data = c.fetchall() # can also use fetchone() to get one row
#    print(data)   
    x = []
    for row in c.fetchall():
        print(row)
        x.append(float(row[0]))

    return x
    
def readFromDBTreatment(treatment):
    # * means everything
    c.execute("SELECT [APD90(ms)] FROM FullSummary WHERE treatment = '%s' AND [APD90(ms)] != '' " % treatment)
#    c.execute("SELECT * FROM stuffToPlot WHERE value > 2 AND value < 8 AND keyword='Python'")
#    data = c.fetchall() # can also use fetchone() to get one row
#    print(data)   
    x = []
    for row in c.fetchall():
        print(row)
        x.append(float(row[0]))

    return x
    
def readQueryFromDB(query):

    #    c.execute("SELECT [APD90(ms)] FROM FullSummary WHERE treatment = '%s' AND pacingfrequency GLOB '1*' AND concentration = 'High' AND [APD90(ms)] != '' " % treatment)
    c.execute(query)
    x = []
    for row in c.fetchall():
        x.append(row)
    return x
#createTable()
#dataEntry()
#    
#for i in range(10):
#    dynamicDataEntry()
#    time.sleep(1) # To make timestamp go up a second
start = time.time()
#x=readFromDB()

#drugs = ('DMSO','Dofetilide','A','B','C','D','E')
#drugs = ('DMSO')

col = "APD50(ms)"
queryM = "SELECT [%s], Donorage FROM FullSummary WHERE treatment = 'DMSO' AND pacingfrequency GLOB '1*' AND [APD90(ms)] != '' AND donorsex = 'M' " % col
queryF = "SELECT [%s], Donorage FROM FullSummary WHERE treatment = 'DMSO' AND pacingfrequency GLOB '1*' AND [APD90(ms)] != '' AND donorsex = 'F' " % col
APDMeans = []
APDStds = []
APDAges = []

#        APD.append(float(row[0]))
#    APDMeans.append(np.mean(APD))
#    APDStds.append(np.std(APD))
    
    # DO GENDER GENDE    

males = readQueryFromDB(queryM)
females = readQueryFromDB(queryF)

c.close()
conn.close()
end = time.time()

M = np.zeros([len(males),2])
F = np.zeros([len(females),2])

for i,line in enumerate(males):
    M[i,0] = float(line[0])
    M[i,1] = float(line[1])
    
for i,line in enumerate(females):
    F[i,0] = float(line[0])
    F[i,1] = float(line[1])    
print "Queries took %f seconds.\n" % (end-start)

plt.subplot(1,2,1)
plt.scatter(F[:,1],F[:,0])
plt.xlabel('Age')
plt.ylabel(col)
plt.xlim([10,70])
plt.ylim([0,500])
plt.title('Women')

plt.subplot(1,2,2)
plt.scatter(M[:,1],M[:,0])
plt.xlabel('Age')
plt.ylabel(col)
plt.xlim([10,70])
plt.ylim([0,500])
plt.title('Men')