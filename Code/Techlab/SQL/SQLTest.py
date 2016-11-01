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

import sqlite3
import time
import datetime
import random

# Define connection and cursor

# If database doesn't exist, sql will create it
conn = sqlite3.connect('tutorial.db')
c = conn.cursor()

# Let's create a table
def createTable():
    # Convention only! SQL is blind to casing!
    # ALL CAPS - Pure SQL
    # Regular - non-SQL commands
    c.execute('CREATE TABLE IF NOT EXISTS stuffToPlot(unix REAL, datestamp TEXT, keyword TEXT,value REAL)')
    
def dataEntry():
    c.execute("INSERT INTO stuffToPlot VALUES(123,'01/01/2016', 'Python', 5)")
    conn.commit() # Need to do to make changes
    c.close()    
    conn.close()
    
def dynamicDataEntry():
    unix = time.time()    
    date = str(datetime.datetime.fromtimestamp(time.time()).strftime('%d-%m-%Y %H:%M:%S'))
    keyword = 'Python'
    value = random.randrange(0,10)    
    # Question marks are sqllite specific, mysql would use %s instead of ?
    # Datatypes: REAL INT BLOB TEXT and one other (NONE?) (for SQLLite, mySQL has more)
    c.execute("INSERT INTO stuffToPlot (unix, datestamp, keyword, value) VALUES (?, ?, ?, ?)", (unix, date, keyword, value))
    conn.commit()    
    
def readFromDB():
    # * means everything
    c.execute("SELECT keyword,value,datestamp FROM stuffToPlot WHERE value > 2 AND value < 8 AND keyword='Python'")
#    c.execute("SELECT * FROM stuffToPlot WHERE value > 2 AND value < 8 AND keyword='Python'")
#    data = c.fetchall() # can also use fetchone() to get one row
#    print(data)    
    for row in c.fetchall():
        print(row)
        print(row[-1])
#    c.execute('SELECT datestamp FROM stuffToPlot')
##    data = c.fetchall() # can also use fetchone() to get one row
##    print(data)    
#    for row in c.fetchall():
#        print row
        
    
#createTable()
#dataEntry()
#    
#for i in range(10):
#    dynamicDataEntry()
#    time.sleep(1) # To make timestamp go up a second
    
readFromDB()
c.close()
conn.close()

