# -*- coding: utf-8 -*-
"""
Test Basic Python
Created on Thu Oct 01 10:34:17 2015

@author: Oliver Britton
"""

import neuron
import nrn

print "Hello"

#%% Lists
my_list = [1, 3, 5, 8, 13]

print my_list

print len(my_list)

print my_list
print my_list[-1]


print my_list
print my_list[2:4]
print my_list[2:-1]
print my_list[:2]
print my_list[2:]

list_a = [1,3,5,8,13]
list_b = list(list_a)
list_b.reverse()
print "List a ", list_a
print "List b ", list_b

mixedList = ['abc', 1.0, 2, "a string"]
print mixedList
print mixedList[0]

#%% Range

print range(10) 
print range(0,-10,-2)

#%% Loops and iterators

someRange = range(10)
for i in someRange:
    print "The value:",i  

for i in xrange(100000): # xrange doesn't construct the list in memory - faster!
    pass

y = ['a','b','c','d','e']
x = range(len(y))
print "x =", x
print "y =", y
print zip(x,y)

for x_val, y_val in zip(x,y):
    print "idx", x_val, "=", y_val
    
my_tuple = (1, 'two', 3)
print my_tuple
print my_tuple[1]

#%% Dictionaries

my_name = "Charles"
my_age = 3.2

about_me = {'name' : my_name, 'age' : my_age, 'height' : "5'8"}

about_me['height'] = 3.2
about_me['phobia'] = "Monkfish"
print about_me['height']

#%%
for k,v in about_me.iteritems():
    print "key=", k, " val = ", v
    
if 'hair_color' in about_me:
    print "Yes."
else:
    print "No."
    

#%% Functions

def PrintHello():
    print "Wassup."
    
PrintHello()

def MyPrint(arg):
    print arg
    print arg
    
MyPrint("Hi")

def MyPrint(arg="Hello"):
    print arg
    
MyPrint()
MyPrint(range(4))

def fib(n=5):
    """Get a Fibonacci series up to n."""
    a, b = 0, 1
    series = [a]
    while b < n:
        a, b = b, a+b
        series.append(a)
    return series
    
print fib(10)

multi_line_str = """This is the first line
This is the second,
and a third.""" 
# The 3 quotation marks allow text to continue
# onto further lines below.

print multi_line_str

#%% Classes
class Contact(object):
    """A given person for my database of dudes."""
    
    def __init__(self,firstName=None, lastName=None, email=None, phone=None):
        self.firstName = firstName
        self.lastName = lastName
        self.email = email
        self.phone = phone
        
    def printInfo(self):
        """Print all of the information of this contact."""
        myStr = "Contact info:"
        if self.firstName:
            myStr += " " + self.firstName
        if self.lastName:
            myStr += " " + self.lastName
        if self.email:
            myStr += " " + self.email
        if self.phone:
            myStr += " " + self.phone
        
        print myStr
        
#%%
    
bob = Contact('Bob','Smith')
joe = Contact(email = 'someone@somewhere.com')
joe.firstName = "Joe"

name = joe.firstName
print name

joe.printInfo()
help(Contact)

#%% Importing modules

from itertools import izip
y = ['a', 'b', 'c', 'd', 'e']
x = range(len(y))
for (x_val,y_val) in izip(x, y):
    print "idx", x_val, "=", y_val

import numpy
myVec= numpy.arange(0,1,0.1)
print myVec

#%% Pickling - save a Python object for later restoration with a load command

import pickle
contacts = [joe,bob]

with open('contacts.p','w') as pickle_file: # Make a new file
    pickle.dump(contacts, pickle_file) # Write contact list
    
with open('contacts.p', 'r') as pickle_file: # Open the file for reading
    contacts2 = pickle.load(pickle_file) # Load the pickled contents
    
for elem in contacts2:
    elem.printInfo()
