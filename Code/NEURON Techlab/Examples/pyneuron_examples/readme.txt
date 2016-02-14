Welcome to the NEURON with Python introductory tutorial.

pyneuron_examples.zip contains code demonstrating basic use of NEURON
with the Python interpreter. 

model*.py was written by Ben Torben-Nielsen for use in the Okinawa Computational Neuroscience
Course, and later modified & presented by Sam Neymotin at the Organization for Computational
Neuroscience 2014 meeting's NEURON tutorial (organized by Bill Lytton). Errors in the final
code can be blamed on Sam.

NOTE: You will need to have NEURON installed with Python enabled ( http://www.neuron.yale.edu/neuron/ ).
For these examples, you also need to have Matplotlib installed ( http://matplotlib.org/ ).
These simulations were run/tested on Linux but may work on Windows.

Brief description of the py files contained, with the important commands/objects used:
 model1.py - import neuron, Section, insert, record, Vector, save output
 model2.py - change hh params, record sodium, potassium currents
 model3.py - IClamp inputs, play random values into IClamp (current clamp)
 hhtest.py - Another test of the hodgkin huxley mechanism with IClamp
 model4.py - Simple use of Synapse (Exp2Syn) with NetStim noisy inputs
 model4B.py - Create & connect 2 cells - drive 1 with NetStim
 model5.py - Connect 2 Sections together into a cell (soma+dendrite), vary position of synapse on
             dendrite to see how EPSPs and spikes change. Command-line arg of 0/1 toggles insertion
	     of an HH channel in the soma.

To run a particular simulation, first unzip pyneuron_examples.zip, then cd to the directory you just created.

You can run directly from the shell with:
 python -i file.py
The -i will keep the interpreter running after the file is executed.

Alternatively, you can start the python interpreter with:   python
Then, from the python interpreter you can load a particular simulation with:
from simfilename import *
Using import with * is generally not considered clean python practice, since you
generally want to keep your namespaces clean (e.g., in case you have two variables
with the same name in two .py files). 
Note that simfilename should not include the .py extension.

An alternative command to run the sim after python starts uses execfile: 
 execfile("simfilename.py")
Execfile executes all the code in the .py file. This allows you to have all the code/variables
in the simfilename.py file exist at the top/global-level namespace. This has some advantages if you
have a simulation organized over many .py files.

for questions/comments email samnemo at gmail dot com 

References:
 This paper is highly recommended for learning the basics of NEURON+Python:
  http://www.frontiersin.org/Neuroinformatics/10.3389/neuro.11.001.2009/full 
 and there's plenty of documentation at the official NEURON site:
  http://www.neuron.yale.edu/neuron/docs
 and info at the user's forum:
  http://www.neuron.yale.edu/phpBB/
