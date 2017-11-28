Hippocampal CA3 pyramidal neuron model from the paper 
M. Migliore, E. Cook, D.B. Jaffe, D.A. Turner and D. Johnston, Computer
simulations of morphologically reconstructed CA3 hippocampal neurons, J.
Neurophysiol. 73, 1157-1168 (1995). 

The paper shows how bursting and non bursting firing modes of CA3 
neurons can depend on the densities of the Ca-independent conductances,
that prevent the build-up of the Ca-dependent depolarizing envelope
(compare bursting and short non bursting simulations).
Calcium accumulation during long current injections results in
spike adaptation modulated by Ca-dependent conductances and KM.

Under unix systems:
to compile the mod files use the command 
nrnivmodl 
and run the simulation hoc file with the command 
nrngui test_a.hoc

This will open a window from which three simulations using the A neuron
could be run to show bursting and non bursting firing of the same cell.


Under Windows:
to compile the mod files use the "mknrndll DOS box" and 
follow on-screen instructions.
A double click on the simulation file
test_a.hoc 
will open the simulation window.

Questions on how to use this model should be directed to
michele.migliore@pa.ibf.cnr.it

July 5th 2007 Model updated to run with Model View (1e-6 replaced 0's
in ca3a.geo) TMM
August 24th 2007 Updated to allow shape plots to be updated via
inclusion of flushPlot() in test_a.hoc TMM
