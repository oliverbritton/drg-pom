TITLE KDR

: v1/2 of minf is from Rola and the k is from Fedulova. mtau from Fedulova. The VC fits to Everill's
: voltage clamp recordings (fig 6/4).

NEURON {
	SUFFIX kdr7
	USEION k READ ek WRITE ik
	RANGE  gbar
	GLOBAL minf, mtau
}

PARAMETER {
	gbar = 0.0026   (mho/cm2)							
	celsius
	ek		(mV)
	v 		(mV)	
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 	  (mA/cm2)
	minf mtau (ms)	 	
}
 

STATE { m }

BREAKPOINT {
        SOLVE states METHOD cnexp
	ik = gbar*m*(v - ek)
} 

INITIAL {
	trates(v)
	m=minf   
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
}

PROCEDURE trates(v) {         
        minf = 1/(1 + exp((v + 13.2)/-12.5))         
	mtau = -23 + 69.46*exp(-0.0142*v) 
	:mtau = 2.04*exp(-0.1519*v) + 18*exp(-0.007095*v) :fast avrg
}





