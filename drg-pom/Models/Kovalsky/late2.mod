TITLE late2.mod   ttx-sensitive late sodium channel
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
NEURON {
        SUFFIX ls :'l' for late, 's' for sodium    
        USEION na READ ena WRITE ina
        RANGE gnabar, gna, ina
        GLOBAL minf, hinf, htau
}
 
PARAMETER {
        gnabar = 0 (S/cm2)	<0,1e9>              
}
 
STATE {
        m h
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        
	gna (S/cm2)	
	ina (mA/cm2)
        minf hinf 
	htau (ms)
}
 
LOCAL mexp, hexp
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        m = minf
        gna = gnabar*m*h
	ina = gna*(v - ena)	
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf	
}

? states
DERIVATIVE states {  
        rates(v)
        h' = (hinf-h)/htau       
}

LOCAL q10

? rates
PROCEDURE rates(v(mV)) {  
        TABLE minf, hinf, htau DEPEND celsius FROM -100 TO 100 WITH 200
UNITSOFF
	        :"m" sodium activation system taken from HH
	minf = 1/(1+exp((-51.8-v)/4.6)) :Baker 2000
                :"h" sodium inactivation system
	htau = 1/(0.04 * exp(v/25.5)) + 63.2            :fig 3C
	hinf = 0.9827/(1 + exp(-((v + 55.67)/-6.552)))  :fig 3E 
}
UNITSON
