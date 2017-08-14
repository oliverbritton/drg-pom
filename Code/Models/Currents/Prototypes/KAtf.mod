TITLE IKA from Tigerholm

COMMENT
IKA from Tigerholm 2014

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX katf
		 USEION k READ ek WRITE ik
		 RANGE gbar, gk, ik	
	     GLOBAL ninf, hinf, ntau, htau
}

PARAMETER {
		 gbar = 0.001275 (S/cm2) <0,1e9>		 
}

STATE {
		n h
}

ASSIGNED {
		 v (mV)
		 celsius (degC)
		 ek (mV)
		 
		 gk (S/cm2)
		 ik (mA/cm2)
		 ninf hinf
		 ntau (ms) htau (ms)
}

LOCAL nexp, hexp

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gk = gbar*n*h
		ik = gk*(v - ek)
}

INITIAL {
	rates(v)
	n = ninf
	h = hinf
}

? states
DERIVATIVE states {
		rates(v)
		n' = (ninf-n)/ntau
		h' = (hinf-h)/htau
}

LOCAL q10

? rates
PROCEDURE rates(v(mV)) { : Computes rate and other constants at current v.
						 : Call once from HOC to initialize inf at resting v.
						 TABLE ninf, ntau, hinf, htau DEPEND celsius FROM -100 TO 100 WITH 200
						 
UNITSOFF
		
		q10 = 3.3
		ninf = (1/(1 + exp(-(v+5.4+15)/16.4)))^4
		ntau = 0.25 + 10.04*exp((-(v+24.67)^2)/(2*34.8^2))*q10
		
		hinf = 1/(1 + exp((v+49.9 + 15)/4.6))
		htau = 20 + 50*exp((-(v+40)^2)/(2*40^2))*q10
		htau = htautrap(htau)
}
UNITSON


FUNCTION htautrap(x) {  : Trap for htau following Sheets /ChoiWaxman/Tigerholm- set it to 5 ms if less than 5 ms
		if (x < 5) {
				htautrap = 5
		}else{
				htautrap = x
		}
}

