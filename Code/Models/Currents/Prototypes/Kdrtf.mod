TITLE IKdr from Tigerholm 2014

COMMENT
IKdr from Tigerholm 2014

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX kdrtf
		 USEION k READ ek WRITE ik
		 RANGE gbar, gk, ik	
	     GLOBAL ninf, ntau
}

PARAMETER {
		 gbar = 0.0018 (S/cm2) <0,1e9>		 
}

STATE {
		n
}

ASSIGNED {
		 v (mV)
		 celsius (degC)
		 ek (mV)
		 
		 gk (S/cm2)
		 ik (mA/cm2)
		 ninf
		 ntau (ms)
}

LOCAL nexp

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gk = gbar*n*n*n*n
		ik = gk*(v - ek)
}

INITIAL {
	rates(v)
	n = ninf
}

? states
DERIVATIVE states {
		rates(v)
		n' = (ninf-n)/ntau
}

LOCAL q10

? rates
PROCEDURE rates(v(mV)) { : Computes rate and other constants at current v.
						 : Call once from HOC to initialize inf at resting v.
						 TABLE ninf, ntau DEPEND celsius FROM -100 TO 100 WITH 200
						 
UNITSOFF
		q10 = 1.0 : 3.3 Disable q10 for now
		ninf = 1/(1 + exp(-(v+45)/15.4))
		ntau = taucalc(v)*q10 : needs an extra factor of ^(deltaT/10) to get right q10
}
UNITSON

FUNCTION taucalc(x) {  : Equation for ntau following Sheets/Tigerholm - use
		if (x > -31) {
				 taucalc =  0.16+0.8*exp(-0.0267*(x+11))
		}else{
				taucalc = 1000*(0.000688 + 1/(exp((x+75.2)/6.5) + exp(-(x-131.5)/(34.8))))
		}
}
