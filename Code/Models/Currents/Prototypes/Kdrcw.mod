TITLE IKdr from Choi/Waxman

COMMENT
IKdr from Choi and Waxman 2011
Recorded at 22

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX kdrcw
		 USEION k READ ek WRITE ik
		 RANGE gbar, gk, ik	
	     GLOBAL ninf, ntau
}

PARAMETER {
		 gbar = 0.0035 (S/cm2) <0,1e9>		 
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
        gk = gbar*n
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

:LOCAL q10

? rates
PROCEDURE rates(v(mV)) { : Computes rate and other constants at current v.
						 : Call once from HOC to initialize inf at resting v.
						 LOCAL alpha_n, beta_n
						 TABLE ninf, ntau DEPEND celsius FROM -100 TO 100 WITH 200
						 
UNITSOFF
		
		alpha_n = 0.001265*(v+14.273)/(1-exp(-(v+14.273)/10))
		beta_n = 0.125*exp(-(v+55)/2.5)
		
		ninf = alpha_n/(alpha_n + beta_n)
		ntau = 1/(alpha_n + beta_n)
}
UNITSON

