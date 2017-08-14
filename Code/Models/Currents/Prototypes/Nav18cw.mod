TITLE Nav 1.8 from Choi/Waxman

COMMENT
INav1.8 from Choi and Waxman 2011

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX nav18cw
		 USEION na READ ena WRITE ina
		 RANGE gbar, gna, ina	
	     GLOBAL minf, hinf, mtau, htau
}

PARAMETER {
		 gbar = 0.026 (S/cm2) <0,1e9>		 
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
		 mtau (ms) htau (ms)
}

LOCAL mexp, hexp

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gbar*m*h
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
		m' = (minf-m)/mtau
		h' = (hinf-h)/htau
}

:LOCAL q10

? rates
PROCEDURE rates(v(mV)) { : Computes rate and other constants at current v.
						 : Call once from HOC to initialize inf at resting v.
						 LOCAL alpha_m, beta_m
						 TABLE minf, mtau, hinf, htau DEPEND celsius FROM -100 TO 100 WITH 200
						 
UNITSOFF
		alpha_m = 2.85 - 2.839/(1 + exp((v-1.159)/13.95))
		beta_m = 7.6205/(1 + exp((v+46.463)/8.8289))
		
		minf = alpha_m/(alpha_m + beta_m)
		mtau = 1/(alpha_m + beta_m)
		
		hinf = 1/(1+exp((v+32.2)/4))
		htau = 1.218 + 42.043*exp(-((v+38.1)^2)/(2*15.19^2))		
}
UNITSON
