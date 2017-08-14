TITLE Nav 1.7 from Choi/Waxman

COMMENT
INav1.7 from Choi and Waxman 2011

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX nav17cw
		 USEION na READ ena WRITE ina
		 RANGE gbar, gna, ina	
	     GLOBAL minf, hinf, sinf, mtau, htau, stau
}

PARAMETER {
		 gbar = 0.018 (S/cm2) <0,1e9>		 
}

STATE {
		 m h s
}

ASSIGNED {
		 v (mV)
		 celsius (degC)
		 ena (mV)
		 
		 gna (S/cm2)
		 ina (mA/cm2)
		 minf hinf sinf
		 mtau (ms) htau (ms) stau (ms)
}

LOCAL mexp, hexp, sexp

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gbar*m*m*m*h*s
		ina = gna*(v - ena)
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
	s = sinf
}

? states
DERIVATIVE states {
		rates(v)
		m' = (minf-m)/mtau
		h' = (hinf-h)/htau
		s' = (sinf-s)/stau
}

:LOCAL q10

? rates
PROCEDURE rates(v(mV)) { : Computes rate and other constants at current v.
						 : Call once from HOC to initialize inf at resting v.
						 LOCAL alpha_m, beta_m, alpha_h, beta_h, alpha_s, beta_s
						 TABLE minf, mtau, hinf, htau, sinf, stau DEPEND celsius FROM -100 TO 100 WITH 200
						 
UNITSOFF
		alpha_m = 15.5/(1 + exp(-(v-5)/(12.08)))
		beta_m = 35.2/(1 + exp((v+72.7)/16.7))
		
		minf = alpha_m/(alpha_m + beta_m)
		mtau = 1/(alpha_m + beta_m)
		
		alpha_h = 0.38685/(1 + exp((v+122.35)/15.29))
		beta_h = -0.00283 + 2.00283/(1 + exp(-(v+5.5266)/12.70195))
		
		hinf = alpha_h/(alpha_h + beta_h)
		htau = 1/(alpha_h + beta_h)
		
		alpha_s = 0.00003 + 0.00092/(1 + exp((v+93.9)/16.6))
		beta_s = 132.05 - 132.05/(1 + exp((v-384.9)/28.5))
		
		sinf = alpha_s/(alpha_s + beta_s)
		stau = 1/(alpha_s + beta_s)
}
UNITSON
