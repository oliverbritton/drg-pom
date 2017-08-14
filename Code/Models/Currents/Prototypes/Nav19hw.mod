TITLE Nav 1.9 from Huang Waxman

COMMENT
INav1.9 from Huang et al. 2014 (Derived from human Nav 1.9 in rat)
Conductance is rounded up from Tigerholm

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX nav19hw
		 USEION na READ ena WRITE ina
		 RANGE gbar, gna, ina	
	     GLOBAL minf, hinf, sinf, mtau, htau, stau
}

PARAMETER {
		 gbar = 0.0001 (S/cm2) <0,1e9>		 
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
        gna = gbar*m*h*s
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
		alpha_m = 0.751/(1 + exp(-(v+32.26)/13.71))
		beta_m = 5.68/(1 + exp((v+123.71)/13.94))
		?hyperpolarisation of activation gate closing

		minf = alpha_m/(alpha_m + beta_m)
		mtau = 1/(alpha_m + beta_m)
		
		alpha_h = 0.082/(1 + exp((v+113.69)/17.4))
		beta_h = 0.24/(1 + exp(-(v-10.1)/17.2))
		
		hinf = alpha_h/(alpha_h + beta_h)
		htau = 1/(alpha_h + beta_h)
		
		alpha_s = 0.019/(1 + exp((v+154.51)/11.46))
		beta_s = 0.000376/(1 + exp(-(v+60.92)/15.79))
		
		sinf = alpha_s/(alpha_s + beta_s)
		stau = 1/(alpha_s + beta_s)
		
		
}
UNITSON
