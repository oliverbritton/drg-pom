TITLE Human Nav 1.8 from Han/Waxman et al.

COMMENT
INav1.8 (human) from Han, Waxman et al. 2015 
Paper: Human Nav1.8: enhanced persistent and ramp currents contribute to distinct
firing properties of human DRG neurons
Chongyang Han, Mark Estacion, Jianying Huang, Dymtro Vasylyev, Peng Zhao,
Sulayman D. Dib-Hajj, and Stephen G. Waxman

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX nav18hw
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
        gna = gbar*m*m*m*h
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
						 LOCAL alpha_m, beta_m, alpha_h, beta_h
						 TABLE minf, mtau, hinf, htau DEPEND celsius FROM -100 TO 100 WITH 200
						 
UNITSOFF
		alpha_m = 7.35 - 7.35/(1 + exp((v+1.38)/10.9))
		beta_m = 5.97/(1 + exp((v+56.43)/18.26))
		
		minf = alpha_m/(alpha_m + beta_m)
		mtau = 1/(alpha_m + beta_m)
		
		alpha_h = 0.011 + 1.39/(1 + exp((v+78.04)/11.32))
		beta_h = 0.56 - 0.56/(1 + exp((v-21.82)/20.03))
		
		hinf = alpha_h/(alpha_h + beta_h)
		htau = 1/(alpha_h + beta_h)		
}
UNITSON
