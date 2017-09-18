TITLE Nav 1.7 from Vasylyev, Waxman et al. 2014

COMMENT
INav1.7 from Vasylyev, Waxman et al. 2014

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX nav17vw
		 USEION na READ ena WRITE ina
		 RANGE gbar, gna, ina	
	     GLOBAL minf, hinf, mtau, htau
}

PARAMETER {
		 gbar = 0.18 (S/cm2) <0,1e9>		 
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
		alpha_m = 10.22 - 10.22/(1 + exp((v+7.19)/15.43))
		beta_m = 23.76/(1 + exp((v+70.37)/14.53))
		
		minf = alpha_m/(alpha_m + beta_m)
		mtau = 1/(alpha_m + beta_m)
		
		alpha_h = 0.0744/(1 + exp((v+99.76)/11.07))
		beta_h = 2.54 - 2.54/(1 + exp((v+7.8)/10.68))
		
		hinf = alpha_h/(alpha_h + beta_h)
		htau = 1/(alpha_h + beta_h)
		
}
UNITSON
