TITLE Nav 1.8 from Tigerholm

COMMENT
INav1.8 from Tigerholm et al. 2014 dervied from Sheets et al. 2007

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX nav18tf
		 USEION na READ ena WRITE ina
		 RANGE gbar, gna, ina	
	     GLOBAL minf, hinf, sinf, uinf, mtau, htau, stau, utau
}

PARAMETER {
		 gbar = 0.2427124 (S/cm2) <0,1e9>		 
}

STATE {
		 m h s u
}

ASSIGNED {
		 v (mV)
		 celsius (degC)
		 ena (mV)
		 
		 gna (S/cm2)
		 ina (mA/cm2)
		 minf hinf sinf uinf
		 mtau (ms) htau (ms) stau (ms) utau (ms)
}

LOCAL mexp, hexp, sexp, uexp

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gbar*m*m*m*h*s*u
		ina = gna*(v - ena)
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
	s = sinf
	u = uinf
}

? states
DERIVATIVE states {
		rates(v)
		m' = (minf-m)/mtau
		h' = (hinf-h)/htau
		s' = (sinf-s)/stau
		u' = (uinf-u)/utau
}

:LOCAL q10

? rates
PROCEDURE rates(v(mV)) { : Computes rate and other constants at current v.
						 : Call once from HOC to initialize inf at resting v.
						 LOCAL alpha_m, beta_m, alpha_s, beta_s, alpha_u, beta_u
						 TABLE minf, mtau, hinf, htau, sinf, stau, uinf, utau DEPEND celsius FROM -100 TO 100 WITH 200
						 
UNITSOFF
		alpha_m = 2.85 - 2.839/(1 + exp((v-1.159)/13.95))
		beta_m = 7.6205/(1 + exp((v+46.463)/8.8289))
		
		minf = alpha_m/(alpha_m + beta_m)
		mtau = 1/(alpha_m + beta_m)
		
		hinf = 1/(1+exp((v+32.2)/4))
		htau = 1.218 + 42.043*exp(-((v+38.1)^2)/(2*15.19^2))	

		alpha_s = 0.001 * 5.4203 / (1 + exp((v+79.816)/16.269))
		beta_s  = 0.001 * 5.0757 / (1 + exp(-(v+15.968)/11.542))
		
		sinf = 1/(1+exp((v+45.0)/8))
		stau = 1/(alpha_s + beta_s)
		
		alpha_u =  0.002 * 2.0434 / (1 + exp((v+67.499)/19.51))
		beta_u =  0.002 * 1.9952 / (1 + exp(-(v+30.963)/14.792))
		
		uinf = 1/(1+exp((v+51.0)/8))
		utau = 1/(alpha_s + beta_s) 
}
UNITSON
