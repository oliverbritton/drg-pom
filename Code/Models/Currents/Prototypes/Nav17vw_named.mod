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
		 SUFFIX nav17vw_named
		 USEION na READ ena WRITE ina
		 RANGE gbar, gna, ina, m_alpha_independentrate, m_alpha_rate,	m_alpha_vhalf, m_alpha_tau, m_beta_independentrate, m_beta_rate, m_beta_vhalf, m_beta_tau, h_alpha_independentrate, h_alpha_rate, h_alpha_vhalf, h_alpha_tau, h_beta_independentrate, h_beta_rate, h_beta_vhalf, h_beta_tau   	
	     GLOBAL minf, hinf, mtau, htau
}

PARAMETER {
		 gbar = 0.18 (S/cm2) <0,1e9>		

		: m gate parameters
		m_alpha_independentrate = 0.0
		m_alpha_rate = 10.22
		m_alpha_vhalf = 7.19
		m_alpha_tau = 15.43
		
		m_beta_independentrate = 0.0
		m_beta_rate = 23.76
		m_beta_vhalf = 70.37
		m_beta_tau = 14.53
		
		: h gate parameters
		h_alpha_independentrate = 0.0
		h_alpha_rate = 0.0744
		h_alpha_vhalf = 99.76
		h_alpha_tau = 11.07
	
		h_beta_independentrate = 0.0
		h_beta_rate = 2.54
		h_beta_vhalf = 7.8
		h_beta_tau = 10.68
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
		alpha_m = m_alpha_rate*(1.0 - 1.0/(1 + exp((v+m_alpha_vhalf)/m_alpha_tau)))
		beta_m = m_beta_rate/(1 + exp((v+m_beta_vhalf)/m_beta_tau))
		
		minf = alpha_m/(alpha_m + beta_m)
		mtau = 1/(alpha_m + beta_m)
		
		alpha_h = h_alpha_rate/(1 + exp((v+h_alpha_vhalf)/h_alpha_tau))
		beta_h = h_beta_rate*(1.0 - 1.0/(1 + exp((v+h_beta_vhalf)/h_beta_tau)))
		
		hinf = alpha_h/(alpha_h + beta_h)
		htau = 1/(alpha_h + beta_h)
		
}
UNITSON
