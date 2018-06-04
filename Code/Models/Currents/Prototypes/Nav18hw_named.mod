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
		 SUFFIX nav18hw_named
		 USEION na READ ena WRITE ina
		 RANGE gbar, gna, ina, m_alpha_independentrate, m_alpha_rate,	m_alpha_vhalf, m_alpha_tau, m_beta_independentrate, m_beta_rate, m_beta_vhalf, m_beta_tau, h_alpha_independentrate, h_alpha_rate, h_alpha_vhalf, h_alpha_tau, h_beta_independentrate, h_beta_rate, h_beta_vhalf, h_beta_tau   
	     GLOBAL minf, hinf, mtau, htau
}

PARAMETER {
		 gbar = 0.026 (S/cm2) <0,1e9>		 
		 
		: m gate parameters
		m_alpha_independentrate = 0.0
		m_alpha_rate = 7.35
		m_alpha_vhalf = 1.38
		m_alpha_tau = 10.9 
		
		m_beta_independentrate = 0.0
		m_beta_rate = 5.97
		m_beta_vhalf = 56.43 
		m_beta_tau = 18.26 
		
		:alpha_m = 7.35 - 7.35/(1 + exp((v+1.38)/10.9))
		:beta_m = 5.97/(1 + exp((v+56.43)/18.26))
		
		:alpha_m = m_alpha_independentrate - m_alpha_rate/(1 + exp((v+m_alpha_vhalf)/m_alpha_tau))
		:beta_m = m_beta_rate/(1 + exp((v+m_beta_vhalf)/m_beta_tau))
		 
		 : h gate parameters
		h_alpha_independentrate = 0.011
		h_alpha_rate = 1.39
		h_alpha_vhalf = 78.04 
		h_alpha_tau = 11.32
	
		h_beta_independentrate = 0.0
		h_beta_rate = 0.56
		h_beta_vhalf = -21.82 
		h_beta_tau = 20.03 
		 
		:alpha_h = 0.011 + 1.39/(1 + exp((v+78.04)/11.32))
		:beta_h = 0.56 - 0.56/(1 + exp((v-21.82)/20.03))
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
						 :TABLE minf, mtau, hinf, htau DEPEND celsius FROM -100 TO 100 WITH 200
						 
UNITSOFF
		alpha_m = m_alpha_rate*(1.0 - 1.0/(1 + exp((v+m_alpha_vhalf)/m_alpha_tau)))
		beta_m = m_beta_rate/(1 + exp((v+m_beta_vhalf)/m_beta_tau))
		
		minf = alpha_m/(alpha_m + beta_m)
		mtau = 1/(alpha_m + beta_m)
		
		alpha_h = h_alpha_independentrate + h_alpha_rate/(1 + exp((v+h_alpha_vhalf)/h_alpha_tau))
		beta_h = h_beta_rate*(1.0 - 1.0/(1 + exp((v+h_beta_vhalf )/h_beta_tau)))
		
		hinf = alpha_h/(alpha_h + beta_h)
		htau = 1/(alpha_h + beta_h)		
}
UNITSON
