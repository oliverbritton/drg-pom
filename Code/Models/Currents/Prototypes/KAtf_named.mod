TITLE IKA from Tigerholm

COMMENT
IKA from Tigerholm 2014 with variable parameters

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX katf_named
		 USEION k READ ek WRITE ik
		 RANGE gbar, gk, ik, q10, n_inf_vhalf, n_inf_tau, n_tau_independentrate, n_tau_rate, n_tau_vhalf, n_tau_tau, h_inf_vhalf, h_inf_tau, h_tau_independentrate, h_tau_rate, h_tau_vhalf, h_tau_tau
	     GLOBAL ninf, hinf, ntau, htau
}

PARAMETER {
		 gbar = 0.001275 (S/cm2) <0,1e9>		 
		 q10 = 3.3
		 
		:n gate parameters
		 n_inf_vhalf = 20.4 :5.4 + 15
		 n_inf_tau = 16.4
		 n_tau_independentrate = 0.25
		 n_tau_rate = 10.04
		 n_tau_vhalf = 24.67
		 n_tau_tau = 2422.08 : 2*34.8^2
		 
		:h gate parameters
		h_inf_vhalf = 64.9 : 49.9 + 15
		h_inf_tau = 4.6
		h_tau_independentrate = 20
		h_tau_rate = 50
		h_tau_vhalf = 40
		h_tau_tau = 3200 :2*40^2
}

STATE {
		n h
}

ASSIGNED {
		 v (mV)
		 celsius (degC)
		 ek (mV)
		 
		 gk (S/cm2)
		 ik (mA/cm2)
		 ninf hinf
		 ntau (ms) htau (ms)
}

LOCAL nexp, hexp

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gk = gbar*n*h
		ik = gk*(v - ek)
}

INITIAL {
	rates(v)
	n = ninf
	h = hinf
}

? states
DERIVATIVE states {
		rates(v)
		n' = (ninf-n)/ntau
		h' = (hinf-h)/htau
}

: Remove any locals if converting future files to named
: Disable TABLE command

? rates
PROCEDURE rates(v(mV)) { : Computes rate and other constants at current v.
						 : Call once from HOC to initialize inf at resting v.
						 : TABLE ninf, ntau, hinf, htau DEPEND celsius FROM -100 TO 100 WITH 200
						 
UNITSOFF

		:Transformed n gate equations (leave n^4 dependency in)
		ninf = (1/(1 + exp(-(v+n_inf_vhalf)/n_inf_tau)))^4
		ntau = n_tau_independentrate + n_tau_rate*exp((-(v+n_tau_vhalf)^2)/(n_tau_tau))*q10


		
		:Transformed h gate equations 
		hinf = 1/(1 + exp((v+h_inf_vhalf)/h_inf_tau))
		htau = h_tau_independentrate + h_tau_rate*exp((-(v+h_tau_vhalf)^2)/(h_tau_tau))*q10
		htau = htautrap(htau)
		

}
UNITSON


FUNCTION htautrap(x) {  : Trap for htau following Sheets /ChoiWaxman/Tigerholm- set it to 5 ms if less than 5 ms
		if (x < 5) {
				htautrap = 5
		}else{
				htautrap = x
		}
}

