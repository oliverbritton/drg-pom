TITLE IKM from Tigerholm 2014

COMMENT
IKM from Tigerholm 2014

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX kmtf
		 USEION k READ ek WRITE ik
		 RANGE gbar, gk, ik	
	     GLOBAL ninf, nsTau, nfTau
}

PARAMETER {
		 gbar = 0.0018 (S/cm2) <0,1e9>		 
}

STATE {
		ns nf
}

ASSIGNED {
		 v (mV)
		 celsius (degC)
		 ek (mV)
		 
		 gk (S/cm2)
		 ik (mA/cm2)
		 ninf
		 nsTau (ms) nfTau (ms)
}

LOCAL nexp

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gk = gbar*(ns/4 + 3*nf/4)
		ik = gk*(v - ek)
}

INITIAL {
	rates(v)
	ns = ninf
	nf = ninf
}

? states
DERIVATIVE states {
		rates(v)
		ns' = (ninf-ns)/nsTau
		nf' = (ninf-nf)/nfTau
}

LOCAL q10

? rates
PROCEDURE rates(v(mV)) { : Computes rate and other constants at current v.
						 : Call once from HOC to initialize inf at resting v.
						 LOCAL nfTau_alpha, nfTau_beta
						 TABLE ninf, nsTau, nfTau DEPEND celsius FROM -100 TO 100 WITH 200
						 
UNITSOFF
		
		
		q10 = 3.3
		ninf = 1/(1 + exp(-(v+30)/6))
		nsTau = nsTauCalc(v,q10)
		
		nfTau_alpha = 0.00395*exp((v+30)/40)
		nfTau_beta = 0.00395*exp(-(v+30)/20)*q10
		nfTau = 1/(nfTau_alpha + nfTau_beta)
		
}
UNITSON

FUNCTION nsTauCalc(x,q10) {  : Equation for nsTau
		if (x < -60) {
				 nsTauCalc =  219*q10
		}else{
				nsTauCalc = 13*x + 1000*q10
		}
}
