TITLE Ih from Kouranova et al. 2008

COMMENT
Ih from Kouranova et al. 2008
Recorded at 22

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX hcnkn
		 NONSPECIFIC_CURRENT i
		 RANGE i, gbar, gs, gf	, eh
		 GLOBAL nsinf, nfinf, nstau, nftau
}

PARAMETER {
		 gbar = 0.001 (S/cm2) <0,1e9>		
		 gs = 2.0
		 gf = 1.0	
		 eh = -30.0 (mV)
}

STATE {
		ns nf
}

ASSIGNED {
		 v (mV)
		 celsius (degC)

		 i (mA/cm2)
		 nsinf
		 nfinf
		 nstau (ms)
		 nftau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
		i = gbar*(gf*nf + gs*ns)*(v - eh)
}

INITIAL {
	rates(v)
	ns = 0
	nf = 0
}

? states
DERIVATIVE states {
		rates(v)
		ns' = (nsinf-ns)/nstau
		nf' = (nfinf-nf)/nftau
}

:LOCAL q10

? rates
PROCEDURE rates(v(mV)) { : Computes rate and other constants at current v.
						 : Call once from HOC to initialize inf at resting v.

						 TABLE nsinf, nstau, nfinf, nftau DEPEND celsius FROM -100 TO 100 WITH 200
						 
UNITSOFF
		
		nsinf = 1/(1 + exp((v+87.2)/9.7))
		nfinf = nsinf
		
		if (v > -70.0 ) {
		nstau = 300.0 + 542.0 * exp((v+25.0)/20.0) 
		nftau = 140.0 + 50.0 * exp(-(v+25.0)/20.0) 
		}else{
			nstau = 2500.0 + 100.0 * exp((v+240.0)/50.0)
			nftau = 250.0 + 12.0 * exp((v+240.0)/50.0)
		}
		
}
UNITSON

