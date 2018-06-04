TITLE Borg-Graham type generic K-A channel for a Sympathetic Preganglionic Neuron

COMMENT
	Description: A-type transient K current for a Sympathetic Preganglionic Neuron.	
	Author: Linford Briant
	
	A-type transient K current = "IA"
	
	Sympathetic Preganglionic Neurones = "SPNs"
	
	Note that this is a modified version of IA found widely on SenseLab. This version has had steady-state kinetics 
	that have been fit to data for the IA in SPN according to Whyment et al. (2011).
	
	Whyment et al. (2011), PMID: 21211550

ENDCOMMENT



UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 		(mV)
	ek 		(mV)
	celsius 	(degC)
	gkabar=0.012 	(mho/cm2)
	vhalfn=-45	(mV)
	vhalfl=-67	(mV)
	vhalfm=-67	(mV)
	vhalfk=-45	(mV)
	a0l=0.023	(/ms)
	a0n=0.04	(/ms)
	zetan=-4	(1)
	zetal=2    	(1)
	gmn=0.45   	(1)
	gml=1 	  	(1)
	zetam=4		(1)
	zetak=-5	(1)
}


NEURON {
	SUFFIX borgka
	USEION k READ ek WRITE ik
        RANGE gkabar,gka,vhalfn,vhalfl,a0l,a0n,zetan,zetal,gmn,gml,zetam,zetak,vhalfm,vhalfk
        GLOBAL ninf,linf,taul,taun
}

STATE {
	n
	l
}

INITIAL {
        rates(v)
        n=ninf
        l=linf
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        linf      
        taul
        taun
        gka
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gka = gkabar*n*l
	ik = gka*(v-ek)
}


FUNCTION alpn(v(mV)) {
  alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alpk(v(mV)) {
  alpk = exp(1.e-3*zetak*(v-vhalfk)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV)) {
  betn = exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alpl(v(mV)) {
  alpl = exp(1.e-3*zetal*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alpm(v(mV)) {
  alpm = exp(1.e-3*zetam*(v-vhalfm)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betl(v(mV)) {
  betl = exp(1.e-3*zetal*gml*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE states { 
        rates(v)
        n' = (ninf - n)/taun
        l' = (linf - l)/taul
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,q10,b
        q10=3^((celsius-30)/10)
        a = alpn(v)
	b = alpk(v)
        ninf = 1/(1 + b)
        taun = betn(v)/(q10*a0n*(1 + a))
        a = alpl(v)
	b = alpm(v)
        linf = 1/(1 + b)
        taul = betl(v)/(q10*a0l*(1 + a))
}
