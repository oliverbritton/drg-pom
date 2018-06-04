: Ca-dependent K channels (BK and SK)


NEURON {
	SUFFIX CadepK
	USEION ca READ ica
	USEION k READ ek WRITE ik
	RANGE gbkbar, gskbar
	GLOBAL ca0, tau, stau
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	B = .26 (mM-cm2/mA-ms)
}

PARAMETER {
	gbkbar = .01	(S/cm2)	: maximum permeability
	gskbar = .01	(S/cm2)	: maximum permeability
	ca0 = .00007	(mM)
	tau = 9		(ms)
	alphar = 7.5	(/ms)
	stau = 10		(ms)
}

ASSIGNED {
	v		(mV)
	ek		(mV)
	ik		(mA/cm2)
	ica		(mA/cm2)
	area		(microm2)
      gbk		(S/cm2)
      gsk		(S/cm2)
}

STATE { ca_i (mM) q r s }

BREAKPOINT {
	SOLVE state METHOD cnexp
	gbk = gbkbar*r*s*s
	gsk = gskbar*q*q
	ik = (gbk+gsk)*(v - ek)
}

DERIVATIVE state {	: exact when v held constant; integrates over dt step
	ca_i' = -B*ica-(ca_i-ca0)/tau
	q' = alphaq(ca_i)*(1-q)-betaq(ca_i)*q
	r' = alphar*(1-r)-betar(v)*r
	s' = (sinf(ca_i)-s)/stau
}

INITIAL {
	ca_i = ca0
	q = alphaq(ca_i)/(alphaq(ca_i)+betaq(ca_i))
	r = alphar/(alphar+betar(v))
      s = sinf(ca_i)
}

FUNCTION exp1(A (/ms), d, k, x (mM)) (/ms) {
	UNITSOFF
	exp1 = A/exp((12*log10(x)+d)/k)
	UNITSON
}

FUNCTION alphaq(x (mM)) (/ms) {
	alphaq = exp1(0.00246,28.48,-4.5,x)
}

FUNCTION betaq(x (mM)) (/ms) {
	betaq = exp1(0.006,60.4,35,x)
}

FUNCTION betar(v (mV)) (/ms) {
	UNITSOFF
	betar = 0.11/exp((v-35)/14.9)
	UNITSON
}

FUNCTION sinf(x (mM)) {
	UNITSOFF
	sinf = 1/(1+4/(1000*x))
	UNITSON
}