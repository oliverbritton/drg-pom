TITLE shift_hh_6.mod  
 
COMMENT
shift HH in 6 mV to the right. 
v1/2 of minf is -34.12 and the slope factor is 9.146. 
v1/2 of hinf is -56.39 and the slope factor is 7.22.
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX Shh_6
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk, ina
        GLOBAL minf, hinf, ninf, mtau, htau, ntau
}
 
PARAMETER {
        gnabar = 0 (S/cm2)	<0,1e9>
        gkbar = 0 (S/cm2)	<0,1e9>
        gl = 0 (S/cm2)	<0,1e9>
        el = 0 (mV)
}
 
STATE {
        m h n
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

	gna (S/cm2)
	gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf
	mtau (ms) htau (ms) ntau (ms)
}
 
LOCAL mexp, hexp, nexp        
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
	ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
	ik = gk*(v - ek)      
        il = gl*(v - el)
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}
 
LOCAL q10


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
        TABLE minf, mtau, hinf, htau, ninf, ntau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)
                :"m" sodium activation system
        alpha = .1 * vtrap(-(v+34),10)
        beta =  4 * exp(-(v+59)/18)
        sum = alpha + beta
	mtau = 1/(q10*sum)
        minf = alpha/sum
                :"h" sodium inactivation system
        alpha = .07 * exp(-(v+59)/20)
        beta = 1 / (exp(-(v+29)/10) + 1)
        sum = alpha + beta
	htau = 1/(q10*sum)
        hinf = alpha/sum
                :"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10) 
        beta = .125*exp(-(v+65)/80)
	sum = alpha + beta
        ntau = 1/(q10*sum)
        ninf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
