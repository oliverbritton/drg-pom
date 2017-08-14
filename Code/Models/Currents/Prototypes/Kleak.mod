TITLE IKleak

COMMENT
Simple potassium leak channel

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX kleak
		 USEION k READ ek WRITE ik
		 RANGE gbar, gk, ik	

}

PARAMETER {
		 gbar = 0.001 (S/cm2) <0,1e9>		 
}

ASSIGNED {
		 v (mV)
		 celsius (degC)
		 ek (mV)
		 
		 ik (mA/cm2)
}

? currents
BREAKPOINT {
		ik = gbar*(v - ek)
}







