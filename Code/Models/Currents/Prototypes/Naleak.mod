TITLE INaleak

COMMENT
Simple sodium leak channel

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX naleak
		 USEION na READ ena WRITE ina
		 RANGE gbar, ina	

}

PARAMETER {
		 gbar = 0.001 (S/cm2) <0,1e9>		 
}

ASSIGNED {
		 v (mV)
		 celsius (degC)
		 ena (mV)
		 
		 ina (mA/cm2)
}

? currents
BREAKPOINT {
		ina = gbar*(v - ena)
}







