TITLE Na dummy

COMMENT
Dummy to set up sodium conc

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX k_dummy
		 USEION k READ ek WRITE ik
}

PARAMETER {
}

ASSIGNED {
		 v (mV)
		 celsius (degC)
		 ek (mV)
		 
		 ik (mA/cm2)
}

? currents
BREAKPOINT {
}







