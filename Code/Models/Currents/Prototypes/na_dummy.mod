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
		 SUFFIX na_dummy
		 USEION na READ ena WRITE ina
}

PARAMETER {
}

ASSIGNED {
		 v (mV)
		 celsius (degC)
		 ena (mV)
		 
		 ina (mA/cm2)
}

? currents
BREAKPOINT {
}







