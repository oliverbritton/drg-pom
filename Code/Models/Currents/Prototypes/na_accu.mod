TITLE Constant sodium accumulator

COMMENT
Adds Na to cell at constant rate

ENDCOMMENT

UNITS {

		 (mA) = (milliamp)
		 (mV) = (millivolt)
		 (S) = (siemens)
}

NEURON {
		 SUFFIX na_accu
		 NONSPECIFIC_CURRENT i
		 USEION na READ ena WRITE ina
		 RANGE gna, ina, i
}

PARAMETER {
		 gna = 0.0001 (S/cm2) <0,1e9>		 
}

ASSIGNED {
		 v (mV)
		 ena (mV)
		 
		 ina (mA/cm2)
		 i (mA/cm2)
}


BREAKPOINT {
		ina = gna*(v - ena) : ena will be around 60 mV so (v-ena) will be -ve as inward current is negative.
		i = -ina : opposite to cancel out effect on vm
}



