
COMMENT
Concentration of potassium (no buffering)
(equivalent modified euler with standard method and
equivalent to diagonalized linear solver with CVODE )
ENDCOMMENT

NEURON {
   SUFFIX k_conc
   USEION k READ ik WRITE ki
}

UNITS {
   (mM) = (milli/liter)
   (um) = (micron)
   FARADAY = (faraday) (coulomb)
   PI = (pi) (1)
}

ASSIGNED {
   ik (milliamp/cm2)
   diam (um)
}

STATE {
   ki (mM)
}

BREAKPOINT {
   SOLVE conc METHOD sparse
}

KINETIC conc {
   COMPARTMENT PI*diam*diam/4 {ki}
   ~ ki << (-ik/(FARADAY)*PI*diam*(1e4))
}