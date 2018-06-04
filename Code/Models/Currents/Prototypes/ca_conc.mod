
COMMENT
Concentration of calcium (no buffering)
(equivalent modified euler with standard method and
equivalent to diagonalized linear solver with CVODE )
ENDCOMMENT

NEURON {
   SUFFIX ca_conc
   USEION ca READ ica WRITE cai
}

UNITS {
   (mM) = (milli/liter)
   (um) = (micron)
   FARADAY = (faraday) (coulomb)
   PI = (pi) (1)
}

ASSIGNED {
   ica (milliamp/cm2)
   diam (um)
}

STATE {
   cai (mM)
}

BREAKPOINT {
   SOLVE conc METHOD sparse
}

KINETIC conc {
   COMPARTMENT PI*diam*diam/4 {cai}
   ~ cai << (-ica/(FARADAY)*PI*diam*(1e4))
}