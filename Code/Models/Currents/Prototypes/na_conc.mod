
COMMENT
Longitudinal diffusion of sodium (no buffering)
(equivalent modified euler with standard method and
equivalent to diagonalized linear solver with CVODE )
ENDCOMMENT

NEURON {
   SUFFIX na_conc
   USEION na READ ina WRITE nai
}

UNITS {
   (mM) = (milli/liter)
   (um) = (micron)
   FARADAY = (faraday) (coulomb)
   PI = (pi) (1)
}

ASSIGNED {
   ina (milliamp/cm2)
   diam (um)
}

STATE {
   nai (mM)
}

BREAKPOINT {
   SOLVE conc METHOD sparse
}

KINETIC conc {
   COMPARTMENT PI*diam*diam/4 {nai}
   ~ nai << (-ina/(FARADAY)*PI*diam*(1e4))
}