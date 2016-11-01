TITLE Test of intracellular sodium accumulation

COMMENT
Test whether we can write to nai

ENDCOMMENT

NEURON {
        SUFFIX naiTest
        USEION na WRITE nai
        RANGE nai0
}

PARAMETER {
        nai0 = 10 (milli/liter)
}

STATE {
        nai (milli/liter)
}


INITIAL {
        nai = nai0
}