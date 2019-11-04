#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _inter_reg();
extern void _kdr7_reg();
extern void _late2_reg();
extern void _shift_hh_6_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," inter.mod");
fprintf(stderr," kdr7.mod");
fprintf(stderr," late2.mod");
fprintf(stderr," shift_hh_6.mod");
fprintf(stderr, "\n");
    }
_inter_reg();
_kdr7_reg();
_late2_reg();
_shift_hh_6_reg();
}
