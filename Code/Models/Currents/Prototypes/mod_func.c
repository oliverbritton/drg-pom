#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _KAtf_reg();
extern void _KMtf_reg();
extern void _Kacw_reg();
extern void _Kdrcw_reg();
extern void _Kdrtf_reg();
extern void _Nav17cw_reg();
extern void _Nav17vw_reg();
extern void _Nav18cw_reg();
extern void _Nav18hw_reg();
extern void _Nav18tf_reg();
extern void _Nav19hw_reg();
extern void _Nav19tf_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," KAtf.mod");
fprintf(stderr," KMtf.mod");
fprintf(stderr," Kacw.mod");
fprintf(stderr," Kdrcw.mod");
fprintf(stderr," Kdrtf.mod");
fprintf(stderr," Nav17cw.mod");
fprintf(stderr," Nav17vw.mod");
fprintf(stderr," Nav18cw.mod");
fprintf(stderr," Nav18hw.mod");
fprintf(stderr," Nav18tf.mod");
fprintf(stderr," Nav19hw.mod");
fprintf(stderr," Nav19tf.mod");
fprintf(stderr, "\n");
    }
_KAtf_reg();
_KMtf_reg();
_Kacw_reg();
_Kdrcw_reg();
_Kdrtf_reg();
_Nav17cw_reg();
_Nav17vw_reg();
_Nav18cw_reg();
_Nav18hw_reg();
_Nav18tf_reg();
_Nav19hw_reg();
_Nav19tf_reg();
}
