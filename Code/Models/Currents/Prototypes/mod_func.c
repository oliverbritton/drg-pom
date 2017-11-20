#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _IRamp_reg();
extern void _KAtf_reg();
extern void _KMtf_reg();
extern void _Kacw_reg();
extern void _Kdrcw_reg();
extern void _Kdrtf_reg();
extern void _Kleak_reg();
extern void _Naleak_reg();
extern void _Nav16bk_reg();
extern void _Nav16hs_reg();
extern void _Nav17cw_reg();
extern void _Nav17vw_reg();
extern void _Nav17vw_named_reg();
extern void _Nav18cw_reg();
extern void _Nav18hw_reg();
extern void _Nav18hw_named_reg();
extern void _Nav18tf_reg();
extern void _Nav19hw_reg();
extern void _Nav19tf_reg();
extern void _hcn_kn_reg();
extern void _hcn_tf_reg();
extern void _k_conc_reg();
extern void _k_dummy_reg();
extern void _na_accu_reg();
extern void _na_conc_reg();
extern void _na_dummy_reg();
extern void _nadifl_reg();
extern void _naiTest_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," IRamp.mod");
fprintf(stderr," KAtf.mod");
fprintf(stderr," KMtf.mod");
fprintf(stderr," Kacw.mod");
fprintf(stderr," Kdrcw.mod");
fprintf(stderr," Kdrtf.mod");
fprintf(stderr," Kleak.mod");
fprintf(stderr," Naleak.mod");
fprintf(stderr," Nav16bk.mod");
fprintf(stderr," Nav16hs.mod");
fprintf(stderr," Nav17cw.mod");
fprintf(stderr," Nav17vw.mod");
fprintf(stderr," Nav17vw_named.mod");
fprintf(stderr," Nav18cw.mod");
fprintf(stderr," Nav18hw.mod");
fprintf(stderr," Nav18hw_named.mod");
fprintf(stderr," Nav18tf.mod");
fprintf(stderr," Nav19hw.mod");
fprintf(stderr," Nav19tf.mod");
fprintf(stderr," hcn_kn.mod");
fprintf(stderr," hcn_tf.mod");
fprintf(stderr," k_conc.mod");
fprintf(stderr," k_dummy.mod");
fprintf(stderr," na_accu.mod");
fprintf(stderr," na_conc.mod");
fprintf(stderr," na_dummy.mod");
fprintf(stderr," nadifl.mod");
fprintf(stderr," naiTest.mod");
fprintf(stderr, "\n");
    }
_IRamp_reg();
_KAtf_reg();
_KMtf_reg();
_Kacw_reg();
_Kdrcw_reg();
_Kdrtf_reg();
_Kleak_reg();
_Naleak_reg();
_Nav16bk_reg();
_Nav16hs_reg();
_Nav17cw_reg();
_Nav17vw_reg();
_Nav17vw_named_reg();
_Nav18cw_reg();
_Nav18hw_reg();
_Nav18hw_named_reg();
_Nav18tf_reg();
_Nav19hw_reg();
_Nav19tf_reg();
_hcn_kn_reg();
_hcn_tf_reg();
_k_conc_reg();
_k_dummy_reg();
_na_accu_reg();
_na_conc_reg();
_na_dummy_reg();
_nadifl_reg();
_naiTest_reg();
}
