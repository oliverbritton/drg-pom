/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__kmtf
#define _nrn_initial _nrn_initial__kmtf
#define nrn_cur _nrn_cur__kmtf
#define _nrn_current _nrn_current__kmtf
#define nrn_jacob _nrn_jacob__kmtf
#define nrn_state _nrn_state__kmtf
#define _net_receive _net_receive__kmtf 
#define _f_rates _f_rates__kmtf 
#define rates rates__kmtf 
#define states states__kmtf 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gbar _p[0]
#define gk _p[1]
#define ik _p[2]
#define ns _p[3]
#define nf _p[4]
#define Dns _p[5]
#define Dnf _p[6]
#define ek _p[7]
#define _g _p[8]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_nsTauCalc(void);
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_kmtf", _hoc_setdata,
 "nsTauCalc_kmtf", _hoc_nsTauCalc,
 "rates_kmtf", _hoc_rates,
 0, 0
};
#define nsTauCalc nsTauCalc_kmtf
 extern double nsTauCalc( double , double );
 /* declare global and static user variables */
#define nfTau nfTau_kmtf
 double nfTau = 0;
#define nsTau nsTau_kmtf
 double nsTau = 0;
#define ninf ninf_kmtf
 double ninf = 0;
#define usetable usetable_kmtf
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gbar_kmtf", 0, 1e+009,
 "usetable_kmtf", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "nsTau_kmtf", "ms",
 "nfTau_kmtf", "ms",
 "gbar_kmtf", "S/cm2",
 "gk_kmtf", "S/cm2",
 "ik_kmtf", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double nf0 = 0;
 static double ns0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ninf_kmtf", &ninf_kmtf,
 "nsTau_kmtf", &nsTau_kmtf,
 "nfTau_kmtf", &nfTau_kmtf,
 "usetable_kmtf", &usetable_kmtf,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"kmtf",
 "gbar_kmtf",
 0,
 "gk_kmtf",
 "ik_kmtf",
 0,
 "ns_kmtf",
 "nf_kmtf",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 9, _prop);
 	/*initialize range parameters*/
 	gbar = 0.0018;
 	_prop->param = _p;
 	_prop->param_size = 9;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _KMtf_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 9, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 kmtf F:/CLPC48/drg-pom/drgpom/models/KMtf.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _znexp ;
 static double _zq10 ;
 static double *_t_ninf;
 static double *_t_nsTau;
 static double *_t_nfTau;
static int _reset;
static char *modelname = "IKM from Tigerholm 2014";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_rates(double);
static int rates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_rates(double);
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Dns = ( ninf - ns ) / nsTau ;
   Dnf = ( ninf - nf ) / nfTau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Dns = Dns  / (1. - dt*( ( ( ( - 1.0 ) ) ) / nsTau )) ;
 Dnf = Dnf  / (1. - dt*( ( ( ( - 1.0 ) ) ) / nfTau )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
    ns = ns + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / nsTau)))*(- ( ( ( ninf ) ) / nsTau ) / ( ( ( ( - 1.0 ) ) ) / nsTau ) - ns) ;
    nf = nf + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / nfTau)))*(- ( ( ( ninf ) ) / nfTau ) / ( ( ( ( - 1.0 ) ) ) / nfTau ) - nf) ;
   }
  return 0;
}
 static double _mfac_rates, _tmin_rates;
 static void _check_rates();
 static void _check_rates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_rates)/200.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 201; _x += _dx, _i++) {
    _f_rates(_x);
    _t_ninf[_i] = ninf;
    _t_nsTau[_i] = nsTau;
    _t_nfTau[_i] = nfTau;
   }
   _sav_celsius = celsius;
  }
 }

 static int rates(double _lv){ _check_rates();
 _n_rates(_lv);
 return 0;
 }

 static void _n_rates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rates(_lv); return; 
}
 _xi = _mfac_rates * (_lv - _tmin_rates);
 if (isnan(_xi)) {
  ninf = _xi;
  nsTau = _xi;
  nfTau = _xi;
  return;
 }
 if (_xi <= 0.) {
 ninf = _t_ninf[0];
 nsTau = _t_nsTau[0];
 nfTau = _t_nfTau[0];
 return; }
 if (_xi >= 200.) {
 ninf = _t_ninf[200];
 nsTau = _t_nsTau[200];
 nfTau = _t_nfTau[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 ninf = _t_ninf[_i] + _theta*(_t_ninf[_i+1] - _t_ninf[_i]);
 nsTau = _t_nsTau[_i] + _theta*(_t_nsTau[_i+1] - _t_nsTau[_i]);
 nfTau = _t_nfTau[_i] + _theta*(_t_nfTau[_i+1] - _t_nfTau[_i]);
 }

 
static int  _f_rates (  double _lv ) {
   double _lnfTau_alpha , _lnfTau_beta ;
  _zq10 = 3.3 ;
   ninf = 1.0 / ( 1.0 + exp ( - ( _lv + 30.0 ) / 6.0 ) ) ;
   nsTau = nsTauCalc ( _threadargscomma_ _lv , _zq10 ) ;
   _lnfTau_alpha = 0.00395 * exp ( ( _lv + 30.0 ) / 40.0 ) ;
   _lnfTau_beta = 0.00395 * exp ( - ( _lv + 30.0 ) / 20.0 ) * _zq10 ;
   nfTau = 1.0 / ( _lnfTau_alpha + _lnfTau_beta ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
    _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double nsTauCalc (  double _lx , double _zq10 ) {
   double _lnsTauCalc;
 if ( _lx < - 60.0 ) {
     _lnsTauCalc = 219.0 * _zq10 ;
     }
   else {
     _lnsTauCalc = 13.0 * _lx + 1000.0 * _zq10 ;
     }
   
return _lnsTauCalc;
 }
 
static void _hoc_nsTauCalc(void) {
  double _r;
   _r =  nsTauCalc (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  nf = nf0;
  ns = ns0;
 {
   rates ( _threadargscomma_ v ) ;
   ns = ninf ;
   nf = ninf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ek = _ion_ek;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gk = gbar * ( ns / 4.0 + 3.0 * nf / 4.0 ) ;
   ik = gk * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ek = _ion_ek;
 { error =  states();
 if(error){fprintf(stderr,"at line 45 in file KMtf.mod:\n        SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(ns) - _p;  _dlist1[0] = &(Dns) - _p;
 _slist1[1] = &(nf) - _p;  _dlist1[1] = &(Dnf) - _p;
   _t_ninf = makevector(201*sizeof(double));
   _t_nsTau = makevector(201*sizeof(double));
   _t_nfTau = makevector(201*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "KMtf.mod";
static const char* nmodl_file_text = 
  "TITLE IKM from Tigerholm 2014\n"
  "\n"
  "COMMENT\n"
  "IKM from Tigerholm 2014\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "UNITS {\n"
  "\n"
  "		 (mA) = (milliamp)\n"
  "		 (mV) = (millivolt)\n"
  "		 (S) = (siemens)\n"
  "}\n"
  "\n"
  "NEURON {\n"
  "		 SUFFIX kmtf\n"
  "		 USEION k READ ek WRITE ik\n"
  "		 RANGE gbar, gk, ik	\n"
  "	     GLOBAL ninf, nsTau, nfTau\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "		 gbar = 0.0018 (S/cm2) <0,1e9>		 \n"
  "}\n"
  "\n"
  "STATE {\n"
  "		ns nf\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "		 v (mV)\n"
  "		 celsius (degC)\n"
  "		 ek (mV)\n"
  "		 \n"
  "		 gk (S/cm2)\n"
  "		 ik (mA/cm2)\n"
  "		 ninf\n"
  "		 nsTau (ms) nfTau (ms)\n"
  "}\n"
  "\n"
  "LOCAL nexp\n"
  "\n"
  "? currents\n"
  "BREAKPOINT {\n"
  "        SOLVE states METHOD cnexp\n"
  "        gk = gbar*(ns/4 + 3*nf/4)\n"
  "		ik = gk*(v - ek)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	rates(v)\n"
  "	ns = ninf\n"
  "	nf = ninf\n"
  "}\n"
  "\n"
  "? states\n"
  "DERIVATIVE states {\n"
  "		rates(v)\n"
  "		ns' = (ninf-ns)/nsTau\n"
  "		nf' = (ninf-nf)/nfTau\n"
  "}\n"
  "\n"
  "LOCAL q10\n"
  "\n"
  "? rates\n"
  "PROCEDURE rates(v(mV)) { : Computes rate and other constants at current v.\n"
  "						 : Call once from HOC to initialize inf at resting v.\n"
  "						 LOCAL nfTau_alpha, nfTau_beta\n"
  "						 TABLE ninf, nsTau, nfTau DEPEND celsius FROM -100 TO 100 WITH 200\n"
  "						 \n"
  "UNITSOFF\n"
  "		\n"
  "		\n"
  "		q10 = 3.3\n"
  "		ninf = 1/(1 + exp(-(v+30)/6))\n"
  "		nsTau = nsTauCalc(v,q10)\n"
  "		\n"
  "		nfTau_alpha = 0.00395*exp((v+30)/40)\n"
  "		nfTau_beta = 0.00395*exp(-(v+30)/20)*q10\n"
  "		nfTau = 1/(nfTau_alpha + nfTau_beta)\n"
  "		\n"
  "}\n"
  "UNITSON\n"
  "\n"
  "FUNCTION nsTauCalc(x,q10) {  : Equation for nsTau\n"
  "		if (x < -60) {\n"
  "				 nsTauCalc =  219*q10\n"
  "		}else{\n"
  "				nsTauCalc = 13*x + 1000*q10\n"
  "		}\n"
  "}\n"
  ;
#endif
