// eNRTL.h : Headerdatei zur Gemischberechnung nach eNRTL(siehe [1])
// Autor: Robert Pack, majo (majo02@avt.rwth-aachen.de), thbe (thbe02@avt.rwth-aachen.de)
// UML-Diagramm:
//   look at "Dokumentation" folder
// Quellen:
// [1] Pack, R., 2011; Development of a Computational Routine for Estimation
//     of Activity Coefficents in Electrolyte Systems, Studienarbeit, AVT.PT


#define ELESU				4.80320425e-10	// electron charge in esu	[statculomb]
#define BOLTZESU			1.3806488e-16	// boltzmann constant in esu[] 
#define ENRTLTREF			298.15			// reference temperature	[K]
#define RALG				8.3144621		// universal gas constant	[J / mok / K]
#define	PI					3.141592654		// pi						[-]
#define AVOG				6.02214129e+23	// Avogadro number			[1 / mol]
#define DEFGMELCCMCAAQ      8.0             // default GMELCC for water - ion pait interaction
#define DEFGMELCCCAMAQ      -4.0            // default GMELCC for ion pair - water interaction
#define DEFGMELCCMCANONAQ   10              // default GMELCC for non-water - ion pair interaction
#define DEFGMELCCCAMNONAQ   -2.0            // default GMELCC for ion pair - non-water interaction
#define DEFGMELCDMCA        0.0             // default GMELCD
#define DEFGMELCDCAM        0.0             // default GMELCD
#define DEFGMELCN           0.2             // default non-randomness facrtor (GMELCN)
#define RKTTRMAX			0.99			// max reduced temperature for rackett model 

#define EPSMAXRECALC		1E-6			// accuracy at wich new calculations for enartl sums are undertaken

#define LSTRING				256
#define INPUTDELI			" \t"

#define ALL					-1				// reference to all components
#define NOCOMP				-2				// reference to no components
#pragma warning(disable : 4996)

#define TRUE 1
#define FALSE 0
#define ARRSIZE(arr) (sizeof(arr)/sizeof(arr[0]))
#define FILEOUT

#include <stdio.h>
#include <ctime>

#include "props.h"
#include "stoffdaten.hpp"


#if !defined(ENRTLSTRUCT)
#define ENRTLSTRUCT


//function for initialize1.cpp
void calc_binary_ms_diffkoeff(
	tpropdata *pd,
	int species_one,
	int species_two,
	double c_binary,
	double *Diff_Koeff
	);
void calculate_prop(
    const char *MethodName,
    const int *OutputLength,
    const int *InputLengths,
    const double *MethodInputs,
    double *MethodOutputs
    );
void calculate_prop_deriv(
    const char *MethodName,
    const int *OutputLength,
    const int *InputLengths,
    const double *MethodInputs,
    const char *DInputName,
    const int *DInputLength,
    double *DerivOutputs 
	);
void calculate_therm_factor(
    double *MethodInputs,
    double *Therm_factors);

extern "C" {	
void calc_therm_factor_C( 
    const int nSpc,
    const double Temp,
    const double Press,
    const double *Mole_frac,
    double *Therm_factors);
}	

extern "C" {
void calc_ms_diff_matrix_from_moledens_C(
    const int nSpc,
    const double Temp,
    const double Press,
    const double *Mole_dens,
    double *D_ij_out);
}

extern "C" {
void calc_ms_diff_matrix_from_molefrac_C(
    const int nSpc,
    const double Temp,
    const double Press,
    const double *Mole_frac,
    double *D_ij_out);
}
void copy_5_3_array(
	double source[5][3], 
	double target[5][3]);
double eval_binary_ms(
	tpropdata *pd, 
	double concentration[], 
	int species_one, 
	int species_two);
char *fgetz(
	char *s, int n, 
	FILE *stream);
int get_komp(
	tpropdata *pd, 
	char *string, 
	int *komp);
tpropdata *getglob_pd();
paramstruct getglob_pm();
extern "C" {
int init_enrtl(
	const char *filename,
        int &nSpc); 
}	
void init_param(
	tpropdata *pd, 
	FILE *input);
void init_paramstruct(
	tpropdata *pd, 
	paramstruct *pm);
void init_plib(
	tpropdata *pd);
void init_plib_data(
	tpropdata *pd);
void init_plib_segment(
	tpropdata *pd);
void init_reac(
	tpropdata *pd); 
void kill_paramstruct(
	paramstruct *pm);
void kill_plib(
	tpropdata *pd);
void killstruct();
static char* mystrtok(
	char *input, 
	char *deli);
static int readline(
	FILE *ein, 
	char *zeile, 
	char *befehl, 
	char **feld);


//function for eNRTL1.cpp
void set_param_default(
	tpropdata *pd);
int update_pm(
	tpropdata *pd, 
	paramstruct *pm, 
	double t, 
	double *x);
void get_index_input(
	tplibdata *ld, 
	int *species, 
	int *index, 
	int size);
void get_lngm_x(
	tpropdata *pd, 
	paramstruct *pm, 
	double *lngm, 
	double* lngm_dt, 
	double* lngm_dx);
int get_index_calc(
	tplibdata *ld, 
	int *species, 
	int *index, 
	int choose);
void eval_sums(
	tpropdata *pd, 
	paramstruct *pm);
void eval_sums_dt(
	tpropdata *pd, 
	paramstruct *pm);
void eval_sums_dx(
	tpropdata *pd, 
	paramstruct *pm);
void eval_sums_mx(
	tpropdata *pd, 
	paramstruct *pm);
void eval_sums_ca(
	tpropdata *pd, 
	paramstruct *pm);
void eval_lngm_sr_x(
	tpropdata *pd, 
	paramstruct *pm, 
	double *lngm_sr, 
	double *lngm_sr_dt, 
	double *lngm_sr_dx);
void eval_lngm_sr_x_ref(
	tpropdata *pd, 
	paramstruct *pm, 
	double *lngm_ref, 
	double *lngm_ref_dt, 
	double *lngm_ref_dx);
void eval_lngm_lr_x(
	tpropdata *pd, 
	paramstruct *pm, 
	double *lngm_lr, 
	double *lngm_lr_dt, 
	double *lngm_lr_dx);
void eval_lngm_born_x(
	tpropdata *pd, 
	double t, 
	double *x, 
	double *lngm_born, 
	double *lngm_born_dt, 
	double *lngm_born_dx);
void get_param(
	tpropdata* pd, 
	paramstruct* pm); 
void get_param_dt(
	tpropdata* pd, 
	paramstruct* pm);
void get_param_dx(
	tpropdata* pd, 
	paramstruct* pm);
double sum_prototype(
	tpropdata *pd, 
	paramstruct *pm, 
	int *lngm_spec, 
	char sum_specs);
int testmolec(
	tplibdata *ld, 
	int komp);
void get_aphi_x(
	tpropdata *pd, 
	double t, 
	double *x, 
	double *aphi, 
	double *aphi_dt, 
	double *aphi_dx);
void get_mxdiec_x(
	tpropdata* pd, 
	double t, 
	double* x, 
	double* mxdiec, 
	double* mxdiec_dt, 
	double* mxdiec_dx);
void get_mxdens_x(
	tpropdata* pd, 
	double t, 
	double* x, 
	double* mxdens, 
	double* mxdens_dt, 
	double* mxdens_dx);


#endif
