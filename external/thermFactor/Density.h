// Density.h : extension to the gFOStoffdaten class
// Autor: majo (majo02@avt.rwth-aachen.de), thbe (thbe02@avt.rwth-aachen.de)
// UML-Diagramm:
//   look at "Dokumentation" folder
// Quellen:
// [1] Pack, R., 2011; Development of a Computational Routine for Estimation
//     of Activity Coefficents in Electrolyte Systems, Studienarbeit, AVT.PT
// [2] Mathias, P.M., 2004; Correlation for the Density of Multicomponent 
//	   Aqueous Electrolytes, Ind. Eng. Chem. Res. 2004, 43, p.: 6247-6252

#include <cmath>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>

#include "props.h"
#include "stoffdaten.hpp"
#include "Matrixoperation.h"

#if !defined(DENSITY)
#define DENSITY

class Density{
public:
	bool DensityInitialized;

	Density();
	Density(bool Anzeige);
	void calc_molar_mass(tpropdata *pd, //calculates the molar mass
		double *x, 
		double *MolarMass);
	double calc_molar_volume(tpropdata *pd, //calculates the molare volume of the electrolyes
		double *V_inf, double *V_zero, double *V_e_inf, double *V_e_zero,
		double *x_e, double x_w,
		double **B_ij, double **C_ij);
	void calc_density(tpropdata *pd, //calculates density
		const double *DenMethodInputs,
		double *Density_measure);
	void calc_density_deriv(tpropdata *pd, //calculates first derivativ of density
		double *density_dx,
		const double *MethodInputs);
	double calc_V_ex(tpropdata *pd, //calculates excess volume 
		double *x_e, double *V_inf, 
		double x_w, double spec_vol);
	void calc_ms_diff_matrix_from_molefrac(
		tpropdata *pd,
		const double *MSMethodInputs,
		double *D_ij_out);
	void calc_ms_diff_matrix_from_moledens(
		tpropdata *pd,
		const double *MSMethodInputs,
		double *D_ij_out);
	double eval_binary_ms(tpropdata *pd, 
		double concentration[], 
		int species_one, int species_two);
	void calculate_binary_X_ij(
		tpropdata *pd,
		double temp,
		int species_one,
		int species_two,
		double *B_ij,
		double *C_ij);
	void calc_binary_ms_diffkoeff(
		tpropdata *pd,
		int species_one,
		int species_two,
		double c_binary,
		double *Diff_Koeff);
	void select_ion_volumes(tpropdata *pd, 
		int species, 
		double *V_0, 
		double *V_inf);

protected:
	Matrixoperation calc_matrix;

	void kill_binary_coefficents(tpropdata *pd, 
		double **B_ij, double **C_ij);
	void fill_binary_coefficents(tpropdata *pd,
		double **B_ij, double **C_ij,
		double temperature);
	void calc_V_e_zero_deriv(tpropdata *pd,
		double *V_e_zero_dx,
		double *x_e, double x_w,
		double **C_ij);
	void calc_V_e_inf_deriv(tpropdata *pd,
		double *V_e_inf_dx,
		double *x_e, double x_w,
		double **B_ij);
	void calc_V_e_deriv(tpropdata *pd,
		double *V_e_dx,
		double *x_e, double x_w, 
		double **B_ij, double **C_ij);
	void calc_V_ex_deriv(tpropdata *pd,
		double *V_ex_dx,
		double *x_e, double x_w, 
		double **B_ij, double **C_ij);
};

#endif
