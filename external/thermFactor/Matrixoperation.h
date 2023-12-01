// Matrixoperation.h : matrixcalculation
// Autor: majo (majo02@avt.rwth-aachen.de), thbe (thbe02@avt.rwth-aachen.de)
// UML-Diagramm:
// -
// Quellen:
// [1] Pack, R., 2011; Development of a Computational Routine for Estimation
//     of Activity Coefficents in Electrolyte Systems, Studienarbeit, AVT.PT
// [2] Mathias, P.M., 2004; Correlation for the Density of Multicomponent 
//	   Aqueous Electrolytes, Ind. Eng. Chem. Res. 2004, 43, p.: 6247-6252

#include <stdio.h>
#include <math.h>
#include <iostream>

#if !defined(MATRIXOPERATION)
#define MATRIXOPERATION

class Matrixoperation{
public:
	bool MatrixoperationInitialized;

	Matrixoperation();
	Matrixoperation(bool Ausgabe); 
	void copy_5_3_array(double source[5][3], double target[5][3]);
	void Eval_4_species(double **DK, double *V, double *z, double *c , double VE, double *EigenValues);
};

#endif
