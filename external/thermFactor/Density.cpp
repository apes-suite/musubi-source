// Density.cpp : extension to the gFOStoffdaten class
// Autor: majo (majo02@avt.rwth-aachen.de), thbe (thbe02@avt.rwth-aachen.de)
// UML-Diagramm:
//   look at "Dokumentation" folder
// Quellen:
// [1] Pack, R., 2011; Development of a Computational Routine for Estimation
//     of Activity Coefficents in Electrolyte Systems, Studienarbeit, AVT.PT
// [2] Mathias, P.M., 2004; Correlation for the Density of Multicomponent 
//	   Aqueous Electrolytes, Ind. Eng. Chem. Res. 2004, 43, p.: 6247-6252

#include "Density.h"

using namespace std;


Density::Density(){};

Density::Density(bool Ausgabe){
	DensityInitialized = true;
	
	if(Ausgabe){
		if(DensityInitialized){
			cout<<">> Initialisierung von Density \t\t [erfolgreich]"<<endl<<endl;
		}else{
			cout<<">> Initialisierung von Density \t\t [FEHLER]"<<endl<<endl;
		}
	}
}

void Density::calc_ms_diff_matrix_from_molefrac(tpropdata *pd,	
	const double *MSMethodInputs,
	double *D_ij_out)
{
	int i, j, index;
	
	// introduce an Array of diffusion coefficients 
	double** D_ij = new double* [pd->nstoff];
	for(i=0; i < pd-> nstoff; i++)
	{
		D_ij[i] = new double[pd->nstoff];
	}


	// introduce mole fractions from the inputs
	double* x		= new double[pd-> nstoff];
	for(i=0; i<pd-> nstoff; ++i){
		x[i] = MSMethodInputs[i+2];
	} 


	// call density calculation
	double *density = new double[3]; //specific density [kmol/m3] is density[0]
	calc_density(pd, MSMethodInputs, density);

	// calculate molar concentrations
	double* c		= new double[pd-> nstoff] ; // in [mol/m3] 
	for(i=0; i<pd-> nstoff; ++i){
		c[i] = x[i] * density[0] * 1000; // !!! Visser model expects [mol/m3] as input
	} 

	double checker;
	// fill the MS diffusion matrix
	for(i=0; i < pd-> nstoff; i++) {
		for(j=0; j < pd-> nstoff; j++) {
			if(i == j) 
			{
				D_ij[i][j] = 1; 
			}
			else
			{
				D_ij[i][j] = eval_binary_ms(pd ,c ,i ,j);
				checker = D_ij[i][j];
			}
		}
	}

	// convert two dimensional array to a vector
	for(i=0; i < pd-> nstoff; i++) {
		for(j=0; j < pd-> nstoff; j++) {
			index = j*pd->nstoff + i;
			D_ij_out[index] = D_ij[i][j];
		}
	}

	// Speicher wieder freigeben
	for(i=0; i < pd-> nstoff; i++)
	{
		delete[]  D_ij[i];
	}
	delete[] D_ij;
	delete[] x;
	delete[] c;
	delete[] density;
}

void Density::calc_ms_diff_matrix_from_moledens(tpropdata *pd,	
	const double *MSMethodInputs,
	double *D_ij_out)
{
	int i, j, index;
	
	// introduce an Array of diffusion coefficients 
	double** D_ij = new double* [pd->nstoff];
	for(i=0; i < pd-> nstoff; i++)
	{
		D_ij[i] = new double[pd->nstoff];
	}


	// introduce molar density from the inputs
	double* c		= new double[pd-> nstoff];
	for(i=0; i<pd-> nstoff; ++i){
		c[i] = MSMethodInputs[i+2];
	} 

	double checker;
	// fill the MS diffusion matrix
	for(i=0; i < pd-> nstoff; i++) {
		for(j=0; j < pd-> nstoff; j++) {
			if(i == j) 
			{
				D_ij[i][j] = 1; 
			}
			else
			{
				D_ij[i][j] = eval_binary_ms(pd ,c ,i ,j);
				checker = D_ij[i][j];
			}
		}
	}

	// convert two dimensional array to a vector
	for(i=0; i < pd-> nstoff; i++) {
		for(j=0; j < pd-> nstoff; j++) {
			index = j*pd->nstoff + i;
			D_ij_out[index] = D_ij[i][j];
		}
	}

	// Speicher wieder freigeben
	for(i=0; i < pd-> nstoff; i++)
	{
		delete[]  D_ij[i];
	}
	delete[] D_ij;
	delete[] c;
}

double Density::eval_binary_ms(tpropdata *pd, double concentration[], int species_one, int species_two)
{
	int i;
	double* D_ij_binary = new double [3];
	int * is_anion  = new int[pd->nstoff];
	int* is_cation  = new int [pd->nstoff];
	int* is_neutral = new int [pd->nstoff]; 

	int ani_counter = 0, cati_counter = 0, neutral_counter = 0 ;

	// get local libdata struct
	tplibdata* ld	= (tplibdata*) pd->libdata;

	// initialize the labels 
	for(i=0; i< pd->nstoff; i++)
	{
		is_anion[i] = 0;
		is_cation[i] = 0;
		is_neutral[i] = 0;
	}

	// first the types of the species are evaluated
	for(i=0; i< pd->nstoff; i++)
	{
		if(ld->charge[i] < 0) 
		{
			is_anion[i] = 1;
			ani_counter += 1;
		}
		else if (ld->charge[i] > 0) 
		{
			is_cation[i] = 1;
			cati_counter += 1; 
		}
		else
		{
			is_neutral[i] = 1;
			neutral_counter += 1; 
		}
	}
	// in the case of a binary anion cation pair we can directly calculate binary DK
	// with the higher concentration
	if(is_anion[species_one] == 1 && is_cation[species_two] == 1  || is_anion[species_two] == 1 && is_cation[species_one])
	{	// binary double cation anion diffusion coefficient
		double D_ca=0;
		if(concentration[species_one] > concentration[species_two])
		{
			calc_binary_ms_diffkoeff(pd,species_one,species_two,concentration[species_one], D_ij_binary);
		}
		else
		{
			calc_binary_ms_diffkoeff(pd,species_one,species_two,concentration[species_two], D_ij_binary);
		}
		D_ca = D_ij_binary[0];
        //Speicher wieder freigeben
        	delete [] D_ij_binary;
        	delete [] is_anion;
        	delete [] is_cation;
        	delete [] is_neutral;
		return D_ca;
	}

	// in case we have two cations
	else if (is_cation[species_one] == 1 && is_cation[species_two] == 1  )
	{	// binary double cation cation diffusion coefficient
		double D_cc=1;
		// diffusion coefficients of cation-(common) anion interactions
		double* D_K1_A = new double [ani_counter];
		double* D_K2_A = new double [ani_counter];
		// counter to fill D_K1_A and D_K2_A
		int j = 0;
		// get all diffusion coefficients of cation-(common) anion interactions

		for (i=0; i < pd->nstoff; i++)
		{
			if(is_anion[i] == 1)
			{	
				// get binary diffusion coefficient of first cation and ith anion
				calc_binary_ms_diffkoeff(pd,species_one,i,concentration[species_one], D_ij_binary);
				D_K1_A[j] = D_ij_binary[0];
				// get binary diffusion coefficient of second cation and ith anion
				calc_binary_ms_diffkoeff(pd,species_one,i,concentration[species_one], D_ij_binary);
				D_K2_A[j] = D_ij_binary[0];
				j++;

			}
		}
		// calculate the cation-cation diffusivity as geometric mean
		for(i=0; i < ani_counter; i++)
		{
			if(ani_counter == 0) 
			{
				printf("\n division by zero! check ani_counter");
				throw 9999;
			}
			D_cc *= pow(D_K1_A[i] * D_K2_A[i], 1/ani_counter);
		}

		delete [] D_K1_A;
		delete [] D_K2_A;
        //Speicher wieder freigeben
        	delete [] D_ij_binary;
        	delete [] is_anion;
        	delete [] is_cation;
        	delete [] is_neutral;
		return -D_cc;
	}

	// in case we have two anions (has no theoretical basis according to Visser!!!) 
	else if (is_anion[species_one] == 1 && is_anion[species_two] == 1  )
	{	// binary double cation cation diffusion coefficient
		double D_aa=1;
		// diffusion coefficients of (common) cation- anion interactions
		double* D_A1_K = new double [cati_counter];
		double* D_A2_K = new double [cati_counter];
		// counter to fill D_A1_K and D_A2_K
		int j = 0;
		// get all diffusion coefficients of (common) cation- anion interactions

		for (i=0; i < pd->nstoff; i++)
		{
			if(is_cation[i] == 1)
			{	
				// get binary diffusion coefficient of first cation and ith anion
				calc_binary_ms_diffkoeff(pd,species_one,i,concentration[species_one], D_ij_binary);
				D_A1_K[j] = D_ij_binary[0];
				// get binary diffusion coefficient of second cation and ith anion
				calc_binary_ms_diffkoeff(pd,species_one,i,concentration[species_one], D_ij_binary);
				D_A2_K[j] = D_ij_binary[0];
				j++;

			}
		}
		// calculate the cation-cation diffusivity as geometric mean
		for(i=0; i < cati_counter; i++)
		{
			if(cati_counter == 0) 
			{
				printf("\n division by zero! check cati_counter");
				throw 9999;
			}
			D_aa *= pow(D_A1_K[i] * D_A2_K[i], 1/cati_counter);
		}

		delete [] D_A1_K;
		delete [] D_A2_K;
        //Speicher wieder freigeben
        	delete [] D_ij_binary;
        	delete [] is_anion;
        	delete [] is_cation;
        	delete [] is_neutral;
		return -D_aa;
	}
	// in case we have a cation - solvent interaction 
	else if (is_cation[species_one] == 1 && is_neutral[species_two] == 1 
		|| is_neutral[species_one] ==1 && is_cation[species_two] ==1 )
	{
		// binary cation water diffusion coefficient
		double D_cw=0;
		double ionic_strenght = 0;
		// diffusion coefficients of (common) cation- anion interactions
		double* D_K_w_A = new double [ani_counter];
		// counter to fill D_A1_K and D_A2_K
		int j = 0;
		// calculate ionic strength of the solution
		for(i=0; i < pd->nstoff; i++)
		{
			ionic_strenght += (ld->charge[i]) * (ld->charge[i]) * concentration[i];
		};

		// get all cation - water diffusion coefficients determined for the different anions
		for (i=0; i < pd->nstoff; i++)
		{
			// if the cation is species one
			if(is_anion[i] == 1 && is_cation[species_one] == 1)
			{	
				calc_binary_ms_diffkoeff(pd,species_one,i,ionic_strenght, D_ij_binary);
				D_K_w_A[j] = D_ij_binary[1]; // D_ij_binary[1] is the cation-water diffusivity
				j++;

			}
			// if the cation is species two
			else if(is_anion[i] == 1 && is_cation[species_two] == 1)
			{	
				calc_binary_ms_diffkoeff(pd,species_two,i,ionic_strenght, D_ij_binary);
				D_K_w_A[j] = D_ij_binary[1]; // D_ij_binary[1] is the cation-water diffusivity
				j++;

			}
		}
		// we need the total concentration of anions
		double sum_c_a = 0;
		for (i=0; i < pd->nstoff; i++)
		{
			if(is_anion[i] == 1 )
			{
				sum_c_a += concentration[i];
			}
		}


		// calculate the cation-water diffusivity by empirical correlation
		j =0; // reset the counter for diffusion coefficient to zero
		for (i=0; i < pd->nstoff; i++)
		{
			if(is_anion[i] == 1 )
			{
				if(sum_c_a == 0) 
				{
					printf("\n division by zero! check sum_c_a");
					throw 9999;
				}
				D_cw += concentration[i] / sum_c_a * D_K_w_A[j];
				j++;
			}
		}
		
		delete [] D_K_w_A;
                //Speicher wieder freigeben
        	delete [] D_ij_binary;
        	delete [] is_anion;
        	delete [] is_cation;
        	delete [] is_neutral;

		return D_cw;
	}
	// in case we have a anion - solvent pair
	else if (is_anion[species_one] == 1 && is_neutral[species_two] == 1 
		|| is_neutral[species_one] ==1 && is_anion[species_two] ==1 )
	{
		// binary cation water diffusion coefficient
		double D_aw=0;
		double ionic_strenght = 0;
		// diffusion coefficients of between anion and water for different 
		// in the presence of different cations
		double* D_A_w_K = new double [cati_counter];
		// counter to fill D_A1_K and D_A2_K
		int j = 0;
		// calculate ionic strength of the solution
		for(i=0; i < pd->nstoff; i++)
		{
			ionic_strenght += (ld->charge[i]) * (ld->charge[i]) * concentration[i];
		};

		// calculate all anion - water diffusion coefficients determined for the different cations
		for (i=0; i < pd->nstoff; i++)
		{
			// if the cation is species one
			if(is_cation[i] == 1 && is_anion[species_one] == 1)
			{	
				calc_binary_ms_diffkoeff(pd,species_one,i,ionic_strenght, D_ij_binary);
				D_A_w_K[j] = D_ij_binary[2]; // D_ij_binary[2] is the anion-water diffusivity
				j++;

			}
			// if the cation is species two
			else if(is_cation[i] == 1 && is_anion[species_two] == 1)
			{	
				calc_binary_ms_diffkoeff(pd,species_two,i,ionic_strenght, D_ij_binary);
				D_A_w_K[j] = D_ij_binary[2]; // D_ij_binary[2] is the anion-water diffusivity
				j++;

			}
		}
		// we need the total concentration of cations
		double sum_c_c = 0;
		for (i=0; i < pd->nstoff; i++)
		{
			if(is_cation[i] == 1 )
			{
				sum_c_c += concentration[i];
			}
		}

		// calculate the cation-water diffusivity by empirical correlation
		j =0; // reset the counter for diffusion coefficient to zero
		for (i=0; i < pd->nstoff; i++)
		{
			if(is_cation[i] == 1 )
			{
				if(sum_c_c == 0) 
				{
					printf("\n division by zero! check sum_c_c");
					throw 9999;
				}
				D_aw += concentration[i] / sum_c_c * D_A_w_K[j];
				j++;
			}
		}

		delete [] D_A_w_K;
        //Speicher wieder freigeben
        	delete [] D_ij_binary;
        	delete [] is_anion;
        	delete [] is_cation;
        	delete [] is_neutral;

		return D_aw;
	}
	else 
	{
		printf(" \n the type of binary pair is not supported by the current implementation! ");
		printf(" \n e.g. only one neutral species is allowed ");
		throw 9999;
	}
//Speicher wieder freigeben
	delete [] D_ij_binary;
	delete [] is_anion;
	delete [] is_cation;
	delete [] is_neutral;

}

void Density::calc_density(tpropdata *pd,	
	const double *DenMethodInputs,
	double *Density_measure)
{   // function to calculate the density of an electrolyte solution according to the model published by Mathias 2004
	// returns different important density measures in the Density_measure array 
	// note that the molar quantities that are returned are given in [kmol/m3] and [m3/kmol] respectively

	tplibdata* ld	= (tplibdata*) pd->libdata; 
	
	// check if water is the first component in the system
	if( strcmpu(pd->stname[0],"water") || !strcmpu(pd->stname[0],"H2O")  
		|| !strcmpu(pd->stname[0],"wasser") || !strcmpu(pd->stname[0],"H20") ) 
	{
		printf("\n first species in the system must be water (note: only one solvent possible in density model)");
		throw 9999;
	}

	// introduce n-1 dimensional (!!no water!!) arrays of binary coefficients 
	double** B_ij = new double* [pd->nstoff - 1];
	for(int i=0; i < pd-> nstoff - 1 ; i++)
	{
		B_ij[i] = new double[pd->nstoff - 1];
	}
	double** C_ij = new double* [pd->nstoff - 1];
	for(int i=0; i < pd-> nstoff - 1; i++)
	{
		C_ij[i] = new double[pd->nstoff - 1];
	}
	
	double temp = DenMethodInputs[0]; // local variable for the temperature
	double* V_inf = new double[pd->nstoff - 1];	 // n-1 dimensional (no water) array for the pure and infinite dilution densities of the ions
	double* V_zero = new double[pd->nstoff - 1];
	double* MolarMass = new double;	//molar mass of the solution
	double V_e = 0;	// total specific volume of the electrolyte
	double c_t = 0; // total molar density
	double spec_vol = 0; // specific volume
	double rho = 0;	// total density of the solution
	double V_ideal = 0; // ideal contribution of the molar volume 
	double V_ex = 0; // excess contribution of the molar volume
	double V_e_inf = 0; 
	double V_e_zero = 0;
	double V_w  = 0.01801528*1000 /1000; // (kg/kmol) / (kg/m3)

	// introduce mole fractions of ions (!!no water!!) from the inputs
	// note that the species list needs to be shifted to exclude water
	double* x_e	= new double[pd-> nstoff - 1];
	for(int i=0; i<pd-> nstoff-1; ++i){
		x_e[i] = DenMethodInputs[i+3];
	} 
	
	double x_w = DenMethodInputs[2]; //mole fraction of water

	double* x	= new double[pd-> nstoff]; //complete vector of mole fractions
	for(int i=0; i<pd-> nstoff; ++i){
		x[i] = DenMethodInputs[i+2];
	}

	fill_binary_coefficents(pd, B_ij, C_ij, temp);

	// calculates molar volume of the electrolyte soluion
	V_e = calc_molar_volume(pd, V_inf, V_zero, &V_e_inf, &V_e_zero, x_e, x_w, B_ij, C_ij);

	// calculates the total concentration of the solution (!!impelemtation is only for 25°C!! 
	// has to be extended to a temperature dependent model for the partial molar volume of water
	// add contribution of water to the specific volume
	spec_vol = x_w*V_w;
	// add the contributions of the different ions
	for(int i=0; i<pd->nstoff-1; i++){
		// add ideal contribution
		spec_vol += x_e[i]*V_inf[i];
		// non ideal contributions
		V_ex += (V_e-V_inf[i])*x_e[i];
	}
	// add nonideal contributions to the specifiv volume  [m3/mol]
	spec_vol += V_ex;

	//if to avoid division by zero
	if(spec_vol == 0) 
		{
			printf("\n division by zero! check spec_vol in the density calculation");
			throw 9999;
		}
	// calculate specific molar density
	c_t = 1 / spec_vol;
	
	// calculate total density
	calc_molar_mass(pd, x, MolarMass);
	rho = c_t * (*MolarMass);

	// different density measures as output in one array
	Density_measure[0] = c_t; // specific molar density [kmol/m3] = [mol/l]
	Density_measure[1] = rho; // density [kg/m3]
	Density_measure[2] = calc_V_ex(pd, x_e, V_inf, x_w, spec_vol); // exzess volume [kmol/m3]

	// delete variables
	delete [] V_inf; 
	delete [] V_zero; ; 
	delete [] MolarMass ;

	delete [] x;
	delete [] x_e;

	kill_binary_coefficents(pd, B_ij, C_ij);
}

void Density::fill_binary_coefficents(tpropdata *pd,
	double **B_ij, double **C_ij,
	double temperature
	)
{
	double *B_binary = new double;
	double *C_binary = new double;
		// fill the B_ij and C_ij matrix
	for(int i=0; i < pd-> nstoff - 1; i++) {
		for(int j=0; j < pd-> nstoff - 1; j++) {
			if(i == j) 
			{
				B_ij[i][j] = 0; 
				C_ij[i][j] = 0; 
			}
			else
			{
				// calculate and fill the binary parameters
				// note that the species list needs to be shifted to exclude water
				calculate_binary_X_ij(pd, temperature, i + 1, j + 1, B_binary, C_binary);
				B_ij[i][j] = *B_binary;
				C_ij[i][j] = *C_binary;
			}
		}
	}

	delete[] B_binary;
	delete[] C_binary;
}

void Density::kill_binary_coefficents(tpropdata *pd, double **B_ij, double **C_ij)
{
	for(int i=0; i < pd-> nstoff - 1 ; i++)
	{
		delete[]  B_ij[i];
		delete[]  C_ij[i];
	}
	delete[] B_ij;
	delete[] C_ij;
}

void Density::calculate_binary_X_ij(
	tpropdata *pd,
	double temp,
	int species_one,
	int species_two,
	double *B_ij,
	double *C_ij
	)
{
	// function to select the binary parameters for the density model
	// binary coefficients are given in [m3/kmol]
	int  i;

	double B0_ij;
	double B1_ij;
	double C0_ij;
	double C1_ij;

	if( !strcmpu(pd->stname[species_one],"H") && !strcmpu(pd->stname[species_two],"Cl")  
		|| !strcmpu(pd->stname[species_one],"Cl") && !strcmpu(pd->stname[species_two],"H") ) 
	{
		B0_ij = -0.00355;
		B1_ij = 0.00549;
		C0_ij = -0.00509;
		C1_ij = 0;


	}
	else 	if( !strcmpu(pd->stname[species_one],"Na") && !strcmpu(pd->stname[species_two],"OH")  
		|| !strcmpu(pd->stname[species_one],"OH") && !strcmpu(pd->stname[species_two],"Na") ) 
	{
		B0_ij = -0.00763;
		B1_ij = 0.00981;
		C0_ij = 0.01006;
		C1_ij = 0;
	}
	else 	if( !strcmpu(pd->stname[species_one],"Na") && !strcmpu(pd->stname[species_two],"SO4")  
		|| !strcmpu(pd->stname[species_one],"SO4") && !strcmpu(pd->stname[species_two],"Na") ) 
	{
		B0_ij = 0;
		B1_ij = 0;
		C0_ij = 0;
		C1_ij = 0;
	}
	else 	if( !strcmpu(pd->stname[species_one],"H") && !strcmpu(pd->stname[species_two],"SO4")  
		|| !strcmpu(pd->stname[species_one],"SO4") && !strcmpu(pd->stname[species_two],"H") ) 
	{
		B0_ij = 0;
		B1_ij = 0;
		C0_ij = 0;
		C1_ij = 0;
	}
	else 	if( !strcmpu(pd->stname[species_one],"Na") && !strcmpu(pd->stname[species_two],"Cl")  
		|| !strcmpu(pd->stname[species_one],"Cl") && !strcmpu(pd->stname[species_two],"Na") ) 
	{
		B0_ij = 0.00155;
		B1_ij = 0.00776;
		C0_ij = 0.01160;
		C1_ij = 0;
	}
	else
	{
		B0_ij = 0;
		B1_ij = 0;
		C0_ij = 0;
		C1_ij = 0;
	//	printf( "\n No binary paramaters in the density model for the pair: %s %s",  pd->stname[species_one] ,pd->stname[species_two] );
	//	printf("\n -->dummy values are used");
	}

	*B_ij = B0_ij + (temp - 298.15)/ 298.15 *  B1_ij;
	*C_ij = C0_ij + (temp - 298.15)/ 298.15 *  C1_ij;
}

void Density::select_ion_volumes(tpropdata *pd, int species, double *V_0, double *V_inf)
	// function that selects volume of an ionic species at infinite dilution and pure state
	// values are given in [m3/kmol]
{
	
	if( !strcmpu(pd->stname[species],"H"))
	{
		*V_0 = 0.01682;
		*V_inf = -0.00318;
	}
	else if( !strcmpu(pd->stname[species],"Na"))
	{
		*V_0 = 0.01148;
		*V_inf = -0.00852;
	}
	else if( !strcmpu(pd->stname[species],"Cl"))
	{
		*V_0 = 0.01469;
		*V_inf = 0.02469;
	}
	else if( !strcmpu(pd->stname[species],"SO4"))
	{
		*V_0 = 0.03235;
		*V_inf = 0.03235;
	}
	else if( !strcmpu(pd->stname[species],"WATER"))
	{
		*V_0 = 0.01806862037;
		*V_inf = 0.01806862037;
	}
	else
	{
		*V_0 = 0.001;
		*V_inf = 0.03;
		printf( "\n No pure densities for ion: %s %s",  pd->stname[species] );
		printf("\n -->dummy values are used");
	}
}

void Density::calc_molar_mass(tpropdata *pd, double *x, double *MolarMass)
{ // function that calculates the molar mass of the solution
	double M_i;
	double M=0;

	for(int i=0; i <= pd-> nstoff - 1 ; i++)
	{
		M_i = 0;
		if( !strcmpu(pd->stname[i],"H")) {M_i = 1.00794;}
		 else if( !strcmpu(pd->stname[i],"H2O") || !strcmpu(pd->stname[i],"water")) {M_i = 18.01528;}
	 	  else if( !strcmpu(pd->stname[i],"Na")) {M_i = 22.989770;}
		   else if( !strcmpu(pd->stname[i],"Cl")) {M_i = 35.4527;}
		    else if( !strcmpu(pd->stname[i],"SO4")) {M_i = 96.061;}
		     else {
				M_i = 0;
				printf( "\n No molar mass for species: %s %s",  pd->stname[i] );
				printf("\n -->calculated molar mass is wrong !!!");
			}
	
		M += (M_i * x[i]);
	}
	*MolarMass = M;
}

double Density::calc_V_ex(tpropdata *pd,
	double *x_e, double *V_inf, 
	double x_w, double spec_vol)
	// method to calculate the molar exzess volume 
{

	// calculate ideal contribution of the specific molar volume
	double V_e_inf_first=0;
	for (int i=0; i < pd->nstoff - 1; i++){
		V_e_inf_first  += x_e[i] * V_inf[i];
	}
	double V_ideal  = x_w * 0.01806862037 + V_e_inf_first;
	
	// calculate the excess contribution to the specific molar volume
	return spec_vol - V_ideal;
}

double Density::calc_molar_volume(tpropdata *pd,
	double *V_inf, double *V_zero, double *V_e_inf, double *V_e_zero,
	double *x_e, double x_w,
	double **B_ij, double **C_ij 
	)
// calculates molar volume proposed by [2] and returns the total specific volume of the electrolyte
//   V_inf: infinite dilution specific volume
//   V_zero: pure specific volume
//   V_e_inf: electrolyte infinite dilution specific volume
//   V_e_zero: pure  specific electrolyte volume
// input parameters:
//   x_e: molefracs of electrolyts
//   x_w: molefrac of water
//   B_ij,C_ij: symmetric binary coefficients
{
	// introduces the infinite dilution specific volume of the solution
	double V_e_inf_first = 0;
	double V_e_inf_second = 0;
	bool calc_V_e_inf = false;

	// introduces the pure specific volume of the solution
	double V_e_zero_first = 0;
	double V_e_zero_second = 0;
	bool calc_V_e_zero = false;

	// scalar quantities which are used to fill the arrays
	double* V_zero_i = new double; 
	double* V_inf_i = new double; 

	// fill the array for the pure and infinite dilution densities of the ions
	if((V_inf!=NULL)||(V_zero!=NULL)){
		for(int i=0; i < pd-> nstoff - 1; i++) {
			select_ion_volumes(pd, i+1, V_zero_i, V_inf_i);
			V_inf[i] = *V_inf_i;
			V_zero[i] = *V_zero_i;
		}
	}

	// calculates infinite dilution and pure specific volume of the electrolyte
	if(V_e_inf!=NULL){
		// calculates the first term infinite dilution specific volume of the electrolyte
		for (int i=0; i < pd->nstoff - 1; i++){
			V_e_inf_first  += x_e[i] * V_inf[i];
		}
		// calculates the second term infinite dilution specific volume of the electrolyte
		for(int i=0; i < pd->nstoff - 1; i++) {
			for(int j=0; j < pd->nstoff - 1; j++) {
				V_e_inf_second += x_e[i] * x_e[j] * B_ij[i][j];
			}
		}
		*V_e_inf = V_e_inf_first + V_e_inf_second;
		calc_V_e_inf = true;
	}

	// calculates pure specific volume of the electrolyte
	if(V_e_zero!=NULL){
		// calculates the first term of pure specific volume of the electrolyte
		for (int i=0; i < pd->nstoff - 1; i++){
			V_e_zero_first += x_e[i] * V_zero[i];
		}
		// calculates the second term of pure specific volume of the electrolyte
		for(int i=0; i < pd->nstoff - 1; i++) {
			for(int j=0; j < pd->nstoff - 1; j++) {
				V_e_zero_second += x_e[i] * x_e[j] * C_ij[i][j];   
			}
		}
		// calculates the total molar volume of the electrolyte solution 
		*V_e_zero = V_e_zero_first + V_e_zero_second;
		calc_V_e_zero = true;
	}

	// clears memory
	delete[] V_zero_i;
	delete[] V_inf_i;

	// calculates the total specific volume of the electrolyte
	return (calc_V_e_inf&&calc_V_e_zero)?((*V_e_inf)+((*V_e_zero)-(*V_e_inf))*(1-x_w)):0;
}

void Density::calc_V_e_zero_deriv(tpropdata *pd,
	double *V_e_zero_dx,
	double *x_e, double x_w,
	double **C_jk
	)
{
	double V_tmp_C_ik;
	double V_tmp_C_ji;

	double *V_zero = new double[pd->nstoff - 1];
	double *V_inf = new double[pd->nstoff - 1];

	//zurücksetzen der Ausgabewerte
	for(int i=0; i<pd->nstoff-1; i++){
		V_e_zero_dx[i] = 0;
	}
	
	calc_molar_volume(pd, V_inf, V_zero, NULL, NULL, x_e, x_w, NULL, NULL);

	for(int count=0; count<pd->nstoff-1; count++){
		V_tmp_C_ik = 0;
		V_tmp_C_ji = 0;
		for(int k=0; k<pd->nstoff-1; k++){
			V_tmp_C_ik += x_e[k]*C_jk[count][k];
		}
		for(int j=0; j<pd->nstoff-1; j++){
			V_tmp_C_ji += x_e[j]*C_jk[j][count];
		}

		V_e_zero_dx[count] = V_zero[count] + V_tmp_C_ik + V_tmp_C_ji;
	}

	delete[] V_zero;
	delete[] V_inf;
}

void Density::calc_V_e_inf_deriv(tpropdata *pd,
	double *V_e_inf_dx, //MethodOutputs
	double *x_e, double x_w,
	double **B_jk
	)
{
	double V_tmp_B_ik;
	double V_tmp_B_ji;

	double *V_zero = new double[pd->nstoff - 1];
	double *V_inf = new double[pd->nstoff - 1];
	
	//zurücksetzen der Ausgabewerte
	for(int i=0; i<pd->nstoff-1; i++){
		V_e_inf_dx[i] = 0;
	}

	calc_molar_volume(pd, V_inf, V_zero, NULL, NULL, x_e, x_w, NULL, NULL);

	for(int count=0; count<pd->nstoff-1; count++){
		V_tmp_B_ik = 0;
		V_tmp_B_ji = 0;
		for(int k=0; k<pd->nstoff-1; k++){
			V_tmp_B_ik += x_e[k]*B_jk[count][k];
		}
		for(int j=0; j<pd->nstoff-1; j++){
			V_tmp_B_ji += x_e[j]*B_jk[j][count];
		}

		V_e_inf_dx[count] = V_inf[count] + V_tmp_B_ik + V_tmp_B_ji;
	}

	delete[] V_zero;
	delete[] V_inf;
}

void Density::calc_V_e_deriv(tpropdata *pd,
	double *V_e_dx, //MethodOutput
	double *x_e, double x_w, 
	double **B_ij, double **C_ij
	)
{
	double *V_zero = new double[pd->nstoff - 1];
	double *V_inf = new double[pd->nstoff - 1];

	double *V_e_inf_dx = new double[pd->nstoff - 1];
	calc_V_e_inf_deriv(pd, V_e_inf_dx, x_e, x_w, B_ij);

	double *V_e_zero_dx = new double[pd->nstoff - 1];
	calc_V_e_zero_deriv(pd, V_e_zero_dx, x_e, x_w, C_ij);

	double V_e_zero = 0;
	double V_e_inf = 0;
	calc_molar_volume(pd, V_inf, V_zero, &V_e_inf, &V_e_zero, x_e, x_w, B_ij, C_ij);
	
	//zurücksetzen der Ausgabewerte
	for(int i=0; i<pd->nstoff-1; i++){
		V_e_dx[i] = 0;
	}

	for(int count=0; count<pd->nstoff-1; count++){
		V_e_dx[count] = V_e_inf_dx[count]+(V_e_zero_dx[count]-V_e_inf_dx[count])*(1-x_w)+(V_e_zero-V_e_inf);
	}

	delete[] V_zero;
	delete[] V_inf;
	delete[] V_e_inf_dx;
	delete[] V_e_zero_dx;
}

void Density::calc_V_ex_deriv(tpropdata *pd,
	double *V_ex_dx, //MethodOutput
	double *x_e, double x_w, 
	double **B_ij, double **C_ij
	)
{
	double *V_zero = new double[pd->nstoff - 1];
	double *V_inf = new double[pd->nstoff - 1];
	double V_e_zero = 0;
	double V_e_inf = 0;
	double V_e = calc_molar_volume(pd, V_inf, V_zero, &V_e_inf, &V_e_zero, x_e, x_w, B_ij, C_ij);

	double *V_e_dx = new double[pd->nstoff - 1];
	calc_V_e_deriv(pd, V_e_dx, x_e, x_w, B_ij, C_ij);

	//zurücksetzen der Ausgabewerte
	for(int i=0; i<pd->nstoff-1; i++){
		V_ex_dx[i] = 0;
	}

	for(int j=0; j<pd->nstoff-1; j++){
		V_ex_dx[j] = V_e_dx[j]*(1-x_w);
		V_ex_dx[j] += V_e-V_inf[j];
	}

	delete[] V_zero;
	delete[] V_inf;
	delete[] V_e_dx;
}

void Density::calc_density_deriv(tpropdata *pd,
	double *density_dx, //MethodOutput
	const double *MethodInputs
	)
{
	double *V_zero = new double[pd->nstoff - 1];
	double *V_inf = new double[pd->nstoff - 1];
	double* x_e	= new double[pd-> nstoff - 1];
	double V_e_zero = 0;
	double V_e_inf = 0;
	double V_e = 0;
	double x_w = 0;
	double temperature = MethodInputs[0];
	double *V_ex_dx = new double[pd->nstoff-1];
	double V_ex = 0;
	double tmp_Vinf = 0;
	double V_w  = 0.01801528*1000 /1000; // (kg/kmol) / (kg/m3)
	double density_dx_water = 0;

	//fills array of binary coefficients
	double** B_ij = new double* [pd->nstoff-1];
	for(int i=0; i < pd->nstoff-1 ; i++)
	{
		B_ij[i] = new double[pd->nstoff-1];
	}
	double** C_ij = new double* [pd->nstoff-1];
	for(int i=0; i < pd->nstoff-1; i++)
	{
		C_ij[i] = new double[pd->nstoff-1];
	}

	fill_binary_coefficents(pd, B_ij, C_ij, temperature);

	//converts molefracs in array
	x_w = MethodInputs[2];
	for(int i=0; i<pd-> nstoff-1; ++i){
		x_e[i] = MethodInputs[i+3];
	} 

	//claculates derivative of exzess volume and molar volume of electrolyte
	calc_V_ex_deriv(pd, V_ex_dx, x_e, x_w, B_ij, C_ij);
	V_e = calc_molar_volume(pd, V_inf, V_zero, &V_e_inf, &V_e_zero, x_e, x_w, B_ij, C_ij);
	
	//calculates exzess volume
	for(int j=0; j<pd->nstoff-1; j++){
		V_ex += (V_e-V_inf[j])*x_e[j];
	}

	//calculates derivative of density 
	for(int count=0; count<pd->nstoff-1; count++){
		tmp_Vinf += x_e[count]*V_inf[count];
	}
	tmp_Vinf += x_w*V_w;

	density_dx[0] = -(1/((tmp_Vinf+V_ex)*(tmp_Vinf+V_ex))) * (V_w); //density derivative of water

	for(int count=1; count<pd->nstoff+1; count++){ //density derivative of electrolyt
		density_dx[count]=-(1/((tmp_Vinf+V_ex)*(tmp_Vinf+V_ex))) * (V_inf[count-1]+V_ex_dx[count-1]);
	}

	//delete memory
	delete[] V_ex_dx;
	delete[] V_zero;
	delete[] V_inf;
	kill_binary_coefficents(pd, B_ij, C_ij);

}

void Density::calc_binary_ms_diffkoeff(
	tpropdata *pd,
	int species_one,
	int species_two,
	double c_binary,
	double *Diff_Koeff
	)
{
	// function to calculate the binary MS diffusion coefficient according to the model published by Visser (phd thesis), p.132ff
	// in Diff_Koeff three coefficients are filled: 0: anion-cation , 1:cation- water and 2: anion-water
	int  i;

	double p_ij[5][3];

	if( !strcmpu(pd->stname[species_one],"H") && !strcmpu(pd->stname[species_two],"Cl")  
		|| !strcmpu(pd->stname[species_one],"Cl") && !strcmpu(pd->stname[species_two],"H") ) 
	{
		double tmp_p_ij [5] [3]  = 
		{
			{0.0, -9.319e-09, -2.042e-09},
			{-1.499E-12, -2.365E-12, -2.667E-13},
			{-3.525E-14, -1.963E-14, -8.223E-15},
			{-2.138E-16, -9.878E-17 -6.180E-17},
			{-7.940E-12, -4.381E-11, -1.367E-11}
		};

		calc_matrix.copy_5_3_array(tmp_p_ij, p_ij);
	}
	else 	if( !strcmpu(pd->stname[species_one],"Na") && !strcmpu(pd->stname[species_two],"OH")  
		|| !strcmpu(pd->stname[species_one],"OH") && !strcmpu(pd->stname[species_two],"Na") ) 
	{
		double tmp_p_ij [5] [3] = 
		{
			{0.0, -1.351E-09, -5.241E-09},
			{-4.062E-13, -8.968E-14, -2.952E-12},
			{-5.383E-15, -1.077E-15, -2.501E-14},
			{-2.153E-19, -2.094E-17, -4.772E-17},
			{-3.318E-12, -7.984E-12, -6.509E-11}
		};
		calc_matrix.copy_5_3_array(tmp_p_ij, p_ij);
	}
	else 	if( !strcmpu(pd->stname[species_one],"Na") && !strcmpu(pd->stname[species_two],"SO4")  
		|| !strcmpu(pd->stname[species_one],"SO4") && !strcmpu(pd->stname[species_two],"Na") ) 
	{
		double tmp_p_ij [5] [3] = 
		{
			{0.0, -1.367E-09, -1.079E-09},
			{-2.943E-14, -1.724E-12, -1.125E-12},
			{-5.839E-16, -6.571E-14, -4.653E-14},
			{-1.592E-19, -6.180E-16, -4.435E-16},
			{-6.080E-13, -2.698E-12, -6.830E-12}
		};
		calc_matrix.copy_5_3_array(tmp_p_ij, p_ij);
	}
	else 	if( !strcmpu(pd->stname[species_one],"H") && !strcmpu(pd->stname[species_two],"SO4")  
		|| !strcmpu(pd->stname[species_one],"SO4") && !strcmpu(pd->stname[species_two],"H") ) 
	{
		double tmp_p_ij [5] [3]=
		{
			{0.0, -9.313E-09 -1.068E-09},
			{-2.545E-13, -1.431E-11, -9.998E-13},
			{-9.089E-16, -2.197E-13, -1.437E-14},
			{-2.240E-17, -1.148E-15, -7.335E-17},
			{0.0, -2.176E-10, -1.525E-11}
		};
		calc_matrix.copy_5_3_array(tmp_p_ij, p_ij);
	}
	else 	if( !strcmpu(pd->stname[species_one],"Na") && !strcmpu(pd->stname[species_two],"Cl")  
		|| !strcmpu(pd->stname[species_one],"Cl") && !strcmpu(pd->stname[species_two],"Na") ) 
	{
		double tmp_p_ij [5] [3]=
		{
			{0.0, 1.336E-09, 2.035E-09},
			{8.018E-14, -3.064E-14, -2.236E-13},
			{-2.090E-16, -3.910E-15, -3.792E-15},
			{-7.028E-18, 3.765E-17, 3.777E-17},
			{2.176E-12, -1.766E-12, 8.317E-12}
		};
		calc_matrix.copy_5_3_array(tmp_p_ij, p_ij);
	}
	else
	{
		double tmp_p_ij [5] [3]=
		{
			{1, 1, 1},
			{1, 1, 1},
			{1, 1, 1},
			{1, 1, 1},
			{1, 1, 1}
		};
		printf( "\n No binary diffusivities for the pair: %s %s",  pd->stname[species_one] ,pd->stname[species_two] );
		printf("\n -->dummy values are used");
		calc_matrix.copy_5_3_array(tmp_p_ij, p_ij);
	}

	for(i=0; i < 3; ++i){
		Diff_Koeff[i] = 
			p_ij[0][i] 
		+  p_ij[1][i] * pow(c_binary, 1) 
			+  p_ij[2][i] * pow(c_binary, 1.5) 
			+  p_ij[3][i] * pow(c_binary, 2)
			+  p_ij[4][i] * pow(c_binary, 0.5);

	}
}

