// initialize1.cpp : fasst die Vorbereitung zur Gemischberechnungen in eNRTL.cpp zusammen
// Autor: Robert Pack, majo (majo02@avt.rwth-aachen.de), thbe (thbe02@avt.rwth-aachen.de)
// UML-Diagramm:
// - 
// Quellen:
// [1] Pack, R., 2011; Development of a Computational Routine for Estimation
//     of Activity Coefficents in Electrolyte Systems, Studienarbeit, AVT.PT
// Vorraussetzungen:
// - Windows (init_enrtl() ermittelt eigenen Ordner)


#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#define CONTROLHEADER
#include "eNRTL.h"
#include "Density.h"

//added by thbe
#include <exception>
//#include <direct.h>
//#include <Windows.h>


tpropdata* glob_pd;
paramstruct* glob_pm;

tpropdata *getglob_pd(){
	return glob_pd;
}

paramstruct getglob_pm(){
	return *glob_pm;
}

extern "C" {
int init_enrtl(const char *filename, int &nSpc){ //Initialisiert die Simulation
// input: Datei im Programmordner mit Gemischparametern
	// open input file
	FILE *input=NULL;
	char* myfile = new char[strlen(filename)+1];
	strcpy(myfile,filename);
 	input = fopen(myfile, "r");

	if(input!=NULL){ //Pruefen ob Datei ladbar, vermeinden des fehlerhaften Dateizugriffs
		tpropdata* pd = new tpropdata; // initialize data structure and read user input 
		init_param(pd, input); //liest benoetigte Gemischparameter ein
		tplibdata* ld = (tplibdata*) pd->libdata;
		paramstruct* pm = new paramstruct; // initialize structure to hold shared calculation between all property methods 
		init_paramstruct(pd, pm); //adressiert Speicher zur Berechnung der weiteren Gemischparameter

		if(glob_pd != NULL){
			delete glob_pd;
		}

		glob_pd = pd;
		glob_pm = pm;

	        nSpc = pd->nstoff;

		fclose(input);

		return 0; //falls Laden der Parameter erfolgreich 
	}else{
		return 1; //falls Fehler auftretten
	}
        delete [] myfile;
}
}

void killstruct()
{
	// kill tpropdata struct
	tpropdata* pd = glob_pd;
	kill_plib(pd);

	// kill paramstruct stuct
	paramstruct* pm = glob_pm; 
	kill_paramstruct(pm);

}

void calculate_therm_factor(
	double *MethodInputs,
	double *Therm_factors
	){
	int i, j, index;

	tpropdata* pd	= glob_pd;
	tplibdata* ld	= (tplibdata*) pd->libdata; 
	paramstruct* pm	= glob_pm;


	// introduce mole fractions from the inputs
	double* x = new double[pd->nstoff];
	for(i=0; i<pd->nstoff; ++i){
		x[i] = MethodInputs[i+2];
	} 

	for(i=0; i<pd->nstoff; ++i){
		for(j=0; j<pd->nstoff; ++j){
			index = j*pd->nstoff + i;
			if(i != j){
				Therm_factors[index] = x[i] * pm->lngm_dx[index];
			}else{
				Therm_factors[index] = 1 + x[i] * pm->lngm_dx[index];
			}
		}
	}
        delete [] x;
}

extern "C" {
void calc_therm_factor_C(
       const int nSpc,
       const double Temp,
       const double Press,
       const double *Mole_frac,
       double *Therm_factors) {
       
	double* MethodInputs;
	tpropdata* pd	= getglob_pd();
	tplibdata* ld	= (tplibdata*) pd->libdata; 
	paramstruct pm	= getglob_pm();

//    printf("in calc_thermo_factor_C \n");
//	printf("nSpc = %d, nstoff = %d \n",nSpc, pd->nstoff);
//	printf("Temp = %f, press = %f \n", Temp, Press);
//	for(int i=0;i<nSpc;i++){
//	        printf("Molefrac, %d, %f \t",i, Mole_frac[i]);
//	}
//	printf("\n");
//        if(nSpc != pd->nstoff) {
//	  printf("ERROR: nspecies defined in fortran and mixture property");
//	  printf("file are different\n");
//          exit(1);
//	}
	MethodInputs = new double[nSpc+2](); // 1 temperature + 1 pressure + number of species
	MethodInputs[0]	= Temp;		// temperature
	MethodInputs[1]	= Press;	// pressure
	for(int i=0;i<nSpc;i++){
		MethodInputs[i+2] = Mole_frac[i];	// mole fractions of the species
        }

        //aktualisiert die Parameterstruktur
	if(update_pm(pd,&pm,MethodInputs[0],&MethodInputs[2])){
	//	printf("Calculate acitity coeff\n");
		get_lngm_x(pd, &pm, pm.lngm, pm.lngm_dt, pm.lngm_dx);
	}

	//for(int i=0;i<nSpc*nSpc;i++){
        //	printf("ACiti, %d, %f \t",i, pm.lngm_dx[i]);
	//}
	//printf("\n");
	// thermodynamic factors
	calculate_therm_factor(MethodInputs, Therm_factors);
//	for(int i=0;i<nSpc*nSpc;i++){
//        	printf("Thermo, %d, %f \t",i, Therm_factors[i]);
//	}
//	printf("\n");
        delete [] MethodInputs;
}
}

extern "C" {
  void calc_ms_diff_matrix_from_molefrac_C(
       const int nSpc,
       const double Temp,
       const double Press,
       const double *Mole_frac,
       double *D_ij_out) {

	double* MethodInputs;
	tpropdata* pd	= getglob_pd();
	tplibdata* ld	= (tplibdata*) pd->libdata; 
	paramstruct pm	= getglob_pm();
	Density calc_edens(false);

	MethodInputs = new double[nSpc+2](); // 1 temperature + 1 pressure + number of species
	MethodInputs[0]	= Temp;		// temperature
	MethodInputs[1]	= Press;	// pressure
	for(int i=0;i<nSpc;i++){
		MethodInputs[i+2] = Mole_frac[i];	// mole fractions of the species
        }

// MS diff coefficients 
	calc_edens.calc_ms_diff_matrix_from_molefrac(pd, MethodInputs, D_ij_out);
        delete [] MethodInputs;
}
}

extern "C" {
  void calc_ms_diff_matrix_from_moledens_C(
       const int nSpc,
       const double Temp,
       const double Press,
       const double *Mole_dens,
       double *D_ij_out) {

	double* MethodInputs;
	tpropdata* pd	= getglob_pd();
	tplibdata* ld	= (tplibdata*) pd->libdata; 
	paramstruct pm	= getglob_pm();
	Density calc_edens(false);

	MethodInputs = new double[nSpc+2](); // 1 temperature + 1 pressure + number of species
	MethodInputs[0]	= Temp;		// temperature
	MethodInputs[1]	= Press;	// pressure
	for(int i=0;i<nSpc;i++){
		MethodInputs[i+2] = Mole_dens[i];	// mole density of the species
        }

// MS diff coefficients 
	calc_edens.calc_ms_diff_matrix_from_moledens(pd, MethodInputs, D_ij_out);
        delete [] MethodInputs;
}
}

void calculate_prop(
	const char *MethodName,	//-in-; gibt die gesuchte Groesse an
	const int *OutputLength, //-in-; Laenge des Ergebnissarrays
	const int *InputLengths, //-in-; Laenge der Randbedingungen MethodInputs
	const double *MethodInputs, //-in-; Randbedingungen; 0: Temperatur; 1: Druck; 2-...: Stoffmengenanteile
	double *MethodOutputs	//-out-; Ergebnissarray
	){

	int i;
	double t = MethodInputs[0];
	double *lngm=NULL, *lngm_dt=NULL, *lngm_dx=NULL;

	// structures to save user input and parameters 
	tpropdata* pd	= glob_pd;
	tplibdata *ld	= (tplibdata*) pd->libdata;
	paramstruct* pm	= glob_pm; 

	double* x = new double[*InputLengths]();
	for(i=0; i<*InputLengths; ++i){ //Uebernimmt die Stoffmengen aus der Eingabe
		x[i] = MethodInputs[i+2];
	}

	//definiert das eNRTL-Modell als symetrisch
	pd->refined_enrtl = FALSE; 
	ld->enrtlref = 0;

	for(i=0; i<*OutputLength; ++i){ //leer den Speicher der Ausgabe
		MethodOutputs[i] = 0.; 
	}

	if (strcmp(MethodName, "LIQUIDFUGACITY") == 0){
		if(update_pm(pd, pm, t, x)){ //aktuallisiert Parameterliste und gibt bei Erfolg true aus
			lngm = pm->lngm;
		}else{
			if(!pm->lngm_curr){ //prueft ob Aktivitaetskoeffizien nicht berechnet wurde
				lngm = pm->lngm;
			}
		}

		get_lngm_x(pd, pm, lngm, lngm_dt, lngm_dx); //Berechnet den Aktivitaetskoeffizienten und die Ableitungen

		for(i=0; i<pd->nstoff; ++i){ //gibt exp() des Aktivitaetskoeffizienten an
			MethodOutputs[i] = exp(pm->lngm[i]);
		}

	}
	/*else if (strcmp(MethodName, "LIQUIDEXESSGIBBSFREEENERGY") == 0)
	{
		if(update_pm(pd, pm, t, x)) lngm = pm->lngm;
		else if(!pm->lngm_curr) lngm = pm->lngm;

		get_lngm_x(pd, pm, lngm, lngm_dt, lngm_dx);

		for(i=0; i<pd->nstoff; ++i){
			*MethodOutputs += x[i] * pm->lngm[i]; 
		}

		*MethodOutputs *= RALG * t; 


	}
	else if (strcmp(MethodName, "LIQUIDEXCESSENTHALPY") == 0)
	{
		if(update_pm(pd, pm, t, x)) lngm_dt = pm->lngm_dt;
		else if(!pm->lngm_dt_curr) lngm_dt = pm->lngm_dt;

		get_lngm_x(pd, pm, lngm, lngm_dt, lngm_dx);

		for(i=0; i<pd->nstoff; ++i){
			*MethodOutputs += x[i] * pm->lngm_dt[i]; 
		}

		*MethodOutputs *= -RALG * pow(t, 2); 


	}
	else if (strcmp(MethodName, "LIQUIDEXCESSENTROPY") == 0)
	{
		if(update_pm(pd, pm, t, x)){
			lngm	= pm->lngm;
			lngm_dt	= pm->lngm_dt;
		}else{
			if(!pm->lngm_curr)		lngm	= pm->lngm;
			if(!pm->lngm_dt_curr)	lngm_dt	= pm->lngm_dt;
		}

		get_lngm_x(pd, pm, lngm, lngm_dt, lngm_dx);

		for(i=0; i<pd->nstoff; ++i){
			*MethodOutputs += x[i] * (pm->lngm[i] + t* pm->lngm_dt[i]); 
		}

		*MethodOutputs *= -RALG; 

	}*/

	delete [] x;
}

void calculate_prop_deriv(
	const char *MethodName,
	const int *OutputLength,
	const int *InputLengths, 
	const double *MethodInputs,
	const char *DInputName,
	const int *DInputLength,
	double *DerivOutputs 
	){
	int i, j, index;
	double t = MethodInputs[0], value=0.; 

	double *lngm=NULL, *lngm_dt=NULL, *lngm_dx=NULL;

	tpropdata* pd	= glob_pd;
	tplibdata* ld	= (tplibdata*) pd->libdata; 
	paramstruct* pm	= glob_pm;


	double* x = new double[*InputLengths]();
	for(i=0; i<*InputLengths; ++i){ //uebernimmt die Stoffmenge aus der Eingabe
		x[i] = MethodInputs[i+2];
	} 

	for(i=0; i<*OutputLength; ++i){ //leer Speicher der Ausgaben
		DerivOutputs[i] = 0.; 
	}

	if (strcmp(MethodName, "LIQUIDFUGACITY") == 0){
		if (strcmp(DInputName, "TEMPERATURE") == 0){
			if(update_pm(pd, pm, t, x)){ //aktuallisiert Parameterliste und gibt bei Erfolg true aus
				lngm_dt = pm->lngm_dt;
			}else{
				if(!pm->lngm_dt_curr){ //prueft ob Aktivitaetskoeffizien nicht berechnet wurde
					lngm_dt = pm->lngm_dt;
				}
			}

			get_lngm_x(pd, pm, lngm, lngm_dt, lngm_dx); //berechnet Aktivitaetskoeffizienten

			for(i=0; i<pd->nstoff; ++i){ //uebergibt Ableitung ueber der Temperatur an die Ausgabe
				DerivOutputs[i] = pm->lngm_dt[i]; 
			}
		}

		if (strcmp(DInputName, "PRESSURE") == 0){
			for(i=0; i<pd->nstoff; ++i){
				DerivOutputs[i] = 0.;
			}
		}

		if (strcmp(DInputName, "COMPOSITION") == 0){
			if(update_pm(pd, pm, t, x)){ //aktuallisiert Parameterliste und gibt bei Erfolg true aus
				lngm_dx = pm->lngm_dx;
			}else{
				if(!pm->lngm_dx_curr){ //prueft ob Aktivitaetskoeffizien nicht berechnet wurde
					lngm_dx = pm->lngm_dx;
				}
			}

			get_lngm_x(pd, pm, lngm, lngm_dt, lngm_dx); //berechnet Aktivitaetskoeffizienten

			for(i=0; i<(pd->nstoff)*(pd->nstoff); ++i){ //uebergibt Ableitung ueber der Stoffmenge an die Ausgabe
				DerivOutputs[i] = pm->lngm_dx[i]; 
			}

		}
	}
	/*else if (strcmp(MethodName, "LiquidExcessGibbsFreeEnergy") == 0)
	{
		if (strcmp(DInputName, "TEMPERATURE") == 0)
		{
			if(update_pm(pd, pm, t, x)){
				lngm	= pm->lngm;
				lngm_dt	= pm->lngm_dt;
			}else{
				if(!pm->lngm_curr)		lngm	= pm->lngm;
				if(!pm->lngm_dt_curr)	lngm_dt	= pm->lngm_dt;
			}

			get_lngm_x(pd, pm, lngm, lngm_dt, lngm_dx);

			for(i=0; i<pd->nstoff; ++i){
				*DerivOutputs	+= x[i] * pm->lngm_dt[i];
				value			+= x[i] * pm->lngm[i];
			}

			*DerivOutputs *= (RALG * t); 
			*DerivOutputs += (RALG * value); 


		}
		else if (strcmp(DInputName, "PRESSURE") == 0)
		{
			DerivOutputs[0] = 0.;

		}
		else if (strcmp(DInputName, "COMPOSITION") == 0)
		{
			if(update_pm(pd, pm, t, x)){
				lngm	= pm->lngm;
				lngm_dx	= pm->lngm_dx;
			}else{
				if(!pm->lngm_curr)		lngm	= pm->lngm;
				if(!pm->lngm_dx_curr)	lngm_dx	= pm->lngm_dx;
			}

			get_lngm_x(pd, pm, lngm, lngm_dt, lngm_dx);

			for(i=0; i<pd->nstoff; ++i){
				for(j=0; j<pd->nstoff; ++j){
					index = j*pd->nstoff + i;
					if(i != j) DerivOutputs[i]	+= x[j] * pm->lngm_dx[index];
					else DerivOutputs[i] += x[j] * pm->lngm_dx[index] + pm->lngm[i];
				}
			}

		}
	}*/
}

void init_plib(tpropdata *pd){ //erstellt die Parameterstruktur und setzt Initialwerte(NULL, nur bei Tref=ENRTLTREF)
	tplibdata *ld;

	pd->libdata = new tplibdata;

	ld = (tplibdata*) pd->libdata;

	pd->stname		= NULL;
	pd->refined_enrtl = FALSE;

	ld->born		= NULL;
	ld->charge		= NULL;
	ld->molec_v		= NULL;
	ld->ani_v		= NULL;
	ld->cati_v		= NULL;
	//ld->solid_v		= NULL;
	ld->enrtlref	= 0;
	ld->henry		= NULL;
	ld->nseg_v		= NULL;
	ld->nseg		= -1;
	ld->mixmolvol	= 1;
	ld->defgmelc	= 1;

	ld->gmelcc		= NULL;
	ld->gmelcd		= NULL;
	ld->gmelce		= NULL;
	ld->gmelcn		= NULL;

	ld->nrtla		= NULL;
	ld->nrtlta		= NULL;
	ld->nrtltb		= NULL;

	ld->cpdiec		= NULL;

	ld->mm			= NULL;
	ld->tc			= NULL;
	ld->pc			= NULL;
	ld->vc			= NULL;
	ld->zc			= NULL;

	ld->rktzra		= NULL;
	ld->vcrkt		= NULL;

	ld->dG_0f		= NULL; 
	ld->dH_0f		= NULL; 
	//	ld->poly_cp_a	= NULL; 
	//	ld->cc_cp_b		= NULL; 
	ld->stoic_coeff	= NULL; 

	ld->nreac		= 0; 
	ld->Tref		= ENRTLTREF; 
}

void init_param(tpropdata *pd, FILE *input){ //liest die benoetigten Parameter aus der uebergebenen Datei
	char row[LSTRING], *field[10], *command=NULL;
	int ident, testcase=FALSE, columns, *species=NULL, i, j, m, index;
	int chargecount=0; //Anzahl der Electrolyte-Molecule Pair Parameter
	tplibdata *ld;

	init_plib(pd); //erstellt die Parameterstruktur und setzt Startwerte
	ld = (tplibdata*) pd->libdata;

	while(!feof(input)){ //Durchlauf bis zum Ende der Datei oder finden von "BEGIN" -> alles vor "BEGIN" wird ignoriert
		fgetz(row, LSTRING, input); //liest naechste Zeile ein
		command = strtok(row, INPUTDELI);

		if(command){
			if(strcmpu(command, "BEGIN") == 0){
				break;  //Schleifenabbruch bei "BEGIN"
			}
			if(strcmpu(command, "TEST") == 0){
				testcase = TRUE;
			}
		}
	}

	while(!feof(input)){ //Durchlauf bis zum Dateiende
		if((columns=readline(input, row, command, field))>0){ //Uebergabe: Anzahl der ? und der eingelesenen Zeile(row), des Befehls(command), der Werte(field)
			ident = FALSE; //Pruefbit, falls Anweisung gefunden gleich true

			if(!ident && !strcmpu(command,"STNAME")){ //Namen der Komponenten setzen
				if(columns>0){
					if(get_komp(pd,field[0],&i)==TRUE){
						delete[] pd->stname[i]; //evtl. vorhandenen Speicher wieder bereitstellen
						pd->stname[i]= new char[strlen(field[1])+1]; //neuen Speicher benutzen
						strcpy(pd->stname[i],field[1]); //Speicher fuellen
					}
				}
				ident=TRUE;
			}

			if(!ident && !strcmpu(command,"NSTOFF")){ //Anzahl der Komponenten setzen
				if(columns>0){
					int n;		
					sscanf(field[0],"%d",&n);
					ident=TRUE;
					if(n>0){
						pd->nstoff=n;
						init_plib_data(pd);	//erzeugt die Parameterstruktur auf Grundlage der Komponentenanzahl
					}
				}
			}

			if(!ident && !strcmpu(command,"NREAC")){ //Anzahl der Reaktionen(?) einlesen
				if(columns>0){
					int n;		
					sscanf(field[0],"%d",&n);
					ident=TRUE;
					if(n>0) {
						ld->nreac=n;
						init_reac(pd); //erstellt die Parameterstruktur auf Grundlage der Reaktionszahlen
					}
				}
			}

			if(!ident && !strcmpu(command,"NSOLID")){ //Anzahl der festen Komponenten einlesen 
					if(columns>0){
						int n;		
						sscanf(field[0],"%d",&n);
						ident=TRUE;
						if(n>0) {
							ld->nsolid=n;
						}
					}
				}

			if(!ident && !strcmpu(command,"SEG")){ //? Bedeutung
				double wert;
				if(columns>0){
					if(get_komp(pd,field[0],&i)==TRUE){
						if(columns>1){
							for (j=0; j<ld->nseg; ++j){
								sscanf(field[j+1],"%lf",&wert);
								ld->nseg_v[i][j] = wert;
							}
							ident = TRUE;
						}
					}	
				}
			}

			if(!ident && !strcmpu(command,"STOIC")){ //liest stoechiometrische Koeffizienten(?)
				double wert;
				if(columns>0){
					sscanf(field[0],"%lf",&i);
					if(columns > pd->nstoff){
						for(j=0; j<pd->nstoff; ++j){
							index = i * pd->nstoff + j; 
							sscanf(field[j+1],"%lf",&ld->stoic_coeff[index]); 
						}
						ident = TRUE;
					}	
				}
			}

			if(!ident && !strcmpu(command,"SOLID")){ //definiert einen Stoff als Feststoff
				double wert;
				if(columns>0){
					for(j=0; j<columns; ++j){
						if(get_komp(pd,field[j],&i)==TRUE){
							ld->solid[i] = TRUE; 
						}
					}
					ident = TRUE;	
				}
			}

			if(!ident && !strcmpu(command,"HENRY")){ //deklariert einzelne Komponente als Henry-Komponente
				double wert;
				if(columns>0){
					for(j=0; j<columns; ++j){
						if(get_komp(pd,field[j],&i)==TRUE){
							ld->henry[i] = TRUE; 
						}
					}
					ident = TRUE;	
				}
			}

			if(!ident && !strcmpu(command,"BORN")){ //setzt die Born-Korrektur
				if (ld->born == NULL){
					init_plib_segment(pd);
				}

				if(columns>0){
					if(sscanf(field[0],"%d",&i)==1){
						if((i<ld->nseg && i>=0) || (i<pd->nstoff && i>=0)){
							if(columns>1){
								sscanf(field[1],"%lf",&ld->born[i]);
							}
						}
					}	
				}
				ident=TRUE;
			}

			if(!ident && !strcmpu(command,"DGFORM")){ //Freie Bildungsenthalpie
				if(columns>0){
					if(sscanf(field[0],"%d",&i)==1){
						if(columns>1) {
							sscanf(field[1],"%lf",&ld->dG_0f[i]);
						}
					}	
				}
				ident=TRUE;
			}

			if(!ident && !strcmpu(command,"DHFORM")){ //Bildungsenthalpie
				if(columns>0){
					if(sscanf(field[0],"%d",&i)==1){
						if(columns>1){
							sscanf(field[1],"%lf",&ld->dH_0f[i]);
						}
					}	
				}
				ident=TRUE;
			}

			if(!ident && !strcmpu(command,"CPDIEC")){ //parameters for dielectric constant of spiecies
				if(columns>0){
					if(get_komp(pd,field[0],&i)==TRUE){
						if(columns>3){
							sscanf(field[1],"%lf",&ld->cpdiec[i*3]);
							sscanf(field[2],"%lf",&ld->cpdiec[i*3+1]);
							sscanf(field[3],"%lf",&ld->cpdiec[i*3+2]);
						}
					}	
				}
				ident=TRUE;
			}

			if(!ident && !strcmpu(command,"CHARGE")){ //setzt die Ladungen der Anionen und Kathionen
				if (ld->charge == NULL){
					init_plib_segment(pd);
				}

				if(columns>0){
					if(sscanf(field[0],"%d",&i)==1){
						if((i<ld->nseg && i>=0) || (i<pd->nstoff && i>=0)){
							if(columns>1){
								sscanf(field[1],"%i",&ld->charge[i]);
								ident=TRUE;
								chargecount += 1; 
							}	
						}
					}	
				}
			}

			if (((chargecount == pd->nstoff) && (ld->nseg == -1)) || ((chargecount == ld->nseg) && (ld->nseg != -1))){
			//prueft ob die Anzahl der Electrolyte-Molecule Pair Parameter gleich der Anzahl der Komponenten und ? (ld->nseg==-1) oder
			// die Anzahl der Electrolyte-Molecule Pair Parameter gleich ? und nicht -1 ist
				ld->molec = 0;
				ld->cati = 0;
				ld->ani = 0;

				for (i=0; i<chargecount; i++){	//ordnet die Komponenten in Kationen, Anionen und neutralen Komponeten ein 
					if(ld->charge[i]==0){ //falls Komponente neutral
						ld->molec += 1; //erhoeht die Summe der neutralen Komponenten um eins
						ld->molec_v = (int*) realloc(ld->molec_v, ld->molec*sizeof(int)); //erweitert den Speicher der neutralen Komponeten
						ld->molec_v[ld->molec-1] = i; //ordnet der neutralen Komponente eine Gemischkomponente zu
					}else{
						if(ld->charge[i]>0){ //falls Kationen
							ld->cati += 1;
							ld->cati_v = (int*) realloc(ld->cati_v, ld->cati*sizeof(int)); 
							ld->cati_v[ld->cati-1] = i; 
						}else{
							if(ld->charge[i]<0){ //falls Anionen
								ld->ani += 1;
								ld->ani_v = (int*) realloc(ld->ani_v, ld->ani*sizeof(int));
								ld->ani_v[ld->ani-1] = i;
							}
						}
					}
				}
					
				chargecount = 2*ld->molec*ld->cati*ld->ani + ld->ani*ld->cati*(ld->cati-1) + ld->cati*ld->ani*(ld->ani-1);

				ld->gmelcc	= new double[chargecount]();
				ld->gmelcd	= new double[chargecount]();
				ld->gmelce	= new double[chargecount]();
				ld->gmelcn	= new double[chargecount]();

				ld->nrtlta	= new double[ld->molec*(ld->molec-1)]();
				ld->nrtltb	= new double[ld->molec*(ld->molec-1)]();
				ld->nrtla	= new double[ld->molec*(ld->molec-1)]();

				for (i=0; i<chargecount; i++){ //setzt die neuerstellten Variablen auf 0
					ld->gmelcc[i] = 0.0;
					ld->gmelcd[i] = 0.0;
					ld->gmelce[i] = 0.0;
					ld->gmelcn[i] = 0.0;
				}

				chargecount = 0;
			}

			if(!ident && !strcmpu(command,"GMELCC")){ //liest einen Electrolyte-Moleecule Pair Parameter
				if(!testcase){
					species = (int*) realloc(species, (columns-1)*sizeof(int)); 
					for(i=0; i<(columns-1); i++){
						sscanf(field[i],"%i",&species[i]);	//liest die Ladungskonfiguration
					}
					get_index_input(ld, species, &m, columns-1); //gibt den Index des zu setzenden Parameters aufgrund der Ladungskonfiguration zurueck
					sscanf(field[columns-1],"%lf",&ld->gmelcc[m]); //setzt den Parameter
					ident=TRUE;
				}
			}

			if(!ident && !strcmpu(command,"GMELCD")){ //sliest einen Electrolyte-Moleecule Pair Parameter
			//Vorgehen aehnlich zu "GMELCC"
				if(!testcase){
					species = (int*) realloc(species, (columns-1)*sizeof(int));
					for(i=0; i<(columns-1); i++){
						sscanf(field[i],"%i",&species[i]);
					}
					get_index_input(ld, species, &m, columns-1);
					sscanf(field[columns-1],"%lf",&ld->gmelcd[m]);
					ident=TRUE;
				}
			}

			if(!ident && !strcmpu(command,"GMELCE")){ //liest einen Electrolyte-Moleecule Pair Parameter
			//Vorgehen aehnlich zu "GMELCC"
				if(!testcase){
					species = (int*) realloc(species, (columns-1)*sizeof(int));
					for(i=0; i<(columns-1); i++){
						sscanf(field[i],"%i",&species[i]);
					}
					get_index_input(ld, species, &m, columns-1);
					sscanf(field[columns-1],"%lf",&ld->gmelce[m]);
					ident=TRUE;
				}
			}	

			if(!ident && !strcmpu(command,"GMELCN")){ //liest einen Electrolyte-Moleecule Pair Parameter
			//Vorgehen aehnlich zu "GMELCC"
				if(!testcase){
					species = (int*) realloc(species, (columns-1)*sizeof(int));
					for(i=0; i<(columns-1); i++){
						sscanf(field[i],"%i",&species[i]);
					}
					get_index_input(ld, species, &m, columns-1);
					sscanf(field[columns-1],"%lf",&ld->gmelcn[m]);
					ident=TRUE;
				}
			}

			if(!ident && !strcmpu(command,"NRTL")){ //? Bedeutung
				if(columns>1){
					if(get_komp(pd,field[0],&i)==TRUE && get_komp(pd,field[1],&j)==TRUE){
						if(i != j){
							species = (int*) realloc(species, 2*sizeof(int));
							species[0] = i;
							species[1] = j;
							get_index_input(ld, species, &index, 2);
							if(columns>2) sscanf(field[2],"%lf",&ld->nrtlta[index]);
							if(columns>3) sscanf(field[3],"%lf",&ld->nrtltb[index]);
							if(columns>4) sscanf(field[4],"%lf",&ld->nrtla[index]);
						}
					}
				}
				ident=TRUE;
			}

			if(!ident && !strcmpu(command,"KRIT")){ //liest kritische Stoffparameter: Temperatur, Druck, Volumen, Komprimierbarkeit
				if(columns>0){
					if(get_komp(pd,field[0],&i)==TRUE){
						if(columns>1) sscanf(field[1],"%lf",&ld->tc[i]);
						if(columns>2) sscanf(field[2],"%lf",&ld->pc[i]);
						if(columns>3) sscanf(field[3],"%lf",&ld->vc[i]);
						if(columns>4) sscanf(field[4],"%lf",&ld->zc[i]);
					}
				}
				ident=TRUE;
			}

			if(!ident && !strcmpu(command,"RKTZRA")){ //liest Komprimierbarkeit nach dem Rackett-Model
				if(columns>0){
					if(get_komp(pd,field[0],&i)==TRUE){
						if(columns>1) sscanf(field[1],"%lf",&ld->rktzra[i]);
					}
				}
				ident=TRUE;
			}

			if(!ident && !strcmpu(command,"VCRKT")){ //liest kritisches Volumen nach dem Rackett-Model
				if(columns>0){
					if(get_komp(pd,field[0],&i)==TRUE){
						if(columns>1) sscanf(field[1],"%lf",&ld->vcrkt[i]);
					}
				}
				ident=TRUE;
			}

			if(!ident && !strcmpu(command, "MM")){ //liest die molaren Massen ein
				if(columns>0){
					if(get_komp(pd, field[0], &i)){
						sscanf(field[1],"%lf",&ld->mm[i]);
					}
				}
			}
		}
	}

	// set default parameters for missing parameters
	if (ld->defgmelc) {
		set_param_default(pd);
	}

        free(species);
}

int get_komp(tpropdata *pd, char *string, int *komp){ //sucht nach einer Komponente in der Parameterstruktur
//gibt das vorhandensein einer Komponente(Name oder Nummer im Komponentenarray) an(true/false), ausser bei "all" und "-1", dann -1
	if(strcmpu(string,"all")==0){
		*komp=ALL;
		return ALL;
	}

	for(int i=0;i<pd->nstoff;i++){
		if(strcmpu(string,pd->stname[i])==0){
			*komp=i;
			return TRUE;
		} 
	}

	if(sscanf(string,"%d",komp)==1){
		if(*komp<pd->nstoff && *komp>=0){
			return TRUE;
		}
		if(*komp==ALL){
			return ALL;
		}
	}

	return FALSE;
}

static int readline(FILE *ein, char *zeile, char *befehl, char **feld){ //liest naechste Zeile aus einer Datei
//ein: Zeiger auf Datei, zeile: ganze eingelesene Zeile, befehl: nur der befehl in der zeile, feld: anweisung in der Zeile
//Rueckgabewert ist die Groesse des Feldes
	char   *rest;
	int    anzahl; 

	fgetz(zeile,LSTRING,ein); //liest naechste Zeile aus der Datei *ein

	befehl=mystrtok(zeile,INPUTDELI); //gibt den Befehl in der Zeile zurueck

	if(!befehl){
		return(-1);
	}

	/* Zeiger auf den rest der Zeile erzeugen */
	rest=mystrtok(NULL,NULL); //??? =NULL, warum nicht sofort rest=""
	if(!rest){
		rest="";
	}

	if(strcmpu(befehl,"MESSAGE")==0){ //Ausgabe einer Nachricht
		printf("\nprops: %s",rest);
		return(0);
	}

	anzahl=0;
	while(1){
		feld[anzahl]=mystrtok(NULL,INPUTDELI); //?
		if(feld[anzahl]){
			anzahl++;
		}else{ 
			break;
		}
		if(anzahl==10){
			break;
		}
	}

	return anzahl;
}

char *fgetz(char *s, int n, FILE *stream){ //liest naechste Zeile aus einer Datei
//s: Speicher fuer eingelesenen Daten; n: Anzahl der eingelesenen Zeichen; stream: Dateiname
	char *result,*ptr,*tmp_ptr;

	*s='\0';
	result=fgets(s,n,stream);

	if(result!=NULL) { //pruefen ob Zeichen vorhanden
		/* remove comments or newline */
		if((ptr=strpbrk(s,";\n\r"))!=NULL){ //suche nach ";", Zeilenumbruch oder Wagenruecklauf
			*ptr='\0'; //ersetzen der Zeichen durch das Ende des Strings
		}
		/* read continuation */
		if((ptr=strstr(s,"..."))!=NULL){ //suchen nach Zeichensatz "..."
			if(!(tmp_ptr=fgetz(ptr,n-(int)(ptr-s),stream))){ //pruefen ob naechste nicht vorhanden ist
				*ptr='\0'; //ersetzen der Zeichen durch das Ende des Strings 
			}
		}
	} 
	return result;
}

static char* mystrtok(char *input, char *deli){
	static char *next=NULL;
	char *res;

	if(!deli){
		return next;
	}

	if(input){ //durchsuchen nach angegebenem Zeichen, falls zu durchsuchender Zeichensatz vorhanden
		next=input;
		while(strchr(deli, *next) && *next!='\0'){ //prueft bis deli ungleich input[x] ist; bei deli="\t" wird ein "\t" am anfang von input zukuenftig ignoriert
			next++;
		}
		if(*next=='\0'){ //pruefen ob Ende des Strings
			next=NULL;
		}
	}

	res=next;
	if(!res){
		return NULL;
	}

	if(*res=='\"'){ //falls  '"'
		res++; 
		next=strchr(res, '\"'); //suche Gegenstueck zum ersten '"'
		if(!next){
			next=NULL;
		}else{
			*next='\0';
			next++;
			while(strchr(deli, *next) && *next!='\0'){ //ueberprueft ob betrachtetes Zeichen gleich zu suchendes Zeichen ist
				next++;
			}
		}
	}else{	
		while(!strchr(deli, *next) && *next!='\0'){ //sucht bis zum ersten deli in next
			next++;
		}
		if(*next!='\0'){ //setzt das Stringende an die gefundene Position
			*next='\0';
			next++;
			while(strchr(deli, *next) && *next!='\0'){ //ueberprueft ob betrachtetes Zeichen gleich zu suchendes Zeichen ist
				next++;
			}
		}
	}

	if(*next=='\0') next=NULL;

	return res;
}

void init_plib_segment(tpropdata *pd){ //erstellt die Parameterstruktur auf Grundlage der Segmentanzahl
//Initialisierung von pd->libdata siehe init_plib() und init_plib_data(..)
	tplibdata *ld = (tplibdata*) pd->libdata;
	int segcount, i, j;

	if(ld->nseg == -1){
		segcount = pd->nstoff;	
	}else{
		segcount = ld->nseg;
	}

	ld->charge	=  new int[segcount](); //(int*) malloc(segcount*sizeof(int));
	//ld->enrtlref=  new int[segcount]; //(int*) malloc(segcount*sizeof(int));
	ld->born	=  new double[segcount](); //(double*) malloc(segcount*sizeof(double));
	ld->nseg_v	=  new double*[pd->nstoff]; 
	//ld->nseg_v	= (double**) get_matrix(pd->nstoff, segcount, sizeof(double));

	for(i=0; i<pd->nstoff; ++i) ld->nseg_v[i] = new double[segcount]();

	for(i=0; i<segcount; i++){
		ld->born[i]		= 0.;
		//ld->enrtlref[i]	= 0;
		ld->charge[i]	= 0;

		for(j=0; j<pd->nstoff; j++){
			ld->nseg_v[j][i]  = 0.;
		}
	}
}

void init_plib_data(tpropdata *pd){ //erstellt die Parameterstruktur auf Grundlage der Komponentenzahl
//Initialisierung von pd->libdata siehe init_plib()
	int i;
	tplibdata *ld = (tplibdata*) pd->libdata;

	ld->tc		=	new double[pd->nstoff]();
	ld->pc		=	new double[pd->nstoff]();
	ld->vc		=	new double[pd->nstoff]();
	ld->zc		=	new double[pd->nstoff]();
	ld->mm		=	new double[pd->nstoff]();
	ld->cpdiec	=	new double[3*pd->nstoff]();
	ld->rktzra	=	new double[pd->nstoff]();
	ld->vcrkt	=	new double[pd->nstoff]();
	ld->henry	=	new int[pd->nstoff]();
	ld->solid	=	new int[pd->nstoff]();
	ld->dG_0f	=	new double[pd->nstoff]();
	ld->dH_0f	=	new double[pd->nstoff]();
	ld->cc_cp_aq	=	new double[2*pd->nstoff]();
	ld->poly_cp_aq	=	new double[pd->nstoff*8]();

	for(i=0; i<pd->nstoff; ++i){
		ld->henry[i] = FALSE;
		ld->solid[i] = FALSE;
	}

	pd->stname	= new char*[pd->nstoff];

	for(i=0; i<pd->nstoff; i++){
		pd->stname[i] = new char[strlen("Stoff")+(int)log10((double)(i+1)+2)];
	}
}

void init_reac(tpropdata *pd){ //erstellt die Parameterstruktur fuer die stoechiometrischen Koeffizienten
	tplibdata *ld = (tplibdata*) pd->libdata; 
	ld->stoic_coeff	= new double[ld->nreac*pd->nstoff]();
}  

void kill_plib(tpropdata *pd){ //loescht die Parameterstruktur
	int i;
	tplibdata *ld = (tplibdata*) pd->libdata;

	if(ld){
		delete [] ld->born;
		delete [] ld->charge;
		delete [] ld->molec_v;
		delete [] ld->ani_v;
		delete [] ld->cati_v;
		//		delete [] ld->solid_v; 

		delete [] ld->henry;
		delete [] ld->gmelcc;
		delete [] ld->gmelcd;
		delete [] ld->gmelce;
		delete [] ld->gmelcn;

		delete [] ld->mm;
		delete [] ld->tc;
		delete [] ld->pc;
		delete [] ld->vc;
		delete [] ld->zc;

		delete [] ld->cpdiec;

		delete [] ld->nrtla;
		delete [] ld->nrtlta;
		delete [] ld->nrtltb;

		delete [] ld->rktzra;
		delete [] ld->vcrkt;

		delete [] ld->dG_0f; 
		delete [] ld->dH_0f; 
		delete [] ld->stoic_coeff; 
		//		delete [] ld->cc_cp_a; 
		//		delete [] ld->cc_cp_b; 

		for(i=0; i<pd->nstoff; ++i) delete [] ld->nseg_v[i];

		delete ld;

		for(i=0; i<pd->nstoff; ++i) delete [] pd->stname[i];

		delete [] pd->stname;
		delete pd;
	}
} 

int strcmpu(const char *string1, const char *string2){ //Vergleicht zwei Strings unabhaengig von gross-/kleinscheibung
// und gibt Anzahl der gleichen Zeichen wieder
	int  res;
	size_t i, l1=strlen(string1),l2=strlen(string2);
	char *hlp1,*hlp2;

	hlp1= (char*) malloc((l1+l2+2)*sizeof(char));
	hlp2=hlp1+l1+1;

	for(i=0;i<=l1;i++){
		hlp1[i]=toupper(string1[i]); /* String mit term-0 */
	}
	for(i=0;i<=l2;i++){
		hlp2[i]=toupper(string2[i]); /* kopieren!         */
	}

	res=strcmp(hlp1,hlp2); 

	free(hlp1);

	return(res);
}

void* get_matrix(int m, int n, int el_size){
	void **ptr;
	int i, offset, size;

	if(m*n <= 0) return NULL;

	if((m*sizeof(void*)) % el_size )	offset = m*sizeof(void*)/el_size +1;
	else								offset = m*sizeof(void*)/el_size;

	size = (offset + m*n) * el_size;
	ptr  = (void**) malloc(size);

	if(!ptr) return(ptr);

	ptr[0]						= (void*) ((char*) ptr		+ offset*el_size);
	for(i=1; i<m; i++) ptr[i]	= (void*) ((char*) ptr[i-1] + n		*el_size);

	return((void*)ptr);
}

void init_paramstruct(tpropdata *pd, paramstruct *pm){ //Initialisiert die Parameterstruktur und setze Startwerte
	int param_count, param_count_dx, i; 
	tplibdata* ld = (tplibdata*) pd->libdata;

	param_count	= ld->molec*((ld->molec -1)+2*(ld->cati+ld->ani)) + 2*(ld->cati+ld->ani);
	param_count_dx	= 4*ld->molec*ld->cati*ld->ani + ld->cati*ld->ani*ld->cati 
		+ ld->ani*ld->cati*ld->ani;
	int param_input	= ld->molec*((ld->molec -1)+2*(ld->cati*ld->ani)) + ld->cati*ld->ani*(ld->ani-1)
		+ ld->ani*ld->cati*(ld->cati-1);

	pm->alpha			= new double[param_count]();
	pm->alpha_dx		= new double[param_count_dx]();
	pm->alpha_ca_dx		= new double[param_count_dx]();
	pm->alpha_input		= new double[param_input]();
	pm->alpha_dt_input	= new double[param_input]();

	pm->gg				= new double[param_count]();
	pm->gg_dt			= new double[param_count]();
	pm->gg_dx			= new double[param_count_dx]();
	pm->gg_dt_dx		= new double[param_count_dx]();
	pm->gg_ca_dx		= new double[param_count_dx]();
	pm->gg_input		= new double[param_input]();
	pm->gg_dt_input		= new double[param_input]();

	pm->tau				= new double[param_count]();
	pm->tau_dt			= new double[param_count]();
	pm->tau_dx			= new double[param_count_dx]();
	pm->tau_dt_dx		= new double[param_count_dx]();
	pm->tau_ca_dx		= new double[param_count_dx]();
	pm->tau_input		= new double[param_input]();
	pm->tau_dt_input	= new double[param_input]();

	for(i=0; i<param_count; ++i){
		pm->alpha[i]	= 0.;
		pm->gg[i]		= 0.;
		pm->gg_dt[i]	= 0.;
		pm->tau[i]		= 0.;
		pm->tau_dt[i]	= 0.;
	}

	for(i=0; i<param_count_dx; ++i){
		pm->alpha_dx[i]		= 0.;
		pm->alpha_ca_dx[i]	= 0.;
		pm->gg_dx[i]		= 0.;
		pm->gg_dt_dx[i]		= 0.;
		pm->gg_ca_dx[i]		= 0.;
		pm->tau_dx[i]		= 0.;
		pm->tau_dt_dx[i]	= 0.;
		pm->tau_ca_dx[i]	= 0.;
	}

	pm->last_x_eff		= new double[ld->molec+ld->cati+ld->ani]();
	pm->last_x			= new double[pd->nstoff]();
	pm->chrfrac			= new double[ld->cati+ld->ani]();
	pm->chrfrac_dx		= new double[ld->cati*ld->cati+ld->ani*ld->ani]();

	for(i=0; i<ld->molec+ld->cati+ld->ani; ++i){
		pm->last_x_eff[i]= 0.;
	}

	param_count_dx		= ld->molec*(ld->molec+ld->cati+ld->ani) 
		+ ld->cati*(ld->molec+ld->ani) + ld->ani*(ld->molec+ld->cati);

	pm->sum_xg			= new double[ld->molec+ld->cati+ld->ani]();
	pm->sum_xg_dt		= new double[ld->molec+ld->cati+ld->ani]();
	pm->sum_xg_dx		= new double[param_count_dx]();
	pm->sum_xg_dt_dx	= new double[param_count_dx]();
	pm->sum_xg_mx		= new double[ld->molec]();
	pm->sum_xg_mx_dt	= new double[ld->molec]();
	pm->sum_xg_ca		= new double[ld->molec+ld->cati+ld->ani]();
	pm->sum_xg_ca_dt	= new double[ld->molec+ld->cati+ld->ani]();

	pm->sum_xgt			= new double[ld->molec+ld->cati+ld->ani]();
	pm->sum_xgt_dt		= new double[ld->molec+ld->cati+ld->ani]();
	pm->sum_xgt_dx		= new double[param_count_dx]();
	pm->sum_xgt_dt_dx	= new double[param_count_dx]();
	pm->sum_xgt_mx		= new double[ld->molec]();
	pm->sum_xgt_mx_dt	= new double[ld->molec]();
	pm->sum_xgt_ca		= new double[ld->molec+ld->cati+ld->ani]();
	pm->sum_xgt_ca_dt	= new double[ld->molec+ld->cati+ld->ani]();

	pm->lngm			= new double[pd->nstoff]();
	pm->lngm_dt			= new double[pd->nstoff]();
	pm->lngm_dx			= new double[pd->nstoff*pd->nstoff]();
	pm->lngm_dt_dx		= new double[pd->nstoff*pd->nstoff]();

	pm->last_t = 0.;

	return;
}

void kill_paramstruct(paramstruct *pm){	//loescht die Parameterstruktur
	delete [] pm->alpha;
	delete [] pm->alpha_dx;
	delete [] pm->alpha_ca_dx;
	delete [] pm->alpha_input;
	delete [] pm->alpha_dt_input;

	delete [] pm->gg;
	delete [] pm->gg_dx;
	delete [] pm->gg_dt;
	delete [] pm->gg_dt_dx;
	delete [] pm->gg_ca_dx;
	delete [] pm->gg_input; 
	delete [] pm->gg_dt_input; 

	delete [] pm->tau;
	delete [] pm->tau_dt;
	delete [] pm->tau_dx;
	delete [] pm->tau_dt_dx;
	delete [] pm->tau_ca_dx;
	delete [] pm->tau_input; 
	delete [] pm->tau_dt_input;

	delete [] pm->last_x_eff; 
	delete [] pm->last_x;
	delete [] pm->chrfrac;
	delete [] pm->chrfrac_dx;

	delete [] pm->sum_xg;
	delete [] pm->sum_xg_dt;
	delete [] pm->sum_xg_dx;
	delete [] pm->sum_xg_dt_dx;
	delete [] pm->sum_xg_mx;
	delete [] pm->sum_xg_mx_dt;
	delete [] pm->sum_xg_ca;
	delete [] pm->sum_xg_ca_dt;

	delete [] pm->sum_xgt;
	delete [] pm->sum_xgt_dt;
	delete [] pm->sum_xgt_dx;
	delete [] pm->sum_xgt_dt_dx;
	delete [] pm->sum_xgt_mx;
	delete [] pm->sum_xgt_mx_dt;
	delete [] pm->sum_xgt_ca;
	delete [] pm->sum_xgt_ca_dt;

	delete [] pm->lngm;
	delete [] pm->lngm_dt;
	delete [] pm->lngm_dx; 
	delete [] pm->lngm_dt_dx; 

	delete pm; 

	return;
}

double eval_binary_ms(tpropdata *pd, double concentration[], int species_one, int species_two)
{
	int i;
	double* D_ij_binary = new double [3]();
	int * is_anion  = new int[pd->nstoff]();
	int* is_cation  = new int [pd->nstoff]();
	int* is_neutral = new int [pd->nstoff]();

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
		return D_ca;
	}

	// in case we have two cations
	else if (is_cation[species_one] == 1 && is_cation[species_two] == 1  )
	{	// binary double cation cation diffusion coefficient
		double D_cc=1;
		// diffusion coefficients of cation-(common) anion interactions
		double* D_K1_A = new double [ani_counter]();
		double* D_K2_A = new double [ani_counter]();
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
		return -D_cc;
	}

	// in case we have two anions (has no theoretical basis according to Visser!!!) 
	else if (is_anion[species_one] == 1 && is_anion[species_two] == 1  )
	{	// binary double cation cation diffusion coefficient
		double D_aa=1;
		// diffusion coefficients of (common) cation- anion interactions
		double* D_A1_K = new double [cati_counter]();
		double* D_A2_K = new double [cati_counter]();
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
		double* D_K_w_A = new double [ani_counter]();
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
		double* D_A_w_K = new double [cati_counter]();
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
		return D_aw;
	}
	else 
	{
		printf(" \n the type of binary pair is not supported by the current implementation! ");
		printf(" \n e.g. only one neutral species is allowed ");
		throw 9999;
	}

}

void calc_binary_ms_diffkoeff(
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

		copy_5_3_array(tmp_p_ij, p_ij);
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
		copy_5_3_array(tmp_p_ij, p_ij);
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
		copy_5_3_array(tmp_p_ij, p_ij);
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
		copy_5_3_array(tmp_p_ij, p_ij);
	}
	else 	if( !strcmpu(pd->stname[species_one],"Na") && !strcmpu(pd->stname[species_two],"Cl")  
		|| !strcmpu(pd->stname[species_one],"Cl") && !strcmpu(pd->stname[species_two],"Na") ) 
	{
		double tmp_p_ij [5] [3]=
		{
			{0.0, -1.336E-09, -2.035E-09},
			{-8.018E-14, -3.064E-14, -2.236E-13},
			{-2.090E-16, -3.910E-15, -3.792E-15},
			{-7.028E-18, -3.765E-17, -3.777E-17},
			{-2.176E-12, -1.766E-12, -8.317E-12}
		};
		copy_5_3_array(tmp_p_ij, p_ij);
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
		copy_5_3_array(tmp_p_ij, p_ij);
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

void copy_5_3_array(double source[5][3], double target[5][3])
{
	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			target[i][j] = source[i][j];
		}
	}
}

