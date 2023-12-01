// props.h : Definition der Datenstruktur
// Autor: majo (majo02@avt.rwth-aachen.de), thbe (thbe02@avt.rwth-aachen.de)
// UML-Diagramm:
// -
// Quellen:
// -

#ifndef PROPS
#define PROPS

typedef struct tplibdata{ //benötigte Gemischparameter
	double *mm;	//Molmasse
	double *tc;	//kritische Temperatur
	double *pc;	//kritischer Druck
	double *vc; //spezifisches Volumen am kritischen Punkt
	double *zc; //komprimierbarkeit am kritischen Punkt
	double *rktzra; //Komprimierbarkeit nach dem Rackett-Model
	double *vcrkt; //kritisches Volumen nach dem Rackett-Model
	int *charge; //gibt die Ladung jeder Komponente an
	int molec; //gibt die Anzahl der Lösungsmittel an
	int cati; //gibt die Anzahl der Kationen an
	int ani; //gibt die Anzahl der Anionen an
	int *molec_v; //listet die Lösungsmittel auf
	int *ani_v; //listet die Anionen auf
	int *cati_v; //listet die Kationen auf 
	int nseg; //Anzahl der Segmente
	int	*henry; //gibt Henry-Komponenten an
	int *solid; //gibt feste Komponete an
	int enrtlref; 
	int defgmelc; 
	int mixmolvol; //Gemischvolumen
	int nreac; //Anzahl der Reaktionen
	int nsolid; //Anzahl fester Komponeten
	double *gmelcc; //electrolyte-molecule pair parameter
	double *gmelcd; //electrolyte-molecule pair parameter
	double *gmelce; //electrolyte-molecule pair parameter
	double *gmelcn; //electrolyte-molecule pair parameter
	double *cpdiec; //dielectric constant
	double *born; //born correction
	double **nseg_v; //amount of each kind of segment contained in a specific species
	double *nrtlta; 
	double *nrtltb; 
	double *nrtla;
	double *dG_0f; //freie Bildungeenthalpie
	double *dH_0f; //Bildungsenthalpie
	double *cc_cp_aq; 
	double *poly_cp_aq; 
	double *stoic_coeff; //stöchiometrische Koeffizienten
	double Tref; //Referenz Temperatur
} tplibdata;

typedef struct paramstruct{ //berechnete Gemischparameter
	double *tau; //binary interaction parameter
	double *tau_dt; //temperature derivative of binary interaction parameter
	double *gg; //binary interaction energy
	double *gg_dt; //temperature derivative of binary interaction energy
	double *alpha; //nonrandom factor parameter
	double *gg_ca_dx;
	double *tau_ca_dx;
	double *alpha_ca_dx; 
	double *tau_input;
	double *gg_input;
	double *alpha_input;
	double *tau_dt_input;
	double *gg_dt_input;
	double *alpha_dt_input; 
	double *last_x;	//letzte Gemischzusammensetzung
	double *last_x_eff;
	double last_t; //letzte Temperatur
	double *lngm; //Aktivitätskoeffizient
	double *lngm_dt; //temperaturliche Änderung des Aktivitätskoefizienten
	double *lngm_dx; //stoffliche Änderung des Aktivitätskoeffizienten
	double *lngm_dt_dx; //stofflich und temperaturliche Änderung des Aktivitätskoeffizienten
	double *chrfrac;
	double *chrfrac_dx;
	double *sum_xg;
	double *sum_xgt;
	double *sum_xg_mx;
	double *sum_xgt_mx; 
	double *sum_xg_mx_dt;
	double *sum_xgt_mx_dt;
	double *sum_xg_ca;
	double *sum_xgt_ca;
	double *sum_xg_ca_dt;
	double *sum_xgt_ca_dt;
	double *sum_xg_dt;
	double *sum_xgt_dt;
	double *sum_xg_dx;
	double *sum_xgt_dx;
	double *gg_dx; //amount of substance derivative of binary interaction energy
	double *tau_dx; //amount of substance derivative of binary interaction parameter
	double *alpha_dx; //amount of substance derivative of nonrandom factor parameter
	double *gg_dt_dx;
	double *tau_dt_dx;
	double *sum_xg_dt_dx;
	double *sum_xgt_dt_dx;
	int lngm_curr; //ist true, falls Aktivitätskoeffizient für aktuelle Systemkonfiguration berechnet wurde
	int lngm_dt_curr;
	int lngm_dx_curr;
} paramstruct; 

int strcmpu(
	const char *string1,
	const char *string2);

#endif //PROPS
