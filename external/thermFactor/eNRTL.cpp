// eNRTL1.cpp : Gemischberechnungen auf Grundlage von [1]
// Autor: Robert Pack, majo (majo02@avt.rwth-aachen.de), thbe (thbe02@avt.rwth-aachen.de)
// UML-Diagramm:
// -
// Quellen:
// [1] Pack, R., 2011; Development of a Computational Routine for Estimation
//     of Activity Coefficents in Electrolyte Systems, Studienarbeit, AVT.PT


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "eNRTL.h"
#include <time.h>


void get_lngm_x(tpropdata *pd, paramstruct *pm, double *lngm, double* lngm_dt, double* lngm_dx)
{
	int i, j, index, count;
	tplibdata *ld=(tplibdata*)pd->libdata;
	double *lngm_sr=NULL, *lngm_sr_dt=NULL, *lngm_sr_dx=NULL, *lngm_sr_ref=NULL; 
	double *lngm_sr_ref_dt=NULL, *lngm_sr_ref_dx=NULL, *lngm_lr=NULL, *lngm_lr_dt=NULL, *lngm_lr_dx=NULL;

	if(lngm){
		count			= ld->molec+ld->cati+ld->ani; 
		lngm_sr			= new double[count];
		lngm_sr_ref		= new double[count];
		lngm_lr			= new double[count];
		for(i=0; i<count; ++i){
			lngm_sr[i]		= 0.;
			lngm_sr_ref[i]	= 0.;
			lngm_lr[i]		= 0.;
		}
	}
	if(lngm_dt){
		count			= ld->molec+ld->cati+ld->ani; 
		lngm_sr_dt		= new double[count];
		lngm_sr_ref_dt	= new double[count];
		lngm_lr_dt		= new double[count]; 
		for(i=0; i<count; ++i){
			lngm_sr_dt[i]		= 0.;
			lngm_sr_ref_dt[i]	= 0.;
			lngm_lr_dt[i]		= 0.;
		}
	}
	if(lngm_dx){
		count			= (ld->molec+ld->cati+ld->ani)*(ld->molec+ld->cati+ld->ani); 
		lngm_sr_dx		= new double[count];
		lngm_sr_ref_dx	= new double[count];
		lngm_lr_dx		= new double[count];
		for(i=0; i<count; ++i){
			lngm_sr_dx[i]	= 0.;
			lngm_lr_dx[i]	= 0.;
		}
	}

	eval_lngm_sr_x(pd, pm, lngm_sr, lngm_sr_dt, lngm_sr_dx);
	eval_lngm_sr_x_ref(pd, pm, lngm_sr_ref, lngm_sr_ref_dt, lngm_sr_ref_dx);
	eval_lngm_lr_x(pd, pm, lngm_lr, lngm_lr_dt, lngm_lr_dx);

	if(ld->nseg == -1){
		for(i=0; i<pd->nstoff; ++i){
			if(lngm){
				if(ld->charge[i] == 0){
					lngm[i] = lngm_sr[i] - lngm_sr_ref[i] + lngm_lr[i];
				}else{
					lngm[i] = abs(ld->charge[i]) * (lngm_sr[i] - lngm_sr_ref[i]) + lngm_lr[i];
				}
			} 
			if(lngm_dt){
				if(ld->charge[i] == 0){
					lngm_dt[i] = lngm_sr_dt[i] - lngm_sr_ref_dt[i] + lngm_lr_dt[i];
				}else{
					lngm_dt[i] = abs(ld->charge[i]) * (lngm_sr_dt[i] - lngm_sr_ref_dt[i]) + lngm_lr_dt[i];
				}
			} 
			if(lngm_dx){
				for(j=0; j<pd->nstoff; ++j){
					index = i*pd->nstoff + j;
					if(ld->charge[i] == 0){
						lngm_dx[index] = lngm_sr_dx[index] - lngm_sr_ref_dx[index] + lngm_lr_dx[index];
					}else{
						lngm_dx[index] = abs(ld->charge[i]) * (lngm_sr_dx[index] - lngm_sr_ref_dx[index]) + lngm_lr_dx[index];
					}
				}
			} 
		}

		if(lngm){
			pm->lngm_curr = TRUE;
		}
		if(lngm_dt){
			pm->lngm_dt_curr = TRUE;
		}
		if(lngm_dx){
			pm->lngm_dx_curr = TRUE;
		}

	}else{
		if(ld->nseg > 0){
			double* x = pm->last_x;
			//get_lngm_fh(pd, lngm_x_fh, x);
			for (i=0; i<pd->nstoff; ++i){
				for (j=0; j<ld->molec+ld->ani+ld->cati; ++j){
					if(lngm){
						lngm[i] += ld->nseg_v[i][j] * (lngm_sr[j] - lngm_sr_ref[j] + lngm_lr[j]);
					}
					if(lngm_dt){
						lngm_dt[i] += ld->nseg_v[i][j] * (lngm_sr_dt[j] - lngm_sr_ref_dt[j] + lngm_lr_dt[j]);
					}
				}
			//lngm_x[i] += lngm_x_fh[i];
			}

			if(lngm){
				pm->lngm_curr = TRUE;
			}
			if(lngm_dt){
				pm->lngm_dt_curr = TRUE;
			}
			if(lngm_dx){
				pm->lngm_dx_curr = TRUE;
			}

			delete [] x;
		}
	}

	delete [] lngm_sr;
	delete [] lngm_sr_dt;
	delete [] lngm_sr_dx;
	delete [] lngm_sr_ref;
	delete [] lngm_sr_ref_dt;
	delete [] lngm_sr_ref_dx;
	delete [] lngm_lr;
	delete [] lngm_lr_dt;
	delete [] lngm_lr_dx;

	return;
}

void set_param_default(tpropdata *pd)
{
	int i, j, k, index, water = -1, species[3];
	tplibdata *ld= (tplibdata*) pd->libdata;
    
	for (i=0; i<ld->molec; i++) {
		if (!strcmpu(pd->stname[ld->molec_v[i]],"H2O") || !strcmpu(pd->stname[ld->molec_v[i]],"WATER") || !strcmpu(pd->stname[ld->molec_v[i]],"WASSER") ) {
			water = i;
		}
	}

    for (i=0; i<ld->molec; i++) {
        for (j=0; j<ld->ani; j++) {
            for (k=0; k<ld->cati; k++) {
                species[0] = ld->molec_v[i];
                species[2] = ld->ani_v[j];
                species[1] = ld->cati_v[k];
                get_index_input(ld, species, &index, 3);
                if (ld->gmelcc[index] == 0.0) {
                    if (i == water) {
                        ld->gmelcc[index] = DEFGMELCCMCAAQ;
                    }else{
                        ld->gmelcc[index] = DEFGMELCCMCANONAQ;
                    }
                    ld->gmelcd[index] = DEFGMELCDMCA;
                }
                if (ld->gmelcc[index+1] == 0.0) {
                    if (i == water) {
                        ld->gmelcc[index+1] = DEFGMELCCCAMAQ;
                    }else{
                        ld->gmelcc[index+1] = DEFGMELCCCAMNONAQ;
                    }
                    ld->gmelcd[index+1] = DEFGMELCDCAM;
                }
                if ((ld->gmelcn[index] == 0.0) & (ld->gmelcn[index+1] == 0.0)) ld->gmelcn[index] = ld->gmelcn[index+1] = DEFGMELCN;
                else if ((ld->gmelcn[index] != 0.0) & (ld->gmelcn[index+1] == 0.0)) ld->gmelcn[index+1] = ld->gmelcn[index];
                else if ((ld->gmelcn[index] == 0.0) & (ld->gmelcn[index+1] != 0.0)) ld->gmelcn[index] = ld->gmelcn[index+1];
            }
        }
    }
}

int update_pm(tpropdata *pd, paramstruct *pm, double t, double *x)
{
	int i, j, k, index, trig=0; 
	double test_recalc=0., TotalSeg=0., sum_xc_2=0., sum_xa_2=0.;
	long double sum_x_mole_a = 0., sum_x_mole_c = 0.;
	tplibdata* ld = (tplibdata*) pd->libdata; 

	for(i=0; i<ld->molec+ld->cati+ld->ani; ++i){
		test_recalc += abs(x[i] - pm->last_x[i]);
	} //berechnet die Prüfsumme test_recalc: es werden die aktuellen Stoffmengenanteile von den alten subtrahiert

	test_recalc += abs(t - pm->last_t) / 100;

	if(test_recalc > 0.){	//EPSMAXRECALC){
		for(i=0; i<pd->nstoff; ++i){
			pm->last_x[i]		= x[i];
		} //ersetzt die alten Stoffmengenanteile durch die neuen
		 
		double* x_eff	= pm->last_x_eff; 
		pm->last_t		= t; 

		for (i=0; i<ld->ani+ld->cati+ld->molec; i++) {
			x_eff[i] = 0.;
		} //
	
		if (ld->nseg == -1) { //für den Fall das nseg=-1
			for(i=0; i<pd->nstoff; ++i){
				for(j=0; j<pd->nstoff; ++j){
					if(i == j) ld->nseg_v[i][j] = 1.;
					else ld->nseg_v[i][j] = 0.;
				}
			} //erzeugt Diagonalmatrix in nseg_v

			for (i=0; i<pd->nstoff; i++) {
				if(ld->charge[i] == 0) x_eff[i] = x[i];
				if(ld->charge[i] < 0) x_eff[i] = (double) -ld->charge[i] * x[i];
				if(ld->charge[i] > 0) x_eff[i] = (double) ld->charge[i] * x[i];
			} //vereinfachte Variante falls nseg=-1, möglich durch die Diagonalmatrix
		}else { //für den Fall das nseg!=-1
			for (i=0; i<ld->ani+ld->cati+ld->molec; i++) {
				for (j=0; j<pd->nstoff; j++) {
					x_eff[i] += x[j] * ld->nseg_v[j][i];
					TotalSeg += x[j] * ld->nseg_v[j][i];
				}
			} //berechnet x_eff mit nseg_v
	
			for (i=0; i<ld->ani+ld->cati+ld->molec; i++) {
				x_eff[i] /= TotalSeg;
				if (ld->charge[i] > 0) x_eff[i] *= ld->charge[i];
				else if (ld->charge[i] < 0) x_eff[i] *= -ld->charge[i];
			} //setzt die effektiven Stoffmengenanteile entsprechend der Ladung der Komponente
		}
	
		double* chrfrac = pm->chrfrac; 
		for(i=0; i<ld->cati+ld->ani; ++i) chrfrac[i] = 0.; //setzt 

		for(i=0; i<ld->molec+ld->ani+ld->cati; i++){
			if(ld->charge[i]<0){
				sum_x_mole_a += x_eff[i];
			}else if (ld->charge[i]>0) {
				sum_x_mole_c += x_eff[i];
			}
		} //berechnet die Summe der Anionen sum_x_mol_a und der Kationen sum_x_mole_c

		for (i=0; i<ld->cati; i++) {
			if (sum_x_mole_c == 0.) chrfrac[i] = 1. / ld->cati;
			else chrfrac[i] = x_eff[ld->cati_v[i]] / sum_x_mole_c;
		} //berechnet chrfrac für die Kationen
	
		for (i=0; i<ld->ani; i++) {
			if (sum_x_mole_a == 0.) chrfrac[i+ld->cati] = 1. / ld->ani;
			else chrfrac[i+ld->cati] = x_eff[ld->ani_v[i]] / sum_x_mole_a;
		} //berechnet chrfrac für die Anionen

		double* chrfrac_dx	= pm->chrfrac_dx;
		for(i=0; i<ld->cati*ld->cati+ld->ani*ld->ani; ++i) chrfrac_dx[i] = 0.; //setzt die stoffliche Änderung von chrfrac auf 0

		for(i=0; i<ld->cati; ++i){
			sum_xc_2 += x_eff[ld->cati_v[i]];
		}	
		sum_xc_2 *= sum_xc_2; //berechnet das Quadrat der Summen der Stoffmengen der Kationen

		for(i=0; i<ld->ani; ++i){
			sum_xa_2 += x_eff[ld->ani_v[i]];
		}
		sum_xa_2 *= sum_xa_2; //berechnet das Quadrat der Summen der Stoffmengen der Anionen

		if(sum_xa_2 == 0 ) sum_xa_2 = 1.;
		if(sum_xc_2 == 0 ) sum_xc_2 = 1.;

		for(i=0; i<ld->cati; ++i){
			for(j=0; j<ld->cati; ++j){
				index = i*ld->cati+j;
				if(i == j){
					for(k=0; k<ld->cati; ++k){
						if(i != k) chrfrac_dx[index] += x_eff[ld->cati_v[k]];
					}
					chrfrac_dx[index] *= ld->charge[ld->cati_v[i]];
				}else{
					chrfrac_dx[index] -= x_eff[ld->cati_v[i]];
					chrfrac_dx[index] *= ld->charge[ld->cati_v[j]];
				}
				chrfrac_dx[index] /= sum_xc_2;
			}
		} 
		for(i=0; i<ld->ani; ++i){
			for(j=0; j<ld->ani; ++j){
				index = ld->cati*ld->cati+i*ld->ani+j;
				chrfrac_dx[index] = 0.; 
				if(i == j){
					for(k=0; k<ld->ani; ++k){
						if(i != k) chrfrac_dx[index] += x_eff[ld->ani_v[k]];
					}
					chrfrac_dx[index] *= -ld->charge[ld->ani_v[i]];
				}else{
					chrfrac_dx[index] -= x_eff[ld->ani_v[i]];
					chrfrac_dx[index] *= -ld->charge[ld->ani_v[j]];
				}
				chrfrac_dx[index] /= sum_xa_2;
			}
		}

		get_param(pd, pm); 
		get_param_dt(pd, pm);
		get_param_dx(pd, pm);

		eval_sums(pd, pm); 
		eval_sums_dt(pd, pm);
		eval_sums_dx(pd, pm);
		eval_sums_mx(pd, pm);
		eval_sums_ca(pd, pm);

		pm->lngm_curr = pm->lngm_dt_curr = pm->lngm_dx_curr = FALSE; //setzt die Aktivitätskoeffizienten auf nicht berechnet

		return TRUE;
	}

	return FALSE;
}

void eval_sums(tpropdata *pd, paramstruct *pm)
{
	int i, j, index, speciesA[2]; 
	tplibdata* ld = (tplibdata*) pd->libdata; 

	double* gg		= pm->gg; 
	double* tau		= pm->tau; 
	double* x_eff	= pm->last_x_eff; 
	double* sum_xg	= pm->sum_xg;
	double* sum_xgt	= pm->sum_xgt;

	for(i=0; i<ld->molec; ++i){
		speciesA[1]	= ld->molec_v[i]; 
		sum_xg[i]	= 0.;
		sum_xgt[i]	= 0.;
		for(j=0; j<ld->molec; ++j){
			speciesA[0]	= ld->molec_v[j];
			if(get_index_calc(ld, speciesA, &index, 2)){
				sum_xg[i]	+= x_eff[ld->molec_v[j]] * gg[index];
				sum_xgt[i]	+= x_eff[ld->molec_v[j]] * gg[index] * tau[index];
			}else{
				sum_xg[i]	+= x_eff[ld->molec_v[j]];
			}
		}
		for(j=0; j<ld->cati; ++j){
			speciesA[0]	= ld->cati_v[j]; 
			if(get_index_calc(ld, speciesA, &index, 2)){
				sum_xg[i]	+= x_eff[ld->cati_v[j]] * gg[index];
				sum_xgt[i]	+= x_eff[ld->cati_v[j]] * gg[index] * tau[index];
			}
		}
		for(j=0; j<ld->ani; ++j){
			speciesA[0]	= ld->ani_v[j]; 
			if(get_index_calc(ld, speciesA, &index, 2)){
				sum_xg[i]	+= x_eff[ld->ani_v[j]] * gg[index];
				sum_xgt[i]	+= x_eff[ld->ani_v[j]] * gg[index] * tau[index];
			}
		}
	}

	for(i=0; i<ld->cati; ++i){
		speciesA[1]	= ld->cati_v[i];
		sum_xg[ld->molec+i]		= 0.; 
		sum_xgt[ld->molec+i]	= 0.; 
		for(j=0; j<ld->molec; ++j){
			speciesA[0]	= ld->molec_v[j];
			if(get_index_calc(ld, speciesA, &index, 2)){
				sum_xg[ld->molec+i]		+= x_eff[ld->molec_v[j]] * gg[index];
				sum_xgt[ld->molec+i]	+= x_eff[ld->molec_v[j]] * gg[index] * tau[index];
			}
		}
		for(j=0; j<ld->ani; ++j){
			speciesA[0]	= ld->ani_v[j];
			if(get_index_calc(ld, speciesA, &index, 2)){
				sum_xg[ld->molec+i]		+= x_eff[ld->ani_v[j]] * gg[index];
				sum_xgt[ld->molec+i]	+= x_eff[ld->ani_v[j]] * gg[index] * tau[index];
			}
		}
	}

	for(i=0; i<ld->ani; ++i){
		speciesA[1]	= ld->ani_v[i];
		pm->sum_xg[ld->molec+ld->cati+i]	= 0.; 
		pm->sum_xgt[ld->molec+ld->cati+i]	= 0.; 
		for(j=0; j<ld->molec; ++j){
			speciesA[0]	= ld->molec_v[j];
			if(get_index_calc(ld, speciesA, &index, 2)){
				sum_xg[ld->molec+ld->cati+i]	+= x_eff[ld->molec_v[j]] * gg[index];
				sum_xgt[ld->molec+ld->cati+i]	+= x_eff[ld->molec_v[j]] * gg[index] * tau[index];
			}
		}
		for(j=0; j<ld->cati; ++j){
			speciesA[0]	= ld->cati_v[j]; 
			if(get_index_calc(ld, speciesA, &index, 2)){
				sum_xg[ld->molec+ld->cati+i]	+= x_eff[ld->cati_v[j]] * gg[index];
				sum_xgt[ld->molec+ld->cati+i]	+= x_eff[ld->cati_v[j]] * gg[index] * tau[index];
			}
		}
	}

	return; 
}
void eval_sums_dt(tpropdata *pd, paramstruct *pm)
{
	int i, j, index, speciesA[2]; 
	tplibdata* ld = (tplibdata*) pd->libdata; 

	double* gg		= pm->gg; 
	double* gg_dt	= pm->gg_dt;
	double* tau		= pm->tau;
	double* tau_dt	= pm->tau_dt;

	for(i=0; i<ld->molec; ++i){
		speciesA[1]	= ld->molec_v[i]; 
		pm->sum_xg_dt[i]	= 0.;
		pm->sum_xgt_dt[i]	= 0.;
		for(j=0; j<ld->molec; ++j){
			speciesA[0]	= ld->molec_v[j];
			if(get_index_calc(ld, speciesA, &index, 2)){
				pm->sum_xg_dt[i]	+= pm->last_x_eff[ld->molec_v[j]] * gg_dt[index];
				pm->sum_xgt_dt[i]	+= pm->last_x_eff[ld->molec_v[j]] 
					* (gg_dt[index] * tau[index] + gg[index] * tau_dt[index]);
			}
		}
		for(j=0; j<ld->cati; ++j){
			speciesA[0]	= ld->cati_v[j]; 
			if(get_index_calc(ld, speciesA, &index, 2)){
				pm->sum_xg_dt[i]	+= pm->last_x_eff[ld->cati_v[j]] * gg_dt[index];
				pm->sum_xgt_dt[i]	+= pm->last_x_eff[ld->cati_v[j]] 
					* (gg_dt[index] * tau[index] + gg[index] * tau_dt[index]);			
			}
		}
		for(j=0; j<ld->ani; ++j){
			speciesA[0]	= ld->ani_v[j]; 
			if(get_index_calc(ld, speciesA, &index, 2)){
				pm->sum_xg_dt[i]	+= pm->last_x_eff[ld->ani_v[j]] * gg_dt[index];
				pm->sum_xgt_dt[i]	+= pm->last_x_eff[ld->ani_v[j]] 
					* (gg_dt[index] * tau[index] + gg[index] * tau_dt[index]);
			}
		}
	}

	for(i=0; i<ld->cati; ++i){
		speciesA[1]	= ld->cati_v[i];
		pm->sum_xg_dt[ld->molec+i]	= 0.; 
		pm->sum_xgt_dt[ld->molec+i]	= 0.; 
		for(j=0; j<ld->molec; ++j){
			speciesA[0]	= ld->molec_v[j];
			if(get_index_calc(ld, speciesA, &index, 2)){
				pm->sum_xg_dt[ld->molec+i]	+= pm->last_x_eff[ld->molec_v[j]] * gg_dt[index];
				pm->sum_xgt_dt[ld->molec+i]	+= pm->last_x_eff[ld->molec_v[j]] 
					* (gg_dt[index] * tau[index] + gg[index] * tau_dt[index]);
			}
		}
		for(j=0; j<ld->ani; ++j){
			speciesA[0]	= ld->ani_v[j];
			if(get_index_calc(ld, speciesA, &index, 2)){
				pm->sum_xg_dt[ld->molec+i]	+= pm->last_x_eff[ld->ani_v[j]] * gg_dt[index];
				pm->sum_xgt_dt[ld->molec+i]	+= pm->last_x_eff[ld->ani_v[j]] 
					* (gg_dt[index] * tau[index] + gg[index] * tau_dt[index]);
			}
		}
	}

	for(i=0; i<ld->ani; ++i){
		speciesA[1]	= ld->ani_v[i];
		pm->sum_xg_dt[ld->molec+ld->cati+i]		= 0.; 
		pm->sum_xgt_dt[ld->molec+ld->cati+i]	= 0.; 
		for(j=0; j<ld->molec; ++j){
			speciesA[0]	= ld->molec_v[j];
			if(get_index_calc(ld, speciesA, &index, 2)){
				pm->sum_xg_dt[ld->molec+ld->cati+i]		+= pm->last_x_eff[ld->molec_v[j]] * gg_dt[index];
				pm->sum_xgt_dt[ld->molec+ld->cati+i]	+= pm->last_x_eff[ld->molec_v[j]] 
					* (gg_dt[index] * tau[index] + gg[index] * tau_dt[index]);
			}
		}
		for(j=0; j<ld->cati; ++j){
			speciesA[0]	= ld->cati_v[j]; 
			if(get_index_calc(ld, speciesA, &index, 2)){
				pm->sum_xg_dt[ld->molec+ld->cati+i]		+= pm->last_x_eff[ld->cati_v[j]] * gg_dt[index];
				pm->sum_xgt_dt[ld->molec+ld->cati+i]	+= pm->last_x_eff[ld->cati_v[j]] 
					* (gg_dt[index] * tau[index] + gg[index] * tau_dt[index]);
			}
		}
	}

	return; 
}
void eval_sums_dx(tpropdata *pd, paramstruct *pm)
{
	int i, j, k, speciesA[2], speciesA_dx[3];
	int index, index_calc, index_dx;
	tplibdata* ld = (tplibdata*) pd->libdata;

	double* tau			= pm->tau;
	double* tau_dx		= pm->tau_dx;
	double* gg			= pm->gg;
	double* gg_dx		= pm->gg_dx;
	double* alpha		= pm->alpha;
	double* alpha_dx	= pm->alpha_dx;

	double* x_eff		= pm->last_x_eff;

	double* sum_xg_dx	= pm->sum_xg_dx;
	double* sum_xgt_dx	= pm->sum_xgt_dx;

	for(i=0; i<ld->molec; ++i){
		speciesA[1]		= ld->molec_v[i]; 
		speciesA_dx[1]	= ld->molec_v[i]; 
		for(j=0; j<ld->molec; ++j){
			index = i*(ld->molec+ld->cati+ld->ani) + j; 
			speciesA[0]	= ld->molec_v[j]; 
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				sum_xg_dx[index]	= gg[index_calc]; 
				sum_xgt_dx[index]	= gg[index_calc] * tau[index_calc]; 
			}else{
				sum_xg_dx[index]	= 1.; 
				sum_xgt_dx[index]	= 0.; 
			}
		}
		for(j=0; j<ld->cati; ++j){
			index = i*(ld->molec+ld->cati+ld->ani) + ld->molec + j; 
			speciesA_dx[2]		= ld->cati_v[j]; 
			sum_xg_dx[index]	= 0.;
			sum_xgt_dx[index]	= 0.;
			for(k=0; k<ld->cati; ++k){
				if(j == k){
					speciesA[0]	= ld->cati_v[k];
					if(get_index_calc(ld, speciesA, &index_calc, 2)){
						sum_xg_dx[index]	+=  ld->charge[ld->cati_v[k]] * gg[index_calc]; 
						sum_xgt_dx[index]	+=  ld->charge[ld->cati_v[k]] * gg[index_calc] * tau[index_calc]; 
					}
				}			
			}
			for(k=0; k<ld->ani; ++k){
				speciesA[0]		= ld->ani_v[k];
				speciesA_dx[0]	= ld->ani_v[k]; 
				if(get_index_calc(ld, speciesA_dx, &index_dx, 3) && get_index_calc(ld, speciesA, &index_calc, 2)){
					sum_xg_dx[index]	+= x_eff[ld->ani_v[k]] * gg_dx[index_dx];
					sum_xgt_dx[index]	+= x_eff[ld->ani_v[k]] * (gg_dx[index_dx]
						* tau[index_calc] + gg[index_calc] * tau_dx[index_dx]); 
				}
			}
		}
		for(j=0; j<ld->ani; ++j){
			index = i*(ld->molec+ld->cati+ld->ani) + ld->molec + ld->cati + j;
			speciesA_dx[2]		= ld->ani_v[j]; 
			sum_xg_dx[index]	= 0.;
			sum_xgt_dx[index]	= 0.;
			for(k=0; k<ld->cati; ++k){
				speciesA[0]		= ld->cati_v[k];  
				speciesA_dx[0]	= ld->cati_v[k]; 
				if(get_index_calc(ld, speciesA_dx, &index_dx, 3) && get_index_calc(ld, speciesA, &index_calc, 2)){
					sum_xg_dx[index]	+= x_eff[ld->cati_v[k]] * gg_dx[index_dx];
					sum_xgt_dx[index]	+= x_eff[ld->cati_v[k]] * (gg_dx[index_dx]
						*tau[index_calc] + gg[index_calc] * tau_dx[index_dx]);
				}
			}
			for(k=0; k<ld->ani; ++k){
				if(j == k){
					speciesA[0]	= ld->ani_v[k]; 
					if(get_index_calc(ld, speciesA, &index_calc, 2)){ 
						sum_xg_dx[index]	+= abs(ld->charge[ld->ani_v[k]]) * gg[index_calc]; 
						sum_xgt_dx[index]	+= abs(ld->charge[ld->ani_v[k]]) * gg[index_calc] * tau[index_calc];
					}
				}
			}
		}
	}

	for(i=0; i<ld->cati; ++i){
		speciesA[1]		= ld->cati_v[i];
		speciesA_dx[1]	= ld->cati_v[i]; 

		for(j=0; j<ld->molec; ++j){
			index = ld->molec*(ld->molec+ld->cati+ld->ani) + i*(ld->molec+ld->ani) + j;
			speciesA[0]	= ld->molec_v[j];
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				sum_xg_dx[index]	= gg[index_calc];
				sum_xgt_dx[index]	= gg[index_calc] * tau[index_calc];
			}
		}

		for(j=0; j<ld->ani; ++j){
			speciesA_dx[2]	= ld->ani_v[j];
			index = ld->molec*(ld->molec+ld->cati+ld->ani) + i*(ld->molec+ld->ani) + ld->molec + j;
			sum_xg_dx[index]	= 0.;
			sum_xgt_dx[index]	= 0.;

			for(k=0; k<ld->molec; ++k){
				speciesA_dx[0]	= ld->molec_v[k];
				speciesA[0]		= ld->molec_v[k];
				if(get_index_calc(ld, speciesA_dx, &index_dx, 3) && get_index_calc(ld, speciesA, &index_calc, 2)){
					sum_xg_dx[index]	+= x_eff[ld->molec_v[k]] * gg_dx[index_dx];
					sum_xgt_dx[index]	+= x_eff[ld->molec_v[k]] * (gg_dx[index_dx] * tau[index_calc]
						+ gg[index_calc] * tau_dx[index_dx]);
 				}
			}

			for(k=0; k<ld->ani; ++k){
				speciesA[0]		= ld->ani_v[k];
				speciesA_dx[0]	= ld->ani_v[k];
				if(get_index_calc(ld, speciesA_dx, &index_dx, 3) && get_index_calc(ld, speciesA, &index_calc, 2)){
					if(j == k){
						sum_xg_dx[index]	+= -ld->charge[ld->ani_v[j]] * gg[index_calc];
						sum_xgt_dx[index]	+= -ld->charge[ld->ani_v[j]] * gg[index_calc] * tau[index_calc];
					}
					sum_xg_dx[index]	+= x_eff[ld->ani_v[k]] * gg_dx[index_dx];
					sum_xgt_dx[index]	+= x_eff[ld->ani_v[k]] * (gg_dx[index_dx] * tau[index_calc]
						+ gg[index_calc] * tau_dx[index_dx]);
				}
			}
		}
	}

	for(i=0; i<ld->ani; ++i){
		speciesA[1]		= ld->ani_v[i];
		speciesA_dx[1]	= ld->ani_v[i];
		for(j=0; j<ld->molec; ++j){
			index = ld->molec*(ld->molec+ld->cati+ld->ani) + ld->cati*(ld->molec+ld->ani) 
				+ i*(ld->molec+ld->cati) + j;
			speciesA[0]	= ld->molec_v[j];
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				sum_xg_dx[index]	= gg[index_calc];
				sum_xgt_dx[index]	= gg[index_calc] * tau[index_calc];
			}
		}

		for(j=0; j<ld->cati; ++j){
			index = ld->molec*(ld->molec+ld->cati+ld->ani) + ld->cati*(ld->molec+ld->ani) 
				+ i*(ld->molec+ld->cati) + ld->molec + j;
			speciesA_dx[2]	= ld->cati_v[j];
			sum_xg_dx[index]	= 0.;
			sum_xgt_dx[index]	= 0.;

			for(k=0; k<ld->molec; ++k){
				speciesA_dx[0]	= ld->molec_v[k];
				speciesA[0]		= ld->molec_v[k];
				if(get_index_calc(ld, speciesA_dx, &index_dx, 3) && get_index_calc(ld, speciesA, &index_calc, 2)){
					sum_xg_dx[index]	+= x_eff[ld->molec_v[k]] * gg_dx[index_dx];
					sum_xgt_dx[index]	+= x_eff[ld->molec_v[k]] * (gg_dx[index_dx] * tau[index_calc]
						+ gg[index_calc] * tau_dx[index_dx]);
 				}
			}

			for(k=0; k<ld->cati; ++k){
				speciesA[0]		= ld->cati_v[k];
				speciesA_dx[0]	= ld->cati_v[k];
				if(get_index_calc(ld, speciesA_dx, &index_dx, 3) && get_index_calc(ld, speciesA, &index_calc, 2)){
					if(j == k){
						sum_xg_dx[index]	+= ld->charge[ld->cati_v[j]] * gg[index_calc];
						sum_xgt_dx[index]	+= ld->charge[ld->cati_v[j]] * gg[index_calc] * tau[index_calc];
					}
					sum_xg_dx[index]	+= x_eff[ld->cati_v[k]] * gg_dx[index_dx];
					sum_xgt_dx[index]	+= x_eff[ld->cati_v[k]] * (gg_dx[index_dx] * tau[index_calc]
						+ gg[index_calc] * tau_dx[index_dx]);
				}
			}
		}
	}

	return;
}

void eval_sums_mx(tpropdata *pd, paramstruct *pm)
{
	int i, j, index, speciesA[2]; 
	double sum_solv=0.; 
	tplibdata* ld = (tplibdata*) pd->libdata; 

	double* x_solv		= new double[ld->molec+ld->cati+ld->ani];
	double* x_eff		= pm->last_x_eff; 

	double* sum_xg_mx	= pm->sum_xg_mx; 
	double* sum_xgt_mx	= pm->sum_xgt_mx; 

	double* gg			= pm->gg; 
	double* tau			= pm->tau; 

	for(i=0; i<ld->molec+ld->cati+ld->ani; ++i){
		if(!ld->henry[i] && testmolec(ld, i)){
			sum_solv += x_eff[i];
			x_solv[i] = x_eff[i]; 
		}else{
			x_solv[i]	= 0.; 
		}
	}

	for(i=0; i<ld->molec+ld->cati+ld->ani; ++i){
		x_solv[i] /= sum_solv; 
	}

	for(i=0; i<ld->molec; ++i){
		speciesA[1]	= ld->molec_v[i]; 
		sum_xg_mx[i]	= 0.;
		sum_xgt_mx[i]	= 0.; 
		for(j=0; j<ld->molec; ++j){
			speciesA[0]	= ld->molec_v[j];
			if(get_index_calc(ld, speciesA, &index, 2) && !ld->henry[ld->molec_v[j]]){
				sum_xg_mx[i]	+= x_solv[ld->molec_v[j]] * gg[index];
				sum_xgt_mx[i]	+= x_solv[ld->molec_v[j]] * gg[index] * tau[index];
			}else if(!ld->henry[j]){
				sum_xg_mx[i]	+= x_solv[ld->molec_v[j]];
			}
		}
	}

	delete [] x_solv; 

	return; 
}

void eval_sums_ca(tpropdata *pd, paramstruct *pm)
{
	int i, j, index, speciesA[2]; 
	tplibdata* ld = (tplibdata*) pd->libdata; 
	
	double* gg			= pm->gg; 
	double* tau			= pm->tau; 
	double* sum_xg_ca	= pm->sum_xg_ca; 
	double* sum_xgt_ca	= pm->sum_xgt_ca; 
	double* x_ca		= new double[ld->molec+ld->cati+ld->ani]; 
	
	double sum_charge=0.; 

	for(i=0; i<ld->molec+ld->cati+ld->ani; ++i){
		sum_charge += (double) abs(ld->charge[i]); 
	}
	
	for(i=0; i<ld->molec+ld->cati+ld->ani; ++i){
		x_ca[i] = (double) abs(ld->charge[i]) / sum_charge; 
	}

	for(i=0; i<ld->molec; ++i){
		speciesA[1]	= ld->molec_v[i]; 

		for(j=0; j<ld->cati; ++j){
			speciesA[0]	= ld->cati_v[j]; 
			if(get_index_calc(ld, speciesA, &index, 2)){
				sum_xg_ca[i]	+= x_ca[ld->cati_v[j]] * gg[index];
				sum_xgt_ca[i]	+= x_ca[ld->cati_v[j]] * gg[index] * tau[index];
			}
		}
		for(j=0; j<ld->ani; ++j){
			speciesA[0]	= ld->ani_v[j]; 
			if(get_index_calc(ld, speciesA, &index, 2)){
				sum_xg_ca[i]	+= x_ca[ld->ani_v[j]] * gg[index];
				sum_xgt_ca[i]	+= x_ca[ld->ani_v[j]] * gg[index] * tau[index];
			}
		}
	}

	for(i=0; i<ld->cati; ++i){
		speciesA[1]	= ld->cati_v[i];
		for(j=0; j<ld->ani; ++j){
			speciesA[0]	= ld->ani_v[j];
			if(get_index_calc(ld, speciesA, &index, 2)){
				sum_xg_ca[ld->molec+i]	+= x_ca[ld->ani_v[j]] * gg[index];
				sum_xgt_ca[ld->molec+i]	+= x_ca[ld->ani_v[j]] * gg[index] * tau[index];
			}
		}
	}

	for(i=0; i<ld->ani; ++i){
		speciesA[1]	= ld->ani_v[i];
		for(j=0; j<ld->cati; ++j){
			speciesA[0]	= ld->cati_v[j]; 
			if(get_index_calc(ld, speciesA, &index, 2)){
				sum_xg_ca[ld->molec+ld->cati+i]		+= x_ca[ld->cati_v[j]] * gg[index];
				sum_xgt_ca[ld->molec+ld->cati+i]	+= x_ca[ld->cati_v[j]] * gg[index] * tau[index];
			}
		}
	}

	delete [] x_ca; 

	return; 
}

int get_index_calc(tplibdata *ld, int *species, int *index, int choose)
{	
	int i, p1, p2, der;
	
	//size_t choose1 = sizeof(species) / sizeof(species[0]);

	switch(choose){
	
		case 2:
			if((ld->charge[species[0]] == 0) && (ld->charge[species[1]] == 0) && (species[0]!=species[1])){
				for(i=0; i<ld->molec; ++i){
					if(ld->molec_v[i] == species[0]){ p1 = i; }
					if(ld->molec_v[i] == species[1]){ p2 = i; }
				}
				if(p1 < p2) *index = p1 * (ld->molec - 1) + p2 - 1;
				else if (p2 < p1) *index = p1 * (ld->molec - 1) + p2;
				return TRUE;
			}else if((ld->charge[species[0]] == 0) && (ld->charge[species[1]] > 0)){
				for(i=0; i<ld->molec; ++i)	if(ld->molec_v[i] == species[0]){ p1 = i; break; }
				for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == species[1])	{ p2 = i; break; }	
				*index = ld->molec*(ld->molec - 1) + p2*ld->molec + p1;
				return TRUE;
			}else if((ld->charge[species[0]] == 0) && (ld->charge[species[1]] < 0)){
				for(i=0; i<ld->molec; ++i)	if(ld->molec_v[i] == species[0]){ p1 = i; break; }
				for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == species[1])	{ p2 = i; break; }
				*index = ld->molec*(ld->molec - 1) + ld->cati*ld->molec + p2*ld->molec + p1;
				return TRUE;
			}else if((ld->charge[species[0]] > 0) && (ld->charge[species[1]] == 0)){
				for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == species[0]) { p1 = i; break; }
				for(i=0; i<ld->molec; ++i)	if(ld->molec_v[i] == species[1]){ p2 = i; break; }
				*index = ld->molec*((ld->molec -1)+ld->cati+ld->ani) + p1*ld->molec + p2;
				return TRUE;
			}else if((ld->charge[species[0]] < 0) && (ld->charge[species[1]] == 0)){
				for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == species[0])	{ p1 = i; break; }
				for(i=0; i<ld->molec; ++i)	if(ld->molec_v[i] == species[1]){ p2 = i; break; }
				*index = ld->molec*((ld->molec -1)+ld->cati+ld->ani)
					+ ld->cati*ld->molec + p1*ld->molec + p2;
				return TRUE;
			}else if((ld->charge[species[0]] > 0) && (ld->charge[species[1]] < 0)){
				for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == species[0]) { p1 = i; break; }
				for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == species[1])	{ p2 = i; break; }
				*index = ld->molec*((ld->molec -1)+2*(ld->cati+ld->ani)) 
					+ p2*ld->cati + p1; 
				return TRUE;
			}else if((ld->charge[species[0]] < 0) && (ld->charge[species[1]] > 0)){
				for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == species[0])	{ p1 = i; break; }
				for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == species[1])	{ p2 = i; break; }
				*index = ld->molec*((ld->molec -1)+2*(ld->cati+ld->ani)) 
					+ ld->ani*ld->cati + p2*ld->ani + p1;
				return TRUE; 
			}
			break;

		case 3:
			if((ld->charge[species[0]] == 0) && (ld->charge[species[1]] > 0) && (ld->charge[species[2]] < 0)){
				for(i=0; i<ld->molec; ++i)	if(ld->molec_v[i] == species[0]){ p1 = i; break; }
				for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == species[1])	{ p2 = i; break; }
				for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == species[2])	{ der= i; break; }
				*index = (p2*ld->molec + p1) * ld->ani + der;
				return TRUE;
			}else if((ld->charge[species[0]] == 0) && (ld->charge[species[1]] < 0) && (ld->charge[species[2]] > 0)){
				for(i=0; i<ld->molec; ++i)	if(ld->molec_v[i] == species[0]){ p1 = i; break; }
				for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == species[1])	{ p2 = i; break; }
				for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == species[2])	{ der= i; break; }
				*index = ld->cati*ld->molec*ld->ani + (p2*ld->molec + p1) * ld->cati + der;
				return TRUE;
			}else if((ld->charge[species[0]] > 0) && (ld->charge[species[1]] == 0) && (ld->charge[species[2]] < 0)){
				for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == species[0]) { p1 = i; break; }
				for(i=0; i<ld->molec; ++i)	if(ld->molec_v[i] == species[1]){ p2 = i; break; }
				for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == species[2])	{ der= i; break; }
				*index = 2*ld->molec*ld->cati*ld->ani + (p1*ld->molec + p2) * ld->ani + der;
				return TRUE;
			}else if((ld->charge[species[0]] < 0) && (ld->charge[species[1]] == 0) && (ld->charge[species[2]] > 0)){
				for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == species[0])	{ p1 = i; break; }
				for(i=0; i<ld->molec; ++i)	if(ld->molec_v[i] == species[1]){ p2 = i; break; }
				for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == species[2])	{ der= i; break; }
				*index = 3*ld->molec*ld->cati*ld->ani + (p1*ld->molec + p2) * ld->cati + der;
				return TRUE;
			}else if((ld->charge[species[0]] > 0) && (ld->charge[species[1]] < 0) && (ld->charge[species[2]] > 0)){
				for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == species[0]) { p1 = i; break; }
				for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == species[1])	{ p2 = i; break; }
				for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == species[2])	{ der= i; break; }
				*index = 4*ld->molec*ld->cati*ld->ani + (p2*ld->cati + p1) * ld->cati + der; 
				return TRUE;
			}else if((ld->charge[species[0]] < 0) && (ld->charge[species[1]] > 0) && (ld->charge[species[2]] < 0)){
				for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == species[0])	{ p1 = i; break; }
				for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == species[1])	{ p2 = i; break; }
				for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == species[2])	{ der= i; break; }
				*index = 4*ld->molec*ld->cati*ld->ani + ld->cati*ld->ani*ld->cati 
					+ (p2*ld->ani + p1) * ld->ani + der;
				return TRUE; 
			}
			break;
	}

	return FALSE;
}

void get_index_input(tplibdata *ld, int *species, int *index, int size){
	
	int trig=0, i, j, k, center, molec, p1, p2;
	
	//added by thbe
	//int ani=0, cati=0;
	//old version: 
	int ani,cati;
	//end adding

	switch(size){	
		case 2:				/* indices for molecular - molecular parameter  */
			for(i=0; i<ld->molec; i++){
				if(ld->molec_v[i] == species[0]){
					p1 = i;
				}
				if(ld->molec_v[i] == species[1]){
					p2 = i;
				}
			}
			if(p1 < p2){ 
				*index = p1 * (ld->molec - 1) + p2 - 1;
			}else{
				if(p2 < p1){
					*index = p1 * (ld->molec - 1) + p2;
				}
			}
			break;


		case 3:				/* indices for molecular - electrolyte parameters */	
			for(i=0; i<3; i++){
				//printf("(%f)-----\n",(float) i);
				//printf("charge: %f\n",(float) ld->charge[species[i]]);
				if(ld->charge[species[i]]==0){
					for(j=0; j<ld->molec; j++){
						if(ld->molec_v[j]==species[i]){
							molec = j;
						}
						if(i==2){
							trig=1;
						}
					}
				}else{
					j = ld->charge[species[i]]/abs(ld->charge[species[i]]);
					//printf("j: %f\n",(float) j);
					switch(j){
						case 1:
							for(j=0; j<ld->cati; j++){
								if(ld->cati_v[j]==species[i]){
									cati = j;
								}
							}
							break;
							
						case -1:
							for(j=0; j<ld->ani; j++){
								//printf("ld->ani_v: %f\n",(float) ld->ani_v[j]);
								//printf("species: %f\n",(float) species[i]);
								if(ld->ani_v[j]==species[i]){
									ani = j;
									//printf("ani/j: %f / %f\n",(float) ani, (float) j);
								}
							}
							break;
							
						default:
							break;
					}
				}
			}
			
			//printf("ani: %f\n\n", (float) ani);
			*index = 2*(molec*ld->cati*ld->ani + ani*ld->cati + cati) + trig;
			break;
			
		case 4:				/* indices for electrolyte - electrolyte parameters */	
			for (i=0; i<2; i++){
				if(species[0]==species[i+2]){
					center = 0;
					break;
				}
				if(species[1]==species[i+2]){
					center = 1;
					break;
				}
			}
			
			j = ld->charge[species[center]]/abs(ld->charge[species[center]]);
			
			switch(j){
					
				case 1:		/* cationic center */
					p1 = 1 - center;
					p2 = 3 - i;
					for(k=0; k<ld->cati; k++){
						if(species[center] == ld->cati_v[k]){
							center = k;
						}
					}
					for(k=0; k<ld->ani; k++) {
						if(species[p1] == ld->ani_v[k]){
							ani = k;
						}
						if(species[p2] == ld->ani_v[k]){
							cati = k;
						}
					}
					
					if(cati > ani){
						*index = 2*ld->molec*ld->ani*ld->cati + ld->ani*ld->cati*(ld->cati-1) + ld->ani*(ld->ani-1)*center + cati*(ld->ani-1) + ani;
					}else{
						if (ani > cati){
							*index = 2*ld->molec*ld->ani*ld->cati + ld->ani*ld->cati*(ld->cati-1) + ld->ani*(ld->ani-1)*center + cati*(ld->ani-1) + ani - 1;
						}
					}
					break;
					
				case -1:	/* anionic center */
					p1 = 1 - center;
					p2 = 3 - i;
					for(k=0; k<ld->ani; k++){
						if(species[center] == ld->ani_v[k]){
							center = k;
						}
					}
					for(k=0; k<ld->cati; k++){
						if(species[p1] == ld->cati_v[k]){
							ani = k;
						}
						if(species[p2] == ld->cati_v[k]){
							cati = k;
						}
					}
					
					if(cati > ani){
						*index = 2*ld->molec*ld->ani*ld->cati + ld->cati*(ld->cati-1)*center + cati*(ld->cati-1) + ani;
					}else{
						if (ani > cati){
							*index = 2*ld->molec*ld->ani*ld->cati + ld->cati*(ld->cati-1)*center + cati*(ld->cati-1) + ani - 1;
					
						}
					}
					break;
					
				default:
					break;
			}
			break;
			
		default:
			break;
	}
}

void get_param(tpropdata* pd, paramstruct* pm)
{
	int i, j, k, speciesA[2], speciesB[3], speciesC[4], index_input, index_calc; 
	tplibdata* ld	= (tplibdata*) pd->libdata; 
	double tau0=0., tau1=0.;

	const double t		= pm->last_t;
	const double value	= ((ENRTLTREF-t)/t + log(t/ENRTLTREF));

	const int param_count	= ld->molec*((ld->molec -1)+2*(ld->cati+ld->ani)) + 2*(ld->cati+ld->ani); 
	const int param_input	= ld->molec*(2*(ld->cati*ld->ani)) + ld->cati*ld->ani*(ld->ani-1)
		+ ld->ani*ld->cati*(ld->cati-1); 

	double* gg			= pm->gg;
	double* tau			= pm->tau; 
	double* alpha		= pm->alpha; 
	double* gg_input	= pm->gg_input;
	double* tau_input	= pm->tau_input; 
	double* alpha_input	= pm->alpha_input;
	double* chrfrac		= pm->chrfrac; 
	
	for(i=0; i<param_count; ++i){
		gg[i]		= 0.; 
		tau[i]		= 0.; 
		alpha[i]	= 0.; 
	}
	
	for(i=0; i<ld->molec; ++i){
		speciesA[0] = ld->molec_v[i];	
		for(j=0; j<ld->molec; ++j){
			speciesA[1] = ld->molec_v[j];
			if (get_index_calc(ld, speciesA, &index_calc, 2)) {
				tau[index_calc]			= ld->nrtlta[index_calc] + ld->nrtltb[index_calc] / t;
				alpha[index_calc]		= ld->nrtla[index_calc];
				gg[index_calc]			= exp(-alpha[index_calc]*tau[index_calc]);
			}
		}
	
		speciesB[0]	= ld->molec_v[i]; 
		for(j=0; j<ld->cati; ++j){
			speciesB[1]	= ld->cati_v[j]; 
			for(k=0; k<ld->ani; ++k){
				speciesB[2]	= ld->ani_v[k]; 
				get_index_input(ld, speciesB, &index_input, 3); 

				tau_input[index_input]	= ld->gmelcc[index_input] + ld->gmelcd[index_input]/t 
					+ ld->gmelce[index_input] * value;
				alpha_input[index_input] = ld->gmelcn[index_input]; 
				gg_input[index_input]	= exp(-alpha_input[index_input]*tau_input[index_input]); 

				tau_input[index_input+1]	= ld->gmelcc[index_input+1] + ld->gmelcd[index_input+1]/t 
					+ ld->gmelce[index_input+1] * value;
				alpha_input[index_input+1] = ld->gmelcn[index_input+1]; 
				gg_input[index_input+1]	= exp(-alpha_input[index_input+1]*tau_input[index_input+1]);

				speciesA[0] = ld->molec_v[i];
				speciesA[1]	= ld->cati_v[j]; 
				if(get_index_calc(ld, speciesA, &index_calc, 2)){
					gg[index_calc]		+= chrfrac[ld->cati+k] * gg_input[index_input];	
					alpha[index_calc]	+= chrfrac[ld->cati+k] * alpha_input[index_input];
				}
				
				speciesA[0]	= ld->cati_v[j]; 
				speciesA[1] = ld->molec_v[i]; 
				if(get_index_calc(ld, speciesA, &index_calc, 2)){
					gg[index_calc]		+= chrfrac[ld->cati+k] * gg_input[index_input+1];	
					alpha[index_calc]	+= chrfrac[ld->cati+k] * alpha_input[index_input+1]; 
				}			
			
				speciesA[0]	= ld->molec_v[i]; 
				speciesA[1]	= ld->ani_v[k]; 
				if(get_index_calc(ld, speciesA, &index_calc, 2)){
					gg[index_calc]		+= chrfrac[j] * gg_input[index_input];
					alpha[index_calc]	+= chrfrac[j] * alpha_input[index_input];
				}
				
				speciesA[0] = ld->ani_v[k]; 
				speciesA[1] = ld->molec_v[i]; 
				if(get_index_calc(ld, speciesA, &index_calc, 2)){
					gg[index_calc]		+= chrfrac[j] * gg_input[index_input+1];
					alpha[index_calc]	+= chrfrac[j] * alpha_input[index_input+1]; 
				}
			}
		}
		
		for(j=0; j<ld->cati; ++j){
			speciesA[0] = ld->molec_v[i]; 
			speciesA[1] = ld->cati_v[j]; 
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				tau[index_calc] = - log(gg[index_calc]) / alpha[index_calc]; 
			}

			speciesA[0]	= ld->cati_v[j]; 
			speciesA[1] = ld->molec_v[i];
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				tau[index_calc] = - log(gg[index_calc]) / alpha[index_calc]; 
			}
		}

		for(j=0; j<ld->ani; ++j){
			speciesA[0] = ld->molec_v[i]; 
			speciesA[1] = ld->ani_v[j]; 
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				tau[index_calc] = - log(gg[index_calc]) / alpha[index_calc]; 
			}

			speciesA[0] = ld->ani_v[j]; 
			speciesA[1] = ld->molec_v[i]; 
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				tau[index_calc] = - log(gg[index_calc]) / alpha[index_calc]; 
			}
 		}
	}

	for(i=0; i<ld->cati; ++i){
		for(j=0; j<ld->ani; ++j){
			speciesC[1] = speciesC[3] = ld->ani_v[j]; 
			speciesC[0] = ld->cati_v[i];

			speciesA[0] = ld->cati_v[i]; 
			speciesA[1] = ld->ani_v[j]; 
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				for(k=0; k<ld->cati; ++k){
					if(i != k){
						speciesC[2] = ld->cati_v[k];
						get_index_input(ld, speciesC, &index_input, 4); 
						tau_input[index_input] = ld->gmelcc[index_input] + ld->gmelcd[index_input]/t 
							+ ld->gmelce[index_input] * value;
						alpha_input[index_input] =  ld->gmelcn[index_input];
						gg_input[index_input] = exp(-alpha_input[index_input]*tau_input[index_input]);

						gg[index_calc]		+= chrfrac[k] * gg_input[index_input]; 
						alpha[index_calc]	+= chrfrac[k] * alpha_input[index_input]; 
					}else{
						gg[index_calc]		+= chrfrac[k]; 
						alpha[index_calc]	+= chrfrac[k] * 0.2; 
					}	 	
				}
				tau[index_calc] = - log(gg[index_calc]) / alpha[index_calc]; 
			}
			
			speciesC[1] = speciesC[3] = ld->cati_v[i];
			speciesC[0] = ld->ani_v[j]; 

			speciesA[0] = ld->ani_v[j]; 
			speciesA[1] = ld->cati_v[i]; 
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				for(k=0; k<ld->ani; ++k){
					if(j != k){
						speciesC[2] = ld->ani_v[k]; 
						get_index_input(ld, speciesC, &index_input, 4); 
						tau_input[index_input] = ld->gmelcc[index_input] + ld->gmelcd[index_input]/t 
							+ ld->gmelce[index_input] * value;
						alpha_input[index_input] =  ld->gmelcn[index_input];
						gg_input[index_input] = exp(-alpha_input[index_input]*tau_input[index_input]);

						gg[index_calc]		+= chrfrac[ld->cati+k] * gg_input[index_input];  
						alpha[index_calc]	+= chrfrac[ld->cati+k] * alpha_input[index_input]; 
					}else{
						gg[index_calc]		+= chrfrac[ld->cati+k]; 
						alpha[index_calc]	+= chrfrac[ld->cati+k] * 0.2; 
					}
				}
				tau[index_calc] = - log(gg[index_calc]) / alpha[index_calc];
			}
		}
	}

	return;
}
void get_param_dt(tpropdata* pd, paramstruct* pm)
{
	int i, j, k, speciesA[2], speciesB[3], speciesC[4], index_input, index_calc; 
	tplibdata* ld	= (tplibdata*) pd->libdata; 
	double tau0=0., tau1=0.;
	double t = pm->last_t;

	int param_count	= ld->molec*((ld->molec -1)+2*(ld->cati+ld->ani)) + 2*(ld->cati+ld->ani); 
	
	double* gg_dt		= pm->gg_dt;
	double* gg_dt_input	= pm->gg_dt_input;
	double* gg			= pm->gg;
	double* gg_input	= pm->gg_input;
	double* tau_dt		= pm->tau_dt;
	double* tau_dt_input= pm->tau_dt_input;
	double* tau			= pm->tau;
	double* alpha		= pm->alpha;	
	double* alpha_input	= pm->alpha_input;
	double* chrfrac		= pm->chrfrac; 
	
	for(i=0; i<param_count; ++i){
		gg_dt[i]		= 0.; 
		tau_dt[i]		= 0.;  
	}
	
	for(i=0; i<ld->molec; ++i){
		speciesA[0] = ld->molec_v[i];	
		for(j=0; j<ld->molec; ++j){
			speciesA[1] = ld->molec_v[j];
			if (get_index_calc(ld, speciesA, &index_calc, 2)) {
				tau_dt[index_calc]	= - ld->nrtltb[index_calc] / pow(t, 2);
				gg_dt[index_calc]	= - alpha[index_calc] * gg[index_calc] * tau_dt[index_calc];
			}
		}
	
		speciesB[0]	= ld->molec_v[i]; 
		for(j=0; j<ld->cati; ++j){
			speciesB[1]	= ld->cati_v[j]; 
			for(k=0; k<ld->ani; ++k){
				speciesB[2]	= ld->ani_v[k]; 
				get_index_input(ld, speciesB, &index_input, 3); 
				tau_dt_input[index_input]	= ld->gmelce[index_input]/t - (ld->gmelce[index_input]
						* ENRTLTREF + ld->gmelcd[index_input]) / pow(t, 2);

				tau_dt_input[index_input+1]	= ld->gmelce[index_input+1]/t - (ld->gmelce[index_input+1]
						* ENRTLTREF + ld->gmelcd[index_input+1]) / pow(t, 2);

				speciesA[0] = ld->molec_v[i];
				speciesA[1]	= ld->cati_v[j]; 
				if(get_index_calc(ld, speciesA, &index_calc, 2)){
					gg_dt[index_calc]	-= chrfrac[ld->cati+k] * alpha_input[index_input] 
						* gg_input[index_input] * tau_dt_input[index_input]; 
				}

				speciesA[0]	= ld->cati_v[j]; 
				speciesA[1] = ld->molec_v[i]; 
				if(get_index_calc(ld, speciesA, &index_calc, 2)){
					gg_dt[index_calc]	-= chrfrac[ld->cati+k] * alpha_input[index_input+1] 
						* gg_input[index_input+1] * tau_dt_input[index_input+1]; 
				}			
			
				speciesA[0]	= ld->molec_v[i]; 
				speciesA[1]	= ld->ani_v[k]; 
				if(get_index_calc(ld, speciesA, &index_calc, 2)){
					gg_dt[index_calc]	-=  chrfrac[j] * alpha_input[index_input] 
						* gg_input[index_input] * tau_dt_input[index_input];
				}
				
				speciesA[0] = ld->ani_v[k]; 
				speciesA[1] = ld->molec_v[i]; 
				if(get_index_calc(ld, speciesA, &index_calc, 2)){
					gg_dt[index_calc]	-=  chrfrac[j] * alpha_input[index_input+1] 
						* gg_input[index_input+1] * tau_dt_input[index_input+1];
				}
			}
		}
		
		for(j=0; j<ld->cati; ++j){
			speciesA[0] = ld->molec_v[i]; 
			speciesA[1] = ld->cati_v[j]; 
			if(get_index_calc(ld, speciesA, &index_calc,2)){
				tau_dt[index_calc] = - gg_dt[index_calc] / (alpha[index_calc] * gg[index_calc]); 
			}

			speciesA[0]	= ld->cati_v[j]; 
			speciesA[1] = ld->molec_v[i];
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				tau_dt[index_calc] = - gg_dt[index_calc] / (alpha[index_calc] * gg[index_calc]);
			}
		}

		for(j=0; j<ld->ani; ++j){
			speciesA[0] = ld->molec_v[i]; 
			speciesA[1] = ld->ani_v[j]; 
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				tau_dt[index_calc] = - gg_dt[index_calc] / (alpha[index_calc] * gg[index_calc]); 
			}

			speciesA[0] = ld->ani_v[j]; 
			speciesA[1] = ld->molec_v[i]; 
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				tau_dt[index_calc] = - gg_dt[index_calc] / (alpha[index_calc] * gg[index_calc]); 
			}
 		}
	}

	for(i=0; i<ld->cati; ++i){
		for(j=0; j<ld->ani; ++j){
			speciesC[1] = speciesC[3] = ld->ani_v[j]; 
			speciesC[0] = ld->cati_v[i];

			speciesA[0] = ld->cati_v[i]; 
			speciesA[1] = ld->ani_v[j]; 
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				for(k=0; k<ld->cati; ++k){
					if(i != k){
						speciesC[2] = ld->cati_v[k];
						get_index_input(ld, speciesC, &index_input, 4); 
						tau_dt_input[index_input]	= ld->gmelce[index_input]/t - (ld->gmelce[index_input]
							* ENRTLTREF + ld->gmelcd[index_input]) / pow(t, 2);
						gg_dt[index_calc]	-= chrfrac[k] * alpha_input[index_input] 
							* gg_input[index_input] * tau_dt_input[index_input]; 
					}	 	
				}
				tau_dt[index_calc] = - gg_dt[index_calc] / (alpha[index_calc] * gg[index_calc]);
			}
			
			speciesC[1] = speciesC[3] = ld->cati_v[i];
			speciesC[0] = ld->ani_v[j]; 

			speciesA[0] = ld->ani_v[j]; 
			speciesA[1] = ld->cati_v[i]; 
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				for(k=0; k<ld->ani; ++k){
					if(j != k){
						speciesC[2] = ld->ani_v[k]; 
						get_index_input(ld, speciesC, &index_input, 4); 
						tau_dt_input[index_input]	= ld->gmelce[index_input]/t - (ld->gmelce[index_input]
							* ENRTLTREF + ld->gmelcd[index_input]) / pow(t, 2);
						gg_dt[index_calc]	-= chrfrac[ld->cati+k] * alpha_input[index_input] 
							* gg_input[index_input] * tau_dt_input[index_input]; 
					}
				}
				tau_dt[index_calc] = - gg_dt[index_calc] / (alpha[index_calc] * gg[index_calc]);
			}
		}
	}

	return;
}
void get_param_dx(tpropdata* pd, paramstruct* pm)
{
	int i, j, k, l, speciesA[2], speciesA_dx[3], speciesB[3], speciesC[4], index_input, index_calc, index_chr, index_input_1; 
	tplibdata* ld	= (tplibdata*) pd->libdata; 
	double tau0=0., tau1=0.;
	double t = pm->last_t;

	int param_count	= 4*ld->molec*ld->cati*ld->ani + ld->cati*ld->ani*ld->cati 
		+ ld->ani*ld->cati*ld->ani; 
	
	double* gg			= pm->gg;
	double* gg_dx		= pm->gg_dx;
	double* gg_input	= pm->gg_input;
	double* tau			= pm->tau;
	double* tau_dx		= pm->tau_dx;
	double* alpha		= pm->alpha;
	double* alpha_dx	= pm->alpha_dx;
	double* alpha_input	= pm->alpha_input;
	double* chrfrac		= pm->chrfrac;
	double* chrfrac_dx	= pm->chrfrac_dx;
	
	for(i=0; i<param_count; ++i){
		gg_dx[i]	= 0.; 
		tau_dx[i]	= 0.; 
		alpha_dx[i]	= 0.; 
	}

	for(i=0; i<ld->molec; ++i){
		speciesB[0]	= ld->molec_v[i]; 
		for(j=0; j<ld->cati; ++j){
			speciesB[1]	= ld->cati_v[j]; 
			for(k=0; k<ld->ani; ++k){
				speciesB[2]	= ld->ani_v[k]; 
				get_index_input(ld, speciesB, &index_input, 3);  

				for(l=0; l<ld->ani; ++l){
					speciesA_dx[2]	= ld->ani_v[l];
					index_chr = ld->cati*ld->cati + k*ld->ani + l;
					
					speciesA_dx[0]	= ld->molec_v[i];	
					speciesA_dx[1]	= ld->cati_v[j]; 
					if(get_index_calc(ld, speciesA_dx, &index_calc, 3)){ // mc
						gg_dx[index_calc]	+= chrfrac_dx[index_chr] * gg_input[index_input];	
						alpha_dx[index_calc]+= chrfrac_dx[index_chr] * alpha_input[index_input]; 
					}

					speciesA_dx[0]	= ld->cati_v[j]; 
					speciesA_dx[1] = ld->molec_v[i]; 
					if(get_index_calc(ld, speciesA_dx, &index_calc, 3)){ // cm
						gg_dx[index_calc]	+= chrfrac_dx[index_chr] * gg_input[index_input+1];	
						alpha_dx[index_calc]+= chrfrac_dx[index_chr] * alpha_input[index_input+1];
					}
				}

				for(l=0; l<ld->cati; ++l){
					speciesA_dx[2] = ld->cati_v[l];

					index_chr = j*ld->cati + l;

					speciesA_dx[0]	= ld->molec_v[i]; 
					speciesA_dx[1]	= ld->ani_v[k]; 
					if(get_index_calc(ld, speciesA_dx, &index_calc, 3)){ // ma
						gg_dx[index_calc]	+= chrfrac_dx[index_chr] * gg_input[index_input];
						alpha_dx[index_calc]+= chrfrac_dx[index_chr] * alpha_input[index_input];
					}

					speciesA_dx[0] = ld->ani_v[k]; 
					speciesA_dx[1] = ld->molec_v[i]; 
					if(get_index_calc(ld, speciesA_dx, &index_calc, 3)){ // am
						gg_dx[index_calc]	+= chrfrac_dx[index_chr] * gg_input[index_input+1];	
						alpha_dx[index_calc]+= chrfrac_dx[index_chr] * alpha_input[index_input+1]; 
					}
				}	
			}
		}
		
		for(j=0; j<ld->cati; ++j){
			for(k=0; k<ld->ani; ++k){
				speciesA_dx[2] = ld->ani_v[k];

				speciesA[0]		= ld->molec_v[i]; 
				speciesA[1]		= ld->cati_v[j]; 
				speciesA_dx[0]	= ld->molec_v[i]; 
				speciesA_dx[1]	= ld->cati_v[j]; 
				if(get_index_calc(ld, speciesA_dx, &index_calc, 3) && get_index_calc(ld, speciesA, &index_chr, 2)){
					tau_dx[index_calc]	= - (gg_dx[index_calc] * alpha[index_chr] / gg[index_chr]
						- alpha_dx[index_calc] * log(gg[index_chr])) / pow(alpha[index_chr], 2);
				}

				speciesA[0]		= ld->cati_v[j]; 
				speciesA[1]		= ld->molec_v[i];
				speciesA_dx[0]	= ld->cati_v[j]; 
				speciesA_dx[1]	= ld->molec_v[i];
				if(get_index_calc(ld, speciesA_dx, &index_calc, 3) && get_index_calc(ld, speciesA, &index_chr, 2)){
					tau_dx[index_calc]	= - (gg_dx[index_calc] * alpha[index_chr] / gg[index_chr]
						- alpha_dx[index_calc] * log(gg[index_chr])) / pow(alpha[index_chr], 2);
				}
			}
		}

		for(j=0; j<ld->ani; ++j){
			for(k=0; k<ld->cati; ++k){
				speciesA_dx[2] = ld->cati_v[k];

				speciesA[0]		= ld->molec_v[i]; 
				speciesA[1]		= ld->ani_v[j]; 
				speciesA_dx[0]	= ld->molec_v[i]; 
				speciesA_dx[1]	= ld->ani_v[j];
				if(get_index_calc(ld, speciesA_dx, &index_calc, 3) && get_index_calc(ld, speciesA, &index_chr, 2)){
					tau_dx[index_calc]	= - (gg_dx[index_calc] * alpha[index_chr] / gg[index_chr]
						- alpha_dx[index_calc] * log(gg[index_chr])) / pow(alpha[index_chr], 2); 
				}

				speciesA[0]		= ld->ani_v[j]; 
				speciesA[1]		= ld->molec_v[i]; 
				speciesA_dx[0] = ld->ani_v[j]; 
				speciesA_dx[1] = ld->molec_v[i];
				if(get_index_calc(ld, speciesA_dx, &index_calc, 3) && get_index_calc(ld, speciesA, &index_chr, 2)){
					tau_dx[index_calc]	= - (gg_dx[index_calc] * alpha[index_chr] / gg[index_chr]
						- alpha_dx[index_calc] * log(gg[index_chr])) / pow(alpha[index_chr], 2);  
				}
			}
 		}
	}

	for(i=0; i<ld->cati; ++i){
		for(j=0; j<ld->ani; ++j){
			speciesC[1] = speciesC[3] = ld->ani_v[j]; 
			speciesC[0] = ld->cati_v[i];
			speciesA_dx[0]	= ld->cati_v[i]; 	
			speciesA_dx[1]	= ld->ani_v[j];
			for(k=0; k<ld->cati; ++k){
				speciesC[2] = ld->cati_v[k];
				get_index_input(ld, speciesC, &index_input, 4); 
				for(l=0; l<ld->cati; ++l){
					speciesA_dx[2]	= ld->cati_v[l];
					index_chr		= k*ld->cati + l;
					if(get_index_calc(ld, speciesA_dx, &index_calc, 3)){
						if(i != k){
							gg_dx[index_calc]	+= chrfrac_dx[index_chr] * gg_input[index_input]; 
							alpha_dx[index_calc]+= chrfrac_dx[index_chr] * alpha_input[index_input]; 
						}else{
							gg_dx[index_calc]	+= chrfrac_dx[index_chr]; 
							alpha_dx[index_calc]+= chrfrac_dx[index_chr] * 0.2; 
						}
					}
				}
			}

			for(k=0; k<ld->cati; ++k){
				speciesA[0]		= ld->cati_v[i]; 	
				speciesA[1]		= ld->ani_v[j];
				speciesA_dx[0]	= ld->cati_v[i]; 	
				speciesA_dx[1]	= ld->ani_v[j];
				speciesA_dx[2]	= ld->cati_v[k];
				if(get_index_calc(ld, speciesA_dx, &index_calc, 3) && get_index_calc(ld, speciesA, &index_chr, 2)){
					tau_dx[index_calc]	= - (gg_dx[index_calc] * alpha[index_chr] / gg[index_chr]
						- alpha_dx[index_calc] * log(gg[index_chr])) / pow(alpha[index_chr], 2); 
				}
			}
			
			speciesC[1] = speciesC[3] = ld->cati_v[i];
			speciesC[0] = ld->ani_v[j]; 
			speciesA_dx[0] = ld->ani_v[j]; 
			speciesA_dx[1] = ld->cati_v[i]; 
			for(k=0; k<ld->ani; ++k){
				speciesC[2] = ld->ani_v[k]; 
				for(l=0; l<ld->ani; ++l){
					speciesA_dx[2]	= ld->ani_v[l];
					index_chr	= ld->cati*ld->cati + k*ld->ani + l;
					if(get_index_calc(ld, speciesA_dx, &index_calc, 3)){
						if(j != k){
							get_index_input(ld, speciesC, &index_input, 4); 
							gg_dx[index_calc]		+= chrfrac_dx[index_chr] * gg_input[index_input]; 
							alpha_dx[index_calc]	+= chrfrac_dx[index_chr] * alpha_input[index_input]; 
						}else{
							gg_dx[index_calc]	+= chrfrac_dx[index_chr]; 
							alpha_dx[index_calc]	+= chrfrac_dx[index_chr] * 0.2; 
						}
					}
				}
			}
			
			for(k=0; k<ld->ani; ++k){
				speciesA[0]		= ld->ani_v[j]; 
				speciesA[1]		= ld->cati_v[i];
				speciesA_dx[0]	= ld->ani_v[j]; 
				speciesA_dx[1]	= ld->cati_v[i];
				speciesA_dx[2]	= ld->ani_v[k];
				if(get_index_calc(ld, speciesA_dx, &index_calc, 3) && get_index_calc(ld, speciesA, &index_chr, 2)){
					tau_dx[index_calc]	= - (gg_dx[index_calc] * alpha[index_chr] / gg[index_chr]
						- alpha_dx[index_calc] * log(gg[index_chr])) / pow(alpha[index_chr], 2); 
				}
			}
		}
	}

	return;
}

void eval_lngm_sr_x(tpropdata *pd, paramstruct *pm, double *lngm_sr, double *lngm_sr_dt, double *lngm_sr_dx)
{
	int i, j, index, index_sum=0, index_sum_dx, lngm_spec[2]; 
	tplibdata* ld = (tplibdata*) pd->libdata; 

	double* sum_xg	= pm->sum_xg; 
	double* sum_xgt	= pm->sum_xgt; 

	if(lngm_sr){

		lngm_spec[1] = -2; 
		for(i=0; i<ld->molec; ++i){
			index			= ld->molec_v[i];
			index_sum		= i; 
			lngm_spec[0]	= ld->molec_v[i]; 
			lngm_sr[index]	= sum_xgt[index_sum] / sum_xg[index_sum]; 
			lngm_sr[index] += sum_prototype(pd, pm, lngm_spec, 'm');
			lngm_sr[index] += sum_prototype(pd, pm, lngm_spec, 'c');
			lngm_sr[index] += sum_prototype(pd, pm, lngm_spec, 'a');
		}

		for(i=0; i<ld->cati; ++i){
			index			= ld->cati_v[i];
			index_sum		= ld->molec + i; 
			lngm_spec[0]	= ld->cati_v[i]; 
			lngm_sr[index]	= sum_xgt[index_sum] / sum_xg[index_sum]; 
			lngm_sr[index] += sum_prototype(pd, pm, lngm_spec, 'm');
			lngm_sr[index] += sum_prototype(pd, pm, lngm_spec, 'a');
		}

		for(i=0; i<ld->ani; ++i){
			index			= ld->ani_v[i]; 
			index_sum		= ld->molec + ld->cati + i; 
			lngm_spec[0]	= ld->ani_v[i]; 
			lngm_sr[index]	= sum_xgt[index_sum] / sum_xg[index_sum]; 
			lngm_sr[index] += sum_prototype(pd, pm, lngm_spec, 'm'); 
			lngm_sr[index] += sum_prototype(pd, pm, lngm_spec, 'c'); 
		}

	}

	if(lngm_sr_dt){

		double* sum_xgt_dt	= pm->sum_xgt_dt;
		double* sum_xg_dt	= pm->sum_xg_dt;

		lngm_spec[1] = -1; 
		for(i=0; i<ld->molec; ++i){
			index				= ld->molec_v[i];
			index_sum			= i; 
			lngm_spec[0]		= ld->molec_v[i]; 
			lngm_sr_dt[index]	= (sum_xgt_dt[index_sum] * sum_xg[index_sum] 
				- sum_xgt[index_sum] * sum_xg_dt[index_sum]) / pow(sum_xg[index_sum], 2); 
			lngm_sr_dt[index]  += sum_prototype(pd, pm, lngm_spec, 'm');
			lngm_sr_dt[index]  += sum_prototype(pd, pm, lngm_spec, 'c');
			lngm_sr_dt[index]  += sum_prototype(pd, pm, lngm_spec, 'a');
		}

		for(i=0; i<ld->cati; ++i){
			index				= ld->cati_v[i];
			index_sum			= ld->molec + i; 
			lngm_spec[0]		= ld->cati_v[i]; 
			lngm_sr_dt[index]	= (sum_xgt_dt[index_sum] * sum_xg[index_sum] 
				- sum_xgt[index_sum] * sum_xg_dt[index_sum]) / pow(sum_xg[index_sum], 2); 
			lngm_sr_dt[index]  += sum_prototype(pd, pm, lngm_spec, 'm');
			lngm_sr_dt[index]  += sum_prototype(pd, pm, lngm_spec, 'a');
		}

		for(i=0; i<ld->ani; ++i){
			index				= ld->ani_v[i]; 
			index_sum			= ld->molec + ld->cati + i; 
			lngm_spec[0]		= ld->ani_v[i]; 
			lngm_sr_dt[index]	= (sum_xgt_dt[index_sum] * sum_xg[index_sum] 
				- sum_xgt[index_sum] * sum_xg_dt[index_sum]) / pow(sum_xg[index_sum], 2);
			lngm_sr_dt[index]  += sum_prototype(pd, pm, lngm_spec, 'm'); 
			lngm_sr_dt[index]  += sum_prototype(pd, pm, lngm_spec, 'c'); 
		}

	}

	if(lngm_sr_dx){

		double* sum_xg_dx	= pm->sum_xg_dx; 
		double* sum_xgt_dx	= pm->sum_xgt_dx; 

		for(i=0; i<ld->molec; ++i){
			lngm_spec[0] = ld->molec_v[i]; 
			for(j=0; j<ld->molec; ++j){
				lngm_spec[1]		= ld->molec_v[j]; 
				index				= i*(ld->molec+ld->ani+ld->cati) + j;  
				index_sum			= i; 
				index_sum_dx		= i*(ld->molec+ld->ani+ld->cati) + j;
				lngm_sr_dx[index]	= (sum_xgt_dx[index_sum_dx] * sum_xg[index_sum] - sum_xgt[index_sum] 
									  * sum_xg_dx[index_sum_dx]) / pow(sum_xg[index_sum], 2); 
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'm');
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'c');
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'a');
			}
			for(j=0; j<ld->cati; ++j){
				lngm_spec[1]		= ld->cati_v[j]; 
				index				= i*(ld->molec+ld->ani+ld->cati) + ld->molec + j;
				index_sum			= i; 
				index_sum_dx		= i*(ld->molec+ld->ani+ld->cati) + ld->molec + j;
				lngm_sr_dx[index]	= (sum_xgt_dx[index_sum_dx] * sum_xg[index_sum] - sum_xgt[index_sum] 
									  * sum_xg_dx[index_sum_dx]) / pow(sum_xg[index_sum], 2); 
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'm');
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'c');
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'a');
 			}
			for(j=0; j<ld->ani; ++j){
				lngm_spec[1]		= ld->ani_v[j]; 
				index				= i*(ld->molec+ld->ani+ld->cati) + ld->molec + ld->cati + j;
				index_sum			= i; 
				index_sum_dx		= i*(ld->molec+ld->ani+ld->cati) + ld->molec + ld->cati + j;
				lngm_sr_dx[index]	= (sum_xgt_dx[index_sum_dx] * sum_xg[index_sum] - sum_xgt[index_sum] 
									  * sum_xg_dx[index_sum_dx]) / pow(sum_xg[index_sum], 2); 
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'm'); 
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'c');
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'a');
 			}
		}
		
		for(i=0; i<ld->cati; ++i){
			lngm_spec[0] = ld->cati_v[i]; 
			for(j=0; j<ld->molec; ++j){
				lngm_spec[1]		= ld->molec_v[j]; 
				index				= ld->molec*(ld->molec+ld->ani+ld->cati) + i*(ld->molec+ld->ani+ld->cati) + j; 
				index_sum			= ld->molec + i; 
				index_sum_dx		= ld->molec*(ld->molec+ld->ani+ld->cati) + i*(ld->molec+ld->ani) + j; 
				lngm_sr_dx[index]	= (sum_xgt_dx[index_sum_dx] * sum_xg[index_sum] - sum_xgt[index_sum] 
									  * sum_xg_dx[index_sum_dx]) / pow(sum_xg[index_sum], 2);  
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'm');
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'a');
			}
			for(j=0; j<ld->cati; ++j){
				lngm_spec[1]		= ld->cati_v[j]; 
				index				= ld->molec*(ld->molec+ld->ani+ld->cati)
									  + i*(ld->molec+ld->ani+ld->cati) + ld->molec + j; 
				lngm_sr_dx[index]	= 0.;  
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'm');
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'a');
			}
			for(j=0; j<ld->ani; ++j){
				lngm_spec[1]		= ld->ani_v[j]; 
				index				= (ld->molec+i)*(ld->molec+ld->ani+ld->cati) + ld->molec + ld->cati + j;  
				index_sum			= ld->molec + i; 
				index_sum_dx		= ld->molec*(ld->molec+ld->ani+ld->cati) + i*(ld->molec+ld->ani) + ld->molec + j;
				lngm_sr_dx[index]	= (sum_xgt_dx[index_sum_dx] * sum_xg[index_sum] - sum_xgt[index_sum] 
									  * sum_xg_dx[index_sum_dx]) / pow(sum_xg[index_sum], 2);  
				lngm_sr_dx[index] += sum_prototype(pd, pm, lngm_spec, 'm');
				lngm_sr_dx[index] += sum_prototype(pd, pm, lngm_spec, 'a');
			}
		}

		for(i=0; i<ld->ani; ++i){
			lngm_spec[0] = ld->ani_v[i];
			for(j=0; j<ld->molec; ++j){
				lngm_spec[1]		= ld->molec_v[j]; 
				index				= (ld->molec+ld->cati+i)*(ld->molec+ld->ani+ld->cati) + j; 
				index_sum			= ld->molec + ld->cati + i; 
				index_sum_dx		= (ld->molec)*(ld->molec+ld->ani+ld->cati) + ld->cati*(ld->molec+ld->ani) 
									  + i*(ld->molec+ld->cati) + j; 
				lngm_sr_dx[index]	= (sum_xgt_dx[index_sum_dx] * sum_xg[index_sum] - sum_xgt[index_sum] 
									  * sum_xg_dx[index_sum_dx]) / pow(sum_xg[index_sum], 2); 
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'm');
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'c');
			}
			for(j=0; j<ld->cati; ++j){
				lngm_spec[1]		= ld->cati_v[j]; 
				index				= (ld->molec+ld->cati+i)*(ld->molec+ld->ani+ld->cati) + ld->molec + j; 
				index_sum			= ld->molec + ld->cati + i; 
				index_sum_dx		= (ld->molec)*(ld->molec+ld->ani+ld->cati) + ld->cati*(ld->molec+ld->ani) 
									  + i*(ld->molec+ld->cati) + ld->molec +j;
				lngm_sr_dx[index]	= (sum_xgt_dx[index_sum_dx] * sum_xg[index_sum] - sum_xgt[index_sum] 
									  * sum_xg_dx[index_sum_dx]) / pow(sum_xg[index_sum], 2); 
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'm');
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'c');
			}
			for(j=0; j<ld->ani; ++j){
				lngm_spec[1]		= ld->ani_v[j]; 
				index				= (ld->molec+ld->cati+i)*(ld->molec+ld->ani+ld->cati) + ld->molec + ld->cati + j;  
				lngm_sr_dx[index]	= 0.;  
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'm');
				lngm_sr_dx[index]  += sum_prototype(pd, pm, lngm_spec, 'c');
			}
		}

	}

	return; 
}

void eval_lngm_sr_x_ref(tpropdata *pd, paramstruct *pm, double *lngm_ref, double *lngm_ref_dt, double *lngm_ref_dx)
{
	int i, j, water, speciesA[2], speciesA_dx[3], index, index_calc, index_dx; 
	double sum_solv=0.; 
	tplibdata *ld = (tplibdata*) pd->libdata;

	double* tau		= pm->tau; 
	double* tau_dt	= pm->tau_dt;
	double* tau_dx	= pm->tau_dx;
	double* alpha	= pm->alpha;
	double* gg		= pm->gg;
	double* gg_dt	= pm->gg_dt;
	double* gg_dx	= pm->gg_dx;
	
	double* x_ref	= new double[ld->molec+ld->cati+ld->ani]; 

	for (i=0; i<pd->nstoff; ++i) {
		if (!strcmpu(pd->stname[i],"H2O") || !strcmpu(pd->stname[i],"WATER") 
			|| !strcmpu(pd->stname[i],"WASSER")) {
			water = i;
			break;
		}
	}

	for(i=0; i<ld->molec+ld->cati+ld->ani; ++i){
		if(!ld->henry[i] && (ld->charge[i] == 0)){
			sum_solv += pm->last_x_eff[ld->molec_v[i]];
			x_ref[i] =  pm->last_x_eff[ld->molec_v[i]]; 
		}else{
			x_ref[i]	= 0.; 
		}
	}

	for(i=0; i<ld->molec+ld->cati+ld->ani; ++i){
		x_ref[i] /= sum_solv; 
	}

	double* sum_xg_mx		= pm->sum_xg_mx; 
	double* sum_xgt_mx		= pm->sum_xgt_mx; 
	double* sum_xg_mx_dt	= pm->sum_xg_mx_dt; 
	double* sum_xgt_mx_dt	= pm->sum_xgt_mx_dt; 

	// molecular species
	for(i=0; i<ld->molec; ++i){
		if(ld->henry[ld->molec_v[i]]){
			if(lngm_ref) lngm_ref[ld->molec_v[i]] = sum_xgt_mx[i] / sum_xg_mx[i];
			if(lngm_ref_dt)	lngm_ref_dt[ld->molec_v[i]] = (sum_xgt_mx_dt[i] * sum_xg_mx[i] 
				- sum_xgt_mx[i] * sum_xg_mx_dt[i]) / pow(sum_xg_mx[i], 2);
			speciesA[0]	= ld->molec_v[i];
			for(j=0; j<ld->molec; ++j){
				speciesA[1] = ld->molec_v[j]; 
				if(get_index_calc(ld, speciesA, &index, 2)){
					if(lngm_ref) lngm_ref[ld->molec_v[i]] += x_ref[ld->molec_v[j]] * gg[index] / sum_xg_mx[j]
						* (tau[index] - sum_xgt_mx[j] / sum_xg_mx[j]); 
					if(lngm_ref_dt) lngm_ref_dt[ld->molec_v[i]] += x_ref[ld->molec_v[j]] * ((gg_dt[index] * sum_xg_mx[j]
						- gg[index] * sum_xg_mx_dt[j]) / pow(sum_xg_mx[j], 2) * (tau[index] - sum_xgt_mx[j] / sum_xg_mx[j])
						+ gg[index] / sum_xg_mx[j] * (tau_dt[index] - (sum_xgt_mx_dt[j] * sum_xg_mx[j] - sum_xgt_mx[j] 
						* sum_xg_mx_dt[j]) / pow(sum_xg_mx[j], 2))); 
				}
			}	 
		}else{
			if(lngm_ref) lngm_ref[ld->molec_v[i]] = 0.;
			if(lngm_ref_dt)	lngm_ref_dt[ld->molec_v[i]] = 0.;
			if(lngm_ref_dx){
				for(j=0; j<ld->molec+ld->cati+ld->ani; ++j){
					index = i*(ld->molec+ld->cati+ld->ani) + j; 
					lngm_ref_dx[index]	= 0.; 
				}
			}
		}
	}

	if(ld->enrtlref == 0){ // infinite dilution aqueous solution reference state for ionic species (unsymmetric)

		for(i=0; i<ld->cati; ++i){
			speciesA[0] = speciesA_dx[0] = water;  
			speciesA[1] = speciesA_dx[1] = ld->cati_v[i]; 
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				if(lngm_ref) lngm_ref[ld->cati_v[i]] = tau[index_calc];
				if(lngm_ref_dt) lngm_ref_dt[ld->cati_v[i]] = tau_dt[index_calc];
				if(lngm_ref_dx){
					for(j=0; j<ld->molec+ld->cati+ld->ani; ++j){
						speciesA_dx[2] = j; 
						index = (ld->molec+i)*(ld->molec+ld->cati+ld->ani) + j; 
						if(get_index_calc(ld, speciesA_dx, &index_dx, 3)){
							lngm_ref_dx[index] = tau_dx[index_dx]; 
						}else{
							lngm_ref_dx[index] = 0.;
						}
					}
				}
			}
			speciesA[0] = speciesA_dx[0] = ld->cati_v[i]; 
			speciesA[1] = speciesA_dx[1] = water;
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				if(lngm_ref) lngm_ref[ld->cati_v[i]] += tau[index_calc] * gg[index_calc]; 
				if(lngm_ref_dt) lngm_ref_dt[ld->cati_v[i]] += tau_dt[index_calc] 
					* gg[index_calc] + tau[index_calc] * gg_dt[index_calc];
				if(lngm_ref_dx){
					for(j=0; j<ld->molec+ld->cati+ld->ani; ++j){
						speciesA_dx[2] = j; 
						index = (ld->molec+i)*(ld->molec+ld->cati+ld->ani) + j; 
						if(get_index_calc(ld, speciesA_dx, &index_dx, 3)){
							lngm_ref_dx[index] += tau_dx[index_dx] * gg[index_calc]
								+ tau[index_calc] * gg_dx[index_dx]; 
						}else{
							lngm_ref_dx[index] += 0.;
						}
					}
				}
			}
		}
		for(i=0; i<ld->ani; ++i){
			speciesA[0] = speciesA_dx[0] = water;  
			speciesA[1] = speciesA_dx[1] = ld->ani_v[i]; 
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				if(lngm_ref) lngm_ref[ld->ani_v[i]] = tau[index_calc]; 
				if(lngm_ref_dt) lngm_ref_dt[ld->ani_v[i]] = tau_dt[index_calc]; 
				if(lngm_ref_dx){
					for(j=0; j<ld->molec+ld->cati+ld->ani; ++j){
						speciesA_dx[2] = j; 
						index = (ld->molec+ld->cati+i)*(ld->molec+ld->cati+ld->ani) + j; 
						if(get_index_calc(ld, speciesA_dx, &index_dx, 3)){
							lngm_ref_dx[index] = tau_dx[index_dx]; 
						}else{
							lngm_ref_dx[index] = 0.;
						}
					}
				}
			}
			speciesA[0] = speciesA_dx[0] = ld->ani_v[i]; 
			speciesA[1] = speciesA_dx[1] = water;
			if(get_index_calc(ld, speciesA, &index_calc, 2)){
				if(lngm_ref) lngm_ref[ld->ani_v[i]] += tau[index_calc] * gg[index_calc]; 
				if(lngm_ref_dt) lngm_ref_dt[ld->ani_v[i]] += tau_dt[index_calc] 
					* gg[index_calc] + tau[index_calc] * gg_dt[index_calc]; 
				if(lngm_ref_dx){
					for(j=0; j<ld->molec+ld->cati+ld->ani; ++j){
						speciesA_dx[2] = j; 
						index = (ld->molec+ld->cati+i)*(ld->molec+ld->cati+ld->ani) + j; 
						if(get_index_calc(ld, speciesA_dx, &index_dx, 3)){
							lngm_ref_dx[index] += tau_dx[index_dx] * gg[index_calc]
								+ tau[index_calc] * gg_dx[index_dx]; 
						}else{
							lngm_ref_dx[index] += 0.;
						}
					}
				}
			}
		}

	}else if(ld->enrtlref == 1){ // pure fused salt reference state for ionic species (symmetric)
		
		double sum_charge=0.; 

		double* sum_xg_ca	= pm->sum_xg_ca; 
		double* sum_xgt_ca	= pm->sum_xgt_ca; 

		for(i=0; i<ld->molec+ld->cati+ld->ani; ++i){
			sum_charge += (double) abs(ld->charge[i]); 
		}

		for(i=0; i<ld->molec+ld->cati+ld->ani; ++i){
			x_ref[i] = (double) abs(ld->charge[i]) / sum_charge; 
		}

		for(i=0; i<ld->cati; ++i){
			speciesA[0] = ld->cati_v[i]; 
			lngm_ref[ld->cati_v[i]] = sum_xgt_ca[ld->molec+i] / sum_xg_ca[ld->molec+i]; 
			for(j=0; j<ld->ani; ++j){
				speciesA[1] = ld->ani_v[j]; 
				if(get_index_calc(ld, speciesA, &index_calc, 2)){
					if(lngm_ref) lngm_ref[ld->cati_v[i]] += x_ref[ld->ani_v[j]] * gg[index_calc] / sum_xg_ca[ld->molec+ld->cati+j]
						* (tau[index_calc] - sum_xgt_ca[ld->molec+ld->cati+j] / sum_xg_ca[ld->molec+ld->cati+j]);
				}
			}
		}

		for(i=0; i<ld->ani; ++i){
			speciesA[0]	= ld->ani_v[i]; 
			lngm_ref[ld->ani_v[i]] = sum_xgt_ca[ld->molec+ld->cati+i] / sum_xg_ca[ld->molec+ld->cati+i]; 
			for(j=0; j<ld->cati; ++j){
				speciesA[1] = ld->cati_v[j]; 
				if(get_index_calc(ld, speciesA, &index_calc, 2)){
					if(lngm_ref) lngm_ref[ld->ani_v[i]] += x_ref[ld->cati_v[j]] * gg[index_calc] / sum_xg_ca[ld->molec+j]
						* (tau[index_calc] - sum_xgt_ca[ld->molec+j] / sum_xg_ca[ld->molec+j]); 
				} 
			}
		}

	}

	delete [] x_ref; 

	return; 
}

void eval_lngm_lr_x(tpropdata *pd, paramstruct *pm, double *lngm_lr, double *lngm_lr_dt, double *lngm_lr_dx)
{
	int i, j, index;
	double aphi, aphi_dt, *aphi_dx=NULL, Ix=0., Ix1=0., Ix2=0., Ix_0=0. ,sum_0=0., value;
	tplibdata *ld=(tplibdata*)pd->libdata;
	double *lngm_born=NULL, *lngm_born_dt=NULL, *lngm_born_dx=NULL;

	static const double cap = 14.9;

	double* x_0				= new double[ld->molec+ld->cati+ld->ani];
	double* Ix_0_n			= new double[ld->cati+ld->ani]; 
	double* x				= pm->last_x;
	double t				= pm->last_t;

	if(!lngm_lr)	aphi	= 0.0;
	if(!lngm_lr_dt) aphi_dt	= 0.0;
	if(lngm_lr_dx)	aphi_dx	= new double[ld->molec+ld->cati+ld->ani];

	get_aphi_x(pd, t, x, &aphi, &aphi_dt, aphi_dx);

	for (i = 0; i < ld->molec+ld->cati+ld->ani; ++i) {
		Ix += x[i] * ld->charge[i] * ld->charge[i];
	} 

	Ix *= 0.5;

	if(Ix<0.) {
		Ix1 = -sqrt(fabs(Ix));
		Ix2 = -pow(fabs(Ix),1.5);
	}else {
		Ix1 = sqrt(Ix);
		Ix2 = pow(Ix,1.5);
	}

	if(lngm_lr){ 
		for(i=0; i<ld->molec; ++i){
			lngm_lr[ld->molec_v[i]] = -aphi * 2 * Ix2 / (1+cap*Ix1);
		}
	}
	if(lngm_lr_dt){
		for(i=0; i<ld->molec; ++i){
			lngm_lr_dt[ld->molec_v[i]] = -aphi_dt * 2 * Ix2 / (1+cap*Ix1);
		}
	}
	if(lngm_lr_dx){
		for(i=0; i<ld->molec; ++i){
			for(j=0; j<ld->molec+ld->cati+ld->ani; ++j){
				index = i*(ld->molec+ld->cati+ld->ani) + j; 
				lngm_lr_dx[index] = -aphi_dx[j] * 2 * Ix2 / (1+cap*Ix1);  
				lngm_lr_dx[index] -= 0.5 * aphi * (2*cap*Ix+3*Ix1) / pow((1+cap*Ix1), 2) * ld->charge[j] * ld->charge[j]; 
			}
		}
	}
	
	if (ld->enrtlref == 0) {

		if(lngm_lr) lngm_born		= new double[ld->molec+ld->cati+ld->ani];
		if(lngm_lr_dt) lngm_born_dt	= new double[ld->molec+ld->cati+ld->ani];
		if(lngm_lr_dx) lngm_born_dx	= new double[(ld->molec+ld->cati+ld->ani)*(ld->molec+ld->cati+ld->ani)];

		eval_lngm_born_x(pd, t, x, lngm_born, lngm_born_dt, lngm_born_dx);

		for (i=0; i<ld->cati; ++i) {
			value = ((2*ld->charge[ld->cati_v[i]]*ld->charge[ld->cati_v[i]])
				/cap*log(fabs(1+cap*Ix1))+((ld->charge[ld->cati_v[i]]*ld->charge[ld->cati_v[i]])*Ix1-2*Ix2)
				/(1+cap*Ix1));
			if(lngm_lr){
				lngm_lr[ld->cati_v[i]] = aphi * value;
				lngm_lr[ld->cati_v[i]] += lngm_born[ld->cati_v[i]];
			} 
			if(lngm_lr_dt){
				lngm_lr_dt[ld->cati_v[i]] = aphi_dt * value;
				lngm_lr_dt[ld->cati_v[i]] += lngm_born_dt[ld->cati_v[i]];
			}
			if(lngm_lr_dx){
				for(j=0; j<ld->molec+ld->cati+ld->ani; ++j){
					index = (ld->molec+i)*(ld->molec+ld->cati+ld->ani) + j; 
					lngm_lr_dx[index] = aphi_dx[j] * value;
					lngm_lr_dx[index] += 0.5 * aphi * (2.*cap*Ix1+3.)*(ld->charge[ld->cati_v[i]]
						* ld->charge[ld->cati_v[i]]-2.*Ix) / (2.*Ix1*pow(1+cap*Ix1,2)) * ld->charge[j] * ld->charge[j];  
					lngm_lr_dx[index] += lngm_born_dx[index];
				}
			}
		}
		for (i=0; i<ld->ani; ++i) {
			value = ((2*ld->charge[ld->ani_v[i]]*ld->charge[ld->ani_v[i]])
				/cap*log(fabs(1+cap*Ix1))+((ld->charge[ld->ani_v[i]]*ld->charge[ld->ani_v[i]])*Ix1-2*Ix2)
				/(1+cap*Ix1));
			if(lngm_lr){
				lngm_lr[ld->ani_v[i]] = aphi * value;
				lngm_lr[ld->ani_v[i]] += lngm_born[ld->ani_v[i]];
			} 
			if(lngm_lr_dt){
				lngm_lr_dt[ld->ani_v[i]] = aphi_dt * value;
				lngm_lr_dt[ld->ani_v[i]] += lngm_born_dt[ld->ani_v[i]];
			}
			if(lngm_lr_dx){
				for(j=0; j<ld->molec+ld->cati+ld->ani; ++j){
					index = (ld->molec+ld->cati+i)*(ld->molec+ld->cati+ld->ani) + j; 
					lngm_lr_dx[index] = aphi_dx[j] * value;
					lngm_lr_dx[index] += 0.5 * aphi * (2.*cap*Ix1+3.)*(ld->charge[ld->ani_v[i]]
						* ld->charge[ld->ani_v[i]]-2.*Ix) / (2.*Ix1*pow(1+cap*Ix1,2)) * ld->charge[j] * ld->charge[j];  
					lngm_lr_dx[index] += lngm_born_dx[index];
				}
			}
		}

	}else if (ld->enrtlref == 1){

		for (i=0; i<ld->ani+ld->cati+ld->molec; ++i) {
			if(ld->charge[i] != 0) sum_0 += x[i];
		}
		for (i=0; i<ld->molec+ld->cati+ld->ani; ++i) {
			if(ld->charge[i] != 0) x_0[i] = x[i] / sum_0;
			else x_0[i] = 0.; 
		}
		for (j = 0; j<ld->molec+ld->cati+ld->ani; ++j) {
			Ix_0 += x_0[i] * ld->charge[i] * ld->charge[i];
		}
		Ix_0 *= 0.5; 
		for(i=0; i<ld->cati; ++i){
			Ix_0_n[i] = 0.; 
			for(j=0; j<ld->cati; ++j){
				if(i != j) Ix_0_n[i] -= ld->charge[ld->cati_v[j]] * ld->charge[ld->cati_v[j]] * x_0[ld->cati_v[j]] / sum_0;
				else Ix_0_n[i] += ld->charge[ld->cati_v[j]] * ld->charge[ld->cati_v[j]] * (1 - x_0[ld->cati_v[j]]) / sum_0;
			}
			for(j=0; j<ld->ani; ++j){
				Ix_0_n[i] -= ld->charge[ld->ani_v[j]] * ld->charge[ld->ani_v[j]] * x_0[ld->ani_v[j]] / sum_0;
			}
		}
		for(i=0; i<ld->ani; ++i){
			for(j=0; j<ld->cati; ++j){
				Ix_0_n[ld->cati+i] += - ld->charge[ld->ani_v[j]] * ld->charge[ld->ani_v[j]] * x_0[ld->ani_v[j]] / sum_0;
			}
			for(j=0; j<ld->ani; ++j){
				if(i != j) Ix_0_n[ld->cati+i] -= ld->charge[ld->ani_v[j]] * ld->charge[ld->ani_v[j]] * x_0[ld->ani_v[j]] / sum_0;
				else Ix_0_n[ld->cati+i] += ld->charge[ld->ani_v[j]] * ld->charge[ld->ani_v[j]] * (1 - x_0[ld->ani_v[j]]) / sum_0;
			}
		}
		for(i=0; i<ld->cati; ++i){
			value = ((2*ld->charge[ld->cati_v[i]]*ld->charge[ld->cati_v[i]])/cap*log(fabs((1+cap*Ix1)
				/(1+cap*sqrt(fabs(Ix_0)))))+((ld->charge[ld->cati_v[i]]*ld->charge[ld->cati_v[i]])*Ix1-2*Ix2)
				/(1+cap*Ix1)-pow(fabs(Ix_0),-0.5)*2*Ix/(1+cap*sqrt(fabs(Ix_0)))*Ix_0_n[i]);
			if(lngm_lr) lngm_lr[ld->cati_v[i]] = aphi * value;
			if(lngm_lr_dt) lngm_lr_dt[ld->cati_v[i]] = aphi_dt * value;
		}
		for(i=0; i<ld->ani; ++i){
			value = ((2*ld->charge[ld->ani_v[i]]*ld->charge[ld->ani_v[i]])/cap*log(fabs((1+cap*Ix1)
				/(1+cap*sqrt(fabs(Ix_0)))))+((ld->charge[ld->ani_v[i]]*ld->charge[ld->ani_v[i]])*Ix1-2*Ix2)
				/(1+cap*Ix1)-pow(fabs(Ix_0),-0.5)*2*Ix/(1+cap*sqrt(fabs(Ix_0)))*Ix_0_n[ld->cati+i]);
			if(lngm_lr) lngm_lr[ld->ani_v[i]] = aphi * value;
			if(lngm_lr_dt) lngm_lr[ld->ani_v[i]] = aphi_dt * value;
		}
		 
	}

	delete [] lngm_born;
	delete [] lngm_born_dt;
	delete [] lngm_born_dx;
	delete [] Ix_0_n; 
	delete [] x_0; 
	delete [] aphi_dx;
}

void eval_lngm_born_x(tpropdata *pd, double t, double *x, double *lngm_born, double *lngm_born_dt, double *lngm_born_dx)
{ 	
	double diec_w, diec_w_dt, prefac, prefac_dt, prefac_dx;
	double mxdiec, mxdiec_dt, *mxdiec_dx=NULL;
	int i, j, index = -1;
	tplibdata *ld=(tplibdata*)pd->libdata;

	if(!lngm_born_dt)	mxdiec_dt	= 0.0;
	if(lngm_born_dx)	mxdiec_dx	= new double[pd->nstoff];

	get_mxdiec_x(pd, t, x, &mxdiec, &mxdiec_dt, mxdiec_dx); 
	
	for (i=0; i<ld->molec; i++) {
		if (!strcmpu(pd->stname[ld->molec_v[i]],"H2O") || !strcmpu(pd->stname[ld->molec_v[i]],"WATER") 
			|| !strcmpu(pd->stname[ld->molec_v[i]],"WASSER") ) {
			index = i;
		}
	}
	
	if(index != -1){

		diec_w		= ld->cpdiec[index*3]+ld->cpdiec[index*3+1]*(1/t-1/ld->cpdiec[index*3+2]);
		if(lngm_born_dt) diec_w_dt	= - ld->cpdiec[index*3+1] / pow(t, 2);

		if(lngm_born) prefac		= ELESU * ELESU / (2*BOLTZESU*t )* (1/mxdiec-1/diec_w) * 0.01;
		if(lngm_born_dt) prefac_dt	=  -ELESU*ELESU / (2*BOLTZESU*t) * ( mxdiec_dt/pow(mxdiec, 2) 
			- diec_w_dt/pow(diec_w, 2) + ( 1 / mxdiec - 1 / diec_w ) / t) * 0.01;
		if(lngm_born_dx) prefac_dx	= ELESU * ELESU / (2*BOLTZESU*t) * 0.01;
	
		for(i=0;i<ld->molec+ld->cati+ld->ani;i++)
		{
			if(ld->charge[i] == 0){
				if(lngm_born)	 lngm_born[i]		= 0.;
				if(lngm_born_dt) lngm_born_dt[i]	= 0.;
				if(lngm_born_dx){
					for(j=0; j<pd->nstoff; ++j){
						index = i*pd->nstoff + j;
						lngm_born_dx[index]	= 0.;
					}
				} 
			}else{
				if(lngm_born) lngm_born[i] = prefac*ld->charge[i]*ld->charge[i]/ld->born[i];
				if(lngm_born_dt) lngm_born_dt[i] = prefac_dt*ld->charge[i]*ld->charge[i]/ld->born[i];
				if(lngm_born_dx){
					for(j=0; j<pd->nstoff; ++j){
						index = i*pd->nstoff+j;
						lngm_born_dx[index] = - prefac_dx * mxdiec_dx[j] / pow(mxdiec, 2)
							* ld->charge[i] * ld->charge[i] / ld->born[i]; 
					}
				}
			}
		}

	}else {

		if(lngm_born) for(i=0;i<ld->molec+ld->cati+ld->ani;i++) lngm_born[i] = 0.;
		if(lngm_born_dt) for(i=0;i<ld->molec+ld->cati+ld->ani;++i) lngm_born_dt[i] = 0.;

	}

	delete [] mxdiec_dx;

	return;
}

int calc_index(tplibdata *ld, int *index, char param_type, int spec0, int spec1, int spec2 = -1)
{
	int i, j, p1, p2; 

	switch(param_type){
		case 'c':
			if(spec2 == -1){
				if((ld->charge[spec0] == 0) && (ld->charge[spec1] == 0) && (spec0 != spec1)){
					for(i=0; i<ld->molec; ++i){
						if(ld->molec_v[i] == spec0){ p1 = i; }
						if(ld->molec_v[i] == spec1){ p2 = i; }
					}
					if(p1 < p2) *index = p1 * (ld->molec - 1) + p2 - 1;
					else if (p2 < p1) *index = p1 * (ld->molec - 1) + p2;
					return TRUE;
				}else if((ld->charge[spec0] == 0) && (ld->charge[spec1] > 0)){
					for(i=0; i<ld->molec; ++i)	if(ld->molec_v[i] == spec0){ p1 = i; break; }
					for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == spec1)	{ p2 = i; break; }	
					*index = ld->molec*(ld->molec - 1) + p2*ld->molec + p1;
					return TRUE;
				}else if((ld->charge[spec0] == 0) && (ld->charge[spec1] < 0)){
					for(i=0; i<ld->molec; ++i)	if(ld->molec_v[i] == spec0){ p1 = i; break; }
					for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == spec1)	{ p2 = i; break; }
					*index = ld->molec*(ld->molec - 1) + ld->cati*ld->molec + p2*ld->molec + p1;
					return TRUE;
				}else if((ld->charge[spec0] > 0) && (ld->charge[spec1] == 0)){
					for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == spec0) { p1 = i; break; }
					for(i=0; i<ld->molec; ++i)	if(ld->molec_v[i] == spec1){ p2 = i; break; }
					*index = ld->molec*((ld->molec -1)+ld->cati+ld->ani) + p1*ld->molec + p2;
					return TRUE;
				}else if((ld->charge[spec0] < 0) && (ld->charge[spec1] == 0)){
					for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == spec0)	{ p1 = i; break; }
					for(i=0; i<ld->molec; ++i)	if(ld->molec_v[i] == spec1){ p2 = i; break; }
					*index = ld->molec*((ld->molec -1)+ld->cati+ld->ani)
						+ ld->cati*ld->molec + p1*ld->molec + p2;
					return TRUE;
				}else if((ld->charge[spec0] > 0) && (ld->charge[spec1] < 0)){
					for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == spec0) { p1 = i; break; }
					for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == spec1)	{ p2 = i; break; }
					*index = ld->molec*((ld->molec -1)+2*(ld->cati+ld->ani)) 
						+ p2*ld->cati + p1; 
					return TRUE;
				}else if((ld->charge[spec0] < 0) && (ld->charge[spec1] > 0)){
					for(i=0; i<ld->ani; ++i)	if(ld->ani_v[i] == spec0)	{ p1 = i; break; }
					for(i=0; i<ld->cati; ++i)	if(ld->cati_v[i] == spec1)	{ p2 = i; break; }
					*index = ld->molec*((ld->molec -1)+2*(ld->cati+ld->ani)) 
						+ ld->ani*ld->cati + p2*ld->ani + p1;
					return TRUE; 
				}
			}else{
			
			}
			break;

		case 'p':

			break;
	}
}

void get_diss_const(tpropdata *pd, paramstruct *pm, double (*int_cp)(double t) , double *ln_K)
{
	int i, j, index; 
	tplibdata *ld = (tplibdata*) pd->libdata; 

	double* dG_0	= new double[ld->nreac];
	double* dH_0	= new double[ld->nreac];
	double* dcp_0	= new double[ld->nreac];

	double t		= pm->last_t; 

	for(i=0; i<ld->nreac; ++i){
		dG_0[i] = 0.; 
		dH_0[i] = 0.; 
		ln_K[i] = 0.; 
		for(j=0; j<ld->molec+ld->cati+ld->ani; ++j){
			index = i*ld->molec+ld->cati+ld->ani + j; 
			dG_0[i] += ld->stoic_coeff[index] * ld->dG_0f[j];  
			dG_0[i] += ld->stoic_coeff[index] * ld->dH_0f[j];
		}
		
		ln_K[i] -= dG_0[i] / (RALG * ld->Tref) + dH_0[i] / RALG * (1 / t - 1 / ld->Tref);  
	}

	

	delete [] dG_0; 
	delete [] dcp_0; 

	return; 
}

int testmolec(tplibdata *ld, int komp)
{
	int i, segcount;

	segcount = ld->molec+ld->cati+ld->ani;

	for (i=0; i<segcount; i++) {
		if (ld->nseg_v[komp][i] != 0  && ld->charge[i] != 0) {
			return 0;
		}
	}
	
	return 1;
}

double sum_prototype(tpropdata *pd, paramstruct *pm, int *lngm_spec, char sum_specs)
{
	int i, sum_count, speciesA[2], speciesA_dx[3], *specs, index_sum=0, index_tg, index_dx;
	double result=0.;
	tplibdata *ld = (tplibdata*) pd->libdata; 

	double* gg		= pm->gg;
	double* tau		= pm->tau;
	double* sum_xg	= pm->sum_xg;
	double* sum_xgt	= pm->sum_xgt;

	double* x_eff	= pm->last_x_eff;

	switch(sum_specs){
		case 'm':
			sum_count	= ld->molec;
			specs		= ld->molec_v; 
			index_sum	= 0; 
			break;

		case 'c':
			sum_count	= ld->cati;
			specs		= ld->cati_v; 
			index_sum	= ld->molec; 
			break;

		case 'a':
			sum_count	= ld->ani;
			specs		= ld->ani_v; 
			index_sum	= ld->molec+ld->cati;
			break; 
	}

	if(lngm_spec[1] == -2){ // regular sum (no derivative)

		speciesA[0]	= lngm_spec[0]; 
		for(i=0; i<sum_count; ++i){
			speciesA[1]	= specs[i]; 
			if(get_index_calc(ld, speciesA, &index_tg, 2)){
				result += x_eff[specs[i]] * gg[index_tg] / sum_xg[index_sum]
					* (tau[index_tg] - sum_xgt[index_sum] / sum_xg[index_sum]);
			}else if(sum_specs == 'm'){
				result -= x_eff[specs[i]] * sum_xgt[index_sum] / pow(sum_xg[index_sum], 2);
			}
			++index_sum; 
		}

	}else if(lngm_spec[1] == -1){ // temperature derivative

		double* gg_dt		= pm->gg_dt;
		double* tau_dt		= pm->tau_dt;
		double* sum_xg_dt	= pm->sum_xg_dt;
		double* sum_xgt_dt	= pm->sum_xgt_dt;

		speciesA[0] = lngm_spec[0];
		for(i=0; i<sum_count; ++i){
			speciesA[1] = specs[i];
			if(get_index_calc(ld, speciesA, &index_tg, 2)){
				result += x_eff[specs[i]] * (gg_dt[index_tg] * sum_xg[index_sum] - gg[index_tg] * sum_xg_dt[index_sum])
					/ pow(sum_xg[index_sum], 2) * (tau[index_tg] - sum_xgt[index_sum] / sum_xg[index_sum])
					+ x_eff[specs[i]] * gg[index_tg] / sum_xg[index_sum] * (tau_dt[index_tg] - (sum_xgt_dt[index_sum]
					* sum_xg[index_sum] - sum_xgt[index_sum] * sum_xg_dt[index_sum]) / pow(sum_xg[index_sum], 2));
			}else if(sum_specs == 'm'){
				result += x_eff[specs[i]] * (2 * sum_xg_dt[index_sum] * sum_xgt[index_sum] 
					- sum_xgt_dt[index_sum] * sum_xg[index_sum]) / pow(sum_xg[index_sum], 3);
			}
			++index_sum;
		}
	
	}else{ // composition derivative 
		
		int index_sum_dx, sum_dx_step, der_spec;
		double* gg_dx		= pm->gg_dx;
		double* tau_dx		= pm->tau_dx;
		double* alpha_dx	= pm->alpha_dx;
		double* sum_xg_dx	= pm->sum_xg_dx;
		double* sum_xgt_dx	= pm->sum_xgt_dx;

		if(ld->charge[lngm_spec[1]] == 0){
			for(i=0; i<ld->molec; ++i){
				if(ld->molec_v[i] == lngm_spec[1]) { der_spec = i; break; }
			}
		}else if(ld->charge[lngm_spec[1]] > 0){
			for(i=0; i<ld->cati; ++i){
				if(ld->cati_v[i] == lngm_spec[1]) { der_spec = i; break; }
			}
		}else if(ld->charge[lngm_spec[1]] < 0){
			for(i=0; i<ld->ani; ++i){
				if(ld->ani_v[i] == lngm_spec[1]) { der_spec = i; break; }
			}
		}
		
		switch(sum_specs){
			case 'm':
				sum_dx_step = ld->molec+ld->cati+ld->ani;
				if(ld->charge[lngm_spec[1]] == 0)	  index_sum_dx	= der_spec;
				else if(ld->charge[lngm_spec[1]] > 0) index_sum_dx	= ld->molec + der_spec;
				else if(ld->charge[lngm_spec[1]] < 0) index_sum_dx	= ld->molec + ld->cati + der_spec;
				break;

			case 'c':
				sum_dx_step = ld->molec+ld->ani;
				if(ld->charge[lngm_spec[1]] == 0) index_sum_dx = ld->molec*(ld->molec+ld->cati+ld->ani) + der_spec;
				else if(ld->charge[lngm_spec[1]] < 0) index_sum_dx = ld->molec*(ld->molec+ld->cati+ld->ani) 
					+ ld->molec + der_spec;
				else index_sum_dx = -1;
				break;

			case 'a':
				sum_dx_step = ld->molec+ld->cati;
				if(ld->charge[lngm_spec[1]] == 0) index_sum_dx = ld->molec*(ld->molec+ld->cati+ld->ani) 
					+ ld->cati*(ld->molec+ld->ani) + der_spec;
				else if(ld->charge[lngm_spec[1]] > 0) index_sum_dx = ld->molec*(ld->molec+ld->cati+ld->ani) 
					+ ld->cati*(ld->molec+ld->ani) + ld->molec + der_spec;
				else index_sum_dx = -1;
				break;
		}

		speciesA[0]	= speciesA_dx[0] = lngm_spec[0];
		speciesA_dx[2] = lngm_spec[1];
		for(i=0; i<sum_count; ++i){
			speciesA[1]	= speciesA_dx[1] = specs[i];
			if((lngm_spec[1] == specs[i]) && get_index_calc(ld, speciesA_dx, &index_dx, 3) 
				 && get_index_calc(ld, speciesA, &index_tg, 2) && (index_sum_dx >= 0)){
				result += ((gg[index_tg] + x_eff[specs[i]] * gg_dx[index_dx]) * sum_xg[index_sum]
					- x_eff[specs[i]] * gg[index_tg] * sum_xg_dx[index_sum_dx]) / pow(sum_xg[index_sum], 2)
					* (tau[index_tg] - sum_xgt[index_sum] / sum_xg[index_sum]);
				result += x_eff[specs[i]] * gg[index_tg] / sum_xg[index_sum] * (tau_dx[index_dx] 
					- (sum_xgt_dx[index_sum_dx] * sum_xg[index_sum] - sum_xgt[index_sum] 
					* sum_xg_dx[index_sum_dx]) / pow(sum_xg[index_sum], 2));
			}else if((lngm_spec[1] == specs[i]) && get_index_calc(ld, speciesA, &index_tg, 2) && (index_sum_dx >= 0)){
				result += (2 * x_eff[specs[i]] * gg[index_tg] * sum_xgt[index_sum] * sum_xg_dx[index_sum_dx] 
					- sum_xg[index_sum] * (gg[index_tg] * sum_xgt[index_sum] + x_eff[specs[i]] * gg[index_tg] * sum_xgt_dx[index_sum_dx]))
					/ pow(sum_xg[index_sum], 3) + tau[index_tg] * gg[index_tg] * (sum_xg[index_sum] - x_eff[specs[i]] 
					* sum_xg_dx[index_sum_dx]) / pow(sum_xg[index_sum], 2);
				/*result += gg[index_tg] * (sum_xg[index_sum] - x_eff[specs[i]] * sum_xg_dx[index_sum_dx])
					/ pow(sum_xg[index_sum], 2) * (tau[index_tg] - sum_xgt[index_sum] / sum_xg[index_sum]);
				result -= x_eff[specs[i]] * gg[index_tg] * (sum_xgt_dx[index_sum_dx] * sum_xg[index_sum]
					- sum_xgt[index_sum] * sum_xgt_dx[index_sum_dx]) / pow(sum_xg[index_sum], 3);*/
			}else if((lngm_spec[1] == specs[i]) && (index_sum_dx >= 0) 
				&& (lngm_spec[0] == specs[i]) && (ld->charge[lngm_spec[0]] == 0)){
				result += x_eff[specs[i]] * (2 * sum_xgt[index_sum] * sum_xg_dx[index_sum_dx] 
					- sum_xgt_dx[index_sum_dx] * sum_xg[index_sum] ) / pow(sum_xg[index_sum], 3);
				result -= sum_xgt[index_sum] / pow(sum_xg[index_sum], 2); 
			}else if((lngm_spec[1] == specs[i]) && get_index_calc(ld, speciesA, &index_tg, 2)){
				if (ld->charge[specs[i]] == 0) result += gg[index_tg] / sum_xg[index_sum] 
					* (tau[index_tg] - sum_xgt[index_sum] / sum_xg[index_sum]); 
				else result += abs(ld->charge[specs[i]]) * gg[index_tg] / sum_xg[index_sum] 
					* (tau[index_tg] - sum_xgt[index_sum] / sum_xg[index_sum]);
			}else if(get_index_calc(ld, speciesA_dx, &index_dx, 3)  
				&& get_index_calc(ld, speciesA, &index_tg, 2) && (index_sum_dx >= 0)){
				result += x_eff[specs[i]] * (gg_dx[index_dx] * sum_xg[index_sum] - gg[index_tg] 
					* sum_xg_dx[index_sum_dx]) / pow(sum_xg[index_sum], 2) * (tau[index_tg] 
					- sum_xgt[index_sum] / sum_xg[index_sum] );
				result += x_eff[specs[i]] * gg[index_tg] / sum_xg[index_sum] * (tau_dx[index_dx] 
					- (sum_xgt_dx[index_sum_dx] * sum_xg[index_sum] - sum_xgt[index_sum] 
					* sum_xg_dx[index_sum_dx]) / pow(sum_xg[index_sum], 2));
			}else if(get_index_calc(ld, speciesA, &index_tg, 2) && (index_sum_dx >= 0)){
				result -= x_eff[specs[i]] * gg[index_tg] * sum_xg_dx[index_sum_dx] / pow(sum_xg[index_sum], 2)
					* (tau[index_tg] - sum_xgt[index_sum] / sum_xg[index_sum]);
				result -=  x_eff[specs[i]] * gg[index_tg] / sum_xg[index_sum] * (sum_xgt_dx[index_sum_dx] 
					* sum_xg[index_sum] - sum_xgt[index_sum] * sum_xg_dx[index_sum_dx]) / pow(sum_xg[index_sum], 2);
			}else if((index_sum_dx >= 0) && (lngm_spec[0] == specs[i]) && (ld->charge[lngm_spec[0]] == 0)){
				result += x_eff[specs[i]] * (2 * sum_xgt[index_sum] * sum_xg_dx[index_sum_dx] 
					- sum_xgt_dx[index_sum_dx] * sum_xg[index_sum] ) / pow(sum_xg[index_sum], 3);
			}
			if(index_sum_dx >= 0) index_sum_dx += sum_dx_step;
			++index_sum;
		}

	} 

	return result;
}

void get_aphi_x(tpropdata *pd, double t, double *x, double *aphi, double *aphi_dt, double *aphi_dx)
{
	int i, j; 
	double mxdens=0., mxdiec=0., value=0., value1=0., mix_mm=0., aphi_reg;
	double pow_mxdiec, sqrt_mxdens; 
	double mxdens_dt, mxdiec_dt;
	double *x_dx=NULL, *mxdens_dx=NULL, *mxdiec_dx=NULL, *xm=NULL; 
	tplibdata *ld = (tplibdata*) pd->libdata;

	if(!aphi_dt){
		mxdens_dt	= 0.0;
		mxdiec_dt	= 0.0;
	}  
	if(aphi_dx){
		mxdens_dx	= new double[pd->nstoff];
		mxdiec_dx	= new double[pd->nstoff];
		x_dx		= new double[pd->nstoff];		
		xm			= new double[pd->nstoff];
	}  

	for (i=0; i<pd->nstoff; ++i) {
		if(aphi_dx) x_dx[i]	= 0.;
		if (testmolec(ld, i) && !ld->henry[i]) {
			value += x[i];
			if(aphi_dx){
				for(j=0; j<pd->nstoff; ++j){
					if (testmolec(ld, j) && !ld->henry[j] && (i != j)) {
						x_dx[i] += x[j];
					}
				}
			}
		}
	}

	for (i=0; i<pd->nstoff; ++i) {
		if (testmolec(ld, i) && !ld->henry[i]) {
			mix_mm		+=  x[i] / value * ld->mm[i];
			if(aphi_dx) x_dx[i]	/= pow(value , 2);
		}
	}

	if(aphi_dx){
		for(i=0; i<pd->nstoff; ++i){
			if(testmolec(ld, i) && !ld->henry[i]){
				xm[i]	= x_dx[i] * ld->mm[i]; 
				for(j=0; j<pd->nstoff; ++j){
					if((i != j) && testmolec(ld, j) && !ld->henry[j]){
						xm[i]	-= x[j] * ld->mm[j] / pow(value, 2);
					} 
				} 
			}
		}
	}

	get_mxdens_x(pd, t, x, &mxdens, &mxdens_dt, mxdens_dx); 
	get_mxdiec_x(pd, t, x, &mxdiec, &mxdiec_dt, mxdiec_dx); 

	if(mxdiec<0.)	pow_mxdiec	= -pow(fabs(ELESU*ELESU/(BOLTZESU*mxdiec*t)),1.5);
	else			pow_mxdiec	= pow(ELESU*ELESU/(BOLTZESU*mxdiec*t),1.5);
		
	if(mxdens<0.)	sqrt_mxdens = -sqrt(fabs(2*PI*AVOG*mxdens/1000));	
	else			sqrt_mxdens = sqrt(2*PI*AVOG*mxdens/1000);
	
	//calculate  Debye-Hueckel parameter
	aphi_reg = - 1. / 3. * sqrt(1000/mix_mm) * sqrt_mxdens * pow_mxdiec;

	if(aphi) *aphi = aphi_reg; 

	if(aphi_dt){
		value		= 2*PI*AVOG*mxdens/1000; 
		value1		= pow(ELESU, 2)/(BOLTZESU*mxdiec*t);

		if((mxdiec*mxdens)<0.)	value = -sqrt(fabs(value1*value));
		else					value = sqrt(value1*value);

		*aphi_dt	= 1. / 3. * PI * AVOG / 1000 * pow_mxdiec / sqrt_mxdens * mxdens_dt; 
		*aphi_dt	-= 0.5 * value *  pow(ELESU/(mxdiec*t),2) / BOLTZESU * (mxdiec_dt*t+mxdiec);
		*aphi_dt	*= - sqrt(1000/mix_mm);
	}
	
	if(aphi_dx){
		value		= 2*PI*AVOG*mxdens/mix_mm; 
		value1		= pow(ELESU, 2)/(BOLTZESU*mxdiec*t);

		if((mxdiec*mxdens)<0.)	value = -sqrt(fabs(value1*value));
		else					value = sqrt(value1*value);

		if(mxdens<0.)	value1 = -sqrt(fabs(2*PI*AVOG*mxdens/mix_mm));
		else			value1 = sqrt(2*PI*AVOG*mxdens/mix_mm);

		for(i=0; i<pd->nstoff; ++i){
			if(testmolec(ld, i) && !ld->henry[i]){
				aphi_dx[i]	= - 1. / 3. * (PI * AVOG  * pow_mxdiec / value1 * (mxdens_dx[i] * mix_mm - xm[i] * mxdens) 
					/ pow(mix_mm, 2) - 3. / 2. * value * pow(ELESU/mxdiec,2) / (BOLTZESU*t) * mxdiec_dx[i]);
			}else{
				aphi_dx[i]	= 0.; 
			}
		}
	}

	delete [] mxdens_dx;
	delete [] mxdiec_dx; 
	delete [] xm; 
	delete [] x_dx;

	return;
}

void get_mxdiec_x(tpropdata* pd, double t, double* x, double* mxdiec, double* mxdiec_dt, double* mxdiec_dx)
{
	int i, j; 
	double	value = 0., value1 = 0., mix_mm = 0., sum_x=0., xme=0., xm=0.; 

	tplibdata* ld = (tplibdata*) pd->libdata; 

	double* x_dx	= new double[pd->nstoff];

	for (i=0; i<pd->nstoff; ++i) {
		x_dx[i] = 0.;
		if (testmolec(ld, i) && !ld->henry[i]) {
			sum_x += x[i];
			for(j=0; j<pd->nstoff; ++j){
				if (testmolec(ld, j) && !ld->henry[j] && (i != j)) {
					x_dx[i] += x[j];
				}
			}
		}
	}
	
	for (i=0; i<pd->nstoff; ++i) {
		if (testmolec(ld, i) && !ld->henry[i]) {
			mix_mm		+= x[i] / sum_x * ld->mm[i];
			x_dx[i]		/= pow(sum_x, 2);
		}
	}

	value = 0.; 
	value1	= 0.; 
	if(mxdiec_dt) *mxdiec_dt = 0.; 

	for(i=0; i<pd->nstoff; ++i){
		if (testmolec(ld, i) && !ld->henry[i]) {
			value1		+= x[i] / sum_x * ld->mm[i] * (ld->cpdiec[i*3] + ld->cpdiec[i*3+1] 
								* (1/t - 1/ld->cpdiec[i*3+2]));
			value		+= x[i] / sum_x * ld->mm[i];
			if(mxdiec_dt) *mxdiec_dt	-= x[i]/sum_x*ld->mm[i]*ld->cpdiec[i*3+1]/pow(t, 2);
		}
	}

	if(mxdiec)	  *mxdiec		= value1 / value; 
	if(mxdiec_dt) *mxdiec_dt	/= value; 

	if(mxdiec_dx){
		for(i=0; i<pd->nstoff; ++i){
			mxdiec_dx[i]	= 0.;  
			if(testmolec(ld, i) && !ld->henry[i]){
				xme	= x_dx[i] * ld->mm[i] * (ld->cpdiec[i*3] + ld->cpdiec[i*3+1] * (1/t - 1/ld->cpdiec[i*3+2]));
				xm	= x_dx[i] * ld->mm[i]; 
				for(j=0; j<pd->nstoff; ++j){
					if((i != j) && testmolec(ld, j) && !ld->henry[j]){
						xme	-= x[j] * ld->mm[j] / pow(sum_x, 2) *  (ld->cpdiec[j*3] + ld->cpdiec[j*3+1] 
							* (1/t - 1/ld->cpdiec[j*3+2]));
						xm	-= x[j] * ld->mm[j] / pow(sum_x, 2);
					} 
				}
				mxdiec_dx[i] = (xme*value - xm*value1) / pow(value, 2); 
			}
		}
	}

	delete [] x_dx; 

	return; 
}

void get_mxdens_x(tpropdata* pd, double t, double* x, double* mxdens, double* mxdens_dt, double* mxdens_dx)
{
	int i, j; 
	double value = 0., value1=0., tr=0., mix_mm=0., sum_x=0.; 
	tplibdata* ld = (tplibdata*) pd->libdata; 

	double* x_dx		= new double[pd->nstoff];

	for (i=0; i<pd->nstoff; ++i) {
		x_dx[i] = 0.;
		if (testmolec(ld, i) && !ld->henry[i]) {
			sum_x += x[i];
			for(j=0; j<pd->nstoff; ++j){
				if (testmolec(ld, j) && !ld->henry[j] && (i != j)) {
					x_dx[i] += x[j];
				}
			}
		}
	}
	
	for (i=0; i<pd->nstoff; ++i) {
		if (testmolec(ld, i) && !ld->henry[i]) {
			mix_mm		+= x[i] / sum_x * ld->mm[i];
			x_dx[i]		/= pow(sum_x, 2);
		}
	}

	value = 0.; 

	for(i=0; i<pd->nstoff; ++i){
		if(testmolec(ld, i) && !ld->henry[i]){
			tr = t / ld->tc[i];
			if(tr>0.99) tr = 1- (1-0.99) * exp((0.99-tr) / (1-0.99));
			value	+= x[i] / sum_x / ld->mm[i] * RALG * ld->tc[i] / ld->pc[i] 
						* pow(ld->rktzra[i], 1 + pow((1-tr),2./7.)) * 1.e6;
			if(mxdens_dt) value1	-= 2. / 7. * x[i] / sum_x / ld->mm[i] * RALG / ld->pc[i] 
						* log(ld->rktzra[i]) * pow(ld->rktzra[i],1+pow((1-tr),2./7.)) 
						* pow(1-tr, -5./7.) * 1.e6;
		}
	}
	
	if(mxdens)    *mxdens		= 1 / value; 
	if(mxdens_dt) *mxdens_dt	= - value1 / pow(value, 2);

	if(mxdens_dx){
		for(i=0; i<pd->nstoff; ++i){
			mxdens_dx[i]	= 0.;
			if(testmolec(ld, i) && !ld->henry[i]){
				tr = t/ld->tc[i];
				if(tr>0.99) tr = 1- (1-0.99) * exp((0.99-tr) / (1-0.99));
				mxdens_dx[i] -= x_dx[i] * RALG / ld->mm[i] * ld->tc[i] / ld->pc[i] 
						* pow(ld->rktzra[i], 1 + pow((1-tr),2./7.)) * 1.e6;
				for(j=0; j<pd->nstoff; ++j){
					if(testmolec(ld, j) && (i != j) && !ld->henry[j]){
						tr = t / ld->tc[j];
						if(tr>0.99) tr = 1- (1-0.99) * exp((0.99-tr) / (1-0.99));
						mxdens_dx[i] += x[j] / pow(sum_x, 2) * RALG / ld->mm[j] * ld->tc[j] / ld->pc[j] 
							* pow(ld->rktzra[j], 1 + pow((1-tr),2./7.)) * 1.e6;
					}		
				}
				mxdens_dx[i] /= pow(value, 2);
			}
		} 
	} 

	delete [] x_dx; 

	return; 
}
