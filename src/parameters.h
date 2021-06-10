/*
 * set_input.h
 *
 *  Created on: Dec 2, 2015
 *      Author: razvan
 */

#ifndef SRC_PARAMETERS_H_
#define SRC_PARAMETERS_H_

#include "IO_classes.h"
#include <string>

namespace POLARES {

class Parameters {

protected:
	int check_input;

public:
//	cuba input parameters
	double SEED;
	long long MINEVAL;
	long long int MAXEVAL_LO;
	long long int MAXEVAL_1st;
	long long int maxeval_1st_aux;
	long long int MAXEVAL_gamma_loop;
	long long int MAXEVAL_2nd;
	long long int MAXEVAL_2nd_add;
	double EPSREL;
	int no_cores;
	int GRIDNO_elastic;
	int GRIDNO_brems;
	int GRIDNO_brems_hadr;
	int GRIDNO_brems_hadr_interf;
	int GRIDNO_brems_1st;
	int GRIDNO_brems_test;
//	int GRIDNO_brems_l1k;
//	int GRIDNO_brems_l2k;
	int GRIDNO_brems_interf;
	int GRIDNO_brems_2nd;
	int GRIDNO_brems_2nd_l1k1;
	int GRIDNO_brems_2nd_l1k2;
	int GRIDNO_brems_2nd_l2k1;
	int GRIDNO_brems_2nd_l2k2;
	int GRIDNO_brems_2nd_add_l1k1;
	int GRIDNO_brems_2nd_add_l1k2;
	int GRIDNO_brems_2nd_add_l2k1;
	int GRIDNO_brems_2nd_add_l2k2;
	int GRIDNO_brems_interf_2nd;
	int GRIDNO_brems_interf_2nd_1;
	int GRIDNO_brems_interf_2nd_2;
//	Vegas
	long long int NSTART, NINCREASE, NBATCH;
//	Suave
	long long int NNEW, NMIN;
	double FLATNESS;

	enum Flags {asymmetry, cuts_born, polarization, ps, form_factors, vac_pol, tpe, echo_input,
		lepton, brems, brems_add, shiftQ2_2nd, LO, order, int_method, kappa_weak, GL, PS,
		brems_hadr, target, int_output, Delta_cut, hadr_corr};
	enum Cuts {E, E_prime, theta_l, theta_l_deg, cos_thl, E_gamma, theta_gamma, cos_thg,
		theta_gamma_deg, phi_gamma, Q2_elastic, E_gamma_prime, E_p};
//	 Flags
	int flag[30];
	double min[20], max[20];
	double P;                            //degree of polarization
	double m, m2, M, M2;			 	 //lepton and nucleon masses
	double sw2;				 			 //the weak mixing angle and the running fine-structure constant
	double Z_lepton, Z_target;
	double en, l1;                       //energy of the incoming electron beam
	double thl, Q2, eps, thg;            //scattering angle and Q^2 (only to be used for diff_thl and shiftQ2)
	double aux, Delta_E;
	double Delta_eps, mu_dim, lambda;
	std::string output_file;

	int read_input(const Input& input);
	int final_param(const Input& input);
	int set_thl(const double thl_deg);

	Parameters();

};

}

#endif /* SRC_PARAMETERS_H_ */
