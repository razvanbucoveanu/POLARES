/*
 * classes.h
 *
 *  Created on: May 13, 2015
 *      Author: razvan
 */

#ifndef SRC_IO_CLASSES_H_
#define SRC_IO_CLASSES_H_

#include <vector>
#include <string>
#include "const.h"

namespace POLARES {

class Final_State {

public:
	double E;
	double E_prime_l, theta_l, phi_l;
	double E_p, theta_p, phi_p;
	double E_gamma, theta_gamma, phi_gamma;
	double E_gamma_prime, theta_gamma_prime, phi_gamma_prime;
	double Q2, weight, avg_weight, sigma_diff;
	double l_2[4], k_1[4], k_2[4], p_2[4], l_1[4], p_1[4];
	long long int seed, event_no, event_type, failed_ev, acc_ev;

	Final_State():
		E(0), E_prime_l(0), theta_l(0), phi_l(0),
		E_p(0), theta_p(0), phi_p(0),
		E_gamma(0), theta_gamma(0), phi_gamma(0),
		E_gamma_prime(0), theta_gamma_prime(0), phi_gamma_prime(0),
		Q2(0), weight(0), avg_weight(0), sigma_diff(0),
		seed(0), event_no(0), event_type(0), failed_ev(0), acc_ev(0) {};
};

class Output {

public:
	double sigma_unpol_born, sigma_pol_born, sigma_unpol_elastic_1st,
	sigma_unpol_elastic_2nd, sigma_pol_elastic_1st, sigma_pol_elastic_2nd,
	sigma_unpol_inelastic_1st, sigma_unpol_inelastic_2nd,
	sigma_unpol_inelastic_loop, sigma_pol_inelastic_loop,
	sigma_pol_inelastic_1st, sigma_pol_inelastic_2nd,
	sigma_unpol_1st, sigma_pol_1st, sigma_unpol_2nd,
	sigma_pol_2nd, sigma_born, sigma_1st, sigma_2nd,
	asymm_born, asymm_1st, asymm_2nd,
	rel_asymm_1st, rel_asymm_2nd, sigma_unpol_2nd_add, sigma_pol_2nd_add,
	sigma_unpol_inelastic_1st_hadr, sigma_unpol_inelastic_1st_hadr_interf;
	std::vector<double> sigma_unpol_1st_vect, sigma_pol_1st_vect, ev_brems_1st,
	sigma_unpol_1st_elastic_vect, sigma_unpol_1st_inelastic_vect,
	sigma_unpol_2nd_vect, sigma_1st_vect, sigma_2nd_vect,
	sigma_pol_2nd_vect, ev_brems_2nd, ev_brems_2nd_l1k1, ev_brems_2nd_l1k2,
	ev_brems_2nd_l2k1, ev_brems_2nd_l2k2;
	double shiftQ2, rel_shiftQ2, Q2;

	Output():
		sigma_unpol_born(0), sigma_pol_born(0), sigma_unpol_elastic_1st(0),
		sigma_unpol_elastic_2nd(0), sigma_pol_elastic_1st(0),
		sigma_pol_elastic_2nd(0), sigma_unpol_inelastic_1st(0),
		sigma_unpol_inelastic_2nd(0), sigma_unpol_inelastic_loop(0),
		sigma_pol_inelastic_loop(0), sigma_pol_inelastic_1st(0),
		sigma_pol_inelastic_2nd(0), sigma_unpol_1st(0), sigma_pol_1st(0),
		sigma_unpol_2nd(0), sigma_pol_2nd(0), sigma_born(0),sigma_1st(0),
		sigma_2nd(0), asymm_born(0), asymm_1st(0), asymm_2nd(0), rel_asymm_1st(0),
		rel_asymm_2nd(0), sigma_unpol_2nd_add(0), sigma_pol_2nd_add(0),
		sigma_unpol_inelastic_1st_hadr(0), sigma_unpol_inelastic_1st_hadr_interf(0),
		shiftQ2(0), rel_shiftQ2(0), Q2(0) {};
};

class Input {

public:
	double thl_min, thl_max, polarization, E_prime_min, E_prime_max, Delta, Delta1,
	thl_deg, thg_deg, E, E_min, E_max, Delta_E, thg_min, thg_max, E_gamma_max,
	Q2min, Q2max, E_p_min, E_p_max;
	double Delta_eps, mu_dim, lambda, sw2;
	double flatness, seed, epsrel;
	int no_cores;
	long long int no_eval_2nd, no_eval_1st, no_eval_gamma_loop, no_eval_2nd_add,
	no_eval_LO, no_min_eval, nstart, nincrease, nbatch, nnew, nmin;
	enum Flags {vac_pol, order, brems, asymmetry, LO, kappa_weak, brems_add,
		cuts_born, ps, form_factors, tpe, echo_input, lepton, int_method, GL, PS,
		brems_hadr, target, int_output, Delta_cut, hadr_corr};
	int flag[30];
	std::string input_file;
	std::string output_file;
//	double me;

	Input()
	:thl_min(25.), thl_max(45.), polarization(1.), E_prime_min(0.045), E_prime_max(1e10),
	 Delta(0.01), Delta1(0.01), thl_deg(35.), thg_deg(35.), E(0.155), E_min(0.155), E_max(0.155),
	 Delta_E(0.001), thg_min(0.), thg_max(180.), E_gamma_max(1e10), Q2min(0.0044), Q2max(0.0134),
	 E_p_min(0.), E_p_max(1e10), Delta_eps(0.), mu_dim(constants::me2), lambda(constants::me2),
	 sw2(constants::sw2_msbar), flatness(5.),
	 seed(1), epsrel(0), no_cores(4), no_eval_2nd(1e9), no_eval_1st(1e8), no_eval_gamma_loop(1e8),
	 no_eval_2nd_add(1e8), no_eval_LO(1e7), no_min_eval(100000), nstart(1000), nincrease(500), nbatch(1000),
	 nnew(100000), nmin(200), flag{3,2,2,1,1,1,1}, input_file(), output_file() {};
};
}
#endif /* SRC_IO_CLASSES_H_ */
