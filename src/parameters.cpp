/*
 * input_set.cpp
 *
 *  Created on: Dec 3, 2015
 *      Author: razvan
 */

#include <cmath>
#include "const.h"
#include "ReadConfigFile.h"    //Code that reads from the input file v2.1 (Copyright (c) 2004 Richard J. Wagner)
#include <iostream>
#include <stdlib.h>
#include "parameters.h"
#include "config.h"

using namespace constants;
using namespace POLARES;

Parameters::Parameters()
: SEED(1), MINEVAL(0), MAXEVAL_LO(1e7), MAXEVAL_1st(1e7), maxeval_1st_aux(0),
MAXEVAL_gamma_loop(1e7), MAXEVAL_2nd(1e8), MAXEVAL_2nd_add(1e8),
EPSREL(1e-100), no_cores(4), GRIDNO_elastic(0), GRIDNO_brems(1),
GRIDNO_brems_hadr(13), GRIDNO_brems_hadr_interf(14),
GRIDNO_brems_1st(2), GRIDNO_brems_test(3), GRIDNO_brems_interf(6),
GRIDNO_brems_2nd(7), GRIDNO_brems_2nd_l1k1(8), GRIDNO_brems_2nd_l1k2(9),
GRIDNO_brems_2nd_l2k1(10), GRIDNO_brems_2nd_l2k2(11),
GRIDNO_brems_interf_2nd(12), GRIDNO_brems_interf_2nd_1(13),
GRIDNO_brems_interf_2nd_2(14), GRIDNO_brems_2nd_add_l1k1(15),
GRIDNO_brems_2nd_add_l1k2(16), GRIDNO_brems_2nd_add_l2k1(17),
GRIDNO_brems_2nd_add_l2k2(18), NSTART(1000),
NINCREASE(500), NBATCH(1000), NNEW(100000), NMIN(200), FLATNESS(5), P(1),
m(me), m2(me2), M(mpr), M2(mpr2), sw2(sw2_msbar), Z_lepton(-1.), Z_target(1.),
en(0.155), l1(0), thl(0), Q2(0), eps(0), thg(0), aux(0), Delta_E(0),
Delta_eps(0), mu_dim(0), lambda(m2), check_input(0) {
}

int Parameters::read_input(const Input& input) {
//		    function for reading the input file

		check_input = 1;

		if (input.input_file.empty()) {

//			basic input
			en = input.E;
			min[theta_l] = input.thl_min * pi/180.;
			min[theta_l_deg] = input.thl_min;
			max[theta_l] = input.thl_max * pi/180.;
			max[theta_l_deg] = input.thl_max;
			min[Q2_elastic] = input.Q2min;
			max[Q2_elastic] = input.Q2max;
			min[E_gamma] = input.Delta;
			min[E_gamma_prime] = min[E_gamma];
			max[E_gamma] = input.E_gamma_max;
			max[E_gamma_prime] = max[E_gamma];
			min[theta_gamma] = input.thg_min * pi/180.;
			max[theta_gamma] = input.thg_max * pi/180.;
			min[E_prime] = input.E_prime_min;
			max[E_prime] = input.E_prime_max;
			P = input.polarization;
			sw2 = input.sw2;
//			if (input.events != 0) flag[ps] = input.events;

//			input for cuba library
			FLATNESS = input.flatness;
			SEED = input.seed;
			EPSREL = input.epsrel;
			no_cores = input.no_cores;
			MAXEVAL_LO = input.no_eval_LO;
			MAXEVAL_1st = input.no_eval_1st;
			MAXEVAL_2nd = input.no_eval_2nd;
			MAXEVAL_gamma_loop = input.no_eval_gamma_loop;
			MAXEVAL_2nd_add = input.no_eval_2nd_add;
			MINEVAL = input.no_min_eval;
			NSTART = input.nstart;
			NBATCH = input.nbatch;
			NINCREASE = input.nincrease;
			NNEW = input.nnew;
			NMIN = input.nmin;

//			input for flags
			flag[form_factors] = input.flag[input.form_factors];
			flag[echo_input] = input.flag[input.echo_input];
			flag[brems] = input.flag[input.brems];
			flag[brems_add] = input.flag[input.brems_add];
			flag[brems_hadr] = input.flag[input.brems_hadr];
			flag[PS] = input.flag[input.PS];
			flag[asymmetry] = input.flag[input.asymmetry];
			flag[LO] = input.flag[input.LO];
			flag[order] = input.flag[input.order];
			flag[hadr_corr] = input.flag[input.hadr_corr];
			flag[kappa_weak] = input.flag[input.kappa_weak];
			flag[vac_pol] = input.flag[input.vac_pol];
			flag[tpe] = input.flag[input.tpe];
			flag[lepton] = input.flag[input.lepton];
			flag[target] = input.flag[input.target];
			flag[int_method] = input.flag[input.int_method];
			flag[GL] = input.flag[input.GL];
			flag[int_output] = input.flag[input.int_output];
			flag[cuts_born] = input.flag[input.cuts_born];

//			Input for testing gamma_loop
			Delta_eps = input.Delta_eps;
			mu_dim = input.mu_dim;
			lambda = input.lambda;

//			Input for event generator
			min[E] = input.E_min;
			max[E] = input.E_max;
			Delta_E = input.Delta_E;

//			input for differential cross sections
			thl = input.thl_deg * pi/180.;
			thg = input.thg_deg * pi/180.;

		}
		else {
		 ConfigFile CFile(input.input_file + ".in");

		 CFile.readInto<double>(en,"Incident Lepton Energy", input.E);

		 CFile.readInto<double>(min[theta_l_deg],"theta_l min", input.thl_min);
		 min[theta_l] = (min[theta_l_deg]*pi)/180.;

		 CFile.readInto<double>(max[theta_l_deg],"theta_l max", input.thl_max);
		 max[theta_l] = (max[theta_l_deg]*pi)/180.;

		 CFile.readInto<double>(min[Q2_elastic],"Q^2 min", input.Q2min);

		 CFile.readInto<double>(max[Q2_elastic],"Q^2 max", input.Q2max);

		 CFile.readInto<double>(min[E_gamma],"Delta", input.Delta);

		 CFile.readInto<double>(min[E_gamma_prime],"Delta", input.Delta);

		 CFile.readInto<double>(max[E_gamma],"E_gamma max", input.E_gamma_max);

		 CFile.readInto<double>(max[E_gamma_prime],"E_gamma max", input.E_gamma_max);

		 CFile.readInto<double>(min[E_p],"E_p min", input.E_p_min);

		 CFile.readInto<double>(min[E],"E min", input.E_min);

		 CFile.readInto<double>(max[E],"E max", input.E_max);

		 CFile.readInto<double>(Delta_E,"Delta E", input.Delta_E);

		 CFile.readInto<double>(min[theta_gamma_deg],"theta_gamma min", input.thg_min);
		 min[theta_gamma] = min[theta_gamma_deg]*pi/180.;

		 CFile.readInto<double>(max[theta_gamma_deg],"theta_gamma max", input.thg_max);
		 max[theta_gamma] = max[theta_gamma_deg]*pi/180.;

		 CFile.readInto<double>(min[E_prime],"E' min", input.E_prime_min);

		 CFile.readInto<double>(max[E_prime],"E' max", input.E_prime_max);

		 CFile.readInto<double>(P,"Polarization",input.polarization);

		 CFile.readInto<double>(sw2,"sin2thetaW",input.sw2);

		 CFile.readInto<double>(FLATNESS,"FLATNESS", input.flatness);

		 CFile.readInto<double>(SEED,"Seed", input.seed);

		 CFile.readInto<double>(EPSREL,"Relative Accuracy", input.epsrel);

		 CFile.readInto<double>(Delta_eps,"Delta_eps", input.Delta_eps);

		 CFile.readInto<double>(mu_dim,"mu_dim", input.mu_dim);

		 CFile.readInto<double>(lambda,"lambda", input.lambda);

		 CFile.readInto<int>(no_cores,"Number of cores", input.no_cores);

		 CFile.readInto<long long int>(MAXEVAL_LO,"Maximum Number of Evaluations LO", input.no_eval_LO);

		 CFile.readInto<long long int>(MAXEVAL_1st,"Maximum Number of Evaluations 1st", input.no_eval_1st);

		 CFile.readInto<long long int>(MAXEVAL_gamma_loop,"Maximum Number of Evaluations gamma_loop", input.no_eval_gamma_loop);

		 CFile.readInto<long long int>(MAXEVAL_2nd,"Maximum Number of Evaluations 2nd", input.no_eval_2nd);

		 CFile.readInto<long long int>(MAXEVAL_2nd_add,"Maximum Number of Evaluations 2nd sg finite", input.no_eval_2nd_add);

		 CFile.readInto<long long int>(MINEVAL,"Minimum Number of Evaluations", input.no_min_eval);

		 CFile.readInto<long long int>(NSTART,"NSTART", input.nstart);

		 CFile.readInto<long long int>(NINCREASE,"NINCREASE", input.nincrease);

		 CFile.readInto<long long int>(NBATCH,"NBATCH", input.nbatch);

		 CFile.readInto<long long int>(NNEW,"NNEW", input.nnew);

		 CFile.readInto<long long int>(NMIN,"NMIN", input.nmin);

		 CFile.readInto<int>(flag[form_factors], "Form Factors", input.flag[input.form_factors]);

		 CFile.readInto<int>(flag[echo_input], "Echo Input", input.flag[input.echo_input]);

//		 CFile.readInto<string>(output_file, "Output", "POLARES");

		 CFile.readInto<int>(flag[brems], "Bremsstrahlung Type", input.flag[input.brems]);

		 CFile.readInto<int>(flag[Delta_cut], "Type of Delta Cut", input.flag[input.Delta_cut]);

		 CFile.readInto<int>(flag[brems_add], "Bremsstrahlung Add", input.flag[input.brems_add]);

		 CFile.readInto<int>(flag[brems_hadr], "Hadronic Radiation", input.flag[input.brems_hadr]);

		 CFile.readInto<int>(flag[hadr_corr], "Hadronic corrections", input.flag[input.hadr_corr]);

		 CFile.readInto<int>(flag[PS], "Phase Space Parametrization", input.flag[input.PS]);

		 CFile.readInto<int>(flag[asymmetry], "Asymmetry", input.flag[input.asymmetry]);

		 CFile.readInto<int>(flag[LO], "Leading Order", input.flag[input.LO]);

		 CFile.readInto<int>(flag[order], "Order SP_loop", input.flag[input.order]);

		 CFile.readInto<int>(flag[kappa_weak], "Kappa Form Factor", input.flag[input.kappa_weak]);

		 CFile.readInto<int>(flag[vac_pol], "Vacuum Polarization", input.flag[input.vac_pol]);

		 CFile.readInto<int>(flag[tpe], "Two-photon exchange", input.flag[input.tpe]);

		 CFile.readInto<int>(flag[lepton], "Incident Lepton", input.flag[input.lepton]);

		 CFile.readInto<int>(flag[target], "Target Particle", input.flag[input.target]);

		 CFile.readInto<int>(flag[int_method], "Integration method", input.flag[input.int_method]);

		 CFile.readInto<int>(flag[GL], "Gamma Loop", input.flag[input.GL]);

		 CFile.readInto<int>(flag[int_output],"Integration Output level", input.flag[input.int_output]);

		 CFile.readInto<int>(flag[cuts_born], "Type of Cuts", input.flag[input.cuts_born]);
		}

		 return 0;
}

int Parameters::final_param(const Input& input) {

	if (check_input == 0) {
		std::cerr << "Warning! Input was not set!\n\n";
		exit (EXIT_FAILURE);
	}

	if (flag[form_factors] < 0 || flag[form_factors] > 6) {
		std::cerr<<"Warning! Wrong input for flag[form_factors]! Value changed to default!\n\n";
		flag[form_factors] = 0;
	}

	if (flag[form_factors] > 4 && flag[target] == 0) {
		std::cerr<<"Warning! Form factors for carbon-12 selected and "
				"proton as target particle! Value changed to default\n\n";
		flag[form_factors] = 0;
	}

	if (flag[form_factors] <= 4 && flag[target] == 1) {
		std::cerr<<"Warning! Form factors for proton selected and "
				"carbon-12 as target particle! Value changed to default\n\n";
		flag[form_factors] = 5;
	}

	if (flag[echo_input] < 0 || flag[echo_input] > 1) {
		std::cerr<<"Warning! Wrong input for flag[echo_input]! Value changed to default!\n\n";
		flag[echo_input] = 0;
	}

	if (flag[brems] < 0 || flag[brems] > 15) {
		std::cerr<<"Warning! Wrong input for flag[brems]! Value changed to default!\n\n";
		flag[brems] = 0;
	}
	else if (flag[brems] > 3) {
		std::cerr<<"Warning! This value inserted for flag[brems] is only for testing. "
				"Use with care!\n\n";
	}

	if (flag[brems_add] < 0 || flag[brems_add] > 1) {
		std::cerr<<"Warning! Wrong input for flag[brems_add]! Value changed to default!\n\n";
		flag[brems_add] = 0;
	}

	if (flag[Delta_cut] < 0 || flag[Delta_cut] > 1) {
		std::cerr<<"Warning! Wrong input for flag[Delta_cut]! Value changed to default!\n\n";
		flag[Delta_cut] = 0;
	}

	if (flag[brems_hadr] < 0 || flag[brems_hadr] > 3) {
		std::cerr<<"Warning! Wrong input for flag[brems_hadr]! Value changed to default!\n\n";
		flag[brems_hadr] = 0;
	}

	if (flag[hadr_corr] < 0 || flag[hadr_corr] > 3) {
		std::cerr<<"Warning! Wrong input for flag[brems_hadr]! Value changed to default!\n\n";
		flag[brems_hadr] = 0;
	}

	if (flag[PS] < 0 || flag[PS] > 1) {
		std::cerr<<"Warning! Wrong input for flag[PS]! Value changed to default!\n\n";
		flag[PS] = 0;
	}

	if (flag[asymmetry] < 0 || flag[asymmetry] > 1) {
		std::cerr<<"Warning! Wrong input for flag[asymmetry]! Value changed to default!\n\n";
		flag[asymmetry] = 0;
	}

	if (flag[LO] < 0 || flag[LO] > 1) {
		std::cerr<<"Warning! Wrong input for flag[LO]! Value changed to default!\n\n";
		flag[LO] = 0;
	}

	if (flag[order] < 0 || flag[order] > 2) {
		std::cerr<<"Warning! Wrong input for flag[order]! Value changed to default!\n\n";
		flag[order] = 0;
	}

	if (flag[kappa_weak] < 0 || flag[kappa_weak] > 1) {
		std::cerr<<"Warning! Wrong input for flag[kappa_weak]! Value changed to default!\n\n";
		flag[kappa_weak] = 0;
	}

	if (flag[vac_pol] < 0 || flag[vac_pol] > 5) {
		std::cerr<<"Warning! Wrong input for flag[vac_pol]! Value changed to default!\n\n";
		flag[vac_pol] = 0;
	}

	if (flag[tpe] < 0 || flag[tpe] > 2) {
		std::cerr<<"Warning! Wrong input for flag[tpe]! Value changed to default!\n\n";
		flag[tpe] = 0;
	}

	if (flag[tpe] == 2 && en != 0.155) {
		std::cerr<<"Warning! Option 2 fo flag[tpe] is valid only for E=155 MeV! "
				"Value changed to default!\n\n";
		flag[tpe] = 0;
	}

	if (flag[lepton] < 0 || flag[lepton] > 3) {
		std::cerr<<"Warning! Wrong input for flag[lepton]! Value changed to default!\n\n";
		flag[lepton] = 0;
	}

	if (flag[target] < 0 || flag[target] > 2) {
		std::cerr<<"Warning! Wrong input for flag[target]! Value changed to default!\n\n";
		flag[target] = 0;
	}

	if (flag[int_method] < 0 || flag[int_method] > 2) {
		std::cerr<<"Warning! Wrong input for flag[int_method]! Value changed to default!\n\n";
		flag[int_method] = 0;
	}

	if (flag[GL] < 0 || flag[GL] > 1) {
		std::cerr<<"Warning! Wrong input for flag[Gl]! Value changed to default!\n\n";
		flag[GL] = 0;
	}

	if (flag[int_output] < 0 || flag[int_output] > 2) {
		std::cerr<<"Warning! Wrong input for flag[int_output]! Value changed to default!\n\n";
		flag[int_output] = 0;
	}

	if (flag[cuts_born] < 0 || flag[cuts_born] > 1) {
		std::cerr<<"Warning! Wrong input for flag[cuts_born]! Value changed to default!\n\n";
		flag[cuts_born] = 0;
	}

	if (flag[lepton] == 0 || flag[lepton] == 1) m = me;
	if (flag[lepton] == 2 || flag[lepton] == 3) m = m_mu;
	if (flag[lepton] == 0 || flag[lepton] == 2) Z_lepton = -1.;
	if (flag[lepton] == 1 || flag[lepton] == 3) Z_lepton = 1.;
	if (flag[target] == 0) {
		M = mpr;
		Z_target = 1.;
	}
	if (flag[target] == 1) {
		M = m_carbon12;
		Z_target = 6.;
	}
	if (flag[target] == 2) {
		M = me;
		Z_target = -1.;
	}

	m2 = pow(m,2.);
	M2 = pow(M,2.);

//	if (en > 10.) std::cerr <<"Warning! Integration not tested for energies above "
//			"10 GeV. Use with care!\n\n";
	if (en <= 0.) {
		std::cerr<<"Warning! Wrong input for E! Value changed to default!\n\n";
		en=0.155;
	}

	if (flag[cuts_born] == 0) {
		if ((min[theta_l_deg] < 0. || min[theta_l_deg] > 180.)
				|| (max[theta_l_deg] < 0. || max[theta_l_deg] > 180.)
				|| (max[theta_l_deg] <= min[theta_l_deg])) {
			std::cerr<<"Warning! Wrong input for thl_min and/or thl_max! "
					"Values changed to default!\n\n";
			min[theta_l_deg] = 25.;
			min[theta_l] = 25.*pi/180.;
			max[theta_l_deg] = 45.;
			max[theta_l] = 45.*pi/180.;
		}
	}
	if (flag[cuts_born] == 1) {
		double Q2max_abs = 4.*M2*(en*en-m2)/(M*(2.*en+M)+m2);
		double Q2min = 4.*en*en*pow(sin(25.*pi/2./180.),2.)/(1.+2.*en*pow(sin(25.*pi/2./180.),2.)/M);
		double Q2max = 4.*en*en*pow(sin(45.*pi/2./180.),2.)/(1.+2.*en*pow(sin(45.*pi/2./180.),2.)/M);
		if (min[Q2_elastic] < 0. || max[Q2_elastic] > Q2max_abs || min[Q2_elastic] > max[Q2_elastic]) {
			std::cerr<<"Warning! Wrong input for Q2_min and/or Q2_max! "
					"Values changed to default!\n\n";
			min[Q2_elastic] = Q2min;
			max[Q2_elastic] = Q2max;
		}
	}

	if (min[E_gamma] <= 0. || min[E_gamma] >= en - min[E_prime]) {
		std::cerr<<"Warning! Wrong input for Delta! Value changed to default!\n\n";
		min[E_gamma] = 0.01;
		min[E_gamma_prime] = 0.01;
	}
	if (min[E_prime] < m || min[E_prime] >= en) {
		std::cerr<<"Warning! Wrong input for E'min! Value changed to default!\n\n";
		min[E_prime] = m;
	}
	if (max[E_prime] > en - min[E_gamma] || max[E_prime] <= 0) {
//		std::cerr<<"Warning! Wrong input for E'max! Value changed to default!";
		max[E_prime] = en;
	}
	if (max[E_gamma] > en - min[E_prime]) max[E_gamma] = en - min[E_prime];
	if (max[E_gamma] < min[E_gamma] || max[E_gamma] <= 0) max[E_gamma] = en - min[E_prime];
	if (max[E_gamma_prime] < min[E_gamma_prime] || max[E_gamma_prime] <= 0) max[E_gamma_prime] = en - min[E_prime];

	if (min[E_p] < M) min[E_p] = M;
	if (max[E_p] > en + M) max [E_p] = en + M;

	GRIDNO_elastic = int(en*10000.);
	GRIDNO_brems_1st = 100000 + int(en*10000.);
	GRIDNO_brems = 4;
	GRIDNO_brems_test = 1;
//	GRIDNO_brems_l1k = 700000 + int(en*10000.);
//	GRIDNO_brems_l2k = 800000 + int(en*10000.);
	GRIDNO_brems_interf = 3;
//	GRIDNO_brems_2nd = 200000 + int(en*10000.);
	GRIDNO_brems_2nd = 5;
	GRIDNO_brems_2nd_l1k1 = 200000 + int(en*10000.);
	GRIDNO_brems_2nd_l1k2 = 300000 + int(en*10000.);
	GRIDNO_brems_2nd_l2k1 = 400000 + int(en*10000.);
	GRIDNO_brems_2nd_l2k2 = 500000 + int(en*10000.);
	GRIDNO_brems_2nd_add_l1k1 = 11;
	GRIDNO_brems_2nd_add_l1k2 = 12;
	GRIDNO_brems_2nd_add_l2k1 = 13;
	GRIDNO_brems_2nd_add_l2k2 = 14;
	GRIDNO_brems_interf_2nd = 2;
	GRIDNO_brems_interf_2nd_1 = 8;
	GRIDNO_brems_interf_2nd_2 = 9;
	GRIDNO_brems_hadr = 6;
	GRIDNO_brems_hadr_interf = 7;

	// HS(17.03.2017)
	 if (P < -1.) P=-1.;
	 if (P >  1.) P= 1.;
	// HS-end

	 if (sw2 <= 0.) {
		 std::cerr<<"Warning! Wrong input for sin2thetaW! Value changed to default!\n\n";
		 sw2 = sw2_msbar;
	 }

	 if (input.Delta_E != 0.0001 && input.Delta_E != 0.001 && input.Delta_E != 0.01) {
		 std::cerr << "Warning! Initialization failed! The increment for the energy "
				 "of the incident electron has to be 0.0001, 0.001 or 0.01 GeV";
		 exit (EXIT_FAILURE);
	 }

//	if (GRIDNO_elastic >= 1602) GRIDNO_elastic = 0;
//	if (GRIDNO_brems >= 3602) GRIDNO_brems = 0;

	l1 = sqrt(pow(en, 2.) - m2);

	if (flag[cuts_born] == 0) {
		if (min[theta_l_deg] > 90)
			min[Q2_elastic] = -((2*(en - m)*(en + m)*M*(-M + en*(-1 + pow(cos(min[theta_l]),2))) -
					2*pow(pow(cos(min[theta_l]),2)*pow(en - m,2)*pow(en + m,2)*M2*
							((-1 + pow(cos(min[theta_l]),2))*m2 + M2),0.5))*pow(2*en*M
									-(-1 + pow(cos(min[theta_l]),2))*pow(en,2) + pow(cos(min[theta_l]),2)*m2 + M2,-1));
		else
			min[Q2_elastic] = -2*((en - m)*(en + m)*M*(-M + en*(-1 + pow(cos(min[theta_l]),2))) +
					pow(pow(cos(min[theta_l]),2)*pow(en - m,2)*pow(en + m,2)*M2*
							((-1 + pow(cos(min[theta_l]),2))*m2 + M2),0.5))*pow(2*en*M
									-(-1 + pow(cos(min[theta_l]),2))*pow(en,2) + pow(cos(min[theta_l]),2)*m2 + M2,-1);

		if (max[theta_l_deg] > 90)
			max[Q2_elastic] = -((2*(en - m)*(en + m)*M*(-M + en*(-1 + pow(cos(max[theta_l]),2))) -
					2*pow(pow(cos(max[theta_l]),2)*pow(en - m,2)*pow(en + m,2)*M2*
							((-1 + pow(cos(max[theta_l]),2))*m2 + M2),0.5))* pow(2*en*M
									-(-1 + pow(cos(max[theta_l]),2))*pow(en,2) + pow(cos(max[theta_l]),2)*m2 + M2,-1));
		else
			max[Q2_elastic] = -2*((en - m)*(en + m)*M*(-M + en*(-1 + pow(cos(max[theta_l]),2))) +
					pow(pow(cos(max[theta_l]),2)*pow(en - m,2)*pow(en + m,2)*M2*
							((-1 + pow(cos(max[theta_l]),2))*m2 + M2),0.5))*pow(2*en*M
									-(-1 + pow(cos(max[theta_l]),2))*pow(en,2) + pow(cos(max[theta_l]),2)*m2 + M2,-1);
	}
	else {
		min[theta_l] = acos((-min[Q2_elastic] * (1. + en/M) + 2.*pow(en,2.) - 2.*m2) /
				(2.*l1*sqrt(pow(en-min[Q2_elastic]/(2.*M),2.)-m2)));
		min[theta_l_deg] = min[theta_l]*180./pi;
		max[theta_l] = acos((-max[Q2_elastic] * (1. + en/M) + 2.*pow(en,2.) - 2.*m2) /
				(2.*l1*sqrt(pow(en-max[Q2_elastic]/(2.*M),2.)-m2)));
		max[theta_l_deg] = max[theta_l]*180./pi;
	}

	if (flag[tpe] != 0) {
		eps = 1. / (1. + 2.*(1. + max[Q2_elastic]/(4.*M2))*pow(tan(max[theta_l]/2.),2.));
		if (eps <= 1e-4) {
			std::cerr<<"Warning! Epsilon is out of allowed range for interpolation. "
					"The two-photon exchange correction is disabled \n\n";
			flag[tpe] = 0;
		}

		eps = 1. / (1. + 2.*(1. + min[Q2_elastic]/(4.*M2))*pow(tan(min[theta_l]/2.),2.));
		if (eps >= 0.9999) {
			std::cerr<<"Warning! Epsilon is out of allowed range for interpolation. "
					"The two-photon exchange correction is disabled \n\n";
			flag[tpe] = 0;
		}

		if (max[Q2_elastic] >= 3 || min[Q2_elastic] <= 1e-7) {
			std::cerr<<"Warning! Q^2 is out of allowed range for interpolation. "
					"The two-photon exchange correction is disabled \n\n";
			flag[tpe] = 0;
		}
	}

	if (flag[vac_pol] == 3) {
		if (min[Q2_elastic] <= 0. || max[Q2_elastic] > 1e6) {
			std::cerr<<"Warning! Q^2 is out of allowed range for interpolation."
					" The hadronic vacuum polarization correction is disabled \n\n";
			flag[vac_pol] = 2;
		}
	}

	if (flag[vac_pol] == 4) {
		if (min[Q2_elastic] <= 0. || max[Q2_elastic] >= 1e6) {
			std::cerr<<"Warning! Q^2 is out of allowed range for interpolation."
					" The hadronic vacuum polarization correction is disabled \n\n";
			flag[vac_pol] = 2;
		}
	}

	if (flag[vac_pol] == 5) {
		if (min[Q2_elastic] <= 0. || max[Q2_elastic] >= 9992.0015999999996) {
			std::cerr<<"Warning! Q^2 is out of allowed range for interpolation."
					" The hadronic vacuum polarization correction is disabled \n\n";
			flag[vac_pol] = 2;
		}
	}

//#ifndef POLARES_USE_LOOPTOOLS
//	if (flag[order] == 2 && flag[GL] == 0) {
//	std::cerr << "Warning! Second order corrections and exact calculation for 1loop*1HP selected, "
//			"but Looptools was not included."
//			" Please check the manual for how to include Looptools or select 'Gamma Loop=1' \n\n";
//	exit (EXIT_FAILURE);
//	}
//#endif

	if (flag[order] == 2 && flag[order] == 2) {
		maxeval_1st_aux = MAXEVAL_1st;
		MAXEVAL_1st = MAXEVAL_gamma_loop;
	}

	min[cos_thl] = cos(max[theta_l]);
	max[cos_thl] = cos(min[theta_l]);
	min[cos_thg] = cos(max[theta_gamma]);
	max[cos_thg] = cos(min[theta_gamma]);

	output_file = input.input_file;

		return 0;
}

int Parameters::set_thl(const double thl_deg) {

	thl = thl_deg * pi / 180.;

	if (thl_deg > 90.)
		Q2 = -((2*(en - m)*(en + m)*M*(-M + en*(-1 + pow(cos(thl),2))) -
				2*pow(pow(cos(thl),2)*pow(en - m,2)*pow(en + m,2)*M2*
						((-1 + pow(cos(thl),2))*m2 + M2),0.5))*pow(2*en*M
								-(-1 + pow(cos(thl),2))*pow(en,2) + pow(cos(thl),2)*m2 + M2,-1));
	else
		Q2 = -2*((en - m)*(en + m)*M*(-M + en*(-1 + pow(cos(thl),2))) +
				pow(pow(cos(thl),2)*pow(en - m,2)*pow(en + m,2)*M2*
						((-1 + pow(cos(thl),2))*m2 + M2),0.5))*pow(2*en*M
								-(-1 + pow(cos(thl),2))*pow(en,2) + pow(cos(thl),2)*m2 + M2,-1);

	eps = 1. / (1. + 2.*(1. + Q2/(4.*M2))*pow(tan(thl/2.),2));

	if (flag[tpe] != 0) {
		eps = 1. / (1. + 2.*(1. + max[Q2_elastic]/(4.*M2))*pow(tan(max[theta_l]/2.),2.));
		if (eps <= 1e-4) {
			std::cerr<<"Warning! Epsilon is out of allowed range for interpolation. The two-photon exchange correction is disabled \n\n";
			flag[tpe] = 0;
		}

		eps = 1. / (1. + 2.*(1. + min[Q2_elastic]/(4.*M2))*pow(tan(min[theta_l]/2.),2.));
		if (eps >= 0.9999) {
			std::cerr<<"Warning! Epsilon is out of allowed range for interpolation. The two-photon exchange correction is disabled \n\n";
			flag[tpe] = 0;
		}

		if (max[Q2_elastic] >= 3 || min[Q2_elastic] <= 1e-7) {
			std::cerr<<"Warning! Q^2 is out of allowed range for interpolation. The two-photon exchange correction is disabled \n\n";
			flag[tpe] = 0;
		}
	}

	if (flag[vac_pol] != 0) {
		if (Q2 <= 4.6e-7 || Q2 > 1e6) {
			std::cerr<<"Warning! Q^2 is out of allowed range for interpolation. The hadronic vacuum polarization correction is disabled \n\n";
			flag[vac_pol] = 2;
		}
	}

	return 0;
}
