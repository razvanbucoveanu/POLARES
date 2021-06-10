/*
 * cross_sections.h
 *
 *  Created on: Dec 4, 2015
 *      Author: razvan
 */

#ifndef SRC_CROSS_SECTIONS_H_
#define SRC_CROSS_SECTIONS_H_

#include "form_factors.h"
#include "melem.h"
#include "melem_pol.h"
#include "virtual_corrections.h"
#include <fstream>
#include "parameters.h"

namespace POLARES {

class Parameters;
class Interpolation;

class Cross_Sections {

protected:
	const Interpolation* interpolation;
	const Parameters* param;
	//variables to calculate
	mutable long double s, melem_interf, melem2, sigma_born;
	//form factors
	mutable double gpe, gpm, gpze, gpzm, gae, f1, f2, f1z, f2z;
	//auxiliary variables
	mutable long double tau, eps, thl, l2, a, b, l1k, cospsi, x, phig;
	Form_factors ff;
	Melem melem_2nd;
	Melem_pol melem_2nd_pol;
//	FILE* fout;
//	mutable int i;
	mutable long double gamma_loop;
	mutable long double m, m2, m4;
	mutable long double m_me, m_mu, m_tau;
	mutable long double M, M2, M4;
	mutable double alpha, pi, pi2, gf;

public:
	Virtual_Corrections VC;

	int set_param(const Parameters* param, const Interpolation* interpolation);

	//For proton target
	long double asymm_born_test (const long double Q2)const;
	long double crsect_born (const long double Q2)const;
	long double crsect_elastic (const long double Q2)const;
	long double crsect_born_thl (const long double Q2)const;
	long double crsect_elastic_thl (const long double Q2)const;
	long double interf_born (const long double Q2)const;
	long double interf_born_thl (const long double Q2)const;
	long double interf_elastic (const long double Q2)const;
	long double interf_elastic_thl (const long double Q2)const;
	long double crsect_brems_1st (const long double en1, const long double thl, const long double eg,
			const long double thg, const long double phig)const;
	long double crsect_brems_1st_test (const long double en1, const long double thl, const long double eg,
			const long double thg, const long double phig)const;
	long double crsect_brems_1st_ps2 (const long double en1, const long double thl, const long double eg,
			const long double thg, const long double phig)const;
	long double crsect_brems_1st_sg_diff (const long double en1, const long double thl, const long double eg,
			const long double thg, const long double phig)const;
	long double crsect_brems_2nd (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double crsect_brems_2nd_phig1 (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double crsect_brems_2nd_l1k1_1 (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double crsect_brems_2nd_l1k2_1 (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double crsect_brems_2nd_l2k1_1 (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double crsect_brems_2nd_l2k2_1 (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double crsect_brems_2nd_l1k1_2 (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double crsect_brems_2nd_l1k2_2 (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double crsect_brems_2nd_l2k1_2 (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double crsect_brems_2nd_l2k2_2 (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double crsect_brems_2nd_add (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double crsect_brems_2nd_add_phig1 (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double interf_brems_1st (const long double en1, const long double thl, const long double eg, const long double thg)const;
	long double interf_brems_1st_ps2 (const long double en1, const long double thl, const long double eg,
			const long double thg, const long double phig)const;
	long double interf_brems_1st_test (const long double en1, const long double thl, const long double eg, const long double thg)const;
	long double interf_brems_2nd (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double crsect_brems_2nd_pol_add (const long double en1, const long double thl, const long double eg, const long double thg,
			const long double phig, const long double eg1, const long double thg1, const long double phig1)const;
	long double crsect_brems_1st_hadr (const long double en1, const long double thl, const long double eg,
			const long double thg, const long double phig)const;
	long double crsect_brems_1st_hadr_interf (const long double en1, const long double thl, const long double eg,
			const long double thg, const long double phig)const;

	//For Carbon-12 target
	long double crsect_born_carbon (const long double Q2)const;
	long double crsect_elastic_carbon (const long double Q2)const;
	long double crsect_brems_1st_carbon (const long double en1, const long double thl, const long double eg,
			const long double thg, const long double phig)const;
	long double interf_born_carbon (const long double Q2)const;
	long double interf_elastic_carbon (const long double Q2)const;
	long double interf_brems_1st_carbon (const long double en1, const long double thl, const long double eg, const long double thg)const;

	long double crsect_born_carbon_thl (const long double Q2)const;
	long double crsect_elastic_carbon_thl (const long double Q2)const;
	long double interf_born_carbon_thl (const long double Q2)const;
	long double interf_elastic_carbon_thl (const long double Q2)const;
	long double asymm_born_carbon_test (const long double Q2)const;

	//For electron target
	long double crsect_born_e(const long double Q2)const;
	long double crsect_brems_1st_e (const long double en1, const long double thl, const long double eg,
				const long double thg, const long double phig)const;

Cross_Sections();

};

}

#endif /* SRC_CROSS_SECTIONS_H_ */
