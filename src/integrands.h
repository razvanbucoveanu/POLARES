/*
 * integrands.h
 *
 *  Created on: Dec 4, 2015
 *      Author: razvan
 */

#ifndef SRC_INTEGRANDS_H_
#define SRC_INTEGRANDS_H_

#include "cross_sections.h"
#include "gsl_rand.h"
//#include "virtual_corrections.h"
//#include <fstream>

namespace POLARES {

class Interpolation;
class Parameters;
class Final_State;

class Integrands {

protected:
	//classes
	const Parameters* param;
	const Interpolation* interpolation;
	Final_State* FS;
	Rand rand;

	//variables for born cross-sections
	//auxiliary
	mutable long double Q2, Jacobian, max_wgt;
	//limits
	long double Q2min, Q2max;

	//variables for hard-photon bremsstrahlung
	//auxiliary
	mutable int events_no, acc_events;
	mutable int acc_events_1, acc_events_2, acc_events_3, acc_events_4, acc_events_5;
	mutable bool event_count;
	mutable long double a, b, c, d, f, n2, l1, l2, a1, a2, y, Jacobian_eg, Jacobian_eg1;
	//calculated limits
	mutable long double thgmin, thgmax, egmin, egmax, thlmin;
	mutable long double eg1min, eg1max, thg1min, thg1max, phig1min, phig1max, phig1min_cos, phig1max_cos;
	//generated values
	mutable long double thl, thg, phig, en1, eg, eg1, thg1, phig1;
	mutable long double m, m2, M, M2, pi;
//	FILE* fout;

public:
	Cross_Sections CS;
//	Virtual_Corrections VC;

	void set_param (const Parameters* input, const Interpolation* interpolation, Final_State* FS);
	void set_param (const Parameters* input, const Interpolation* interpolation);

//	Leading order integrands
	int integrand_born (const double xx[], double ff[], const double weight[])const;
	int integrand_elastic (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_born (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_elastic (const double xx[], double ff[], const double weight[])const;
//	First order integrands
	int integrand_brems_1st (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_test (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_hadr (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_hadr_interf (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_l1k (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_l2k (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_2diff_v1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_2diff_v2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_ps2_2diff (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_ps2_3diff (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_ps2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_sg_diff (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_brems_1st (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_thl (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_thl_ps2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_hadr_thl (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_hadr_interf_thl (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_brems_1st_thl (const double xx[], double ff[], const double weight[])const;
//	Second order integrands
	int integrand_brems_2nd (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l1k1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l1k2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l2k1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l2k2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l1k1_1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l1k2_1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l2k1_1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l2k2_1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l1k1_2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l1k2_2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l2k1_2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l2k2_2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l1k1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l1k2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l2k1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l2k2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l1k1_1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l1k2_1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l2k1_1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l2k2_1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l1k1_2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l1k2_2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l2k1_2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l2k2_2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_thl (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_thl_1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_thl_2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_thl (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_thl_1 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_thl_2 (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_2diff (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_1diff (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_2diff (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l1k1_2diff (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l1k2_2diff (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l2k1_2diff (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_add_l2k2_2diff (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l1k1_2diff (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l1k2_2diff (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l2k1_2diff (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_l2k2_2diff (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_pol_add (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_2nd_pol_add_thl (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_brems_2nd (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_brems_2nd_thl (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_brems_2nd_1 (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_brems_2nd_2 (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_brems_2nd_thl_1 (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_brems_2nd_thl_2 (const double xx[], double ff[], const double weight[])const;
	//	Integrands for the shift in Q2
	int integrand_shiftQ2 (const double xx[], double ff[], const double weight[])const;

//	Integrands for Carbon-12 target
	int integrand_born_carbon (const double xx[], double ff[], const double weight[])const;
	int integrand_elastic_carbon (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_born_carbon (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_elastic_carbon (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_carbon (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_brems_1st_carbon (const double xx[], double ff[], const double weight[])const;

	int integrand_brems_1st_carbon_thl (const double xx[], double ff[], const double weight[])const;
	int integrand_interf_brems_1st_carbon_thl (const double xx[], double ff[], const double weight[])const;

	//	Integrands for electron target
	int integrand_born_e (const double xx[], double ff[], const double weight[])const;
	int integrand_brems_1st_e (const double xx[], double ff[], const double weight[])const;

	static int cuba_integrand_born (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_born(xx,ff,weight);
	}

	static int cuba_integrand_elastic (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_elastic(xx,ff,weight);
	}

	static int cuba_integrand_interf_born (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_born(xx,ff,weight);
	}

	static int cuba_integrand_interf_elastic (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_elastic(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_test (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_test(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_hadr (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_hadr(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_hadr_interf (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_hadr_interf(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_l1k (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_l1k(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_l2k (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_l2k(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_2diff_v1 (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_2diff_v1(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_2diff_v2 (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_2diff_v2(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_ps2_2diff (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_ps2_2diff(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_ps2_3diff (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_ps2_3diff(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_ps2 (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_ps2(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_sg_diff (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_sg_diff(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l1k1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l1k1(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l1k2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l1k2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l2k1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l2k1(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l2k2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l2k2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l1k1_1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l1k1_1(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l1k2_1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l1k2_1(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l2k1_1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l2k1_1(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l2k2_1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l2k2_1(xx,ff,weight);
	}


	static int cuba_integrand_brems_2nd_add_l1k1_2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l1k1_2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l1k2_2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l1k2_2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l2k1_2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l2k1_2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l2k2_2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l2k2_2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l1k1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l1k1(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l1k2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l1k2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l2k1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l2k1(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l2k2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l2k2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l1k1_1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l1k1_1(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l1k2_1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l1k2_1(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l2k1_1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l2k1_1(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l2k2_1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l2k2_1(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l1k1_2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l1k1_2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l1k2_2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l1k2_2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l2k1_2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l2k1_2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l2k2_2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l2k2_2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_thl (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_thl(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_thl_1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_thl_1(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_thl_2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_thl_2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_thl (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_thl(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_thl_1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_thl_1(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_thl_2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_thl_2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_2diff (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_2diff(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_1diff (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_1diff(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_2diff (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_2diff(xx,ff,weight);
	}


	static int cuba_integrand_brems_2nd_add_l1k1_2diff (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l1k1_2diff(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l1k2_2diff (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l1k2_2diff(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l2k1_2diff (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l2k1_2diff(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_add_l2k2_2diff (const int *ndim, const double xx[],
				   const int *ncomp, double ff[], void *userdata,
				   const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_add_l2k2_2diff(xx,ff,weight);
	}


	static int cuba_integrand_brems_2nd_l1k1_2diff (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l1k1_2diff(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l1k2_2diff (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l1k2_2diff(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l2k1_2diff (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l2k1_2diff(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_l2k2_2diff (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_l2k2_2diff(xx,ff,weight);
	}

	static int cuba_integrand_interf_brems_1st (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_brems_1st(xx,ff,weight);
	}

	static int cuba_integrand_interf_brems_2nd (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_brems_2nd(xx,ff,weight);
	}

	static int cuba_integrand_interf_brems_2nd_1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_brems_2nd_1(xx,ff,weight);
	}

	static int cuba_integrand_interf_brems_2nd_2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_brems_2nd_2(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_pol_add (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_pol_add(xx,ff,weight);
	}

	static int cuba_integrand_brems_2nd_pol_add_thl (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_2nd_pol_add_thl(xx,ff,weight);
	}

	static int cuba_integrand_interf_brems_2nd_thl (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_brems_2nd_thl(xx,ff,weight);
	}

	static int cuba_integrand_interf_brems_2nd_thl_1 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_brems_2nd_thl_1(xx,ff,weight);
	}

	static int cuba_integrand_interf_brems_2nd_thl_2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_brems_2nd_thl_2(xx,ff,weight);
	}

	static int cuba_integrand_shiftQ2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_shiftQ2(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_thl (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_thl(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_thl_ps2 (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_thl_ps2(xx,ff,weight);
	}

	static int cuba_integrand_interf_brems_1st_thl (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_brems_1st_thl(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_hadr_thl (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_hadr_thl(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_hadr_interf_thl (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_hadr_interf_thl(xx,ff,weight);
	}

	static int cuba_integrand_born_carbon (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_born_carbon(xx,ff,weight);
	}

	static int cuba_integrand_elastic_carbon (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_elastic_carbon(xx,ff,weight);
	}

	static int cuba_integrand_interf_born_carbon (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_born_carbon(xx,ff,weight);
	}

	static int cuba_integrand_interf_elastic_carbon (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_elastic_carbon(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_carbon (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_carbon(xx,ff,weight);
	}

	static int cuba_integrand_interf_brems_1st_carbon (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_brems_1st_carbon(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_carbon_thl (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_carbon_thl(xx,ff,weight);
	}

	static int cuba_integrand_interf_brems_1st_carbon_thl (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_interf_brems_1st_carbon_thl(xx,ff,weight);
	}

	static int cuba_integrand_born_e (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_born_e(xx,ff,weight);
	}

	static int cuba_integrand_brems_1st_e (const int *ndim, const double xx[],
			const int *ncomp, double ff[], void *userdata,
			const int *nvec, const int *core, const double weight[], const int *iter) {
		const Integrands* integrand (static_cast<const Integrands*>(userdata));
		return integrand->integrand_brems_1st_e(xx,ff,weight);
	}

	Integrands ();

};

}

#endif /* SRC_INTEGRANDS_H_ */
