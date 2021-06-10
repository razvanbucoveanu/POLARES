/*
 * virtual_corrections.h
 *
 *  Created on: Dec 1, 2015
 *      Author: razvan
 */

#ifndef SRC_VIRTUAL_CORRECTIONS_H_
#define SRC_VIRTUAL_CORRECTIONS_H_

#include "scalar_integrals.h"
#include "config.h"

#include "gamma_loop.h"

namespace POLARES{

class Parameters;
class Interpolation;

class Virtual_Corrections {

private:
	const Interpolation* interpolation;
	const Parameters* param;
//	double en, delta;
	mutable double v, S, f1e, f2e;
	mutable double p1p2, G1, G2, G3;
	mutable double m, m2, m4;
	mutable double M, M2, pi, pi2;
	mutable double m_me, m_mu, m_tau;
	mutable double alpha, sw2, mw2, mz2;
	//leptonic contribution
	double Jkapflkl(const double T3l, const double Ql, const double Q2, const double ml)const;
	double Jkapflgr(const double T3l, const double Ql, const double Q2, const double ml)const;
	//light quarks contribution
	double Jkapfqkl(const double T3q, const double Qq, const double Q2, const double mq)const;
	double Jkapfqgr(const double T3q, const double Qq, const double Q2, const double mq)const;
	double f1e_1st(const double Q2)const;
	double f2e_1st(const double Q2)const;
	double f1e_2nd(const double Q2)const;

	Scalar_Integrals SI;
public:
	Gamma_Loop GL;
	//vacuum polarization
	double d_vac_1st(const double Q2, const double m_lepton)const;
	double d_vac_2nd(const double Q2, const double m_lepton)const;
	double d_box_Feshbach(const double Q2)const;
	double d_box_MTj(const double Q2)const;

	//vertex correction
	double d_vert(const double Q2, const double f1, const double f2)const;
	double d_vert_pol(const double Q2, const double f1, const double f2,
			const double f1z, const double f2z, const double gae, const double ga, const double gv)const;
	double d_vert_pol_g(const double Q2, const double f1, const double f2,
			const double f1z, const double f2z, const double gae, const double ga, const double gv)const;
	double d_vert_pol_Z(const double Q2, const double f1, const double f2,
			const double f1z, const double f2z, const double gae, const double ga, const double gv)const;
	double d_vert_pol_quad_gZ(const double Q2, const double f1, const double f2,
			const double f1z, const double f2z, const double gae, const double ga, const double gv)const;
	double d_vert_quad(const double Q2, const double f1, const double f2)const;
//	double d_vert_test(const double Q2)const;
	double d_vert_carbon(const double Q2)const;

	//soft-bremsstrahlung
	double d_brems_ee(const double Q2)const;
	double d_brems_ee(const double Q2, const double en1)const;
	double d_brems_hadr(const double Q2)const;
	double d_brems_hadr_interf(const double Q2)const;
	double d_brems_ee_test(const double Q2)const;

	//one hard-photon + one loop corrections
	long double d_gamma_loop(const long double Q2h, const long double Q2e, const long double l1k,
    		const long double S, const long double U, const long double f1, const long double f2)const;

	long double d_gamma_loop_pol(const long double Q2h, const long double Q2e, const long double l1k,
    		const long double S1, const long double S2, const long double f1, const long double f2)const;

	//total 2nd order correction
	double d_2nd_total(const double Q2, const double f1, const double f2)const;

	double d_2nd_pol_total(const double Q2, const double f1, const double f2,
			const double f1z, const double f2z, const double gae, const double ga, const double gv)const;

	/**
	 * weak corrections
	 * arXiv : 1107.4683, formulas extracted from Jegerlehner's code alphaQED
	 * using effective quark masses, A. Weber, H. Spiesberger
	 */
	double kappa_weak(const double Q2)const;
	double running_alpha(const double Q2)const;

	int set_param(const Parameters* param, const Interpolation* interpolation);

	Virtual_Corrections();

};

} //namespace POLARES

#endif /* SRC_VIRTUAL_CORRECTIONS_H_ */
