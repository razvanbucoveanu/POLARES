/*
 * scalar_integrals.h
 *
 *  Created on: Mar 27, 2017
 *      Author: razvan
 */

#ifndef SCALAR_INTEGRALS_H_
#define SCALAR_INTEGRALS_H_

#include "config.h"

#ifdef POLARES_USE_LOOPTOOLS
# include "looptools_interface.h"
#endif

namespace POLARES {

class Parameters;

class Scalar_Integrals {

private:
	const Parameters* param;
	mutable long double Delta_eps, mu_dim2, lambda2;
	mutable long double x_Q2e, x_Q2h, beta_Q2e, beta_Q2h;
#ifdef POLARES_USE_LOOPTOOLS
    LoopTools LT;
#endif
public:
//	Scalar_Integrals();

	/**
	 * 1-point function
	 * @return A0(m2)
	 */
	long double A0_m(const long double m2)const;
	/**
	 * 2-point function
	 * @return B0(0,m2,m2)
	 */
	long double B0_0mm(const long double m2)const;
	/**
	 * 2-point function
	 * @return B0(m2,0,m2)
	 */
	long double B0_m0m(const long double m2)const;
	/**
	 * 2-point function
	 * @return B0(0,0,m2)
	 */
	long double B0_00m(const long double m2)const;
	/**
	 * 2-point function
	 * @return B0(-Q2e,m2,m2)
	 */
	long double B0_qmm(const long double Q2e, const long double m2)const;
	/**
	 * 2-point function
	 * @return B0(M2,0,m2)
	 */
	long double B0_M0m(const long double sll, const long double m2)const;
	/**
	 * IR divergent 3-point function
	 * regularized with a photon mass lambda
	 * @return C0(m2,m2,-Q2e,m2,0,m2)
	 */
	long double C0_mmqm0m(const long double Q2e, const long double m2)const;
	/**
	 * IR divergent 3-point function
	 * regularized with a photon mass lambda
	 * @return C0(m2,m2,0,m2,0,m2)
	 */
	long double C0_mm0m0m(const long double m)const;
	/**
	 * 3-point function
	 * @return C0(m2,0,M2,0,m2,m2)
	 */
	long double C0_m0M0mm(const long double sll, const long double m2)const;
	/**
	 * 3-point function
	 * @return C0(0,-Q2e,-Q2h,m2,m2,m2)
	 */
	long double C0_0qQmmm(const long double Q2e, const long double Q2h, const long double m2)const;
	/**
	 * 3-point function calculated by LoopTools
	 * to be used only if LoopTools is included
	 * @return C0(m2,M2,-Q2h,m2,0,m2)
	 */
	long double C0_mMQm0m(const long double sll, const long double Q2h, const long double m2)const;
	/**
	 * IR divergent 4-point function calculated by LoopTools
	 * to be used only if LoopTools is included
	 * regularized with a photon mass lambda
	 * @return D0(m2,m2,0,-Q2h,-Q2e,M2,m2,0,m2,m2)
	 */
	long double D0_mm0QqMm0mm(const long double Q2h, const long double Q2e, const long double sll, const long double m2)const;

	int set_param(const Parameters* param);
};

}

#endif /* SCALAR_INTEGRALS_H_ */
