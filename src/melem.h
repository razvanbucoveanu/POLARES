/*
 * melem.h
 *
 *  Created on: Apr 18, 2016
 *      Author: razvan
 */

#ifndef MELEM_H_
#define MELEM_H_

#include "parameters.h"

namespace POLARES {

class Melem {

protected:
	const Parameters* param;
	mutable long double m, m2, m4, m6;
	mutable long double M, M2, M4, alpha;

public:
	int set_param(const Parameters* param);
	long double melem2(const long double l1k1, const long double l1k2, const long double k1k2, const long double Q2e,
			const long double Q2h, const long double Q2k, const long double S, const long double Sk, const long double Sq2,
			const long double f1, const long double f2)const;
	long double melem2_R2(const long double p1l1, const long double p1l2, const long double p1k1, const long double p1k2,
			const long double p1p2, const long double l1p2, const long double l2p2, const long double k1p2, const long double p2k2,
			const long double l1l2, const long double l1k1, const long double l2k1, const long double l1k2, const long double l2k2,
			const long double k1k2, const long double eg1, const long double f1, const long double f2)const;
	long double melem2_Rinterf(const long double l1k1, const long double l1k2, const long double l2k1, const long double l2k2,
			const long double Q2e, const long double Q2h, const long double S, const long double Sk, const long double Sq2,
			const long double f1, const long double f2)const;
	long double melem2_test(const long double l1k1, const long double l1k2, const long double k1k2, const long double Q2e,
			const long double Q2h, const long double Q2k, const long double S, const long double Sk, const long double Sq2,
			const long double f1, const long double f2)const;
	long double melem2_l1k1_1(const long double l1k1, const long double l1k2, const long double k1k2, const long double Q2e,
			const long double Q2h, const long double Q2k, const long double S, const long double Sk, const long double Sq2,
			const long double f1, const long double f2)const;
	long double melem2_l1k2_1(const long double l1k1, const long double l1k2, const long double k1k2, const long double Q2e,
			const long double Q2h, const long double Q2k, const long double S, const long double Sk, const long double Sq2,
			const long double f1, const long double f2)const;
	long double melem2_l2k1_1(const long double l1k1, const long double l1k2, const long double k1k2, const long double Q2e,
			const long double Q2h, const long double Q2k, const long double S, const long double Sk, const long double Sq2,
			const long double f1, const long double f2)const;
	long double melem2_l2k2_1(const long double l1k1, const long double l1k2, const long double k1k2, const long double Q2e,
			const long double Q2h, const long double Q2k, const long double S, const long double Sk, const long double Sq2,
			const long double f1, const long double f2)const;
	long double melem2_l1k1_2(const long double l1k1, const long double l1k2, const long double k1k2, const long double Q2e,
			const long double Q2h, const long double Q2k, const long double S, const long double Sk, const long double Sq2,
			const long double f1, const long double f2)const;
	long double melem2_l1k2_2(const long double l1k1, const long double l1k2, const long double k1k2, const long double Q2e,
			const long double Q2h, const long double Q2k, const long double S, const long double Sk, const long double Sq2,
			const long double f1, const long double f2)const;
	long double melem2_l2k1_2(const long double l1k1, const long double l1k2, const long double k1k2, const long double Q2e,
			const long double Q2h, const long double Q2k, const long double S, const long double Sk, const long double Sq2,
			const long double f1, const long double f2)const;
	long double melem2_l2k2_2(const long double l1k1, const long double l1k2, const long double k1k2, const long double Q2e,
			const long double Q2h, const long double Q2k, const long double S, const long double Sk, const long double Sq2,
			const long double f1, const long double f2)const;

	Melem();
};

}

#endif /* MELEM_H_ */
