/*
 * melem_pol.h
 *
 *  Created on: May 23, 2019
 *      Author: razvan
 */

#ifndef SRC_MELEM_POL_H_
#define SRC_MELEM_POL_H_


#include "parameters.h"

namespace POLARES {

class Melem_pol {

protected:
	const Parameters* param;
	mutable double m, m2, m4, m6;
	mutable double M, M2, M4, alpha;

public:
	int set_param(const Parameters* param);
	double melem_interf(const double l1k1, const double l1k2, const double k1k2, const double Q2e,
			const double Q2h, const double Q2k, const double S, const double Sk, const double Sq2,
			const double f1, const double f2, const double f1z, const double f2z, const double gae,
			const double ga, const double gv)const;
	double melem2_pol_add(const double l1k1, const double l1k2, const double l2k1, const double l2k2,
			const double Q2e, const double Q2h, const double S, const double Sk, const double Sq2,
			const double f1, const double f2, const double f1z, const double f2z, const double gae,
			const double ga, const double gv)const;
	double melem2_pol_add_interf_g(const double l1k1, const double l1k2, const double l2k1, const double l2k2,
			const double Q2e, const double Q2h, const double S, const double Sk, const double Sq2,
			const double f1, const double f2, const double f1z, const double f2z, const double gae,
			const double ga, const double gv)const;
	double melem2_pol_add_interf_Z(const double l1k1, const double l1k2, const double l2k1, const double l2k2,
			const double Q2e, const double Q2h, const double S, const double Sk, const double Sq2,
			const double f1, const double f2, const double f1z, const double f2z, const double gae,
			const double ga, const double gv)const;

	Melem_pol();
};

}


#endif /* SRC_MELEM_POL_H_ */
