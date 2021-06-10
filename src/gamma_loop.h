/*
 * gamma_loop.h
 *
 *  Created on: Mar 28, 2017
 *      Author: razvan
 */

#ifndef GAMMA_LOOP_H_
#define GAMMA_LOOP_H_

#include "scalar_integrals.h"

namespace POLARES {

class Gamma_Loop {

private:
	mutable long double deltaZ1, deltam, lambda2;
	mutable long double melem_interf;
	const Parameters* param;
	mutable long double m, m2, m4, m6, m8, m10, m12, m14;
	mutable long double M, M2, c;

public:
	long double M1gamma(const long double Q2h, const long double Q2e, const long double l1k,
			const long double S, const long double U, const long double G1, const long double G23)const;
    Scalar_Integrals SI;
    /**
     * Self-energy with a hard-photon emitted from the off-shell line
     */
    long double se_off_shell(const long double Q2h, const long double Q2e, const long double l1k,
    		const long double S, const long double U, const long double G1, const long double G23)const;
    /**
     * Self-energy with a hard-photon emitted from the on-shell line
     */
    long double se_on_shell(const long double Q2h, const long double Q2e, const long double l1k,
    		const long double S, const long double U, const long double G1, const long double G23)const;
    /**
     * Vertex correction with a hard-photon emitted from the off-shell line
     */
    long double vb_off_shell(const long double Q2h, const long double Q2e, const long double l1k,
    		const long double S, const long double U, const long double G1, const long double G23)const;

    /**
     * Vertex correction with a hard-photon emitted from the on-shell line
     */
    long double vb_on_shell(const long double Q2h, const long double Q2e, const long double l1k,
    		const long double S, const long double U, const long double G1, const long double G23)const;

    Gamma_Loop();
	int set_param(const Parameters* param);

};

}

#endif /* GAMMA_LOOP_H_ */
