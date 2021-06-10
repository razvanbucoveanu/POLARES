/*
 * looptools_interface.h
 *
 *  Created on: Mar 10, 2017
 *      Author: razvan
 */

#ifndef LOOPTOOLS_INTERFACE_H_
#define LOOPTOOLS_INTERFACE_H_

namespace POLARES {

class LoopTools {

private:

public:
	double C0_looptools(const double p1, const double p2, const double p3,
			const double m1, const double m2, double const m3)const;

	double D0_looptools(const double p1, const double p2, const double p3,
				const double p4, const double Q12, const double Q23, const double m1,
				const double m2, double const m3, const double m4)const;

	void set_Delta_eps(const double Delta_eps);
	void set_mudim2(const double mu_dim2);
	void set_lambda2(const double lambda2);

	LoopTools();
	~LoopTools();

};


}

#endif /* LOOPTOOLS_INTERFACE_H_ */
