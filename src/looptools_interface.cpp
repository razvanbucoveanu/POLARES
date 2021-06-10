/*
 * looptools.cpp
 *
 *  Created on: Mar 10, 2017
 *      Author: razvan
 */

#include "config.h"

#ifdef POLARES_USE_LOOPTOOLS

#include <clooptools.h>
#include <iostream>

#include "looptools_interface.h"
#include "const.h"

using namespace constants;
using namespace POLARES;

LoopTools::LoopTools() {
	ltini();
	setversionkey(1*KeyC0);
	setversionkey(1*KeyD0);
//	setlambda(me2);
//	setdelta(0.);
//	setmudim(1.);
}

LoopTools::~LoopTools() {
	clearcache();
//	ltexi();
}

double LoopTools::C0_looptools(const double p1, const double p2, const double p3,
			const double m1, const double m2, double const m3)const {

	clearcache();

	return Re(C0(p1,p2,p3,m1,m2,m3));
}

double LoopTools::D0_looptools(const double p1, const double p2, const double p3,
			const double p4, const double Q12, const double Q23, const double m1,
			const double m2, double const m3, const double m4)const {

	clearcache();

	return Re(D0(p1,p2,p3,p4,Q12,Q23,m1,m2,m3,m4));
}

void LoopTools::set_Delta_eps(const double Delta_eps) {

	setdelta(Delta_eps);
}

void LoopTools::set_mudim2(const double mu_dim2) {

	setmudim(mu_dim2);
}

void LoopTools::set_lambda2(const double lambda2) {

	setlambda(lambda2);
}

#endif /* POLARES_USE_LOOPTOOLS */
