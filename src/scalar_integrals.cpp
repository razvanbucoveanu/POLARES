/*
 * scalar_integrals.cpp
 *
 *  Created on: Mar 27, 2017
 *      Author: razvan
 */

#include "scalar_integrals.h"
#include "parameters.h"
#include <cmath>
#include "const.h"
#include <gsl/gsl_sf_dilog.h>
#include <stdlib.h>
#include <iostream>

using namespace constants;
using namespace POLARES;

//Scalar_Integrals::Scalar_Integrals()
//: Delta_eps(0.), mu_dim2(me2), lambda2(me2),
//  x_Q2e(0.), x_Q2h(0.), beta_Q2e(0.), beta_Q2h(0.), param(0) {
//}

int Scalar_Integrals::set_param(const Parameters* param) {
	this->param = param;
	Delta_eps = param->Delta_eps;
	mu_dim2 = param->mu_dim;
	lambda2 = param->lambda;
#ifdef POLARES_USE_LOOPTOOLS
	LT.set_Delta_eps(Delta_eps);
	LT.set_mudim2(mu_dim2);
	LT.set_lambda2(lambda2);
#endif

	return 0;
}

long double Scalar_Integrals::A0_m(const long double m2)const {

	return m2*(Delta_eps + 1. - log(m2/mu_dim2));
}

long double Scalar_Integrals::B0_0mm(const long double m2)const {

	return Delta_eps - log(m2/mu_dim2);
}

long double Scalar_Integrals::B0_m0m(const long double m2)const {

	return Delta_eps + 2. - log(m2/mu_dim2);
}

long double Scalar_Integrals::B0_00m(const long double m2)const {

	return Delta_eps + 1. - log(m2/mu_dim2);
}

long double Scalar_Integrals::B0_qmm(const long double Q2, const long double m2)const {

	beta_Q2e = sqrt(1. + (4.*m2)/Q2);
	x_Q2e = (beta_Q2e - 1.)/(beta_Q2e + 1.);

	return Delta_eps + 2. - beta_Q2e*log(1./x_Q2e) - log(m2/mu_dim2);
}

long double Scalar_Integrals::B0_M0m(const long double sll, const long double m2)const {

	long double res;
	if (sll > m2) res = Delta_eps + 2. + (m2 - sll)*log(sll/m2-1.)/sll - log(m2/mu_dim2);
	else if (sll < m2) res = Delta_eps + 2. + (m2 - sll)*log(1.-sll/m2)/sll - log(m2/mu_dim2);
	return res;
}

long double Scalar_Integrals::C0_mmqm0m(const long double Q2e, const long double m2)const {

	beta_Q2e = sqrt(1. + (4.*m2)/Q2e);
	x_Q2e = (beta_Q2e - 1.)/(beta_Q2e + 1.);

	return -1./Q2e/beta_Q2e*(log(lambda2/m2)*log(x_Q2e) - 2.*gsl_sf_dilog(-x_Q2e)
			- 2.*log(x_Q2e)*log(1.+x_Q2e) + 1./2.*pow(log(x_Q2e),2.) - pi*pi/6.);
}

long double Scalar_Integrals::C0_mm0m0m(const long double m2)const {

	return log(lambda2/m2)/2./m2;
}

long double Scalar_Integrals::C0_m0M0mm(const long double sll, const long double m2)const {

	return -1./(m2 - sll)*(pow(pi,2.)/6. - gsl_sf_dilog(sll/m2));
}

long double Scalar_Integrals::C0_0qQmmm(const long double Q2e, const long double Q2h, const long double m2)const {

	beta_Q2e = sqrt(1. + (4.*m2)/Q2e);
	beta_Q2h = sqrt(1. + (4.*m2)/Q2h);
	x_Q2e = (beta_Q2e - 1.)/(beta_Q2e + 1.);
	x_Q2h = (beta_Q2h - 1.)/(beta_Q2h + 1.);

	return 1./2. * 1./(Q2h - Q2e) * (pow(log(x_Q2e),2.) - pow(log(x_Q2h),2.));
}

long double Scalar_Integrals::C0_mMQm0m(const long double sll, const long double Q2h,
		const long double m2)const {

	long double t = -Q2h;

	long double lambda = sqrt(pow(m2 - t,2) - 2*(m2 + t)*sll + pow(sll,2));
	long double beta = sqrt(1 - (4*m2)/t);
	long double result;

	if (sll < 0)
	result = pow(lambda,-1)*(log(2)*log(2*(lambda - m2 + t - sll)*pow(lambda + m2 - t - sll,-1)*
        pow(sll,2)*pow(lambda + m2 - t + sll,-2)) -
     log((m2 - sll)*pow(m2,-1))*log((lambda - m2 + t - sll)*
        pow(lambda + m2 - t + sll,-1)) +
     gsl_sf_dilog(((lambda + m2 - t + sll)*pow(sll,-1))/2.) +
     gsl_sf_dilog(((lambda - m2 + t + sll)*pow(sll,-1))/2.) -
     gsl_sf_dilog((lambda*m2 + m2*t + m2*sll - t*sll - beta*t*sll - pow(m2,2))*
       pow(2*m2 - t - beta*t,-1)*pow(sll,-1)) -
     gsl_sf_dilog((lambda*m2 + m2*t + m2*sll - t*sll + beta*t*sll - pow(m2,2))*
       pow(2*m2 - t + beta*t,-1)*pow(sll,-1)) -
     gsl_sf_dilog((2*m2 - t - beta*t)*sll*
       pow(-(lambda*m2) + m2*t + m2*sll - t*sll - beta*t*sll - pow(m2,2),-1)) -
     gsl_sf_dilog((2*m2 - t + beta*t)*sll*
       pow(-(lambda*m2) + m2*t + m2*sll - t*sll + beta*t*sll - pow(m2,2),-1)) +
     pow(log(-(sll*pow(lambda + m2 - t - sll,-1))),2)/2. -
     pow(log(-(sll*pow(lambda - m2 + t - sll,-1))),2)/2. -
     pow(log(-((lambda*m2 - m2*t - m2*sll + t*sll + beta*t*sll + pow(m2,2))*
           pow(2*m2 - t - beta*t,-1)*pow(sll,-1))),2)/2. -
     pow(log(-((lambda*m2 - m2*t - m2*sll + t*sll - beta*t*sll + pow(m2,2))*
           pow(2*m2 - t + beta*t,-1)*pow(sll,-1))),2)/2. -
     pow(log((lambda + m2 - t - sll)*pow(lambda + m2 - t + sll,-1)),2)/2. +
     pow(log(-(sll*pow(lambda + m2 - t + sll,-1))),2));

	if (sll > m2)
	result = (pow(lambda,-1)*(2*log(2)*log(-2*(lambda - m2 - sll + t)*pow(sll,2)*
          pow(lambda + m2 - sll - t,-1)*pow(lambda + m2 + sll - t,-2)) -
       2*log(-((m2 - sll)*pow(m2,-1)))*
        log(-((lambda - m2 - sll + t)*pow(lambda + m2 + sll - t,-1))) +
       2*gsl_sf_dilog(((lambda - m2 + sll + t)*pow(sll,-1))/2.) -
       2*gsl_sf_dilog(2*sll*pow(lambda + m2 + sll - t,-1)) +
       2*gsl_sf_dilog((-(lambda*m2) + m2*sll + m2*t - sll*t - beta*sll*t - pow(m2,2))*
          pow(sll,-1)*pow(2*m2 - t - beta*t,-1)) -
       2*gsl_sf_dilog((lambda*m2 + m2*sll + m2*t - sll*t - beta*sll*t - pow(m2,2))*
          pow(sll,-1)*pow(2*m2 - t - beta*t,-1)) -
       2*gsl_sf_dilog((lambda*m2 + m2*sll + m2*t - sll*t + beta*sll*t - pow(m2,2))*
          pow(sll,-1)*pow(2*m2 - t + beta*t,-1)) -
       2*gsl_sf_dilog(sll*(2*m2 - t + beta*t)*
          pow(-(lambda*m2) + m2*sll + m2*t - sll*t + beta*sll*t - pow(m2,2),-1)) -
       pow(log(((lambda + m2 + sll - t)*pow(sll,-1))/2.),2) +
       pow(log(sll*pow(lambda + m2 - sll - t,-1)),2) +
       2*pow(log(sll*pow(lambda + m2 + sll - t,-1)),2) -
       pow(log((lambda + m2 - sll - t)*pow(lambda + m2 + sll - t,-1)),2) -
       pow(log(-(sll*pow(lambda - m2 - sll + t,-1))),2) -
       pow(log(-((lambda*m2 - m2*sll - m2*t + sll*t + beta*sll*t + pow(m2,2))*
            pow(sll,-1)*pow(2*m2 - t - beta*t,-1))),2) -
       pow(log((lambda*m2 - m2*sll - m2*t + sll*t - beta*sll*t + pow(m2,2))*
          pow(sll,-1)*pow(2*m2 - t + beta*t,-1)),2) +
       pow(log(sll*(2*m2 - t - beta*t)*
          pow(-(lambda*m2) + m2*sll + m2*t - sll*t - beta*sll*t - pow(m2,2),-1)),2)))/2.;

	if (sll < m2 && sll > 0)
	result = (pow(lambda,-1)*(2*log(2)*log(-2*(lambda - m2 - sll + t)*pow(sll,2)*
      pow(lambda + m2 - sll - t,-1)*pow(lambda + m2 + sll - t,-2)) -
   2*log((m2 - sll)*pow(m2,-1))*
    log(-((lambda - m2 - sll + t)*pow(lambda + m2 + sll - t,-1))) +
   2*gsl_sf_dilog(((lambda - m2 + sll + t)*pow(sll,-1))/2.) -
   2*gsl_sf_dilog(2*sll*pow(lambda + m2 + sll - t,-1)) +
   2*gsl_sf_dilog((-(lambda*m2) + m2*sll + m2*t - sll*t - beta*sll*t - pow(m2,2))*
      pow(sll,-1)*pow(2*m2 - t - beta*t,-1)) -
   2*gsl_sf_dilog((lambda*m2 + m2*sll + m2*t - sll*t - beta*sll*t - pow(m2,2))*
      pow(sll,-1)*pow(2*m2 - t - beta*t,-1)) -
   2*gsl_sf_dilog((lambda*m2 + m2*sll + m2*t - sll*t + beta*sll*t - pow(m2,2))*
      pow(sll,-1)*pow(2*m2 - t + beta*t,-1)) -
   2*gsl_sf_dilog(sll*(2*m2 - t + beta*t)*
      pow(-(lambda*m2) + m2*sll + m2*t - sll*t + beta*sll*t - pow(m2,2),-1)) -
   pow(log(((lambda + m2 + sll - t)*pow(sll,-1))/2.),2) +
   pow(log(sll*pow(lambda + m2 - sll - t,-1)),2) +
   2*pow(log(sll*pow(lambda + m2 + sll - t,-1)),2) -
   pow(log((lambda + m2 - sll - t)*pow(lambda + m2 + sll - t,-1)),2) -
   pow(log(-(sll*pow(lambda - m2 - sll + t,-1))),2) -
   pow(log((lambda*m2 - m2*sll - m2*t + sll*t + beta*sll*t + pow(m2,2))*
      pow(sll,-1)*pow(2*m2 - t - beta*t,-1)),2) -
   pow(log((lambda*m2 - m2*sll - m2*t + sll*t - beta*sll*t + pow(m2,2))*
      pow(sll,-1)*pow(2*m2 - t + beta*t,-1)),2) +
   pow(log(-(sll*(2*m2 - t - beta*t)*
        pow(-(lambda*m2) + m2*sll + m2*t - sll*t - beta*sll*t - pow(m2,2),-1))),
    2)))/2.;

	return result;
}

long double Scalar_Integrals::D0_mm0QqMm0mm(const long double Q2h, const long double Q2e,
		const long double sll, const long double m2)const {

	long double t = -Q2h;
	long double ull = -Q2e;

	long double betat = sqrt(1 - (4*m2)/t);
	long double betau = sqrt(1 - (4*m2)/ull);
	long double result;

	if (sll > m2)
		result = -(betau*pow(m2 - sll,-1)*pow(4*m2 - ull,-1)*
      (pi2 + 12*log((-1 + betat)*pow(1 + betat,-1))*
         log((-1 + betau)*pow(1 + betau,-1)) +
        12*log(4*betau*pow(1 + betau,-2))*log((-1 + betau)*pow(1 + betau,-1)) -
        12*log((-1 + betat)*pow(1 + betat,-1))*
         log(2*(betat + betau)*pow(1 + betat,-1)*pow(1 + betau,-1)) -
        12*log((-1 + betau)*pow(1 + betau,-1))*
         log(2*(betat + betau)*pow(1 + betat,-1)*pow(1 + betau,-1)) +
        12*log((-1 + betat)*(-1 + betau)*pow(1 + betat,-1)*pow(1 + betau,-1))*
         log(2*(betat + betau)*pow(1 + betat,-1)*pow(1 + betau,-1)) +
        12*log(((-1 + betat)*(-1 + betau)*pow(betat + betau,-1))/2.)*
         log(((1 + betat)*(1 + betau)*pow(betat + betau,-1))/2.) +
        6*log((-1 + betau)*pow(1 + betau,-1))*log(m2*pow(lambda2,-1)) +
        12*log((-1 + betau)*pow(1 + betau,-1))*log(-((m2 - sll)*pow(m2,-1))) -
        12*gsl_sf_dilog(-2*(betat - betau)*pow(1 + betat,-1)*pow(-1 + betau,-1)) +
        6*gsl_sf_dilog(pow(-1 + betau,2)*pow(1 + betau,-2)) +
        12*gsl_sf_dilog(-((-1 + betat)*(-1 + betau)*pow(betat + betau,-1))/2.) -
        6*pow(log((-1 + betau)*pow(1 + betau,-1)),2) -
        6*pow(log(2*(betat + betau)*pow(1 + betat,-1)*pow(1 + betau,-1)),2)))/6.;

	if (sll < m2)
	result = -(betau*pow(m2 - sll,-1)*pow(4*m2 - ull,-1)*
      (pi2 + 12*log((-1 + betat)*pow(1 + betat,-1))*
         log((-1 + betau)*pow(1 + betau,-1)) +
        12*log(4*betau*pow(1 + betau,-2))*log((-1 + betau)*pow(1 + betau,-1)) -
        12*log((-1 + betat)*pow(1 + betat,-1))*
         log(2*(betat + betau)*pow(1 + betat,-1)*pow(1 + betau,-1)) -
        12*log((-1 + betau)*pow(1 + betau,-1))*
         log(2*(betat + betau)*pow(1 + betat,-1)*pow(1 + betau,-1)) +
        12*log((-1 + betat)*(-1 + betau)*pow(1 + betat,-1)*pow(1 + betau,-1))*
         log(2*(betat + betau)*pow(1 + betat,-1)*pow(1 + betau,-1)) -
        12*log(2*(betat + betau)*pow(1 + betat,-1)*pow(1 + betau,-1))*
         log(((-1 + betat)*(-1 + betau)*pow(betat + betau,-1))/2.) +
        6*log((-1 + betau)*pow(1 + betau,-1))*log(m2*pow(lambda2,-1)) +
        12*log((-1 + betau)*pow(1 + betau,-1))*log((m2 - sll)*pow(m2,-1)) -
        12*gsl_sf_dilog(-2*(betat - betau)*pow(1 + betat,-1)*pow(-1 + betau,-1)) +
        6*gsl_sf_dilog(pow(-1 + betau,2)*pow(1 + betau,-2)) +
        12*gsl_sf_dilog(-((-1 + betat)*(-1 + betau)*pow(betat + betau,-1))/2.) -
        6*pow(log((-1 + betau)*pow(1 + betau,-1)),2) -
        6*pow(log(2*(betat + betau)*pow(1 + betat,-1)*pow(1 + betau,-1)),2)))/6.;

	return result;
}
