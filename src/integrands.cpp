/*
 * integrands.cpp
 *
 *  Created on: Dec 4, 2015
 *      Author: razvan
 */

#include "const.h"
#include "integrands.h"
#include "interpolation.h"
#include "parameters.h"
#include <cmath>
#include <iostream>

//using namespace constants;
using namespace POLARES;

Integrands::Integrands ()
:a(0), b(0), c(0), d(0), f(0), n2(0), l1(0), l2(0), a1(0), a2(0), y(0), Jacobian_eg(0),Jacobian_eg1(0),
 thgmin(0), thgmax(0), egmin(0), egmax(0), thlmin(0), phig1min(0),phig1max(0),
 phig1min_cos(0),phig1max_cos(0), thg1min(0), thg1max(0), eg1min(0), eg1max(0),
 thl(0), thg(0), phig(0), en1(0), eg(0), eg1(0), thg1(0), phig1(0),
 Q2(0), Jacobian(0), max_wgt(0), Q2min(0), Q2max(0), events_no(0), acc_events(0),
 acc_events_1(0), acc_events_2(0), acc_events_3(0), acc_events_4(0), acc_events_5(0),
 event_count(false), param(NULL), interpolation(NULL), FS(NULL), m(constants::me),
 m2(constants::me2), M(constants::mpr), M2(constants::mpr2), pi(constants::pi) {}

void Integrands::set_param (const Parameters* param, const Interpolation* interpolation, Final_State* FS) {
	this->interpolation = interpolation;
	this->param = param;
	this->FS = FS;
	m = param->m;
	m2 = pow(m,2.);
	M = param->M;
	M2 = pow(M,2.);
	Q2min = param->min[param->Q2_elastic];
	Q2max = param->max[param->Q2_elastic];
	CS.set_param (param, interpolation);
//	VC.set_param(param);
	l1 = param->l1;
	std::cout.precision(6);
}

void Integrands::set_param (const Parameters* param, const Interpolation* interpolation) {
	this->interpolation = interpolation;
	this->param = param;
//	this->FS = FS;
	m = param->m;
	m2 = pow(m,2.);
	M = param->M;
	M2 = pow(M,2.);
	Q2min = param->min[param->Q2_elastic];
	Q2max = param->max[param->Q2_elastic];
	CS.set_param (param, interpolation);
//	VC.set_param(param);
	l1 = param->l1;
	std::cout.precision(6);
//	fout = fopen("events_test.txt", "w");
}

int Integrands::integrand_born (const double xx[], double ff[], const double weight[])const {

	Q2 = 1./Q2min+(1./Q2max-1./Q2min)*xx[0];
	Q2 = 1./Q2;
	Jacobian = - pow(Q2, 2.);
	ff[0] = CS.crsect_born(Q2) * (1./Q2max-1./Q2min) * Jacobian;

	return 0;
}

int Integrands::integrand_elastic (const double xx[], double ff[], const double weight[])const {

	Q2 = 1./Q2min+(1./Q2max-1./Q2min)*xx[0];
	Q2 = 1./Q2;
	Jacobian = - pow(Q2, 2);

	ff[0] = CS.crsect_elastic(Q2) * (1./Q2max-1./Q2min) * Jacobian;

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = (param->en - Q2 / (2.*M));
		FS->theta_l = acos((-Q2 * (1. + param->en/M) +
				2.*pow(param->en,2.) - 2.*m2)/(2.*l1*sqrt(pow(param->en-Q2/(2.*M),2.)-m2)));
		FS->E_gamma = 0.;
		FS->theta_gamma = 0.;
		FS->phi_gamma = 0.;
		FS->E_gamma_prime = 0.;
		FS->theta_gamma_prime = 0.;
		FS->phi_gamma_prime = 0.;
		FS->weight = weight[0]*weight[2];
		FS->sigma_diff = ff[0];
		if (param->flag[param->asymmetry] == 1)
			FS->weight *= 1. + param->P*CS.interf_elastic(Q2)/CS.crsect_elastic(Q2);


	}
	return 0;
}

int Integrands::integrand_interf_born (const double xx[], double ff[], const double weight[])const {

		Q2 = Q2min + (Q2max - Q2min) * xx[0];
		ff[0] = CS.interf_born(Q2) * (Q2max - Q2min);
		return 0;
	}

int Integrands::integrand_interf_elastic (const double xx[], double ff[], const double weight[])const {

	Q2 = Q2min+(Q2max-Q2min)*xx[0];

//Function to be integrated
	ff[0] = CS.interf_elastic(Q2)*(Q2max-Q2min);

	return 0;
}

// Integrand for Sigma | E_gamma > Delta
int Integrands::integrand_brems_1st (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	if (param->min[param->E_prime] >= param->en) return 0;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thlmin = (param->en*en1 - M*(param->en-en1) - m2)/(l1*l2);

	//check if the lower limit calculated is lower than the maximum limit given
	if (param->max[param->cos_thl] < thlmin) return 0;

	//if the lower limit calculated is lower than the given limit or -1 resets the lower limit
	if (thlmin < param->min[param->cos_thl]) thlmin = param->min[param->cos_thl];

	//generates number for cos(theta_l)
	thl = thlmin+(param->max[param->cos_thl]-thlmin)*xx[1];

	//check if -1<cos(thl)<1
	if (pow(thl, 2.) > 1.) return 0;

	thl = acos(thl);

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin >= thgmax) return 0;

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[3];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;

	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);
//	if (rand.uniform() < 0.5)
//	phig = 2.*pi - phig;

	if (std::abs(sin(thl)*sin(thg)*sin(phig)) < 1e-40) return 0;

//	long double test = m2 + M*(param->en-en1-eg) - param->en*(en1+eg) + l1*l2*cos(thl) +
//			l1*eg*cos(thg) + en1*eg - l2*eg*sin(thl)*sin(thg)*cos(phig) - l2*eg*cos(thl)*cos(thg);
//
//	if (std::abs(test) > 1e-16) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st(en1,thl,eg,thg,phig)*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-thlmin)*(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	if (param->flag[param->ps] == 1) {

		FS->weight = weight[0]*weight[5];
		FS->E_prime_l = en1;
		FS->theta_l = thl;
		FS->E_gamma = eg;
		FS->theta_gamma = thg;
		FS->phi_gamma = phig;
		FS->E_gamma_prime = 0.;
		FS->theta_gamma_prime = 0.;
		FS->phi_gamma_prime = 0.;
		FS->sigma_diff = ff[0];
		if (param->flag[param->asymmetry] == 1)
			FS->weight *= 1. + param->P*CS.interf_brems_1st(en1,thl,eg,thg)
									/CS.crsect_brems_1st(en1,thl,eg,thg,phig);
//		}
	};

	return 0;
}

// Integrand for Sigma | E_gamma > Delta
int Integrands::integrand_brems_1st_test (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	if (param->min[param->E_prime] >= param->en) return 0;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thlmin = (param->en*en1 - M*(param->en-en1) - m2)/(l1*l2);

	//check if the lower limit calculated is lower than the maximum limit given
	if (param->max[param->cos_thl] < thlmin) return 0;

	//if the lower limit calculated is lower than the given limit or -1 resets the lower limit
	if (thlmin < param->min[param->cos_thl]) thlmin = param->min[param->cos_thl];

	//generates number for cos(theta_l)
	thl = thlmin+(param->max[param->cos_thl]-thlmin)*xx[1];

	//check if -1<cos(thl)<1
	if (pow(thl, 2.) > 1.) return 0;

	thl = acos(thl);

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = -eg*l1*sin(thl);
	d = eg*(l1*cos(thl)-l2);
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

//	std::cout << thgmin <<"\t"<< thgmax <<"\n";

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin >= thgmax) return 0;

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[3];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)*cos(thl)+en1*eg-eg*l2*cos(thg)+m2;
	b = -eg*l1*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;

	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);
//	if (rand.uniform() < 0.5)
//	phig = 2.*pi - phig;

	if (std::abs(sin(thl)*sin(thg)*sin(phig)) < 1e-40) return 0;

//	long double test = m2 + M*(param->en-en1-eg) - param->en*(en1+eg) + l1*l2*cos(thl) -
//			l2*eg*cos(thg) + en1*eg + l1*eg*sin(thl)*sin(thg)*cos(phig) + l1*eg*cos(thl)*cos(thg);
//
//	std::cout << test << "\n";
//
//	if (std::abs(test) > 1e-16) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st_test(en1,thl,eg,thg,phig)*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-thlmin)*(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	return 0;
}

int Integrands::integrand_brems_1st_hadr (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	if (param->min[param->E_prime] >= param->en) return 0;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thlmin = (param->en*en1 - M*(param->en-en1) - m2)/(l1*l2);

	//check if the lower limit calculated is lower than the maximum limit given
	if (param->max[param->cos_thl] < thlmin) return 0;

	//if the lower limit calculated is lower than the given limit or -1 resets the lower limit
	if (thlmin < param->min[param->cos_thl]) thlmin = param->min[param->cos_thl];

	//generates number for cos(theta_l)
	thl = thlmin+(param->max[param->cos_thl]-thlmin)*xx[1];

	//check if -1<cos(thl)<1
	if (pow(thl, 2.) > 1.) return 0;

	thl = acos(thl);

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin >= thgmax) return 0;

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[3];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;

	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st_hadr(en1,thl,eg,thg,phig)*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-thlmin)*(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	if (param->flag[param->ps] == 1) {

		FS->weight = weight[0]*weight[5];
			FS->E_prime_l = en1;
			FS->theta_l = thl;
			FS->E_gamma = eg;
			FS->theta_gamma = thg;
			FS->phi_gamma = phig;
			FS->E_gamma_prime = 0.;
			FS->theta_gamma_prime = 0.;
			FS->phi_gamma_prime = 0.;
			FS->sigma_diff = ff[0];

//		}
	};

	return 0;
}

int Integrands::integrand_brems_1st_hadr_interf (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	if (param->min[param->E_prime] >= param->en) return 0;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thlmin = (param->en*en1 - M*(param->en-en1) - m2)/(l1*l2);

	//check if the lower limit calculated is lower than the maximum limit given
	if (param->max[param->cos_thl] < thlmin) return 0;

	//if the lower limit calculated is lower than the given limit or -1 resets the lower limit
	if (thlmin < param->min[param->cos_thl]) thlmin = param->min[param->cos_thl];

	//generates number for cos(theta_l)
	thl = thlmin+(param->max[param->cos_thl]-thlmin)*xx[1];

	//check if -1<cos(thl)<1
	if (pow(thl, 2.) > 1.) return 0;

	thl = acos(thl);

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin >= thgmax) return 0;

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[3];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;

	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st_hadr_interf(en1,thl,eg,thg,phig)*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-thlmin)*(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	if (param->flag[param->ps] == 1) {

		FS->weight = weight[0]*weight[5];
			FS->E_prime_l = en1;
			FS->theta_l = thl;
			FS->E_gamma = eg;
			FS->theta_gamma = thg;
			FS->phi_gamma = phig;
			FS->E_gamma_prime = 0.;
			FS->theta_gamma_prime = 0.;
			FS->phi_gamma_prime = 0.;
			FS->sigma_diff = ff[0];

//		}
	};

	return 0;
}



int Integrands::integrand_brems_1st_l1k (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	if (param->min[param->E_prime] >= param->en) return 0;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thlmin = (param->en*en1 - M*(param->en-en1) - m2)/(l1*l2);

	//check if the lower limit calculated is lower than the maximum limit given
	if (param->max[param->cos_thl] < thlmin) return 0;

	//if the lower limit calculated is lower than the given limit or -1 resets the lower limit
	if (thlmin < param->min[param->cos_thl]) thlmin = param->min[param->cos_thl];

	//generates number for cos(theta_l)
	thl = thlmin+(param->max[param->cos_thl]-thlmin)*xx[1];

	//check if -1<cos(thl)<1
	if (pow(thl, 2.) > 1.) return 0;

	thl = acos(thl);

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin >= thgmax) return 0;

	long double l1k_min = param->en * eg - l1 * eg * thgmax;
	long double l1k_max = param->en * eg - l1 * eg * thgmin;

	long double l1k = log(l1k_min) + (log(l1k_max) - log(l1k_min)) * xx[3];

	l1k = exp(l1k);
	long double Jacobian_thg = l1k / (l1 * eg);

	thg = (param->en * eg - l1k) / (l1 * eg);

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;

	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st(en1,thl,eg,thg,phig)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-thlmin)*(log(egmax)-log(egmin))*Jacobian_eg
			*(log(l1k_max)-log(l1k_min))*Jacobian_thg;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double l2k = en1*eg - eg*l2*cospsi;

	long double d1 = pow(l1k,2.);
	long double d2 = pow(l2k,2.);

	ff[0] *= d2/(d1+d2);

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = en1;
		FS->theta_l = thl;
		FS->E_gamma = eg;
		FS->theta_gamma = thg;
		FS->phi_gamma = phig;
		FS->E_gamma_prime = 0.;
		FS->theta_gamma_prime = 0.;
		FS->phi_gamma_prime = 0.;
		FS->weight = weight[0]*weight[5];

	};

	return 0;
}

int Integrands::integrand_brems_1st_l2k (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	if (param->min[param->E_prime] >= param->en) return 0;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thlmin = (param->en*en1 - M*(param->en-en1) - m2)/(l1*l2);

	//check if the lower limit calculated is lower than the maximum limit given
	if (param->max[param->cos_thl] < thlmin) return 0;

	//if the lower limit calculated is lower than the given limit or -1 resets the lower limit
	if (thlmin < param->min[param->cos_thl]) thlmin = param->min[param->cos_thl];

	//generates number for cos(theta_l)
	thl = thlmin+(param->max[param->cos_thl]-thlmin)*xx[1];

	//check if -1<cos(thl)<1
	if (pow(thl, 2.) > 1.) return 0;

	thl = acos(thl);

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;

	egmax = -y/a2;
//	std::cout << "test: " << param->min[param->E_gamma] <<"\t" << param->max[param->E_gamma] << "\n";
	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin >= thgmax) return 0;

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[3];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;


	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st(en1,thl,eg,thg,phig)*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-thlmin)*(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double l1k = param->en*eg - l1*eg*cos(thg);
	long double l2k = en1*eg - eg*l2*cospsi;

	long double d1 = pow(l1k,2.);
	long double d2 = pow(l2k,2.);

	ff[0] *= d1/(d1+d2);

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = en1;
		FS->theta_l = thl;
		FS->E_gamma = eg;
		FS->theta_gamma = thg;
		FS->phi_gamma = phig;
		FS->E_gamma_prime = 0.;
		FS->theta_gamma_prime = 0.;
		FS->phi_gamma_prime = 0.;
		FS->weight = weight[0]*weight[5];

	};

	return 0;
}

// Integrand for Sigma | E_gamma > Delta
int Integrands::integrand_brems_1st_2diff_v1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	//check if cos(theta_gamma) is inside the cut-off limits
	if (cos(param->thg) < thgmin || cos(param->thg) > thgmax) return 0;

	thg = param->thg;

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;

	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);

	std::cout.precision(14);
//	std::cout << en1 <<"\t"<< thl*180./pi <<"\t"<< eg <<"\t"<< thg*180./pi <<"\n";

	//function to be calculated
	ff[0] = CS.crsect_brems_1st(en1,thl,eg,thg,phig)*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg;

	return 0;
}

// Integrand for Sigma | E_gamma > Delta
int Integrands::integrand_brems_1st_2diff_v2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);


	egmin = y/a1;

	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[0];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;


	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[1];

	//check if cos(theta_gamma) is inside the cut-off limits
	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;


	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);

	if (acos(cospsi) * 180./pi <= 20.) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st(en1,thl,eg,thg,phig)*(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	return 0;
}

// Integrand for Sigma | E_gamma > Delta
// Second phase space parametrization
int Integrands::integrand_brems_1st_ps2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

	thg = thgmin+(thgmax-thgmin)*xx[0];
	thg = acos(thg);
//	long double Jacobian_thg = sin(thg);

	long double phigmin = 0.;
	long double phigmax = 2.*pi;

	phig = phigmin+(phigmax-phigmin)*xx[1];

	egmin = param->min[param->E_gamma];
	egmax = param->max[param->E_gamma];

	if (egmax > M*(param->en-m)/(M+param->en-l1*cos(thg)))
		egmax = M*(param->en-m)/(M+param->en-l1*cos(thg));
	if (egmin > egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	thlmin = param->min[param->cos_thl];
	long double thlmax = param->max[param->cos_thl];

	//generates number for cos(theta_l)
	thl = thlmin+(thlmax-thlmin)*xx[3];
	thl = acos(thl);
//	long double Jacobian_thl = sin(thl);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double A = l1*cos(thl)-eg*cospsi;
	long double B = param->en+M-eg;
	long double C = eg*(param->en+M-l1*cos(thg))-M*param->en-m2;

	if (m2*(pow(A,2.)-pow(B,2.))+pow(C,2.) <= 0.) return 0;

	if (std::abs(pow(A,2.)-pow(B,2.)) < 1.e-12) return 0;

	en1 = (B*C - A*sqrt(m2*(pow(A,2.)-pow(B,2.))+pow(C,2.)))/(pow(A,2.)-pow(B,2.));

	if (en1 < param->min[param->E_prime] || en1 > param->en - eg) return 0;

	l2 = sqrt(pow(en1,2.)-m2);

	long double test = m2 + M*(param->en-en1-eg) - param->en*(en1+eg) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*eg - l2*eg*sin(thl)*sin(thg)*cos(phig) - l2*eg*cos(thl)*cos(thg);

	if (std::abs(test) > 1e-16) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st_ps2(en1,thl,eg,thg,phig)*Jacobian_eg
			*(thgmax-thgmin)*(phigmax-phigmin)*(log(egmax)-log(egmin))*(thlmax-thlmin);

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = en1;
		FS->theta_l = thl;
		FS->E_gamma = eg;
		FS->theta_gamma = thg;
		FS->phi_gamma = phig;
		FS->E_gamma_prime = 0.;
		FS->theta_gamma_prime = 0.;
		FS->phi_gamma_prime = 0.;
		FS->weight = weight[0]*weight[5];
		FS->sigma_diff = ff[0];
		if (param->flag[param->asymmetry] == 1)
			FS->weight *= 1. + param->P*CS.interf_brems_1st_ps2(en1,thl,eg,thg,phig)
									/CS.crsect_brems_1st_ps2(en1,thl,eg,thg,phig);

	};

	return 0;
}

int Integrands::integrand_brems_1st_sg_diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

	thg = thgmin+(thgmax-thgmin)*xx[0];
	thg = acos(thg);

	long double phigmin = 0.;
	long double phigmax = 2.*pi;

	phig = phigmin+(phigmax-phigmin)*xx[1];

	egmin = 0.;
	egmax = param->min[param->E_gamma];

	//generates photon energy
	eg = egmax*xx[2];

	thlmin = param->min[param->cos_thl];
	long double thlmax = param->max[param->cos_thl];

	//generates number for cos(theta_l)
	thl = thlmin+(thlmax-thlmin)*xx[3];
	thl = acos(thl);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double A = l1*cos(thl)-eg*cospsi;
	long double B = param->en+M-eg;
	long double C = eg*(param->en+M-l1*cos(thg))-M*param->en-m2;

	if (m2*(pow(A,2.)-pow(B,2.))+pow(C,2.) <= 0.) return 0;

	if (std::abs(pow(A,2.)-pow(B,2.)) < 1.e-12) return 0;

	en1 = (B*C - A*sqrt(m2*(pow(A,2.)-pow(B,2.))+pow(C,2.)))/(pow(A,2.)-pow(B,2.));

	if (en1 < param->min[param->E_prime] || en1 > param->en - eg) return 0;

	l2 = sqrt(pow(en1,2.)-m2);

	long double test = m2 + M*(param->en-en1-eg) - param->en*(en1+eg) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*eg - l2*eg*sin(thl)*sin(thg)*cos(phig) - l2*eg*cos(thl)*cos(thg);

	if (std::abs(test) > 1e-16) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st_sg_diff(en1,thl,eg,thg,phig)
			*(thgmax-thgmin)*(phigmax-phigmin)*egmax*(thlmax-thlmin);

	return 0;
}

// Integrand for Sigma | E_gamma > Delta
int Integrands::integrand_brems_1st_ps2_2diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	thg = param->thg;

	long double phigmin = 0.;
	long double phigmax = 2.*pi;

	phig = phigmin+(phigmax-phigmin)*xx[0];

	egmin = param->min[param->E_gamma];
	egmax = param->max[param->E_gamma];

	//generates photon energy
	eg = egmin+(egmax-egmin)*xx[1];
	thl = param->thl;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double A = l1*cos(thl)-eg*cospsi;
	long double B = param->en+M-eg;
	long double C = eg*(param->en+M-l1*cos(thg))-M*param->en-m2;

	if (m2*(pow(A,2.)-pow(B,2.))+pow(C,2.) <= 0.) return 0;

	if (std::abs(pow(A,2.)-pow(B,2.)) < 1.e-12) return 0;

	en1 = (B*C - A*sqrt(m2*(pow(A,2.)-pow(B,2.))+pow(C,2.)))/(pow(A,2.)-pow(B,2.));

	if (en1 < m || en1 > param->en - eg) return 0;

	l2 = sqrt(pow(en1,2.)-m2);

	long double test = m2 + M*(param->en-en1-eg) - param->en*(en1+eg) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*eg - l2*eg*sin(thl)*sin(thg)*cos(phig) - l2*eg*cos(thl)*cos(thg);

	if (std::abs(test) > 1e-16) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st_ps2(en1,thl,eg,thg,phig)*(phigmax-phigmin)*(egmax-egmin);

	return 0;
}

int Integrands::integrand_brems_1st_ps2_3diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	thg = param->thg;
	phig = 0.;

	egmin = param->min[param->E_gamma];
	egmax = param->max[param->E_gamma];

	//generates photon energy
	eg = egmin+(egmax-egmin)*xx[0];


	thl = param->thl;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double A = l1*cos(thl)-eg*cospsi;
	long double B = param->en+M-eg;
	long double C = eg*(param->en+M-l1*cos(thg))-M*param->en-m2;

	if (m2*(pow(A,2.)-pow(B,2.))+pow(C,2.) <= 0.) return 0;

	if (std::abs(pow(A,2.)-pow(B,2.)) < 1.e-12) return 0;

	en1 = (B*C - A*sqrt(m2*(pow(A,2.)-pow(B,2.))+pow(C,2.)))/(pow(A,2.)-pow(B,2.));

	if (en1 < m || en1 > param->en - eg) return 0;

	l2 = sqrt(pow(en1,2.)-m2);

	long double test = m2 + M*(param->en-en1-eg) - param->en*(en1+eg) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*eg - l2*eg*sin(thl)*sin(thg)*cos(phig) - l2*eg*cos(thl)*cos(thg);

	if (std::abs(test) > 1e-16) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st_ps2(en1,thl,eg,thg,phig)*(egmax-egmin);

	return 0;
}

int Integrands::integrand_brems_2nd (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[3];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

	//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin > thgmax) return 0;

	//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
					*(param->max[param->E_prime]-param->min[param->E_prime])
					*(param->max[param->cos_thl]-param->min[param->cos_thl])
					*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
					*(thgmax-thgmin)*(thg1max-thg1min)*2.*(phig1max-phig1min);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
					+ cos(thg1)*cos(thg);
	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);
	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	if (std::abs(l1k1) < 1e-60 || std::abs(l1k2) < 1e-60
			|| std::abs(l2k1) < 1e-60 || std::abs(l2k2) < 1e-60) return 0;
	if (std::abs(l1k1+l1k2-k1k2) < 1e-60 || std::abs(l2k1+l2k2+k1k2) < 1e-60) return 0;

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = en1;
		FS->theta_l = thl;
		FS->E_gamma = eg;
		FS->theta_gamma = thg;
		FS->phi_gamma = phig;
		FS->E_gamma_prime = eg1;
		FS->theta_gamma_prime = thg1;
		FS->phi_gamma_prime = phig1;
		FS->weight = weight[0]*weight[8];
		if (param->flag[param->asymmetry] == 1)
			FS->weight *= 1. + param->P*CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			/CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1);
	};

	return 0;
}

int Integrands::integrand_brems_2nd_l1k1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[3];
	eg = exp(eg);
	Jacobian_eg = eg;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

	long double l1k1_min = param->en * eg - l1 * eg * thgmax;
	long double l1k1_max = param->en * eg - l1 * eg * thgmin;

	long double l1k1 = log(l1k1_min) + (log(l1k1_max) - log(l1k1_min)) * xx[4];

	l1k1 = exp(l1k1);
	long double Jacobian_thg = l1k1 / (l1 * eg);

	thg = (param->en * eg - l1k1) / (l1 * eg);

//	generates number for cos(theta'_gamma)
//	thg = thgmin+(thgmax-thgmin)*xx[4];

	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double phigmax = acos(phigmin_cos);
	long double phigmin = acos(phigmax_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	phig = phigmin+(phigmax-phigmin)*xx[5];
	if (rand.uniform() < 0.5) phig = 2.*pi - phig;
//	phig = acos(phig);

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}
	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_phig1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(log(l1k1_max)-log(l1k1_min))*Jacobian_thg*(thg1max-thg1min)
			*2.*(phigmax-phigmin);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

//	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);

	ff[0] *= d2*d3*d4/(d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4));

//	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
//			*((1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))/d1
//			+ (1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))/d2
//			+ (1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))/d3
//			+ (1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))/d4);

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = en1;
		FS->theta_l = thl;
		FS->E_gamma = eg;
		FS->theta_gamma = thg;
		FS->phi_gamma = phig;
		FS->E_gamma_prime = eg1;
		FS->theta_gamma_prime = thg1;
		FS->phi_gamma_prime = phig1;
		FS->weight = weight[0]*weight[8];
		if (param->flag[param->asymmetry] == 1)
			FS->weight *= 1. + param->P*CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
									/CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1);

	};

	return 0;
}

int Integrands::integrand_brems_2nd_l1k2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[3];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

//	long double gamma_test = pow(param->en-en1-eg-eg1+M,2.) - mpr2;

//	if (gamma < 0.) return 0.;

	thg1min = -1.;
	thg1max = 1.;

//	generates number for cos(theta'_gamma)
//	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	long double l1k2_min = param->en * eg1 - l1 * eg1 * thg1max;
	long double l1k2_max = param->en * eg1 - l1 * eg1 * thg1min;

	long double l1k2 = log(l1k2_min) + (log(l1k2_max) - log(l1k2_min)) * xx[4];

	l1k2 = exp(l1k2);
	long double Jacobian_thg1 = l1k2 / (l1 * eg1);

	thg1 = (param->en * eg1 - l1k2) / (l1 * eg1);

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

	if ( gamma <= 0.) return 0;

//	std::cout<<sqrt(gamma)<<"\t"<<eg*g3*sqrt(gamma_test)<<"\n";

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	 phig1max = acos(phig1min_cos);
	 phig1min = acos(phig1max_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;
//	phig1 = acos(phig1);

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(log(l1k2_max)-log(l1k2_min))*Jacobian_thg1
			*2.*(phig1max-phig1min);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
//	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);

	ff[0] *= d1*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4));

//	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
//			*((1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))/d1
//			+ (1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))/d2
//			+ (1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))/d3
//			+ (1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))/d4);

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = en1;
		FS->theta_l = thl;
		FS->E_gamma = eg;
		FS->theta_gamma = thg;
		FS->phi_gamma = phig;
		FS->E_gamma_prime = eg1;
		FS->theta_gamma_prime = thg1;
		FS->phi_gamma_prime = phig1;
		FS->weight = weight[0]*weight[8];
		if (param->flag[param->asymmetry] == 1)
			FS->weight *= 1. + param->P*CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
									/CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1);

	};

	return 0;
}

int Integrands::integrand_brems_2nd_l2k2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[3];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

//	 phig1max = acos(phig1min_cos);
//	 phig1min = acos(phig1max_cos);

	long double l2k2_min = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1max_cos
			- l2*eg1*cos(thl)*cos(thg1);
	long double l2k2_max = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1min_cos
			- l2*eg1*cos(thl)*cos(thg1);

	long double l2k2 = log(l2k2_min) + (log(l2k2_max) - log(l2k2_min)) * xx[5];

	l2k2 = exp(l2k2);

	phig1 = (en1*eg1 - l2*eg1*cos(thl)*cos(thg1) - l2k2)/(l2*eg1*sin(thl)*sin(thg1));

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	if (pow(phig1, 2.) >= 1.) return 0;

//	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	phig1 = acos(phig1);
	long double Jacobian_phig1 = l2k2 / (l2*eg1*sin(thl)*sin(thg1)) / sin(phig1);
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(thg1max-thg1min)
			*2.*(log(l2k2_max)-log(l2k2_min))*Jacobian_phig1;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
//	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);

	ff[0] *= d1*d2*d3 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3));

//	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
//			*((1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))/d1
//			+ (1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))/d2
//			+ (1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))/d3
//			+ (1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))/d4);

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = en1;
		FS->theta_l = thl;
		FS->E_gamma = eg;
		FS->theta_gamma = thg;
		FS->phi_gamma = phig;
		FS->E_gamma_prime = eg1;
		FS->theta_gamma_prime = thg1;
		FS->phi_gamma_prime = phig1;
		FS->weight = weight[0]*weight[8];
		if (param->flag[param->asymmetry] == 1)
			FS->weight *= 1. + param->P*CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
									/CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1);

	};

	return 0;
}

int Integrands::integrand_brems_2nd_l2k1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[3];
	eg = exp(eg);
	Jacobian_eg = eg;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

//	generates number for cos(theta'_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[4];

	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double l2k1_min = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmax_cos
			- l2*eg*cos(thl)*cos(thg);
	long double l2k1_max = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmin_cos
			- l2*eg*cos(thl)*cos(thg);

	long double l2k1 = log(l2k1_min) + (log(l2k1_max) - log(l2k1_min)) * xx[5];

	l2k1 = exp(l2k1);

	phig = (en1*eg - l2*eg*cos(thl)*cos(thg) - l2k1)/(l2*eg*sin(thl)*sin(thg));

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	if (pow(phig, 2.) >= 1.) return 0;

//	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	phig = acos(phig);
	long double Jacobian_phig = l2k1 / (l2*eg*sin(thl)*sin(thg)) / sin(phig);
	if (rand.uniform() < 0.5) phig = 2.*pi - phig;

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

//function to be calculated

	ff[0] = CS.crsect_brems_2nd_phig1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(thg1max-thg1min)
			*2.*(log(l2k1_max)-log(l2k1_min))*Jacobian_phig;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
//	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);

	ff[0] *= d1*d2*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4));

//	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
//			*((1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))/d1
//			+ (1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))/d2
//			+ (1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))/d3
//			+ (1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))/d4);

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = en1;
		FS->theta_l = thl;
		FS->E_gamma = eg;
		FS->theta_gamma = thg;
		FS->phi_gamma = phig;
		FS->E_gamma_prime = eg1;
		FS->theta_gamma_prime = thg1;
		FS->phi_gamma_prime = phig1;
		FS->weight = weight[0]*weight[8];
		if (param->flag[param->asymmetry] == 1)
			FS->weight *= 1. + param->P*CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
									/CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1);

	};

	return 0;
}

int Integrands::integrand_brems_2nd_l1k1_1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[3];
	eg = exp(eg);
	Jacobian_eg = eg;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

	long double l1k1_min = param->en * eg - l1 * eg * thgmax;
	long double l1k1_max = param->en * eg - l1 * eg * thgmin;

	long double l1k1 = log(l1k1_min) + (log(l1k1_max) - log(l1k1_min)) * xx[4];

	l1k1 = exp(l1k1);
	long double Jacobian_thg = l1k1 / (l1 * eg);

	thg = (param->en * eg - l1k1) / (l1 * eg);

//	generates number for cos(theta'_gamma)
//	thg = thgmin+(thgmax-thgmin)*xx[4];

	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double phigmax = acos(phigmin_cos);
	long double phigmin = acos(phigmax_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	phig = phigmin+(phigmax-phigmin)*xx[5];
	if (rand.uniform() < 0.5) phig = 2.*pi - phig;
//	phig = acos(phig);

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}
	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_l1k1_1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(log(l1k1_max)-log(l1k1_min))*Jacobian_thg*(thg1max-thg1min)
			*2.*(phigmax-phigmin);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

//	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= (1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))
			/(d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4);

	return 0;
}

int Integrands::integrand_brems_2nd_l1k2_1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[3];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

//	long double gamma_test = pow(param->en-en1-eg-eg1+M,2.) - mpr2;

//	if (gamma < 0.) return 0.;

	thg1min = -1.;
	thg1max = 1.;

//	generates number for cos(theta'_gamma)
//	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	long double l1k2_min = param->en * eg1 - l1 * eg1 * thg1max;
	long double l1k2_max = param->en * eg1 - l1 * eg1 * thg1min;

	long double l1k2 = log(l1k2_min) + (log(l1k2_max) - log(l1k2_min)) * xx[4];

	l1k2 = exp(l1k2);
	long double Jacobian_thg1 = l1k2 / (l1 * eg1);

	thg1 = (param->en * eg1 - l1k2) / (l1 * eg1);

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

	if ( gamma <= 0.) return 0;

//	std::cout<<sqrt(gamma)<<"\t"<<eg*g3*sqrt(gamma_test)<<"\n";

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	 phig1max = acos(phig1min_cos);
	 phig1min = acos(phig1max_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;
//	phig1 = acos(phig1);

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_2nd_l1k2_1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(log(l1k2_max)-log(l1k2_min))*Jacobian_thg1
			*2.*(phig1max-phig1min);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
//	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= (1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))
		/(d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4);

	return 0;
}

int Integrands::integrand_brems_2nd_l2k2_1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[3];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

//	 phig1max = acos(phig1min_cos);
//	 phig1min = acos(phig1max_cos);

	long double l2k2_min = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1max_cos
			- l2*eg1*cos(thl)*cos(thg1);
	long double l2k2_max = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1min_cos
			- l2*eg1*cos(thl)*cos(thg1);

	long double l2k2 = log(l2k2_min) + (log(l2k2_max) - log(l2k2_min)) * xx[5];

	l2k2 = exp(l2k2);

	phig1 = (en1*eg1 - l2*eg1*cos(thl)*cos(thg1) - l2k2)/(l2*eg1*sin(thl)*sin(thg1));

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	if (pow(phig1, 2.) >= 1.) return 0;

//	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	phig1 = acos(phig1);
	long double Jacobian_phig1 = l2k2 / (l2*eg1*sin(thl)*sin(thg1)) / sin(phig1);
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_l2k2_1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(thg1max-thg1min)
			*2.*(log(l2k2_max)-log(l2k2_min))*Jacobian_phig1;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
//	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= (1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))
		/(d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4);

	return 0;
}

//int Integrands::integrand_brems_2nd_l2k2_1 (const double xx[], double ff[], const double weight[])const {
//
//	ff[0] = 0.;
//	//xmin and xmax are limits which are defined in the param
//	//generates number for the energy of the scattered electron in the limits fixed by the param
//	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
//	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2
//
//	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
//	thl = acos(thl);
//
//	egmin = param->min[param->E_gamma];
//	egmax = param->en - en1;
//
//	//generates photon energy
//	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
//	eg = exp(eg);
//	Jacobian_eg = eg;
//
//	eg1min = param->min[param->E_gamma_prime];
//	eg1max = param->en - en1 - eg;
//
//	if (eg1max <= eg1min) return 0;
//
//	//generates photon energy
//	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[3];
//	eg1 = exp(eg1);
//	Jacobian_eg1 = eg1;
//
//	thg1min = -1.;
//	thg1max = 1.;
//
////	generates number for cos(theta'_gamma)
////	thg1 = thg1min+(thg1max-thg1min)*xx[4];
////	std::cout << en1 <<"\t"<< eg1 <<"\n";
//
//	long double l2k2_min = en1 * eg1 - l2 * eg1 * thg1max;
//	long double l2k2_max = en1 * eg1 - l2 * eg1 * thg1min;
//
//	long double l2k2 = log(l2k2_min) + (log(l2k2_max) - log(l2k2_min)) * xx[4];
//
//	l2k2 = exp(l2k2);
//	long double Jacobian_thg1 = l2k2 / (l2 * eg1);
//
//	thg1 = (en1 * eg1 - l2k2) / (l2 * eg1);
//
//	if (pow(thg1, 2.) >= 1.) return 0;
//
//	thg1 = acos(thg1);
//
////	std::cout << thg1 << "\n";
//
//	long double g1 = -pow(eg,2.)*l1*eg1*sin(thl)*sin(thg1);
//	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l1*eg*sin(thl),2.);
//	long double g3 = - l1 * eg1 * sin(thl) * sin(thg1);
//	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
//			+ eg*eg1 + l1*l2*cos(thl) - l2*eg1*cos(thg1) + l1*eg1*cos(thg1)*cos(thl);
//	long double beta2 = -l2*eg + l1*eg*cos(thl) - eg*eg1*cos(thg1);
//
//	long double gamma2 = g1 + g3*g4;
//	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
//	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
////	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);
//
////	std::cout<<gamma<<"\t"<<gamma_test<<"\n";
//
//	if ( gamma <= 0.) return 0;
//
////	std::cout << gamma << "\n";
//
//	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
//	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);
//
//
//	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;
//
//	if (phig1max_cos > 1.) phig1max_cos = 1.;
//	if (phig1min_cos < -1.) phig1min_cos = -1.;
//
//	 phig1max = acos(phig1min_cos);
//	 phig1min = acos(phig1max_cos);
//
////	std::cout << phig1min <<"\t"<< phig1max << "\n";
//
//	phig1 = phig1min+(phig1max-phig1min)*xx[5];
//	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;
//
//	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;
//
//	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
//	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) -
//			l2*eg1*cos(thg1) + en1*(eg+eg1) + l1*eg1*sin(thl)*sin(thg1)*cos(phig1) +
//			l1*eg1*cos(thl)*cos(thg1) + eg*eg1;
//	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);
//
//	if (beta <= 0.) return 0;
//
////	std::cout<<beta<<"\n";
//
//	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
////	std::cout << thgmin <<"\t"<< thgmax << "\n";
//
//	if (thgmin < -1.) thgmin = -1.;
//	if (thgmax > 1.) thgmax = 1.;
//	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
//	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
//	if (thgmin > thgmax) return 0;
//
////	generates number for cos(theta_gamma)
//	thg = thgmin+(thgmax-thgmin)*xx[6];
//
//	//check if -1<cos(thg)<1
//	if (pow(thg, 2.) >= 1.) return 0;
//
////	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;
//
//	thg = acos(thg);
//
//	//check if the generated values are real numbers
//	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;
//
//	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
//	long double alpha2 = -l1*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
//	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1) + l1*l2*cos(thl) +
//			l1*eg*cos(thg)*cos(thl) + l1*eg1*cos(thg1)*cos(thl) + l1*eg1*sin(thl)*sin(thg1)*cos(phig1) -
//			l2*eg*cos(thg) - l2*eg1*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);
//
////	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;
//
//	if ( pow(alpha3,2.) / pow(alpha1,2.)+pow(alpha2,2.) > 1.) return 0;
//
//	if (alpha1 > 0.)
//		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
//	else {
//		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
//		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
//	}
//
//	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) -
//			l2*eg*cos(thg) - l2*eg1*cos(thg1) + en1*(eg+eg1) + l1*eg*sin(thl)*sin(thg)*cos(phig) +
//			l1*eg1*sin(thl)*sin(thg1)*cos(phig1) + l1*eg*cos(thl)*cos(thg) + l1*eg1*cos(thl)*cos(thg1) +
//			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
//					cos(thg)*cos(thg1));
//
////	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
////	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
////	std::cout<<test<<"\n";
//
//	if (std::abs(test) > 1e-15) return 0;
//
////	if (phig < 0.) phig = 2.*pi + phig;
//
//	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
//	if (x < 1e-60) return 0;
//
//	//function to be calculated
//
//	ff[0] = CS.crsect_brems_2nd_l2k2_1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
//			*(param->max[param->E_prime]-param->min[param->E_prime])
//			*(param->max[param->cos_thl]-param->min[param->cos_thl])
//			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
//			*(thgmax-thgmin)*(l2k2_max-l2k2_min)*Jacobian_thg1*2.*(phig1max-phig1min);
//
//	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
//	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);
//
//	long double l2k1 = en1*eg - l2*eg*cos(thg);
//	long double l1k1 = param->en*eg - eg*param->l1*cospsi;
//	long double l1k2 = param->en*eg1 - eg1*param->l1*cospsi1;
//
//	long double d1 = l1k1;
//	long double d2 = l1k2;
//	long double d3 = l2k1;
//	long double d4 = l2k2;
//
//	ff[0] *= (1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))
//					/(d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4);
//
//	return 0;
//}

int Integrands::integrand_brems_2nd_l2k1_1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[3];
	eg = exp(eg);
	Jacobian_eg = eg;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

//	generates number for cos(theta'_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[4];

	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double l2k1_min = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmax_cos
			- l2*eg*cos(thl)*cos(thg);
	long double l2k1_max = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmin_cos
			- l2*eg*cos(thl)*cos(thg);

	long double l2k1 = log(l2k1_min) + (log(l2k1_max) - log(l2k1_min)) * xx[5];

	l2k1 = exp(l2k1);

	phig = (en1*eg - l2*eg*cos(thl)*cos(thg) - l2k1)/(l2*eg*sin(thl)*sin(thg));

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	if (pow(phig, 2.) >= 1.) return 0;

//	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	phig = acos(phig);
	long double Jacobian_phig = l2k1 / (l2*eg*sin(thl)*sin(thg)) / sin(phig);
	if (rand.uniform() < 0.5) phig = 2.*pi - phig;

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

//function to be calculated

	ff[0] = CS.crsect_brems_2nd_l2k1_1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(thg1max-thg1min)
			*2.*(log(l2k1_max)-log(l2k1_min))*Jacobian_phig;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
//	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= (1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))
		/(d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4);

	return 0;
}


int Integrands::integrand_brems_2nd_l1k1_2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[3];
	eg = exp(eg);
	Jacobian_eg = eg;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

	long double l1k1_min = 1./(param->en * eg - l1 * eg * thgmin);
	long double l1k1_max = 1./(param->en * eg - l1 * eg * thgmax);

	long double l1k1 = l1k1_min + (l1k1_max - l1k1_min) * xx[4];

	l1k1 = 1./l1k1;
	long double Jacobian_thg = pow(l1k1,2.) / (l1 * eg);

	thg = (param->en * eg - l1k1) / (l1 * eg);

//	generates number for cos(theta'_gamma)
//	thg = thgmin+(thgmax-thgmin)*xx[4];

	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double phigmax = acos(phigmin_cos);
	long double phigmin = acos(phigmax_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	phig = phigmin+(phigmax-phigmin)*xx[5];
	if (rand.uniform() < 0.5) phig = 2.*pi - phig;
//	phig = acos(phig);

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}
	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_l1k1_2(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(l1k1_max-l1k1_min)*Jacobian_thg*(thg1max-thg1min)
			*2.*(phigmax-phigmin);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

//	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

//	long double d1 = pow(l1k1,2.);
//	long double d2 = pow(l1k2,2.);
//	long double d3 = pow(l2k1,2.);
//	long double d4 = pow(l2k2,2.);

//	ff[0] *= (1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))
//			/(d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4);

	return 0;
}

int Integrands::integrand_brems_2nd_l1k2_2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[3];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

//	long double gamma_test = pow(param->en-en1-eg-eg1+M,2.) - mpr2;

//	if (gamma < 0.) return 0.;

	thg1min = -1.;
	thg1max = 1.;

//	generates number for cos(theta'_gamma)
//	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	long double l1k2_min = 1./(param->en * eg1 - l1 * eg1 * thg1min);
	long double l1k2_max = 1./(param->en * eg1 - l1 * eg1 * thg1max);

	long double l1k2 = l1k2_min + (l1k2_max - l1k2_min) * xx[4];

	l1k2 = 1./l1k2;
	long double Jacobian_thg1 = pow(l1k2,2.) / (l1 * eg1);

	thg1 = (param->en * eg1 - l1k2) / (l1 * eg1);

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

	if ( gamma <= 0.) return 0;

//	std::cout<<sqrt(gamma)<<"\t"<<eg*g3*sqrt(gamma_test)<<"\n";

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	 phig1max = acos(phig1min_cos);
	 phig1min = acos(phig1max_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;
//	phig1 = acos(phig1);

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_2nd_l1k2_2(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(l1k2_max-l1k2_min)*Jacobian_thg1
			*2.*(phig1max-phig1min);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
//	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

//	long double d1 = pow(l1k1,2.);
//	long double d2 = pow(l1k2,2.);
//	long double d3 = pow(l2k1,2.);
//	long double d4 = pow(l2k2,2.);

//	ff[0] *= (1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))
//		/(d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4);

	return 0;
}

int Integrands::integrand_brems_2nd_l2k2_2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[3];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

//	 phig1max = acos(phig1min_cos);
//	 phig1min = acos(phig1max_cos);

	long double l2k2_min = 1./(en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1min_cos
			- l2*eg1*cos(thl)*cos(thg1));
	long double l2k2_max = 1./(en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1max_cos
			- l2*eg1*cos(thl)*cos(thg1));

	long double l2k2 = l2k2_min + (l2k2_max - l2k2_min) * xx[5];

	l2k2 = 1./l2k2;

	phig1 = (en1*eg1 - l2*eg1*cos(thl)*cos(thg1) - l2k2)/(l2*eg1*sin(thl)*sin(thg1));

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	if (pow(phig1, 2.) >= 1.) return 0;

//	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	phig1 = acos(phig1);
	long double Jacobian_phig1 = pow(l2k2,2.) / (l2*eg1*sin(thl)*sin(thg1)) / sin(phig1);
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_l2k2_2(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(thg1max-thg1min)
			*2.*(l2k2_max-l2k2_min)*Jacobian_phig1;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
//	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

//	long double d1 = pow(l1k1,2.);
//	long double d2 = pow(l1k2,2.);
//	long double d3 = pow(l2k1,2.);
//	long double d4 = pow(l2k2,2.);

//	ff[0] *= (1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))
//		/(d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4);

	return 0;
}

//int Integrands::integrand_brems_2nd_l2k2_2 (const double xx[], double ff[], const double weight[])const {
//
//	ff[0] = 0.;
//	//xmin and xmax are limits which are defined in the param
//	//generates number for the energy of the scattered electron in the limits fixed by the param
//	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
//	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2
//
//	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
//	thl = acos(thl);
//
//	egmin = param->min[param->E_gamma];
//	egmax = param->en - en1;
//
//	//generates photon energy
//	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
//	eg = exp(eg);
//	Jacobian_eg = eg;
//
//	eg1min = param->min[param->E_gamma_prime];
//	eg1max = param->en - en1 - eg;
//
//	if (eg1max <= eg1min) return 0;
//
//	//generates photon energy
//	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[3];
//	eg1 = exp(eg1);
//	Jacobian_eg1 = eg1;
//
//	thg1min = -1.;
//	thg1max = 1.;
//
////	generates number for cos(theta'_gamma)
////	thg1 = thg1min+(thg1max-thg1min)*xx[4];
////	std::cout << en1 <<"\t"<< eg1 <<"\n";
//
//	long double l2k2_min = 1./(en1 * eg1 - l2 * eg1 * thg1min);
//	long double l2k2_max = 1./(en1 * eg1 - l2 * eg1 * thg1max);
//
//	long double l2k2 = l2k2_min + (l2k2_max - l2k2_min) * xx[4];
//
//	l2k2 = 1./l2k2;
//	long double Jacobian_thg1 = pow(l2k2,2.) / (l2 * eg1);
//
//	thg1 = (en1*eg1 - l2k2) / (l2 * eg1);
//
//	if (pow(thg1, 2.) >= 1.) return 0;
//
//	thg1 = acos(thg1);
//
////	std::cout << thg1 << "\n";
//
//	long double g1 = -pow(eg,2.)*l1*eg1*sin(thl)*sin(thg1);
//	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l1*eg*sin(thl),2.);
//	long double g3 = - l1 * eg1 * sin(thl) * sin(thg1);
//	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
//			+ eg*eg1 + l1*l2*cos(thl) - l2*eg1*cos(thg1) + l1*eg1*cos(thg1)*cos(thl);
//	long double beta2 = -l2*eg + l1*eg*cos(thl) - eg*eg1*cos(thg1);
//
//	long double gamma2 = g1 + g3*g4;
//	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
//	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
////	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);
//
////	std::cout<<gamma<<"\t"<<gamma_test<<"\n";
//
//	if ( gamma <= 0.) return 0;
//
////	std::cout << gamma << "\n";
//
//	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
//	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);
//
//
//	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;
//
//	if (phig1max_cos > 1.) phig1max_cos = 1.;
//	if (phig1min_cos < -1.) phig1min_cos = -1.;
//
//	 phig1max = acos(phig1min_cos);
//	 phig1min = acos(phig1max_cos);
//
////	std::cout << phig1min <<"\t"<< phig1max << "\n";
//
//	phig1 = phig1min+(phig1max-phig1min)*xx[5];
//	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;
//
//	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;
//
//	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
//	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) -
//			l2*eg1*cos(thg1) + en1*(eg+eg1) + l1*eg1*sin(thl)*sin(thg1)*cos(phig1) +
//			l1*eg1*cos(thl)*cos(thg1) + eg*eg1;
//	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);
//
//	if (beta <= 0.) return 0;
//
////	std::cout<<beta<<"\n";
//
//	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
////	std::cout << thgmin <<"\t"<< thgmax << "\n";
//
//	if (thgmin < -1.) thgmin = -1.;
//	if (thgmax > 1.) thgmax = 1.;
//	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
//	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
//	if (thgmin > thgmax) return 0;
//
////	generates number for cos(theta_gamma)
//	thg = thgmin+(thgmax-thgmin)*xx[6];
//
//	//check if -1<cos(thg)<1
//	if (pow(thg, 2.) >= 1.) return 0;
//
////	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;
//
//	thg = acos(thg);
//
//	//check if the generated values are real numbers
//	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;
//
//	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
//	long double alpha2 = -l1*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
//	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1) + l1*l2*cos(thl) +
//			l1*eg*cos(thg)*cos(thl) + l1*eg1*cos(thg1)*cos(thl) + l1*eg1*sin(thl)*sin(thg1)*cos(phig1) -
//			l2*eg*cos(thg) - l2*eg1*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);
//
////	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;
//
//	if ( pow(alpha3,2.) / pow(alpha1,2.)+pow(alpha2,2.) > 1.) return 0;
//
//	if (alpha1 > 0.)
//		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
//	else {
//		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
//		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
//	}
//
//	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) -
//			l2*eg*cos(thg) - l2*eg1*cos(thg1) + en1*(eg+eg1) + l1*eg*sin(thl)*sin(thg)*cos(phig) +
//			l1*eg1*sin(thl)*sin(thg1)*cos(phig1) + l1*eg*cos(thl)*cos(thg) + l1*eg1*cos(thl)*cos(thg1) +
//			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
//					cos(thg)*cos(thg1));
//
////	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
////	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
////	std::cout<<test<<"\n";
//
//	if (std::abs(test) > 1e-15) return 0;
//
////	if (phig < 0.) phig = 2.*pi + phig;
//
//	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
//	if (x < 1e-60) return 0;
//
//	//function to be calculated
//
//	ff[0] = CS.crsect_brems_2nd_l2k2_2(en1,thl,eg,thg,phig,eg1,thg1,phig1)
//			*(param->max[param->E_prime]-param->min[param->E_prime])
//			*(param->max[param->cos_thl]-param->min[param->cos_thl])
//			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
//			*(thgmax-thgmin)*(l2k2_max-l2k2_min)*Jacobian_thg1*2.*(phig1max-phig1min);
//
////	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
////	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);
//
////	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
////	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
////	long double l2k1 = en1*eg - eg*l2*cospsi;
////	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;
//
////	long double d1 = pow(l1k1,2.);
////	long double d2 = pow(l1k2,2.);
////	long double d3 = pow(l2k1,2.);
////	long double d4 = pow(l2k2,2.);
//
////	ff[0] *= (1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))
////		/(d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4);
//
//	return 0;
//}

int Integrands::integrand_brems_2nd_l2k1_2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[3];
	eg = exp(eg);
	Jacobian_eg = eg;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

//	generates number for cos(theta'_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[4];

	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double l2k1_min = 1./(en1*eg - l2*eg*sin(thl)*sin(thg)*phigmin_cos
			- l2*eg*cos(thl)*cos(thg));
	long double l2k1_max = 1./(en1*eg - l2*eg*sin(thl)*sin(thg)*phigmax_cos
			- l2*eg*cos(thl)*cos(thg));

	long double l2k1 = l2k1_min + (l2k1_max - l2k1_min) * xx[5];

	l2k1 = 1./l2k1;

	phig = (en1*eg - l2*eg*cos(thl)*cos(thg) - l2k1)/(l2*eg*sin(thl)*sin(thg));

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	if (pow(phig, 2.) >= 1.) return 0;

//	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	phig = acos(phig);
	long double Jacobian_phig = pow(l2k1,2.) / (l2*eg*sin(thl)*sin(thg)) / sin(phig);
	if (rand.uniform() < 0.5) phig = 2.*pi - phig;

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

//function to be calculated

	ff[0] = CS.crsect_brems_2nd_l2k1_2(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(thg1max-thg1min)
			*2.*(l2k1_max-l2k1_min)*Jacobian_phig;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
//	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

//	long double d1 = pow(l1k1,2.);
//	long double d2 = pow(l1k2,2.);
//	long double d3 = pow(l2k1,2.);
//	long double d4 = pow(l2k2,2.);

//	ff[0] *= (1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))
//		/(d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4);

	return 0;
}

//int Integrands::integrand_brems_2nd_l2k1_2 (const double xx[], double ff[], const double weight[])const {
//
//	ff[0] = 0.;
//	//xmin and xmax are limits which are defined in the param
//	//generates number for the energy of the scattered electron in the limits fixed by the param
//	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
//	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2
//
//	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
//	thl = acos(thl);
//
//	eg1min = param->min[param->E_gamma_prime];
//	eg1max = param->en - en1;
//
//	//generates photon energy
//	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
//	eg1 = exp(eg1);
//	Jacobian_eg1 = eg1;
//
//	egmin = param->min[param->E_gamma];
//	egmax = param->en - en1 - eg1;
//
//	//generates photon energy
//	eg = log(egmin)+(log(egmax)-log(egmin))*xx[3];
//	eg = exp(eg);
//	Jacobian_eg = eg;
//
//	thgmin = param->min[param->cos_thg];
//	thgmax = param->max[param->cos_thg];
//
////	generates number for cos(theta'_gamma)
//	thg = thgmin+(thgmax-thgmin)*xx[4];
//
//	if (pow(thg, 2.) >= 1.) return 0;
//
////	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;
//
//	thg = acos(thg);
//
//	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
//	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
//	long double g3 = l2 * eg * sin(thl) * sin(thg);
//	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
//			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
//	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);
//
//	long double gamma2 = g1 + g3*g4;
//	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
//	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
////	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);
//
////	std::cout<<gamma<<"\t"<<gamma_test<<"\n";
//
//	if ( gamma <= 0.) return 0;
//
////	if (param->flag[param->ps] == 1) FS->acc_ev ++;
//
//	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
//	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);
//
////	std::cout<<phig1min<<"\t"<<phig1max<<"\n";
//
//	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;
//
//	if (phigmax_cos > 1.) phigmax_cos = 1.;
//	if (phigmin_cos < -1.) phigmin_cos = -1.;
//
//	long double l2k1_min = 1./(en1*eg - l2*eg*sin(thl)*sin(thg)*phigmin_cos
//			- l2*eg*cos(thl)*cos(thg));
//	long double l2k1_max = 1./(en1*eg - l2*eg*sin(thl)*sin(thg)*phigmax_cos
//			- l2*eg*cos(thl)*cos(thg));
//
//	long double l2k1 = l2k1_min + (l2k1_max - l2k1_min) * xx[5];
//
//	l2k1 = 1./l2k1;
//
//	phig = (en1*eg - l2*eg*cos(thl)*cos(thg) - l2k1)/(l2*eg*sin(thl)*sin(thg));
//
////	std::cout << phig1min <<"\t"<< phig1max << "\n";
//
//	if (pow(phig, 2.) >= 1.) return 0;
//
////	phig1 = phig1min+(phig1max-phig1min)*xx[5];
//	phig = acos(phig);
//	long double Jacobian_phig = pow(l2k1,2.) / (l2*eg*sin(thl)*sin(thg)) / sin(phig);
//	if (rand.uniform() < 0.5) phig = 2.*pi - phig;
//
//	if (2.*g1 * cos(phig) + g2 < 0.) return 0;
//
//	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
//	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
//			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
//			l2*eg*cos(thl)*cos(thg) + eg*eg1;
//	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);
//
//	if (beta <= 0.) return 0;
//
////	std::cout<<beta<<"\n";
//
//	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
////
////	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
////
////	if (thgmax < -1. || thgmin > 1.) return 0;
////
//	if (thg1min < -1.) thg1min = -1.;
//	if (thg1max > 1.) thg1max = 1.;
//	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
//	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
//	if (thg1min > thg1max) return 0;
//
////	generates number for cos(theta_gamma)
//	thg1 = thg1min+(thg1max-thg1min)*xx[6];
//
//	//check if -1<cos(thg)<1
//	if (pow(thg1, 2.) >= 1.) return 0;
//
//	thg1 = acos(thg1);
//
//	//check if the generated values are real numbers
//	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;
//
////	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;
//
//	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
//	long double alpha2 = -l1*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
//	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1) + l1*l2*cos(thl)
//			+ l1*eg*cos(thg)*cos(thl) + l1*eg1*cos(thg1)*cos(thl) + l1*eg*sin(thl)*sin(thg)*cos(phig)
//			- l2*eg*cos(thg) - l2*eg1*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);
//
////	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;
//
//	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;
//
//	if (alpha1 > 0.)
//		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
//	else {
//		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
//		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
//	}
//
//	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) -
//			l2*eg*cos(thg) - l2*eg1*cos(thg1) + en1*(eg+eg1) + l1*eg*sin(thl)*sin(thg)*cos(phig) +
//			l1*eg1*sin(thl)*sin(thg1)*cos(phig1) + l1*eg*cos(thl)*cos(thg) + l1*eg1*cos(thl)*cos(thg1) +
//			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
//					cos(thg)*cos(thg1));
//
////	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
////	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
////	std::cout<<test<<"\n";
//
//	if (std::abs(test) > 1e-15) return 0;
//
////	if (phig < 0.) phig = 2.*pi + phig;
//
//	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
//	if (x < 1e-60) return 0;
//
////function to be calculated
//
//	ff[0] = CS.crsect_brems_2nd_l2k1_2(en1,thl,eg,thg,phig,eg1,thg1,phig1)
//			*(param->max[param->E_prime]-param->min[param->E_prime])
//			*(param->max[param->cos_thl]-param->min[param->cos_thl])
//			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
//			*(thgmax-thgmin)*(thg1max-thg1min)
//			*2.*(l2k1_max-l2k1_min)*Jacobian_phig;
//
//	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
//	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);
//
//	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
//	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
////	long double l2k1 = en1*eg - eg*l2*cospsi;
//	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;
//
////	long double d1 = pow(l1k1,2.);
////	long double d2 = pow(l1k2,2.);
////	long double d3 = pow(l2k1,2.);
////	long double d4 = pow(l2k2,2.);
//
////	ff[0] *= (1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))
////		/(d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4);
//
//	return 0;
//}

int Integrands::integrand_brems_2nd_add (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma];

	//generates photon energy
	eg1 = eg1max*xx[3];

	thg1min = -1.;
	thg1max = 1.;

	//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin > thgmax) return 0;

	//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
					+ cos(thg1)*cos(thg);
	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);
	long double d5 = pow(l1k1 + l1k2 - k1k2,2.);
	long double d6 = pow(l2k1 + l2k2 + k1k2,2.);

	if (d1*d3*d5 < 1e-60) return 0;
	if (d2*d4*d6 < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
					*(param->max[param->E_prime]-param->min[param->E_prime])
					*(param->max[param->cos_thl]-param->min[param->cos_thl])
					*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
					*(thgmax-thgmin)*(thg1max-thg1min)*2.*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_brems_2nd_add_l1k1_1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[2];

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[3];
	eg = exp(eg);
	Jacobian_eg = eg;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

	long double l1k1_min = param->en * eg - l1 * eg * thgmax;
	long double l1k1_max = param->en * eg - l1 * eg * thgmin;

	long double l1k1 = log(l1k1_min) + (log(l1k1_max) - log(l1k1_min)) * xx[4];

	l1k1 = exp(l1k1);
	long double Jacobian_thg = l1k1 / (l1 * eg);

	thg = (param->en * eg - l1k1) / (l1 * eg);

	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double phigmax = acos(phigmin_cos);
	long double phigmin = acos(phigmax_cos);

	phig = phigmin+(phigmax-phigmin)*xx[5];

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add_phig1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(log(l1k1_max)-log(l1k1_min))*Jacobian_thg*(thg1max-thg1min)
			*(phigmax-phigmin);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))/d1;

	return 0;
}

int Integrands::integrand_brems_2nd_add_l1k1_2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[2];

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[3];
	eg = exp(eg);
	Jacobian_eg = eg;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

	long double l1k1_min = param->en * eg - l1 * eg * thgmax;
	long double l1k1_max = param->en * eg - l1 * eg * thgmin;

	long double l1k1 = log(l1k1_min) + (log(l1k1_max) - log(l1k1_min)) * xx[4];

	l1k1 = exp(l1k1);
	long double Jacobian_thg = l1k1 / (l1 * eg);

	thg = (param->en * eg - l1k1) / (l1 * eg);

	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double phigmax = acos(phigmin_cos);
	long double phigmin = acos(phigmax_cos);

	phig = phigmin+(phigmax-phigmin)*xx[5];
	phig = 2.*pi - phig;

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add_phig1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(log(l1k1_max)-log(l1k1_min))*Jacobian_thg*(thg1max-thg1min)
			*(phigmax-phigmin);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))/d1;

	return 0;
}

int Integrands::integrand_brems_2nd_add_l1k2_1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[3];

	thg1min = -1.;
	thg1max = 1.;

	//	generates number for cos(theta'_gamma)
	//	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	long double l1k2_min = param->en * eg1 - l1 * eg1 * thg1max;
	long double l1k2_max = param->en * eg1 - l1 * eg1 * thg1min;

	long double l1k2 = log(l1k2_min) + (log(l1k2_max) - log(l1k2_min)) * xx[4];

	l1k2 = exp(l1k2);
	long double Jacobian_thg1 = l1k2 / (l1 * eg1);

	thg1 = (param->en * eg1 - l1k2) / (l1 * eg1);

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
	if (thgmin > thgmax) return 0;

	//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));


	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
					*(param->max[param->E_prime]-param->min[param->E_prime])
					*(param->max[param->cos_thl]-param->min[param->cos_thl])
					*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
					*(thgmax-thgmin)*(log(l1k2_max)-log(l1k2_min))*Jacobian_thg1
					*(phig1max-phig1min);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
					*(1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))/d2;

	return 0;
}

int Integrands::integrand_brems_2nd_add_l1k2_2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[3];

	thg1min = -1.;
	thg1max = 1.;

	//	generates number for cos(theta'_gamma)
	//	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	long double l1k2_min = param->en * eg1 - l1 * eg1 * thg1max;
	long double l1k2_max = param->en * eg1 - l1 * eg1 * thg1min;

	long double l1k2 = log(l1k2_min) + (log(l1k2_max) - log(l1k2_min)) * xx[4];

	l1k2 = exp(l1k2);
	long double Jacobian_thg1 = l1k2 / (l1 * eg1);

	thg1 = (param->en * eg1 - l1k2) / (l1 * eg1);

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
	if (thgmin > thgmax) return 0;

	//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
					*(param->max[param->E_prime]-param->min[param->E_prime])
					*(param->max[param->cos_thl]-param->min[param->cos_thl])
					*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
					*(thgmax-thgmin)*(log(l1k2_max)-log(l1k2_min))*Jacobian_thg1
					*(phig1max-phig1min);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
					*(1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))/d2;

	return 0;
}



int Integrands::integrand_brems_2nd_add_l2k2_1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[3];

	thg1min = -1.;
	thg1max = 1.;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	long double l2k2_min = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1max_cos
			- l2*eg1*cos(thl)*cos(thg1);
	long double l2k2_max = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1min_cos
			- l2*eg1*cos(thl)*cos(thg1);

	long double l2k2 = log(l2k2_min) + (log(l2k2_max) - log(l2k2_min)) * xx[5];

	l2k2 = exp(l2k2);

	phig1 = (en1*eg1 - l2*eg1*cos(thl)*cos(thg1) - l2k2)/(l2*eg1*sin(thl)*sin(thg1));

	if (pow(phig1, 2.) >= 1.) return 0;

	phig1 = acos(phig1);
	long double Jacobian_phig1 = l2k2 / (l2*eg1*sin(thl)*sin(thg1)) / sin(phig1);

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thgmax-thgmin)*(thg1max-thg1min)
			*(log(l2k2_max)-log(l2k2_min))*Jacobian_phig1;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
//	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))/d4;

	return 0;
}

int Integrands::integrand_brems_2nd_add_l2k2_2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[3];

	thg1min = -1.;
	thg1max = 1.;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	long double l2k2_min = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1max_cos
			- l2*eg1*cos(thl)*cos(thg1);
	long double l2k2_max = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1min_cos
			- l2*eg1*cos(thl)*cos(thg1);

	long double l2k2 = log(l2k2_min) + (log(l2k2_max) - log(l2k2_min)) * xx[5];

	l2k2 = exp(l2k2);

	phig1 = (en1*eg1 - l2*eg1*cos(thl)*cos(thg1) - l2k2)/(l2*eg1*sin(thl)*sin(thg1));

	if (pow(phig1, 2.) >= 1.) return 0;

	phig1 = acos(phig1);
	long double Jacobian_phig1 = l2k2 / (l2*eg1*sin(thl)*sin(thg1)) / sin(phig1);
	phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thgmax-thgmin)*(thg1max-thg1min)
			*(log(l2k2_max)-log(l2k2_min))*Jacobian_phig1;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
//	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))/d4;

	return 0;
}

int Integrands::integrand_brems_2nd_add_l2k1_1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[2];

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[3];
	eg = exp(eg);
	Jacobian_eg = eg;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

//	generates number for cos(theta'_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[4];

	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double l2k1_min = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmax_cos
			- l2*eg*cos(thl)*cos(thg);
	long double l2k1_max = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmin_cos
			- l2*eg*cos(thl)*cos(thg);

	long double l2k1 = log(l2k1_min) + (log(l2k1_max) - log(l2k1_min)) * xx[5];

	l2k1 = exp(l2k1);

	phig = (en1*eg - l2*eg*cos(thl)*cos(thg) - l2k1)/(l2*eg*sin(thl)*sin(thg));

	if (pow(phig, 2.) >= 1.) return 0;

	phig = acos(phig);
	long double Jacobian_phig = l2k1 / (l2*eg*sin(thl)*sin(thg)) / sin(phig);

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add_phig1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thgmax-thgmin)*(thg1max-thg1min)
			*(log(l2k1_max)-log(l2k1_min))*Jacobian_phig;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
//	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))/d3;

	return 0;
}

int Integrands::integrand_brems_2nd_add_l2k1_2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[2];

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[3];
	eg = exp(eg);
	Jacobian_eg = eg;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

//	generates number for cos(theta'_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[4];

	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double l2k1_min = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmax_cos
			- l2*eg*cos(thl)*cos(thg);
	long double l2k1_max = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmin_cos
			- l2*eg*cos(thl)*cos(thg);

	long double l2k1 = log(l2k1_min) + (log(l2k1_max) - log(l2k1_min)) * xx[5];

	l2k1 = exp(l2k1);

	phig = (en1*eg - l2*eg*cos(thl)*cos(thg) - l2k1)/(l2*eg*sin(thl)*sin(thg));

	if (pow(phig, 2.) >= 1.) return 0;

	phig = acos(phig);
	long double Jacobian_phig = l2k1 / (l2*eg*sin(thl)*sin(thg)) / sin(phig);
	phig = 2.*pi - phig;

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add_phig1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thgmax-thgmin)*(thg1max-thg1min)
			*(log(l2k1_max)-log(l2k1_min))*Jacobian_phig;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
//	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))/d3;

	return 0;
}

int Integrands::integrand_brems_2nd_add_l1k1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[2];

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[3];
	eg = exp(eg);
	Jacobian_eg = eg;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

	long double l1k1_min = param->en * eg - l1 * eg * thgmax;
	long double l1k1_max = param->en * eg - l1 * eg * thgmin;

	long double l1k1 = log(l1k1_min) + (log(l1k1_max) - log(l1k1_min)) * xx[4];

	l1k1 = exp(l1k1);
	long double Jacobian_thg = l1k1 / (l1 * eg);

	thg = (param->en * eg - l1k1) / (l1 * eg);

	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double phigmax = acos(phigmin_cos);
	long double phigmin = acos(phigmax_cos);

	phig = phigmin+(phigmax-phigmin)*xx[5];
	if (rand.uniform() < 0.5) phig = 2.*pi - phig;

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add_phig1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(log(l1k1_max)-log(l1k1_min))*Jacobian_thg*(thg1max-thg1min)
			*2.*(phigmax-phigmin);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))/d1;

	return 0;
}

int Integrands::integrand_brems_2nd_add_l1k2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[3];

	thg1min = -1.;
	thg1max = 1.;

	//	generates number for cos(theta'_gamma)
	//	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	long double l1k2_min = param->en * eg1 - l1 * eg1 * thg1max;
	long double l1k2_max = param->en * eg1 - l1 * eg1 * thg1min;

	long double l1k2 = log(l1k2_min) + (log(l1k2_max) - log(l1k2_min)) * xx[4];

	l1k2 = exp(l1k2);
	long double Jacobian_thg1 = l1k2 / (l1 * eg1);

	thg1 = (param->en * eg1 - l1k2) / (l1 * eg1);

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
	if (thgmin > thgmax) return 0;

	//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
					*(param->max[param->E_prime]-param->min[param->E_prime])
					*(param->max[param->cos_thl]-param->min[param->cos_thl])
					*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
					*(thgmax-thgmin)*(log(l1k2_max)-log(l1k2_min))*Jacobian_thg1
					*2.*(phig1max-phig1min);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
					*(1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))/d2;

	return 0;
}

int Integrands::integrand_brems_2nd_add_l2k1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[2];

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[3];
	eg = exp(eg);
	Jacobian_eg = eg;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

//	generates number for cos(theta'_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[4];

	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double l2k1_min = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmax_cos
			- l2*eg*cos(thl)*cos(thg);
	long double l2k1_max = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmin_cos
			- l2*eg*cos(thl)*cos(thg);

	long double l2k1 = log(l2k1_min) + (log(l2k1_max) - log(l2k1_min)) * xx[5];

	l2k1 = exp(l2k1);

	phig = (en1*eg - l2*eg*cos(thl)*cos(thg) - l2k1)/(l2*eg*sin(thl)*sin(thg));

	if (pow(phig, 2.) >= 1.) return 0;

	phig = acos(phig);
	long double Jacobian_phig = l2k1 / (l2*eg*sin(thl)*sin(thg)) / sin(phig);
	if (rand.uniform() < 0.5) phig = 2.*pi - phig;

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add_phig1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thgmax-thgmin)*(thg1max-thg1min)
			*2.*(log(l2k1_max)-log(l2k1_min))*Jacobian_phig;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
//	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))/d3;

	return 0;
}

int Integrands::integrand_brems_2nd_add_l2k2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[3];

	thg1min = -1.;
	thg1max = 1.;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	long double l2k2_min = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1max_cos
			- l2*eg1*cos(thl)*cos(thg1);
	long double l2k2_max = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1min_cos
			- l2*eg1*cos(thl)*cos(thg1);

	long double l2k2 = log(l2k2_min) + (log(l2k2_max) - log(l2k2_min)) * xx[5];

	l2k2 = exp(l2k2);

	phig1 = (en1*eg1 - l2*eg1*cos(thl)*cos(thg1) - l2k2)/(l2*eg1*sin(thl)*sin(thg1));

	if (pow(phig1, 2.) >= 1.) return 0;

	phig1 = acos(phig1);
	long double Jacobian_phig1 = l2k2 / (l2*eg1*sin(thl)*sin(thg1)) / sin(phig1);
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (param->min[param->cos_thg] > thgmin) thgmin = param->min[param->cos_thg];
	if (param->max[param->cos_thg] < thgmax) thgmax = param->max[param->cos_thg];
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thgmax-thgmin)*(thg1max-thg1min)
			*2.*(log(l2k2_max)-log(l2k2_min))*Jacobian_phig1;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
//	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = l1k1;
	long double d2 = l1k2;
	long double d3 = l2k1;
	long double d4 = l2k2;

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))/d4;

	return 0;
}

int Integrands::integrand_brems_2nd_thl (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

	//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[3];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[4];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin > thgmax) return 0;

	//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[5];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
					*(param->max[param->E_prime]-param->min[param->E_prime])
					*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
					*(thgmax-thgmin)*(thg1max-thg1min)*2.*(phig1max-phig1min);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
					+ cos(thg1)*cos(thg);
	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);
	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	if (std::abs(l1k1) < 1e-60 || std::abs(l1k2) < 1e-60
			|| std::abs(l2k1) < 1e-60 || std::abs(l2k2) < 1e-60) return 0;
	if (std::abs(l1k1+l1k2-k1k2) < 1e-60 || std::abs(l2k1+l2k2+k1k2) < 1e-60) return 0;

	return 0;
}

int Integrands::integrand_brems_2nd_thl_1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

	//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[3];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[4];

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin > thgmax) return 0;

	//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[5];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
					*(param->max[param->E_prime]-param->min[param->E_prime])
					*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
					*(thgmax-thgmin)*(thg1max-thg1min)*(phig1max-phig1min);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
					+ cos(thg1)*cos(thg);
	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);
	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	if (std::abs(l1k1) < 1e-60 || std::abs(l1k2) < 1e-60
			|| std::abs(l2k1) < 1e-60 || std::abs(l2k2) < 1e-60) return 0;
	if (std::abs(l1k1+l1k2-k1k2) < 1e-60 || std::abs(l2k1+l2k2+k1k2) < 1e-60) return 0;

	return 0;
}

int Integrands::integrand_brems_2nd_thl_2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

	//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[3];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[4];
	phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin > thgmax) return 0;

	//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[5];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
					*(param->max[param->E_prime]-param->min[param->E_prime])
					*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
					*(thgmax-thgmin)*(thg1max-thg1min)*(phig1max-phig1min);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
					+ cos(thg1)*cos(thg);
	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);
	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	if (std::abs(l1k1) < 1e-60 || std::abs(l1k2) < 1e-60
			|| std::abs(l2k1) < 1e-60 || std::abs(l2k2) < 1e-60) return 0;
	if (std::abs(l1k1+l1k2-k1k2) < 1e-60 || std::abs(l2k1+l2k2+k1k2) < 1e-60) return 0;

	return 0;
}

int Integrands::integrand_brems_2nd_add_thl (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma];

	//generates photon energy
	eg1 = eg1max*xx[2];

	thg1min = -1.;
	thg1max = 1.;

//	thg1min = -1.;
//	thg1max = 1.;

	long double h1 = pow(eg,2.) * l2 * eg1 * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg * l2 * sin(thl);
	long double h4 = l2 * eg1 * sin(thl);
	long double h5 = l1 * eg1 - l2 * eg1 * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg - l2 * eg * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[3];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	 phig1max = acos(phig1min_cos);
	 phig1min = acos(phig1max_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	phig1 = phig1min+(phig1max-phig1min)*xx[4];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[5];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);
	long double d5 = pow(l1k1 + l1k2 - k1k2,2.);
	long double d6 = pow(l2k1 + l2k2 + k1k2,2.);

	if (d1*d3*d5 < 1e-60) return 0;
	if (d2*d4*d6 < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thgmax-thgmin)*(thg1max-thg1min)*2.*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_brems_2nd_add_thl_1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma];

	//generates photon energy
	eg1 = eg1max*xx[2];

	thg1min = -1.;
	thg1max = 1.;

//	thg1min = -1.;
//	thg1max = 1.;

	long double h1 = pow(eg,2.) * l2 * eg1 * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg * l2 * sin(thl);
	long double h4 = l2 * eg1 * sin(thl);
	long double h5 = l1 * eg1 - l2 * eg1 * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg - l2 * eg * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[3];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	 phig1max = acos(phig1min_cos);
	 phig1min = acos(phig1max_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	phig1 = phig1min+(phig1max-phig1min)*xx[4];

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[5];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);
	long double d5 = pow(l1k1 + l1k2 - k1k2,2.);
	long double d6 = pow(l2k1 + l2k2 + k1k2,2.);

	if (d1*d3*d5 < 1e-60) return 0;
	if (d2*d4*d6 < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thgmax-thgmin)*(thg1max-thg1min)*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_brems_2nd_add_thl_2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma];

	//generates photon energy
	eg1 = eg1max*xx[2];

	thg1min = -1.;
	thg1max = 1.;

//	thg1min = -1.;
//	thg1max = 1.;

	long double h1 = pow(eg,2.) * l2 * eg1 * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg * l2 * sin(thl);
	long double h4 = l2 * eg1 * sin(thl);
	long double h5 = l1 * eg1 - l2 * eg1 * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg - l2 * eg * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[3];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	 phig1max = acos(phig1min_cos);
	 phig1min = acos(phig1max_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	phig1 = phig1min+(phig1max-phig1min)*xx[4];
	phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[5];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);
	long double d5 = pow(l1k1 + l1k2 - k1k2,2.);
	long double d6 = pow(l2k1 + l2k2 + k1k2,2.);

	if (d1*d3*d5 < 1e-60) return 0;
	if (d2*d4*d6 < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thgmax-thgmin)*(thg1max-thg1min)*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_brems_2nd_add_2diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma];

	//generates photon energy
	eg1 = eg1max*xx[2];

	thg1min = -1.;
	thg1max = 1.;

	thg1 = thg1min+(thg1max-thg1min)*xx[3];
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	 phig1max = acos(phig1min_cos);
	 phig1min = acos(phig1max_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	phig1 = phig1min+(phig1max-phig1min)*xx[4];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	if (thgmin > thgmax) return 0;

	thg = param->thg;
	if (cos(thg) < thgmin || cos(thg) > thgmax) return 0;

//	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
//	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
//	if (thgmin < -1.) thgmin = -1.;
//	if (thgmax > 1.) thgmax = 1.;

//	generates number for cos(theta_gamma)
//	thg = thgmin+(thgmax-thgmin)*xx[4];

	//check if -1<cos(thg)<1
//	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);
	long double d5 = pow(l1k1 + l1k2 - k1k2,2.);
	long double d6 = pow(l2k1 + l2k2 + k1k2,2.);

	if (d1*d3*d5 < 1e-60) return 0;
	if (d2*d4*d6 < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thg1max-thg1min)*2.*(phig1max-phig1min);

//	std::cout << (log(eg1max)-log(eg1min))*Jacobian_eg1 << "\n";

//	(param->max[param->E_prime]-param->min[param->E_prime])
//						*(log(eg1max)-log(eg1min))*Jacobian_eg1*egmax
//						*(thgmax-thgmin)*2.*(phig1max-phig1min)

	return 0;
}

int Integrands::integrand_brems_2nd_add_l1k1_2diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[1];

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

//	thgmin = param->min[param->cos_thg];
//	thgmax = param->max[param->cos_thg];
//
//	long double l1k1_min = param->en * eg - l1 * eg * thgmax;
//	long double l1k1_max = param->en * eg - l1 * eg * thgmin;
//
//	long double l1k1 = log(l1k1_min) + (log(l1k1_max) - log(l1k1_min)) * xx[4];
//
//	l1k1 = exp(l1k1);
//	long double Jacobian_thg = l1k1 / (l1 * eg);
//
//	thg = (param->en * eg - l1k1) / (l1 * eg);
//
//	if (pow(thg, 2.) >= 1.) return 0;

	thg = param->thg;
	long double l1k1 = param->en * eg - l1 * eg * cos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double phigmax = acos(phigmin_cos);
	long double phigmin = acos(phigmax_cos);

	phig = phigmin+(phigmax-phigmin)*xx[3];
	if (rand.uniform() < 0.5) phig = 2.*pi - phig;

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add_phig1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thg1max-thg1min)*2.*(phigmax-phigmin);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))/d1;
//	ff[0] /= 4.;

	return 0;
}

int Integrands::integrand_brems_2nd_add_l1k2_2diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[2];

	thg1min = -1.;
	thg1max = 1.;

	//	generates number for cos(theta'_gamma)
	//	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	long double l1k2_min = param->en * eg1 - l1 * eg1 * thg1max;
	long double l1k2_max = param->en * eg1 - l1 * eg1 * thg1min;

	long double l1k2 = log(l1k2_min) + (log(l1k2_max) - log(l1k2_min)) * xx[3];

	l1k2 = exp(l1k2);
	long double Jacobian_thg1 = l1k2 / (l1 * eg1);

	thg1 = (param->en * eg1 - l1k2) / (l1 * eg1);

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[4];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	if (thgmin > thgmax) return 0;

	thg = param->thg;
	if (cos(thg) < thgmin || cos(thg) > thgmax) return 0;

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
					*(param->max[param->E_prime]-param->min[param->E_prime])
					*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
					*(log(l1k2_max)-log(l1k2_min))*Jacobian_thg1
					*2.*(phig1max-phig1min);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
					*(1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))/d2;
//	ff[0] /= 4.;

	return 0;
}

int Integrands::integrand_brems_2nd_add_l2k1_2diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = param->thl;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[1];

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

//	thgmin = param->min[param->cos_thg];
//	thgmax = param->max[param->cos_thg];

//	generates number for cos(theta'_gamma)
//	thg = thgmin+(thgmax-thgmin)*xx[4];

	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = param->thg;

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double l2k1_min = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmax_cos
			- l2*eg*cos(thl)*cos(thg);
	long double l2k1_max = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmin_cos
			- l2*eg*cos(thl)*cos(thg);

	long double l2k1 = log(l2k1_min) + (log(l2k1_max) - log(l2k1_min)) * xx[3];

	l2k1 = exp(l2k1);

	phig = (en1*eg - l2*eg*cos(thl)*cos(thg) - l2k1)/(l2*eg*sin(thl)*sin(thg));

	if (pow(phig, 2.) >= 1.) return 0;

	phig = acos(phig);
	long double Jacobian_phig = l2k1 / (l2*eg*sin(thl)*sin(thg)) / sin(phig);
	if (rand.uniform() < 0.5) phig = 2.*pi - phig;

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (param->min[param->cos_thg] > thg1min) thg1min = -1.;
	if (param->max[param->cos_thg] < thg1max) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add_phig1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thg1max-thg1min)*2.*(log(l2k1_max)-log(l2k1_min))*Jacobian_phig;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
//	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))/d3;
//	ff[0] /= 4.;

	return 0;
}

int Integrands::integrand_brems_2nd_add_l2k2_2diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];
	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma_prime];

	//generates photon energy
	eg1 = eg1max*xx[2];

	thg1min = -1.;
	thg1max = 1.;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[3];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	long double l2k2_min = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1max_cos
			- l2*eg1*cos(thl)*cos(thg1);
	long double l2k2_max = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1min_cos
			- l2*eg1*cos(thl)*cos(thg1);

	long double l2k2 = log(l2k2_min) + (log(l2k2_max) - log(l2k2_min)) * xx[4];

	l2k2 = exp(l2k2);

	phig1 = (en1*eg1 - l2*eg1*cos(thl)*cos(thg1) - l2k2)/(l2*eg1*sin(thl)*sin(thg1));

	if (pow(phig1, 2.) >= 1.) return 0;

	phig1 = acos(phig1);
	long double Jacobian_phig1 = l2k2 / (l2*eg1*sin(thl)*sin(thg1)) / sin(phig1);
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	if (thgmin > thgmax) return 0;

	thg = param->thg;
	if (cos(thg) < thgmin || cos(thg) > thgmax) return 0;

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thg1max-thg1min)*2.*(log(l2k2_max)-log(l2k2_min))*Jacobian_phig1;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
//	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))/d4;
//	ff[0] /= 4.;

	return 0;
}

int Integrands::integrand_brems_2nd_1diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[3];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1 = param->thg;

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(gamma2,2.) + pow(g3,2.)*gamma3;

	if ( gamma < 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

//	std::cout<<phig1 * 180. /pi<<"\n";

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta < 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
//
//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[5];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*2.*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_brems_2nd_2diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;


	egmin = param->min[param->E_gamma];

	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;


	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

//	thg1min = -1.;
//	thg1max = 1.;

//	generates number for cos(theta'_gamma)
//	thg1 = thg1min+(thg1max-thg1min)*xx[3];

//	if (pow(thg1, 2.) >= 1.) return 0;

//	thg1 = acos(thg1);
	thg1 = param->thg;

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(gamma2,2.) + pow(g3,2.)*gamma3;

	if ( gamma < 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

//	std::cout<<phig1 * 180. /pi<<"\n";

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta < 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
//
//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[4];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*2.*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_brems_2nd_l1k1_2diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[1];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

//	thgmin = param->min[param->cos_thg];
//	thgmax = param->max[param->cos_thg];

	thg = param->thg;

	long double h1 = pow(eg1,2.) * l2 * eg * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg1 * l2 * sin(thl);
	long double h4 = l2 * eg * sin(thl);
	long double h5 = l1 * eg - l2 * eg * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg1 - l2 * eg1 * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	long double l1k1_min = param->en * eg - l1 * eg * thgmax;
//	long double l1k1_max = param->en * eg - l1 * eg * thgmin;
//
//	long double l1k1 = log(l1k1_min) + (log(l1k1_max) - log(l1k1_min)) * xx[4];
//
//	l1k1 = exp(l1k1);
//	long double Jacobian_thg = l1k1 / (l1 * eg);
//
//	thg = (param->en * eg - l1k1) / (l1 * eg);

//	generates number for cos(theta'_gamma)
//	thg = thgmin+(thgmax-thgmin)*xx[4];

//	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

//	thg = acos(thg);

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double phigmax = acos(phigmin_cos);
	long double phigmin = acos(phigmax_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	phig = phigmin+(phigmax-phigmin)*xx[3];
	if (rand.uniform() < 0.5) phig = 2.*pi - phig;

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_phig1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thg1max-thg1min)*2.*(phigmax-phigmin);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))/d1;

//	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
//			*((1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))/d1
//			+ (1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))/d2
//			+ (1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))/d3
//			+ (1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))/d4);

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = en1;
		FS->theta_l = thl;
		FS->E_gamma = eg;
		FS->theta_gamma = thg;
		FS->phi_gamma = phig;
		FS->E_gamma_prime = eg1;
		FS->theta_gamma_prime = thg1;
		FS->phi_gamma_prime = phig1;
		FS->weight = weight[0]*weight[8];
		if (param->flag[param->asymmetry] == 1)
			FS->weight *= 1. + param->P*CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
									/CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1);

	};

	return 0;
}

int Integrands::integrand_brems_2nd_l1k2_2diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

	long double h1 = pow(eg,2.) * l2 * eg1 * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg * l2 * sin(thl);
	long double h4 = l2 * eg1 * sin(thl);
	long double h5 = l1 * eg1 - l2 * eg1 * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg - l2 * eg * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	generates number for cos(theta'_gamma)
//	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	long double l1k2_min = param->en * eg1 - l1 * eg1 * thg1max;
	long double l1k2_max = param->en * eg1 - l1 * eg1 * thg1min;

	long double l1k2 = log(l1k2_min) + (log(l1k2_max) - log(l1k2_min)) * xx[3];

	l1k2 = exp(l1k2);
	long double Jacobian_thg1 = l1k2 / (l1 * eg1);

	thg1 = (param->en * eg1 - l1k2) / (l1 * eg1);

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	 phig1max = acos(phig1min_cos);
	 phig1min = acos(phig1max_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	phig1 = phig1min+(phig1max-phig1min)*xx[4];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;
//	phig1 = acos(phig1);

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

//	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
//	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
//	thg = thgmin+(thgmax-thgmin)*xx[6];
//
//	//check if -1<cos(thg)<1
//	if (pow(thg, 2.) >= 1.) return 0;
//
//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;
//
//	thg = acos(thg);

	thg = param->thg;

	if (thgmax < cos(thg)) return 0;
	if (thgmin > cos(thg)) return 0;

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))
			*Jacobian_eg1*(log(l1k2_max)-log(l1k2_min))*Jacobian_thg1
			*2.*(phig1max-phig1min);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
//	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))/d2;

//	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
//			*((1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))/d1
//			+ (1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))/d2
//			+ (1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))/d3
//			+ (1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))/d4);

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = en1;
		FS->theta_l = thl;
		FS->E_gamma = eg;
		FS->theta_gamma = thg;
		FS->phi_gamma = phig;
		FS->E_gamma_prime = eg1;
		FS->theta_gamma_prime = thg1;
		FS->phi_gamma_prime = phig1;
		FS->weight = weight[0]*weight[8];
		if (param->flag[param->asymmetry] == 1)
			FS->weight *= 1. + param->P*CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
									/CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1);

	};

	return 0;
}

int Integrands::integrand_brems_2nd_l2k2_2diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

	long double h1 = pow(eg,2.) * l2 * eg1 * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg * l2 * sin(thl);
	long double h4 = l2 * eg1 * sin(thl);
	long double h5 = l1 * eg1 - l2 * eg1 * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg - l2 * eg * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[3];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

//	 phig1max = acos(phig1min_cos);
//	 phig1min = acos(phig1max_cos);

	long double l2k2_min = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1max_cos
			- l2*eg1*cos(thl)*cos(thg1);
	long double l2k2_max = en1*eg1 - l2*eg1*sin(thl)*sin(thg1)*phig1min_cos
			- l2*eg1*cos(thl)*cos(thg1);

	long double l2k2 = log(l2k2_min) + (log(l2k2_max) - log(l2k2_min)) * xx[4];

	l2k2 = exp(l2k2);

	phig1 = (en1*eg1 - l2*eg1*cos(thl)*cos(thg1) - l2k2)/(l2*eg1*sin(thl)*sin(thg1));

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	if (pow(phig1, 2.) >= 1.) return 0;

//	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	phig1 = acos(phig1);
	long double Jacobian_phig1 = l2k2 / (l2*eg1*sin(thl)*sin(thg1)) / sin(phig1);
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

//	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
//	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
//	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
//	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = param->thg;

	if (thgmax < cos(thg)) return 0;
	if (thgmin > cos(thg)) return 0;

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thg1max-thg1min)*2.*(log(l2k2_max)-log(l2k2_min))*Jacobian_phig1;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double l2k1 = en1*eg - eg*l2*cospsi;
//	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))/d4;

//	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
//			*((1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))/d1
//			+ (1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))/d2
//			+ (1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))/d3
//			+ (1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))/d4);

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = en1;
		FS->theta_l = thl;
		FS->E_gamma = eg;
		FS->theta_gamma = thg;
		FS->phi_gamma = phig;
		FS->E_gamma_prime = eg1;
		FS->theta_gamma_prime = thg1;
		FS->phi_gamma_prime = phig1;
		FS->weight = weight[0]*weight[8];
		if (param->flag[param->asymmetry] == 1)
			FS->weight *= 1. + param->P*CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
									/CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1);

	};

	return 0;
}

int Integrands::integrand_brems_2nd_l2k1_2diff (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = param->thl;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[1];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1 - eg1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

//	thgmin = param->min[param->cos_thg];
//	thgmax = param->max[param->cos_thg];

	long double h1 = pow(eg1,2.) * l2 * eg * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg1 * l2 * sin(thl);
	long double h4 = l2 * eg * sin(thl);
	long double h5 = l1 * eg - l2 * eg * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg1 - l2 * eg1 * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	generates number for cos(theta'_gamma)
//	thg = thgmin+(thgmax-thgmin)*xx[4];

//	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = param->thg;

	long double g1 = pow(eg1,2.)*l2*eg*sin(thl)*sin(thg);
	long double g2 = pow(eg*eg1*sin(thg),2.) + pow(l2*eg1*sin(thl),2.);
	long double g3 = l2 * eg * sin(thl) * sin(thg);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg*cos(thg) - l2*eg*cos(thg)*cos(thl);
	long double beta2 = l1*eg1 - l2*eg1*cos(thl) - eg*eg1*cos(thg);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	long double phigmin_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	long double phigmax_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phigmax_cos < -1. || phigmin_cos > 1.) return 0;

	if (phigmax_cos > 1.) phigmax_cos = 1.;
	if (phigmin_cos < -1.) phigmin_cos = -1.;

	long double l2k1_min = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmax_cos
			- l2*eg*cos(thl)*cos(thg);
	long double l2k1_max = en1*eg - l2*eg*sin(thl)*sin(thg)*phigmin_cos
			- l2*eg*cos(thl)*cos(thg);

	long double l2k1 = log(l2k1_min) + (log(l2k1_max) - log(l2k1_min)) * xx[3];

	l2k1 = exp(l2k1);

	phig = (en1*eg - l2*eg*cos(thl)*cos(thg) - l2k1)/(l2*eg*sin(thl)*sin(thg));

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	if (pow(phig, 2.) >= 1.) return 0;

//	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	phig = acos(phig);
	long double Jacobian_phig = l2k1 / (l2*eg*sin(thl)*sin(thg)) / sin(phig);
	if (rand.uniform() < 0.5) phig = 2.*pi - phig;

//	long double phigmax = acos(phigmin_cos);
//	long double phigmin = acos(phigmax_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

//	phig = phigmin+(phigmax-phigmin)*xx[5];
//	phig1 = acos(phig1);

	if (2.*g1 * cos(phig) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg*cos(thl)*cos(thg) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thg1min = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thg1max = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thg1min < -1.) thg1min = -1.;
	if (thg1max > 1.) thg1max = 1.;
	if (thg1min > thg1max) return 0;

//	generates number for cos(theta_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	//check if -1<cos(thg)<1
	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl)
			+ l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig)
			- l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig1 = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}
	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));
	if (x < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_phig1(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thg1max-thg1min)*2.*(log(l2k1_max)-log(l2k1_min))*Jacobian_phig;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg);
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
//	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);

	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
			*(1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))/d3;

//	ff[0] *= d1*d2*d3*d4 / (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)
//			*((1./(d1+d2) + 1./(d1+d3) + 1./(d1+d4))/d1
//			+ (1./(d2+d1) + 1./(d2+d3) + 1./(d2+d4))/d2
//			+ (1./(d3+d1) + 1./(d3+d2) + 1./(d3+d4))/d3
//			+ (1./(d4+d1) + 1./(d4+d2) + 1./(d4+d3))/d4);

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = en1;
		FS->theta_l = thl;
		FS->E_gamma = eg;
		FS->theta_gamma = thg;
		FS->phi_gamma = phig;
		FS->E_gamma_prime = eg1;
		FS->theta_gamma_prime = thg1;
		FS->phi_gamma_prime = phig1;
		FS->weight = weight[0]*weight[8];
		if (param->flag[param->asymmetry] == 1)
			FS->weight *= 1. + param->P*CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
									/CS.crsect_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1);

	};

	return 0;
}

int Integrands::integrand_interf_brems_1st (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//calculates the lower limit for cos(theta_l) with formula 4.5.1
	thlmin = (param->en*en1 - M*(param->en-en1) - m2)/(l1*l2);

	//check if the lower limit calculated is lower than the maximum limit given
	if (param->max[param->cos_thl] < thlmin) return 0;

	//if the lower limit calculated is lower than the given limit or -1 resets the lower limit
	if (thlmin < param->min[param->cos_thl]) thlmin = param->min[param->cos_thl];

	//generates number for cos(theta_l)
	thl = thlmin+(param->max[param->cos_thl]-thlmin)*xx[1];

	//check if -1<cos(thl)<1
	if (pow(thl, 2.) > 1.) return 0;

	thl = acos(thl);

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;


	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[3];

	//check if cos(theta_gamma) is inside the cut-off limits
	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;


	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.interf_brems_1st(en1,thl,eg,thg)*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-thlmin)*(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	return 0;
}

int Integrands::integrand_brems_2nd_pol_add (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma];

	//generates photon energy
	eg1 = eg1max*xx[3];

	thg1min = -1.;
	thg1max = 1.;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);

	if ( gamma <= 0.) return 0;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	 phig1max = acos(phig1min_cos);
	 phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}
	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);
	long double d5 = pow(l1k1 + l1k2 - k1k2,2.);
	long double d6 = pow(l2k1 + l2k2 + k1k2,2.);

	if (d1*d3*d5 < 1e-60) return 0;
	if (d2*d4*d6 < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_pol_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thgmax-thgmin)*(thg1max-thg1min)*2.*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_brems_2nd_pol_add_thl (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = 0.;
	eg1max = param->min[param->E_gamma];

	//generates photon energy
	eg1 = eg1max*xx[2];

	thg1min = -1.;
	thg1max = 1.;

//	thg1min = -1.;
//	thg1max = 1.;

	long double h1 = pow(eg,2.) * l2 * eg1 * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg * l2 * sin(thl);
	long double h4 = l2 * eg1 * sin(thl);
	long double h5 = l1 * eg1 - l2 * eg1 * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg - l2 * eg * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[3];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma <= 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

//	std::cout<<phig1min<<"\t"<<phig1max<<"\n";

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	 phig1max = acos(phig1min_cos);
	 phig1min = acos(phig1max_cos);

//	std::cout << phig1min <<"\t"<< phig1max << "\n";

	phig1 = phig1min+(phig1max-phig1min)*xx[4];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;
//	phig1 = acos(phig1);

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta <= 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin > thgmax) return 0;

//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[5];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

//	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}
	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

//	if (phig < 0.) phig = 2.*pi + phig;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-60) return 0;

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double d1 = pow(l1k1,2.);
	long double d2 = pow(l1k2,2.);
	long double d3 = pow(l2k1,2.);
	long double d4 = pow(l2k2,2.);
	long double d5 = pow(l1k1 + l1k2 - k1k2,2.);
	long double d6 = pow(l2k1 + l2k2 + k1k2,2.);

	if (d1*d3*d5 < 1e-60) return 0;
	if (d2*d4*d6 < 1e-60) return 0;

	//function to be calculated

	ff[0] = CS.crsect_brems_2nd_pol_add(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*eg1max
			*(thgmax-thgmin)*(thg1max-thg1min)*2.*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_shiftQ2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//calculates the lower limit for cos(theta_l) with formula 4.5.1
	thlmin = (param->en*en1 - M*(param->en-en1) - m2)/(l1*l2);

	//check if the lower limit calculated is lower than the maximum limit given
	if (cos(param->thl) < thlmin) return 0;

	thl = param->thl;

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < 0.) egmin = 0.;
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = egmin+(egmax-egmin)*xx[1];

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;


	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[2];

	//check if cos(theta_gamma) is inside the cut-off limits
	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;


	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double Q12 = 2.*M*(param->en-en1-eg);
	long double deltaQ2 = Q12 - param->Q2;

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st(en1,thl,eg,thg,phig)*(param->max[param->E_prime]-param->min[param->E_prime])
		    *(egmax-egmin)*(thgmax-thgmin)*deltaQ2;

	return 0;
}

int Integrands::integrand_brems_1st_thl (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
//	events_no++;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//calculates the lower limit for cos(theta_l) with formula 4.5.1
	thlmin = (param->en*en1 - M*(param->en-en1) - m2)/(l1*l2);

	//check if the lower limit calculated is lower than the maximum limit given
	if (cos(param->thl) < thlmin) return 0;

	thl = param->thl;

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;


	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[2];

	//check if cos(theta_gamma) is inside the cut-off limits
	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;


	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

//	acc_events++;
//	if (events_no > param->MAXEVAL_2nd - 10) std::cout << acc_events <<"\t"<< events_no <<"\n";

	//function to be calculated
	ff[0] = CS.crsect_brems_1st(en1,thl,eg,thg,phig)*(param->max[param->E_prime]-param->min[param->E_prime])
		    *(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	return 0;
}


int Integrands::integrand_brems_1st_hadr_thl (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	if (param->min[param->E_prime] >= param->en) return 0;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin >= thgmax) return 0;

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[2];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;

	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st_hadr(en1,thl,eg,thg,phig)*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	return 0;
}

int Integrands::integrand_brems_1st_hadr_interf_thl (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	if (param->min[param->E_prime] >= param->en) return 0;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin >= thgmax) return 0;

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[2];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;

	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st_hadr_interf(en1,thl,eg,thg,phig)*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	return 0;
}

int Integrands::integrand_brems_1st_thl_ps2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	events_no++;

	thgmin = param->min[param->cos_thg];
	thgmax = param->max[param->cos_thg];

	thg = thgmin+(thgmax-thgmin)*xx[0];
	thg = acos(thg);
//	long double Jacobian_thg = sin(thg);

	long double phigmin = 0.;
	long double phigmax = 2.*pi;

	phig = phigmin+(phigmax-phigmin)*xx[1];

	egmin = param->min[param->E_gamma];
	egmax = param->max[param->E_gamma];

//	if (egmax > M*(param->en-m)/(M+param->en-l1*cos(thg)))
//		egmax = M*(param->en-m)/(M+param->en-l1*cos(thg));
//	if (egmin > egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	thl = param->thl;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double A = l1*cos(thl)-eg*cospsi;
	long double B = param->en+M-eg;
	long double C = eg*(param->en+M-l1*cos(thg))-M*param->en-m2;

	if (m2*(pow(A,2.)-pow(B,2.))+pow(C,2.) <= 0.) return 0;

	if (std::abs(pow(A,2.)-pow(B,2.)) < 1.e-12) return 0;

	en1 = (B*C - A*sqrt(m2*(pow(A,2.)-pow(B,2.))+pow(C,2.)))/(pow(A,2.)-pow(B,2.));

	if (en1 < param->min[param->E_prime] || en1 > param->en - eg) return 0;

	l2 = sqrt(pow(en1,2.)-m2);

	long double test = m2 + M*(param->en-en1-eg) - param->en*(en1+eg) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + en1*eg - l2*eg*sin(thl)*sin(thg)*cos(phig) - l2*eg*cos(thl)*cos(thg);

	if (std::abs(test) > 1e-16) return 0;

//	acc_events++;
//	if (events_no > param->MAXEVAL_2nd - 10) std::cout << acc_events <<"\t"<< events_no <<"\n";

	//function to be calculated
	ff[0] = CS.crsect_brems_1st_ps2(en1,thl,eg,thg,phig)*Jacobian_eg
			*(thgmax-thgmin)*(phigmax-phigmin)*(log(egmax)-log(egmin));

	return 0;
}

int Integrands::integrand_interf_brems_1st_thl (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//calculates the lower limit for cos(theta_l) with formula 4.5.1
	thlmin = (param->en*en1 - M*(param->en-en1) - m2)/(l1*l2);

	//check if the lower limit calculated is lower than the maximum limit given
	if (cos(param->thl) < thlmin) return 0;

	thl = param->thl;

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[2];

	//check if cos(theta_gamma) is inside the cut-off limits
	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;

	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.interf_brems_1st(en1,thl,eg,thg)*(param->max[param->E_prime]-param->min[param->E_prime])
		    *(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	return 0;
}

int Integrands::integrand_interf_brems_2nd (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[3];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

	long double h1 = pow(eg,2.) * l2 * eg1 * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg * l2 * sin(thl);
	long double h4 = l2 * eg1 * sin(thl);
	long double h5 = l1 * eg1 - l2 * eg1 * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg - l2 * eg * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma < 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

//	std::cout<<phig1 * 180. /pi<<"\n";

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta < 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin >= thgmax) return 0;
//
//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}
	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));

	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(thg1max-thg1min)*2.*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_interf_brems_2nd_1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[3];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

	long double h1 = pow(eg,2.) * l2 * eg1 * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg * l2 * sin(thl);
	long double h4 = l2 * eg1 * sin(thl);
	long double h5 = l1 * eg1 - l2 * eg1 * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg - l2 * eg * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma < 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[5];

//	std::cout<<phig1 * 180. /pi<<"\n";

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta < 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin >= thgmax) return 0;
//
//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}
	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));

	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(thg1max-thg1min)*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_interf_brems_2nd_2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = acos(thl);

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[3];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

	long double h1 = pow(eg,2.) * l2 * eg1 * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg * l2 * sin(thl);
	long double h4 = l2 * eg1 * sin(thl);
	long double h5 = l1 * eg1 - l2 * eg1 * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg - l2 * eg * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[4];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma < 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[5];
	phig1 = 2.*pi - phig1;

//	std::cout<<phig1 * 180. /pi<<"\n";

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta < 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin >= thgmax) return 0;
//
//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[6];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}
	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));

	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-param->min[param->cos_thl])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(thg1max-thg1min)*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_interf_brems_2nd_thl (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

	long double h1 = pow(eg,2.) * l2 * eg1 * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg * l2 * sin(thl);
	long double h4 = l2 * eg1 * sin(thl);
	long double h5 = l1 * eg1 - l2 * eg1 * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg - l2 * eg * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[3];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma < 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[4];
	if (rand.uniform() < 0.5) phig1 = 2.*pi - phig1;

//	std::cout<<phig1 * 180. /pi<<"\n";

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta < 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin >= thgmax) return 0;
//
//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[5];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(thg1max-thg1min)*2.*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_interf_brems_2nd_thl_1 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

	long double h1 = pow(eg,2.) * l2 * eg1 * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg * l2 * sin(thl);
	long double h4 = l2 * eg1 * sin(thl);
	long double h5 = l1 * eg1 - l2 * eg1 * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg - l2 * eg * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[3];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma < 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[4];

//	std::cout<<phig1 * 180. /pi<<"\n";

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta < 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin >= thgmax) return 0;
//
//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[5];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(thg1max-thg1min)*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_interf_brems_2nd_thl_2 (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//generates number for cos(theta_l)
//	thl = param->min[param->cos_thl]+(param->max[param->cos_thl]-param->min[param->cos_thl])*xx[1];

	thl = param->thl;

	egmin = param->min[param->E_gamma];
	egmax = param->en - en1;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	eg1min = param->min[param->E_gamma_prime];
	eg1max = param->en - en1 - eg;

	if (eg1max <= eg1min) return 0;

	//generates photon energy
	eg1 = log(eg1min)+(log(eg1max)-log(eg1min))*xx[2];
	eg1 = exp(eg1);
	Jacobian_eg1 = eg1;

	thg1min = -1.;
	thg1max = 1.;

	long double h1 = pow(eg,2.) * l2 * eg1 * sin(thl);
	long double h2 = eg * eg1;
	long double h3 = eg * l2 * sin(thl);
	long double h4 = l2 * eg1 * sin(thl);
	long double h5 = l1 * eg1 - l2 * eg1 * cos(thl);
	long double h6 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
					+ eg*eg1 + l1*l2*cos(thl);
	long double h8 = l1 * eg - l2 * eg * cos(thl);

	long double eta2 = h4 * (h4*h2*h8 - h1*h5);
	long double eta3 = pow(h1,2.) + 2.*h1*h4*h6 + pow(h4,2.) * (pow(h8,2.) + pow(h2,2.) + pow(h3,2.));

	if (eta3 <= 0.) return 0;

//	generates number for cos(theta'_gamma)
	thg1 = thg1min+(thg1max-thg1min)*xx[3];

	if (pow(thg1, 2.) >= 1.) return 0;

	thg1 = acos(thg1);

	long double g1 = pow(eg,2.)*l2*eg1*sin(thl)*sin(thg1);
	long double g2 = pow(eg*eg1*sin(thg1),2.) + pow(l2*eg*sin(thl),2.);
	long double g3 = l2 * eg1 * sin(thl) * sin(thg1);
	long double g4 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + en1*(eg+eg1)
			+ eg*eg1 + l1*l2*cos(thl) + l1*eg1*cos(thg1) - l2*eg1*cos(thg1)*cos(thl);
	long double beta2 = l1*eg - l2*eg*cos(thl) - eg*eg1*cos(thg1);

	long double gamma2 = g1 + g3*g4;
	long double gamma3 = - pow(g4,2.) + pow(beta2,2.) + g2;
	long double gamma = pow(g1,2.) + 2.*g1*g3*g4 + pow(g3,2.)*(pow(beta2,2.)+g2);
//	long double gamma_test = pow(sin(thg1),2.)*(-2.*eta2*cos(thg1) + eta3);

//	std::cout<<gamma<<"\t"<<gamma_test<<"\n";

	if ( gamma < 0.) return 0;

//	if (param->flag[param->ps] == 1) FS->acc_ev ++;

	phig1min_cos = (gamma2 - sqrt(gamma)) / pow(g3,2.);
	phig1max_cos = (gamma2 + sqrt(gamma)) / pow(g3,2.);

	if (phig1max_cos < -1. || phig1min_cos > 1.) return 0;

	if (phig1max_cos > 1.) phig1max_cos = 1.;
	if (phig1min_cos < -1.) phig1min_cos = -1.;

	phig1max = acos(phig1min_cos);
	phig1min = acos(phig1max_cos);

	phig1 = phig1min+(phig1max-phig1min)*xx[4];
	phig1 = 2.*pi - phig1;

//	std::cout<<phig1 * 180. /pi<<"\n";

	if (2.*g1 * cos(phig1) + g2 < 0.) return 0;

	long double beta1 = sqrt(2.*g1 * cos(phig1) + g2);
	long double beta3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg1*cos(thl)*cos(thg1) + eg*eg1;
	long double beta = pow(beta1,2.) + pow(beta2,2.) - pow(beta3,2.);

	if (beta < 0.) return 0;

//	std::cout<<beta<<"\n";

	thgmin = (- beta2 * beta3 - beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
	thgmax = (- beta2 * beta3 + beta1 * sqrt(beta)) / (pow(beta1,2.) + pow(beta2,2.));
//
//	std::cout<<thgmin<<"\t"<<thgmax<<"\n";
//
//	if (thgmax < -1. || thgmin > 1.) return 0;
//
	if (thgmin < -1.) thgmin = -1.;
	if (thgmax > 1.) thgmax = 1.;
	if (thgmin >= thgmax) return 0;
//
//	generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[5];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg || eg1 != eg1 || thg1 != thg1 || phig1 != phig1) return 0;

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
	long double alpha3 = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg1*sin(thl)*sin(thg1)*cos(phig1) -
			l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) + eg*eg1 - eg*eg1*cos(thg)*cos(thg1);

//	if (pow(alpha1,2.) + pow(alpha2,2.) - pow(alpha3,2.) <= 0.) return 0;

	if ( alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) < -1. || alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.)) > 1.) return 0;

	if (alpha1 > 0.)
		phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1);
	else {
		if (alpha2 > 0.) phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) - pi;
		else phig = asin(alpha3 / sqrt(pow(alpha1,2.)+pow(alpha2,2.))) - atan(alpha2/alpha1) + pi;
	}

	long double test = m2 + M*(param->en-en1-eg-eg1) - param->en*(en1+eg+eg1) + l1*l2*cos(thl) +
			l1*eg*cos(thg) + l1*eg1*cos(thg1) + en1*(eg+eg1) - l2*eg*sin(thl)*sin(thg)*cos(phig) -
			l2*eg1*sin(thl)*sin(thg1)*cos(phig1) - l2*eg*cos(thl)*cos(thg) - l2*eg1*cos(thl)*cos(thg1) +
			eg*eg1*(1. - sin(thg)*cos(phig)*sin(thg1)*cos(phig1) - sin(thg)*sin(phig)*sin(thg1)*sin(phig1) -
					cos(thg)*cos(thg1));

//	if (std::abs(test) > 1e-14) std::cout<<test<<"\t"<<phig * 180./pi<<"\n";
//	if (phig * 180./pi > 90.) std::cout<<phig * 180./pi<<"\n";
//	std::cout<<test<<"\n";

	if (std::abs(test) > 1e-15) return 0;

	long double x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.interf_brems_2nd(en1,thl,eg,thg,phig,eg1,thg1,phig1)
			*(param->max[param->E_prime]-param->min[param->E_prime])
			*(log(egmax)-log(egmin))*Jacobian_eg*(log(eg1max)-log(eg1min))*Jacobian_eg1
			*(thgmax-thgmin)*(thg1max-thg1min)*(phig1max-phig1min);

	return 0;
}

int Integrands::integrand_born_carbon (const double xx[], double ff[], const double weight[])const {

	Q2 = 1./Q2min+(1./Q2max-1./Q2min)*xx[0];
	Q2 = 1./Q2;
	Jacobian = - pow(Q2, 2.);
	ff[0] = CS.crsect_born_carbon(Q2) * (1./Q2max-1./Q2min) * Jacobian;

	return 0;
}

int Integrands::integrand_elastic_carbon (const double xx[], double ff[], const double weight[])const {

	Q2 = 1./Q2min+(1./Q2max-1./Q2min)*xx[0];
	Q2 = 1./Q2;
	Jacobian = - pow(Q2, 2);

	ff[0] = CS.crsect_elastic_carbon(Q2) * (1./Q2max-1./Q2min) * Jacobian;

	if (param->flag[param->ps] == 1) {

		FS->E_prime_l = (param->en - Q2 / (2.*M));
		FS->theta_l = acos((-Q2 * (1. + param->en/M) +
				2.*pow(param->en,2.) - 2.*m2)/(2.*l1*sqrt(pow(param->en-Q2/(2.*M),2.)-m2)));
		FS->E_gamma = 0.;
		FS->theta_gamma = 0.;
		FS->phi_gamma = 0.;
		FS->E_gamma_prime = 0.;
		FS->theta_gamma_prime = 0.;
		FS->phi_gamma_prime = 0.;
		FS->weight = weight[0]*weight[2];
		FS->sigma_diff = ff[0];
		if (param->flag[param->asymmetry] == 1)
			FS->weight *= 1. + param->P*CS.interf_elastic_carbon(Q2)
									/CS.crsect_elastic_carbon(Q2);

	}
	return 0;
}

int Integrands::integrand_interf_born_carbon (const double xx[], double ff[], const double weight[])const {

		Q2 = Q2min + (Q2max - Q2min) * xx[0];
		ff[0] = CS.interf_born_carbon(Q2) * (Q2max - Q2min);
		return 0;
	}

int Integrands::integrand_interf_elastic_carbon (const double xx[], double ff[], const double weight[])const {

	Q2 = Q2min+(Q2max-Q2min)*xx[0];

//Function to be integrated
	ff[0] = CS.interf_elastic_carbon(Q2)*(Q2max-Q2min);

	return 0;
}

// Integrand for Sigma | E_gamma > Delta
int Integrands::integrand_brems_1st_carbon (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	if (param->min[param->E_prime] >= param->en) return 0;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thlmin = (param->en*en1 - M*(param->en-en1) - m2)/(l1*l2);

	//check if the lower limit calculated is lower than the maximum limit given
	if (param->max[param->cos_thl] < thlmin) return 0;

	//if the lower limit calculated is lower than the given limit or -1 resets the lower limit
	if (thlmin < param->min[param->cos_thl]) thlmin = param->min[param->cos_thl];

	//generates number for cos(theta_l)
	thl = thlmin+(param->max[param->cos_thl]-thlmin)*xx[1];

	//check if -1<cos(thl)<1
	if (pow(thl, 2.) > 1.) return 0;

	thl = acos(thl);

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin >= thgmax) return 0;

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[3];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;

	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st_carbon(en1,thl,eg,thg,phig)*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-thlmin)*(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	if (param->flag[param->ps] == 1) {

		FS->weight = weight[0]*weight[5];
			FS->E_prime_l = en1;
			FS->theta_l = thl;
			FS->E_gamma = eg;
			FS->theta_gamma = thg;
			FS->phi_gamma = phig;
			FS->E_gamma_prime = 0.;
			FS->theta_gamma_prime = 0.;
			FS->phi_gamma_prime = 0.;
			FS->sigma_diff = ff[0];
			if (param->flag[param->asymmetry] == 1)
				FS->weight *= 1. + param->P*CS.interf_brems_1st_carbon(en1,thl,eg,thg)
										/CS.crsect_brems_1st_carbon(en1,thl,eg,thg,phig);

//		}
	};

	return 0;
}

int Integrands::integrand_interf_brems_1st_carbon (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	//calculates the lower limit for cos(theta_l) with formula 4.5.1
	thlmin = (param->en*en1 - M*(param->en-en1) - m2)/(l1*l2);

	//check if the lower limit calculated is lower than the maximum limit given
	if (param->max[param->cos_thl] < thlmin) return 0;

	//if the lower limit calculated is lower than the given limit or -1 resets the lower limit
	if (thlmin < param->min[param->cos_thl]) thlmin = param->min[param->cos_thl];

	//generates number for cos(theta_l)
	thl = thlmin+(param->max[param->cos_thl]-thlmin)*xx[1];

	//check if -1<cos(thl)<1
	if (pow(thl, 2.) > 1.) return 0;

	thl = acos(thl);

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;


	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[3];

	//check if cos(theta_gamma) is inside the cut-off limits
	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;


	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.interf_brems_1st_carbon(en1,thl,eg,thg)*(param->max[param->E_prime]-param->min[param->E_prime])
			*(param->max[param->cos_thl]-thlmin)*(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	return 0;
}

int Integrands::integrand_brems_1st_carbon_thl (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
//	events_no++;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[2];

	//check if cos(theta_gamma) is inside the cut-off limits
	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;


	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st_carbon(en1,thl,eg,thg,phig)*(param->max[param->E_prime]-param->min[param->E_prime])
		    *(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	return 0;
}

int Integrands::integrand_interf_brems_1st_carbon_thl (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;
	//xmin and xmax are limits which are defined in the param
	//generates number for the energy of the scattered electron in the limits fixed by the param
	en1 = param->min[param->E_prime]+(param->max[param->E_prime]-param->min[param->E_prime])*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thl = param->thl;

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[1];
	eg = exp(eg);
	Jacobian_eg = eg;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[2];

	//check if cos(theta_gamma) is inside the cut-off limits
	if (thg < param->min[param->cos_thg] || thg > param->max[param->cos_thg]) return 0;

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;

	if (a/b >= 1.) return 0;
	if (a/b <= -1.) return 0;

	phig = acos(a/b);

	long double x = sin(thl) * sin(thg) * sin(phig);
	if (x < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.interf_brems_1st_carbon(en1,thl,eg,thg)*(param->max[param->E_prime]-param->min[param->E_prime])
		    *(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	return 0;
}

int Integrands::integrand_born_e (const double xx[], double ff[], const double weight[])const {

	Q2 = 1./Q2min+(1./Q2max-1./Q2min)*xx[0];
	Q2 = 1./Q2;
	Jacobian = - pow(Q2, 2.);
	ff[0] = CS.crsect_born_e(Q2) * (1./Q2max-1./Q2min) * Jacobian;

	return 0;
}

int Integrands::integrand_brems_1st_e (const double xx[], double ff[], const double weight[])const {

	ff[0] = 0.;

	//generates number for the energy of the scattered electron in the limits fixed by the param
	long double en1min = param->min[param->E_prime];
	long double en1max = param->max[param->E_prime];

	en1 = en1min+(en1max-en1min)*xx[0];
	l2 = sqrt(pow(en1, 2.) - m2);         //absolute value of the vector l2

	thlmin = (param->en*en1 - M*(param->en-en1) - m2)/(l1*l2);

	//check if the lower limit calculated is lower than the maximum limit given
	if (param->max[param->cos_thl] < thlmin) return 0;

	//if the lower limit calculated is lower than the given limit or -1 resets the lower limit
	if (param->min[param->cos_thl] > thlmin)
		thlmin = param->min[param->cos_thl];

	//generates number for cos(theta_l)
	thl = thlmin+(param->max[param->cos_thl]-thlmin)*xx[1];

	//check if -1<cos(thl)<1
	if (pow(thl, 2.) > 1.) return 0;

	thl = acos(thl);

	y = M*(param->en-en1)-param->en*en1+l1*l2*cos(thl)+m2;
	a1 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))+(param->en-en1+M);
	a2 = sqrt(pow(l1,2.)+pow(l2,2.)-2.*l1*l2*cos(thl))-(param->en-en1+M);

	egmin = y/a1;
	egmax = -y/a2;

	if (egmax > param->en + M - en1 - param->min[param->E_p])
		egmax = param->en + M - en1 - param->min[param->E_p];

	// if (egmax > param->max[param->E_gamma]) egmax = param->max[param->E_gamma];
	if (param->flag[param->Delta_cut] == 0) {
	if (egmin < param->min[param->E_gamma]) egmin = param->min[param->E_gamma];
	if (egmin >= egmax) return 0;
	}

	//generates photon energy
	eg = log(egmin)+(log(egmax)-log(egmin))*xx[2];
	eg = exp(eg);
	Jacobian_eg = eg;

  // long double Ee = param->en + M - en1 - eg;
  // if (Ee < 1.) return 0;

	f = eg*l2*sin(thl);
	d = eg*(l1-l2*cos(thl));
	c = M*(param->en-en1-eg)-param->en*en1-param->en*eg+en1*eg+l1*l2*cos(thl)+m2;
	n2 = pow(f,2.)+pow(d,2.)-pow(c,2.);

	//check if the square root of n2 is defined
	if (n2 <= 0.) return 0;

	thgmin = (-f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));
	thgmax = (f*sqrt(n2)-d*c)/(pow(d,2)+pow(f,2));

	if (thgmax > param->max[param->cos_thg]) thgmax = param->max[param->cos_thg];
	if (thgmin < param->min[param->cos_thg]) thgmin = param->min[param->cos_thg];
	if (thgmin >= thgmax) return 0;

	//generates number for cos(theta_gamma)
	thg = thgmin+(thgmax-thgmin)*xx[3];

	//check if -1<cos(thg)<1
	if (pow(thg, 2.) >= 1.) return 0;

	thg = acos(thg);

	if (param->flag[param->Delta_cut] == 1) {

		long double E_cm = sqrt(m2 + M2 + 2.*param->en*M);
		long double en_cm = (m2-M2+E_cm*E_cm)/E_cm/2.;
		long double l1_cm = sqrt(en_cm*en_cm - m2);
		long double Ep_cm = E_cm - en_cm;
		long double p1_cm = -l1_cm;
		long double eg_cm = eg*(l1*p1_cm*cos(thg) - param->en*p1_cm+l1_cm*M)
                        		  /(Ep_cm*l1_cm-en_cm*p1_cm);

		if (eg_cm < param->min[param->E_gamma] || eg_cm > param->max[param->E_gamma])
			return 0;
	}

	a = M*(param->en-en1-eg)-param->en*en1+l1*l2*cos(thl)-
			param->en*eg+l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);

	//check if the generated values are real numbers
	if (en1 != en1 || thl != thl || eg != eg || thg != thg) return 0;

	if (a/b > 1.) return 0;
	if (a/b < -1.) return 0;

	phig = acos(a/b);

	if (std::abs(sin(thl)*sin(thg)*sin(phig)) < 1e-40) return 0;

	//function to be calculated
	ff[0] = CS.crsect_brems_1st_e(en1,thl,eg,thg,phig)*(en1max-en1min)
			*(param->max[param->cos_thl]-thlmin)*(log(egmax)-log(egmin))*Jacobian_eg*(thgmax-thgmin);

	return 0;
}
