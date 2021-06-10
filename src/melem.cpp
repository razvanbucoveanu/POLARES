#include "melem.h"
#include <cmath>
#include "const.h"
#include <iostream>

//using namespace constants;
using namespace POLARES;

Melem::Melem():param(0), m(constants::me), m2(constants::me2),
		m4(constants::me4), m6(constants::me6), M(constants::mpr),
		M2(constants::mpr2), M4(constants::mpr4), alpha(constants::alpha) {

}

int Melem::set_param(const Parameters* param) {

	this->param = param;
	m = param->m;
	m2 = pow(m,2.);
	m4 = pow(m,4.);
	m6 = pow(m,6.);
	M = param->M;
	M2 = pow(M,2.);
	M4 = pow(M,4.);

	return 0;
}

long double Melem::melem2(const long double l1k1,const long double l1k2,const long double k1k2,const long double Q2e,
		const long double Q2h,const long double Q2k,const long double S,const long double Sk,const long double Sq2,
		const long double f1, const long double f2)const {

	long double p1p2 = M2 + Q2h/2.;
	long double G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
	long double G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
	long double G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;

	return 16*G1*m2*pow(l1k1,-2) + (16*G2*m2*M2 + 16*G3*m2*M2)*pow(l1k1,-2) -
			  16*G1*pow(l1k1,-1) + (-16*G2*M2 - 16*G3*M2)*pow(l1k1,-1) +
			  (8*G1*m6 - 4*G1*m4*Q2h)*pow(l1k1,-2)*pow(k1k2 - l1k1 - l1k2,-2) +
			  (-8*G1*l1k2*m2 - 16*G1*m4 + 4*G1*l1k2*Q2h + 8*G1*m2*Q2h)*pow(l1k1,-1)*
			   pow(k1k2 - l1k1 - l1k2,-2) +
			  (16*G1*m4 + 4*G1*m2*Q2h + 4*G1*m2*Q2k)*pow(l1k1,-2)*
			   pow(k1k2 - l1k1 - l1k2,-1) +
			  (8*G1*l1k2 - 48*G1*m2 - 4*G1*Q2e - 4*G1*Q2h - 4*G1*Q2k)*pow(l1k1,-1)*
			   pow(k1k2 - l1k1 - l1k2,-1) +
			  (-16*G2*m4 - 16*G3*m4 + 8*G2*l1k2*M2 + 8*G3*l1k2*M2 -
			     32*G2*m2*M2 - 32*G3*m2*M2 - 4*G2*M2*Q2e - 4*G3*M2*Q2e -
			     16*G2*m2*Q2h - 16*G3*m2*Q2h - 4*G2*M2*Q2h - 4*G3*M2*Q2h +
			     4*G2*m2*Q2k + 4*G3*m2*Q2k - 4*G2*Q2h*Q2k - 4*G3*Q2h*Q2k -
			     16*G2*m2*S - 16*G3*m2*S + 8*G2*Q2e*S + 8*G3*Q2e*S - 8*G2*Q2k*S -
			     8*G3*Q2k*S + 40*G2*m2*Sk + 40*G3*m2*Sk - 8*G2*Q2h*Sk - 8*G3*Q2h*Sk -
			     16*G2*S*Sk - 16*G3*S*Sk - 16*G2*m2*Sq2 - 16*G3*m2*Sq2)*pow(l1k1,-1)*
			   pow(k1k2 - l1k1 - l1k2,-1) + 8*G1*m2*pow(l1k2,-2) +
			  (8*G2*m2*M2 + 8*G3*m2*M2)*pow(l1k2,-2) +
			  (8*G1*m6 - 4*G1*m4*Q2h)*pow(k1k2 - l1k1 - l1k2,-2)*pow(l1k2,-2) +
			  (-8*G1*l1k1*m2 + 16*G1*m4 - 4*G1*m2*Q2e - 8*G1*m2*Q2h + 4*G1*m2*Q2k)*
			   pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-2) +
			  (8*G1*k1k2 - 16*G1*m2 + 4*G1*Q2h + 4*G1*Q2k)*pow(l1k1,-1)*pow(l1k2,-1) +
			  (8*G2*k1k2*M2 + 8*G3*k1k2*M2 - 8*G2*M2*Q2e - 8*G3*M2*Q2e +
			     4*G2*M2*Q2h + 4*G3*M2*Q2h - 4*G2*M2*Q2k - 4*G3*M2*Q2k -
			     8*G2*Q2h*S - 8*G3*Q2h*S - 16*G2*Q2k*S - 16*G3*Q2k*S - 32*G2*S*Sk -
			     32*G3*S*Sk - 16*G2*Q2k*Sq2 - 16*G3*Q2k*Sq2 - 32*G2*Sk*Sq2 -
			     32*G3*Sk*Sq2)*pow(l1k1,-1)*pow(l1k2,-1) +
			  (-8*G1*l1k1*m2 - 16*G1*m4 + 4*G1*l1k1*Q2h + 8*G1*m2*Q2h)*
			   pow(k1k2 - l1k1 - l1k2,-2)*pow(l1k2,-1) +
			  (16*G1*m6 - 8*G1*m4*Q2h)*pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-2)*
			   pow(l1k2,-1) + (-64*G1*m2 - 8*G1*Q2e + 4*G1*Q2h)*
			   pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1) +
			  (32*G1*m4 + 4*G1*m2*Q2e - 4*G1*m2*Q2h)*pow(l1k1,-1)*
			   pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1) +
			  (32*G1*m2 - 16*G1*Q2h)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2) +
			  (16*G1*k1k2 + 32*G1*m2 + 8*G1*Q2e - 16*G1*Q2h - 8*G1*Q2k)*pow(l1k1,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1) +
			  (16*G2*k1k2*M2 + 16*G3*k1k2*M2 + 24*G2*M2*Q2e + 24*G3*M2*Q2e -
			     16*G2*M2*Q2h - 16*G3*M2*Q2h + 8*G2*M2*Q2k + 8*G3*M2*Q2k -
			     8*G2*Q2h*Q2k - 8*G3*Q2h*Q2k - 16*G2*Q2h*S - 16*G3*Q2h*S +
			     32*G2*Q2k*S + 32*G3*Q2k*S - 16*G2*Q2h*Sk - 16*G3*Q2h*Sk + 64*G2*S*Sk +
			     64*G3*S*Sk + 32*G2*Q2k*Sq2 + 32*G3*Q2k*Sq2 + 64*G2*Sk*Sq2 +
			     64*G3*Sk*Sq2)*pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1) +
			  (64*G1*m2 + 32*G1*Q2e - 64*G1*Q2h)*pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1) +
			  16*G1*m2*pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-2)*
			   (-8*G2*l1k1*m2*M2 - 8*G3*l1k1*m2*M2 - 4*G2*m2*M2*Q2e -
			     4*G3*m2*M2*Q2e - 8*G2*m2*M2*Q2h - 8*G3*m2*M2*Q2h +
			     4*G2*m2*M2*Q2k + 4*G3*m2*M2*Q2k + 4*G2*m2*Q2h*Q2k +
			     4*G3*m2*Q2h*Q2k - 8*G2*m2*Q2h*S - 8*G3*m2*Q2h*S + 8*G2*m2*Q2h*Sk +
			     8*G3*m2*Q2h*Sk - 8*G2*m2*Q2h*Sq2 - 8*G3*m2*Q2h*Sq2 -
			     4*G2*m2*pow(Q2h,2) - 4*G3*m2*pow(Q2h,2)) +
			  64*G1*m2*pow(2*l1k1 + Q2e - Q2k,-2) +
			  (64*G2*m2*M2 + 64*G3*m2*M2)*pow(2*l1k1 + Q2e - Q2k,-2) +
			  (32*G1*m6 - 16*G1*m4*Q2h)*pow(l1k2,-2)*pow(2*l1k1 + Q2e - Q2k,-2) +
			  (-64*G1*m4 - 16*G1*m2*Q2h - 16*G1*m2*Q2k)*pow(l1k2,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-2) +
			  (128*G1*m6 - 64*G1*m4*Q2h)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
			   pow(2*l1k1 + Q2e - Q2k,-2) +
			  (128*G1*m4 + 32*G1*m2*Q2h + 32*G1*m2*Q2k)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-2) +
			  32*G1*pow(2*l1k1 + Q2e - Q2k,-1) +
			  (32*G2*M2 + 32*G3*M2)*pow(2*l1k1 + Q2e - Q2k,-1) +
			  (-64*G1*m2 - 16*G1*Q2k)*pow(l1k1,-1)*pow(2*l1k1 + Q2e - Q2k,-1) +
			  (-64*G2*m2*M2 - 64*G3*m2*M2 - 16*G2*M2*Q2k - 16*G3*M2*Q2k)*
			   pow(l1k1,-1)*pow(2*l1k1 + Q2e - Q2k,-1) +
			  (-16*G1*l1k2 - 32*G1*m2 + 16*G1*Q2h + 8*G1*Q2k)*
			   pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-1) +
			  (-16*G2*l1k2*M2 - 16*G3*l1k2*M2 - 16*G2*M2*Q2e - 16*G3*M2*Q2e +
			     16*G2*M2*Q2h + 16*G3*M2*Q2h - 8*G2*M2*Q2k - 8*G3*M2*Q2k +
			     8*G2*Q2h*Q2k + 8*G3*Q2h*Q2k - 16*G2*Q2h*S - 16*G3*Q2h*S +
			     32*G2*Q2k*S + 32*G3*Q2k*S + 16*G2*Q2h*Sk + 16*G3*Q2h*Sk + 64*G2*S*Sk +
			     64*G3*S*Sk - 16*G2*Q2h*Sq2 - 16*G3*Q2h*Sq2)*pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-1) +
			  (-16*G1*k1k2*m2 + 32*G1*m4 - 8*G1*m2*Q2e - 16*G1*m2*Q2h +
			     8*G1*m2*Q2k)*pow(l1k2,-2)*pow(2*l1k1 + Q2e - Q2k,-1) +
			  (-16*G2*k1k2*m2*M2 - 16*G3*k1k2*m2*M2 - 8*G2*m2*M2*Q2e -
			     8*G3*m2*M2*Q2e - 16*G2*m2*M2*Q2h - 16*G3*m2*M2*Q2h +
			     8*G2*m2*M2*Q2k + 8*G3*m2*M2*Q2k + 16*G2*m2*Q2h*S +
			     16*G3*m2*Q2h*S)*pow(l1k2,-2)*pow(2*l1k1 + Q2e - Q2k,-1) +
			  (-16*G1*k1k2 - 80*G1*m2 - 16*G1*Q2e - 8*G1*Q2h - 8*G1*Q2k)*pow(l1k2,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-1) +
			  (32*G1*m6 - 16*G1*m4*Q2h + 16*G1*m4*Q2k - 8*G1*m2*Q2h*Q2k)*
			   pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-1) +
			  (64*G1*k1k2*m2 + 128*G1*m4 - 32*G1*k1k2*Q2h - 64*G1*m2*Q2h)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*pow(2*l1k1 + Q2e - Q2k,-1) +
			  (-32*G1*k1k2 + 160*G1*m2 + 16*G1*Q2e + 32*G1*Q2h + 16*G1*Q2k)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-1) +
			  16*G1*m2*Q2e*pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-1) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-2)*(32*G1*m6 - 8*G1*m2*pow(Q2h,2))*
			   pow(2*l1k1 + Q2e - Q2k,-1) +
			  64*G1*m2*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) +
			  (64*G2*m2*M2 + 64*G3*m2*M2)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) +
			  (32*G1*m6 - 16*G1*m4*Q2h)*pow(l1k1,-2)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) +
			  (32*G1*l1k2*m2 - 64*G1*m4 + 16*G1*m2*Q2e + 16*G1*m2*Q2h)*pow(l1k1,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) +
			  (32*G2*l1k2*m2*M2 + 32*G3*l1k2*m2*M2 + 16*G2*m2*M2*Q2e +
			     16*G3*m2*M2*Q2e + 16*G2*m2*M2*Q2h + 16*G3*m2*M2*Q2h +
			     32*G2*m2*Q2h*S + 32*G3*m2*Q2h*S + 32*G2*m2*Q2h*Sq2 +
			     32*G3*m2*Q2h*Sq2)*pow(l1k1,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) +
			  (128*G1*m6 - 64*G1*m4*Q2h)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) +
			  (-64*G1*l1k2*m2 + 128*G1*m4 - 32*G1*m2*Q2e - 32*G1*m2*Q2h)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) +
			  pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   (-128*G1*m6 + 32*G1*m2*pow(Q2h,2))*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) +
			  pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   (-64*G2*l1k2*m2*M2 - 64*G3*l1k2*m2*M2 - 32*G2*m2*M2*Q2e -
			     32*G3*m2*M2*Q2e - 32*G2*m2*M2*Q2h - 32*G3*m2*M2*Q2h +
			     32*G2*m2*Q2h*Q2k + 32*G3*m2*Q2h*Q2k + 64*G2*m2*Q2h*S +
			     64*G3*m2*Q2h*S + 64*G2*m2*Q2h*Sk + 64*G3*m2*Q2h*Sk -
			     32*G2*m2*pow(Q2h,2) - 32*G3*m2*pow(Q2h,2))*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) +
			  (-32*G1*m4 - 8*G1*m2*Q2h - 8*G1*m2*Q2k)*pow(l1k1,-2)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  (16*G1*l1k2 + 96*G1*m2 + 16*G1*Q2e + 16*G1*Q2k)*pow(l1k1,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  (32*G2*m4 + 32*G3*m4 + 16*G2*l1k2*M2 + 16*G3*l1k2*M2 +
			     64*G2*m2*M2 + 64*G3*m2*M2 + 16*G2*M2*Q2e + 16*G3*M2*Q2e +
			     32*G2*m2*Q2h + 32*G3*m2*Q2h - 8*G2*m2*Q2k - 8*G3*m2*Q2k +
			     8*G2*M2*Q2k + 8*G3*M2*Q2k + 32*G2*m2*S + 32*G3*m2*S -
			     16*G2*Q2e*S - 16*G3*Q2e*S + 16*G2*Q2h*S + 16*G3*Q2h*S + 16*G2*Q2k*S +
			     16*G3*Q2k*S - 80*G2*m2*Sk - 80*G3*m2*Sk + 32*G2*S*Sk + 32*G3*S*Sk +
			     32*G2*m2*Sq2 + 32*G3*m2*Sq2 + 16*G2*Q2h*Sq2 + 16*G3*Q2h*Sq2)*
			   pow(l1k1,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  16*G1*m2*pow(k1k2 - l1k1 - l1k2,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) -
			  8*G1*m2*Q2e*pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  (64*G1*m2 + 32*G1*Q2e)*pow(l1k2,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  (64*G2*m2*M2 + 64*G3*m2*M2 + 32*G2*M2*Q2e + 32*G3*M2*Q2e)*
			   pow(l1k2,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  (-32*G1*m6 - 16*G1*m4*Q2e + 16*G1*m4*Q2k + 8*G1*m2*Q2e*Q2k)*
			   pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  (-224*G1*m2 - 32*G1*Q2e)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
			   (-64*G1*l1k2*m2 - 128*G1*m4 + 32*G1*l1k2*Q2h + 96*G1*m2*Q2h -
			     32*G1*m2*Q2k + 16*G1*Q2h*Q2k - 16*G1*pow(Q2h,2))*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  (-32*G1*l1k2 + 64*G1*m2 - 32*G1*Q2k)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  (-32*G2*l1k2*M2 - 32*G3*l1k2*M2 + 32*G2*M2*Q2e + 32*G3*M2*Q2e -
			     32*G2*Q2h*S - 32*G3*Q2h*S - 64*G2*Q2k*S - 64*G3*Q2k*S - 128*G2*S*Sk -
			     128*G3*S*Sk - 32*G2*Q2h*Sq2 - 32*G3*Q2h*Sq2)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  (-256*G1*m6 + 128*G1*m4*Q2h)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
			   pow(2*l1k1 + Q2e - Q2k,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  (-256*G1*m4 - 32*G1*m2*Q2e + 32*G1*m2*Q2h)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  (128*G1*m6 - 64*G1*m4*Q2h + 64*G1*m4*Q2k - 32*G1*m2*Q2h*Q2k)*
			   pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  (128*G1*m6 + 64*G1*m4*Q2e - 64*G1*m4*Q2k - 32*G1*m2*Q2e*Q2k)*
			   pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
			  pow(l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (192*G1*m4 + 48*G1*m2*Q2e - 16*G1*m2*Q2h - 8*G1*Q2e*Q2k -
			     8*G1*pow(Q2e,2) - 16*G1*pow(Q2k,2)) +
			  pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (64*G1*m4 - 16*G1*m2*Q2e - 16*G1*m2*Q2h - 16*G1*Q2e*Q2h +
			     8*G1*Q2e*Q2k + 64*G1*Q2h*Q2k - 8*G1*pow(Q2e,2) - 64*G1*pow(Q2h,2) -
			     16*G1*pow(Q2k,2)) + pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (128*G1*m4 + 64*G1*m2*Q2e - 48*G1*m2*Q2h - 16*G1*m2*Q2k -
			     8*G1*Q2e*Q2k - 8*G1*pow(Q2k,2)) +
			  pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-1)*(64*G1*m4 - 32*G1*m2*Q2h + 16*G1*k1k2*Q2k -
			     16*G1*m2*Q2k - 8*G1*Q2h*Q2k - 8*G1*pow(Q2k,2)) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (32*G1*m4 - 8*G1*m2*Q2e - 8*G1*m2*Q2h - 8*G1*Q2e*Q2h + 4*G1*Q2e*Q2k +
			     32*G1*Q2h*Q2k - 4*G1*pow(Q2e,2) - 32*G1*pow(Q2h,2) - 8*G1*pow(Q2k,2)) +
			  pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (32*G1*m4 - 16*G1*m2*Q2h + 8*G1*l1k2*Q2k - 8*G1*m2*Q2k -
			     4*G1*Q2h*Q2k - 4*G1*pow(Q2k,2)) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (-32*G1*m4 - 16*G1*m2*Q2e - 8*G1*m2*Q2h + 8*G1*m2*Q2k -
			     4*G1*Q2h*Q2k - 4*G1*pow(Q2e,2) - 4*G1*pow(Q2h,2) - 4*G1*pow(Q2k,2)) +
			  pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   (-64*G1*m4 - 32*G1*m2*Q2e + 24*G1*m2*Q2h + 8*G1*m2*Q2k +
			     4*G1*Q2e*Q2k + 4*G1*pow(Q2k,2)) +
			  pow(l1k1,-1)*pow(l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (-32*G1*m4 + 16*G1*m2*Q2h + 8*G1*k1k2*Q2k + 8*G1*m2*Q2k +
			     4*G1*Q2h*Q2k + 4*G1*pow(Q2k,2)) +
			  pow(l1k1,-1)*pow(l1k2,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-96*G1*m4 - 24*G1*m2*Q2e + 8*G1*m2*Q2h + 4*G1*Q2e*Q2k +
			     4*G1*pow(Q2e,2) + 8*G1*pow(Q2k,2)) +
			  pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (64*G1*m4 + 32*G1*m2*Q2e + 16*G1*m2*Q2h - 16*G1*m2*Q2k +
			     8*G1*Q2h*Q2k + 8*G1*pow(Q2e,2) + 8*G1*pow(Q2h,2) + 8*G1*pow(Q2k,2)) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   (-64*G1*m4 + 16*G1*m2*Q2e + 8*G1*m2*Q2h + 8*G1*Q2e*Q2h -
			     8*G1*m2*Q2k - 4*G1*Q2e*Q2k - 32*G1*Q2h*Q2k + 8*G1*pow(Q2e,2) +
			     48*G1*pow(Q2h,2) + 12*G1*pow(Q2k,2)) +
			  pow(l1k1,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-64*G1*m4 + 32*G1*m2*Q2h + 16*G1*l1k2*Q2k + 16*G1*m2*Q2k +
			     16*G1*pow(Q2k,2)) + pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-128*G1*m4 + 32*G1*m2*Q2e + 16*G1*m2*Q2h + 16*G1*Q2e*Q2h -
			     16*G1*m2*Q2k - 8*G1*Q2e*Q2k - 64*G1*Q2h*Q2k + 16*G1*pow(Q2e,2) +
			     96*G1*pow(Q2h,2) + 24*G1*pow(Q2k,2)) +
			  pow(l1k1,-2)*pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*(-32*G1*m6 + 8*G1*m2*pow(Q2k,2)) +
			  pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-2)*(-128*G1*m6 + 32*G1*m2*pow(Q2k,2)) +
			  pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (64*G1*m6 + 32*G1*m4*Q2e - 24*G1*m2*Q2e*Q2k - 4*G1*Q2k*pow(Q2e,2) +
			     8*G1*m2*pow(Q2k,2) - 4*G1*pow(Q2k,3)) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (-64*G1*m6 - 32*G1*m4*Q2e + 24*G1*m2*Q2e*Q2k + 4*G1*Q2k*pow(Q2e,2) -
			     8*G1*m2*pow(Q2k,2) + 4*G1*pow(Q2k,3)) +
			  pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (-64*G1*m6 + 64*G1*m4*Q2h - 96*G1*m4*Q2k + 24*G1*m2*Q2h*Q2k -
			     8*G1*m2*pow(Q2h,2) + 8*G1*pow(Q2k,3)) +
			  pow(l1k1,-1)*pow(l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-64*G1*m6 + 64*G1*m4*Q2h - 96*G1*m4*Q2k + 24*G1*m2*Q2h*Q2k -
			     8*G1*m2*pow(Q2h,2) + 8*G1*pow(Q2k,3)) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (64*G1*m6 + 128*G1*m4*Q2h - 96*G1*m4*Q2k + 24*G1*m2*Q2h*Q2k -
			     40*G1*m2*pow(Q2h,2) + 96*G1*Q2k*pow(Q2h,2) - 64*G1*pow(Q2h,3) -
			     48*G1*Q2h*pow(Q2k,2) + 8*G1*pow(Q2k,3)) +
			  pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-96*G2*m2*M2 - 96*G3*m2*M2 - 48*G2*M2*Q2e - 48*G3*M2*Q2e +
			     16*G2*M2*Q2h + 16*G3*M2*Q2h - 32*G2*Q2h*Sq2 - 32*G3*Q2h*Sq2 -
			     128*G2*S*Sq2 - 128*G3*S*Sq2 - 128*G2*pow(S,2) - 128*G3*pow(S,2)) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (32*G2*m4 + 32*G3*m4 + 16*G2*m2*M2 + 16*G3*m2*M2 +
			     8*G2*M2*Q2e + 8*G3*M2*Q2e + 32*G2*m2*Q2h + 32*G3*m2*Q2h +
			     8*G2*M2*Q2h + 8*G3*M2*Q2h - 8*G2*m2*Q2k - 8*G3*m2*Q2k +
			     8*G2*M2*Q2k + 8*G3*M2*Q2k + 32*G2*m2*S + 32*G3*m2*S -
			     16*G2*Q2e*S - 16*G3*Q2e*S - 16*G2*Q2k*S - 16*G3*Q2k*S - 80*G2*m2*Sk -
			     80*G3*m2*Sk - 32*G2*S*Sk - 32*G3*S*Sk + 32*G2*m2*Sq2 +
			     32*G3*m2*Sq2 - 64*G2*S*Sq2 - 64*G3*S*Sq2 - 64*G2*pow(S,2) -
			     64*G3*pow(S,2)) + pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   (-24*G2*m2*M2 - 24*G3*m2*M2 - 12*G2*M2*Q2e - 12*G3*M2*Q2e +
			     8*G2*M2*Q2h + 8*G3*M2*Q2h - 8*G2*Q2h*S - 8*G3*Q2h*S -
			     8*G2*Q2h*Sq2 - 8*G3*Q2h*Sq2 - 32*G2*S*Sq2 - 32*G3*S*Sq2 -
			     48*G2*pow(S,2) - 48*G3*pow(S,2)) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   (-32*G2*M2*Q2h - 32*G3*M2*Q2h + 32*G2*Q2h*Sq2 + 32*G3*Q2h*Sq2 +
			     128*G2*S*Sq2 + 128*G3*S*Sq2 - 16*G2*pow(Q2h,2) - 16*G3*pow(Q2h,2) +
			     128*G2*pow(S,2) + 128*G3*pow(S,2)) +
			  pow(k1k2 - l1k1 - l1k2,-2)*pow(l1k2,-1)*
			   (4*G2*l1k1*M2*Q2h + 4*G3*l1k1*M2*Q2h + 8*G2*m2*M2*Q2h +
			     8*G3*m2*M2*Q2h - 8*G2*l1k1*Q2h*S - 8*G3*l1k1*Q2h*S -
			     16*G2*m2*Q2h*S - 16*G3*m2*Q2h*S - 16*G2*l1k1*pow(S,2) -
			     16*G3*l1k1*pow(S,2) - 32*G2*m2*pow(S,2) - 32*G3*m2*pow(S,2)) +
			  pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-2)*
			   (4*G2*l1k2*M2*Q2h + 4*G3*l1k2*M2*Q2h + 8*G2*m2*M2*Q2h +
			     8*G3*m2*M2*Q2h - 8*G2*l1k2*Q2h*S - 8*G3*l1k2*Q2h*S -
			     16*G2*m2*Q2h*S - 16*G3*m2*Q2h*S - 16*G2*l1k2*pow(S,2) -
			     16*G3*l1k2*pow(S,2) - 32*G2*m2*pow(S,2) - 32*G3*m2*pow(S,2)) +
			  pow(l1k1,-2)*pow(k1k2 - l1k1 - l1k2,-2)*
			   (-4*G2*m4*M2*Q2h - 4*G3*m4*M2*Q2h + 8*G2*m4*Q2h*S +
			     8*G3*m4*Q2h*S + 16*G2*m4*pow(S,2) + 16*G3*m4*pow(S,2)) +
			  pow(k1k2 - l1k1 - l1k2,-2)*pow(l1k2,-2)*
			   (-4*G2*m4*M2*Q2h - 4*G3*m4*M2*Q2h + 8*G2*m4*Q2h*S +
			     8*G3*m4*Q2h*S + 16*G2*m4*pow(S,2) + 16*G3*m4*pow(S,2)) +
			  pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-2)*pow(l1k2,-1)*
			   (-8*G2*m4*M2*Q2h - 8*G3*m4*M2*Q2h + 16*G2*m4*Q2h*S +
			     16*G3*m4*Q2h*S + 32*G2*m4*pow(S,2) + 32*G3*m4*pow(S,2)) +
			  pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-128*G2*m4*M2 - 128*G3*m4*M2 - 32*G2*m2*M2*Q2e -
			     32*G3*m2*M2*Q2e - 64*G2*m4*Q2h - 64*G3*m4*Q2h +
			     32*G2*m2*M2*Q2h + 32*G3*m2*M2*Q2h + 16*G2*m2*Q2e*Q2h +
			     16*G3*m2*Q2e*Q2h - 48*G2*m2*Q2e*Q2k - 48*G3*m2*Q2e*Q2k +
			     48*G2*m2*Q2h*Q2k + 48*G3*m2*Q2h*Q2k - 128*G2*m2*Q2k*S -
			     128*G3*m2*Q2k*S + 128*G2*m4*Sk + 128*G3*m4*Sk - 128*G2*m2*Q2e*Sk -
			     128*G3*m2*Q2e*Sk + 128*G2*m2*Q2h*Sk + 128*G3*m2*Q2h*Sk -
			     192*G2*m2*Q2k*Sk - 192*G3*m2*Q2k*Sk - 256*G2*m2*S*Sk -
			     256*G3*m2*S*Sk - 64*G2*m2*Q2h*Sq2 - 64*G3*m2*Q2h*Sq2 +
			     64*G2*m2*Q2k*Sq2 + 64*G3*m2*Q2k*Sq2 + 256*G2*m2*S*Sq2 +
			     256*G3*m2*S*Sq2 + 128*G2*m2*Sk*Sq2 + 128*G3*m2*Sk*Sq2 -
			     48*G2*m2*pow(Q2k,2) - 48*G3*m2*pow(Q2k,2) + 256*G2*m2*pow(S,2) +
			     256*G3*m2*pow(S,2) - 256*G2*m2*pow(Sk,2) - 256*G3*m2*pow(Sk,2)) +
			  pow(l1k1,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (32*G2*m2*M2*Q2h + 32*G3*m2*M2*Q2h + 16*G2*l1k2*M2*Q2k +
			     16*G3*l1k2*M2*Q2k - 32*G2*m2*M2*Q2k - 32*G3*m2*M2*Q2k +
			     24*G2*m2*Q2h*Q2k + 24*G3*m2*Q2h*Q2k + 32*G2*m2*Q2h*S +
			     32*G3*m2*Q2h*S - 64*G2*m2*Q2k*S - 64*G3*m2*Q2k*S + 16*G2*Q2h*Q2k*S +
			     16*G3*Q2h*Q2k*S + 48*G2*m2*Q2h*Sk + 48*G3*m2*Q2h*Sk -
			     96*G2*m2*Q2k*Sk - 96*G3*m2*Q2k*Sk - 128*G2*m2*S*Sk -
			     128*G3*m2*S*Sk + 64*G2*Q2k*S*Sk + 64*G3*Q2k*S*Sk + 16*G2*Q2h*Q2k*Sq2 +
			     16*G3*Q2h*Q2k*Sq2 - 24*G2*m2*pow(Q2k,2) - 24*G3*m2*pow(Q2k,2) +
			     32*G2*S*pow(Q2k,2) + 32*G3*S*pow(Q2k,2) - 128*G2*m2*pow(S,2) -
			     128*G3*m2*pow(S,2) + 64*G2*Q2k*pow(S,2) + 64*G3*Q2k*pow(S,2) -
			     96*G2*m2*pow(Sk,2) - 96*G3*m2*pow(Sk,2)) +
			  pow(l1k1,-2)*pow(k1k2 - l1k1 - l1k2,-1)*
			   (16*G2*m4*M2 + 16*G3*m4*M2 + 8*G2*m4*Q2h + 8*G3*m4*Q2h +
			     4*G2*m2*M2*Q2h + 4*G3*m2*M2*Q2h + 4*G2*m2*M2*Q2k +
			     4*G3*m2*M2*Q2k - 8*G2*m2*Q2h*S - 8*G3*m2*Q2h*S - 16*G2*m2*Q2k*S -
			     16*G3*m2*Q2k*S - 16*G2*m4*Sk - 16*G3*m4*Sk - 8*G2*m2*Q2k*Sk -
			     8*G3*m2*Q2k*Sk - 32*G2*m2*S*Sk - 32*G3*m2*S*Sk -
			     32*G2*m2*pow(S,2) - 32*G3*m2*pow(S,2) - 16*G2*m2*pow(Sk,2) -
			     16*G3*m2*pow(Sk,2)) + pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-16*G2*l1k2*m2*M2 - 16*G3*l1k2*m2*M2 - 32*G2*m4*M2 -
			     32*G3*m4*M2 - 8*G2*m2*M2*Q2e - 8*G3*m2*M2*Q2e -
			     16*G2*m4*Q2k - 16*G3*m4*Q2k - 40*G2*m2*M2*Q2k -
			     40*G3*m2*M2*Q2k - 4*G2*M2*Q2e*Q2k - 4*G3*M2*Q2e*Q2k +
			     16*G2*m2*Q2k*S + 16*G3*m2*Q2k*S + 32*G2*m4*Sk + 32*G3*m4*Sk +
			     8*G2*m2*Q2k*Sk + 8*G3*m2*Q2k*Sk + 64*G2*m2*S*Sk + 64*G3*m2*S*Sk -
			     16*G2*Q2e*S*Sk - 16*G3*Q2e*S*Sk + 16*G2*Q2k*S*Sk + 16*G3*Q2k*S*Sk -
			     16*G2*m2*Q2k*Sq2 - 16*G3*m2*Q2k*Sq2 + 32*G2*Q2k*S*Sq2 +
			     32*G3*Q2k*S*Sq2 - 8*G2*m2*pow(Q2k,2) - 8*G3*m2*pow(Q2k,2) -
			     4*G2*M2*pow(Q2k,2) - 4*G3*M2*pow(Q2k,2) + 8*G2*S*pow(Q2k,2) +
			     8*G3*S*pow(Q2k,2) + 32*G2*l1k2*pow(S,2) + 32*G3*l1k2*pow(S,2) +
			     64*G2*m2*pow(S,2) + 64*G3*m2*pow(S,2) + 64*G2*Q2k*pow(S,2) +
			     64*G3*Q2k*pow(S,2) - 16*G2*m2*pow(Sk,2) - 16*G3*m2*pow(Sk,2)) +
			  pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   (16*G2*m4*M2 + 16*G3*m4*M2 + 4*G2*m2*M2*Q2e +
			     4*G3*m2*M2*Q2e + 8*G2*m4*Q2h + 8*G3*m4*Q2h - 4*G2*m2*M2*Q2h -
			     4*G3*m2*M2*Q2h - 2*G2*m2*Q2e*Q2h - 2*G3*m2*Q2e*Q2h +
			     6*G2*m2*Q2e*Q2k + 6*G3*m2*Q2e*Q2k - 6*G2*m2*Q2h*Q2k -
			     6*G3*m2*Q2h*Q2k - 16*G2*m2*Q2k*S - 16*G3*m2*Q2k*S - 16*G2*m4*Sk -
			     16*G3*m4*Sk + 16*G2*m2*Q2e*Sk + 16*G3*m2*Q2e*Sk - 16*G2*m2*Q2h*Sk -
			     16*G3*m2*Q2h*Sk + 24*G2*m2*Q2k*Sk + 24*G3*m2*Q2k*Sk -
			     32*G2*m2*S*Sk - 32*G3*m2*S*Sk + 8*G2*m2*Q2h*Sq2 + 8*G3*m2*Q2h*Sq2 -
			     24*G2*m2*Q2k*Sq2 - 24*G3*m2*Q2k*Sq2 - 32*G2*m2*S*Sq2 -
			     32*G3*m2*S*Sq2 - 48*G2*m2*Sk*Sq2 - 48*G3*m2*Sk*Sq2 +
			     6*G2*m2*pow(Q2k,2) + 6*G3*m2*pow(Q2k,2) - 32*G2*m2*pow(S,2) -
			     32*G3*m2*pow(S,2) + 32*G2*m2*pow(Sk,2) + 32*G3*m2*pow(Sk,2)) +
			  pow(l1k1,-2)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-32*G2*m4*M2 - 32*G3*m4*M2 - 16*G2*m4*Q2h - 16*G3*m4*Q2h -
			     8*G2*m2*M2*Q2h - 8*G3*m2*M2*Q2h - 8*G2*m2*M2*Q2k -
			     8*G3*m2*M2*Q2k - 8*G2*m2*Q2h*Q2k - 8*G3*m2*Q2h*Q2k -
			     16*G2*m2*Q2h*S - 16*G3*m2*Q2h*S + 32*G2*m2*Q2k*S + 32*G3*m2*Q2k*S +
			     32*G2*m4*Sk + 32*G3*m4*Sk - 16*G2*m2*Q2h*Sk - 16*G3*m2*Q2h*Sk +
			     16*G2*m2*Q2k*Sk + 16*G3*m2*Q2k*Sk + 64*G2*m2*S*Sk + 64*G3*m2*S*Sk +
			     64*G2*m2*pow(S,2) + 64*G3*m2*pow(S,2) + 32*G2*m2*pow(Sk,2) +
			     32*G3*m2*pow(Sk,2)) + pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-1)*(-16*G2*m2*M2*Q2h - 16*G3*m2*M2*Q2h +
			     8*G2*l1k2*M2*Q2k + 8*G3*l1k2*M2*Q2k + 16*G2*m2*M2*Q2k +
			     16*G3*m2*M2*Q2k - 4*G2*m2*Q2h*Q2k - 4*G3*m2*Q2h*Q2k -
			     4*G2*M2*Q2h*Q2k - 4*G3*M2*Q2h*Q2k + 16*G2*m2*Q2h*S +
			     16*G3*m2*Q2h*S + 32*G2*m2*Q2k*S + 32*G3*m2*Q2k*S + 8*G2*Q2h*Q2k*S +
			     8*G3*Q2h*Q2k*S - 8*G2*m2*Q2h*Sk - 8*G3*m2*Q2h*Sk + 48*G2*m2*Q2k*Sk +
			     48*G3*m2*Q2k*Sk - 8*G2*Q2h*Q2k*Sk - 8*G3*Q2h*Q2k*Sk + 64*G2*m2*S*Sk +
			     64*G3*m2*S*Sk - 32*G2*Q2k*S*Sk - 32*G3*Q2k*S*Sk + 8*G2*Q2h*Q2k*Sq2 +
			     8*G3*Q2h*Q2k*Sq2 + 12*G2*m2*pow(Q2k,2) + 12*G3*m2*pow(Q2k,2) +
			     4*G2*M2*pow(Q2k,2) + 4*G3*M2*pow(Q2k,2) - 4*G2*Q2h*pow(Q2k,2) -
			     4*G3*Q2h*pow(Q2k,2) - 16*G2*S*pow(Q2k,2) - 16*G3*S*pow(Q2k,2) +
			     64*G2*m2*pow(S,2) + 64*G3*m2*pow(S,2) - 32*G2*Q2k*pow(S,2) -
			     32*G3*Q2k*pow(S,2) + 48*G2*m2*pow(Sk,2) + 48*G3*m2*pow(Sk,2)) +
			  pow(l1k1,-2)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2)*
			   (-16*G2*m4*M2*Q2h - 16*G3*m4*M2*Q2h - 16*G2*m4*Q2h*Q2k -
			     16*G3*m4*Q2h*Q2k - 32*G2*m4*Q2h*S - 32*G3*m4*Q2h*S +
			     64*G2*m4*Q2k*S + 64*G3*m4*Q2k*S - 32*G2*m4*Q2h*Sk -
			     32*G3*m4*Q2h*Sk + 64*G2*m4*Q2k*Sk + 64*G3*m4*Q2k*Sk +
			     128*G2*m4*S*Sk + 128*G3*m4*S*Sk + 16*G2*m4*pow(Q2k,2) +
			     16*G3*m4*pow(Q2k,2) + 64*G2*m4*pow(S,2) + 64*G3*m4*pow(S,2) +
			     64*G2*m4*pow(Sk,2) + 64*G3*m4*pow(Sk,2)) +
			  pow(l1k1,-2)*pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (16*G2*m4*M2*Q2k + 16*G3*m4*M2*Q2k - 32*G2*m4*Q2k*S -
			     32*G3*m4*Q2k*S - 16*G2*m4*Q2k*Sk - 16*G3*m4*Q2k*Sk -
			     64*G2*m4*S*Sk - 64*G3*m4*S*Sk - 32*G2*m2*Q2k*S*Sk -
			     32*G3*m2*Q2k*S*Sk + 8*G2*m2*M2*pow(Q2k,2) +
			     8*G3*m2*M2*pow(Q2k,2) - 16*G2*m2*S*pow(Q2k,2) -
			     16*G3*m2*S*pow(Q2k,2) - 8*G2*m2*Sk*pow(Q2k,2) -
			     8*G3*m2*Sk*pow(Q2k,2) - 64*G2*m4*pow(S,2) - 64*G3*m4*pow(S,2) -
			     32*G2*m2*Q2k*pow(S,2) - 32*G3*m2*Q2k*pow(S,2) -
			     16*G2*m2*Q2k*pow(Sk,2) - 16*G3*m2*Q2k*pow(Sk,2)) +
			  pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (16*G2*m4*M2*Q2k + 16*G3*m4*M2*Q2k + 8*G2*m2*M2*Q2e*Q2k +
			     8*G3*m2*M2*Q2e*Q2k - 32*G2*m4*Q2k*S - 32*G3*m4*Q2k*S -
			     16*G2*m2*Q2e*Q2k*S - 16*G3*m2*Q2e*Q2k*S - 16*G2*m4*Q2k*Sk -
			     16*G3*m4*Q2k*Sk - 8*G2*m2*Q2e*Q2k*Sk - 8*G3*m2*Q2e*Q2k*Sk -
			     64*G2*m4*S*Sk - 64*G3*m4*S*Sk - 32*G2*m2*Q2e*S*Sk -
			     32*G3*m2*Q2e*S*Sk - 64*G2*m4*pow(S,2) - 64*G3*m4*pow(S,2) -
			     32*G2*m2*Q2e*pow(S,2) - 32*G3*m2*Q2e*pow(S,2) -
			     24*G2*m2*Q2e*pow(Sk,2) - 24*G3*m2*Q2e*pow(Sk,2) +
			     8*G2*m2*Q2k*pow(Sk,2) + 8*G3*m2*Q2k*pow(Sk,2)) +
			  pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (64*G2*m4 + 64*G3*m4 - 32*G2*k1k2*M2 - 32*G3*k1k2*M2 +
			     128*G2*m2*M2 + 128*G3*m2*M2 + 16*G2*M2*Q2e + 16*G3*M2*Q2e +
			     64*G2*m2*Q2h + 64*G3*m2*Q2h + 32*G2*M2*Q2h + 32*G3*M2*Q2h -
			     16*G2*m2*Q2k - 16*G3*m2*Q2k + 16*G2*Q2h*Q2k + 16*G3*Q2h*Q2k -
			     64*G2*m2*S - 64*G3*m2*S + 32*G2*Q2e*S + 32*G3*Q2e*S + 32*G2*Q2h*S +
			     32*G3*Q2h*S - 32*G2*Q2k*S - 32*G3*Q2k*S - 160*G2*m2*Sk -
			     160*G3*m2*Sk + 32*G2*Q2h*Sk + 32*G3*Q2h*Sk - 64*G2*S*Sk - 64*G3*S*Sk +
			     32*G2*Q2e*Sq2 + 32*G3*Q2e*Sq2 + 32*G2*Q2h*Sq2 + 32*G3*Q2h*Sq2 -
			     32*G2*Q2k*Sq2 - 32*G3*Q2k*Sq2 - 128*G2*S*Sq2 - 128*G3*S*Sq2 -
			     64*G2*Sk*Sq2 - 64*G3*Sk*Sq2 - 64*G2*pow(S,2) - 64*G3*pow(S,2) -
			     64*G2*pow(Sq2,2) - 64*G3*pow(Sq2,2)) +
			  pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   (32*G2*m4 + 32*G3*m4 + 32*G2*m2*M2 + 32*G3*m2*M2 +
			     8*G2*M2*Q2e + 8*G3*M2*Q2e + 32*G2*m2*Q2h + 32*G3*m2*Q2h +
			     8*G2*M2*Q2h + 8*G3*M2*Q2h - 8*G2*m2*Q2k - 8*G3*m2*Q2k +
			     8*G2*M2*Q2k + 8*G3*M2*Q2k - 32*G2*m2*S - 32*G3*m2*S +
			     16*G2*Q2e*S + 16*G3*Q2e*S + 16*G2*Q2k*S + 16*G3*Q2k*S - 80*G2*m2*Sk -
			     80*G3*m2*Sk + 32*G2*S*Sk + 32*G3*S*Sk + 16*G2*Q2e*Sq2 +
			     16*G3*Q2e*Sq2 + 16*G2*Q2k*Sq2 + 16*G3*Q2k*Sq2 - 128*G2*S*Sq2 -
			     128*G3*S*Sq2 + 32*G2*Sk*Sq2 + 32*G3*Sk*Sq2 - 96*G2*pow(S,2) -
			     96*G3*pow(S,2) - 32*G2*pow(Sq2,2) - 32*G3*pow(Sq2,2)) +
			  pow(l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (-32*G2*m4 - 32*G3*m4 - 16*G2*k1k2*M2 - 16*G3*k1k2*M2 -
			     64*G2*m2*M2 - 64*G3*m2*M2 - 16*G2*M2*Q2e - 16*G3*M2*Q2e -
			     32*G2*m2*Q2h - 32*G3*m2*Q2h - 8*G2*M2*Q2h - 8*G3*M2*Q2h +
			     8*G2*m2*Q2k + 8*G3*m2*Q2k + 32*G2*m2*S + 32*G3*m2*S - 16*G2*Q2e*S -
			     16*G3*Q2e*S + 16*G2*Q2h*S + 16*G3*Q2h*S + 16*G2*Q2k*S + 16*G3*Q2k*S +
			     80*G2*m2*Sk + 80*G3*m2*Sk + 32*G2*S*Sk + 32*G3*S*Sk - 16*G2*Q2e*Sq2 -
			     16*G3*Q2e*Sq2 + 16*G2*Q2k*Sq2 + 16*G3*Q2k*Sq2 + 64*G2*S*Sq2 +
			     64*G3*S*Sq2 + 32*G2*Sk*Sq2 + 32*G3*Sk*Sq2 + 32*G2*pow(S,2) +
			     32*G3*pow(S,2) + 32*G2*pow(Sq2,2) + 32*G3*pow(Sq2,2)) +
			  pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
			   (-16*G2*M2*Q2h - 16*G3*M2*Q2h - 32*G2*Q2h*S - 32*G3*Q2h*S -
			     32*G2*Q2h*Sq2 - 32*G3*Q2h*Sq2 + 128*G2*S*Sq2 + 128*G3*S*Sq2 +
			     64*G2*pow(S,2) + 64*G3*pow(S,2) + 64*G2*pow(Sq2,2) + 64*G3*pow(Sq2,2)) +
			  pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-2)*
			   (128*G2*m4*M2 + 128*G3*m4*M2 + 64*G2*m4*Q2h + 64*G3*m4*Q2h +
			     32*G2*m2*M2*Q2h + 32*G3*m2*M2*Q2h + 32*G2*m2*M2*Q2k +
			     32*G3*m2*M2*Q2k + 64*G2*m2*Q2h*S + 64*G3*m2*Q2h*S +
			     128*G2*m2*Q2k*S + 128*G3*m2*Q2k*S - 128*G2*m4*Sk - 128*G3*m4*Sk -
			     64*G2*m2*Q2k*Sk - 64*G3*m2*Q2k*Sk + 256*G2*m2*S*Sk +
			     256*G3*m2*S*Sk + 64*G2*m2*Q2h*Sq2 + 64*G3*m2*Q2h*Sq2 +
			     128*G2*m2*Q2k*Sq2 + 128*G3*m2*Q2k*Sq2 - 512*G2*m2*S*Sq2 -
			     512*G3*m2*S*Sq2 + 256*G2*m2*Sk*Sq2 + 256*G3*m2*Sk*Sq2 -
			     256*G2*m2*pow(S,2) - 256*G3*m2*pow(S,2) - 128*G2*m2*pow(Sk,2) -
			     128*G3*m2*pow(Sk,2) - 256*G2*m2*pow(Sq2,2) - 256*G3*m2*pow(Sq2,2)) +
			  pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (64*G2*m4*M2 + 64*G3*m4*M2 + 32*G2*m4*Q2h + 32*G3*m4*Q2h +
			     16*G2*m2*M2*Q2h + 16*G3*m2*M2*Q2h - 8*G2*M2*Q2e*Q2h -
			     8*G3*M2*Q2e*Q2h + 32*G2*m4*Q2k + 32*G3*m4*Q2k +
			     48*G2*m2*M2*Q2k + 48*G3*m2*M2*Q2k - 24*G2*m2*Q2h*Q2k -
			     24*G3*m2*Q2h*Q2k + 16*G2*M2*Q2h*Q2k + 16*G3*M2*Q2h*Q2k -
			     8*G2*Q2e*Q2h*Q2k - 8*G3*Q2e*Q2h*Q2k - 32*G2*m2*Q2h*S -
			     32*G3*m2*Q2h*S - 32*G2*Q2e*Q2h*S - 32*G3*Q2e*Q2h*S + 64*G2*m2*Q2k*S +
			     64*G3*m2*Q2k*S + 16*G2*Q2e*Q2k*S + 16*G3*Q2e*Q2k*S - 16*G2*Q2h*Q2k*S -
			     16*G3*Q2h*Q2k*S - 80*G2*m2*Q2h*Sk - 80*G3*m2*Q2h*Sk -
			     16*G2*Q2e*Q2h*Sk - 16*G3*Q2e*Q2h*Sk + 32*G2*m2*Q2k*Sk +
			     32*G3*m2*Q2k*Sk + 16*G2*Q2e*Q2k*Sk + 16*G3*Q2e*Q2k*Sk -
			     16*G2*Q2h*Q2k*Sk - 16*G3*Q2h*Q2k*Sk + 128*G2*m2*S*Sk +
			     128*G3*m2*S*Sk + 96*G2*Q2e*S*Sk + 96*G3*Q2e*S*Sk - 96*G2*Q2h*S*Sk -
			     96*G3*Q2h*S*Sk - 128*G2*Q2k*S*Sk - 128*G3*Q2k*S*Sk - 64*G2*m4*Sq2 -
			     64*G3*m4*Sq2 + 16*G2*m2*Q2h*Sq2 + 16*G3*m2*Q2h*Sq2 +
			     64*G2*m2*Q2k*Sq2 + 64*G3*m2*Q2k*Sq2 + 16*G2*Q2h*Q2k*Sq2 +
			     16*G3*Q2h*Q2k*Sq2 - 128*G2*m2*S*Sq2 - 128*G3*m2*S*Sq2 +
			     32*G2*Q2e*S*Sq2 + 32*G3*Q2e*S*Sq2 + 32*G2*Q2h*S*Sq2 + 32*G3*Q2h*S*Sq2 -
			     64*G2*Q2k*S*Sq2 - 64*G3*Q2k*S*Sq2 + 64*G2*m2*Sk*Sq2 +
			     64*G3*m2*Sk*Sq2 - 32*G2*Q2k*Sk*Sq2 - 32*G3*Q2k*Sk*Sq2 +
			     16*G2*m2*pow(Q2h,2) + 16*G3*m2*pow(Q2h,2) + 8*G2*M2*pow(Q2h,2) +
			     8*G3*M2*pow(Q2h,2) + 8*G2*Q2k*pow(Q2h,2) + 8*G3*Q2k*pow(Q2h,2) +
			     16*G2*S*pow(Q2h,2) + 16*G3*S*pow(Q2h,2) + 16*G2*Sk*pow(Q2h,2) +
			     16*G3*Sk*pow(Q2h,2) + 16*G2*m2*pow(Q2k,2) + 16*G3*m2*pow(Q2k,2) +
			     16*G2*M2*pow(Q2k,2) + 16*G3*M2*pow(Q2k,2) - 32*G2*S*pow(Q2k,2) -
			     32*G3*S*pow(Q2k,2) - 16*G2*Sk*pow(Q2k,2) - 16*G3*Sk*pow(Q2k,2) -
			     128*G2*m2*pow(S,2) - 128*G3*m2*pow(S,2) + 96*G2*Q2e*pow(S,2) +
			     96*G3*Q2e*pow(S,2) - 64*G2*Q2h*pow(S,2) - 64*G3*Q2h*pow(S,2) -
			     160*G2*Q2k*pow(S,2) - 160*G3*Q2k*pow(S,2) + 64*G2*m2*pow(Sk,2) +
			     64*G3*m2*pow(Sk,2) + 32*G2*Q2e*pow(Sk,2) + 32*G3*Q2e*pow(Sk,2) -
			     32*G2*Q2h*pow(Sk,2) - 32*G3*Q2h*pow(Sk,2) - 32*G2*Q2k*pow(Sk,2) -
			     32*G3*Q2k*pow(Sk,2) - 64*G2*m2*pow(Sq2,2) - 64*G3*m2*pow(Sq2,2)) +
			  pow(l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-2)*
			   (-64*G2*m4*M2 - 64*G3*m4*M2 - 32*G2*m4*Q2h - 32*G3*m4*Q2h -
			     16*G2*m2*M2*Q2h - 16*G3*m2*M2*Q2h - 16*G2*m2*M2*Q2k -
			     16*G3*m2*M2*Q2k - 16*G2*m2*Q2h*Q2k - 16*G3*m2*Q2h*Q2k +
			     32*G2*m2*Q2h*S + 32*G3*m2*Q2h*S - 64*G2*m2*Q2k*S - 64*G3*m2*Q2k*S +
			     64*G2*m4*Sk + 64*G3*m4*Sk - 32*G2*m2*Q2h*Sk - 32*G3*m2*Q2h*Sk +
			     32*G2*m2*Q2k*Sk + 32*G3*m2*Q2k*Sk - 128*G2*m2*S*Sk -
			     128*G3*m2*S*Sk + 32*G2*m2*Q2h*Sq2 + 32*G3*m2*Q2h*Sq2 -
			     64*G2*m2*Q2k*Sq2 - 64*G3*m2*Q2k*Sq2 + 256*G2*m2*S*Sq2 +
			     256*G3*m2*S*Sq2 - 128*G2*m2*Sk*Sq2 - 128*G3*m2*Sk*Sq2 +
			     128*G2*m2*pow(S,2) + 128*G3*m2*pow(S,2) + 64*G2*m2*pow(Sk,2) +
			     64*G3*m2*pow(Sk,2) + 128*G2*m2*pow(Sq2,2) + 128*G3*m2*pow(Sq2,2)) +
			  pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (-32*G2*k1k2*M2*Q2h - 32*G3*k1k2*M2*Q2h - 64*G2*m2*M2*Q2h -
			     64*G3*m2*M2*Q2h - 64*G2*k1k2*Q2h*S - 64*G3*k1k2*Q2h*S -
			     128*G2*m2*Q2h*S - 128*G3*m2*Q2h*S - 64*G2*k1k2*Q2h*Sq2 -
			     64*G3*k1k2*Q2h*Sq2 - 128*G2*m2*Q2h*Sq2 - 128*G3*m2*Q2h*Sq2 +
			     256*G2*k1k2*S*Sq2 + 256*G3*k1k2*S*Sq2 + 512*G2*m2*S*Sq2 +
			     512*G3*m2*S*Sq2 + 128*G2*k1k2*pow(S,2) + 128*G3*k1k2*pow(S,2) +
			     256*G2*m2*pow(S,2) + 256*G3*m2*pow(S,2) + 128*G2*k1k2*pow(Sq2,2) +
			     128*G3*k1k2*pow(Sq2,2) + 256*G2*m2*pow(Sq2,2) + 256*G3*m2*pow(Sq2,2)) \
			+ pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (128*G2*m4*M2*Q2h + 128*G3*m4*M2*Q2h + 256*G2*m4*Q2h*S +
			     256*G3*m4*Q2h*S + 256*G2*m4*Q2h*Sq2 + 256*G3*m4*Q2h*Sq2 -
			     1024*G2*m4*S*Sq2 - 1024*G3*m4*S*Sq2 - 512*G2*m4*pow(S,2) -
			     512*G3*m4*pow(S,2) - 512*G2*m4*pow(Sq2,2) - 512*G3*m4*pow(Sq2,2)) +
			  pow(l1k2,-2)*pow(2*l1k1 + Q2e - Q2k,-2)*
			   (-16*G2*m4*M2*Q2h - 16*G3*m4*M2*Q2h - 16*G2*m4*Q2h*Q2k -
			     16*G3*m4*Q2h*Q2k + 32*G2*m4*Q2h*S + 32*G3*m4*Q2h*S -
			     64*G2*m4*Q2k*S - 64*G3*m4*Q2k*S - 32*G2*m4*Q2h*Sk -
			     32*G3*m4*Q2h*Sk + 64*G2*m4*Q2k*Sk + 64*G3*m4*Q2k*Sk -
			     128*G2*m4*S*Sk - 128*G3*m4*S*Sk + 32*G2*m4*Q2h*Sq2 +
			     32*G3*m4*Q2h*Sq2 - 64*G2*m4*Q2k*Sq2 - 64*G3*m4*Q2k*Sq2 +
			     128*G2*m4*S*Sq2 + 128*G3*m4*S*Sq2 - 128*G2*m4*Sk*Sq2 -
			     128*G3*m4*Sk*Sq2 + 16*G2*m4*pow(Q2k,2) + 16*G3*m4*pow(Q2k,2) +
			     64*G2*m4*pow(S,2) + 64*G3*m4*pow(S,2) + 64*G2*m4*pow(Sk,2) +
			     64*G3*m4*pow(Sk,2) + 64*G2*m4*pow(Sq2,2) + 64*G3*m4*pow(Sq2,2)) +
			  pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*pow(2*l1k1 + Q2e - Q2k,-2)*
			   (-64*G2*m4*M2*Q2h - 64*G3*m4*M2*Q2h - 128*G2*m4*Q2h*S -
			     128*G3*m4*Q2h*S - 128*G2*m4*Q2h*Sq2 - 128*G3*m4*Q2h*Sq2 +
			     512*G2*m4*S*Sq2 + 512*G3*m4*S*Sq2 + 256*G2*m4*pow(S,2) +
			     256*G3*m4*pow(S,2) + 256*G2*m4*pow(Sq2,2) + 256*G3*m4*pow(Sq2,2)) +
			  pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2)*
			   (-64*G2*m4*M2*Q2h - 64*G3*m4*M2*Q2h - 128*G2*m4*Q2h*S -
			     128*G3*m4*Q2h*S - 128*G2*m4*Q2h*Sq2 - 128*G3*m4*Q2h*Sq2 +
			     512*G2*m4*S*Sq2 + 512*G3*m4*S*Sq2 + 256*G2*m4*pow(S,2) +
			     256*G3*m4*pow(S,2) + 256*G2*m4*pow(Sq2,2) + 256*G3*m4*pow(Sq2,2)) +
			  pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-64*G2*m4*M2*Q2k - 64*G3*m4*M2*Q2k - 32*G2*m2*M2*Q2e*Q2k -
			     32*G3*m2*M2*Q2e*Q2k - 128*G2*m4*Q2k*S - 128*G3*m4*Q2k*S -
			     64*G2*m2*Q2e*Q2k*S - 64*G3*m2*Q2e*Q2k*S + 64*G2*m4*Q2k*Sk +
			     64*G3*m4*Q2k*Sk + 32*G2*m2*Q2e*Q2k*Sk + 32*G3*m2*Q2e*Q2k*Sk -
			     256*G2*m4*S*Sk - 256*G3*m4*S*Sk - 128*G2*m2*Q2e*S*Sk -
			     128*G3*m2*Q2e*S*Sk - 128*G2*m4*Q2k*Sq2 - 128*G3*m4*Q2k*Sq2 -
			     64*G2*m2*Q2e*Q2k*Sq2 - 64*G3*m2*Q2e*Q2k*Sq2 + 512*G2*m4*S*Sq2 +
			     512*G3*m4*S*Sq2 + 256*G2*m2*Q2e*S*Sq2 + 256*G3*m2*Q2e*S*Sq2 -
			     256*G2*m4*Sk*Sq2 - 256*G3*m4*Sk*Sq2 - 128*G2*m2*Q2e*Sk*Sq2 -
			     128*G3*m2*Q2e*Sk*Sq2 + 256*G2*m4*pow(S,2) + 256*G3*m4*pow(S,2) +
			     128*G2*m2*Q2e*pow(S,2) + 128*G3*m2*Q2e*pow(S,2) +
			     96*G2*m2*Q2e*pow(Sk,2) + 96*G3*m2*Q2e*pow(Sk,2) -
			     32*G2*m2*Q2k*pow(Sk,2) - 32*G3*m2*Q2k*pow(Sk,2) +
			     256*G2*m4*pow(Sq2,2) + 256*G3*m4*pow(Sq2,2) +
			     128*G2*m2*Q2e*pow(Sq2,2) + 128*G3*m2*Q2e*pow(Sq2,2)) +
			  pow(l1k1,-1)*pow(l1k2,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-32*G2*m4*M2 - 32*G3*m4*M2 - 8*G2*m2*M2*Q2e -
			     8*G3*m2*M2*Q2e - 16*G2*m4*Q2h - 16*G3*m4*Q2h +
			     32*G2*m2*M2*Q2h + 32*G3*m2*M2*Q2h + 4*G2*m2*Q2e*Q2h +
			     4*G3*m2*Q2e*Q2h + 8*G2*M2*Q2e*Q2h + 8*G3*M2*Q2e*Q2h -
			     8*G2*m2*M2*Q2k - 8*G3*m2*M2*Q2k - 12*G2*m2*Q2e*Q2k -
			     12*G3*m2*Q2e*Q2k - 4*G2*M2*Q2e*Q2k - 4*G3*M2*Q2e*Q2k +
			     8*G2*m2*Q2h*Q2k + 8*G3*m2*Q2h*Q2k + 8*G2*M2*Q2h*Q2k +
			     8*G3*M2*Q2h*Q2k + 8*G2*Q2e*Q2h*Q2k + 8*G3*Q2e*Q2h*Q2k -
			     16*G2*m2*Q2h*S - 16*G3*m2*Q2h*S + 24*G2*Q2e*Q2h*S + 24*G3*Q2e*Q2h*S +
			     32*G2*m2*Q2k*S + 32*G3*m2*Q2k*S - 16*G2*Q2e*Q2k*S - 16*G3*Q2e*Q2k*S +
			     32*G2*m4*Sk + 32*G3*m4*Sk - 32*G2*m2*Q2e*Sk - 32*G3*m2*Q2e*Sk +
			     24*G2*m2*Q2h*Sk + 24*G3*m2*Q2h*Sk + 16*G2*Q2e*Q2h*Sk +
			     16*G3*Q2e*Q2h*Sk - 16*G2*m2*Q2k*Sk - 16*G3*m2*Q2k*Sk -
			     32*G2*Q2e*Q2k*Sk - 32*G3*Q2e*Q2k*Sk + 16*G2*Q2h*Q2k*Sk +
			     16*G3*Q2h*Q2k*Sk + 64*G2*m2*S*Sk + 64*G3*m2*S*Sk - 32*G2*Q2e*S*Sk -
			     32*G3*Q2e*S*Sk - 16*G2*m2*Q2h*Sq2 - 16*G3*m2*Q2h*Sq2 +
			     8*G2*Q2e*Q2h*Sq2 + 8*G3*Q2e*Q2h*Sq2 + 16*G2*m2*Q2k*Sq2 +
			     16*G3*m2*Q2k*Sq2 - 8*G2*Q2h*Q2k*Sq2 - 8*G3*Q2h*Q2k*Sq2 -
			     64*G2*m2*S*Sq2 - 64*G3*m2*S*Sq2 - 64*G2*Q2h*S*Sq2 - 64*G3*Q2h*S*Sq2 +
			     32*G2*Q2k*S*Sq2 + 32*G3*Q2k*S*Sq2 + 32*G2*m2*Sk*Sq2 +
			     32*G3*m2*Sk*Sq2 + 32*G2*Q2k*Sk*Sq2 + 32*G3*Q2k*Sk*Sq2 +
			     4*G2*M2*pow(Q2e,2) + 4*G3*M2*pow(Q2e,2) + 4*G2*m2*pow(Q2h,2) +
			     4*G3*m2*pow(Q2h,2) - 4*G2*m2*pow(Q2k,2) - 4*G3*m2*pow(Q2k,2) -
			     8*G2*Q2e*pow(Q2k,2) - 8*G3*Q2e*pow(Q2k,2) + 8*G2*Q2h*pow(Q2k,2) +
			     8*G3*Q2h*pow(Q2k,2) - 32*G2*Sk*pow(Q2k,2) - 32*G3*Sk*pow(Q2k,2) +
			     16*G2*Sq2*pow(Q2k,2) + 16*G3*Sq2*pow(Q2k,2) - 8*G2*pow(Q2k,3) -
			     8*G3*pow(Q2k,3) - 64*G2*m2*pow(S,2) - 64*G3*m2*pow(S,2) -
			     64*G2*Q2h*pow(S,2) - 64*G3*Q2h*pow(S,2) + 32*G2*Q2k*pow(S,2) +
			     32*G3*Q2k*pow(S,2) - 32*G2*m2*pow(Sk,2) - 32*G3*m2*pow(Sk,2) -
			     32*G2*Q2e*pow(Sk,2) - 32*G3*Q2e*pow(Sk,2) - 32*G2*Q2k*pow(Sk,2) -
			     32*G3*Q2k*pow(Sk,2) - 16*G2*Q2h*pow(Sq2,2) - 16*G3*Q2h*pow(Sq2,2)) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-32*G2*m4*M2 - 32*G3*m4*M2 - 24*G2*m2*M2*Q2e -
			     24*G3*m2*M2*Q2e - 16*G2*m4*Q2h - 16*G3*m4*Q2h -
			     16*G2*m2*M2*Q2h - 16*G3*m2*M2*Q2h + 4*G2*m2*Q2e*Q2h +
			     4*G3*m2*Q2e*Q2h - 8*G2*m2*M2*Q2k - 8*G3*m2*M2*Q2k -
			     12*G2*m2*Q2e*Q2k - 12*G3*m2*Q2e*Q2k - 4*G2*M2*Q2e*Q2k -
			     4*G3*M2*Q2e*Q2k + 16*G2*m2*Q2h*Q2k + 16*G3*m2*Q2h*Q2k +
			     8*G2*M2*Q2h*Q2k + 8*G3*M2*Q2h*Q2k - 4*G2*Q2e*Q2h*Q2k -
			     4*G3*Q2e*Q2h*Q2k - 48*G2*m2*Q2h*S - 48*G3*m2*Q2h*S + 8*G2*Q2e*Q2h*S +
			     8*G3*Q2e*Q2h*S + 32*G2*m2*Q2k*S + 32*G3*m2*Q2k*S - 16*G2*Q2e*Q2k*S -
			     16*G3*Q2e*Q2k*S + 32*G2*m4*Sk + 32*G3*m4*Sk - 32*G2*m2*Q2e*Sk -
			     32*G3*m2*Q2e*Sk + 40*G2*m2*Q2h*Sk + 40*G3*m2*Q2h*Sk -
			     8*G2*Q2e*Q2h*Sk - 8*G3*Q2e*Q2h*Sk - 48*G2*m2*Q2k*Sk -
			     48*G3*m2*Q2k*Sk + 64*G2*m2*S*Sk + 64*G3*m2*S*Sk - 32*G2*Q2e*S*Sk -
			     32*G3*Q2e*S*Sk - 48*G2*m2*Q2h*Sq2 - 48*G3*m2*Q2h*Sq2 -
			     8*G2*Q2e*Q2h*Sq2 - 8*G3*Q2e*Q2h*Sq2 + 48*G2*m2*Q2k*Sq2 +
			     48*G3*m2*Q2k*Sq2 + 8*G2*Q2h*Q2k*Sq2 + 8*G3*Q2h*Q2k*Sq2 +
			     64*G2*m2*S*Sq2 + 64*G3*m2*S*Sq2 + 32*G2*Q2k*S*Sq2 + 32*G3*Q2k*S*Sq2 +
			     96*G2*m2*Sk*Sq2 + 96*G3*m2*Sk*Sq2 - 4*G2*M2*pow(Q2e,2) -
			     4*G3*M2*pow(Q2e,2) - 12*G2*m2*pow(Q2h,2) - 12*G3*m2*pow(Q2h,2) -
			     16*G2*M2*pow(Q2h,2) - 16*G3*M2*pow(Q2h,2) + 4*G2*Q2e*pow(Q2h,2) +
			     4*G3*Q2e*pow(Q2h,2) - 16*G2*Sq2*pow(Q2h,2) - 16*G3*Sq2*pow(Q2h,2) -
			     12*G2*m2*pow(Q2k,2) - 12*G3*m2*pow(Q2k,2) + 64*G2*m2*pow(S,2) +
			     64*G3*m2*pow(S,2) + 32*G2*Q2k*pow(S,2) + 32*G3*Q2k*pow(S,2) -
			     64*G2*m2*pow(Sk,2) - 64*G3*m2*pow(Sk,2) + 16*G2*Q2h*pow(Sq2,2) +
			     16*G3*Q2h*pow(Sq2,2)) + pow(l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-64*G2*m4*M2 - 64*G3*m4*M2 - 48*G2*m2*M2*Q2e -
			     48*G3*m2*M2*Q2e - 32*G2*m4*Q2h - 32*G3*m4*Q2h -
			     32*G2*m2*M2*Q2h - 32*G3*m2*M2*Q2h + 8*G2*m2*Q2e*Q2h +
			     8*G3*m2*Q2e*Q2h - 16*G2*m2*M2*Q2k - 16*G3*m2*M2*Q2k -
			     24*G2*m2*Q2e*Q2k - 24*G3*m2*Q2e*Q2k - 8*G2*M2*Q2e*Q2k -
			     8*G3*M2*Q2e*Q2k + 32*G2*m2*Q2h*Q2k + 32*G3*m2*Q2h*Q2k +
			     16*G2*M2*Q2h*Q2k + 16*G3*M2*Q2h*Q2k - 8*G2*Q2e*Q2h*Q2k -
			     8*G3*Q2e*Q2h*Q2k + 96*G2*m2*Q2h*S + 96*G3*m2*Q2h*S -
			     16*G2*Q2e*Q2h*S - 16*G3*Q2e*Q2h*S - 64*G2*m2*Q2k*S - 64*G3*m2*Q2k*S +
			     32*G2*Q2e*Q2k*S + 32*G3*Q2e*Q2k*S + 64*G2*m4*Sk + 64*G3*m4*Sk -
			     64*G2*m2*Q2e*Sk - 64*G3*m2*Q2e*Sk + 80*G2*m2*Q2h*Sk +
			     80*G3*m2*Q2h*Sk - 16*G2*Q2e*Q2h*Sk - 16*G3*Q2e*Q2h*Sk -
			     96*G2*m2*Q2k*Sk - 96*G3*m2*Q2k*Sk - 128*G2*m2*S*Sk -
			     128*G3*m2*S*Sk + 64*G2*Q2e*S*Sk + 64*G3*Q2e*S*Sk - 32*G2*Q2e*Q2h*Sq2 -
			     32*G3*Q2e*Q2h*Sq2 + 32*G2*m2*Q2k*Sq2 + 32*G3*m2*Q2k*Sq2 +
			     32*G2*Q2e*Q2k*Sq2 + 32*G3*Q2e*Q2k*Sq2 + 16*G2*Q2h*Q2k*Sq2 +
			     16*G3*Q2h*Q2k*Sq2 + 128*G2*m2*S*Sq2 + 128*G3*m2*S*Sq2 +
			     64*G2*Q2k*S*Sq2 + 64*G3*Q2k*S*Sq2 + 64*G2*m2*Sk*Sq2 +
			     64*G3*m2*Sk*Sq2 + 64*G2*Q2e*Sk*Sq2 + 64*G3*Q2e*Sk*Sq2 -
			     8*G2*M2*pow(Q2e,2) - 8*G3*M2*pow(Q2e,2) - 24*G2*m2*pow(Q2h,2) -
			     24*G3*m2*pow(Q2h,2) - 32*G2*M2*pow(Q2h,2) - 32*G3*M2*pow(Q2h,2) +
			     8*G2*Q2e*pow(Q2h,2) + 8*G3*Q2e*pow(Q2h,2) - 32*G2*Sq2*pow(Q2h,2) -
			     32*G3*Sq2*pow(Q2h,2) - 24*G2*m2*pow(Q2k,2) - 24*G3*m2*pow(Q2k,2) +
			     128*G2*m2*pow(S,2) + 128*G3*m2*pow(S,2) + 64*G2*Q2k*pow(S,2) +
			     64*G3*Q2k*pow(S,2) - 128*G2*m2*pow(Sk,2) - 128*G3*m2*pow(Sk,2) +
			     32*G2*Q2h*pow(Sq2,2) + 32*G3*Q2h*pow(Sq2,2)) +
			  pow(l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (64*G2*m4*M2 + 64*G3*m4*M2 + 16*G2*m2*M2*Q2e +
			     16*G3*m2*M2*Q2e + 32*G2*m4*Q2h + 32*G3*m4*Q2h -
			     64*G2*m2*M2*Q2h - 64*G3*m2*M2*Q2h - 8*G2*m2*Q2e*Q2h -
			     8*G3*m2*Q2e*Q2h - 16*G2*M2*Q2e*Q2h - 16*G3*M2*Q2e*Q2h +
			     16*G2*m2*M2*Q2k + 16*G3*m2*M2*Q2k + 24*G2*m2*Q2e*Q2k +
			     24*G3*m2*Q2e*Q2k + 8*G2*M2*Q2e*Q2k + 8*G3*M2*Q2e*Q2k -
			     16*G2*m2*Q2h*Q2k - 16*G3*m2*Q2h*Q2k - 16*G2*M2*Q2h*Q2k -
			     16*G3*M2*Q2h*Q2k - 16*G2*Q2e*Q2h*Q2k - 16*G3*Q2e*Q2h*Q2k -
			     32*G2*m2*Q2h*S - 32*G3*m2*Q2h*S + 48*G2*Q2e*Q2h*S + 48*G3*Q2e*Q2h*S +
			     64*G2*m2*Q2k*S + 64*G3*m2*Q2k*S - 32*G2*Q2e*Q2k*S - 32*G3*Q2e*Q2k*S -
			     64*G2*m4*Sk - 64*G3*m4*Sk + 64*G2*m2*Q2e*Sk + 64*G3*m2*Q2e*Sk -
			     48*G2*m2*Q2h*Sk - 48*G3*m2*Q2h*Sk - 32*G2*Q2e*Q2h*Sk -
			     32*G3*Q2e*Q2h*Sk + 32*G2*m2*Q2k*Sk + 32*G3*m2*Q2k*Sk +
			     64*G2*Q2e*Q2k*Sk + 64*G3*Q2e*Q2k*Sk - 32*G2*Q2h*Q2k*Sk -
			     32*G3*Q2h*Q2k*Sk + 128*G2*m2*S*Sk + 128*G3*m2*S*Sk - 64*G2*Q2e*S*Sk -
			     64*G3*Q2e*S*Sk + 32*G2*Q2e*Q2h*Sq2 + 32*G3*Q2e*Q2h*Sq2 +
			     32*G2*m2*Q2k*Sq2 + 32*G3*m2*Q2k*Sq2 - 32*G2*Q2e*Q2k*Sq2 -
			     32*G3*Q2e*Q2k*Sq2 + 16*G2*Q2h*Q2k*Sq2 + 16*G3*Q2h*Q2k*Sq2 +
			     128*G2*m2*S*Sq2 + 128*G3*m2*S*Sq2 + 128*G2*Q2h*S*Sq2 +
			     128*G3*Q2h*S*Sq2 - 64*G2*Q2k*S*Sq2 - 64*G3*Q2k*S*Sq2 +
			     64*G2*m2*Sk*Sq2 + 64*G3*m2*Sk*Sq2 - 64*G2*Q2e*Sk*Sq2 -
			     64*G3*Q2e*Sk*Sq2 - 64*G2*Q2k*Sk*Sq2 - 64*G3*Q2k*Sk*Sq2 -
			     8*G2*M2*pow(Q2e,2) - 8*G3*M2*pow(Q2e,2) - 8*G2*m2*pow(Q2h,2) -
			     8*G3*m2*pow(Q2h,2) + 8*G2*m2*pow(Q2k,2) + 8*G3*m2*pow(Q2k,2) +
			     16*G2*Q2e*pow(Q2k,2) + 16*G3*Q2e*pow(Q2k,2) - 16*G2*Q2h*pow(Q2k,2) -
			     16*G3*Q2h*pow(Q2k,2) + 64*G2*Sk*pow(Q2k,2) + 64*G3*Sk*pow(Q2k,2) -
			     32*G2*Sq2*pow(Q2k,2) - 32*G3*Sq2*pow(Q2k,2) + 16*G2*pow(Q2k,3) +
			     16*G3*pow(Q2k,3) + 128*G2*m2*pow(S,2) + 128*G3*m2*pow(S,2) +
			     128*G2*Q2h*pow(S,2) + 128*G3*Q2h*pow(S,2) - 64*G2*Q2k*pow(S,2) -
			     64*G3*Q2k*pow(S,2) + 64*G2*m2*pow(Sk,2) + 64*G3*m2*pow(Sk,2) +
			     64*G2*Q2e*pow(Sk,2) + 64*G3*Q2e*pow(Sk,2) + 64*G2*Q2k*pow(Sk,2) +
			     64*G3*Q2k*pow(Sk,2) + 32*G2*Q2h*pow(Sq2,2) + 32*G3*Q2h*pow(Sq2,2)) +
			  pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2)*
			   (64*G2*m4*M2*Q2h + 64*G3*m4*M2*Q2h + 32*G2*m4*Q2h*Q2k +
			     32*G3*m4*Q2h*Q2k + 128*G2*m4*Q2h*S + 128*G3*m4*Q2h*S -
			     128*G2*m4*Q2k*S - 128*G3*m4*Q2k*S - 64*G2*m2*Q2h*Q2k*S -
			     64*G3*m2*Q2h*Q2k*S + 64*G2*m4*Q2h*Sk + 64*G3*m4*Q2h*Sk -
			     64*G2*m2*Q2h*Q2k*Sk - 64*G3*m2*Q2h*Q2k*Sk - 256*G2*m4*S*Sk -
			     256*G3*m4*S*Sk - 128*G2*m2*Q2h*S*Sk - 128*G3*m2*Q2h*S*Sk +
			     64*G2*m4*Q2h*Sq2 + 64*G3*m4*Q2h*Sq2 - 128*G2*m4*Q2k*Sq2 -
			     128*G3*m4*Q2k*Sq2 - 256*G2*m4*S*Sq2 - 256*G3*m4*S*Sq2 -
			     128*G2*m2*Q2h*S*Sq2 - 128*G3*m2*Q2h*S*Sq2 - 256*G2*m4*Sk*Sq2 -
			     256*G3*m4*Sk*Sq2 + 32*G2*m2*M2*pow(Q2h,2) +
			     32*G3*m2*M2*pow(Q2h,2) + 16*G2*m2*Q2k*pow(Q2h,2) +
			     16*G3*m2*Q2k*pow(Q2h,2) + 64*G2*m2*S*pow(Q2h,2) +
			     64*G3*m2*S*pow(Q2h,2) + 32*G2*m2*Sk*pow(Q2h,2) +
			     32*G3*m2*Sk*pow(Q2h,2) + 32*G2*m2*Sq2*pow(Q2h,2) +
			     32*G3*m2*Sq2*pow(Q2h,2) - 16*G2*m2*Q2h*pow(Q2k,2) -
			     16*G3*m2*Q2h*pow(Q2k,2) - 256*G2*m4*pow(S,2) - 256*G3*m4*pow(S,2) -
			     128*G2*m2*Q2h*pow(S,2) - 128*G3*m2*Q2h*pow(S,2) -
			     64*G2*m2*Q2h*pow(Sk,2) - 64*G3*m2*Q2h*pow(Sk,2) -
			     64*G2*m2*Q2h*pow(Sq2,2) - 64*G3*m2*Q2h*pow(Sq2,2)) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-2)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (-16*G2*m4*M2*Q2h - 16*G3*m4*M2*Q2h - 8*G2*m4*Q2h*Q2k -
			     8*G3*m4*Q2h*Q2k + 32*G2*m4*Q2h*S + 32*G3*m4*Q2h*S -
			     32*G2*m4*Q2k*S - 32*G3*m4*Q2k*S - 16*G2*m2*Q2h*Q2k*S -
			     16*G3*m2*Q2h*Q2k*S - 16*G2*m4*Q2h*Sk - 16*G3*m4*Q2h*Sk +
			     16*G2*m2*Q2h*Q2k*Sk + 16*G3*m2*Q2h*Q2k*Sk - 64*G2*m4*S*Sk -
			     64*G3*m4*S*Sk - 32*G2*m2*Q2h*S*Sk - 32*G3*m2*Q2h*S*Sk +
			     16*G2*m4*Q2h*Sq2 + 16*G3*m4*Q2h*Sq2 - 16*G2*m2*Q2h*Q2k*Sq2 -
			     16*G3*m2*Q2h*Q2k*Sq2 + 64*G2*m4*S*Sq2 + 64*G3*m4*S*Sq2 +
			     32*G2*m2*Q2h*S*Sq2 + 32*G3*m2*Q2h*S*Sq2 - 32*G2*m2*Q2h*Sk*Sq2 -
			     32*G3*m2*Q2h*Sk*Sq2 - 8*G2*m2*M2*pow(Q2h,2) -
			     8*G3*m2*M2*pow(Q2h,2) - 4*G2*m2*Q2k*pow(Q2h,2) -
			     4*G3*m2*Q2k*pow(Q2h,2) + 16*G2*m2*S*pow(Q2h,2) +
			     16*G3*m2*S*pow(Q2h,2) - 8*G2*m2*Sk*pow(Q2h,2) -
			     8*G3*m2*Sk*pow(Q2h,2) + 8*G2*m2*Sq2*pow(Q2h,2) +
			     8*G3*m2*Sq2*pow(Q2h,2) + 4*G2*m2*Q2h*pow(Q2k,2) +
			     4*G3*m2*Q2h*pow(Q2k,2) + 64*G2*m4*pow(S,2) + 64*G3*m4*pow(S,2) +
			     32*G2*m2*Q2h*pow(S,2) + 32*G3*m2*Q2h*pow(S,2) +
			     16*G2*m2*Q2h*pow(Sk,2) + 16*G3*m2*Q2h*pow(Sk,2) +
			     16*G2*m2*Q2h*pow(Sq2,2) + 16*G3*m2*Q2h*pow(Sq2,2)) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-1)*(64*G2*m4*M2 + 64*G3*m4*M2 +
			     32*G2*m4*Q2h + 32*G3*m4*Q2h + 16*G2*m2*M2*Q2h +
			     16*G3*m2*M2*Q2h + 32*G2*m4*Q2k + 32*G3*m4*Q2k +
			     80*G2*m2*M2*Q2k + 80*G3*m2*M2*Q2k + 8*G2*M2*Q2e*Q2k +
			     8*G3*M2*Q2e*Q2k + 16*G2*m2*Q2h*Q2k + 16*G3*m2*Q2h*Q2k -
			     8*G2*M2*Q2h*Q2k - 8*G3*M2*Q2h*Q2k - 32*G2*Q2e*Q2k*S -
			     32*G3*Q2e*Q2k*S + 32*G2*Q2h*Q2k*S + 32*G3*Q2h*Q2k*S +
			     32*G2*m2*Q2k*Sk + 32*G3*m2*Q2k*Sk + 16*G2*Q2e*Q2k*Sk +
			     16*G3*Q2e*Q2k*Sk - 16*G2*Q2h*Q2k*Sk - 16*G3*Q2h*Q2k*Sk -
			     128*G2*Q2e*S*Sk - 128*G3*Q2e*S*Sk + 128*G2*Q2h*S*Sk + 128*G3*Q2h*S*Sk +
			     128*G2*Q2k*S*Sk + 128*G3*Q2k*S*Sk - 64*G2*m4*Sq2 - 64*G3*m4*Sq2 +
			     48*G2*m2*Q2h*Sq2 + 48*G3*m2*Q2h*Sq2 - 16*G2*Q2e*Q2k*Sq2 -
			     16*G3*Q2e*Q2k*Sq2 + 16*G2*Q2h*Q2k*Sq2 + 16*G3*Q2h*Q2k*Sq2 -
			     128*G2*m2*S*Sq2 - 128*G3*m2*S*Sq2 + 160*G2*Q2e*S*Sq2 +
			     160*G3*Q2e*S*Sq2 - 160*G2*Q2h*S*Sq2 - 160*G3*Q2h*S*Sq2 -
			     320*G2*Q2k*S*Sq2 - 320*G3*Q2k*S*Sq2 - 64*G2*m2*Sk*Sq2 -
			     64*G3*m2*Sk*Sq2 - 96*G2*Q2e*Sk*Sq2 - 96*G3*Q2e*Sk*Sq2 +
			     96*G2*Q2h*Sk*Sq2 + 96*G3*Q2h*Sk*Sq2 + 96*G2*Q2k*Sk*Sq2 +
			     96*G3*Q2k*Sk*Sq2 - 8*G2*m2*pow(Q2h,2) - 8*G3*m2*pow(Q2h,2) +
			     16*G2*m2*pow(Q2k,2) + 16*G3*m2*pow(Q2k,2) + 16*G2*M2*pow(Q2k,2) +
			     16*G3*M2*pow(Q2k,2) + 32*G2*S*pow(Q2k,2) + 32*G3*S*pow(Q2k,2) -
			     16*G2*Sk*pow(Q2k,2) - 16*G3*Sk*pow(Q2k,2) + 32*G2*Sq2*pow(Q2k,2) +
			     32*G3*Sq2*pow(Q2k,2) - 128*G2*m2*pow(S,2) - 128*G3*m2*pow(S,2) +
			     96*G2*Q2e*pow(S,2) + 96*G3*Q2e*pow(S,2) - 96*G2*Q2h*pow(S,2) -
			     96*G3*Q2h*pow(S,2) - 256*G2*Q2k*pow(S,2) - 256*G3*Q2k*pow(S,2) +
			     64*G2*m2*pow(Sk,2) + 64*G3*m2*pow(Sk,2) + 32*G2*Q2e*pow(Sk,2) +
			     32*G3*Q2e*pow(Sk,2) - 32*G2*Q2h*pow(Sk,2) - 32*G3*Q2h*pow(Sk,2) -
			     32*G2*Q2k*pow(Sk,2) - 32*G3*Q2k*pow(Sk,2) - 64*G2*m2*pow(Sq2,2) -
			     64*G3*m2*pow(Sq2,2) + 64*G2*Q2e*pow(Sq2,2) + 64*G3*Q2e*pow(Sq2,2) -
			     64*G2*Q2h*pow(Sq2,2) - 64*G3*Q2h*pow(Sq2,2) - 96*G2*Q2k*pow(Sq2,2) -
			     96*G3*Q2k*pow(Sq2,2)) + pow(l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (32*G2*k1k2*m2*M2 + 32*G3*k1k2*m2*M2 + 64*G2*m4*M2 +
			     64*G3*m4*M2 + 16*G2*m2*M2*Q2e + 16*G3*m2*M2*Q2e +
			     32*G2*m4*Q2k + 32*G3*m4*Q2k + 80*G2*m2*M2*Q2k +
			     80*G3*m2*M2*Q2k + 8*G2*M2*Q2e*Q2k + 8*G3*M2*Q2e*Q2k +
			     32*G2*m2*Q2k*S + 32*G3*m2*Q2k*S - 64*G2*m4*Sk - 64*G3*m4*Sk -
			     16*G2*m2*Q2k*Sk - 16*G3*m2*Q2k*Sk + 128*G2*m2*S*Sk +
			     128*G3*m2*S*Sk - 32*G2*Q2e*S*Sk - 32*G3*Q2e*S*Sk + 32*G2*Q2k*S*Sk +
			     32*G3*Q2k*S*Sk + 64*G2*m2*Q2k*Sq2 + 64*G3*m2*Q2k*Sq2 -
			     128*G2*k1k2*S*Sq2 - 128*G3*k1k2*S*Sq2 - 256*G2*m2*S*Sq2 -
			     256*G3*m2*S*Sq2 - 192*G2*Q2k*S*Sq2 - 192*G3*Q2k*S*Sq2 +
			     128*G2*m2*Sk*Sq2 + 128*G3*m2*Sk*Sq2 - 32*G2*Q2e*Sk*Sq2 -
			     32*G3*Q2e*Sk*Sq2 + 32*G2*Q2k*Sk*Sq2 + 32*G3*Q2k*Sk*Sq2 +
			     16*G2*m2*pow(Q2k,2) + 16*G3*m2*pow(Q2k,2) + 8*G2*M2*pow(Q2k,2) +
			     8*G3*M2*pow(Q2k,2) + 16*G2*S*pow(Q2k,2) + 16*G3*S*pow(Q2k,2) +
			     16*G2*Sq2*pow(Q2k,2) + 16*G3*Sq2*pow(Q2k,2) - 64*G2*k1k2*pow(S,2) -
			     64*G3*k1k2*pow(S,2) - 128*G2*m2*pow(S,2) - 128*G3*m2*pow(S,2) -
			     128*G2*Q2k*pow(S,2) - 128*G3*Q2k*pow(S,2) + 32*G2*m2*pow(Sk,2) +
			     32*G3*m2*pow(Sk,2) - 64*G2*k1k2*pow(Sq2,2) - 64*G3*k1k2*pow(Sq2,2) -
			     128*G2*m2*pow(Sq2,2) - 128*G3*m2*pow(Sq2,2) - 64*G2*Q2k*pow(Sq2,2) -
			     64*G3*Q2k*pow(Sq2,2)) + pow(l1k1,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (-32*G2*m2*M2*Q2h - 32*G3*m2*M2*Q2h + 16*G2*k1k2*M2*Q2k +
			     16*G3*k1k2*M2*Q2k + 32*G2*m2*M2*Q2k + 32*G3*m2*M2*Q2k -
			     8*G2*m2*Q2h*Q2k - 8*G3*m2*Q2h*Q2k - 8*G2*M2*Q2h*Q2k -
			     8*G3*M2*Q2h*Q2k - 32*G2*m2*Q2h*S - 32*G3*m2*Q2h*S -
			     64*G2*m2*Q2k*S - 64*G3*m2*Q2k*S - 16*G2*Q2h*Q2k*S - 16*G3*Q2h*Q2k*S -
			     16*G2*m2*Q2h*Sk - 16*G3*m2*Q2h*Sk + 96*G2*m2*Q2k*Sk +
			     96*G3*m2*Q2k*Sk - 16*G2*Q2h*Q2k*Sk - 16*G3*Q2h*Q2k*Sk -
			     128*G2*m2*S*Sk - 128*G3*m2*S*Sk + 64*G2*Q2k*S*Sk + 64*G3*Q2k*S*Sk -
			     32*G2*m2*Q2h*Sq2 - 32*G3*m2*Q2h*Sq2 - 64*G2*m2*Q2k*Sq2 -
			     64*G3*m2*Q2k*Sq2 + 256*G2*m2*S*Sq2 + 256*G3*m2*S*Sq2 -
			     128*G2*Q2k*S*Sq2 - 128*G3*Q2k*S*Sq2 - 128*G2*m2*Sk*Sq2 -
			     128*G3*m2*Sk*Sq2 + 64*G2*Q2k*Sk*Sq2 + 64*G3*Q2k*Sk*Sq2 +
			     24*G2*m2*pow(Q2k,2) + 24*G3*m2*pow(Q2k,2) + 8*G2*M2*pow(Q2k,2) +
			     8*G3*M2*pow(Q2k,2) - 8*G2*Q2h*pow(Q2k,2) - 8*G3*Q2h*pow(Q2k,2) +
			     32*G2*S*pow(Q2k,2) + 32*G3*S*pow(Q2k,2) + 32*G2*Sq2*pow(Q2k,2) +
			     32*G3*Sq2*pow(Q2k,2) + 128*G2*m2*pow(S,2) + 128*G3*m2*pow(S,2) -
			     64*G2*Q2k*pow(S,2) - 64*G3*Q2k*pow(S,2) + 96*G2*m2*pow(Sk,2) +
			     96*G3*m2*pow(Sk,2) + 128*G2*m2*pow(Sq2,2) + 128*G3*m2*pow(Sq2,2) -
			     64*G2*Q2k*pow(Sq2,2) - 64*G3*Q2k*pow(Sq2,2)) +
			  pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (32*G2*l1k2*M2*Q2h + 32*G3*l1k2*M2*Q2h + 64*G2*m2*M2*Q2h +
			     64*G3*m2*M2*Q2h + 16*G2*M2*Q2h*Q2k + 16*G3*M2*Q2h*Q2k +
			     64*G2*l1k2*Q2h*S + 64*G3*l1k2*Q2h*S + 128*G2*m2*Q2h*S +
			     128*G3*m2*Q2h*S + 32*G2*Q2h*Q2k*S + 32*G3*Q2h*Q2k*S +
			     64*G2*l1k2*Q2h*Sq2 + 64*G3*l1k2*Q2h*Sq2 + 128*G2*m2*Q2h*Sq2 +
			     128*G3*m2*Q2h*Sq2 + 32*G2*Q2h*Q2k*Sq2 + 32*G3*Q2h*Q2k*Sq2 -
			     256*G2*l1k2*S*Sq2 - 256*G3*l1k2*S*Sq2 - 512*G2*m2*S*Sq2 -
			     512*G3*m2*S*Sq2 + 128*G2*Q2h*S*Sq2 + 128*G3*Q2h*S*Sq2 -
			     128*G2*Q2k*S*Sq2 - 128*G3*Q2k*S*Sq2 - 16*G2*M2*pow(Q2h,2) -
			     16*G3*M2*pow(Q2h,2) - 32*G2*S*pow(Q2h,2) - 32*G3*S*pow(Q2h,2) -
			     32*G2*Sq2*pow(Q2h,2) - 32*G3*Sq2*pow(Q2h,2) - 128*G2*l1k2*pow(S,2) -
			     128*G3*l1k2*pow(S,2) - 256*G2*m2*pow(S,2) - 256*G3*m2*pow(S,2) +
			     64*G2*Q2h*pow(S,2) + 64*G3*Q2h*pow(S,2) - 64*G2*Q2k*pow(S,2) -
			     64*G3*Q2k*pow(S,2) - 128*G2*l1k2*pow(Sq2,2) - 128*G3*l1k2*pow(Sq2,2) -
			     256*G2*m2*pow(Sq2,2) - 256*G3*m2*pow(Sq2,2) + 64*G2*Q2h*pow(Sq2,2) +
			     64*G3*Q2h*pow(Sq2,2) - 64*G2*Q2k*pow(Sq2,2) - 64*G3*Q2k*pow(Sq2,2)) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   (32*G2*m4*M2 + 32*G3*m4*M2 + 16*G2*m4*Q2h + 16*G3*m4*Q2h +
			     40*G2*m2*M2*Q2h + 40*G3*m2*M2*Q2h - 8*G2*M2*Q2e*Q2h -
			     8*G3*M2*Q2e*Q2h + 16*G2*m4*Q2k + 16*G3*m4*Q2k +
			     40*G2*m2*M2*Q2k + 40*G3*m2*M2*Q2k + 4*G2*M2*Q2e*Q2k +
			     4*G3*M2*Q2e*Q2k + 8*G2*m2*Q2h*Q2k + 8*G3*m2*Q2h*Q2k -
			     4*G2*M2*Q2h*Q2k - 4*G3*M2*Q2h*Q2k - 16*G2*Q2e*Q2k*S -
			     16*G3*Q2e*Q2k*S + 16*G2*Q2h*Q2k*S + 16*G3*Q2h*Q2k*S +
			     16*G2*m2*Q2k*Sk + 16*G3*m2*Q2k*Sk + 8*G2*Q2e*Q2k*Sk +
			     8*G3*Q2e*Q2k*Sk - 8*G2*Q2h*Q2k*Sk - 8*G3*Q2h*Q2k*Sk - 64*G2*Q2e*S*Sk -
			     64*G3*Q2e*S*Sk + 64*G2*Q2h*S*Sk + 64*G3*Q2h*S*Sk + 64*G2*Q2k*S*Sk +
			     64*G3*Q2k*S*Sk - 32*G2*m4*Sq2 - 32*G3*m4*Sq2 - 8*G2*m2*Q2h*Sq2 -
			     8*G3*m2*Q2h*Sq2 + 8*G2*Q2e*Q2h*Sq2 + 8*G3*Q2e*Q2h*Sq2 -
			     8*G2*Q2e*Q2k*Sq2 - 8*G3*Q2e*Q2k*Sq2 + 8*G2*Q2h*Q2k*Sq2 +
			     8*G3*Q2h*Q2k*Sq2 - 192*G2*m2*S*Sq2 - 192*G3*m2*S*Sq2 +
			     112*G2*Q2e*S*Sq2 + 112*G3*Q2e*S*Sq2 - 112*G2*Q2h*S*Sq2 -
			     112*G3*Q2h*S*Sq2 - 160*G2*Q2k*S*Sq2 - 160*G3*Q2k*S*Sq2 -
			     32*G2*m2*Sk*Sq2 - 32*G3*m2*Sk*Sq2 - 48*G2*Q2e*Sk*Sq2 -
			     48*G3*Q2e*Sk*Sq2 + 48*G2*Q2h*Sk*Sq2 + 48*G3*Q2h*Sk*Sq2 +
			     48*G2*Q2k*Sk*Sq2 + 48*G3*Q2k*Sk*Sq2 + 12*G2*m2*pow(Q2h,2) +
			     12*G3*m2*pow(Q2h,2) + 40*G2*M2*pow(Q2h,2) + 40*G3*M2*pow(Q2h,2) -
			     4*G2*Q2e*pow(Q2h,2) - 4*G3*Q2e*pow(Q2h,2) + 24*G2*Sq2*pow(Q2h,2) +
			     24*G3*Sq2*pow(Q2h,2) + 4*G2*pow(Q2h,3) + 4*G3*pow(Q2h,3) +
			     8*G2*m2*pow(Q2k,2) + 8*G3*m2*pow(Q2k,2) + 8*G2*M2*pow(Q2k,2) +
			     8*G3*M2*pow(Q2k,2) + 16*G2*S*pow(Q2k,2) + 16*G3*S*pow(Q2k,2) -
			     8*G2*Sk*pow(Q2k,2) - 8*G3*Sk*pow(Q2k,2) + 16*G2*Sq2*pow(Q2k,2) +
			     16*G3*Sq2*pow(Q2k,2) - 192*G2*m2*pow(S,2) - 192*G3*m2*pow(S,2) +
			     80*G2*Q2e*pow(S,2) + 80*G3*Q2e*pow(S,2) - 80*G2*Q2h*pow(S,2) -
			     80*G3*Q2h*pow(S,2) - 128*G2*Q2k*pow(S,2) - 128*G3*Q2k*pow(S,2) +
			     32*G2*m2*pow(Sk,2) + 32*G3*m2*pow(Sk,2) + 16*G2*Q2e*pow(Sk,2) +
			     16*G3*Q2e*pow(Sk,2) - 16*G2*Q2h*pow(Sk,2) - 16*G3*Q2h*pow(Sk,2) -
			     16*G2*Q2k*pow(Sk,2) - 16*G3*Q2k*pow(Sk,2) - 32*G2*m2*pow(Sq2,2) -
			     32*G3*m2*pow(Sq2,2) + 32*G2*Q2e*pow(Sq2,2) + 32*G3*Q2e*pow(Sq2,2) -
			     64*G2*Q2h*pow(Sq2,2) - 64*G3*Q2h*pow(Sq2,2) - 48*G2*Q2k*pow(Sq2,2) -
			     48*G3*Q2k*pow(Sq2,2)) + pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (64*G2*m4*M2 + 64*G3*m4*M2 + 32*G2*m4*Q2h + 32*G3*m4*Q2h +
			     80*G2*m2*M2*Q2h + 80*G3*m2*M2*Q2h - 16*G2*M2*Q2e*Q2h -
			     16*G3*M2*Q2e*Q2h + 32*G2*m4*Q2k + 32*G3*m4*Q2k +
			     80*G2*m2*M2*Q2k + 80*G3*m2*M2*Q2k + 8*G2*M2*Q2e*Q2k +
			     8*G3*M2*Q2e*Q2k + 16*G2*m2*Q2h*Q2k + 16*G3*m2*Q2h*Q2k -
			     8*G2*M2*Q2h*Q2k - 8*G3*M2*Q2h*Q2k + 32*G2*Q2e*Q2k*S +
			     32*G3*Q2e*Q2k*S - 32*G2*Q2h*Q2k*S - 32*G3*Q2h*Q2k*S +
			     32*G2*m2*Q2k*Sk + 32*G3*m2*Q2k*Sk + 16*G2*Q2e*Q2k*Sk +
			     16*G3*Q2e*Q2k*Sk - 16*G2*Q2h*Q2k*Sk - 16*G3*Q2h*Q2k*Sk +
			     128*G2*Q2e*S*Sk + 128*G3*Q2e*S*Sk - 128*G2*Q2h*S*Sk - 128*G3*Q2h*S*Sk -
			     128*G2*Q2k*S*Sk - 128*G3*Q2k*S*Sk - 64*G2*m4*Sq2 - 64*G3*m4*Sq2 -
			     16*G2*m2*Q2h*Sq2 - 16*G3*m2*Q2h*Sq2 + 16*G2*Q2e*Q2h*Sq2 +
			     16*G3*Q2e*Q2h*Sq2 + 16*G2*Q2e*Q2k*Sq2 + 16*G3*Q2e*Q2k*Sq2 -
			     16*G2*Q2h*Q2k*Sq2 - 16*G3*Q2h*Q2k*Sq2 - 384*G2*m2*S*Sq2 -
			     384*G3*m2*S*Sq2 + 96*G2*Q2e*S*Sq2 + 96*G3*Q2e*S*Sq2 -
			     96*G2*Q2h*S*Sq2 - 96*G3*Q2h*S*Sq2 - 192*G2*Q2k*S*Sq2 -
			     192*G3*Q2k*S*Sq2 - 64*G2*m2*Sk*Sq2 - 64*G3*m2*Sk*Sq2 +
			     32*G2*Q2e*Sk*Sq2 + 32*G3*Q2e*Sk*Sq2 - 32*G2*Q2h*Sk*Sq2 -
			     32*G3*Q2h*Sk*Sq2 - 32*G2*Q2k*Sk*Sq2 - 32*G3*Q2k*Sk*Sq2 +
			     24*G2*m2*pow(Q2h,2) + 24*G3*m2*pow(Q2h,2) + 80*G2*M2*pow(Q2h,2) +
			     80*G3*M2*pow(Q2h,2) - 8*G2*Q2e*pow(Q2h,2) - 8*G3*Q2e*pow(Q2h,2) +
			     48*G2*Sq2*pow(Q2h,2) + 48*G3*Sq2*pow(Q2h,2) + 8*G2*pow(Q2h,3) +
			     8*G3*pow(Q2h,3) + 16*G2*m2*pow(Q2k,2) + 16*G3*m2*pow(Q2k,2) +
			     16*G2*M2*pow(Q2k,2) + 16*G3*M2*pow(Q2k,2) - 32*G2*S*pow(Q2k,2) -
			     32*G3*S*pow(Q2k,2) - 16*G2*Sk*pow(Q2k,2) - 16*G3*Sk*pow(Q2k,2) -
			     384*G2*m2*pow(S,2) - 384*G3*m2*pow(S,2) + 160*G2*Q2e*pow(S,2) +
			     160*G3*Q2e*pow(S,2) - 160*G2*Q2h*pow(S,2) - 160*G3*Q2h*pow(S,2) -
			     256*G2*Q2k*pow(S,2) - 256*G3*Q2k*pow(S,2) + 64*G2*m2*pow(Sk,2) +
			     64*G3*m2*pow(Sk,2) + 32*G2*Q2e*pow(Sk,2) + 32*G3*Q2e*pow(Sk,2) -
			     32*G2*Q2h*pow(Sk,2) - 32*G3*Q2h*pow(Sk,2) - 32*G2*Q2k*pow(Sk,2) -
			     32*G3*Q2k*pow(Sk,2) - 64*G2*m2*pow(Sq2,2) - 64*G3*m2*pow(Sq2,2) -
			     64*G2*Q2h*pow(Sq2,2) - 64*G3*Q2h*pow(Sq2,2) - 32*G2*Q2k*pow(Sq2,2) -
			     32*G3*Q2k*pow(Sq2,2)) + pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   (-32*G2*m4*M2 - 32*G3*m4*M2 - 16*G2*m4*Q2h - 16*G3*m4*Q2h -
			     8*G2*m2*M2*Q2h - 8*G3*m2*M2*Q2h - 16*G2*m4*Q2k - 16*G3*m4*Q2k -
			     40*G2*m2*M2*Q2k - 40*G3*m2*M2*Q2k - 4*G2*M2*Q2e*Q2k -
			     4*G3*M2*Q2e*Q2k - 8*G2*m2*Q2h*Q2k - 8*G3*m2*Q2h*Q2k +
			     4*G2*M2*Q2h*Q2k + 4*G3*M2*Q2h*Q2k - 16*G2*Q2e*Q2k*S -
			     16*G3*Q2e*Q2k*S + 16*G2*Q2h*Q2k*S + 16*G3*Q2h*Q2k*S -
			     16*G2*m2*Q2k*Sk - 16*G3*m2*Q2k*Sk - 8*G2*Q2e*Q2k*Sk -
			     8*G3*Q2e*Q2k*Sk + 8*G2*Q2h*Q2k*Sk + 8*G3*Q2h*Q2k*Sk - 64*G2*Q2e*S*Sk -
			     64*G3*Q2e*S*Sk + 64*G2*Q2h*S*Sk + 64*G3*Q2h*S*Sk + 64*G2*Q2k*S*Sk +
			     64*G3*Q2k*S*Sk + 32*G2*m4*Sq2 + 32*G3*m4*Sq2 - 24*G2*m2*Q2h*Sq2 -
			     24*G3*m2*Q2h*Sq2 - 8*G2*Q2e*Q2k*Sq2 - 8*G3*Q2e*Q2k*Sq2 +
			     8*G2*Q2h*Q2k*Sq2 + 8*G3*Q2h*Q2k*Sq2 + 64*G2*m2*S*Sq2 +
			     64*G3*m2*S*Sq2 - 16*G2*Q2e*S*Sq2 - 16*G3*Q2e*S*Sq2 + 16*G2*Q2h*S*Sq2 +
			     16*G3*Q2h*S*Sq2 + 96*G2*Q2k*S*Sq2 + 96*G3*Q2k*S*Sq2 +
			     32*G2*m2*Sk*Sq2 + 32*G3*m2*Sk*Sq2 - 16*G2*Q2e*Sk*Sq2 -
			     16*G3*Q2e*Sk*Sq2 + 16*G2*Q2h*Sk*Sq2 + 16*G3*Q2h*Sk*Sq2 +
			     16*G2*Q2k*Sk*Sq2 + 16*G3*Q2k*Sk*Sq2 + 4*G2*m2*pow(Q2h,2) +
			     4*G3*m2*pow(Q2h,2) - 8*G2*m2*pow(Q2k,2) - 8*G3*m2*pow(Q2k,2) -
			     8*G2*M2*pow(Q2k,2) - 8*G3*M2*pow(Q2k,2) + 16*G2*S*pow(Q2k,2) +
			     16*G3*S*pow(Q2k,2) + 8*G2*Sk*pow(Q2k,2) + 8*G3*Sk*pow(Q2k,2) +
			     64*G2*m2*pow(S,2) + 64*G3*m2*pow(S,2) - 48*G2*Q2e*pow(S,2) -
			     48*G3*Q2e*pow(S,2) + 48*G2*Q2h*pow(S,2) + 48*G3*Q2h*pow(S,2) +
			     128*G2*Q2k*pow(S,2) + 128*G3*Q2k*pow(S,2) - 32*G2*m2*pow(Sk,2) -
			     32*G3*m2*pow(Sk,2) - 16*G2*Q2e*pow(Sk,2) - 16*G3*Q2e*pow(Sk,2) +
			     16*G2*Q2h*pow(Sk,2) + 16*G3*Q2h*pow(Sk,2) + 16*G2*Q2k*pow(Sk,2) +
			     16*G3*Q2k*pow(Sk,2) + 32*G2*m2*pow(Sq2,2) + 32*G3*m2*pow(Sq2,2) +
			     16*G2*Q2k*pow(Sq2,2) + 16*G3*Q2k*pow(Sq2,2)) +
			  pow(l1k1,-1)*pow(l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (16*G2*m2*M2*Q2h + 16*G3*m2*M2*Q2h + 8*G2*k1k2*M2*Q2k +
			     8*G3*k1k2*M2*Q2k - 16*G2*m2*M2*Q2k - 16*G3*m2*M2*Q2k +
			     12*G2*m2*Q2h*Q2k + 12*G3*m2*Q2h*Q2k + 4*G2*M2*Q2h*Q2k +
			     4*G3*M2*Q2h*Q2k - 16*G2*m2*Q2h*S - 16*G3*m2*Q2h*S +
			     32*G2*m2*Q2k*S + 32*G3*m2*Q2k*S - 8*G2*Q2h*Q2k*S - 8*G3*Q2h*Q2k*S +
			     24*G2*m2*Q2h*Sk + 24*G3*m2*Q2h*Sk - 48*G2*m2*Q2k*Sk -
			     48*G3*m2*Q2k*Sk + 64*G2*m2*S*Sk + 64*G3*m2*S*Sk - 32*G2*Q2k*S*Sk -
			     32*G3*Q2k*S*Sk - 16*G2*m2*Q2h*Sq2 - 16*G3*m2*Q2h*Sq2 +
			     32*G2*m2*Q2k*Sq2 + 32*G3*m2*Q2k*Sq2 - 128*G2*m2*S*Sq2 -
			     128*G3*m2*S*Sq2 + 64*G2*Q2k*S*Sq2 + 64*G3*Q2k*S*Sq2 +
			     64*G2*m2*Sk*Sq2 + 64*G3*m2*Sk*Sq2 - 32*G2*Q2k*Sk*Sq2 -
			     32*G3*Q2k*Sk*Sq2 - 12*G2*m2*pow(Q2k,2) - 12*G3*m2*pow(Q2k,2) -
			     4*G2*M2*pow(Q2k,2) - 4*G3*M2*pow(Q2k,2) - 16*G2*S*pow(Q2k,2) -
			     16*G3*S*pow(Q2k,2) - 16*G2*Sq2*pow(Q2k,2) - 16*G3*Sq2*pow(Q2k,2) -
			     64*G2*m2*pow(S,2) - 64*G3*m2*pow(S,2) + 32*G2*Q2k*pow(S,2) +
			     32*G3*Q2k*pow(S,2) - 48*G2*m2*pow(Sk,2) - 48*G3*m2*pow(Sk,2) -
			     64*G2*m2*pow(Sq2,2) - 64*G3*m2*pow(Sq2,2) + 32*G2*Q2k*pow(Sq2,2) +
			     32*G3*Q2k*pow(Sq2,2)) + pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-1)*(-32*G2*m4*M2 - 32*G3*m4*M2 -
			     16*G2*m4*Q2h - 16*G3*m4*Q2h - 8*G2*m2*M2*Q2h - 8*G3*m2*M2*Q2h +
			     4*G2*M2*Q2e*Q2h + 4*G3*M2*Q2e*Q2h - 16*G2*m4*Q2k - 16*G3*m4*Q2k -
			     24*G2*m2*M2*Q2k - 24*G3*m2*M2*Q2k + 12*G2*m2*Q2h*Q2k +
			     12*G3*m2*Q2h*Q2k - 8*G2*M2*Q2h*Q2k - 8*G3*M2*Q2h*Q2k +
			     4*G2*Q2e*Q2h*Q2k + 4*G3*Q2e*Q2h*Q2k - 16*G2*m2*Q2h*S -
			     16*G3*m2*Q2h*S - 16*G2*Q2e*Q2h*S - 16*G3*Q2e*Q2h*S + 32*G2*m2*Q2k*S +
			     32*G3*m2*Q2k*S + 8*G2*Q2e*Q2k*S + 8*G3*Q2e*Q2k*S - 8*G2*Q2h*Q2k*S -
			     8*G3*Q2h*Q2k*S + 40*G2*m2*Q2h*Sk + 40*G3*m2*Q2h*Sk +
			     8*G2*Q2e*Q2h*Sk + 8*G3*Q2e*Q2h*Sk - 16*G2*m2*Q2k*Sk -
			     16*G3*m2*Q2k*Sk - 8*G2*Q2e*Q2k*Sk - 8*G3*Q2e*Q2k*Sk +
			     8*G2*Q2h*Q2k*Sk + 8*G3*Q2h*Q2k*Sk + 64*G2*m2*S*Sk + 64*G3*m2*S*Sk +
			     48*G2*Q2e*S*Sk + 48*G3*Q2e*S*Sk - 48*G2*Q2h*S*Sk - 48*G3*Q2h*S*Sk -
			     64*G2*Q2k*S*Sk - 64*G3*Q2k*S*Sk + 32*G2*m4*Sq2 + 32*G3*m4*Sq2 -
			     24*G2*m2*Q2h*Sq2 - 24*G3*m2*Q2h*Sq2 - 16*G2*Q2e*Q2h*Sq2 -
			     16*G3*Q2e*Q2h*Sq2 + 8*G2*Q2e*Q2k*Sq2 + 8*G3*Q2e*Q2k*Sq2 -
			     16*G2*Q2h*Q2k*Sq2 - 16*G3*Q2h*Q2k*Sq2 + 64*G2*m2*S*Sq2 +
			     64*G3*m2*S*Sq2 - 80*G2*Q2e*S*Sq2 - 80*G3*Q2e*S*Sq2 + 80*G2*Q2h*S*Sq2 +
			     80*G3*Q2h*S*Sq2 + 128*G2*Q2k*S*Sq2 + 128*G3*Q2k*S*Sq2 +
			     32*G2*m2*Sk*Sq2 + 32*G3*m2*Sk*Sq2 + 48*G2*Q2e*Sk*Sq2 +
			     48*G3*Q2e*Sk*Sq2 - 48*G2*Q2h*Sk*Sq2 - 48*G3*Q2h*Sk*Sq2 -
			     48*G2*Q2k*Sk*Sq2 - 48*G3*Q2k*Sk*Sq2 - 8*G2*m2*pow(Q2h,2) -
			     8*G3*m2*pow(Q2h,2) - 4*G2*M2*pow(Q2h,2) - 4*G3*M2*pow(Q2h,2) -
			     4*G2*Q2k*pow(Q2h,2) - 4*G3*Q2k*pow(Q2h,2) + 8*G2*S*pow(Q2h,2) +
			     8*G3*S*pow(Q2h,2) - 8*G2*Sk*pow(Q2h,2) - 8*G3*Sk*pow(Q2h,2) +
			     8*G2*Sq2*pow(Q2h,2) + 8*G3*Sq2*pow(Q2h,2) - 8*G2*m2*pow(Q2k,2) -
			     8*G3*m2*pow(Q2k,2) - 8*G2*M2*pow(Q2k,2) - 8*G3*M2*pow(Q2k,2) -
			     16*G2*S*pow(Q2k,2) - 16*G3*S*pow(Q2k,2) + 8*G2*Sk*pow(Q2k,2) +
			     8*G3*Sk*pow(Q2k,2) - 16*G2*Sq2*pow(Q2k,2) - 16*G3*Sq2*pow(Q2k,2) +
			     64*G2*m2*pow(S,2) + 64*G3*m2*pow(S,2) - 48*G2*Q2e*pow(S,2) -
			     48*G3*Q2e*pow(S,2) + 32*G2*Q2h*pow(S,2) + 32*G3*Q2h*pow(S,2) +
			     80*G2*Q2k*pow(S,2) + 80*G3*Q2k*pow(S,2) - 32*G2*m2*pow(Sk,2) -
			     32*G3*m2*pow(Sk,2) - 16*G2*Q2e*pow(Sk,2) - 16*G3*Q2e*pow(Sk,2) +
			     16*G2*Q2h*pow(Sk,2) + 16*G3*Q2h*pow(Sk,2) + 16*G2*Q2k*pow(Sk,2) +
			     16*G3*Q2k*pow(Sk,2) + 32*G2*m2*pow(Sq2,2) + 32*G3*m2*pow(Sq2,2) -
			     32*G2*Q2e*pow(Sq2,2) - 32*G3*Q2e*pow(Sq2,2) + 48*G2*Q2h*pow(Sq2,2) +
			     48*G3*Q2h*pow(Sq2,2) + 48*G2*Q2k*pow(Sq2,2) + 48*G3*Q2k*pow(Sq2,2)) +
			  pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-2)*(64*G2*m4*M2*Q2k + 64*G3*m4*M2*Q2k +
			     128*G2*m4*Q2k*S + 128*G3*m4*Q2k*S - 64*G2*m4*Q2k*Sk -
			     64*G3*m4*Q2k*Sk + 256*G2*m4*S*Sk + 256*G3*m4*S*Sk +
			     128*G2*m2*Q2k*S*Sk + 128*G3*m2*Q2k*S*Sk + 128*G2*m4*Q2k*Sq2 +
			     128*G3*m4*Q2k*Sq2 - 512*G2*m4*S*Sq2 - 512*G3*m4*S*Sq2 -
			     256*G2*m2*Q2k*S*Sq2 - 256*G3*m2*Q2k*S*Sq2 + 256*G2*m4*Sk*Sq2 +
			     256*G3*m4*Sk*Sq2 + 128*G2*m2*Q2k*Sk*Sq2 + 128*G3*m2*Q2k*Sk*Sq2 +
			     32*G2*m2*M2*pow(Q2k,2) + 32*G3*m2*M2*pow(Q2k,2) +
			     64*G2*m2*S*pow(Q2k,2) + 64*G3*m2*S*pow(Q2k,2) -
			     32*G2*m2*Sk*pow(Q2k,2) - 32*G3*m2*Sk*pow(Q2k,2) +
			     64*G2*m2*Sq2*pow(Q2k,2) + 64*G3*m2*Sq2*pow(Q2k,2) -
			     256*G2*m4*pow(S,2) - 256*G3*m4*pow(S,2) - 128*G2*m2*Q2k*pow(S,2) -
			     128*G3*m2*Q2k*pow(S,2) - 64*G2*m2*Q2k*pow(Sk,2) -
			     64*G3*m2*Q2k*pow(Sk,2) - 256*G2*m4*pow(Sq2,2) -
			     256*G3*m4*pow(Sq2,2) - 128*G2*m2*Q2k*pow(Sq2,2) -
			     128*G3*m2*Q2k*pow(Sq2,2)) +
			  pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-1)*(-16*G2*m4*M2*Q2h - 16*G3*m4*M2*Q2h -
			     8*G2*m4*Q2h*Q2k - 8*G3*m4*Q2h*Q2k - 8*G2*m2*M2*Q2h*Q2k -
			     8*G3*m2*M2*Q2h*Q2k + 32*G2*m4*Q2h*S + 32*G3*m4*Q2h*S -
			     32*G2*m4*Q2k*S - 32*G3*m4*Q2k*S + 16*G2*m2*Q2h*Q2k*S +
			     16*G3*m2*Q2h*Q2k*S - 16*G2*m4*Q2h*Sk - 16*G3*m4*Q2h*Sk -
			     16*G2*m2*Q2h*Q2k*Sk - 16*G3*m2*Q2h*Q2k*Sk - 64*G2*m4*S*Sk -
			     64*G3*m4*S*Sk - 32*G2*m2*Q2k*S*Sk - 32*G3*m2*Q2k*S*Sk +
			     16*G2*m4*Q2h*Sq2 + 16*G3*m4*Q2h*Sq2 + 16*G2*m2*Q2h*Q2k*Sq2 +
			     16*G3*m2*Q2h*Q2k*Sq2 + 64*G2*m4*S*Sq2 + 64*G3*m4*S*Sq2 +
			     32*G2*m2*Q2k*S*Sq2 + 32*G3*m2*Q2k*S*Sq2 + 16*G2*m2*Q2h*Sk*Sq2 +
			     16*G3*m2*Q2h*Sk*Sq2 - 48*G2*m2*Q2k*Sk*Sq2 - 48*G3*m2*Q2k*Sk*Sq2 -
			     6*G2*m2*Q2h*pow(Q2k,2) - 6*G3*m2*Q2h*pow(Q2k,2) -
			     16*G2*m2*S*pow(Q2k,2) - 16*G3*m2*S*pow(Q2k,2) +
			     24*G2*m2*Sk*pow(Q2k,2) + 24*G3*m2*Sk*pow(Q2k,2) -
			     24*G2*m2*Sq2*pow(Q2k,2) - 24*G3*m2*Sq2*pow(Q2k,2) +
			     6*G2*m2*pow(Q2k,3) + 6*G3*m2*pow(Q2k,3) + 64*G2*m4*pow(S,2) +
			     64*G3*m4*pow(S,2) + 32*G2*m2*Q2k*pow(S,2) + 32*G3*m2*Q2k*pow(S,2) -
			     8*G2*m2*Q2h*pow(Sk,2) - 8*G3*m2*Q2h*pow(Sk,2) +
			     24*G2*m2*Q2k*pow(Sk,2) + 24*G3*m2*Q2k*pow(Sk,2) -
			     8*G2*m2*Q2h*pow(Sq2,2) - 8*G3*m2*Q2h*pow(Sq2,2) +
			     24*G2*m2*Q2k*pow(Sq2,2) + 24*G3*m2*Q2k*pow(Sq2,2)) +
			  pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*l1k1 + Q2e - Q2k,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-64*G2*m4*M2*Q2h - 64*G3*m4*M2*Q2h - 32*G2*m4*Q2h*Q2k -
			     32*G3*m4*Q2h*Q2k - 32*G2*m2*M2*Q2h*Q2k - 32*G3*m2*M2*Q2h*Q2k -
			     128*G2*m4*Q2h*S - 128*G3*m4*Q2h*S + 128*G2*m4*Q2k*S +
			     128*G3*m4*Q2k*S - 64*G2*m2*Q2h*Q2k*S - 64*G3*m2*Q2h*Q2k*S -
			     64*G2*m4*Q2h*Sk - 64*G3*m4*Q2h*Sk - 64*G2*m2*Q2h*Q2k*Sk -
			     64*G3*m2*Q2h*Q2k*Sk + 256*G2*m4*S*Sk + 256*G3*m4*S*Sk +
			     128*G2*m2*Q2k*S*Sk + 128*G3*m2*Q2k*S*Sk - 64*G2*m4*Q2h*Sq2 -
			     64*G3*m4*Q2h*Sq2 + 128*G2*m4*Q2k*Sq2 + 128*G3*m4*Q2k*Sq2 +
			     256*G2*m4*S*Sq2 + 256*G3*m4*S*Sq2 + 128*G2*m2*Q2k*S*Sq2 +
			     128*G3*m2*Q2k*S*Sq2 + 256*G2*m4*Sk*Sq2 + 256*G3*m4*Sk*Sq2 +
			     64*G2*m2*Q2h*Sk*Sq2 + 64*G3*m2*Q2h*Sk*Sq2 - 64*G2*m2*Q2k*Sk*Sq2 -
			     64*G3*m2*Q2k*Sk*Sq2 - 24*G2*m2*Q2h*pow(Q2k,2) -
			     24*G3*m2*Q2h*pow(Q2k,2) + 64*G2*m2*S*pow(Q2k,2) +
			     64*G3*m2*S*pow(Q2k,2) + 96*G2*m2*Sk*pow(Q2k,2) +
			     96*G3*m2*Sk*pow(Q2k,2) - 32*G2*m2*Sq2*pow(Q2k,2) -
			     32*G3*m2*Sq2*pow(Q2k,2) + 24*G2*m2*pow(Q2k,3) +
			     24*G3*m2*pow(Q2k,3) + 256*G2*m4*pow(S,2) + 256*G3*m4*pow(S,2) +
			     128*G2*m2*Q2k*pow(S,2) + 128*G3*m2*Q2k*pow(S,2) -
			     32*G2*m2*Q2h*pow(Sk,2) - 32*G3*m2*Q2h*pow(Sk,2) +
			     96*G2*m2*Q2k*pow(Sk,2) + 96*G3*m2*Q2k*pow(Sk,2) -
			     32*G2*m2*Q2h*pow(Sq2,2) - 32*G3*m2*Q2h*pow(Sq2,2) +
			     96*G2*m2*Q2k*pow(Sq2,2) + 96*G3*m2*Q2k*pow(Sq2,2)) +
			  pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (32*G2*m4*M2*Q2h + 32*G3*m4*M2*Q2h + 48*G2*m2*M2*Q2h*Q2k +
			     48*G3*m2*M2*Q2h*Q2k + 16*G2*m2*Q2h*Q2k*Sk + 16*G3*m2*Q2h*Q2k*Sk -
			     32*G2*m4*Q2h*Sq2 - 32*G3*m4*Q2h*Sq2 + 8*G2*m2*Q2h*Q2k*Sq2 +
			     8*G3*m2*Q2h*Q2k*Sq2 - 128*G2*m4*S*Sq2 - 128*G3*m4*S*Sq2 -
			     128*G2*m2*Q2k*S*Sq2 - 128*G3*m2*Q2k*S*Sq2 - 64*G2*Q2h*Q2k*S*Sq2 -
			     64*G3*Q2h*Q2k*S*Sq2 - 32*G2*m2*Q2h*Sk*Sq2 - 32*G3*m2*Q2h*Sk*Sq2 +
			     16*G2*m4*pow(Q2h,2) + 16*G3*m4*pow(Q2h,2) -
			     8*G2*m2*M2*pow(Q2h,2) - 8*G3*m2*M2*pow(Q2h,2) +
			     8*G2*m2*Q2k*pow(Q2h,2) + 8*G3*m2*Q2k*pow(Q2h,2) -
			     8*G2*m2*Sq2*pow(Q2h,2) - 8*G3*m2*Sq2*pow(Q2h,2) -
			     8*G2*m2*M2*pow(Q2k,2) - 8*G3*m2*M2*pow(Q2k,2) +
			     8*G2*M2*Q2h*pow(Q2k,2) + 8*G3*M2*Q2h*pow(Q2k,2) +
			     8*G2*Q2h*Sq2*pow(Q2k,2) + 8*G3*Q2h*Sq2*pow(Q2k,2) +
			     32*G2*S*Sq2*pow(Q2k,2) + 32*G3*S*Sq2*pow(Q2k,2) - 128*G2*m4*pow(S,2) -
			     128*G3*m4*pow(S,2) - 128*G2*m2*Q2k*pow(S,2) -
			     128*G3*m2*Q2k*pow(S,2) - 64*G2*Q2h*Q2k*pow(S,2) -
			     64*G3*Q2h*Q2k*pow(S,2) + 32*G2*pow(Q2k,2)*pow(S,2) +
			     32*G3*pow(Q2k,2)*pow(S,2) + 32*G2*m2*Q2h*pow(Sk,2) +
			     32*G3*m2*Q2h*pow(Sk,2) + 16*G2*m2*Q2h*pow(Sq2,2) +
			     16*G3*m2*Q2h*pow(Sq2,2) - 48*G2*m2*Q2k*pow(Sq2,2) -
			     48*G3*m2*Q2k*pow(Sq2,2) - 16*G2*Q2h*Q2k*pow(Sq2,2) -
			     16*G3*Q2h*Q2k*pow(Sq2,2)) +
			  pow(l1k1,-1)*pow(l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (32*G2*m4*M2*Q2h + 32*G3*m4*M2*Q2h - 32*G2*m4*Q2h*Q2k -
			     32*G3*m4*Q2h*Q2k + 48*G2*m2*M2*Q2h*Q2k + 48*G3*m2*M2*Q2h*Q2k -
			     64*G2*m4*Q2h*Sk - 64*G3*m4*Q2h*Sk + 128*G2*m4*Q2k*Sk +
			     128*G3*m4*Q2k*Sk + 64*G2*m2*Q2h*Q2k*Sk + 64*G3*m2*Q2h*Q2k*Sk +
			     32*G2*m4*Q2h*Sq2 + 32*G3*m4*Q2h*Sq2 - 64*G2*m4*Q2k*Sq2 -
			     64*G3*m4*Q2k*Sq2 - 40*G2*m2*Q2h*Q2k*Sq2 - 40*G3*m2*Q2h*Q2k*Sq2 -
			     128*G2*m4*S*Sq2 - 128*G3*m4*S*Sq2 - 128*G2*m2*Q2k*S*Sq2 -
			     128*G3*m2*Q2k*S*Sq2 - 64*G2*Q2h*Q2k*S*Sq2 - 64*G3*Q2h*Q2k*S*Sq2 -
			     128*G2*m4*Sk*Sq2 - 128*G3*m4*Sk*Sq2 - 32*G2*m2*Q2h*Sk*Sq2 -
			     32*G3*m2*Q2h*Sk*Sq2 + 64*G2*m2*Q2k*Sk*Sq2 + 64*G3*m2*Q2k*Sk*Sq2 +
			     16*G2*m4*pow(Q2h,2) + 16*G3*m4*pow(Q2h,2) -
			     8*G2*m2*M2*pow(Q2h,2) - 8*G3*m2*M2*pow(Q2h,2) -
			     16*G2*m2*Sk*pow(Q2h,2) - 16*G3*m2*Sk*pow(Q2h,2) +
			     8*G2*m2*Sq2*pow(Q2h,2) + 8*G3*m2*Sq2*pow(Q2h,2) +
			     32*G2*m4*pow(Q2k,2) + 32*G3*m4*pow(Q2k,2) -
			     8*G2*m2*M2*pow(Q2k,2) - 8*G3*m2*M2*pow(Q2k,2) +
			     24*G2*m2*Q2h*pow(Q2k,2) + 24*G3*m2*Q2h*pow(Q2k,2) +
			     8*G2*M2*Q2h*pow(Q2k,2) + 8*G3*M2*Q2h*pow(Q2k,2) -
			     64*G2*m2*Sk*pow(Q2k,2) - 64*G3*m2*Sk*pow(Q2k,2) +
			     16*G2*Q2h*Sk*pow(Q2k,2) + 16*G3*Q2h*Sk*pow(Q2k,2) +
			     32*G2*m2*Sq2*pow(Q2k,2) + 32*G3*m2*Sq2*pow(Q2k,2) -
			     8*G2*Q2h*Sq2*pow(Q2k,2) - 8*G3*Q2h*Sq2*pow(Q2k,2) +
			     32*G2*S*Sq2*pow(Q2k,2) + 32*G3*S*Sq2*pow(Q2k,2) +
			     32*G2*Sk*Sq2*pow(Q2k,2) + 32*G3*Sk*Sq2*pow(Q2k,2) -
			     16*G2*m2*pow(Q2k,3) - 16*G3*m2*pow(Q2k,3) + 8*G2*Q2h*pow(Q2k,3) +
			     8*G3*Q2h*pow(Q2k,3) - 32*G2*Sk*pow(Q2k,3) - 32*G3*Sk*pow(Q2k,3) +
			     16*G2*Sq2*pow(Q2k,3) + 16*G3*Sq2*pow(Q2k,3) - 8*G2*pow(Q2k,4) -
			     8*G3*pow(Q2k,4) - 128*G2*m4*pow(S,2) - 128*G3*m4*pow(S,2) -
			     128*G2*m2*Q2k*pow(S,2) - 128*G3*m2*Q2k*pow(S,2) -
			     64*G2*Q2h*Q2k*pow(S,2) - 64*G3*Q2h*Q2k*pow(S,2) +
			     32*G2*pow(Q2k,2)*pow(S,2) + 32*G3*pow(Q2k,2)*pow(S,2) +
			     128*G2*m4*pow(Sk,2) + 128*G3*m4*pow(Sk,2) + 32*G2*m2*Q2h*pow(Sk,2) +
			     32*G3*m2*Q2h*pow(Sk,2) - 64*G2*m2*Q2k*pow(Sk,2) -
			     64*G3*m2*Q2k*pow(Sk,2) - 32*G2*pow(Q2k,2)*pow(Sk,2) -
			     32*G3*pow(Q2k,2)*pow(Sk,2) + 16*G2*m2*Q2h*pow(Sq2,2) +
			     16*G3*m2*Q2h*pow(Sq2,2) - 48*G2*m2*Q2k*pow(Sq2,2) -
			     48*G3*m2*Q2k*pow(Sq2,2) - 16*G2*Q2h*Q2k*pow(Sq2,2) -
			     16*G3*Q2h*Q2k*pow(Sq2,2)) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-32*G2*m4*M2*Q2h - 32*G3*m4*M2*Q2h + 16*G2*m2*M2*Q2h*Q2k +
			     16*G3*m2*M2*Q2h*Q2k - 16*G2*m2*Q2h*Q2k*Sk - 16*G3*m2*Q2h*Q2k*Sk +
			     32*G2*m4*Q2h*Sq2 + 32*G3*m4*Q2h*Sq2 + 24*G2*m2*Q2h*Q2k*Sq2 +
			     24*G3*m2*Q2h*Q2k*Sq2 + 128*G2*m4*S*Sq2 + 128*G3*m4*S*Sq2 +
			     256*G2*m2*Q2h*S*Sq2 + 256*G3*m2*Q2h*S*Sq2 - 128*G2*m2*Q2k*S*Sq2 -
			     128*G3*m2*Q2k*S*Sq2 + 64*G2*Q2h*Q2k*S*Sq2 + 64*G3*Q2h*Q2k*S*Sq2 +
			     32*G2*m2*Q2h*Sk*Sq2 + 32*G3*m2*Q2h*Sk*Sq2 - 16*G2*m4*pow(Q2h,2) -
			     16*G3*m4*pow(Q2h,2) - 56*G2*m2*M2*pow(Q2h,2) -
			     56*G3*m2*M2*pow(Q2h,2) + 32*G2*M2*Q2k*pow(Q2h,2) +
			     32*G3*M2*Q2k*pow(Q2h,2) - 24*G2*m2*Sq2*pow(Q2h,2) -
			     24*G3*m2*Sq2*pow(Q2h,2) + 32*G2*Q2k*Sq2*pow(Q2h,2) +
			     32*G3*Q2k*Sq2*pow(Q2h,2) - 8*G2*m2*pow(Q2h,3) - 8*G3*m2*pow(Q2h,3) -
			     32*G2*M2*pow(Q2h,3) - 32*G3*M2*pow(Q2h,3) - 32*G2*Sq2*pow(Q2h,3) -
			     32*G3*Sq2*pow(Q2h,3) + 8*G2*m2*M2*pow(Q2k,2) +
			     8*G3*m2*M2*pow(Q2k,2) - 8*G2*M2*Q2h*pow(Q2k,2) -
			     8*G3*M2*Q2h*pow(Q2k,2) - 8*G2*Q2h*Sq2*pow(Q2k,2) -
			     8*G3*Q2h*Sq2*pow(Q2k,2) - 32*G2*S*Sq2*pow(Q2k,2) -
			     32*G3*S*Sq2*pow(Q2k,2) + 128*G2*m4*pow(S,2) + 128*G3*m4*pow(S,2) +
			     256*G2*m2*Q2h*pow(S,2) + 256*G3*m2*Q2h*pow(S,2) -
			     128*G2*m2*Q2k*pow(S,2) - 128*G3*m2*Q2k*pow(S,2) +
			     64*G2*Q2h*Q2k*pow(S,2) + 64*G3*Q2h*Q2k*pow(S,2) -
			     32*G2*pow(Q2k,2)*pow(S,2) - 32*G3*pow(Q2k,2)*pow(S,2) -
			     32*G2*m2*Q2h*pow(Sk,2) - 32*G3*m2*Q2h*pow(Sk,2) +
			     80*G2*m2*Q2h*pow(Sq2,2) + 80*G3*m2*Q2h*pow(Sq2,2) -
			     48*G2*m2*Q2k*pow(Sq2,2) - 48*G3*m2*Q2k*pow(Sq2,2) -
			     16*G2*Q2h*Q2k*pow(Sq2,2) - 16*G3*Q2h*Q2k*pow(Sq2,2) +
			     32*G2*pow(Q2h,2)*pow(Sq2,2) + 32*G3*pow(Q2h,2)*pow(Sq2,2)) +
			  pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
			   pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
			   (32*G2*m4*M2*Q2k + 32*G3*m4*M2*Q2k - 64*G2*Q2e*Q2k*S*Sk -
			     64*G3*Q2e*Q2k*S*Sk - 32*G2*m4*Q2k*Sq2 - 32*G3*m4*Q2k*Sq2 -
			     128*G2*m4*S*Sq2 - 128*G3*m4*S*Sq2 - 128*G2*m2*Q2k*S*Sq2 -
			     128*G3*m2*Q2k*S*Sq2 + 96*G2*Q2e*Q2k*S*Sq2 + 96*G3*Q2e*Q2k*S*Sq2 -
			     32*G2*m2*Q2k*Sk*Sq2 - 32*G3*m2*Q2k*Sk*Sq2 - 48*G2*Q2e*Q2k*Sk*Sq2 -
			     48*G3*Q2e*Q2k*Sk*Sq2 + 16*G2*m4*pow(Q2k,2) + 16*G3*m4*pow(Q2k,2) +
			     32*G2*m2*M2*pow(Q2k,2) + 32*G3*m2*M2*pow(Q2k,2) -
			     16*G2*Q2e*S*pow(Q2k,2) - 16*G3*Q2e*S*pow(Q2k,2) +
			     16*G2*m2*Sk*pow(Q2k,2) + 16*G3*m2*Sk*pow(Q2k,2) +
			     8*G2*Q2e*Sk*pow(Q2k,2) + 8*G3*Q2e*Sk*pow(Q2k,2) +
			     64*G2*S*Sk*pow(Q2k,2) + 64*G3*S*Sk*pow(Q2k,2) -
			     8*G2*Q2e*Sq2*pow(Q2k,2) - 8*G3*Q2e*Sq2*pow(Q2k,2) -
			     128*G2*S*Sq2*pow(Q2k,2) - 128*G3*S*Sq2*pow(Q2k,2) +
			     48*G2*Sk*Sq2*pow(Q2k,2) + 48*G3*Sk*Sq2*pow(Q2k,2) +
			     8*G2*m2*pow(Q2k,3) + 8*G3*m2*pow(Q2k,3) + 8*G2*M2*pow(Q2k,3) +
			     8*G3*M2*pow(Q2k,3) + 16*G2*S*pow(Q2k,3) + 16*G3*S*pow(Q2k,3) -
			     8*G2*Sk*pow(Q2k,3) - 8*G3*Sk*pow(Q2k,3) + 16*G2*Sq2*pow(Q2k,3) +
			     16*G3*Sq2*pow(Q2k,3) - 128*G2*m4*pow(S,2) - 128*G3*m4*pow(S,2) -
			     128*G2*m2*Q2k*pow(S,2) - 128*G3*m2*Q2k*pow(S,2) +
			     64*G2*Q2e*Q2k*pow(S,2) + 64*G3*Q2e*Q2k*pow(S,2) -
			     96*G2*pow(Q2k,2)*pow(S,2) - 96*G3*pow(Q2k,2)*pow(S,2) +
			     32*G2*m2*Q2k*pow(Sk,2) + 32*G3*m2*Q2k*pow(Sk,2) +
			     16*G2*Q2e*Q2k*pow(Sk,2) + 16*G3*Q2e*Q2k*pow(Sk,2) -
			     16*G2*pow(Q2k,2)*pow(Sk,2) - 16*G3*pow(Q2k,2)*pow(Sk,2) -
			     32*G2*m2*Q2k*pow(Sq2,2) - 32*G3*m2*Q2k*pow(Sq2,2) +
			     32*G2*Q2e*Q2k*pow(Sq2,2) + 32*G3*Q2e*Q2k*pow(Sq2,2) -
			     48*G2*pow(Q2k,2)*pow(Sq2,2) - 48*G3*pow(Q2k,2)*pow(Sq2,2)) +
			  pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
			   pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
			   (-32*G2*m4*M2*Q2k - 32*G3*m4*M2*Q2k - 64*G2*Q2e*Q2k*S*Sk -
			     64*G3*Q2e*Q2k*S*Sk + 32*G2*m4*Q2k*Sq2 + 32*G3*m4*Q2k*Sq2 +
			     128*G2*m4*S*Sq2 + 128*G3*m4*S*Sq2 + 128*G2*m2*Q2k*S*Sq2 +
			     128*G3*m2*Q2k*S*Sq2 - 32*G2*Q2e*Q2k*S*Sq2 - 32*G3*Q2e*Q2k*S*Sq2 +
			     32*G2*m2*Q2k*Sk*Sq2 + 32*G3*m2*Q2k*Sk*Sq2 - 16*G2*Q2e*Q2k*Sk*Sq2 -
			     16*G3*Q2e*Q2k*Sk*Sq2 - 16*G2*m4*pow(Q2k,2) - 16*G3*m4*pow(Q2k,2) -
			     32*G2*m2*M2*pow(Q2k,2) - 32*G3*m2*M2*pow(Q2k,2) -
			     16*G2*Q2e*S*pow(Q2k,2) - 16*G3*Q2e*S*pow(Q2k,2) -
			     16*G2*m2*Sk*pow(Q2k,2) - 16*G3*m2*Sk*pow(Q2k,2) -
			     8*G2*Q2e*Sk*pow(Q2k,2) - 8*G3*Q2e*Sk*pow(Q2k,2) +
			     64*G2*S*Sk*pow(Q2k,2) + 64*G3*S*Sk*pow(Q2k,2) -
			     8*G2*Q2e*Sq2*pow(Q2k,2) - 8*G3*Q2e*Sq2*pow(Q2k,2) +
			     64*G2*S*Sq2*pow(Q2k,2) + 64*G3*S*Sq2*pow(Q2k,2) +
			     16*G2*Sk*Sq2*pow(Q2k,2) + 16*G3*Sk*Sq2*pow(Q2k,2) -
			     8*G2*m2*pow(Q2k,3) - 8*G3*m2*pow(Q2k,3) - 8*G2*M2*pow(Q2k,3) -
			     8*G3*M2*pow(Q2k,3) + 16*G2*S*pow(Q2k,3) + 16*G3*S*pow(Q2k,3) +
			     8*G2*Sk*pow(Q2k,3) + 8*G3*Sk*pow(Q2k,3) + 128*G2*m4*pow(S,2) +
			     128*G3*m4*pow(S,2) + 128*G2*m2*Q2k*pow(S,2) +
			     128*G3*m2*Q2k*pow(S,2) - 64*G2*Q2e*Q2k*pow(S,2) -
			     64*G3*Q2e*Q2k*pow(S,2) + 96*G2*pow(Q2k,2)*pow(S,2) +
			     96*G3*pow(Q2k,2)*pow(S,2) - 32*G2*m2*Q2k*pow(Sk,2) -
			     32*G3*m2*Q2k*pow(Sk,2) - 16*G2*Q2e*Q2k*pow(Sk,2) -
			     16*G3*Q2e*Q2k*pow(Sk,2) + 16*G2*pow(Q2k,2)*pow(Sk,2) +
			     16*G3*pow(Q2k,2)*pow(Sk,2) + 32*G2*m2*Q2k*pow(Sq2,2) +
			     32*G3*m2*Q2k*pow(Sq2,2) + 16*G2*pow(Q2k,2)*pow(Sq2,2) +
			     16*G3*pow(Q2k,2)*pow(Sq2,2));
}

long double Melem::melem2_R2(const long double p1l1, const long double p1l2, const long double p1k1, const long double p1k2,
		const long double p1p2, const long double l1p2, const long double l2p2, const long double k1p2, const long double p2k2,
		const long double l1l2, const long double l1k1, const long double l2k1, const long double l1k2, const long double l2k2,
		const long double k1k2, const long double eg1, const long double f1, const long double f2)const {

	return -8*eg1*pow(l1k2,-1)*pow(l2k1,-2)/M2*
			   (4*f1*f2*M2*(k1k2*(l2k1 - m2)*(M2 - p1p2) +
			        l2k2*m2*(-M2 + p1p2) +
			        (l2p2*m2 + k1p2*(-l2k1 + m2) + l2k1*p1k1 - m2*p1k1 - m2*p1l2)*
			         (p1k2 - p2k2)) + 4*M2*
			      (k1k2*(l2k1 - m2)*M2 - l2k2*m2*M2 - k1p2*l2k1*p1k2 +
			        k1p2*m2*p1k2 + l2p2*m2*p1k2 - l2k1*p1k1*p2k2 + m2*p1k1*p2k2 +
			        m2*p1l2*p2k2)*pow(f1,2) +
			     pow(f2,2)*(-(k1p2*l2k1*M2*p1k2) + k1p2*m2*M2*p1k2 +
			        l2p2*m2*M2*p1k2 + 3*l2k1*M2*p1k1*p1k2 - 3*m2*M2*p1k1*p1k2 -
			        3*m2*M2*p1k2*p1l2 - k1p2*l2k1*p1k2*p1p2 + k1p2*m2*p1k2*p1p2 +
			        l2p2*m2*p1k2*p1p2 - l2k1*p1k1*p1k2*p1p2 + m2*p1k1*p1k2*p1p2 +
			        m2*p1k2*p1l2*p1p2 + 3*k1p2*l2k1*M2*p2k2 - 3*k1p2*m2*M2*p2k2 -
			        3*l2p2*m2*M2*p2k2 - l2k1*M2*p1k1*p2k2 + m2*M2*p1k1*p2k2 +
			        m2*M2*p1l2*p2k2 - k1p2*l2k1*p1p2*p2k2 + k1p2*m2*p1p2*p2k2 +
			        l2p2*m2*p1p2*p2k2 - l2k1*p1k1*p1p2*p2k2 + m2*p1k1*p1p2*p2k2 +
			        m2*p1l2*p1p2*p2k2 + k1k2*(l2k1 - m2)*pow(M2 - p1p2,2) -
			        l2k2*m2*pow(M2 - p1p2,2))) -
			  8*eg1*pow(l1k1,-2)*pow(l2k2,-1)/M2*
			   (4*f1*f2*M2*(k1k2*(l1k1 + m2)*(M2 - p1p2) +
			        l1k2*m2*(-M2 + p1p2) +
			        (l1p2*m2 - k1p2*(l1k1 + m2) + l1k1*p1k1 + m2*p1k1 - m2*p1l1)*
			         (p1k2 - p2k2)) + 4*M2*
			      (-(l1k2*m2*M2) + k1k2*(l1k1 + m2)*M2 - k1p2*l1k1*p1k2 -
			        k1p2*m2*p1k2 + l1p2*m2*p1k2 - l1k1*p1k1*p2k2 - m2*p1k1*p2k2 +
			        m2*p1l1*p2k2)*pow(f1,2) +
			     pow(f2,2)*(-(k1p2*l1k1*M2*p1k2) - k1p2*m2*M2*p1k2 +
			        l1p2*m2*M2*p1k2 + 3*l1k1*M2*p1k1*p1k2 + 3*m2*M2*p1k1*p1k2 -
			        3*m2*M2*p1k2*p1l1 - k1p2*l1k1*p1k2*p1p2 - k1p2*m2*p1k2*p1p2 +
			        l1p2*m2*p1k2*p1p2 - l1k1*p1k1*p1k2*p1p2 - m2*p1k1*p1k2*p1p2 +
			        m2*p1k2*p1l1*p1p2 + 3*k1p2*l1k1*M2*p2k2 + 3*k1p2*m2*M2*p2k2 -
			        3*l1p2*m2*M2*p2k2 - l1k1*M2*p1k1*p2k2 - m2*M2*p1k1*p2k2 +
			        m2*M2*p1l1*p2k2 - k1p2*l1k1*p1p2*p2k2 - k1p2*m2*p1p2*p2k2 +
			        l1p2*m2*p1p2*p2k2 - l1k1*p1k1*p1p2*p2k2 - m2*p1k1*p1p2*p2k2 +
			        m2*p1l1*p1p2*p2k2 - l1k2*m2*pow(M2 - p1p2,2) +
			        k1k2*(l1k1 + m2)*pow(M2 - p1p2,2))) +
			  8*eg1*pow(l2k2,-1)*pow(l2k1 + eg1*(k1k2 + l2k2),-2)/M2*
			   (-4*f1*f2*M2*(l1p2*l2k2*l2p2 - l1k2*l2k1*M2 + l1k1*l2k2*M2 +
			        l1l2*l2k2*M2 - 3*l2k2*m2*M2 - l1p2*l2k2*p1k1 +
			        l1p2*l2k1*p1k2 + k1p2*l2k2*(l1p2 - p1l1) - l2k2*l2p2*p1l1 +
			        l2k2*p1k1*p1l1 - l2k1*p1k2*p1l1 - l1p2*l2k2*p1l2 + l2k2*p1l1*p1l2 +
			        l1k2*l2k1*p1p2 - l1k1*l2k2*p1p2 - l1l2*l2k2*p1p2 +
			        3*l2k2*m2*p1p2 + k1k2*
			         (l1k1*M2 + l1l2*M2 - 3*m2*M2 + k1p2*(l1p2 - p1l1) -
			           l2p2*p1l1 + p1k1*p1l1 + l1p2*(l2p2 - p1k1 - p1l2) + p1l1*p1l2 -
			           l1k1*p1p2 - l1l2*p1p2 + 3*m2*p1p2) - l1p2*l2k1*p2k2 +
			        l2k1*p1l1*p2k2) + 4*M2*
			      (l1k2*l2k1*M2 - l1k1*l2k2*M2 - l1l2*l2k2*M2 +
			        2*l2k2*m2*M2 + l1p2*l2k2*p1k1 - l1p2*l2k1*p1k2 +
			        k1p2*l2k2*p1l1 + l2k2*l2p2*p1l1 + l1p2*l2k2*p1l2 - l2k2*m2*p1p2 +
			        k1k2*(-(l1k1*M2) - l1l2*M2 + 2*m2*M2 + l1p2*p1k1 +
			           k1p2*p1l1 + l2p2*p1l1 + l1p2*p1l2 - m2*p1p2) - l2k1*p1l1*p2k2)*
			      pow(f1,2) + pow(f2,2)*(-3*k1p2*l1p2*l2k2*M2 -
			        3*l1p2*l2k2*l2p2*M2 + l1p2*l2k2*M2*p1k1 - l1p2*l2k1*M2*p1k2 +
			        k1p2*l2k2*M2*p1l1 + l2k2*l2p2*M2*p1l1 - 3*l2k2*M2*p1k1*p1l1 +
			        3*l2k1*M2*p1k2*p1l1 + l1p2*l2k2*M2*p1l2 -
			        3*l2k2*M2*p1l1*p1l2 + k1p2*l1p2*l2k2*p1p2 + l1p2*l2k2*l2p2*p1p2 -
			        2*l1k2*l2k1*M2*p1p2 + 2*l1k1*l2k2*M2*p1p2 +
			        2*l1l2*l2k2*M2*p1p2 - 6*l2k2*m2*M2*p1p2 + l1p2*l2k2*p1k1*p1p2 -
			        l1p2*l2k1*p1k2*p1p2 + k1p2*l2k2*p1l1*p1p2 + l2k2*l2p2*p1l1*p1p2 +
			        l2k2*p1k1*p1l1*p1p2 - l2k1*p1k2*p1l1*p1p2 + l1p2*l2k2*p1l2*p1p2 +
			        l2k2*p1l1*p1l2*p1p2 + 3*l1p2*l2k1*M2*p2k2 - l2k1*M2*p1l1*p2k2 -
			        l1p2*l2k1*p1p2*p2k2 - l2k1*p1l1*p1p2*p2k2 + l1k2*l2k1*M4 -
			        l1k1*l2k2*M4 - l1l2*l2k2*M4 +
			        5*l2k2*m2*M4 + l1k2*l2k1*pow(p1p2,2) -
			        l1k1*l2k2*pow(p1p2,2) - l1l2*l2k2*pow(p1p2,2) +
			        l2k2*m2*pow(p1p2,2) +
			        k1k2*(l2p2*M2*p1l1 - 3*M2*p1k1*p1l1 - 3*M2*p1l1*p1l2 +
			           2*l1k1*M2*p1p2 + 2*l1l2*M2*p1p2 - 6*m2*M2*p1p2 +
			           l2p2*p1l1*p1p2 + p1k1*p1l1*p1p2 + p1l1*p1l2*p1p2 +
			           k1p2*(l1p2*(-3*M2 + p1p2) + p1l1*(M2 + p1p2)) +
			           l1p2*(l2p2*(-3*M2 + p1p2) + (p1k1 + p1l2)*(M2 + p1p2)) -
			           l1k1*M4 - l1l2*M4 + 5*m2*M4 -
			           l1k1*pow(p1p2,2) - l1l2*pow(p1p2,2) + m2*pow(p1p2,2)))) +
			  8*eg1*pow(l1k2,-1)*pow(l1k1 + eg1*(-k1k2 + l1k2),-2)/M2*
			   (-4*f1*f2*M2*(l1k2*l1p2*l2p2 + l1k2*l1l2*M2 - l1k2*l2k1*M2 +
			        l1k1*l2k2*M2 - 3*l1k2*m2*M2 + l1k2*l2p2*p1k1 -
			        l1k1*l2p2*p1k2 - l1k2*l2p2*p1l1 - l1k2*l1p2*p1l2 - l1k2*p1k1*p1l2 +
			        l1k1*p1k2*p1l2 + l1k2*p1l1*p1l2 + k1p2*l1k2*(-l2p2 + p1l2) -
			        l1k2*l1l2*p1p2 + l1k2*l2k1*p1p2 - l1k1*l2k2*p1p2 +
			        3*l1k2*m2*p1p2 + k1k2*
			         (-(l1p2*l2p2) - l1l2*M2 + l2k1*M2 + 3*m2*M2 - l2p2*p1k1 +
			           l2p2*p1l1 + k1p2*(l2p2 - p1l2) + l1p2*p1l2 + p1k1*p1l2 -
			           p1l1*p1l2 + l1l2*p1p2 - l2k1*p1p2 - 3*m2*p1p2) +
			        l1k1*l2p2*p2k2 - l1k1*p1l2*p2k2) +
			     4*M2*(l1k2*(-(l1l2*M2) + l2k1*M2 + 2*m2*M2 - l2p2*p1k1 +
			           l2p2*p1l1 - k1p2*p1l2 + l1p2*p1l2 - m2*p1p2) +
			        k1k2*(l1l2*M2 - l2k1*M2 - 2*m2*M2 + l2p2*p1k1 -
			           l2p2*p1l1 + k1p2*p1l2 - l1p2*p1l2 + m2*p1p2) +
			        l1k1*(-(l2k2*M2) + l2p2*p1k2 + p1l2*p2k2))*pow(f1,2) +
			     pow(f2,2)*(-3*l1k2*l1p2*l2p2*M2 - l1k2*l2p2*M2*p1k1 +
			        l1k1*l2p2*M2*p1k2 + l1k2*l2p2*M2*p1l1 + l1k2*l1p2*M2*p1l2 +
			        3*l1k2*M2*p1k1*p1l2 - 3*l1k1*M2*p1k2*p1l2 -
			        3*l1k2*M2*p1l1*p1l2 + l1k2*l1p2*l2p2*p1p2 +
			        2*l1k2*l1l2*M2*p1p2 - 2*l1k2*l2k1*M2*p1p2 +
			        2*l1k1*l2k2*M2*p1p2 - 6*l1k2*m2*M2*p1p2 - l1k2*l2p2*p1k1*p1p2 +
			        l1k1*l2p2*p1k2*p1p2 + l1k2*l2p2*p1l1*p1p2 + l1k2*l1p2*p1l2*p1p2 -
			        l1k2*p1k1*p1l2*p1p2 + l1k1*p1k2*p1l2*p1p2 + l1k2*p1l1*p1l2*p1p2 +
			        k1p2*l1k2*(l2p2*(3*M2 - p1p2) - p1l2*(M2 + p1p2)) -
			        3*l1k1*l2p2*M2*p2k2 + l1k1*M2*p1l2*p2k2 + l1k1*l2p2*p1p2*p2k2 +
			        l1k1*p1l2*p1p2*p2k2 - l1k2*l1l2*M4 + l1k2*l2k1*M4 -
			        l1k1*l2k2*M4 + 5*l1k2*m2*M4 -
			        l1k2*l1l2*pow(p1p2,2) + l1k2*l2k1*pow(p1p2,2) -
			        l1k1*l2k2*pow(p1p2,2) + l1k2*m2*pow(p1p2,2) -
			        k1k2*(-(l2p2*M2*p1k1) + l2p2*M2*p1l1 + 3*M2*p1k1*p1l2 -
			           3*M2*p1l1*p1l2 + 2*l1l2*M2*p1p2 - 2*l2k1*M2*p1p2 -
			           6*m2*M2*p1p2 - l2p2*p1k1*p1p2 + l2p2*p1l1*p1p2 -
			           p1k1*p1l2*p1p2 + p1l1*p1l2*p1p2 +
			           k1p2*(l2p2*(3*M2 - p1p2) - p1l2*(M2 + p1p2)) +
			           l1p2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) -
			           l1l2*M4 + l2k1*M4 + 5*m2*M4 -
			           l1l2*pow(p1p2,2) + l2k1*pow(p1p2,2) + m2*pow(p1p2,2)))) -
			  4*eg1*pow(l1k1,-1)*pow(l2k1,-1)*pow(l2k2,-1)*
			   pow(l2k1 + eg1*(k1k2 + l2k2),-1)/M2*
			   (-4*M2*pow(f1,2)*(-2*l1k2*l1l2*l2k1*M2 - 4*eg1*l1k2*l1l2*l2k2*M2 +
			        2*eg1*l1k2*l2k1*l2k2*M2 + 2*l1k2*l1l2*m2*M2 -
			        3*l1k2*l2k1*m2*M2 - 4*l1l2*l2k2*m2*M2 +
			        2*l2k1*l2k2*m2*M2 - 2*l1k2*l1p2*l2k1*p1k1 -
			        4*k1p2*l1l2*l2k2*p1k1 + 2*l1l2*l1p2*l2k2*p1k1 -
			        2*l1k2*l2k1*l2p2*p1k1 + 2*k1p2*l1k2*m2*p1k1 - l1k2*l1p2*m2*p1k1 +
			        4*k1p2*l2k2*m2*p1k1 + l2k2*l2p2*m2*p1k1 -
			        2*eg1*l1k2*l1p2*l2k1*p1k2 + 2*l1l2*l1p2*l2k1*p1k2 -
			        2*eg1*k1p2*l1l2*l2k2*p1k2 + 4*eg1*l1l2*l1p2*l2k2*p1k2 -
			        2*eg1*l1p2*l2k1*l2k2*p1k2 + k1p2*l1l2*m2*p1k2 -
			        2*l1l2*l1p2*m2*p1k2 - 2*k1p2*l2k1*m2*p1k2 + l1p2*l2k1*m2*p1k2 +
			        l2k1*l2p2*m2*p1k2 - k1p2*m4*p1k2 - 2*k1p2*l1k2*l2k1*p1l1 +
			        2*k1p2*l1l2*l2k2*p1l1 - k1p2*l1k2*m2*p1l1 -
			        2*k1p2*l1k2*l2k1*p1l2 + k1p2*l2k2*m2*p1l2 + l1k2*l2k1*m2*p1p2 +
			        2*l1l2*l2k2*m2*p1p2 - 2*l2k1*l2k2*m2*p1p2 -
			        2*eg1*l1l2*l2k2*p1k1*p2k2 + l1l2*m2*p1k1*p2k2 -
			        2*l2k1*m2*p1k1*p2k2 - m4*p1k1*p2k2 - 2*eg1*l1k2*l2k1*p1l1*p2k2 +
			        2*l1l2*l2k1*p1l1*p2k2 + 4*eg1*l1l2*l2k2*p1l1*p2k2 -
			        2*eg1*l2k1*l2k2*p1l1*p2k2 - 2*l1l2*m2*p1l1*p2k2 +
			        l2k1*m2*p1l1*p2k2 + l2k1*m2*p1l2*p2k2 - 2*l2k2*M2*pow(l1k1,2) +
			        2*eg1*l2k1*M2*pow(l1k2,2) +
			        k1k2*(-7*l1l2*m2*M2 + 2*l2k1*m2*M2 - m4*M2 +
			           2*l1l2*l2p2*p1k1 - l1p2*m2*p1k1 - 2*l2p2*m2*p1k1 +
			           4*l1p2*l2k1*p1l1 - 2*l1l2*l2p2*p1l1 + 2*l2k1*l2p2*p1l1 -
			           k1p2*m2*p1l1 + 2*l1p2*m2*p1l1 + 2*k1p2*l1l2*p1l2 -
			           2*l1l2*l1p2*p1l2 + 2*l1p2*l2k1*p1l2 - 2*k1p2*m2*p1l2 -
			           2*l2p2*m2*p1l2 + 2*l1k1*
			            (l1l2*M2 - l2p2*p1l1 - l1p2*p1l2 - 2*l2p2*p1l2 -
			              2*m2*(M2 - p1p2)) - 2*l1l2*l2k1*p1p2 + 3*l1l2*m2*p1p2 +
			           m4*p1p2 + 2*eg1*l1l2*
			            (-(l1k2*M2) + l2k2*M2 + l1p2*p1k2 + p1l1*p2k2) +
			           2*M2*pow(l1l2,2)) - 2*l1p2*p1k2*pow(l2k1,2) +
			        2*l1k2*p1p2*pow(l2k1,2) - 2*p1l1*p2k2*pow(l2k1,2) +
			        l1k1*(-2*l1l2*l2k2*M2 - l2k2*m2*M2 +
			           2*l1k2*(l2k1 - eg1*l2k2 + m2)*M2 + 2*l1p2*l2k2*p1k1 +
			           2*l2k2*l2p2*p1k1 + 2*eg1*l1p2*l2k2*p1k2 + 2*l2k1*l2p2*p1k2 +
			           2*eg1*l2k2*l2p2*p1k2 - k1p2*m2*p1k2 - l1p2*m2*p1k2 -
			           l2p2*m2*p1k2 + 2*k1p2*l2k2*p1l1 + 2*k1p2*l2k2*p1l2 -
			           2*l2k1*l2k2*p1p2 + l2k2*m2*p1p2 - m2*p1k1*p2k2 +
			           2*eg1*l2k2*p1l1*p2k2 - m2*p1l1*p2k2 + 2*l2k1*p1l2*p2k2 +
			           2*eg1*l2k2*p1l2*p2k2 - m2*p1l2*p2k2 - 2*eg1*M2*pow(l2k2,2))) +
			     4*f1*f2*M2*(-2*l1k1*l1k2*l2k1*M2 + 2*l1k2*l1l2*l2k1*M2 +
			        2*eg1*l1k1*l1k2*l2k2*M2 + 2*l1k1*l1l2*l2k2*M2 +
			        4*eg1*l1k2*l1l2*l2k2*M2 - 2*l1k1*l2k1*l2k2*M2 -
			        2*eg1*l1k2*l2k1*l2k2*M2 - 2*l1k1*l1k2*m2*M2 -
			        2*l1k2*l1l2*m2*M2 + 4*l1k2*l2k1*m2*M2 +
			        2*l1k1*l2k2*m2*M2 + 6*l1l2*l2k2*m2*M2 -
			        4*l2k1*l2k2*m2*M2 + 2*l1k2*l1p2*l2k1*p1k1 -
			        2*l1k1*l1p2*l2k2*p1k1 - 2*l1l2*l1p2*l2k2*p1k1 +
			        2*l1k2*l2k1*l2p2*p1k1 - 2*l1k1*l2k2*l2p2*p1k1 +
			        l1k2*l1p2*m2*p1k1 - l2k2*l2p2*m2*p1k1 +
			        2*eg1*l1k2*l1p2*l2k1*p1k2 - 2*l1l2*l1p2*l2k1*p1k2 -
			        2*eg1*l1k1*l1p2*l2k2*p1k2 - 4*eg1*l1l2*l1p2*l2k2*p1k2 +
			        2*eg1*l1p2*l2k1*l2k2*p1k2 - 2*l1k1*l2k1*l2p2*p1k2 -
			        2*eg1*l1k1*l2k2*l2p2*p1k2 + l1k1*l1p2*m2*p1k2 +
			        2*l1l2*l1p2*m2*p1k2 - l1p2*l2k1*m2*p1k2 + l1k1*l2p2*m2*p1k2 -
			        l2k1*l2p2*m2*p1k2 - 2*eg1*l1l2*l2k2*p1k1*p1k2 -
			        l1k1*m2*p1k1*p1k2 + l1l2*m2*p1k1*p1k2 - 2*l2k1*m2*p1k1*p1k2 -
			        m4*p1k1*p1k2 - 2*l1k2*l2k1*p1k1*p1l1 + 2*l1k1*l2k2*p1k1*p1l1 +
			        2*l1l2*l2k2*p1k1*p1l1 - l1k2*m2*p1k1*p1l1 -
			        2*eg1*l1k2*l2k1*p1k2*p1l1 + 2*l1l2*l2k1*p1k2*p1l1 +
			        2*eg1*l1k1*l2k2*p1k2*p1l1 + 4*eg1*l1l2*l2k2*p1k2*p1l1 -
			        2*eg1*l2k1*l2k2*p1k2*p1l1 - l1k1*m2*p1k2*p1l1 -
			        2*l1l2*m2*p1k2*p1l1 + l2k1*m2*p1k2*p1l1 - 2*l1k2*l2k1*p1k1*p1l2 +
			        2*l1k1*l2k2*p1k1*p1l2 + l2k2*m2*p1k1*p1l2 +
			        2*l1k1*l2k1*p1k2*p1l2 + 2*eg1*l1k1*l2k2*p1k2*p1l2 -
			        l1k1*m2*p1k2*p1l2 + l2k1*m2*p1k2*p1l2 + 2*l1k1*l1k2*l2k1*p1p2 -
			        2*l1k2*l1l2*l2k1*p1p2 - 2*eg1*l1k1*l1k2*l2k2*p1p2 -
			        2*l1k1*l1l2*l2k2*p1p2 - 4*eg1*l1k2*l1l2*l2k2*p1p2 +
			        2*l1k1*l2k1*l2k2*p1p2 + 2*eg1*l1k2*l2k1*l2k2*p1p2 +
			        2*l1k1*l1k2*m2*p1p2 + 2*l1k2*l1l2*m2*p1p2 -
			        4*l1k2*l2k1*m2*p1p2 - 2*l1k1*l2k2*m2*p1p2 -
			        6*l1l2*l2k2*m2*p1p2 + 4*l2k1*l2k2*m2*p1p2 -
			        2*eg1*l1k2*l1p2*l2k1*p2k2 + 2*l1l2*l1p2*l2k1*p2k2 +
			        2*eg1*l1k1*l1p2*l2k2*p2k2 + 4*eg1*l1l2*l1p2*l2k2*p2k2 -
			        2*eg1*l1p2*l2k1*l2k2*p2k2 + 2*l1k1*l2k1*l2p2*p2k2 +
			        2*eg1*l1k1*l2k2*l2p2*p2k2 - l1k1*l1p2*m2*p2k2 -
			        2*l1l2*l1p2*m2*p2k2 + l1p2*l2k1*m2*p2k2 - l1k1*l2p2*m2*p2k2 +
			        l2k1*l2p2*m2*p2k2 + 2*eg1*l1l2*l2k2*p1k1*p2k2 +
			        l1k1*m2*p1k1*p2k2 - l1l2*m2*p1k1*p2k2 + 2*l2k1*m2*p1k1*p2k2 +
			        m4*p1k1*p2k2 + 2*eg1*l1k2*l2k1*p1l1*p2k2 - 2*l1l2*l2k1*p1l1*p2k2 -
			        2*eg1*l1k1*l2k2*p1l1*p2k2 - 4*eg1*l1l2*l2k2*p1l1*p2k2 +
			        2*eg1*l2k1*l2k2*p1l1*p2k2 + l1k1*m2*p1l1*p2k2 +
			        2*l1l2*m2*p1l1*p2k2 - l2k1*m2*p1l1*p2k2 - 2*l1k1*l2k1*p1l2*p2k2 -
			        2*eg1*l1k1*l2k2*p1l2*p2k2 + l1k1*m2*p1l2*p2k2 -
			        l2k1*m2*p1l2*p2k2 + k1p2*
			         (2*l1k1*l1p2*l2k2 + 2*l1l2*l1p2*l2k2 + 2*k1k2*l1l2*l2p2 +
			           2*l1k1*l2k2*l2p2 - k1k2*l1p2*m2 - 2*k1k2*l2p2*m2 +
			           l2k2*l2p2*m2 + 4*l1l2*l2k2*p1k1 - 4*l2k2*m2*p1k1 +
			           2*eg1*l1l2*l2k2*p1k2 + l1k1*m2*p1k2 - l1l2*m2*p1k2 +
			           2*l2k1*m2*p1k2 + m4*p1k2 - 2*l1k1*l2k2*p1l1 -
			           2*l1l2*l2k2*p1l1 + k1k2*m2*p1l1 - 2*k1k2*l1l2*p1l2 -
			           2*l1k1*l2k2*p1l2 + 2*k1k2*m2*p1l2 - l2k2*m2*p1l2 +
			           l1k2*(-(l1p2*(2*l2k1 + m2)) + m2*(-2*p1k1 + p1l1) +
			              2*l2k1*(-l2p2 + p1l1 + p1l2)) - 2*eg1*l1l2*l2k2*p2k2 -
			           l1k1*m2*p2k2 + l1l2*m2*p2k2 - 2*l2k1*m2*p2k2 - m4*p2k2) +
			        (-2*l1l2*l2k2 + (l1k2 + 2*l2k2)*m2)*pow(k1p2,2) +
			        2*l2k2*M2*pow(l1k1,2) - 2*l2k2*p1p2*pow(l1k1,2) -
			        2*eg1*l2k1*M2*pow(l1k2,2) + 2*eg1*l2k1*p1p2*pow(l1k2,2) +
			        2*l1k2*M2*pow(l2k1,2) + 2*l1p2*p1k2*pow(l2k1,2) -
			        2*p1k2*p1l1*pow(l2k1,2) - 2*l1k2*p1p2*pow(l2k1,2) -
			        2*l1p2*p2k2*pow(l2k1,2) + 2*p1l1*p2k2*pow(l2k1,2) +
			        2*eg1*l1k1*M2*pow(l2k2,2) - 2*eg1*l1k1*p1p2*pow(l2k2,2) -
			        2*l1l2*l2k2*pow(p1k1,2) + l1k2*m2*pow(p1k1,2) +
			        2*l2k2*m2*pow(p1k1,2) +
			        k1k2*(2*eg1*l1k2*l1l2*M2 - 2*l1l2*l2k1*M2 -
			           2*eg1*l1l2*l2k2*M2 + 10*l1l2*m2*M2 - 2*l2k1*m2*M2 +
			           2*m4*M2 - 2*l1l2*l2p2*p1k1 + 2*l2p2*m2*p1k1 +
			           2*l1l2*l2p2*p1l1 - 2*l2k1*l2p2*p1l1 - m2*p1k1*p1l1 +
			           2*eg1*l1l2*p1k2*p1l1 + 2*l2p2*m2*p1l2 + 2*l1l2*p1k1*p1l2 -
			           2*m2*p1k1*p1l2 - 2*l1l2*p1l1*p1l2 + 2*l2k1*p1l1*p1l2 -
			           2*eg1*l1k2*l1l2*p1p2 + 2*l1l2*l2k1*p1p2 + 2*eg1*l1l2*l2k2*p1p2 -
			           10*l1l2*m2*p1p2 + 2*l2k1*m2*p1p2 - 2*m4*p1p2 -
			           2*eg1*l1l2*p1l1*p2k2 +
			           l1p2*(2*l2k1*l2p2 + m2*p1k1 - 4*l2k1*p1l1 - 2*m2*p1l1 -
			              2*l1k1*(l2p2 - p1l2) - 2*l2k1*p1l2 -
			              2*l1l2*(l2p2 + eg1*p1k2 - p1l2 - eg1*p2k2)) -
			           2*M2*pow(l1l2,2) + 2*p1p2*pow(l1l2,2) +
			           (2*l2k1 + m2)*pow(l1p2,2) - m2*pow(l2p2,2) +
			           2*l2k1*pow(p1l1,2) + m2*pow(p1l1,2) - m2*pow(p1l2,2) -
			           2*l1k1*(-4*m2*M2 + p1l1*p1l2 - l2p2*(p1l1 + 2*p1l2) +
			              l1l2*(M2 - p1p2) + 4*m2*p1p2 + pow(l2p2,2) + pow(p1l2,2)))) \
			+ pow(f2,2)*(2*l1k2*l1p2*l2k1*M2*p1k1 - 2*l1k1*l1p2*l2k2*M2*p1k1 -
			        2*l1l2*l1p2*l2k2*M2*p1k1 + 2*l1k2*l2k1*l2p2*M2*p1k1 -
			        2*l1k1*l2k2*l2p2*M2*p1k1 + l1k2*l1p2*m2*M2*p1k1 -
			        l2k2*l2p2*m2*M2*p1k1 + 2*eg1*l1k2*l1p2*l2k1*M2*p1k2 -
			        2*l1l2*l1p2*l2k1*M2*p1k2 - 2*eg1*l1k1*l1p2*l2k2*M2*p1k2 -
			        4*eg1*l1l2*l1p2*l2k2*M2*p1k2 + 2*eg1*l1p2*l2k1*l2k2*M2*p1k2 -
			        2*l1k1*l2k1*l2p2*M2*p1k2 - 2*eg1*l1k1*l2k2*l2p2*M2*p1k2 +
			        l1k1*l1p2*m2*M2*p1k2 + 2*l1l2*l1p2*m2*M2*p1k2 -
			        l1p2*l2k1*m2*M2*p1k2 + l1k1*l2p2*m2*M2*p1k2 -
			        l2k1*l2p2*m2*M2*p1k2 - 6*eg1*l1l2*l2k2*M2*p1k1*p1k2 -
			        3*l1k1*m2*M2*p1k1*p1k2 + 3*l1l2*m2*M2*p1k1*p1k2 -
			        6*l2k1*m2*M2*p1k1*p1k2 - 3*m4*M2*p1k1*p1k2 -
			        6*l1k2*l2k1*M2*p1k1*p1l1 + 6*l1k1*l2k2*M2*p1k1*p1l1 +
			        6*l1l2*l2k2*M2*p1k1*p1l1 - 3*l1k2*m2*M2*p1k1*p1l1 -
			        6*eg1*l1k2*l2k1*M2*p1k2*p1l1 + 6*l1l2*l2k1*M2*p1k2*p1l1 +
			        6*eg1*l1k1*l2k2*M2*p1k2*p1l1 + 12*eg1*l1l2*l2k2*M2*p1k2*p1l1 -
			        6*eg1*l2k1*l2k2*M2*p1k2*p1l1 - 3*l1k1*m2*M2*p1k2*p1l1 -
			        6*l1l2*m2*M2*p1k2*p1l1 + 3*l2k1*m2*M2*p1k2*p1l1 -
			        6*l1k2*l2k1*M2*p1k1*p1l2 + 6*l1k1*l2k2*M2*p1k1*p1l2 +
			        3*l2k2*m2*M2*p1k1*p1l2 + 6*l1k1*l2k1*M2*p1k2*p1l2 +
			        6*eg1*l1k1*l2k2*M2*p1k2*p1l2 - 3*l1k1*m2*M2*p1k2*p1l2 +
			        3*l2k1*m2*M2*p1k2*p1l2 + 4*l1k1*l1k2*l2k1*M2*p1p2 -
			        4*l1k2*l1l2*l2k1*M2*p1p2 - 4*eg1*l1k1*l1k2*l2k2*M2*p1p2 -
			        4*l1k1*l1l2*l2k2*M2*p1p2 - 8*eg1*l1k2*l1l2*l2k2*M2*p1p2 +
			        4*l1k1*l2k1*l2k2*M2*p1p2 + 4*eg1*l1k2*l2k1*l2k2*M2*p1p2 +
			        4*l1k1*l1k2*m2*M2*p1p2 + 4*l1k2*l1l2*m2*M2*p1p2 -
			        8*l1k2*l2k1*m2*M2*p1p2 - 4*l1k1*l2k2*m2*M2*p1p2 -
			        12*l1l2*l2k2*m2*M2*p1p2 + 8*l2k1*l2k2*m2*M2*p1p2 +
			        2*l1k2*l1p2*l2k1*p1k1*p1p2 - 2*l1k1*l1p2*l2k2*p1k1*p1p2 -
			        2*l1l2*l1p2*l2k2*p1k1*p1p2 + 2*l1k2*l2k1*l2p2*p1k1*p1p2 -
			        2*l1k1*l2k2*l2p2*p1k1*p1p2 + l1k2*l1p2*m2*p1k1*p1p2 -
			        l2k2*l2p2*m2*p1k1*p1p2 + 2*eg1*l1k2*l1p2*l2k1*p1k2*p1p2 -
			        2*l1l2*l1p2*l2k1*p1k2*p1p2 - 2*eg1*l1k1*l1p2*l2k2*p1k2*p1p2 -
			        4*eg1*l1l2*l1p2*l2k2*p1k2*p1p2 + 2*eg1*l1p2*l2k1*l2k2*p1k2*p1p2 -
			        2*l1k1*l2k1*l2p2*p1k2*p1p2 - 2*eg1*l1k1*l2k2*l2p2*p1k2*p1p2 +
			        l1k1*l1p2*m2*p1k2*p1p2 + 2*l1l2*l1p2*m2*p1k2*p1p2 -
			        l1p2*l2k1*m2*p1k2*p1p2 + l1k1*l2p2*m2*p1k2*p1p2 -
			        l2k1*l2p2*m2*p1k2*p1p2 + 2*eg1*l1l2*l2k2*p1k1*p1k2*p1p2 +
			        l1k1*m2*p1k1*p1k2*p1p2 - l1l2*m2*p1k1*p1k2*p1p2 +
			        2*l2k1*m2*p1k1*p1k2*p1p2 + m4*p1k1*p1k2*p1p2 +
			        2*l1k2*l2k1*p1k1*p1l1*p1p2 - 2*l1k1*l2k2*p1k1*p1l1*p1p2 -
			        2*l1l2*l2k2*p1k1*p1l1*p1p2 + l1k2*m2*p1k1*p1l1*p1p2 +
			        2*eg1*l1k2*l2k1*p1k2*p1l1*p1p2 - 2*l1l2*l2k1*p1k2*p1l1*p1p2 -
			        2*eg1*l1k1*l2k2*p1k2*p1l1*p1p2 - 4*eg1*l1l2*l2k2*p1k2*p1l1*p1p2 +
			        2*eg1*l2k1*l2k2*p1k2*p1l1*p1p2 + l1k1*m2*p1k2*p1l1*p1p2 +
			        2*l1l2*m2*p1k2*p1l1*p1p2 - l2k1*m2*p1k2*p1l1*p1p2 +
			        2*l1k2*l2k1*p1k1*p1l2*p1p2 - 2*l1k1*l2k2*p1k1*p1l2*p1p2 -
			        l2k2*m2*p1k1*p1l2*p1p2 - 2*l1k1*l2k1*p1k2*p1l2*p1p2 -
			        2*eg1*l1k1*l2k2*p1k2*p1l2*p1p2 + l1k1*m2*p1k2*p1l2*p1p2 -
			        l2k1*m2*p1k2*p1l2*p1p2 - 6*eg1*l1k2*l1p2*l2k1*M2*p2k2 +
			        6*l1l2*l1p2*l2k1*M2*p2k2 + 6*eg1*l1k1*l1p2*l2k2*M2*p2k2 +
			        12*eg1*l1l2*l1p2*l2k2*M2*p2k2 - 6*eg1*l1p2*l2k1*l2k2*M2*p2k2 +
			        6*l1k1*l2k1*l2p2*M2*p2k2 + 6*eg1*l1k1*l2k2*l2p2*M2*p2k2 -
			        3*l1k1*l1p2*m2*M2*p2k2 - 6*l1l2*l1p2*m2*M2*p2k2 +
			        3*l1p2*l2k1*m2*M2*p2k2 - 3*l1k1*l2p2*m2*M2*p2k2 +
			        3*l2k1*l2p2*m2*M2*p2k2 + 2*eg1*l1l2*l2k2*M2*p1k1*p2k2 +
			        l1k1*m2*M2*p1k1*p2k2 - l1l2*m2*M2*p1k1*p2k2 +
			        2*l2k1*m2*M2*p1k1*p2k2 + m4*M2*p1k1*p2k2 +
			        2*eg1*l1k2*l2k1*M2*p1l1*p2k2 - 2*l1l2*l2k1*M2*p1l1*p2k2 -
			        2*eg1*l1k1*l2k2*M2*p1l1*p2k2 - 4*eg1*l1l2*l2k2*M2*p1l1*p2k2 +
			        2*eg1*l2k1*l2k2*M2*p1l1*p2k2 + l1k1*m2*M2*p1l1*p2k2 +
			        2*l1l2*m2*M2*p1l1*p2k2 - l2k1*m2*M2*p1l1*p2k2 -
			        2*l1k1*l2k1*M2*p1l2*p2k2 - 2*eg1*l1k1*l2k2*M2*p1l2*p2k2 +
			        l1k1*m2*M2*p1l2*p2k2 - l2k1*m2*M2*p1l2*p2k2 +
			        2*eg1*l1k2*l1p2*l2k1*p1p2*p2k2 - 2*l1l2*l1p2*l2k1*p1p2*p2k2 -
			        2*eg1*l1k1*l1p2*l2k2*p1p2*p2k2 - 4*eg1*l1l2*l1p2*l2k2*p1p2*p2k2 +
			        2*eg1*l1p2*l2k1*l2k2*p1p2*p2k2 - 2*l1k1*l2k1*l2p2*p1p2*p2k2 -
			        2*eg1*l1k1*l2k2*l2p2*p1p2*p2k2 + l1k1*l1p2*m2*p1p2*p2k2 +
			        2*l1l2*l1p2*m2*p1p2*p2k2 - l1p2*l2k1*m2*p1p2*p2k2 +
			        l1k1*l2p2*m2*p1p2*p2k2 - l2k1*l2p2*m2*p1p2*p2k2 +
			        2*eg1*l1l2*l2k2*p1k1*p1p2*p2k2 + l1k1*m2*p1k1*p1p2*p2k2 -
			        l1l2*m2*p1k1*p1p2*p2k2 + 2*l2k1*m2*p1k1*p1p2*p2k2 +
			        m4*p1k1*p1p2*p2k2 + 2*eg1*l1k2*l2k1*p1l1*p1p2*p2k2 -
			        2*l1l2*l2k1*p1l1*p1p2*p2k2 - 2*eg1*l1k1*l2k2*p1l1*p1p2*p2k2 -
			        4*eg1*l1l2*l2k2*p1l1*p1p2*p2k2 + 2*eg1*l2k1*l2k2*p1l1*p1p2*p2k2 +
			        l1k1*m2*p1l1*p1p2*p2k2 + 2*l1l2*m2*p1l1*p1p2*p2k2 -
			        l2k1*m2*p1l1*p1p2*p2k2 - 2*l1k1*l2k1*p1l2*p1p2*p2k2 -
			        2*eg1*l1k1*l2k2*p1l2*p1p2*p2k2 + l1k1*m2*p1l2*p1p2*p2k2 -
			        l2k1*m2*p1l2*p1p2*p2k2 +
			        k1p2*(6*l1k1*l1p2*l2k2*M2 + 6*l1l2*l1p2*l2k2*M2 +
			           6*k1k2*l1l2*l2p2*M2 + 6*l1k1*l2k2*l2p2*M2 -
			           3*k1k2*l1p2*m2*M2 - 6*k1k2*l2p2*m2*M2 +
			           3*l2k2*l2p2*m2*M2 + 4*l1l2*l2k2*M2*p1k1 -
			           4*l2k2*m2*M2*p1k1 + 2*eg1*l1l2*l2k2*M2*p1k2 +
			           l1k1*m2*M2*p1k2 - l1l2*m2*M2*p1k2 + 2*l2k1*m2*M2*p1k2 +
			           m4*M2*p1k2 - 2*l1k1*l2k2*M2*p1l1 - 2*l1l2*l2k2*M2*p1l1 +
			           k1k2*m2*M2*p1l1 - 2*k1k2*l1l2*M2*p1l2 -
			           2*l1k1*l2k2*M2*p1l2 + 2*k1k2*m2*M2*p1l2 -
			           l2k2*m2*M2*p1l2 - 2*l1k1*l1p2*l2k2*p1p2 -
			           2*l1l2*l1p2*l2k2*p1p2 - 2*k1k2*l1l2*l2p2*p1p2 -
			           2*l1k1*l2k2*l2p2*p1p2 + k1k2*l1p2*m2*p1p2 +
			           2*k1k2*l2p2*m2*p1p2 - l2k2*l2p2*m2*p1p2 +
			           4*l1l2*l2k2*p1k1*p1p2 - 4*l2k2*m2*p1k1*p1p2 +
			           2*eg1*l1l2*l2k2*p1k2*p1p2 + l1k1*m2*p1k2*p1p2 -
			           l1l2*m2*p1k2*p1p2 + 2*l2k1*m2*p1k2*p1p2 + m4*p1k2*p1p2 -
			           2*l1k1*l2k2*p1l1*p1p2 - 2*l1l2*l2k2*p1l1*p1p2 +
			           k1k2*m2*p1l1*p1p2 - 2*k1k2*l1l2*p1l2*p1p2 -
			           2*l1k1*l2k2*p1l2*p1p2 + 2*k1k2*m2*p1l2*p1p2 -
			           l2k2*m2*p1l2*p1p2 +
			           l1k2*(-(l1p2*(2*l2k1 + m2)*(3*M2 - p1p2)) -
			              m2*(2*p1k1 - p1l1)*(M2 + p1p2) +
			              2*l2k1*(l2p2*(-3*M2 + p1p2) + (p1l1 + p1l2)*(M2 + p1p2))) \
			- 6*eg1*l1l2*l2k2*M2*p2k2 - 3*l1k1*m2*M2*p2k2 + 3*l1l2*m2*M2*p2k2 -
			           6*l2k1*m2*M2*p2k2 - 3*m4*M2*p2k2 +
			           2*eg1*l1l2*l2k2*p1p2*p2k2 + l1k1*m2*p1p2*p2k2 -
			           l1l2*m2*p1p2*p2k2 + 2*l2k1*m2*p1p2*p2k2 + m4*p1p2*p2k2) -
			        (2*l1l2*l2k2 - (l1k2 + 2*l2k2)*m2)*(3*M2 - p1p2)*pow(k1p2,2) -
			        4*l2k2*M2*p1p2*pow(l1k1,2) + 4*eg1*l2k1*M2*p1p2*pow(l1k2,2) +
			        2*l1p2*M2*p1k2*pow(l2k1,2) - 6*M2*p1k2*p1l1*pow(l2k1,2) -
			        4*l1k2*M2*p1p2*pow(l2k1,2) + 2*l1p2*p1k2*p1p2*pow(l2k1,2) +
			        2*p1k2*p1l1*p1p2*pow(l2k1,2) - 6*l1p2*M2*p2k2*pow(l2k1,2) +
			        2*M2*p1l1*p2k2*pow(l2k1,2) + 2*l1p2*p1p2*p2k2*pow(l2k1,2) +
			        2*p1l1*p1p2*p2k2*pow(l2k1,2) - 4*eg1*l1k1*M2*p1p2*pow(l2k2,2) -
			        2*l1k1*l1k2*l2k1*M4 + 2*l1k2*l1l2*l2k1*M4 +
			        2*eg1*l1k1*l1k2*l2k2*M4 + 2*l1k1*l1l2*l2k2*M4 +
			        4*eg1*l1k2*l1l2*l2k2*M4 - 6*l1k1*l2k1*l2k2*M4 -
			        2*eg1*l1k2*l2k1*l2k2*M4 - 2*l1k1*l1k2*m2*M4 -
			        2*l1k2*l1l2*m2*M4 + 6*l1k2*l2k1*m2*M4 +
			        4*l1k1*l2k2*m2*M4 + 10*l1l2*l2k2*m2*M4 -
			        8*l2k1*l2k2*m2*M4 + 2*l2k2*pow(l1k1,2)*M4 -
			        2*eg1*l2k1*pow(l1k2,2)*M4 + 6*l1k2*pow(l2k1,2)*M4 +
			        2*eg1*l1k1*pow(l2k2,2)*M4 - 6*l1l2*l2k2*M2*pow(p1k1,2) +
			        3*l1k2*m2*M2*pow(p1k1,2) + 6*l2k2*m2*M2*pow(p1k1,2) +
			        2*l1l2*l2k2*p1p2*pow(p1k1,2) - l1k2*m2*p1p2*pow(p1k1,2) -
			        2*l2k2*m2*p1p2*pow(p1k1,2) - 2*l1k1*l1k2*l2k1*pow(p1p2,2) +
			        2*l1k2*l1l2*l2k1*pow(p1p2,2) + 2*eg1*l1k1*l1k2*l2k2*pow(p1p2,2) +
			        2*l1k1*l1l2*l2k2*pow(p1p2,2) + 4*eg1*l1k2*l1l2*l2k2*pow(p1p2,2) +
			        2*l1k1*l2k1*l2k2*pow(p1p2,2) - 2*eg1*l1k2*l2k1*l2k2*pow(p1p2,2) -
			        2*l1k1*l1k2*m2*pow(p1p2,2) - 2*l1k2*l1l2*m2*pow(p1p2,2) +
			        2*l1k2*l2k1*m2*pow(p1p2,2) + 2*l1l2*l2k2*m2*pow(p1p2,2) +
			        2*l2k2*pow(l1k1,2)*pow(p1p2,2) -
			        2*eg1*l2k1*pow(l1k2,2)*pow(p1p2,2) -
			        2*l1k2*pow(l2k1,2)*pow(p1p2,2) +
			        2*eg1*l1k1*pow(l2k2,2)*pow(p1p2,2) +
			        k1k2*(-2*l1l2*l2p2*M2*p1k1 + 2*l2p2*m2*M2*p1k1 +
			           2*l1l2*l2p2*M2*p1l1 - 2*l2k1*l2p2*M2*p1l1 -
			           3*m2*M2*p1k1*p1l1 + 6*eg1*l1l2*M2*p1k2*p1l1 +
			           2*l2p2*m2*M2*p1l2 + 6*l1l2*M2*p1k1*p1l2 -
			           6*m2*M2*p1k1*p1l2 - 6*l1l2*M2*p1l1*p1l2 +
			           6*l2k1*M2*p1l1*p1l2 - 4*eg1*l1k2*l1l2*M2*p1p2 +
			           4*l1l2*l2k1*M2*p1p2 + 4*eg1*l1l2*l2k2*M2*p1p2 -
			           20*l1l2*m2*M2*p1p2 + 4*l2k1*m2*M2*p1p2 - 4*m4*M2*p1p2 -
			           2*l1l2*l2p2*p1k1*p1p2 + 2*l2p2*m2*p1k1*p1p2 +
			           2*l1l2*l2p2*p1l1*p1p2 - 2*l2k1*l2p2*p1l1*p1p2 +
			           m2*p1k1*p1l1*p1p2 - 2*eg1*l1l2*p1k2*p1l1*p1p2 +
			           2*l2p2*m2*p1l2*p1p2 - 2*l1l2*p1k1*p1l2*p1p2 +
			           2*m2*p1k1*p1l2*p1p2 + 2*l1l2*p1l1*p1l2*p1p2 -
			           2*l2k1*p1l1*p1l2*p1p2 - 2*eg1*l1l2*M2*p1l1*p2k2 -
			           2*eg1*l1l2*p1l1*p1p2*p2k2 +
			           l1p2*(6*l2k1*l2p2*M2 + m2*M2*p1k1 - 4*l2k1*M2*p1l1 -
			              2*m2*M2*p1l1 - 2*l2k1*M2*p1l2 - 2*l2k1*l2p2*p1p2 +
			              m2*p1k1*p1p2 - 4*l2k1*p1l1*p1p2 - 2*m2*p1l1*p1p2 -
			              2*l2k1*p1l2*p1p2 +
			              2*l1k1*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) -
			              2*l1l2*(l2p2*(3*M2 - p1p2) - p1l2*(M2 + p1p2) +
			                 eg1*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)))) +
			           4*M2*p1p2*pow(l1l2,2) +
			           (2*l2k1 + m2)*(3*M2 - p1p2)*pow(l1p2,2) -
			           3*m2*M2*pow(l2p2,2) + m2*p1p2*pow(l2p2,2) +
			           2*eg1*l1k2*l1l2*M4 - 6*l1l2*l2k1*M4 -
			           2*eg1*l1l2*l2k2*M4 + 16*l1l2*m2*M4 -
			           2*l2k1*m2*M4 + 4*m4*M4 -
			           2*pow(l1l2,2)*M4 + 6*l2k1*M2*pow(p1l1,2) +
			           3*m2*M2*pow(p1l1,2) - 2*l2k1*p1p2*pow(p1l1,2) -
			           m2*p1p2*pow(p1l1,2) - 3*m2*M2*pow(p1l2,2) +
			           m2*p1p2*pow(p1l2,2) +
			           2*l1k1*(-3*M2*p1l1*p1l2 - 8*m2*M2*p1p2 + p1l1*p1l2*p1p2 +
			              l2p2*(p1l1 + 2*p1l2)*(M2 + p1p2) +
			              (-3*M2 + p1p2)*pow(l2p2,2) + 8*m2*M4 -
			              3*M2*pow(p1l2,2) + p1p2*pow(p1l2,2) -
			              l1l2*pow(M2 - p1p2,2)) + 2*eg1*l1k2*l1l2*pow(p1p2,2) +
			           2*l1l2*l2k1*pow(p1p2,2) - 2*eg1*l1l2*l2k2*pow(p1p2,2) +
			           4*l1l2*m2*pow(p1p2,2) - 2*l2k1*m2*pow(p1p2,2) -
			           2*pow(l1l2,2)*pow(p1p2,2)))) -
			  4*eg1*pow(l1k1,-1)*pow(l1k2,-1)*pow(l1k1 + eg1*(-k1k2 + l1k2),-1)*
			   pow(l2k1,-1)/M2*(-4*M2*pow(f1,2)*
			      (2*l1k2*l1l2*l2k1*M2 + 2*l1k1*l1l2*l2k2*M2 +
			        2*l1k1*l2k1*l2k2*M2 - 2*l1k1*l1k2*m2*M2 -
			        4*l1k2*l1l2*m2*M2 + l1k2*l2k1*m2*M2 + 3*l1k1*l2k2*m2*M2 +
			        2*l1l2*l2k2*m2*M2 - 2*l2k1*l2k2*m2*M2 -
			        4*k1p2*l1k2*l1l2*p1k1 + 2*l1k2*l1p2*l2k1*p1k1 -
			        2*l1k1*l1p2*l2k2*p1k1 - 2*l1k2*l1l2*l2p2*p1k1 +
			        2*l1k2*l2k1*l2p2*p1k1 - 2*l1k1*l2k2*l2p2*p1k1 +
			        4*k1p2*l1k2*m2*p1k1 - l1k2*l1p2*m2*p1k1 + 2*k1p2*l2k2*m2*p1k1 +
			        l2k2*l2p2*m2*p1k1 + 2*l1k1*l1p2*l2k1*p1k2 -
			        2*l1k1*l1l2*l2p2*p1k2 - 2*k1p2*l1k1*m2*p1k2 - k1p2*l1l2*m2*p1k2 -
			        l1k1*l1p2*m2*p1k2 - k1p2*l2k1*m2*p1k2 + l1p2*l2k1*m2*p1k2 -
			        l1k1*l2p2*m2*p1k2 - 2*l1l2*l2p2*m2*p1k2 + l2k1*l2p2*m2*p1k2 +
			        k1p2*m4*p1k2 + 2*k1p2*l1k2*l2k1*p1l1 - 2*k1p2*l1k1*l2k2*p1l1 -
			        k1p2*l1k2*m2*p1l1 - 2*k1p2*l1k2*l1l2*p1l2 +
			        2*k1p2*l1k2*l2k1*p1l2 - 2*k1p2*l1k1*l2k2*p1l2 +
			        k1p2*l2k2*m2*p1l2 - 2*l1k1*l1k2*l2k1*p1p2 + 2*l1k1*l1k2*m2*p1p2 +
			        2*l1k2*l1l2*m2*p1p2 - l1k2*l2k1*m2*p1p2 - l1k1*l2k2*m2*p1p2 -
			        2*l1k1*m2*p1k1*p2k2 - l1l2*m2*p1k1*p2k2 - l2k1*m2*p1k1*p2k2 +
			        m4*p1k1*p2k2 + 2*l1k1*l2k1*p1l1*p2k2 - l1k1*m2*p1l1*p2k2 +
			        l2k1*m2*p1l1*p2k2 - 2*l1k1*l1l2*p1l2*p2k2 - l1k1*m2*p1l2*p2k2 -
			        2*l1l2*m2*p1l2*p2k2 + l2k1*m2*p1l2*p2k2 -
			        2*l2p2*p1k2*pow(l1k1,2) + 2*l2k2*p1p2*pow(l1k1,2) -
			        2*p1l2*p2k2*pow(l1k1,2) +
			        2*eg1*(l1k1*l2k2*(l2k2*M2 - l2p2*p1k2 - p1l2*p2k2) +
			           k1k2*l1l2*(l1k2*M2 - l2k2*M2 + l2p2*p1k2 + p1l2*p2k2) +
			           l1k2*(l1l2*(2*l2k2*M2 - k1p2*p1k2 - 2*l2p2*p1k2 -
			                 p1k1*p2k2 - 2*p1l2*p2k2) +
			              l1k1*(l2k2*M2 - l2p2*p1k2 - p1l2*p2k2) +
			              l2k1*(-(l2k2*M2) + l1p2*p1k2 + l2p2*p1k2 + p1l1*p2k2 +
			                 p1l2*p2k2)) - l2k1*M2*pow(l1k2,2)) -
			        k1k2*(4*l2k1*m2*M2 - m4*M2 + 2*l1p2*m2*p1k1 +
			           l2p2*m2*p1k1 + 4*l1p2*l2k1*p1l1 + 2*l2k1*l2p2*p1l1 +
			           2*k1p2*m2*p1l1 - 2*l1p2*m2*p1l1 + 2*l1p2*l2k1*p1l2 +
			           k1p2*m2*p1l2 + 2*l2p2*m2*p1l2 -
			           2*l1k1*(m2*M2 + l1p2*p1l2 + l2p2*(p1l1 + 2*p1l2)) -
			           4*l2k1*m2*p1p2 + m4*p1p2 -
			           l1l2*(2*l2k1*M2 + m2*(7*M2 - 3*p1p2) +
			              2*(k1p2*p1l1 + l2p2*p1l1 + l1p2*(p1k1 + p1l2) - l1k1*p1p2)) \
			+ 2*M2*pow(l1l2,2)) - 2*l1k2*M2*pow(l2k1,2)) +
			     4*f1*f2*M2*(-2*l1k1*l1k2*l2k1*M2 - 2*l1k2*l1l2*l2k1*M2 -
			        2*eg1*l1k1*l1k2*l2k2*M2 - 2*l1k1*l1l2*l2k2*M2 -
			        4*eg1*l1k2*l1l2*l2k2*M2 - 2*l1k1*l2k1*l2k2*M2 +
			        2*eg1*l1k2*l2k1*l2k2*M2 + 4*l1k1*l1k2*m2*M2 +
			        6*l1k2*l1l2*m2*M2 - 2*l1k2*l2k1*m2*M2 -
			        4*l1k1*l2k2*m2*M2 - 2*l1l2*l2k2*m2*M2 +
			        2*l2k1*l2k2*m2*M2 - 2*l1k2*l1p2*l2k1*p1k1 +
			        2*l1k1*l1p2*l2k2*p1k1 + 2*l1k2*l1l2*l2p2*p1k1 -
			        2*l1k2*l2k1*l2p2*p1k1 + 2*l1k1*l2k2*l2p2*p1k1 +
			        l1k2*l1p2*m2*p1k1 - l2k2*l2p2*m2*p1k1 - 2*l1k1*l1p2*l2k1*p1k2 -
			        2*eg1*l1k2*l1p2*l2k1*p1k2 + 2*eg1*l1k1*l1k2*l2p2*p1k2 +
			        2*l1k1*l1l2*l2p2*p1k2 + 4*eg1*l1k2*l1l2*l2p2*p1k2 -
			        2*eg1*l1k2*l2k1*l2p2*p1k2 + 2*eg1*l1k1*l2k2*l2p2*p1k2 +
			        l1k1*l1p2*m2*p1k2 - l1p2*l2k1*m2*p1k2 + l1k1*l2p2*m2*p1k2 +
			        2*l1l2*l2p2*m2*p1k2 - l2k1*l2p2*m2*p1k2 -
			        2*eg1*l1k2*l1l2*p1k1*p1k2 - 2*l1k1*m2*p1k1*p1k2 -
			        l1l2*m2*p1k1*p1k2 - l2k1*m2*p1k1*p1k2 + m4*p1k1*p1k2 +
			        2*l1k2*l2k1*p1k1*p1l1 - 2*l1k1*l2k2*p1k1*p1l1 -
			        l1k2*m2*p1k1*p1l1 + 2*l1k1*l2k1*p1k2*p1l1 +
			        2*eg1*l1k2*l2k1*p1k2*p1l1 - l1k1*m2*p1k2*p1l1 +
			        l2k1*m2*p1k2*p1l1 - 2*l1k2*l1l2*p1k1*p1l2 +
			        2*l1k2*l2k1*p1k1*p1l2 - 2*l1k1*l2k2*p1k1*p1l2 +
			        l2k2*m2*p1k1*p1l2 - 2*eg1*l1k1*l1k2*p1k2*p1l2 -
			        2*l1k1*l1l2*p1k2*p1l2 - 4*eg1*l1k2*l1l2*p1k2*p1l2 +
			        2*eg1*l1k2*l2k1*p1k2*p1l2 - 2*eg1*l1k1*l2k2*p1k2*p1l2 -
			        l1k1*m2*p1k2*p1l2 - 2*l1l2*m2*p1k2*p1l2 + l2k1*m2*p1k2*p1l2 +
			        2*l1k1*l1k2*l2k1*p1p2 + 2*l1k2*l1l2*l2k1*p1p2 +
			        2*eg1*l1k1*l1k2*l2k2*p1p2 + 2*l1k1*l1l2*l2k2*p1p2 +
			        4*eg1*l1k2*l1l2*l2k2*p1p2 + 2*l1k1*l2k1*l2k2*p1p2 -
			        2*eg1*l1k2*l2k1*l2k2*p1p2 - 4*l1k1*l1k2*m2*p1p2 -
			        6*l1k2*l1l2*m2*p1p2 + 2*l1k2*l2k1*m2*p1p2 +
			        4*l1k1*l2k2*m2*p1p2 + 2*l1l2*l2k2*m2*p1p2 -
			        2*l2k1*l2k2*m2*p1p2 + 2*l1k1*l1p2*l2k1*p2k2 +
			        2*eg1*l1k2*l1p2*l2k1*p2k2 - 2*eg1*l1k1*l1k2*l2p2*p2k2 -
			        2*l1k1*l1l2*l2p2*p2k2 - 4*eg1*l1k2*l1l2*l2p2*p2k2 +
			        2*eg1*l1k2*l2k1*l2p2*p2k2 - 2*eg1*l1k1*l2k2*l2p2*p2k2 -
			        l1k1*l1p2*m2*p2k2 + l1p2*l2k1*m2*p2k2 - l1k1*l2p2*m2*p2k2 -
			        2*l1l2*l2p2*m2*p2k2 + l2k1*l2p2*m2*p2k2 +
			        2*eg1*l1k2*l1l2*p1k1*p2k2 + 2*l1k1*m2*p1k1*p2k2 +
			        l1l2*m2*p1k1*p2k2 + l2k1*m2*p1k1*p2k2 - m4*p1k1*p2k2 -
			        2*l1k1*l2k1*p1l1*p2k2 - 2*eg1*l1k2*l2k1*p1l1*p2k2 +
			        l1k1*m2*p1l1*p2k2 - l2k1*m2*p1l1*p2k2 +
			        2*eg1*l1k1*l1k2*p1l2*p2k2 + 2*l1k1*l1l2*p1l2*p2k2 +
			        4*eg1*l1k2*l1l2*p1l2*p2k2 - 2*eg1*l1k2*l2k1*p1l2*p2k2 +
			        2*eg1*l1k1*l2k2*p1l2*p2k2 + l1k1*m2*p1l2*p2k2 +
			        2*l1l2*m2*p1l2*p2k2 - l2k1*m2*p1l2*p2k2 +
			        k1p2*(-2*l1k1*l1p2*l2k2 - 2*l1k1*l2k2*l2p2 + l2k2*l2p2*m2 -
			           2*l2k2*m2*p1k1 + 2*l1k1*m2*p1k2 + l1l2*m2*p1k2 +
			           l2k1*m2*p1k2 - m4*p1k2 + 2*l1k1*l2k2*p1l1 +
			           2*l1k1*l2k2*p1l2 - l2k2*m2*p1l2 +
			           k1k2*(2*l1l2*(l1p2 - p1l1) +
			              m2*(-2*l1p2 - l2p2 + 2*p1l1 + p1l2)) - 2*l1k1*m2*p2k2 -
			           l1l2*m2*p2k2 - l2k1*m2*p2k2 + m4*p2k2 +
			           l1k2*(2*l2k1*l2p2 + l1p2*(2*l2k1 - m2) - 4*m2*p1k1 -
			              2*l2k1*p1l1 + m2*p1l1 - 2*l2k1*p1l2 -
			              2*l1l2*(l2p2 - 2*p1k1 - eg1*p1k2 - p1l2 + eg1*p2k2))) +
			        (-2*l1k2*l1l2 + 2*l1k2*m2 + l2k2*m2)*pow(k1p2,2) +
			        2*l2k2*M2*pow(l1k1,2) + 2*l2p2*p1k2*pow(l1k1,2) -
			        2*p1k2*p1l2*pow(l1k1,2) - 2*l2k2*p1p2*pow(l1k1,2) -
			        2*l2p2*p2k2*pow(l1k1,2) + 2*p1l2*p2k2*pow(l1k1,2) +
			        2*eg1*l2k1*M2*pow(l1k2,2) - 2*eg1*l2k1*p1p2*pow(l1k2,2) +
			        2*l1k2*M2*pow(l2k1,2) - 2*l1k2*p1p2*pow(l2k1,2) -
			        2*eg1*l1k1*M2*pow(l2k2,2) + 2*eg1*l1k1*p1p2*pow(l2k2,2) -
			        2*l1k2*l1l2*pow(p1k1,2) + 2*l1k2*m2*pow(p1k1,2) +
			        l2k2*m2*pow(p1k1,2) +
			        k1k2*(-2*eg1*l1k2*l1l2*M2 - 2*l1l2*l2k1*M2 +
			           2*eg1*l1l2*l2k2*M2 - 10*l1l2*m2*M2 + 8*l2k1*m2*M2 -
			           2*m4*M2 + l2p2*m2*p1k1 - 2*eg1*l1l2*l2p2*p1k2 -
			           2*l1l2*l2p2*p1l1 + 2*l2k1*l2p2*p1l1 + 2*l1l2*p1k1*p1l1 -
			           2*m2*p1k1*p1l1 + 2*l2p2*m2*p1l2 - m2*p1k1*p1l2 +
			           2*eg1*l1l2*p1k2*p1l2 + 2*l1l2*p1l1*p1l2 - 2*l2k1*p1l1*p1l2 +
			           2*l1p2*(-(l2k1*l2p2) + m2*p1k1 + 2*l2k1*p1l1 - m2*p1l1 +
			              l1k1*(l2p2 - p1l2) + l1l2*(l2p2 - p1k1 - p1l2) + l2k1*p1l2) +
			           2*eg1*l1k2*l1l2*p1p2 + 2*l1l2*l2k1*p1p2 - 2*eg1*l1l2*l2k2*p1p2 +
			           10*l1l2*m2*p1p2 - 8*l2k1*m2*p1p2 + 2*m4*p1p2 +
			           2*eg1*l1l2*l2p2*p2k2 - 2*eg1*l1l2*p1l2*p2k2 +
			           2*M2*pow(l1l2,2) - 2*p1p2*pow(l1l2,2) +
			           (-2*l2k1 + m2)*pow(l1p2,2) - m2*pow(l2p2,2) -
			           2*l2k1*pow(p1l1,2) + m2*pow(p1l1,2) - m2*pow(p1l2,2) +
			           2*l1k1*(-(m2*M2) + p1l1*p1l2 - l2p2*(p1l1 + 2*p1l2) +
			              m2*p1p2 + l1l2*(-M2 + p1p2) + pow(l2p2,2) + pow(p1l2,2)))) \
			+ pow(f2,2)*(-2*l1k2*l1p2*l2k1*M2*p1k1 + 2*l1k1*l1p2*l2k2*M2*p1k1 +
			        2*l1k2*l1l2*l2p2*M2*p1k1 - 2*l1k2*l2k1*l2p2*M2*p1k1 +
			        2*l1k1*l2k2*l2p2*M2*p1k1 + l1k2*l1p2*m2*M2*p1k1 -
			        l2k2*l2p2*m2*M2*p1k1 - 2*l1k1*l1p2*l2k1*M2*p1k2 -
			        2*eg1*l1k2*l1p2*l2k1*M2*p1k2 + 2*eg1*l1k1*l1k2*l2p2*M2*p1k2 +
			        2*l1k1*l1l2*l2p2*M2*p1k2 + 4*eg1*l1k2*l1l2*l2p2*M2*p1k2 -
			        2*eg1*l1k2*l2k1*l2p2*M2*p1k2 + 2*eg1*l1k1*l2k2*l2p2*M2*p1k2 +
			        l1k1*l1p2*m2*M2*p1k2 - l1p2*l2k1*m2*M2*p1k2 +
			        l1k1*l2p2*m2*M2*p1k2 + 2*l1l2*l2p2*m2*M2*p1k2 -
			        l2k1*l2p2*m2*M2*p1k2 - 6*eg1*l1k2*l1l2*M2*p1k1*p1k2 -
			        6*l1k1*m2*M2*p1k1*p1k2 - 3*l1l2*m2*M2*p1k1*p1k2 -
			        3*l2k1*m2*M2*p1k1*p1k2 + 3*m4*M2*p1k1*p1k2 +
			        6*l1k2*l2k1*M2*p1k1*p1l1 - 6*l1k1*l2k2*M2*p1k1*p1l1 -
			        3*l1k2*m2*M2*p1k1*p1l1 + 6*l1k1*l2k1*M2*p1k2*p1l1 +
			        6*eg1*l1k2*l2k1*M2*p1k2*p1l1 - 3*l1k1*m2*M2*p1k2*p1l1 +
			        3*l2k1*m2*M2*p1k2*p1l1 - 6*l1k2*l1l2*M2*p1k1*p1l2 +
			        6*l1k2*l2k1*M2*p1k1*p1l2 - 6*l1k1*l2k2*M2*p1k1*p1l2 +
			        3*l2k2*m2*M2*p1k1*p1l2 - 6*eg1*l1k1*l1k2*M2*p1k2*p1l2 -
			        6*l1k1*l1l2*M2*p1k2*p1l2 - 12*eg1*l1k2*l1l2*M2*p1k2*p1l2 +
			        6*eg1*l1k2*l2k1*M2*p1k2*p1l2 - 6*eg1*l1k1*l2k2*M2*p1k2*p1l2 -
			        3*l1k1*m2*M2*p1k2*p1l2 - 6*l1l2*m2*M2*p1k2*p1l2 +
			        3*l2k1*m2*M2*p1k2*p1l2 + 4*l1k1*l1k2*l2k1*M2*p1p2 +
			        4*l1k2*l1l2*l2k1*M2*p1p2 + 4*eg1*l1k1*l1k2*l2k2*M2*p1p2 +
			        4*l1k1*l1l2*l2k2*M2*p1p2 + 8*eg1*l1k2*l1l2*l2k2*M2*p1p2 +
			        4*l1k1*l2k1*l2k2*M2*p1p2 - 4*eg1*l1k2*l2k1*l2k2*M2*p1p2 -
			        8*l1k1*l1k2*m2*M2*p1p2 - 12*l1k2*l1l2*m2*M2*p1p2 +
			        4*l1k2*l2k1*m2*M2*p1p2 + 8*l1k1*l2k2*m2*M2*p1p2 +
			        4*l1l2*l2k2*m2*M2*p1p2 - 4*l2k1*l2k2*m2*M2*p1p2 -
			        2*l1k2*l1p2*l2k1*p1k1*p1p2 + 2*l1k1*l1p2*l2k2*p1k1*p1p2 +
			        2*l1k2*l1l2*l2p2*p1k1*p1p2 - 2*l1k2*l2k1*l2p2*p1k1*p1p2 +
			        2*l1k1*l2k2*l2p2*p1k1*p1p2 + l1k2*l1p2*m2*p1k1*p1p2 -
			        l2k2*l2p2*m2*p1k1*p1p2 - 2*l1k1*l1p2*l2k1*p1k2*p1p2 -
			        2*eg1*l1k2*l1p2*l2k1*p1k2*p1p2 + 2*eg1*l1k1*l1k2*l2p2*p1k2*p1p2 +
			        2*l1k1*l1l2*l2p2*p1k2*p1p2 + 4*eg1*l1k2*l1l2*l2p2*p1k2*p1p2 -
			        2*eg1*l1k2*l2k1*l2p2*p1k2*p1p2 + 2*eg1*l1k1*l2k2*l2p2*p1k2*p1p2 +
			        l1k1*l1p2*m2*p1k2*p1p2 - l1p2*l2k1*m2*p1k2*p1p2 +
			        l1k1*l2p2*m2*p1k2*p1p2 + 2*l1l2*l2p2*m2*p1k2*p1p2 -
			        l2k1*l2p2*m2*p1k2*p1p2 + 2*eg1*l1k2*l1l2*p1k1*p1k2*p1p2 +
			        2*l1k1*m2*p1k1*p1k2*p1p2 + l1l2*m2*p1k1*p1k2*p1p2 +
			        l2k1*m2*p1k1*p1k2*p1p2 - m4*p1k1*p1k2*p1p2 -
			        2*l1k2*l2k1*p1k1*p1l1*p1p2 + 2*l1k1*l2k2*p1k1*p1l1*p1p2 +
			        l1k2*m2*p1k1*p1l1*p1p2 - 2*l1k1*l2k1*p1k2*p1l1*p1p2 -
			        2*eg1*l1k2*l2k1*p1k2*p1l1*p1p2 + l1k1*m2*p1k2*p1l1*p1p2 -
			        l2k1*m2*p1k2*p1l1*p1p2 + 2*l1k2*l1l2*p1k1*p1l2*p1p2 -
			        2*l1k2*l2k1*p1k1*p1l2*p1p2 + 2*l1k1*l2k2*p1k1*p1l2*p1p2 -
			        l2k2*m2*p1k1*p1l2*p1p2 + 2*eg1*l1k1*l1k2*p1k2*p1l2*p1p2 +
			        2*l1k1*l1l2*p1k2*p1l2*p1p2 + 4*eg1*l1k2*l1l2*p1k2*p1l2*p1p2 -
			        2*eg1*l1k2*l2k1*p1k2*p1l2*p1p2 + 2*eg1*l1k1*l2k2*p1k2*p1l2*p1p2 +
			        l1k1*m2*p1k2*p1l2*p1p2 + 2*l1l2*m2*p1k2*p1l2*p1p2 -
			        l2k1*m2*p1k2*p1l2*p1p2 + 6*l1k1*l1p2*l2k1*M2*p2k2 +
			        6*eg1*l1k2*l1p2*l2k1*M2*p2k2 - 6*eg1*l1k1*l1k2*l2p2*M2*p2k2 -
			        6*l1k1*l1l2*l2p2*M2*p2k2 - 12*eg1*l1k2*l1l2*l2p2*M2*p2k2 +
			        6*eg1*l1k2*l2k1*l2p2*M2*p2k2 - 6*eg1*l1k1*l2k2*l2p2*M2*p2k2 -
			        3*l1k1*l1p2*m2*M2*p2k2 + 3*l1p2*l2k1*m2*M2*p2k2 -
			        3*l1k1*l2p2*m2*M2*p2k2 - 6*l1l2*l2p2*m2*M2*p2k2 +
			        3*l2k1*l2p2*m2*M2*p2k2 + 2*eg1*l1k2*l1l2*M2*p1k1*p2k2 +
			        2*l1k1*m2*M2*p1k1*p2k2 + l1l2*m2*M2*p1k1*p2k2 +
			        l2k1*m2*M2*p1k1*p2k2 - m4*M2*p1k1*p2k2 -
			        2*l1k1*l2k1*M2*p1l1*p2k2 - 2*eg1*l1k2*l2k1*M2*p1l1*p2k2 +
			        l1k1*m2*M2*p1l1*p2k2 - l2k1*m2*M2*p1l1*p2k2 +
			        2*eg1*l1k1*l1k2*M2*p1l2*p2k2 + 2*l1k1*l1l2*M2*p1l2*p2k2 +
			        4*eg1*l1k2*l1l2*M2*p1l2*p2k2 - 2*eg1*l1k2*l2k1*M2*p1l2*p2k2 +
			        2*eg1*l1k1*l2k2*M2*p1l2*p2k2 + l1k1*m2*M2*p1l2*p2k2 +
			        2*l1l2*m2*M2*p1l2*p2k2 - l2k1*m2*M2*p1l2*p2k2 -
			        2*l1k1*l1p2*l2k1*p1p2*p2k2 - 2*eg1*l1k2*l1p2*l2k1*p1p2*p2k2 +
			        2*eg1*l1k1*l1k2*l2p2*p1p2*p2k2 + 2*l1k1*l1l2*l2p2*p1p2*p2k2 +
			        4*eg1*l1k2*l1l2*l2p2*p1p2*p2k2 - 2*eg1*l1k2*l2k1*l2p2*p1p2*p2k2 +
			        2*eg1*l1k1*l2k2*l2p2*p1p2*p2k2 + l1k1*l1p2*m2*p1p2*p2k2 -
			        l1p2*l2k1*m2*p1p2*p2k2 + l1k1*l2p2*m2*p1p2*p2k2 +
			        2*l1l2*l2p2*m2*p1p2*p2k2 - l2k1*l2p2*m2*p1p2*p2k2 +
			        2*eg1*l1k2*l1l2*p1k1*p1p2*p2k2 + 2*l1k1*m2*p1k1*p1p2*p2k2 +
			        l1l2*m2*p1k1*p1p2*p2k2 + l2k1*m2*p1k1*p1p2*p2k2 -
			        m4*p1k1*p1p2*p2k2 - 2*l1k1*l2k1*p1l1*p1p2*p2k2 -
			        2*eg1*l1k2*l2k1*p1l1*p1p2*p2k2 + l1k1*m2*p1l1*p1p2*p2k2 -
			        l2k1*m2*p1l1*p1p2*p2k2 + 2*eg1*l1k1*l1k2*p1l2*p1p2*p2k2 +
			        2*l1k1*l1l2*p1l2*p1p2*p2k2 + 4*eg1*l1k2*l1l2*p1l2*p1p2*p2k2 -
			        2*eg1*l1k2*l2k1*p1l2*p1p2*p2k2 + 2*eg1*l1k1*l2k2*p1l2*p1p2*p2k2 +
			        l1k1*m2*p1l2*p1p2*p2k2 + 2*l1l2*m2*p1l2*p1p2*p2k2 -
			        l2k1*m2*p1l2*p1p2*p2k2 +
			        k1p2*(-6*l1k1*l1p2*l2k2*M2 - 6*l1k1*l2k2*l2p2*M2 +
			           3*l2k2*l2p2*m2*M2 - 2*l2k2*m2*M2*p1k1 +
			           2*l1k1*m2*M2*p1k2 + l1l2*m2*M2*p1k2 + l2k1*m2*M2*p1k2 -
			           m4*M2*p1k2 + 2*l1k1*l2k2*M2*p1l1 + 2*l1k1*l2k2*M2*p1l2 -
			           l2k2*m2*M2*p1l2 + 2*l1k1*l1p2*l2k2*p1p2 +
			           2*l1k1*l2k2*l2p2*p1p2 - l2k2*l2p2*m2*p1p2 -
			           2*l2k2*m2*p1k1*p1p2 + 2*l1k1*m2*p1k2*p1p2 +
			           l1l2*m2*p1k2*p1p2 + l2k1*m2*p1k2*p1p2 - m4*p1k2*p1p2 +
			           2*l1k1*l2k2*p1l1*p1p2 + 2*l1k1*l2k2*p1l2*p1p2 -
			           l2k2*m2*p1l2*p1p2 +
			           k1k2*(l1l2*(l1p2*(6*M2 - 2*p1p2) - 2*p1l1*(M2 + p1p2)) +
			              m2*(l2p2*(-3*M2 + p1p2) + (2*p1l1 + p1l2)*(M2 + p1p2) +
			                 l1p2*(-6*M2 + 2*p1p2))) - 6*l1k1*m2*M2*p2k2 -
			           3*l1l2*m2*M2*p2k2 - 3*l2k1*m2*M2*p2k2 + 3*m4*M2*p2k2 +
			           2*l1k1*m2*p1p2*p2k2 + l1l2*m2*p1p2*p2k2 + l2k1*m2*p1p2*p2k2 -
			           m4*p1p2*p2k2 + l1k2*
			            (6*l2k1*l2p2*M2 - 4*m2*M2*p1k1 - 2*l2k1*M2*p1l1 +
			              m2*M2*p1l1 - 2*l2k1*M2*p1l2 +
			              l1p2*(2*l2k1 - m2)*(3*M2 - p1p2) - 2*l2k1*l2p2*p1p2 -
			              4*m2*p1k1*p1p2 - 2*l2k1*p1l1*p1p2 + m2*p1l1*p1p2 -
			              2*l2k1*p1l2*p1p2 +
			              2*l1l2*(l2p2*(-3*M2 + p1p2) +
			                 M2*(2*p1k1 + p1l2 + eg1*(p1k2 - 3*p2k2)) +
			                 p1p2*(2*p1k1 + p1l2 + eg1*(p1k2 + p2k2))))) +
			        (-2*l1k2*l1l2 + 2*l1k2*m2 + l2k2*m2)*(3*M2 - p1p2)*pow(k1p2,2) +
			        2*l2p2*M2*p1k2*pow(l1k1,2) - 6*M2*p1k2*p1l2*pow(l1k1,2) -
			        4*l2k2*M2*p1p2*pow(l1k1,2) + 2*l2p2*p1k2*p1p2*pow(l1k1,2) +
			        2*p1k2*p1l2*p1p2*pow(l1k1,2) - 6*l2p2*M2*p2k2*pow(l1k1,2) +
			        2*M2*p1l2*p2k2*pow(l1k1,2) + 2*l2p2*p1p2*p2k2*pow(l1k1,2) +
			        2*p1l2*p1p2*p2k2*pow(l1k1,2) - 4*eg1*l2k1*M2*p1p2*pow(l1k2,2) -
			        4*l1k2*M2*p1p2*pow(l2k1,2) + 4*eg1*l1k1*M2*p1p2*pow(l2k2,2) -
			        6*l1k1*l1k2*l2k1*M4 - 2*l1k2*l1l2*l2k1*M4 -
			        2*eg1*l1k1*l1k2*l2k2*M4 - 2*l1k1*l1l2*l2k2*M4 -
			        4*eg1*l1k2*l1l2*l2k2*M4 - 2*l1k1*l2k1*l2k2*M4 +
			        2*eg1*l1k2*l2k1*l2k2*M4 + 8*l1k1*l1k2*m2*M4 +
			        10*l1k2*l1l2*m2*M4 - 4*l1k2*l2k1*m2*M4 -
			        6*l1k1*l2k2*m2*M4 - 2*l1l2*l2k2*m2*M4 +
			        2*l2k1*l2k2*m2*M4 + 6*l2k2*pow(l1k1,2)*M4 +
			        2*eg1*l2k1*pow(l1k2,2)*M4 + 2*l1k2*pow(l2k1,2)*M4 -
			        2*eg1*l1k1*pow(l2k2,2)*M4 - 6*l1k2*l1l2*M2*pow(p1k1,2) +
			        6*l1k2*m2*M2*pow(p1k1,2) + 3*l2k2*m2*M2*pow(p1k1,2) +
			        2*l1k2*l1l2*p1p2*pow(p1k1,2) - 2*l1k2*m2*p1p2*pow(p1k1,2) -
			        l2k2*m2*p1p2*pow(p1k1,2) + 2*l1k1*l1k2*l2k1*pow(p1p2,2) -
			        2*l1k2*l1l2*l2k1*pow(p1p2,2) - 2*eg1*l1k1*l1k2*l2k2*pow(p1p2,2) -
			        2*l1k1*l1l2*l2k2*pow(p1p2,2) - 4*eg1*l1k2*l1l2*l2k2*pow(p1p2,2) -
			        2*l1k1*l2k1*l2k2*pow(p1p2,2) + 2*eg1*l1k2*l2k1*l2k2*pow(p1p2,2) +
			        2*l1k2*l1l2*m2*pow(p1p2,2) - 2*l1k1*l2k2*m2*pow(p1p2,2) -
			        2*l1l2*l2k2*m2*pow(p1p2,2) + 2*l2k1*l2k2*m2*pow(p1p2,2) -
			        2*l2k2*pow(l1k1,2)*pow(p1p2,2) +
			        2*eg1*l2k1*pow(l1k2,2)*pow(p1p2,2) +
			        2*l1k2*pow(l2k1,2)*pow(p1p2,2) -
			        2*eg1*l1k1*pow(l2k2,2)*pow(p1p2,2) -
			        k1k2*(-(l2p2*m2*M2*p1k1) + 2*eg1*l1l2*l2p2*M2*p1k2 +
			           2*l1l2*l2p2*M2*p1l1 - 2*l2k1*l2p2*M2*p1l1 -
			           6*l1l2*M2*p1k1*p1l1 + 6*m2*M2*p1k1*p1l1 -
			           2*l2p2*m2*M2*p1l2 + 3*m2*M2*p1k1*p1l2 -
			           6*eg1*l1l2*M2*p1k2*p1l2 - 6*l1l2*M2*p1l1*p1l2 +
			           6*l2k1*M2*p1l1*p1l2 - 4*eg1*l1k2*l1l2*M2*p1p2 -
			           4*l1l2*l2k1*M2*p1p2 + 4*eg1*l1l2*l2k2*M2*p1p2 -
			           20*l1l2*m2*M2*p1p2 + 16*l2k1*m2*M2*p1p2 - 4*m4*M2*p1p2 -
			           l2p2*m2*p1k1*p1p2 + 2*eg1*l1l2*l2p2*p1k2*p1p2 +
			           2*l1l2*l2p2*p1l1*p1p2 - 2*l2k1*l2p2*p1l1*p1p2 +
			           2*l1l2*p1k1*p1l1*p1p2 - 2*m2*p1k1*p1l1*p1p2 -
			           2*l2p2*m2*p1l2*p1p2 - m2*p1k1*p1l2*p1p2 +
			           2*eg1*l1l2*p1k2*p1l2*p1p2 + 2*l1l2*p1l1*p1l2*p1p2 -
			           2*l2k1*p1l1*p1l2*p1p2 -
			           2*l1p2*(-3*l2k1*l2p2*M2 + m2*M2*p1k1 + 2*l2k1*M2*p1l1 -
			              m2*M2*p1l1 + l2k1*M2*p1l2 + l2k1*l2p2*p1p2 +
			              m2*p1k1*p1p2 + 2*l2k1*p1l1*p1p2 - m2*p1l1*p1p2 +
			              l2k1*p1l2*p1p2 +
			              l1k1*(l2p2*(3*M2 - p1p2) - p1l2*(M2 + p1p2)) +
			              l1l2*(l2p2*(3*M2 - p1p2) - (p1k1 + p1l2)*(M2 + p1p2))) -
			           6*eg1*l1l2*l2p2*M2*p2k2 + 2*eg1*l1l2*M2*p1l2*p2k2 +
			           2*eg1*l1l2*l2p2*p1p2*p2k2 + 2*eg1*l1l2*p1l2*p1p2*p2k2 +
			           4*M2*p1p2*pow(l1l2,2) +
			           (2*l2k1 - m2)*(3*M2 - p1p2)*pow(l1p2,2) +
			           3*m2*M2*pow(l2p2,2) - m2*p1p2*pow(l2p2,2) +
			           2*eg1*l1k2*l1l2*M4 + 2*l1l2*l2k1*M4 -
			           2*eg1*l1l2*l2k2*M4 + 16*l1l2*m2*M4 -
			           16*l2k1*m2*M4 + 4*m4*M4 -
			           2*pow(l1l2,2)*M4 + 6*l2k1*M2*pow(p1l1,2) -
			           3*m2*M2*pow(p1l1,2) - 2*l2k1*p1p2*pow(p1l1,2) +
			           m2*p1p2*pow(p1l1,2) + 3*m2*M2*pow(p1l2,2) -
			           m2*p1p2*pow(p1l2,2) + 2*eg1*l1k2*l1l2*pow(p1p2,2) +
			           2*l1l2*l2k1*pow(p1p2,2) - 2*eg1*l1l2*l2k2*pow(p1p2,2) +
			           4*l1l2*m2*pow(p1p2,2) - 2*pow(l1l2,2)*pow(p1p2,2) +
			           2*l1k1*(-3*M2*p1l1*p1l2 - 2*m2*M2*p1p2 + p1l1*p1l2*p1p2 +
			              l2p2*(p1l1 + 2*p1l2)*(M2 + p1p2) +
			              (-3*M2 + p1p2)*pow(l2p2,2) + m2*M4 -
			              3*M2*pow(p1l2,2) + p1p2*pow(p1l2,2) +
			              l1l2*(-2*M2*p1p2 + 3*M4 - pow(p1p2,2)) +
			              m2*pow(p1p2,2))))) +
			  8*eg1*pow(l2k1,-2)*pow(l2k1 + eg1*(k1k2 + l2k2),-2)/M2*
			   (-4*M2*pow(f1,2)*(l1k1*m2*(l2k1 + m2)*M2 - 4*l2k1*m4*M2 -
			        2*m6*M2 - l1p2*l2k1*m2*p1k1 - l1p2*m4*p1k1 -
			        k1p2*l2k1*m2*p1l1 - l2k1*l2p2*m2*p1l1 - k1p2*m4*p1l1 -
			        l2p2*m4*p1l1 - l1p2*l2k1*m2*p1l2 - l1p2*m4*p1l2 +
			        2*l2k1*m4*p1p2 + m6*p1p2 - 6*m2*M2*pow(l2k1,2) -
			        l2p2*p1l1*pow(l2k1,2) - l1p2*p1l2*pow(l2k1,2) +
			        3*m2*p1p2*pow(l2k1,2) + l1l2*M2*(l2k1*m2 + m4 + pow(l2k1,2))) +
			     4*f1*f2*M2*(-(l1k1*l2k1*m2*M2) - l1k1*m4*M2 +
			        6*l2k1*m4*M2 + 3*m6*M2 + k1p2*l2k1*m2*p1l1 +
			        l2k1*l2p2*m2*p1l1 + k1p2*m4*p1l1 + l2p2*m4*p1l1 -
			        l2k1*m2*p1k1*p1l1 - m4*p1k1*p1l1 - l2k1*m2*p1l1*p1l2 -
			        m4*p1l1*p1l2 + l1k1*l2k1*m2*p1p2 + l1k1*m4*p1p2 -
			        6*l2k1*m4*p1p2 - 3*m6*p1p2 + 9*m2*M2*pow(l2k1,2) +
			        l2p2*p1l1*pow(l2k1,2) - p1l1*p1l2*pow(l2k1,2) -
			        9*m2*p1p2*pow(l2k1,2) -
			        l1l2*(M2 - p1p2)*(l2k1*m2 + m4 + pow(l2k1,2)) +
			        l1p2*(l2k1*m2*(-k1p2 - l2p2 + p1k1 + p1l2) +
			           m4*(-k1p2 - l2p2 + p1k1 + p1l2) + (-l2p2 + p1l2)*pow(l2k1,2))) -
			     (k1k2*(l2k1 - m2) - l2k2*m2)*pow(eg1,2)*
			      (4*f1*f2*M2*(l1k2*(M2 - p1p2) - (l1p2 - p1l1)*(p1k2 - p2k2)) +
			        4*M2*(l1k2*M2 - l1p2*p1k2 - p1l1*p2k2)*pow(f1,2) +
			        pow(f2,2)*(p1l1*(M2*(3*p1k2 - p2k2) - p1p2*(p1k2 + p2k2)) -
			           l1p2*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)) +
			           l1k2*pow(M2 - p1p2,2))) +
			     pow(f2,2)*(k1p2*l2k1*m2*M2*p1l1 + l2k1*l2p2*m2*M2*p1l1 +
			        k1p2*m4*M2*p1l1 + l2p2*m4*M2*p1l1 -
			        3*l2k1*m2*M2*p1k1*p1l1 - 3*m4*M2*p1k1*p1l1 -
			        3*l2k1*m2*M2*p1l1*p1l2 - 3*m4*M2*p1l1*p1l2 +
			        2*l1k1*l2k1*m2*M2*p1p2 + 2*l1k1*m4*M2*p1p2 -
			        12*l2k1*m4*M2*p1p2 - 6*m6*M2*p1p2 + k1p2*l2k1*m2*p1l1*p1p2 +
			        l2k1*l2p2*m2*p1l1*p1p2 + k1p2*m4*p1l1*p1p2 + l2p2*m4*p1l1*p1p2 +
			        l2k1*m2*p1k1*p1l1*p1p2 + m4*p1k1*p1l1*p1p2 +
			        l2k1*m2*p1l1*p1l2*p1p2 + m4*p1l1*p1l2*p1p2 +
			        l2p2*M2*p1l1*pow(l2k1,2) - 3*M2*p1l1*p1l2*pow(l2k1,2) -
			        18*m2*M2*p1p2*pow(l2k1,2) + l2p2*p1l1*p1p2*pow(l2k1,2) +
			        p1l1*p1l2*p1p2*pow(l2k1,2) +
			        l1p2*(l2k1*m2*(k1p2*(-3*M2 + p1p2) + l2p2*(-3*M2 + p1p2) +
			              (p1k1 + p1l2)*(M2 + p1p2)) +
			           m4*(k1p2*(-3*M2 + p1p2) + l2p2*(-3*M2 + p1p2) +
			              (p1k1 + p1l2)*(M2 + p1p2)) +
			           (l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2))*pow(l2k1,2)) -
			        l1k1*l2k1*m2*M4 - l1k1*m4*M4 +
			        10*l2k1*m4*M4 + 5*m6*M4 +
			        15*m2*pow(l2k1,2)*M4 -
			        l1l2*(l2k1*m2 + m4 + pow(l2k1,2))*pow(M2 - p1p2,2) -
			        l1k1*l2k1*m2*pow(p1p2,2) - l1k1*m4*pow(p1p2,2) +
			        2*l2k1*m4*pow(p1p2,2) + m6*pow(p1p2,2) +
			        3*m2*pow(l2k1,2)*pow(p1p2,2)) +
			     eg1*(-4*M2*pow(f1,2)*(-(l1k1*l2k1*l2k2*M2) -
			           4*l2k1*l2k2*m2*M2 + l1k2*m4*M2 - 2*l2k2*m4*M2 +
			           l1p2*l2k1*l2k2*p1k1 - l1p2*m4*p1k2 + k1p2*l2k1*l2k2*p1l1 +
			           2*l2k1*l2k2*m2*p1p2 + l2k2*m4*p1p2 +
			           k1k2*(l1l2*l2k1*M2 -
			              l2k1*(l2p2*p1l1 + l1p2*p1l2 + m2*(6*M2 - 3*p1p2)) +
			              m4*(-2*M2 + p1p2)) - m4*p1l1*p2k2 +
			           l1k2*M2*pow(l2k1,2) - l1p2*p1k2*pow(l2k1,2) -
			           p1l1*p2k2*pow(l2k1,2)) +
			        4*f1*f2*M2*(l1k1*l2k1*l2k2*M2 + 6*l2k1*l2k2*m2*M2 -
			           l1k2*m4*M2 + 3*l2k2*m4*M2 - l1p2*l2k1*l2k2*p1k1 +
			           l1p2*m4*p1k2 + k1p2*l2k1*l2k2*(l1p2 - p1l1) +
			           l2k1*l2k2*p1k1*p1l1 - m4*p1k2*p1l1 - l1k1*l2k1*l2k2*p1p2 -
			           6*l2k1*l2k2*m2*p1p2 + l1k2*m4*p1p2 - 3*l2k2*m4*p1p2 +
			           k1k2*(9*l2k1*m2*M2 + 3*m4*M2 + l2k1*l2p2*p1l1 -
			              l2k1*p1l1*p1l2 + l1p2*l2k1*(-l2p2 + p1l2) -
			              9*l2k1*m2*p1p2 - 3*m4*p1p2 + l1l2*l2k1*(-M2 + p1p2)) -
			           l1p2*m4*p2k2 + m4*p1l1*p2k2 - l1k2*M2*pow(l2k1,2) +
			           l1p2*p1k2*pow(l2k1,2) - p1k2*p1l1*pow(l2k1,2) +
			           l1k2*p1p2*pow(l2k1,2) - l1p2*p2k2*pow(l2k1,2) +
			           p1l1*p2k2*pow(l2k1,2)) +
			        pow(f2,2)*(-(l1p2*l2k1*l2k2*M2*p1k1) + l1p2*m4*M2*p1k2 +
			           3*l2k1*l2k2*M2*p1k1*p1l1 - 3*m4*M2*p1k2*p1l1 -
			           2*l1k1*l2k1*l2k2*M2*p1p2 - 12*l2k1*l2k2*m2*M2*p1p2 +
			           2*l1k2*m4*M2*p1p2 - 6*l2k2*m4*M2*p1p2 -
			           l1p2*l2k1*l2k2*p1k1*p1p2 + l1p2*m4*p1k2*p1p2 -
			           l2k1*l2k2*p1k1*p1l1*p1p2 + m4*p1k2*p1l1*p1p2 +
			           k1p2*l2k1*l2k2*(l1p2*(3*M2 - p1p2) - p1l1*(M2 + p1p2)) -
			           3*l1p2*m4*M2*p2k2 + m4*M2*p1l1*p2k2 + l1p2*m4*p1p2*p2k2 +
			           m4*p1l1*p1p2*p2k2 + l1p2*M2*p1k2*pow(l2k1,2) -
			           3*M2*p1k2*p1l1*pow(l2k1,2) + 2*l1k2*M2*p1p2*pow(l2k1,2) +
			           l1p2*p1k2*p1p2*pow(l2k1,2) + p1k2*p1l1*p1p2*pow(l2k1,2) -
			           3*l1p2*M2*p2k2*pow(l2k1,2) + M2*p1l1*p2k2*pow(l2k1,2) +
			           l1p2*p1p2*p2k2*pow(l2k1,2) + p1l1*p1p2*p2k2*pow(l2k1,2) +
			           l1k1*l2k1*l2k2*M4 + 10*l2k1*l2k2*m2*M4 -
			           l1k2*m4*M4 + 5*l2k2*m4*M4 -
			           l1k2*pow(l2k1,2)*M4 + l1k1*l2k1*l2k2*pow(p1p2,2) +
			           2*l2k1*l2k2*m2*pow(p1p2,2) - l1k2*m4*pow(p1p2,2) +
			           l2k2*m4*pow(p1p2,2) - l1k2*pow(l2k1,2)*pow(p1p2,2) +
			           k1k2*(l2k1*l2p2*M2*p1l1 - 3*l2k1*M2*p1l1*p1l2 -
			              18*l2k1*m2*M2*p1p2 - 6*m4*M2*p1p2 +
			              l2k1*l2p2*p1l1*p1p2 + l2k1*p1l1*p1l2*p1p2 +
			              l1p2*l2k1*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) +
			              15*l2k1*m2*M4 + 5*m4*M4 -
			              l1l2*l2k1*pow(M2 - p1p2,2) + 3*l2k1*m2*pow(p1p2,2) +
			              m4*pow(p1p2,2))))) -
			  8*eg1*pow(l2k1,-1)*pow(l2k2,-1)*pow(l2k1 + eg1*(k1k2 + l2k2),-2)*
			   pow(M,-2)*(4*M2*pow(f1,2)*
			      (-2*l1k1*l2k1*l2k2*M2 - l1k2*l2k1*m2*M2 + l1k1*l2k2*m2*M2 +
			        l1l2*l2k2*m2*M2 - 2*l2k2*m4*M2 + 2*l1p2*l2k1*l2k2*p1k1 -
			        l1p2*l2k2*m2*p1k1 + l1p2*l2k1*m2*p1k2 + 2*k1p2*l2k1*l2k2*p1l1 -
			        k1p2*l2k2*m2*p1l1 - l2k2*l2p2*m2*p1l1 - l1p2*l2k2*m2*p1l2 +
			        k1k2*l2k1*(-2*l1k1*M2 + 2*(l1p2*p1k1 + k1p2*p1l1) +
			           m2*(2*M2 - p1p2)) + l2k2*m4*p1p2 + l2k1*m2*p1l1*p2k2 -
			        2*l2k2*(k1k2 + l2k2)*(l1k2*M2 - l1p2*p1k2 - p1l1*p2k2)*
			         pow(eg1,2) + eg1*(2*l2k2*
			            (-(l1k2*l2k1*M2) - l1k1*l2k2*M2 + l1p2*l2k2*p1k1 +
			              l1p2*l2k1*p1k2 + k1p2*l2k2*p1l1 + l2k1*p1l1*p2k2) +
			           k1k2*(-2*l1k1*l2k2*M2 + 2*l1l2*l2k2*M2 - 2*l2k2*m2*M2 -
			              l1k2*(2*l2k1 + m2)*M2 + 2*l1p2*l2k2*p1k1 +
			              2*l1p2*l2k1*p1k2 + l1p2*m2*p1k2 + 2*k1p2*l2k2*p1l1 -
			              2*l2k2*l2p2*p1l1 - 2*l1p2*l2k2*p1l2 + l2k2*m2*p1p2 +
			              2*l2k1*p1l1*p2k2 + m2*p1l1*p2k2) +
			           (2*l1l2*M2 - 2*(l2p2*p1l1 + l1p2*p1l2) + m2*(-2*M2 + p1p2))*
			            pow(k1k2,2))) + 4*f1*f2*M2*
			      (-2*k1p2*l1p2*l2k1*l2k2 + k1p2*l1p2*l2k2*m2 + l1p2*l2k2*l2p2*m2 -
			        2*l1k1*l2k1*l2k2*M2 - 2*eg1*l1k2*l2k1*l2k2*M2 -
			        l1k2*l2k1*m2*M2 + l1k1*l2k2*m2*M2 + l1l2*l2k2*m2*M2 -
			        3*l2k2*m4*M2 + 2*l1p2*l2k1*l2k2*p1k1 - l1p2*l2k2*m2*p1k1 +
			        2*eg1*l1p2*l2k1*l2k2*p1k2 + l1p2*l2k1*m2*p1k2 +
			        2*k1p2*l2k1*l2k2*p1l1 - k1p2*l2k2*m2*p1l1 - l2k2*l2p2*m2*p1l1 -
			        2*l2k1*l2k2*p1k1*p1l1 + l2k2*m2*p1k1*p1l1 -
			        2*eg1*l2k1*l2k2*p1k2*p1l1 - l2k1*m2*p1k2*p1l1 -
			        l1p2*l2k2*m2*p1l2 + l2k2*m2*p1l1*p1l2 + 2*l1k1*l2k1*l2k2*p1p2 +
			        2*eg1*l1k2*l2k1*l2k2*p1p2 + l1k2*l2k1*m2*p1p2 -
			        l1k1*l2k2*m2*p1p2 - l1l2*l2k2*m2*p1p2 + 3*l2k2*m4*p1p2 -
			        2*eg1*l1p2*l2k1*l2k2*p2k2 - l1p2*l2k1*m2*p2k2 +
			        2*eg1*l2k1*l2k2*p1l1*p2k2 + l2k1*m2*p1l1*p2k2 +
			        k1k2*(-2*k1p2*(l2k1 + eg1*l2k2)*(l1p2 - p1l1) +
			           l2k1*(2*p1k1*(l1p2 - p1l1) - 2*l1k1*(M2 - p1p2) +
			              3*m2*(M2 - p1p2)) +
			           eg1*(-2*l1k1*l2k2*M2 + 2*l1l2*l2k2*M2 - 3*l2k2*m2*M2 -
			              2*l2k2*l2p2*p1l1 - 2*l2k2*p1k1*p1l1 - 2*l2k1*p1k2*p1l1 -
			              m2*p1k2*p1l1 + 2*l2k2*p1l1*p1l2 -
			              l1k2*(2*l2k1 + m2)*(M2 - p1p2) + 2*l1k1*l2k2*p1p2 -
			              2*l1l2*l2k2*p1p2 + 3*l2k2*m2*p1p2 +
			              l1p2*(2*l2k2*(l2p2 + p1k1 - p1l2) +
			                 (2*l2k1 + m2)*(p1k2 - p2k2)) + 2*l2k1*p1l1*p2k2 +
			              m2*p1l1*p2k2) -
			           2*l2k2*(l1k2*(M2 - p1p2) - (l1p2 - p1l1)*(p1k2 - p2k2))*
			            pow(eg1,2)) + eg1*(-3*m2*M2 - 2*l2p2*p1l1 +
			           2*l1p2*(l2p2 - p1l2) + 2*p1l1*p1l2 + 2*l1l2*(M2 - p1p2) +
			           3*m2*p1p2)*pow(k1k2,2) - 2*eg1*k1p2*l1p2*pow(l2k2,2) -
			        2*eg1*l1k1*M2*pow(l2k2,2) + 2*eg1*l1p2*p1k1*pow(l2k2,2) +
			        2*eg1*k1p2*p1l1*pow(l2k2,2) - 2*eg1*p1k1*p1l1*pow(l2k2,2) +
			        2*eg1*l1k1*p1p2*pow(l2k2,2) - 2*l1k2*M2*pow(eg1,2)*pow(l2k2,2) +
			        2*l1p2*p1k2*pow(eg1,2)*pow(l2k2,2) -
			        2*p1k2*p1l1*pow(eg1,2)*pow(l2k2,2) +
			        2*l1k2*p1p2*pow(eg1,2)*pow(l2k2,2) -
			        2*l1p2*p2k2*pow(eg1,2)*pow(l2k2,2) +
			        2*p1l1*p2k2*pow(eg1,2)*pow(l2k2,2)) +
			     pow(f2,2)*(-6*k1p2*l1p2*l2k1*l2k2*M2 + 3*k1p2*l1p2*l2k2*m2*M2 +
			        3*l1p2*l2k2*l2p2*m2*M2 + 2*l1p2*l2k1*l2k2*M2*p1k1 -
			        l1p2*l2k2*m2*M2*p1k1 + 2*eg1*l1p2*l2k1*l2k2*M2*p1k2 +
			        l1p2*l2k1*m2*M2*p1k2 + 2*k1p2*l2k1*l2k2*M2*p1l1 -
			        k1p2*l2k2*m2*M2*p1l1 - l2k2*l2p2*m2*M2*p1l1 -
			        6*l2k1*l2k2*M2*p1k1*p1l1 + 3*l2k2*m2*M2*p1k1*p1l1 -
			        6*eg1*l2k1*l2k2*M2*p1k2*p1l1 - 3*l2k1*m2*M2*p1k2*p1l1 -
			        l1p2*l2k2*m2*M2*p1l2 + 3*l2k2*m2*M2*p1l1*p1l2 +
			        2*k1p2*l1p2*l2k1*l2k2*p1p2 - k1p2*l1p2*l2k2*m2*p1p2 -
			        l1p2*l2k2*l2p2*m2*p1p2 + 4*l1k1*l2k1*l2k2*M2*p1p2 +
			        4*eg1*l1k2*l2k1*l2k2*M2*p1p2 + 2*l1k2*l2k1*m2*M2*p1p2 -
			        2*l1k1*l2k2*m2*M2*p1p2 - 2*l1l2*l2k2*m2*M2*p1p2 +
			        6*l2k2*m4*M2*p1p2 + 2*l1p2*l2k1*l2k2*p1k1*p1p2 -
			        l1p2*l2k2*m2*p1k1*p1p2 + 2*eg1*l1p2*l2k1*l2k2*p1k2*p1p2 +
			        l1p2*l2k1*m2*p1k2*p1p2 + 2*k1p2*l2k1*l2k2*p1l1*p1p2 -
			        k1p2*l2k2*m2*p1l1*p1p2 - l2k2*l2p2*m2*p1l1*p1p2 +
			        2*l2k1*l2k2*p1k1*p1l1*p1p2 - l2k2*m2*p1k1*p1l1*p1p2 +
			        2*eg1*l2k1*l2k2*p1k2*p1l1*p1p2 + l2k1*m2*p1k2*p1l1*p1p2 -
			        l1p2*l2k2*m2*p1l2*p1p2 - l2k2*m2*p1l1*p1l2*p1p2 -
			        6*eg1*l1p2*l2k1*l2k2*M2*p2k2 - 3*l1p2*l2k1*m2*M2*p2k2 +
			        2*eg1*l2k1*l2k2*M2*p1l1*p2k2 + l2k1*m2*M2*p1l1*p2k2 +
			        2*eg1*l1p2*l2k1*l2k2*p1p2*p2k2 + l1p2*l2k1*m2*p1p2*p2k2 +
			        2*eg1*l2k1*l2k2*p1l1*p1p2*p2k2 + l2k1*m2*p1l1*p1p2*p2k2 -
			        6*eg1*k1p2*l1p2*M2*pow(l2k2,2) +
			        2*eg1*l1p2*M2*p1k1*pow(l2k2,2) +
			        2*eg1*k1p2*M2*p1l1*pow(l2k2,2) -
			        6*eg1*M2*p1k1*p1l1*pow(l2k2,2) +
			        2*eg1*k1p2*l1p2*p1p2*pow(l2k2,2) +
			        4*eg1*l1k1*M2*p1p2*pow(l2k2,2) +
			        2*eg1*l1p2*p1k1*p1p2*pow(l2k2,2) +
			        2*eg1*k1p2*p1l1*p1p2*pow(l2k2,2) +
			        2*eg1*p1k1*p1l1*p1p2*pow(l2k2,2) +
			        2*l1p2*M2*p1k2*pow(eg1,2)*pow(l2k2,2) -
			        6*M2*p1k2*p1l1*pow(eg1,2)*pow(l2k2,2) +
			        4*l1k2*M2*p1p2*pow(eg1,2)*pow(l2k2,2) +
			        2*l1p2*p1k2*p1p2*pow(eg1,2)*pow(l2k2,2) +
			        2*p1k2*p1l1*p1p2*pow(eg1,2)*pow(l2k2,2) -
			        6*l1p2*M2*p2k2*pow(eg1,2)*pow(l2k2,2) +
			        2*M2*p1l1*p2k2*pow(eg1,2)*pow(l2k2,2) +
			        2*l1p2*p1p2*p2k2*pow(eg1,2)*pow(l2k2,2) +
			        2*p1l1*p1p2*p2k2*pow(eg1,2)*pow(l2k2,2) -
			        2*l1k1*l2k1*l2k2*M4 - 2*eg1*l1k2*l2k1*l2k2*M4 -
			        l1k2*l2k1*m2*M4 + l1k1*l2k2*m2*M4 +
			        l1l2*l2k2*m2*M4 - 5*l2k2*m4*M4 -
			        2*eg1*l1k1*pow(l2k2,2)*M4 -
			        2*l1k2*pow(eg1,2)*pow(l2k2,2)*M4 -
			        2*l1k1*l2k1*l2k2*pow(p1p2,2) - 2*eg1*l1k2*l2k1*l2k2*pow(p1p2,2) -
			        l1k2*l2k1*m2*pow(p1p2,2) + l1k1*l2k2*m2*pow(p1p2,2) +
			        l1l2*l2k2*m2*pow(p1p2,2) - l2k2*m4*pow(p1p2,2) -
			        2*eg1*l1k1*pow(l2k2,2)*pow(p1p2,2) -
			        2*l1k2*pow(eg1,2)*pow(l2k2,2)*pow(p1p2,2) +
			        eg1*pow(k1k2,2)*(-2*l2p2*M2*p1l1 + 6*M2*p1l1*p1l2 +
			           6*m2*M2*p1p2 - 2*l2p2*p1l1*p1p2 - 2*p1l1*p1l2*p1p2 +
			           l1p2*(l2p2*(6*M2 - 2*p1p2) - 2*p1l2*(M2 + p1p2)) -
			           5*m2*M4 + 2*l1l2*pow(M2 - p1p2,2) - m2*pow(p1p2,2)) +
			        k1k2*(-2*k1p2*(l2k1 + eg1*l2k2)*
			            (l1p2*(3*M2 - p1p2) - p1l1*(M2 + p1p2)) +
			           2*l2k2*pow(eg1,2)*(l1p2*
			               (M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)) +
			              p1l1*(M2*(-3*p1k2 + p2k2) + p1p2*(p1k2 + p2k2)) -
			              l1k2*pow(M2 - p1p2,2)) +
			           eg1*(-2*l2k2*l2p2*M2*p1l1 - 6*l2k2*M2*p1k1*p1l1 -
			              6*l2k1*M2*p1k2*p1l1 - 3*m2*M2*p1k2*p1l1 +
			              6*l2k2*M2*p1l1*p1l2 + 4*l1k1*l2k2*M2*p1p2 -
			              4*l1l2*l2k2*M2*p1p2 + 6*l2k2*m2*M2*p1p2 -
			              2*l2k2*l2p2*p1l1*p1p2 + 2*l2k2*p1k1*p1l1*p1p2 +
			              2*l2k1*p1k2*p1l1*p1p2 + m2*p1k2*p1l1*p1p2 -
			              2*l2k2*p1l1*p1l2*p1p2 + 2*l2k1*M2*p1l1*p2k2 +
			              m2*M2*p1l1*p2k2 + 2*l2k1*p1l1*p1p2*p2k2 +
			              m2*p1l1*p1p2*p2k2 +
			              l1p2*(2*l2k2*(l2p2*(3*M2 - p1p2) +
			                    (p1k1 - p1l2)*(M2 + p1p2)) +
			                 (2*l2k1 + m2)*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2))\
			) - 2*l1k1*l2k2*M4 + 2*l1l2*l2k2*M4 -
			              5*l2k2*m2*M4 -
			              l1k2*(2*l2k1 + m2)*pow(M2 - p1p2,2) -
			              2*l1k1*l2k2*pow(p1p2,2) + 2*l1l2*l2k2*pow(p1p2,2) -
			              l2k2*m2*pow(p1p2,2)) +
			           l2k1*(2*p1k1*(p1l1*(-3*M2 + p1p2) + l1p2*(M2 + p1p2)) -
			              2*l1k1*pow(M2 - p1p2,2) +
			              m2*(-6*M2*p1p2 + 5*M4 + pow(p1p2,2)))))) -
			  8*eg1*pow(l1k1,-1)*pow(l1k2,-1)*pow(l1k1 + eg1*(-k1k2 + l1k2),-2)*
			   pow(M,-2)*(-4*M2*pow(f1,2)*
			      (-2*l1k1*l1k2*l2k1*M2 + l1k2*l1l2*m2*M2 - l1k2*l2k1*m2*M2 +
			        l1k1*l2k2*m2*M2 - 2*l1k2*m4*M2 + 2*l1k1*l1k2*l2p2*p1k1 +
			        l1k2*l2p2*m2*p1k1 - l1k1*l2p2*m2*p1k2 - l1k2*l2p2*m2*p1l1 +
			        2*k1p2*l1k1*l1k2*p1l2 + k1p2*l1k2*m2*p1l2 - l1k2*l1p2*m2*p1l2 +
			        k1k2*l1k1*(2*l2k1*M2 - 2*(l2p2*p1k1 + k1p2*p1l2) +
			           m2*(2*M2 - p1p2)) + l1k2*m4*p1p2 - l1k1*m2*p1l2*p2k2 +
			        2*(k1k2 - l1k2)*l1k2*(l2k2*M2 - l2p2*p1k2 - p1l2*p2k2)*
			         pow(eg1,2) + eg1*(k1k2*
			            (l1k2*(2*l1l2*M2 + 2*l2k1*M2 - 2*m2*M2 - 2*l2p2*p1k1 -
			                 2*l2p2*p1l1 - 2*k1p2*p1l2 - 2*l1p2*p1l2 + m2*p1p2) +
			              (2*l1k1 - m2)*(l2k2*M2 - l2p2*p1k2 - p1l2*p2k2)) +
			           2*l1k2*(l1k2*(-(l2k1*M2) + l2p2*p1k1 + k1p2*p1l2) +
			              l1k1*(-(l2k2*M2) + l2p2*p1k2 + p1l2*p2k2)) +
			           (-2*l1l2*M2 + 2*m2*M2 + 2*l2p2*p1l1 + 2*l1p2*p1l2 -
			              m2*p1p2)*pow(k1k2,2))) +
			     4*f1*f2*M2*(2*k1p2*l1k1*l1k2*l2p2 + k1p2*l1k2*l2p2*m2 -
			        l1k2*l1p2*l2p2*m2 + 2*l1k1*l1k2*l2k1*M2 +
			        2*eg1*l1k1*l1k2*l2k2*M2 - l1k2*l1l2*m2*M2 +
			        l1k2*l2k1*m2*M2 - l1k1*l2k2*m2*M2 + 3*l1k2*m4*M2 -
			        2*l1k1*l1k2*l2p2*p1k1 - l1k2*l2p2*m2*p1k1 -
			        2*eg1*l1k1*l1k2*l2p2*p1k2 + l1k1*l2p2*m2*p1k2 +
			        l1k2*l2p2*m2*p1l1 - 2*k1p2*l1k1*l1k2*p1l2 - k1p2*l1k2*m2*p1l2 +
			        l1k2*l1p2*m2*p1l2 + 2*l1k1*l1k2*p1k1*p1l2 + l1k2*m2*p1k1*p1l2 +
			        2*eg1*l1k1*l1k2*p1k2*p1l2 - l1k1*m2*p1k2*p1l2 -
			        l1k2*m2*p1l1*p1l2 - 2*l1k1*l1k2*l2k1*p1p2 -
			        2*eg1*l1k1*l1k2*l2k2*p1p2 + l1k2*l1l2*m2*p1p2 -
			        l1k2*l2k1*m2*p1p2 + l1k1*l2k2*m2*p1p2 - 3*l1k2*m4*p1p2 +
			        2*eg1*l1k1*l1k2*l2p2*p2k2 - l1k1*l2p2*m2*p2k2 -
			        2*eg1*l1k1*l1k2*p1l2*p2k2 + l1k1*m2*p1l2*p2k2 +
			        k1k2*(-2*k1p2*(l1k1 + eg1*l1k2)*(l2p2 - p1l2) +
			           l1k1*(-2*l2k1*M2 - 3*m2*M2 + 2*l2p2*p1k1 - 2*p1k1*p1l2 +
			              2*l2k1*p1p2 + 3*m2*p1p2) +
			           eg1*(l1k2*(-2*l1p2*l2p2 - 2*l1l2*M2 - 2*l2k1*M2 +
			                 3*m2*M2 + 2*l2p2*p1k1 + 2*l2p2*p1l1 + 2*l1p2*p1l2 -
			                 2*p1k1*p1l2 - 2*p1l1*p1l2 + 2*l1l2*p1p2 + 2*l2k1*p1p2 -
			                 3*m2*p1p2) -
			              (2*l1k1 - m2)*(l2k2*(M2 - p1p2) -
			                 (l2p2 - p1l2)*(p1k2 - p2k2))) -
			           2*l1k2*(l2k2*(M2 - p1p2) - (l2p2 - p1l2)*(p1k2 - p2k2))*
			            pow(eg1,2)) + eg1*(-3*m2*M2 - 2*l2p2*p1l1 +
			           2*l1p2*(l2p2 - p1l2) + 2*p1l1*p1l2 + 2*l1l2*(M2 - p1p2) +
			           3*m2*p1p2)*pow(k1k2,2) + 2*eg1*k1p2*l2p2*pow(l1k2,2) +
			        2*eg1*l2k1*M2*pow(l1k2,2) - 2*eg1*l2p2*p1k1*pow(l1k2,2) -
			        2*eg1*k1p2*p1l2*pow(l1k2,2) + 2*eg1*p1k1*p1l2*pow(l1k2,2) -
			        2*eg1*l2k1*p1p2*pow(l1k2,2) + 2*l2k2*M2*pow(eg1,2)*pow(l1k2,2) -
			        2*l2p2*p1k2*pow(eg1,2)*pow(l1k2,2) +
			        2*p1k2*p1l2*pow(eg1,2)*pow(l1k2,2) -
			        2*l2k2*p1p2*pow(eg1,2)*pow(l1k2,2) +
			        2*l2p2*p2k2*pow(eg1,2)*pow(l1k2,2) -
			        2*p1l2*p2k2*pow(eg1,2)*pow(l1k2,2)) +
			     pow(f2,2)*(6*k1p2*l1k1*l1k2*l2p2*M2 + 3*k1p2*l1k2*l2p2*m2*M2 -
			        3*l1k2*l1p2*l2p2*m2*M2 - 2*l1k1*l1k2*l2p2*M2*p1k1 -
			        l1k2*l2p2*m2*M2*p1k1 - 2*eg1*l1k1*l1k2*l2p2*M2*p1k2 +
			        l1k1*l2p2*m2*M2*p1k2 + l1k2*l2p2*m2*M2*p1l1 -
			        2*k1p2*l1k1*l1k2*M2*p1l2 - k1p2*l1k2*m2*M2*p1l2 +
			        l1k2*l1p2*m2*M2*p1l2 + 6*l1k1*l1k2*M2*p1k1*p1l2 +
			        3*l1k2*m2*M2*p1k1*p1l2 + 6*eg1*l1k1*l1k2*M2*p1k2*p1l2 -
			        3*l1k1*m2*M2*p1k2*p1l2 - 3*l1k2*m2*M2*p1l1*p1l2 -
			        2*k1p2*l1k1*l1k2*l2p2*p1p2 - k1p2*l1k2*l2p2*m2*p1p2 +
			        l1k2*l1p2*l2p2*m2*p1p2 - 4*l1k1*l1k2*l2k1*M2*p1p2 -
			        4*eg1*l1k1*l1k2*l2k2*M2*p1p2 + 2*l1k2*l1l2*m2*M2*p1p2 -
			        2*l1k2*l2k1*m2*M2*p1p2 + 2*l1k1*l2k2*m2*M2*p1p2 -
			        6*l1k2*m4*M2*p1p2 - 2*l1k1*l1k2*l2p2*p1k1*p1p2 -
			        l1k2*l2p2*m2*p1k1*p1p2 - 2*eg1*l1k1*l1k2*l2p2*p1k2*p1p2 +
			        l1k1*l2p2*m2*p1k2*p1p2 + l1k2*l2p2*m2*p1l1*p1p2 -
			        2*k1p2*l1k1*l1k2*p1l2*p1p2 - k1p2*l1k2*m2*p1l2*p1p2 +
			        l1k2*l1p2*m2*p1l2*p1p2 - 2*l1k1*l1k2*p1k1*p1l2*p1p2 -
			        l1k2*m2*p1k1*p1l2*p1p2 - 2*eg1*l1k1*l1k2*p1k2*p1l2*p1p2 +
			        l1k1*m2*p1k2*p1l2*p1p2 + l1k2*m2*p1l1*p1l2*p1p2 +
			        6*eg1*l1k1*l1k2*l2p2*M2*p2k2 - 3*l1k1*l2p2*m2*M2*p2k2 -
			        2*eg1*l1k1*l1k2*M2*p1l2*p2k2 + l1k1*m2*M2*p1l2*p2k2 -
			        2*eg1*l1k1*l1k2*l2p2*p1p2*p2k2 + l1k1*l2p2*m2*p1p2*p2k2 -
			        2*eg1*l1k1*l1k2*p1l2*p1p2*p2k2 + l1k1*m2*p1l2*p1p2*p2k2 +
			        6*eg1*k1p2*l2p2*M2*pow(l1k2,2) -
			        2*eg1*l2p2*M2*p1k1*pow(l1k2,2) -
			        2*eg1*k1p2*M2*p1l2*pow(l1k2,2) +
			        6*eg1*M2*p1k1*p1l2*pow(l1k2,2) -
			        2*eg1*k1p2*l2p2*p1p2*pow(l1k2,2) -
			        4*eg1*l2k1*M2*p1p2*pow(l1k2,2) -
			        2*eg1*l2p2*p1k1*p1p2*pow(l1k2,2) -
			        2*eg1*k1p2*p1l2*p1p2*pow(l1k2,2) -
			        2*eg1*p1k1*p1l2*p1p2*pow(l1k2,2) -
			        2*l2p2*M2*p1k2*pow(eg1,2)*pow(l1k2,2) +
			        6*M2*p1k2*p1l2*pow(eg1,2)*pow(l1k2,2) -
			        4*l2k2*M2*p1p2*pow(eg1,2)*pow(l1k2,2) -
			        2*l2p2*p1k2*p1p2*pow(eg1,2)*pow(l1k2,2) -
			        2*p1k2*p1l2*p1p2*pow(eg1,2)*pow(l1k2,2) +
			        6*l2p2*M2*p2k2*pow(eg1,2)*pow(l1k2,2) -
			        2*M2*p1l2*p2k2*pow(eg1,2)*pow(l1k2,2) -
			        2*l2p2*p1p2*p2k2*pow(eg1,2)*pow(l1k2,2) -
			        2*p1l2*p1p2*p2k2*pow(eg1,2)*pow(l1k2,2) +
			        2*l1k1*l1k2*l2k1*M4 + 2*eg1*l1k1*l1k2*l2k2*M4 -
			        l1k2*l1l2*m2*M4 + l1k2*l2k1*m2*M4 -
			        l1k1*l2k2*m2*M4 + 5*l1k2*m4*M4 +
			        2*eg1*l2k1*pow(l1k2,2)*M4 +
			        2*l2k2*pow(eg1,2)*pow(l1k2,2)*M4 +
			        2*l1k1*l1k2*l2k1*pow(p1p2,2) + 2*eg1*l1k1*l1k2*l2k2*pow(p1p2,2) -
			        l1k2*l1l2*m2*pow(p1p2,2) + l1k2*l2k1*m2*pow(p1p2,2) -
			        l1k1*l2k2*m2*pow(p1p2,2) + l1k2*m4*pow(p1p2,2) +
			        2*eg1*l2k1*pow(l1k2,2)*pow(p1p2,2) +
			        2*l2k2*pow(eg1,2)*pow(l1k2,2)*pow(p1p2,2) +
			        eg1*pow(k1k2,2)*(-2*l2p2*M2*p1l1 + 6*M2*p1l1*p1l2 +
			           6*m2*M2*p1p2 - 2*l2p2*p1l1*p1p2 - 2*p1l1*p1l2*p1p2 +
			           l1p2*(l2p2*(6*M2 - 2*p1p2) - 2*p1l2*(M2 + p1p2)) -
			           5*m2*M4 + 2*l1l2*pow(M2 - p1p2,2) - m2*pow(p1p2,2)) +
			        k1k2*(-2*k1p2*(l1k1 + eg1*l1k2)*
			            (l2p2*(3*M2 - p1p2) - p1l2*(M2 + p1p2)) +
			           2*l1k2*pow(eg1,2)*(l2p2*
			               (M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)) +
			              p1l2*(M2*(-3*p1k2 + p2k2) + p1p2*(p1k2 + p2k2)) -
			              l2k2*pow(M2 - p1p2,2)) +
			           l1k1*(2*p1k1*(p1l2*(-3*M2 + p1p2) + l2p2*(M2 + p1p2)) -
			              2*l2k1*pow(M2 - p1p2,2) -
			              m2*(-6*M2*p1p2 + 5*M4 + pow(p1p2,2))) +
			           eg1*((2*l1k1 - m2)*
			               (l2p2*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)) +
			                 p1l2*(M2*(-3*p1k2 + p2k2) + p1p2*(p1k2 + p2k2)) -
			                 l2k2*pow(M2 - p1p2,2)) +
			              l1k2*(2*l2p2*M2*p1k1 + 2*l2p2*M2*p1l1 - 6*M2*p1k1*p1l2 -
			                 6*M2*p1l1*p1l2 + 4*l2k1*M2*p1p2 - 6*m2*M2*p1p2 +
			                 2*l2p2*p1k1*p1p2 + 2*l2p2*p1l1*p1p2 + 2*p1k1*p1l2*p1p2 +
			                 2*p1l1*p1l2*p1p2 +
			                 2*l1p2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) -
			                 2*l2k1*M4 + 5*m2*M4 -
			                 2*l1l2*pow(M2 - p1p2,2) - 2*l2k1*pow(p1p2,2) +
			                 m2*pow(p1p2,2)))))) +
			  8*eg1*pow(l1k1,-2)*pow(l1k1 + eg1*(-k1k2 + l1k2),-2)/M2*
			   (4*f1*f2*M2*(m4*(-(l1p2*l2p2) - l1l2*M2 + l2k1*M2 + 3*m2*M2 -
			           l2p2*p1k1 + l2p2*p1l1 + k1p2*(l2p2 - p1l2) + l1p2*p1l2 +
			           p1k1*p1l2 - p1l1*p1l2 + l1l2*p1p2 - l2k1*p1p2 - 3*m2*p1p2) +
			        l1k1*m2*(l1p2*l2p2 + l1l2*M2 - l2k1*M2 - 6*m2*M2 +
			           l2p2*p1k1 - l2p2*p1l1 - l1p2*p1l2 - p1k1*p1l2 + p1l1*p1l2 +
			           k1p2*(-l2p2 + p1l2) - l1l2*p1p2 + l2k1*p1p2 + 6*m2*p1p2) +
			        (-(l1l2*M2) + 9*m2*M2 + l2p2*p1l1 - p1l1*p1l2 +
			           l1p2*(-l2p2 + p1l2) + l1l2*p1p2 - 9*m2*p1p2)*pow(l1k1,2)) -
			     4*M2*pow(f1,2)*(l1k1*m2*
			         (-(l1l2*M2) + l2k1*M2 + 4*m2*M2 - l2p2*p1k1 + l2p2*p1l1 -
			           k1p2*p1l2 + l1p2*p1l2 - 2*m2*p1p2) +
			        m4*(l1l2*M2 - l2k1*M2 - 2*m2*M2 + l2p2*p1k1 - l2p2*p1l1 +
			           k1p2*p1l2 - l1p2*p1l2 + m2*p1p2) +
			        (l1l2*M2 - l2p2*p1l1 - l1p2*p1l2 + m2*(-6*M2 + 3*p1p2))*
			         pow(l1k1,2)) - (-(l1k2*m2) + k1k2*(l1k1 + m2))*pow(eg1,2)*
			      (4*f1*f2*M2*(l2k2*(M2 - p1p2) - (l2p2 - p1l2)*(p1k2 - p2k2)) +
			        4*M2*(l2k2*M2 - l2p2*p1k2 - p1l2*p2k2)*pow(f1,2) +
			        pow(f2,2)*(p1l2*(M2*(3*p1k2 - p2k2) - p1p2*(p1k2 + p2k2)) -
			           l2p2*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)) +
			           l2k2*pow(M2 - p1p2,2))) +
			     pow(f2,2)*(l1k1*m2*(l2p2*M2*p1k1 - l2p2*M2*p1l1 -
			           3*M2*p1k1*p1l2 + 3*M2*p1l1*p1l2 - 2*l1l2*M2*p1p2 +
			           2*l2k1*M2*p1p2 + 12*m2*M2*p1p2 + l2p2*p1k1*p1p2 -
			           l2p2*p1l1*p1p2 + p1k1*p1l2*p1p2 - p1l1*p1l2*p1p2 +
			           l1p2*(l2p2*(3*M2 - p1p2) - p1l2*(M2 + p1p2)) +
			           k1p2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) +
			           l1l2*M4 - l2k1*M4 - 10*m2*M4 +
			           l1l2*pow(p1p2,2) - l2k1*pow(p1p2,2) - 2*m2*pow(p1p2,2)) +
			        m4*(-(l2p2*M2*p1k1) + l2p2*M2*p1l1 + 3*M2*p1k1*p1l2 -
			           3*M2*p1l1*p1l2 + 2*l1l2*M2*p1p2 - 2*l2k1*M2*p1p2 -
			           6*m2*M2*p1p2 - l2p2*p1k1*p1p2 + l2p2*p1l1*p1p2 -
			           p1k1*p1l2*p1p2 + p1l1*p1l2*p1p2 +
			           k1p2*(l2p2*(3*M2 - p1p2) - p1l2*(M2 + p1p2)) +
			           l1p2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) -
			           l1l2*M4 + l2k1*M4 + 5*m2*M4 -
			           l1l2*pow(p1p2,2) + l2k1*pow(p1p2,2) + m2*pow(p1p2,2)) +
			        pow(l1k1,2)*(l2p2*M2*p1l1 - 3*M2*p1l1*p1l2 - 18*m2*M2*p1p2 +
			           l2p2*p1l1*p1p2 + p1l1*p1l2*p1p2 +
			           l1p2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) +
			           15*m2*M4 - l1l2*pow(M2 - p1p2,2) + 3*m2*pow(p1p2,2))) \
			+ eg1*(4*f1*f2*M2*(-(l1k1*l1k2*l2k1*M2) + 6*l1k1*l1k2*m2*M2 -
			           3*l1k2*m4*M2 + l2k2*m4*M2 + l1k1*l1k2*l2p2*p1k1 -
			           l2p2*m4*p1k2 - l1k1*l1k2*p1k1*p1l2 + m4*p1k2*p1l2 +
			           k1p2*l1k1*l1k2*(-l2p2 + p1l2) + l1k1*l1k2*l2k1*p1p2 -
			           6*l1k1*l1k2*m2*p1p2 + 3*l1k2*m4*p1p2 - l2k2*m4*p1p2 +
			           k1k2*(3*m4*(M2 - p1p2) +
			              l1k1*(l1p2*l2p2 + l1l2*M2 - 9*m2*M2 - l2p2*p1l1 -
			                 l1p2*p1l2 + p1l1*p1l2 - l1l2*p1p2 + 9*m2*p1p2)) +
			           l2p2*m4*p2k2 - m4*p1l2*p2k2 + l2k2*M2*pow(l1k1,2) -
			           l2p2*p1k2*pow(l1k1,2) + p1k2*p1l2*pow(l1k1,2) -
			           l2k2*p1p2*pow(l1k1,2) + l2p2*p2k2*pow(l1k1,2) -
			           p1l2*p2k2*pow(l1k1,2)) +
			        4*M2*pow(f1,2)*(l1k1*l1k2*
			            (-(l2k1*M2) + l2p2*p1k1 + k1p2*p1l2 + m2*(4*M2 - 2*p1p2)) \
			+ k1k2*(m4*(2*M2 - p1p2) + l1k1*
			               (l1l2*M2 - l2p2*p1l1 - l1p2*p1l2 + m2*(-6*M2 + 3*p1p2))\
			) + m4*(l2k2*M2 - l2p2*p1k2 + l1k2*(-2*M2 + p1p2) - p1l2*p2k2) +
			           (l2k2*M2 - l2p2*p1k2 - p1l2*p2k2)*pow(l1k1,2)) -
			        pow(f2,2)*(-(l1k1*l1k2*l2p2*M2*p1k1) + l2p2*m4*M2*p1k2 +
			           3*l1k1*l1k2*M2*p1k1*p1l2 - 3*m4*M2*p1k2*p1l2 -
			           2*l1k1*l1k2*l2k1*M2*p1p2 + 12*l1k1*l1k2*m2*M2*p1p2 -
			           6*l1k2*m4*M2*p1p2 + 2*l2k2*m4*M2*p1p2 -
			           l1k1*l1k2*l2p2*p1k1*p1p2 + l2p2*m4*p1k2*p1p2 -
			           l1k1*l1k2*p1k1*p1l2*p1p2 + m4*p1k2*p1l2*p1p2 +
			           k1p2*l1k1*l1k2*(l2p2*(3*M2 - p1p2) - p1l2*(M2 + p1p2)) -
			           3*l2p2*m4*M2*p2k2 + m4*M2*p1l2*p2k2 + l2p2*m4*p1p2*p2k2 +
			           m4*p1l2*p1p2*p2k2 + l2p2*M2*p1k2*pow(l1k1,2) -
			           3*M2*p1k2*p1l2*pow(l1k1,2) + 2*l2k2*M2*p1p2*pow(l1k1,2) +
			           l2p2*p1k2*p1p2*pow(l1k1,2) + p1k2*p1l2*p1p2*pow(l1k1,2) -
			           3*l2p2*M2*p2k2*pow(l1k1,2) + M2*p1l2*p2k2*pow(l1k1,2) +
			           l2p2*p1p2*p2k2*pow(l1k1,2) + p1l2*p1p2*p2k2*pow(l1k1,2) +
			           l1k1*l1k2*l2k1*M4 - 10*l1k1*l1k2*m2*M4 +
			           5*l1k2*m4*M4 - l2k2*m4*M4 -
			           l2k2*pow(l1k1,2)*M4 + l1k1*l1k2*l2k1*pow(p1p2,2) -
			           2*l1k1*l1k2*m2*pow(p1p2,2) + l1k2*m4*pow(p1p2,2) -
			           l2k2*m4*pow(p1p2,2) - l2k2*pow(l1k1,2)*pow(p1p2,2) +
			           k1k2*(-(m4*(-6*M2*p1p2 + 5*M4 + pow(p1p2,2))) +
			              l1k1*(l2p2*M2*p1l1 - 3*M2*p1l1*p1l2 - 18*m2*M2*p1p2 +
			                 l2p2*p1l1*p1p2 + p1l1*p1l2*p1p2 +
			                 l1p2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) +
			                 15*m2*M4 - l1l2*pow(M2 - p1p2,2) +
			                 3*m2*pow(p1p2,2)))))) +
			  8*eg1*pow(l1k1,-1)*pow(l2k2,-1)*pow(l2k1 + eg1*(k1k2 + l2k2),-1)*
			   pow(M,-2)*(4*M2*pow(f1,2)*
			      (-2*l1k2*l1l2*M2 + 2*l1k2*l2k1*M2 + 2*l1k2*m2*M2 +
			        l2k2*m2*M2 + l1k2*l1p2*p1k1 + l1p2*l2k2*p1k1 - k1p2*l1l2*p1k2 +
			        l1l2*l1p2*p1k2 + k1p2*m2*p1k2 + k1p2*l1k2*p1l1 -
			        2*k1k2*l1p2*p1l1 + k1p2*l2k2*p1l1 - 2*l1p2*l2k2*p1l1 -
			        k1k2*l2p2*p1l1 + l1k2*l2p2*p1l1 - k1k2*l1p2*p1l2 + l1k2*l1p2*p1l2 +
			        k1k2*l1l2*p1p2 - l1k2*l2k1*p1p2 - l1k2*m2*p1p2 - l1l2*p1k1*p2k2 +
			        m2*p1k1*p2k2 + l1l2*p1l1*p2k2 +
			        l1k1*(-2*l1k2*M2 + l1p2*p1k2 + l2p2*p1k2 - l2k2*p1p2 +
			           p1l1*p2k2 + p1l2*p2k2) -
			        eg1*(-2*l1k1*p1k2*p2k2 +
			           l1k2*(k1p2*p1k2 - 2*l1p2*p1k2 + p1k1*p2k2 - 2*p1l1*p2k2) +
			           k1k2*(-2*l1k2*M2 + l1p2*p1k2 + p1l1*p2k2) + 2*M2*pow(l1k2,2))\
			) - 4*f1*f2*M2*(l1k2*l1p2*l2p2 + 2*l1k1*l1k2*M2 + 2*l1k2*l1l2*M2 -
			        3*l1k2*l2k1*M2 - l1k1*l2k2*M2 - 3*l1k2*m2*M2 -
			        l2k2*m2*M2 - l1k2*l1p2*p1k1 - l1p2*l2k2*p1k1 - l1k1*l1p2*p1k2 -
			        2*eg1*l1k2*l1p2*p1k2 - l1l2*l1p2*p1k2 - l1k1*l2p2*p1k2 -
			        eg1*l1k2*p1k1*p1k2 - l1l2*p1k1*p1k2 + m2*p1k1*p1k2 +
			        2*l1p2*l2k2*p1l1 - l1k2*l2p2*p1l1 + l1k2*p1k1*p1l1 +
			        l2k2*p1k1*p1l1 + l1k1*p1k2*p1l1 + 2*eg1*l1k2*p1k2*p1l1 +
			        l1l2*p1k2*p1l1 - l1k2*l1p2*p1l2 + l1k1*p1k2*p1l2 + l1k2*p1l1*p1l2 -
			        2*l1k1*l1k2*p1p2 - 2*l1k2*l1l2*p1p2 + 3*l1k2*l2k1*p1p2 +
			        l1k1*l2k2*p1p2 + 3*l1k2*m2*p1p2 + l2k2*m2*p1p2 + l1k1*l1p2*p2k2 +
			        2*eg1*l1k2*l1p2*p2k2 + l1l2*l1p2*p2k2 + l1k1*l2p2*p2k2 +
			        eg1*l1k2*p1k1*p2k2 + l1l2*p1k1*p2k2 - m2*p1k1*p2k2 -
			        2*eg1*l1k1*p1k2*p2k2 - l1k1*p1l1*p2k2 - 2*eg1*l1k2*p1l1*p2k2 -
			        l1l2*p1l1*p2k2 - l1k1*p1l2*p2k2 +
			        k1p2*(l1p2*l2k2 + l1l2*p1k2 - m2*p1k2 - l2k2*p1l1 - l1l2*p2k2 +
			           m2*p2k2 + l1k2*(l1p2 + eg1*p1k2 - p1l1 - eg1*p2k2)) +
			        2*eg1*M2*pow(l1k2,2) - 2*eg1*p1p2*pow(l1k2,2) -
			        l2k2*pow(l1p2,2) + eg1*l1k1*pow(p1k2,2) - l2k2*pow(p1l1,2) -
			        k1k2*(-(l1l2*M2) - l2p2*p1l1 + p1l1*p1l2 + l1l2*p1p2 +
			           l1p2*(l2p2 - eg1*p1k2 - 2*p1l1 - p1l2 + eg1*p2k2) +
			           eg1*(2*l1k2*M2 + p1k2*p1l1 - 2*l1k2*p1p2 - p1l1*p2k2) +
			           pow(l1p2,2) + pow(p1l1,2)) + eg1*l1k1*pow(p2k2,2)) +
			     pow(f2,2)*(-3*l1k2*l1p2*l2p2*M2 + l1k2*l1p2*M2*p1k1 +
			        l1p2*l2k2*M2*p1k1 + l1k1*l1p2*M2*p1k2 +
			        2*eg1*l1k2*l1p2*M2*p1k2 + l1l2*l1p2*M2*p1k2 +
			        l1k1*l2p2*M2*p1k2 + 3*eg1*l1k2*M2*p1k1*p1k2 +
			        3*l1l2*M2*p1k1*p1k2 - 3*m2*M2*p1k1*p1k2 -
			        2*l1p2*l2k2*M2*p1l1 + l1k2*l2p2*M2*p1l1 -
			        3*l1k2*M2*p1k1*p1l1 - 3*l2k2*M2*p1k1*p1l1 -
			        3*l1k1*M2*p1k2*p1l1 - 6*eg1*l1k2*M2*p1k2*p1l1 -
			        3*l1l2*M2*p1k2*p1l1 + l1k2*l1p2*M2*p1l2 -
			        3*l1k1*M2*p1k2*p1l2 - 3*l1k2*M2*p1l1*p1l2 +
			        l1k2*l1p2*l2p2*p1p2 + 4*l1k1*l1k2*M2*p1p2 +
			        4*l1k2*l1l2*M2*p1p2 - 6*l1k2*l2k1*M2*p1p2 -
			        2*l1k1*l2k2*M2*p1p2 - 6*l1k2*m2*M2*p1p2 -
			        2*l2k2*m2*M2*p1p2 + l1k2*l1p2*p1k1*p1p2 + l1p2*l2k2*p1k1*p1p2 +
			        l1k1*l1p2*p1k2*p1p2 + 2*eg1*l1k2*l1p2*p1k2*p1p2 +
			        l1l2*l1p2*p1k2*p1p2 + l1k1*l2p2*p1k2*p1p2 -
			        eg1*l1k2*p1k1*p1k2*p1p2 - l1l2*p1k1*p1k2*p1p2 + m2*p1k1*p1k2*p1p2 -
			        2*l1p2*l2k2*p1l1*p1p2 + l1k2*l2p2*p1l1*p1p2 + l1k2*p1k1*p1l1*p1p2 +
			        l2k2*p1k1*p1l1*p1p2 + l1k1*p1k2*p1l1*p1p2 +
			        2*eg1*l1k2*p1k2*p1l1*p1p2 + l1l2*p1k2*p1l1*p1p2 +
			        l1k2*l1p2*p1l2*p1p2 + l1k1*p1k2*p1l2*p1p2 + l1k2*p1l1*p1l2*p1p2 -
			        3*l1k1*l1p2*M2*p2k2 - 6*eg1*l1k2*l1p2*M2*p2k2 -
			        3*l1l2*l1p2*M2*p2k2 - 3*l1k1*l2p2*M2*p2k2 -
			        eg1*l1k2*M2*p1k1*p2k2 - l1l2*M2*p1k1*p2k2 + m2*M2*p1k1*p2k2 +
			        2*eg1*l1k1*M2*p1k2*p2k2 + l1k1*M2*p1l1*p2k2 +
			        2*eg1*l1k2*M2*p1l1*p2k2 + l1l2*M2*p1l1*p2k2 +
			        l1k1*M2*p1l2*p2k2 + l1k1*l1p2*p1p2*p2k2 +
			        2*eg1*l1k2*l1p2*p1p2*p2k2 + l1l2*l1p2*p1p2*p2k2 +
			        l1k1*l2p2*p1p2*p2k2 - eg1*l1k2*p1k1*p1p2*p2k2 -
			        l1l2*p1k1*p1p2*p2k2 + m2*p1k1*p1p2*p2k2 +
			        2*eg1*l1k1*p1k2*p1p2*p2k2 + l1k1*p1l1*p1p2*p2k2 +
			        2*eg1*l1k2*p1l1*p1p2*p2k2 + l1l2*p1l1*p1p2*p2k2 +
			        l1k1*p1l2*p1p2*p2k2 + k1p2*
			         (-3*l1p2*l2k2*M2 - l1l2*M2*p1k2 + m2*M2*p1k2 +
			           l2k2*M2*p1l1 + l1p2*l2k2*p1p2 - l1l2*p1k2*p1p2 +
			           m2*p1k2*p1p2 + l2k2*p1l1*p1p2 + 3*l1l2*M2*p2k2 -
			           3*m2*M2*p2k2 - l1l2*p1p2*p2k2 + m2*p1p2*p2k2 +
			           l1k2*(l1p2*(-3*M2 + p1p2) + p1l1*(M2 + p1p2) -
			              eg1*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)))) +
			        4*eg1*M2*p1p2*pow(l1k2,2) + 3*l2k2*M2*pow(l1p2,2) -
			        l2k2*p1p2*pow(l1p2,2) - 2*l1k1*l1k2*M4 -
			        2*l1k2*l1l2*M4 + 5*l1k2*l2k1*M4 +
			        3*l1k1*l2k2*M4 + 5*l1k2*m2*M4 +
			        l2k2*m2*M4 - 2*eg1*pow(l1k2,2)*M4 -
			        3*eg1*l1k1*M2*pow(p1k2,2) + eg1*l1k1*p1p2*pow(p1k2,2) +
			        3*l2k2*M2*pow(p1l1,2) - l2k2*p1p2*pow(p1l1,2) -
			        2*l1k1*l1k2*pow(p1p2,2) - 2*l1k2*l1l2*pow(p1p2,2) +
			        l1k2*l2k1*pow(p1p2,2) - l1k1*l2k2*pow(p1p2,2) +
			        l1k2*m2*pow(p1p2,2) + l2k2*m2*pow(p1p2,2) -
			        2*eg1*pow(l1k2,2)*pow(p1p2,2) +
			        k1k2*(-(l2p2*M2*p1l1) + 3*eg1*M2*p1k2*p1l1 + 3*M2*p1l1*p1l2 -
			           4*eg1*l1k2*M2*p1p2 + 2*l1l2*M2*p1p2 - l2p2*p1l1*p1p2 -
			           eg1*p1k2*p1l1*p1p2 - p1l1*p1l2*p1p2 - eg1*M2*p1l1*p2k2 -
			           eg1*p1l1*p1p2*p2k2 +
			           l1p2*(l2p2*(3*M2 - p1p2) - (2*p1l1 + p1l2)*(M2 + p1p2) -
			              eg1*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2))) +
			           (3*M2 - p1p2)*pow(l1p2,2) + 2*eg1*l1k2*M4 -
			           3*l1l2*M4 + 3*M2*pow(p1l1,2) - p1p2*pow(p1l1,2) +
			           2*eg1*l1k2*pow(p1p2,2) + l1l2*pow(p1p2,2)) -
			        3*eg1*l1k1*M2*pow(p2k2,2) + eg1*l1k1*p1p2*pow(p2k2,2))) +
			  4*eg1*m2*pow(l1k1,-1)*pow(l1k2,-1)*pow(l2k1,-1)*pow(l2k2,-1)/M2*
			   (-4*M2*pow(f1,2)*(k1p2*(-l1k2 + l2k2)*p1k2 +
			        (-(l1k2*p1k1) + l2k2*p1k1 - 2*(l1k1 + 2*l1l2 - l2k1)*p1k2)*p2k2 -
			        k1k2*(2*k1p2*p1k2 - l1p2*p1k2 + l2p2*p1k2 + 2*p1k1*p2k2 -
			           p1l1*p2k2 + p1l2*p2k2) + 2*M2*pow(k1k2,2)) -
			     4*f1*f2*M2*(-(k1k2*(2*k1p2 - l1p2 + l2p2 - 2*p1k1 + p1l1 - p1l2)*
			           (p1k2 - p2k2)) + (p1k2 - p2k2)*
			         (k1p2*(-l1k2 + l2k2) + l1k2*p1k1 - l2k2*p1k1 + l1k1*p1k2 +
			           2*l1l2*p1k2 - l2k1*p1k2 - l1k1*p2k2 - 2*l1l2*p2k2 + l2k1*p2k2) +
			        2*(M2 - p1p2)*pow(k1k2,2)) +
			     pow(f2,2)*(-3*l1k2*M2*p1k1*p1k2 + 3*l2k2*M2*p1k1*p1k2 +
			        l1k2*p1k1*p1k2*p1p2 - l2k2*p1k1*p1k2*p1p2 + l1k2*M2*p1k1*p2k2 -
			        l2k2*M2*p1k1*p2k2 + 2*l1k1*M2*p1k2*p2k2 +
			        4*l1l2*M2*p1k2*p2k2 - 2*l2k1*M2*p1k2*p2k2 +
			        l1k2*p1k1*p1p2*p2k2 - l2k2*p1k1*p1p2*p2k2 + 2*l1k1*p1k2*p1p2*p2k2 +
			        4*l1l2*p1k2*p1p2*p2k2 - 2*l2k1*p1k2*p1p2*p2k2 +
			        k1p2*(l1k2 - l2k2)*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)) +
			        k1k2*(l2p2*M2*p1k2 - 6*M2*p1k1*p1k2 + 3*M2*p1k2*p1l1 -
			           3*M2*p1k2*p1l2 + l2p2*p1k2*p1p2 + 2*p1k1*p1k2*p1p2 -
			           p1k2*p1l1*p1p2 + p1k2*p1l2*p1p2 - 3*l2p2*M2*p2k2 +
			           2*M2*p1k1*p2k2 - M2*p1l1*p2k2 + M2*p1l2*p2k2 +
			           l2p2*p1p2*p2k2 + 2*p1k1*p1p2*p2k2 - p1l1*p1p2*p2k2 +
			           p1l2*p1p2*p2k2 + 2*k1p2*
			            (M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)) -
			           l1p2*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2))) -
			        3*l1k1*M2*pow(p1k2,2) - 6*l1l2*M2*pow(p1k2,2) +
			        3*l2k1*M2*pow(p1k2,2) + l1k1*p1p2*pow(p1k2,2) +
			        2*l1l2*p1p2*pow(p1k2,2) - l2k1*p1p2*pow(p1k2,2) -
			        2*pow(k1k2,2)*pow(M2 - p1p2,2) - 3*l1k1*M2*pow(p2k2,2) -
			        6*l1l2*M2*pow(p2k2,2) + 3*l2k1*M2*pow(p2k2,2) +
			        l1k1*p1p2*pow(p2k2,2) + 2*l1l2*p1p2*pow(p2k2,2) -
			        l2k1*p1p2*pow(p2k2,2))) +
			  8*eg1*pow(l1k2,-1)*pow(l1k1 + eg1*(-k1k2 + l1k2),-1)*pow(l2k1,-1)*
			   pow(M,-2)*(4*M2*pow(f1,2)*
			      (2*l1k1*l2k2*M2 + 2*l1l2*l2k2*M2 - 2*l2k1*l2k2*M2 -
			        l1k2*m2*M2 - 2*l2k2*m2*M2 + l1k2*l2p2*p1k1 + l2k2*l2p2*p1k1 -
			        k1p2*l1l2*p1k2 + l1p2*l2k1*p1k2 - l1l2*l2p2*p1k2 + l2k1*l2p2*p1k2 +
			        k1p2*m2*p1k2 - k1k2*l2p2*p1l1 - l2k2*l2p2*p1l1 + k1p2*l1k2*p1l2 -
			        k1k2*l1p2*p1l2 + k1p2*l2k2*p1l2 - l1p2*l2k2*p1l2 -
			        2*k1k2*l2p2*p1l2 + 2*l1k2*l2p2*p1l2 + k1k2*l1l2*p1p2 -
			        l1k2*l2k1*p1p2 - l1k1*l2k2*p1p2 + l2k2*m2*p1p2 - l1l2*p1k1*p2k2 +
			        m2*p1k1*p2k2 + l2k1*p1l1*p2k2 - l1l2*p1l2*p2k2 + l2k1*p1l2*p2k2 +
			        eg1*(-2*l2k1*p1k2*p2k2 +
			           k1k2*(-2*l2k2*M2 + l2p2*p1k2 + p1l2*p2k2) +
			           l2k2*(k1p2*p1k2 + 2*l2p2*p1k2 + p1k1*p2k2 + 2*p1l2*p2k2) -
			           2*M2*pow(l2k2,2))) +
			     4*f1*f2*M2*(l1p2*l2k2*l2p2 + l1k2*l2k1*M2 + 3*l1k1*l2k2*M2 +
			        2*l1l2*l2k2*M2 - 2*l2k1*l2k2*M2 - l1k2*m2*M2 -
			        3*l2k2*m2*M2 + l1k2*l2p2*p1k1 + l2k2*l2p2*p1k1 +
			        l1p2*l2k1*p1k2 - l1l2*l2p2*p1k2 + l2k1*l2p2*p1k2 +
			        2*eg1*l2k2*l2p2*p1k2 + l1l2*p1k1*p1k2 - eg1*l2k2*p1k1*p1k2 -
			        m2*p1k1*p1k2 - l2k2*l2p2*p1l1 - l2k1*p1k2*p1l1 - l1p2*l2k2*p1l2 +
			        2*l1k2*l2p2*p1l2 - l1k2*p1k1*p1l2 - l2k2*p1k1*p1l2 +
			        l1l2*p1k2*p1l2 - l2k1*p1k2*p1l2 - 2*eg1*l2k2*p1k2*p1l2 +
			        l2k2*p1l1*p1l2 - l1k2*l2k1*p1p2 - 3*l1k1*l2k2*p1p2 -
			        2*l1l2*l2k2*p1p2 + 2*l2k1*l2k2*p1p2 + l1k2*m2*p1p2 +
			        3*l2k2*m2*p1p2 - l1p2*l2k1*p2k2 + l1l2*l2p2*p2k2 -
			        l2k1*l2p2*p2k2 - 2*eg1*l2k2*l2p2*p2k2 - l1l2*p1k1*p2k2 +
			        eg1*l2k2*p1k1*p2k2 + m2*p1k1*p2k2 - 2*eg1*l2k1*p1k2*p2k2 +
			        l2k1*p1l1*p2k2 - l1l2*p1l2*p2k2 + l2k1*p1l2*p2k2 +
			        2*eg1*l2k2*p1l2*p2k2 +
			        k1p2*(l1k2*(-l2p2 + p1l2) - (l1l2 - m2)*(p1k2 - p2k2) +
			           l2k2*(-l2p2 + eg1*p1k2 + p1l2 - eg1*p2k2)) -
			        2*eg1*M2*pow(l2k2,2) + 2*eg1*p1p2*pow(l2k2,2) -
			        l1k2*pow(l2p2,2) + eg1*l2k1*pow(p1k2,2) - l1k2*pow(p1l2,2) +
			        k1k2*(-(l1l2*M2) - 2*eg1*l2k2*M2 + l1p2*(l2p2 - p1l2) -
			           eg1*p1k2*p1l2 + p1l1*p1l2 + l1l2*p1p2 + 2*eg1*l2k2*p1p2 +
			           eg1*p1l2*p2k2 - l2p2*(p1l1 + 2*p1l2 + eg1*(-p1k2 + p2k2)) +
			           pow(l2p2,2) + pow(p1l2,2)) + eg1*l2k1*pow(p2k2,2)) +
			     pow(f2,2)*(3*k1k2*l1p2*l2p2*M2 + 3*l1p2*l2k2*l2p2*M2 +
			        l1k2*l2p2*M2*p1k1 + l2k2*l2p2*M2*p1k1 + l1p2*l2k1*M2*p1k2 +
			        eg1*k1k2*l2p2*M2*p1k2 - l1l2*l2p2*M2*p1k2 +
			        l2k1*l2p2*M2*p1k2 + 2*eg1*l2k2*l2p2*M2*p1k2 +
			        3*l1l2*M2*p1k1*p1k2 - 3*eg1*l2k2*M2*p1k1*p1k2 -
			        3*m2*M2*p1k1*p1k2 - k1k2*l2p2*M2*p1l1 - l2k2*l2p2*M2*p1l1 -
			        3*l2k1*M2*p1k2*p1l1 - k1k2*l1p2*M2*p1l2 - l1p2*l2k2*M2*p1l2 -
			        2*k1k2*l2p2*M2*p1l2 + 2*l1k2*l2p2*M2*p1l2 -
			        3*l1k2*M2*p1k1*p1l2 - 3*l2k2*M2*p1k1*p1l2 -
			        3*eg1*k1k2*M2*p1k2*p1l2 + 3*l1l2*M2*p1k2*p1l2 -
			        3*l2k1*M2*p1k2*p1l2 - 6*eg1*l2k2*M2*p1k2*p1l2 +
			        3*k1k2*M2*p1l1*p1l2 + 3*l2k2*M2*p1l1*p1l2 -
			        k1k2*l1p2*l2p2*p1p2 - l1p2*l2k2*l2p2*p1p2 + 2*k1k2*l1l2*M2*p1p2 -
			        2*l1k2*l2k1*M2*p1p2 + 4*eg1*k1k2*l2k2*M2*p1p2 -
			        6*l1k1*l2k2*M2*p1p2 - 4*l1l2*l2k2*M2*p1p2 +
			        4*l2k1*l2k2*M2*p1p2 + 2*l1k2*m2*M2*p1p2 +
			        6*l2k2*m2*M2*p1p2 + l1k2*l2p2*p1k1*p1p2 + l2k2*l2p2*p1k1*p1p2 +
			        l1p2*l2k1*p1k2*p1p2 + eg1*k1k2*l2p2*p1k2*p1p2 -
			        l1l2*l2p2*p1k2*p1p2 + l2k1*l2p2*p1k2*p1p2 +
			        2*eg1*l2k2*l2p2*p1k2*p1p2 - l1l2*p1k1*p1k2*p1p2 +
			        eg1*l2k2*p1k1*p1k2*p1p2 + m2*p1k1*p1k2*p1p2 - k1k2*l2p2*p1l1*p1p2 -
			        l2k2*l2p2*p1l1*p1p2 + l2k1*p1k2*p1l1*p1p2 - k1k2*l1p2*p1l2*p1p2 -
			        l1p2*l2k2*p1l2*p1p2 - 2*k1k2*l2p2*p1l2*p1p2 +
			        2*l1k2*l2p2*p1l2*p1p2 + l1k2*p1k1*p1l2*p1p2 + l2k2*p1k1*p1l2*p1p2 +
			        eg1*k1k2*p1k2*p1l2*p1p2 - l1l2*p1k2*p1l2*p1p2 +
			        l2k1*p1k2*p1l2*p1p2 + 2*eg1*l2k2*p1k2*p1l2*p1p2 -
			        k1k2*p1l1*p1l2*p1p2 - l2k2*p1l1*p1l2*p1p2 - 3*l1p2*l2k1*M2*p2k2 -
			        3*eg1*k1k2*l2p2*M2*p2k2 + 3*l1l2*l2p2*M2*p2k2 -
			        3*l2k1*l2p2*M2*p2k2 - 6*eg1*l2k2*l2p2*M2*p2k2 -
			        l1l2*M2*p1k1*p2k2 + eg1*l2k2*M2*p1k1*p2k2 + m2*M2*p1k1*p2k2 -
			        2*eg1*l2k1*M2*p1k2*p2k2 + l2k1*M2*p1l1*p2k2 +
			        eg1*k1k2*M2*p1l2*p2k2 - l1l2*M2*p1l2*p2k2 +
			        l2k1*M2*p1l2*p2k2 + 2*eg1*l2k2*M2*p1l2*p2k2 +
			        l1p2*l2k1*p1p2*p2k2 + eg1*k1k2*l2p2*p1p2*p2k2 -
			        l1l2*l2p2*p1p2*p2k2 + l2k1*l2p2*p1p2*p2k2 +
			        2*eg1*l2k2*l2p2*p1p2*p2k2 - l1l2*p1k1*p1p2*p2k2 +
			        eg1*l2k2*p1k1*p1p2*p2k2 + m2*p1k1*p1p2*p2k2 -
			        2*eg1*l2k1*p1k2*p1p2*p2k2 + l2k1*p1l1*p1p2*p2k2 +
			        eg1*k1k2*p1l2*p1p2*p2k2 - l1l2*p1l2*p1p2*p2k2 +
			        l2k1*p1l2*p1p2*p2k2 + 2*eg1*l2k2*p1l2*p1p2*p2k2 +
			        k1p2*(l1k2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) -
			           (l1l2 - m2)*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)) +
			           l2k2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2) +
			              eg1*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)))) +
			        4*eg1*M2*p1p2*pow(l2k2,2) + 3*k1k2*M2*pow(l2p2,2) -
			        3*l1k2*M2*pow(l2p2,2) - k1k2*p1p2*pow(l2p2,2) +
			        l1k2*p1p2*pow(l2p2,2) - 3*k1k2*l1l2*M4 +
			        3*l1k2*l2k1*M4 - 2*eg1*k1k2*l2k2*M4 +
			        5*l1k1*l2k2*M4 + 2*l1l2*l2k2*M4 -
			        2*l2k1*l2k2*M4 - l1k2*m2*M4 -
			        5*l2k2*m2*M4 - 2*eg1*pow(l2k2,2)*M4 +
			        3*eg1*l2k1*M2*pow(p1k2,2) - eg1*l2k1*p1p2*pow(p1k2,2) +
			        3*k1k2*M2*pow(p1l2,2) - 3*l1k2*M2*pow(p1l2,2) -
			        k1k2*p1p2*pow(p1l2,2) + l1k2*p1p2*pow(p1l2,2) +
			        k1k2*l1l2*pow(p1p2,2) - l1k2*l2k1*pow(p1p2,2) -
			        2*eg1*k1k2*l2k2*pow(p1p2,2) + l1k1*l2k2*pow(p1p2,2) +
			        2*l1l2*l2k2*pow(p1p2,2) - 2*l2k1*l2k2*pow(p1p2,2) -
			        l1k2*m2*pow(p1p2,2) - l2k2*m2*pow(p1p2,2) -
			        2*eg1*pow(l2k2,2)*pow(p1p2,2) + 3*eg1*l2k1*M2*pow(p2k2,2) -
			        eg1*l2k1*p1p2*pow(p2k2,2))) +
			  8*eg1*m2*pow(l1k2,-1)*pow(l1k1 + eg1*(-k1k2 + l1k2),-1)*pow(l2k2,-1)*
			   pow(l2k1 + eg1*(k1k2 + l2k2),-1)/M2*
			   (4*M2*pow(f1,2)*(-(k1p2*l2k2*p1k2) + l1p2*l2k2*p1k2 - l2k2*l2p2*p1k2 -
			        l2k2*p1k1*p2k2 - 2*l1k1*p1k2*p2k2 - 2*l1l2*p1k2*p2k2 +
			        2*l2k1*p1k2*p2k2 + 2*m2*p1k2*p2k2 + l2k2*p1l1*p2k2 -
			        l2k2*p1l2*p2k2 + l1k2*(-2*l2k2*M2 + k1p2*p1k2 - l1p2*p1k2 +
			           l2p2*p1k2 + p1k1*p2k2 - p1l1*p2k2 + p1l2*p2k2) -
			        k1k2*(2*l1k2*M2 - 2*l2k2*M2 + k1p2*p1k2 - l1p2*p1k2 +
			           l2p2*p1k2 + p1k1*p2k2 - p1l1*p2k2 + p1l2*p2k2) +
			        M2*pow(k1k2,2) + M2*pow(l1k2,2) + M2*pow(l2k2,2)) +
			     4*f1*f2*M2*(-(k1p2*l2k2*p1k2) + l1p2*l2k2*p1k2 - l2k2*l2p2*p1k2 +
			        l2k2*p1k1*p1k2 - l2k2*p1k2*p1l1 + l2k2*p1k2*p1l2 +
			        k1k2*(-2*l1k2*(M2 - p1p2) + 2*l2k2*(M2 - p1p2) -
			           (k1p2 - l1p2 + l2p2 - p1k1 + p1l1 - p1l2)*(p1k2 - p2k2)) +
			        l1k2*(-2*l2k2*(M2 - p1p2) +
			           (k1p2 - l1p2 + l2p2 - p1k1 + p1l1 - p1l2)*(p1k2 - p2k2)) +
			        k1p2*l2k2*p2k2 - l1p2*l2k2*p2k2 + l2k2*l2p2*p2k2 - l2k2*p1k1*p2k2 -
			        2*l1k1*p1k2*p2k2 - 2*l1l2*p1k2*p2k2 + 2*l2k1*p1k2*p2k2 +
			        2*m2*p1k2*p2k2 + l2k2*p1l1*p2k2 - l2k2*p1l2*p2k2 +
			        (M2 - p1p2)*pow(k1k2,2) + (M2 - p1p2)*pow(l1k2,2) +
			        M2*pow(l2k2,2) - p1p2*pow(l2k2,2) + l1k1*pow(p1k2,2) +
			        l1l2*pow(p1k2,2) - l2k1*pow(p1k2,2) - m2*pow(p1k2,2) +
			        l1k1*pow(p2k2,2) + l1l2*pow(p2k2,2) - l2k1*pow(p2k2,2) -
			        m2*pow(p2k2,2)) + pow(f2,2)*
			      (-(k1p2*l2k2*M2*p1k2) + l1p2*l2k2*M2*p1k2 - l2k2*l2p2*M2*p1k2 +
			        3*l2k2*M2*p1k1*p1k2 - 3*l2k2*M2*p1k2*p1l1 +
			        3*l2k2*M2*p1k2*p1l2 - k1p2*l2k2*p1k2*p1p2 + l1p2*l2k2*p1k2*p1p2 -
			        l2k2*l2p2*p1k2*p1p2 - l2k2*p1k1*p1k2*p1p2 + l2k2*p1k2*p1l1*p1p2 -
			        l2k2*p1k2*p1l2*p1p2 + 3*k1p2*l2k2*M2*p2k2 -
			        3*l1p2*l2k2*M2*p2k2 + 3*l2k2*l2p2*M2*p2k2 -
			        l2k2*M2*p1k1*p2k2 - 2*l1k1*M2*p1k2*p2k2 -
			        2*l1l2*M2*p1k2*p2k2 + 2*l2k1*M2*p1k2*p2k2 +
			        2*m2*M2*p1k2*p2k2 + l2k2*M2*p1l1*p2k2 - l2k2*M2*p1l2*p2k2 -
			        k1p2*l2k2*p1p2*p2k2 + l1p2*l2k2*p1p2*p2k2 - l2k2*l2p2*p1p2*p2k2 -
			        l2k2*p1k1*p1p2*p2k2 - 2*l1k1*p1k2*p1p2*p2k2 -
			        2*l1l2*p1k2*p1p2*p2k2 + 2*l2k1*p1k2*p1p2*p2k2 +
			        2*m2*p1k2*p1p2*p2k2 + l2k2*p1l1*p1p2*p2k2 - l2k2*p1l2*p1p2*p2k2 -
			        2*M2*p1p2*pow(l2k2,2) + pow(l2k2,2)*M4 +
			        3*l1k1*M2*pow(p1k2,2) + 3*l1l2*M2*pow(p1k2,2) -
			        3*l2k1*M2*pow(p1k2,2) - 3*m2*M2*pow(p1k2,2) -
			        l1k1*p1p2*pow(p1k2,2) - l1l2*p1p2*pow(p1k2,2) +
			        l2k1*p1p2*pow(p1k2,2) + m2*p1p2*pow(p1k2,2) +
			        pow(k1k2,2)*pow(M2 - p1p2,2) + pow(l1k2,2)*pow(M2 - p1p2,2) +
			        l1k2*(-(l1p2*M2*p1k2) + l2p2*M2*p1k2 - 3*M2*p1k1*p1k2 +
			           3*M2*p1k2*p1l1 - 3*M2*p1k2*p1l2 - l1p2*p1k2*p1p2 +
			           l2p2*p1k2*p1p2 + p1k1*p1k2*p1p2 - p1k2*p1l1*p1p2 +
			           p1k2*p1l2*p1p2 + 3*l1p2*M2*p2k2 - 3*l2p2*M2*p2k2 +
			           M2*p1k1*p2k2 - M2*p1l1*p2k2 + M2*p1l2*p2k2 -
			           l1p2*p1p2*p2k2 + l2p2*p1p2*p2k2 + p1k1*p1p2*p2k2 -
			           p1l1*p1p2*p2k2 + p1l2*p1p2*p2k2 +
			           k1p2*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)) -
			           2*l2k2*pow(M2 - p1p2,2)) -
			        k1k2*(k1p2*M2*p1k2 - l1p2*M2*p1k2 + l2p2*M2*p1k2 -
			           3*M2*p1k1*p1k2 + 3*M2*p1k2*p1l1 - 3*M2*p1k2*p1l2 +
			           k1p2*p1k2*p1p2 - l1p2*p1k2*p1p2 + l2p2*p1k2*p1p2 +
			           p1k1*p1k2*p1p2 - p1k2*p1l1*p1p2 + p1k2*p1l2*p1p2 -
			           3*k1p2*M2*p2k2 + 3*l1p2*M2*p2k2 - 3*l2p2*M2*p2k2 +
			           M2*p1k1*p2k2 - M2*p1l1*p2k2 + M2*p1l2*p2k2 +
			           k1p2*p1p2*p2k2 - l1p2*p1p2*p2k2 + l2p2*p1p2*p2k2 +
			           p1k1*p1p2*p2k2 - p1l1*p1p2*p2k2 + p1l2*p1p2*p2k2 +
			           2*l1k2*pow(M2 - p1p2,2) - 2*l2k2*pow(M2 - p1p2,2)) +
			        pow(l2k2,2)*pow(p1p2,2) + 3*l1k1*M2*pow(p2k2,2) +
			        3*l1l2*M2*pow(p2k2,2) - 3*l2k1*M2*pow(p2k2,2) -
			        3*m2*M2*pow(p2k2,2) - l1k1*p1p2*pow(p2k2,2) -
			        l1l2*p1p2*pow(p2k2,2) + l2k1*p1p2*pow(p2k2,2) + m2*p1p2*pow(p2k2,2))\
			) - 4*eg1*m2*pow(l1k2,-1)*pow(l2k1,-1)*pow(l2k2,-1)*
			   pow(l2k1 + eg1*(k1k2 + l2k2),-1)/M2*
			   (4*f1*f2*M2*(l1k2*(4*l2k2*(M2 - p1p2) +
			           (l2p2 - p1l2)*(p1k2 - p2k2)) +
			        (p1k2 - p2k2)*(k1p2*l2k2 + l1p2*l2k2 - k1k2*l2p2 - l2k2*p1k1 +
			           l1l2*p1k2 - l2k1*p1k2 - 2*eg1*l2k2*p1k2 - m2*p1k2 - l2k2*p1l1 +
			           k1k2*p1l2 - l1l2*p2k2 + l2k1*p2k2 + 2*eg1*l2k2*p2k2 + m2*p2k2)) \
			+ 4*M2*(k1p2*l2k2*p1k2 + l1p2*l2k2*p1k2 - k1k2*l2p2*p1k2 +
			        l2k2*p1k1*p2k2 - 2*l1l2*p1k2*p2k2 + 2*l2k1*p1k2*p2k2 +
			        4*eg1*l2k2*p1k2*p2k2 + 2*m2*p1k2*p2k2 + l2k2*p1l1*p2k2 -
			        k1k2*p1l2*p2k2 + l1k2*(l2p2*p1k2 + 2*l2k2*(M2 - p1p2) +
			           p1l2*p2k2))*pow(f1,2) +
			     pow(f2,2)*(k1p2*l2k2*M2*p1k2 + l1p2*l2k2*M2*p1k2 -
			        k1k2*l2p2*M2*p1k2 - 3*l2k2*M2*p1k1*p1k2 -
			        3*l2k2*M2*p1k2*p1l1 + 3*k1k2*M2*p1k2*p1l2 +
			        k1p2*l2k2*p1k2*p1p2 + l1p2*l2k2*p1k2*p1p2 - k1k2*l2p2*p1k2*p1p2 +
			        l2k2*p1k1*p1k2*p1p2 + l2k2*p1k2*p1l1*p1p2 - k1k2*p1k2*p1l2*p1p2 -
			        3*k1p2*l2k2*M2*p2k2 - 3*l1p2*l2k2*M2*p2k2 +
			        3*k1k2*l2p2*M2*p2k2 + l2k2*M2*p1k1*p2k2 -
			        2*l1l2*M2*p1k2*p2k2 + 2*l2k1*M2*p1k2*p2k2 +
			        4*eg1*l2k2*M2*p1k2*p2k2 + 2*m2*M2*p1k2*p2k2 +
			        l2k2*M2*p1l1*p2k2 - k1k2*M2*p1l2*p2k2 + k1p2*l2k2*p1p2*p2k2 +
			        l1p2*l2k2*p1p2*p2k2 - k1k2*l2p2*p1p2*p2k2 + l2k2*p1k1*p1p2*p2k2 -
			        2*l1l2*p1k2*p1p2*p2k2 + 2*l2k1*p1k2*p1p2*p2k2 +
			        4*eg1*l2k2*p1k2*p1p2*p2k2 + 2*m2*p1k2*p1p2*p2k2 +
			        l2k2*p1l1*p1p2*p2k2 - k1k2*p1l2*p1p2*p2k2 +
			        l1k2*(l2p2*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)) +
			           p1l2*(M2*(-3*p1k2 + p2k2) + p1p2*(p1k2 + p2k2)) +
			           8*l2k2*(-(M2*p1p2) + M4)) + 3*l1l2*M2*pow(p1k2,2) -
			        3*l2k1*M2*pow(p1k2,2) - 6*eg1*l2k2*M2*pow(p1k2,2) -
			        3*m2*M2*pow(p1k2,2) - l1l2*p1p2*pow(p1k2,2) +
			        l2k1*p1p2*pow(p1k2,2) + 2*eg1*l2k2*p1p2*pow(p1k2,2) +
			        m2*p1p2*pow(p1k2,2) + 3*l1l2*M2*pow(p2k2,2) -
			        3*l2k1*M2*pow(p2k2,2) - 6*eg1*l2k2*M2*pow(p2k2,2) -
			        3*m2*M2*pow(p2k2,2) - l1l2*p1p2*pow(p2k2,2) +
			        l2k1*p1p2*pow(p2k2,2) + 2*eg1*l2k2*p1p2*pow(p2k2,2) +
			        m2*p1p2*pow(p2k2,2))) -
			  8*eg1*pow(l1k1,-2)*pow(l1k1 + eg1*(-k1k2 + l1k2),-1)*pow(l2k2,-1)*
			   pow(M,-2)*(-4*M2*pow(f1,2)*
			      (k1k2*(m2*(-2*l1l2*M2 + 2*l2k1*M2 + 2*m2*M2 - l2p2*p1k1 +
			              l2p2*p1l1 - k1p2*p1l2 + l1p2*p1l2 - m2*p1p2) +
			           l1k1*(2*l2k1*M2 - l2p2*p1k1 + l2p2*p1l1 - k1p2*p1l2 +
			              l1p2*p1l2 - l1l2*p1p2 + m2*p1p2)) +
			        m2*(k1p2*l1l2*p1k2 - l1l2*l1p2*p1k2 - k1p2*l2k1*p1k2 +
			           l1p2*l2k1*p1k2 - l2k2*
			            (m2*M2 - 2*(k1p2 - l1p2)*(p1k1 - p1l1)) +
			           l1k2*(2*l1l2*M2 - 2*l2k1*M2 - 2*m2*M2 + l2p2*p1k1 -
			              l2p2*p1l1 + k1p2*p1l2 - l1p2*p1l2 + m2*p1p2) +
			           l1l2*p1k1*p2k2 - l2k1*p1k1*p2k2 + 2*eg1*m2*p1k2*p2k2 -
			           l1l2*p1l1*p2k2 + l2k1*p1l1*p2k2) +
			        l1k1*(k1p2*l1l2*p1k2 - k1p2*l2k1*p1k2 - 2*k1p2*m2*p1k2 +
			           l1p2*m2*p1k2 - l2p2*m2*p1k2 +
			           l1k2*(m2*(M2 - p1p2) + l2k1*(-2*M2 + p1p2)) +
			           l2k2*(2*k1p2*p1k1 - l1p2*p1k1 - k1p2*p1l1 +
			              m2*(M2 + p1p2)) + l1l2*p1k1*p2k2 - l2k1*p1k1*p2k2 -
			           2*m2*p1k1*p2k2 - 2*eg1*m2*p1k2*p2k2 + m2*p1l1*p2k2 -
			           m2*p1l2*p2k2) - (l2p2*p1k2 - l2k2*p1p2 + p1l2*p2k2)*pow(l1k1,2)) \
			+ 4*f1*f2*M2*(-(l1k2*l1p2*l2p2*m2) + 3*l1k1*l1k2*l2k1*M2 -
			        2*l1k1*l1k2*m2*M2 - 2*l1k2*l1l2*m2*M2 +
			        2*l1k2*l2k1*m2*M2 + 3*l1k2*m4*M2 + l2k2*m4*M2 +
			        l1k1*l1p2*l2k2*p1k1 + 2*l1p2*l2k2*m2*p1k1 - l1k2*l2p2*m2*p1k1 -
			        l1k1*l1p2*m2*p1k2 + l1l2*l1p2*m2*p1k2 - l1p2*l2k1*m2*p1k2 +
			        l1k1*l2p2*m2*p1k2 + l1k1*l1l2*p1k1*p1k2 - l1k1*l2k1*p1k1*p1k2 -
			        2*l1k1*m2*p1k1*p1k2 + l1l2*m2*p1k1*p1k2 - l2k1*m2*p1k1*p1k2 -
			        2*l1p2*l2k2*m2*p1l1 + l1k2*l2p2*m2*p1l1 - l1k1*l2k2*p1k1*p1l1 -
			        2*l2k2*m2*p1k1*p1l1 + l1k1*m2*p1k2*p1l1 - l1l2*m2*p1k2*p1l1 +
			        l2k1*m2*p1k2*p1l1 + l1k2*l1p2*m2*p1l2 + l1k2*m2*p1k1*p1l2 -
			        l1k1*m2*p1k2*p1l2 - l1k2*m2*p1l1*p1l2 - 3*l1k1*l1k2*l2k1*p1p2 +
			        2*l1k1*l1k2*m2*p1p2 + 2*l1k2*l1l2*m2*p1p2 -
			        2*l1k2*l2k1*m2*p1p2 - 3*l1k2*m4*p1p2 - l2k2*m4*p1p2 +
			        k1k2*(l1k1*(l1p2*l2p2 - l1l2*M2 - 2*l2k1*M2 + m2*M2 +
			              l2p2*p1k1 - l2p2*p1l1 - l1p2*p1l2 - p1k1*p1l2 + p1l1*p1l2 +
			              l1l2*p1p2 + 2*l2k1*p1p2 - m2*p1p2) +
			           m2*(l1p2*l2p2 + 2*l1l2*M2 - 2*l2k1*M2 - 3*m2*M2 +
			              l2p2*p1k1 - l2p2*p1l1 - l1p2*p1l2 - p1k1*p1l2 + p1l1*p1l2 -
			              2*l1l2*p1p2 + 2*l2k1*p1p2 + 3*m2*p1p2)) +
			        l1k1*l1p2*m2*p2k2 - l1l2*l1p2*m2*p2k2 + l1p2*l2k1*m2*p2k2 -
			        l1k1*l2p2*m2*p2k2 - l1k1*l1l2*p1k1*p2k2 + l1k1*l2k1*p1k1*p2k2 +
			        2*l1k1*m2*p1k1*p2k2 - l1l2*m2*p1k1*p2k2 + l2k1*m2*p1k1*p2k2 +
			        2*eg1*l1k1*m2*p1k2*p2k2 - 2*eg1*m4*p1k2*p2k2 -
			        l1k1*m2*p1l1*p2k2 + l1l2*m2*p1l1*p2k2 - l2k1*m2*p1l1*p2k2 +
			        l1k1*m2*p1l2*p2k2 - k1p2*
			         (m2*(2*l1p2*l2k2 + k1k2*l2p2 - l1k2*l2p2 + 2*l2k2*p1k1 +
			              l1l2*p1k2 - l2k1*p1k2 - 2*l2k2*p1l1 - k1k2*p1l2 +
			              l1k2*p1l2 - l1l2*p2k2 + l2k1*p2k2) +
			           l1k1*(l1p2*l2k2 + k1k2*l2p2 + 2*l2k2*p1k1 + l1l2*p1k2 -
			              l2k1*p1k2 - 2*m2*p1k2 - l2k2*p1l1 - k1k2*p1l2 - l1l2*p2k2 +
			              l2k1*p2k2 + 2*m2*p2k2)) + l2k2*(l1k1 + m2)*pow(k1p2,2) +
			        l2k2*M2*pow(l1k1,2) + l2p2*p1k2*pow(l1k1,2) -
			        p1k2*p1l2*pow(l1k1,2) - l2k2*p1p2*pow(l1k1,2) -
			        l2p2*p2k2*pow(l1k1,2) + p1l2*p2k2*pow(l1k1,2) +
			        l2k2*m2*pow(l1p2,2) + l1k1*l2k2*pow(p1k1,2) +
			        l2k2*m2*pow(p1k1,2) - eg1*l1k1*m2*pow(p1k2,2) +
			        eg1*m4*pow(p1k2,2) + l2k2*m2*pow(p1l1,2) -
			        eg1*l1k1*m2*pow(p2k2,2) + eg1*m4*pow(p2k2,2)) +
			     pow(f2,2)*(-3*l1k2*l1p2*l2p2*m2*M2 + l1k1*l1p2*l2k2*M2*p1k1 +
			        2*l1p2*l2k2*m2*M2*p1k1 - l1k2*l2p2*m2*M2*p1k1 -
			        l1k1*l1p2*m2*M2*p1k2 + l1l2*l1p2*m2*M2*p1k2 -
			        l1p2*l2k1*m2*M2*p1k2 + l1k1*l2p2*m2*M2*p1k2 +
			        3*l1k1*l1l2*M2*p1k1*p1k2 - 3*l1k1*l2k1*M2*p1k1*p1k2 -
			        6*l1k1*m2*M2*p1k1*p1k2 + 3*l1l2*m2*M2*p1k1*p1k2 -
			        3*l2k1*m2*M2*p1k1*p1k2 - 2*l1p2*l2k2*m2*M2*p1l1 +
			        l1k2*l2p2*m2*M2*p1l1 - 3*l1k1*l2k2*M2*p1k1*p1l1 -
			        6*l2k2*m2*M2*p1k1*p1l1 + 3*l1k1*m2*M2*p1k2*p1l1 -
			        3*l1l2*m2*M2*p1k2*p1l1 + 3*l2k1*m2*M2*p1k2*p1l1 +
			        l1k2*l1p2*m2*M2*p1l2 + 3*l1k2*m2*M2*p1k1*p1l2 -
			        3*l1k1*m2*M2*p1k2*p1l2 - 3*l1k2*m2*M2*p1l1*p1l2 +
			        l1k2*l1p2*l2p2*m2*p1p2 - 6*l1k1*l1k2*l2k1*M2*p1p2 +
			        4*l1k1*l1k2*m2*M2*p1p2 + 4*l1k2*l1l2*m2*M2*p1p2 -
			        4*l1k2*l2k1*m2*M2*p1p2 - 6*l1k2*m4*M2*p1p2 -
			        2*l2k2*m4*M2*p1p2 + l1k1*l1p2*l2k2*p1k1*p1p2 +
			        2*l1p2*l2k2*m2*p1k1*p1p2 - l1k2*l2p2*m2*p1k1*p1p2 -
			        l1k1*l1p2*m2*p1k2*p1p2 + l1l2*l1p2*m2*p1k2*p1p2 -
			        l1p2*l2k1*m2*p1k2*p1p2 + l1k1*l2p2*m2*p1k2*p1p2 -
			        l1k1*l1l2*p1k1*p1k2*p1p2 + l1k1*l2k1*p1k1*p1k2*p1p2 +
			        2*l1k1*m2*p1k1*p1k2*p1p2 - l1l2*m2*p1k1*p1k2*p1p2 +
			        l2k1*m2*p1k1*p1k2*p1p2 - 2*l1p2*l2k2*m2*p1l1*p1p2 +
			        l1k2*l2p2*m2*p1l1*p1p2 + l1k1*l2k2*p1k1*p1l1*p1p2 +
			        2*l2k2*m2*p1k1*p1l1*p1p2 - l1k1*m2*p1k2*p1l1*p1p2 +
			        l1l2*m2*p1k2*p1l1*p1p2 - l2k1*m2*p1k2*p1l1*p1p2 +
			        l1k2*l1p2*m2*p1l2*p1p2 - l1k2*m2*p1k1*p1l2*p1p2 +
			        l1k1*m2*p1k2*p1l2*p1p2 + l1k2*m2*p1l1*p1l2*p1p2 +
			        3*l1k1*l1p2*m2*M2*p2k2 - 3*l1l2*l1p2*m2*M2*p2k2 +
			        3*l1p2*l2k1*m2*M2*p2k2 - 3*l1k1*l2p2*m2*M2*p2k2 -
			        l1k1*l1l2*M2*p1k1*p2k2 + l1k1*l2k1*M2*p1k1*p2k2 +
			        2*l1k1*m2*M2*p1k1*p2k2 - l1l2*m2*M2*p1k1*p2k2 +
			        l2k1*m2*M2*p1k1*p2k2 + 2*eg1*l1k1*m2*M2*p1k2*p2k2 -
			        2*eg1*m4*M2*p1k2*p2k2 - l1k1*m2*M2*p1l1*p2k2 +
			        l1l2*m2*M2*p1l1*p2k2 - l2k1*m2*M2*p1l1*p2k2 +
			        l1k1*m2*M2*p1l2*p2k2 - l1k1*l1p2*m2*p1p2*p2k2 +
			        l1l2*l1p2*m2*p1p2*p2k2 - l1p2*l2k1*m2*p1p2*p2k2 +
			        l1k1*l2p2*m2*p1p2*p2k2 - l1k1*l1l2*p1k1*p1p2*p2k2 +
			        l1k1*l2k1*p1k1*p1p2*p2k2 + 2*l1k1*m2*p1k1*p1p2*p2k2 -
			        l1l2*m2*p1k1*p1p2*p2k2 + l2k1*m2*p1k1*p1p2*p2k2 +
			        2*eg1*l1k1*m2*p1k2*p1p2*p2k2 - 2*eg1*m4*p1k2*p1p2*p2k2 -
			        l1k1*m2*p1l1*p1p2*p2k2 + l1l2*m2*p1l1*p1p2*p2k2 -
			        l2k1*m2*p1l1*p1p2*p2k2 + l1k1*m2*p1l2*p1p2*p2k2 +
			        k1p2*(m2*(3*l1k2*l2p2*M2 - 2*l2k2*M2*p1k1 - l1l2*M2*p1k2 +
			              l2k1*M2*p1k2 + 2*l2k2*M2*p1l1 - l1k2*M2*p1l2 -
			              l1k2*l2p2*p1p2 - 2*l2k2*p1k1*p1p2 - l1l2*p1k2*p1p2 +
			              l2k1*p1k2*p1p2 + 2*l2k2*p1l1*p1p2 - l1k2*p1l2*p1p2 +
			              2*l1p2*l2k2*(-3*M2 + p1p2) +
			              k1k2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) +
			              3*l1l2*M2*p2k2 - 3*l2k1*M2*p2k2 - l1l2*p1p2*p2k2 +
			              l2k1*p1p2*p2k2) +
			           l1k1*(-2*l2k2*M2*p1k1 - l1l2*M2*p1k2 + l2k1*M2*p1k2 +
			              2*m2*M2*p1k2 + l2k2*M2*p1l1 - 2*l2k2*p1k1*p1p2 -
			              l1l2*p1k2*p1p2 + l2k1*p1k2*p1p2 + 2*m2*p1k2*p1p2 +
			              l2k2*p1l1*p1p2 + l1p2*l2k2*(-3*M2 + p1p2) +
			              k1k2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) +
			              3*l1l2*M2*p2k2 - 3*l2k1*M2*p2k2 - 6*m2*M2*p2k2 -
			              l1l2*p1p2*p2k2 + l2k1*p1p2*p2k2 + 2*m2*p1p2*p2k2)) +
			        l2k2*(l1k1 + m2)*(3*M2 - p1p2)*pow(k1p2,2) +
			        l2p2*M2*p1k2*pow(l1k1,2) - 3*M2*p1k2*p1l2*pow(l1k1,2) -
			        2*l2k2*M2*p1p2*pow(l1k1,2) + l2p2*p1k2*p1p2*pow(l1k1,2) +
			        p1k2*p1l2*p1p2*pow(l1k1,2) - 3*l2p2*M2*p2k2*pow(l1k1,2) +
			        M2*p1l2*p2k2*pow(l1k1,2) + l2p2*p1p2*p2k2*pow(l1k1,2) +
			        p1l2*p1p2*p2k2*pow(l1k1,2) + 3*l2k2*m2*M2*pow(l1p2,2) -
			        l2k2*m2*p1p2*pow(l1p2,2) + 5*l1k1*l1k2*l2k1*M4 -
			        4*l1k1*l1k2*m2*M4 - 2*l1k2*l1l2*m2*M4 +
			        2*l1k2*l2k1*m2*M4 + 2*l1k1*l2k2*m2*M4 +
			        5*l1k2*m4*M4 + l2k2*m4*M4 +
			        3*l2k2*pow(l1k1,2)*M4 + 3*l1k1*l2k2*M2*pow(p1k1,2) +
			        3*l2k2*m2*M2*pow(p1k1,2) - l1k1*l2k2*p1p2*pow(p1k1,2) -
			        l2k2*m2*p1p2*pow(p1k1,2) - 3*eg1*l1k1*m2*M2*pow(p1k2,2) +
			        3*eg1*m4*M2*pow(p1k2,2) + eg1*l1k1*m2*p1p2*pow(p1k2,2) -
			        eg1*m4*p1p2*pow(p1k2,2) + 3*l2k2*m2*M2*pow(p1l1,2) -
			        l2k2*m2*p1p2*pow(p1l1,2) + l1k1*l1k2*l2k1*pow(p1p2,2) -
			        2*l1k2*l1l2*m2*pow(p1p2,2) + 2*l1k2*l2k1*m2*pow(p1p2,2) -
			        2*l1k1*l2k2*m2*pow(p1p2,2) + l1k2*m4*pow(p1p2,2) +
			        l2k2*m4*pow(p1p2,2) - l2k2*pow(l1k1,2)*pow(p1p2,2) +
			        k1k2*(-(m2*(-(l2p2*M2*p1k1) + l2p2*M2*p1l1 +
			                3*M2*p1k1*p1l2 - 3*M2*p1l1*p1l2 - 4*l2k1*M2*p1p2 -
			                6*m2*M2*p1p2 - l2p2*p1k1*p1p2 + l2p2*p1l1*p1p2 -
			                p1k1*p1l2*p1p2 + p1l1*p1l2*p1p2 +
			                l1p2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) +
			                2*l2k1*M4 + 5*m2*M4 -
			                2*l1l2*pow(M2 - p1p2,2) + 2*l2k1*pow(p1p2,2) +
			                m2*pow(p1p2,2))) -
			           l1k1*(-(l2p2*M2*p1k1) + l2p2*M2*p1l1 + 3*M2*p1k1*p1l2 -
			              3*M2*p1l1*p1l2 - 4*l2k1*M2*p1p2 + 2*m2*M2*p1p2 -
			              l2p2*p1k1*p1p2 + l2p2*p1l1*p1p2 - p1k1*p1l2*p1p2 +
			              p1l1*p1l2*p1p2 +
			              l1p2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) +
			              2*l2k1*M4 - 3*m2*M4 +
			              l1l2*(-2*M2*p1p2 + 3*M4 - pow(p1p2,2)) +
			              2*l2k1*pow(p1p2,2) + m2*pow(p1p2,2))) -
			        3*eg1*l1k1*m2*M2*pow(p2k2,2) + 3*eg1*m4*M2*pow(p2k2,2) +
			        eg1*l1k1*m2*p1p2*pow(p2k2,2) - eg1*m4*p1p2*pow(p2k2,2))) -
			  8*eg1*pow(l1k2,-1)*pow(l2k1,-2)*pow(l2k1 + eg1*(k1k2 + l2k2),-1)*
			   pow(M,-2)*(4*M2*pow(f1,2)*
			      (l1k2*l2k1*m2*M2 - 2*l1l2*l2k2*m2*M2 + l2k1*l2k2*m2*M2 +
			        l1k2*m4*M2 + 2*l2k2*m4*M2 + 2*k1p2*l1k2*l2k1*p1k1 +
			        l1k2*l2k1*l2p2*p1k1 - 2*k1p2*l1k2*m2*p1k1 + l1p2*l2k2*m2*p1k1 -
			        2*l1k2*l2p2*m2*p1k1 - k1p2*l1l2*l2k1*p1k2 + k1p2*l1l2*m2*p1k2 +
			        2*k1p2*l2k1*m2*p1k2 - l1p2*l2k1*m2*p1k2 + l1l2*l2p2*m2*p1k2 +
			        l2k1*l2p2*m2*p1k2 + k1p2*l2k2*m2*p1l1 + l2k2*l2p2*m2*p1l1 +
			        k1p2*l1k2*l2k1*p1l2 - 2*k1p2*l1k2*m2*p1l2 + l1p2*l2k2*m2*p1l2 -
			        2*l1k2*l2p2*m2*p1l2 + l1k2*l2k1*m2*p1p2 - l2k1*l2k2*m2*p1p2 -
			        l2k2*m4*p1p2 + k1k2*(2*l1k1*(l2k1 - m2)*M2 - 2*l1l2*m2*M2 +
			           2*m4*M2 - l1p2*l2k1*p1k1 + l1p2*m2*p1k1 - k1p2*l2k1*p1l1 -
			           l2k1*l2p2*p1l1 + k1p2*m2*p1l1 + l2p2*m2*p1l1 -
			           l1p2*l2k1*p1l2 + l1p2*m2*p1l2 + l1l2*l2k1*p1p2 -
			           l2k1*m2*p1p2 - m4*p1p2) - l1l2*l2k1*p1k1*p2k2 +
			        l1l2*m2*p1k1*p2k2 + 2*l2k1*m2*p1k1*p2k2 +
			        2*eg1*l2k1*m2*p1k2*p2k2 + 2*eg1*m4*p1k2*p2k2 -
			        l2k1*m2*p1l1*p2k2 + l1l2*m2*p1l2*p2k2 + l2k1*m2*p1l2*p2k2 +
			        l1k1*(l2k1*(2*l2k2*M2 - k1p2*p1k2 - l2k2*p1p2 - p1k1*p2k2) +
			           m2*(-2*l2k2*M2 + k1p2*p1k2 + l2p2*p1k2 + p1k1*p2k2 +
			              p1l2*p2k2)) + l1p2*p1k2*pow(l2k1,2) - l1k2*p1p2*pow(l2k1,2) +
			        p1l1*p2k2*pow(l2k1,2)) -
			     4*f1*f2*M2*(l1p2*l2k2*l2p2*m2 - 3*l1k1*l2k1*l2k2*M2 +
			        2*l1k1*l2k2*m2*M2 + 2*l1l2*l2k2*m2*M2 -
			        2*l2k1*l2k2*m2*M2 - l1k2*m4*M2 - 3*l2k2*m4*M2 -
			        l1k2*l2k1*l2p2*p1k1 - l1p2*l2k2*m2*p1k1 + 2*l1k2*l2p2*m2*p1k1 +
			        l1p2*l2k1*m2*p1k2 - l1k1*l2p2*m2*p1k2 - l1l2*l2p2*m2*p1k2 -
			        l2k1*l2p2*m2*p1k2 - l1k1*l2k1*p1k1*p1k2 - l1l2*l2k1*p1k1*p1k2 +
			        l1k1*m2*p1k1*p1k2 + l1l2*m2*p1k1*p1k2 + 2*l2k1*m2*p1k1*p1k2 -
			        l2k2*l2p2*m2*p1l1 + l2k2*m2*p1k1*p1l1 - l2k1*m2*p1k2*p1l1 -
			        l1p2*l2k2*m2*p1l2 + 2*l1k2*l2p2*m2*p1l2 + l1k2*l2k1*p1k1*p1l2 -
			        2*l1k2*m2*p1k1*p1l2 + l1k1*m2*p1k2*p1l2 + l1l2*m2*p1k2*p1l2 +
			        l2k1*m2*p1k2*p1l2 + l2k2*m2*p1l1*p1l2 + 3*l1k1*l2k1*l2k2*p1p2 -
			        2*l1k1*l2k2*m2*p1p2 - 2*l1l2*l2k2*m2*p1p2 +
			        2*l2k1*l2k2*m2*p1p2 + l1k2*m4*p1p2 + 3*l2k2*m4*p1p2 +
			        k1k2*(l1l2*l2k1*M2 + 2*l1l2*m2*M2 - l2k1*m2*M2 -
			           3*m4*M2 + l2k1*l2p2*p1l1 - l2p2*m2*p1l1 - l2k1*p1k1*p1l1 +
			           m2*p1k1*p1l1 - l1p2*(l2k1 - m2)*(l2p2 - p1k1 - p1l2) -
			           l2k1*p1l1*p1l2 + m2*p1l1*p1l2 -
			           2*l1k1*(l2k1 - m2)*(M2 - p1p2) - l1l2*l2k1*p1p2 -
			           2*l1l2*m2*p1p2 + l2k1*m2*p1p2 + 3*m4*p1p2) -
			        l1p2*l2k1*m2*p2k2 + l1k1*l2p2*m2*p2k2 + l1l2*l2p2*m2*p2k2 +
			        l2k1*l2p2*m2*p2k2 + l1k1*l2k1*p1k1*p2k2 + l1l2*l2k1*p1k1*p2k2 -
			        l1k1*m2*p1k1*p2k2 - l1l2*m2*p1k1*p2k2 - 2*l2k1*m2*p1k1*p2k2 -
			        2*eg1*l2k1*m2*p1k2*p2k2 - 2*eg1*m4*p1k2*p2k2 +
			        l2k1*m2*p1l1*p2k2 - l1k1*m2*p1l2*p2k2 - l1l2*m2*p1l2*p2k2 -
			        l2k1*m2*p1l2*p2k2 + k1p2*
			         (l1p2*l2k2*m2 + l1k1*l2k1*p1k2 + l1l2*l2k1*p1k2 -
			           l1k1*m2*p1k2 - l1l2*m2*p1k2 - 2*l2k1*m2*p1k2 -
			           k1k2*(l2k1 - m2)*(l1p2 - p1l1) - l2k2*m2*p1l1 +
			           l1k2*(l2k1*(l2p2 - 2*p1k1 - p1l2) +
			              2*m2*(-l2p2 + p1k1 + p1l2)) - l1k1*l2k1*p2k2 -
			           l1l2*l2k1*p2k2 + l1k1*m2*p2k2 + l1l2*m2*p2k2 + 2*l2k1*m2*p2k2\
			) + l1k2*(l2k1 - m2)*pow(k1p2,2) - l1k2*M2*pow(l2k1,2) -
			        l1p2*p1k2*pow(l2k1,2) + p1k2*p1l1*pow(l2k1,2) +
			        l1k2*p1p2*pow(l2k1,2) + l1p2*p2k2*pow(l2k1,2) -
			        p1l1*p2k2*pow(l2k1,2) - l1k2*m2*pow(l2p2,2) +
			        l1k2*l2k1*pow(p1k1,2) - l1k2*m2*pow(p1k1,2) +
			        eg1*l2k1*m2*pow(p1k2,2) + eg1*m4*pow(p1k2,2) -
			        l1k2*m2*pow(p1l2,2) + eg1*l2k1*m2*pow(p2k2,2) +
			        eg1*m4*pow(p2k2,2)) +
			     pow(f2,2)*(-3*l1p2*l2k2*l2p2*m2*M2 + l1k2*l2k1*l2p2*M2*p1k1 +
			        l1p2*l2k2*m2*M2*p1k1 - 2*l1k2*l2p2*m2*M2*p1k1 -
			        l1p2*l2k1*m2*M2*p1k2 + l1k1*l2p2*m2*M2*p1k2 +
			        l1l2*l2p2*m2*M2*p1k2 + l2k1*l2p2*m2*M2*p1k2 +
			        3*l1k1*l2k1*M2*p1k1*p1k2 + 3*l1l2*l2k1*M2*p1k1*p1k2 -
			        3*l1k1*m2*M2*p1k1*p1k2 - 3*l1l2*m2*M2*p1k1*p1k2 -
			        6*l2k1*m2*M2*p1k1*p1k2 + l2k2*l2p2*m2*M2*p1l1 -
			        3*l2k2*m2*M2*p1k1*p1l1 + 3*l2k1*m2*M2*p1k2*p1l1 +
			        l1p2*l2k2*m2*M2*p1l2 - 2*l1k2*l2p2*m2*M2*p1l2 -
			        3*l1k2*l2k1*M2*p1k1*p1l2 + 6*l1k2*m2*M2*p1k1*p1l2 -
			        3*l1k1*m2*M2*p1k2*p1l2 - 3*l1l2*m2*M2*p1k2*p1l2 -
			        3*l2k1*m2*M2*p1k2*p1l2 - 3*l2k2*m2*M2*p1l1*p1l2 +
			        l1p2*l2k2*l2p2*m2*p1p2 - 6*l1k1*l2k1*l2k2*M2*p1p2 +
			        4*l1k1*l2k2*m2*M2*p1p2 + 4*l1l2*l2k2*m2*M2*p1p2 -
			        4*l2k1*l2k2*m2*M2*p1p2 - 2*l1k2*m4*M2*p1p2 -
			        6*l2k2*m4*M2*p1p2 + l1k2*l2k1*l2p2*p1k1*p1p2 +
			        l1p2*l2k2*m2*p1k1*p1p2 - 2*l1k2*l2p2*m2*p1k1*p1p2 -
			        l1p2*l2k1*m2*p1k2*p1p2 + l1k1*l2p2*m2*p1k2*p1p2 +
			        l1l2*l2p2*m2*p1k2*p1p2 + l2k1*l2p2*m2*p1k2*p1p2 -
			        l1k1*l2k1*p1k1*p1k2*p1p2 - l1l2*l2k1*p1k1*p1k2*p1p2 +
			        l1k1*m2*p1k1*p1k2*p1p2 + l1l2*m2*p1k1*p1k2*p1p2 +
			        2*l2k1*m2*p1k1*p1k2*p1p2 + l2k2*l2p2*m2*p1l1*p1p2 +
			        l2k2*m2*p1k1*p1l1*p1p2 - l2k1*m2*p1k2*p1l1*p1p2 +
			        l1p2*l2k2*m2*p1l2*p1p2 - 2*l1k2*l2p2*m2*p1l2*p1p2 +
			        l1k2*l2k1*p1k1*p1l2*p1p2 - 2*l1k2*m2*p1k1*p1l2*p1p2 +
			        l1k1*m2*p1k2*p1l2*p1p2 + l1l2*m2*p1k2*p1l2*p1p2 +
			        l2k1*m2*p1k2*p1l2*p1p2 + l2k2*m2*p1l1*p1l2*p1p2 +
			        3*l1p2*l2k1*m2*M2*p2k2 - 3*l1k1*l2p2*m2*M2*p2k2 -
			        3*l1l2*l2p2*m2*M2*p2k2 - 3*l2k1*l2p2*m2*M2*p2k2 -
			        l1k1*l2k1*M2*p1k1*p2k2 - l1l2*l2k1*M2*p1k1*p2k2 +
			        l1k1*m2*M2*p1k1*p2k2 + l1l2*m2*M2*p1k1*p2k2 +
			        2*l2k1*m2*M2*p1k1*p2k2 + 2*eg1*l2k1*m2*M2*p1k2*p2k2 +
			        2*eg1*m4*M2*p1k2*p2k2 - l2k1*m2*M2*p1l1*p2k2 +
			        l1k1*m2*M2*p1l2*p2k2 + l1l2*m2*M2*p1l2*p2k2 +
			        l2k1*m2*M2*p1l2*p2k2 - l1p2*l2k1*m2*p1p2*p2k2 +
			        l1k1*l2p2*m2*p1p2*p2k2 + l1l2*l2p2*m2*p1p2*p2k2 +
			        l2k1*l2p2*m2*p1p2*p2k2 - l1k1*l2k1*p1k1*p1p2*p2k2 -
			        l1l2*l2k1*p1k1*p1p2*p2k2 + l1k1*m2*p1k1*p1p2*p2k2 +
			        l1l2*m2*p1k1*p1p2*p2k2 + 2*l2k1*m2*p1k1*p1p2*p2k2 +
			        2*eg1*l2k1*m2*p1k2*p1p2*p2k2 + 2*eg1*m4*p1k2*p1p2*p2k2 -
			        l2k1*m2*p1l1*p1p2*p2k2 + l1k1*m2*p1l2*p1p2*p2k2 +
			        l1l2*m2*p1l2*p1p2*p2k2 + l2k1*m2*p1l2*p1p2*p2k2 +
			        k1p2*(-3*l1p2*l2k2*m2*M2 - l1k1*l2k1*M2*p1k2 -
			           l1l2*l2k1*M2*p1k2 + l1k1*m2*M2*p1k2 + l1l2*m2*M2*p1k2 +
			           2*l2k1*m2*M2*p1k2 + l2k2*m2*M2*p1l1 + l1p2*l2k2*m2*p1p2 -
			           l1k1*l2k1*p1k2*p1p2 - l1l2*l2k1*p1k2*p1p2 + l1k1*m2*p1k2*p1p2 +
			           l1l2*m2*p1k2*p1p2 + 2*l2k1*m2*p1k2*p1p2 + l2k2*m2*p1l1*p1p2 +
			           k1k2*(l2k1 - m2)*(l1p2*(3*M2 - p1p2) - p1l1*(M2 + p1p2)) +
			           l1k2*(2*m2*(l2p2*(3*M2 - p1p2) -
			                 (p1k1 + p1l2)*(M2 + p1p2)) +
			              l2k1*(l2p2*(-3*M2 + p1p2) + (2*p1k1 + p1l2)*(M2 + p1p2))) \
			+ 3*l1k1*l2k1*M2*p2k2 + 3*l1l2*l2k1*M2*p2k2 - 3*l1k1*m2*M2*p2k2 -
			           3*l1l2*m2*M2*p2k2 - 6*l2k1*m2*M2*p2k2 -
			           l1k1*l2k1*p1p2*p2k2 - l1l2*l2k1*p1p2*p2k2 + l1k1*m2*p1p2*p2k2 +
			           l1l2*m2*p1p2*p2k2 + 2*l2k1*m2*p1p2*p2k2) -
			        l1k2*(l2k1 - m2)*(3*M2 - p1p2)*pow(k1p2,2) +
			        l1p2*M2*p1k2*pow(l2k1,2) - 3*M2*p1k2*p1l1*pow(l2k1,2) -
			        2*l1k2*M2*p1p2*pow(l2k1,2) + l1p2*p1k2*p1p2*pow(l2k1,2) +
			        p1k2*p1l1*p1p2*pow(l2k1,2) - 3*l1p2*M2*p2k2*pow(l2k1,2) +
			        M2*p1l1*p2k2*pow(l2k1,2) + l1p2*p1p2*p2k2*pow(l2k1,2) +
			        p1l1*p1p2*p2k2*pow(l2k1,2) + 3*l1k2*m2*M2*pow(l2p2,2) -
			        l1k2*m2*p1p2*pow(l2p2,2) + 5*l1k1*l2k1*l2k2*M4 -
			        2*l1k2*l2k1*m2*M4 - 2*l1k1*l2k2*m2*M4 -
			        2*l1l2*l2k2*m2*M4 + 4*l2k1*l2k2*m2*M4 +
			        l1k2*m4*M4 + 5*l2k2*m4*M4 +
			        3*l1k2*pow(l2k1,2)*M4 - 3*l1k2*l2k1*M2*pow(p1k1,2) +
			        3*l1k2*m2*M2*pow(p1k1,2) + l1k2*l2k1*p1p2*pow(p1k1,2) -
			        l1k2*m2*p1p2*pow(p1k1,2) - 3*eg1*l2k1*m2*M2*pow(p1k2,2) -
			        3*eg1*m4*M2*pow(p1k2,2) + eg1*l2k1*m2*p1p2*pow(p1k2,2) +
			        eg1*m4*p1p2*pow(p1k2,2) + 3*l1k2*m2*M2*pow(p1l2,2) -
			        l1k2*m2*p1p2*pow(p1l2,2) + l1k1*l2k1*l2k2*pow(p1p2,2) +
			        2*l1k2*l2k1*m2*pow(p1p2,2) - 2*l1k1*l2k2*m2*pow(p1p2,2) -
			        2*l1l2*l2k2*m2*pow(p1p2,2) + l1k2*m4*pow(p1p2,2) +
			        l2k2*m4*pow(p1p2,2) - l1k2*pow(l2k1,2)*pow(p1p2,2) +
			        k1k2*(-(l2k1*l2p2*M2*p1l1) + l2p2*m2*M2*p1l1 +
			           3*l2k1*M2*p1k1*p1l1 - 3*m2*M2*p1k1*p1l1 +
			           3*l2k1*M2*p1l1*p1l2 - 3*m2*M2*p1l1*p1l2 +
			           2*l1l2*l2k1*M2*p1p2 + 4*l1l2*m2*M2*p1p2 -
			           2*l2k1*m2*M2*p1p2 - 6*m4*M2*p1p2 - l2k1*l2p2*p1l1*p1p2 +
			           l2p2*m2*p1l1*p1p2 - l2k1*p1k1*p1l1*p1p2 + m2*p1k1*p1l1*p1p2 -
			           l2k1*p1l1*p1l2*p1p2 + m2*p1l1*p1l2*p1p2 +
			           l1p2*(l2k1 - m2)*(l2p2*(3*M2 - p1p2) -
			              (p1k1 + p1l2)*(M2 + p1p2)) - 3*l1l2*l2k1*M4 -
			           2*l1l2*m2*M4 + 3*l2k1*m2*M4 +
			           5*m4*M4 + 2*l1k1*(l2k1 - m2)*pow(M2 - p1p2,2) +
			           l1l2*l2k1*pow(p1p2,2) - 2*l1l2*m2*pow(p1p2,2) -
			           l2k1*m2*pow(p1p2,2) + m4*pow(p1p2,2)) -
			        3*eg1*l2k1*m2*M2*pow(p2k2,2) - 3*eg1*m4*M2*pow(p2k2,2) +
			        eg1*l2k1*m2*p1p2*pow(p2k2,2) + eg1*m4*p1p2*pow(p2k2,2))) -
			  4*eg1*pow(l1k1,-1)*pow(l1k1 + eg1*(-k1k2 + l1k2),-1)*pow(l2k1,-1)*
			   pow(l2k1 + eg1*(k1k2 + l2k2),-1)/M2*
			   (-4*M2*pow(f1,2)*(m2*pow(eg1,2)*
			         (k1p2*(l1k2 - l2k2)*p1k2 +
			           (l1k2*p1k1 - l2k2*p1k1 - 2*(l1k1 + 2*l1l2 - l2k1)*p1k2)*p2k2 -
			           k1k2*(2*k1p2*p1k2 + l1p2*p1k2 - l2p2*p1k2 + 2*p1k1*p2k2 +
			              p1l1*p2k2 - p1l2*p2k2) + 2*M2*pow(k1k2,2)) +
			        eg1*(-4*l1k1*l1k2*l1l2*M2 + 4*l1k2*l1l2*l2k1*M2 +
			           4*l1k1*l1l2*l2k2*M2 - 4*l1l2*l2k1*l2k2*M2 +
			           2*l1k1*l1k2*m2*M2 + 2*l1k2*l1l2*m2*M2 +
			           2*l1k2*l2k1*m2*M2 + 2*l1k1*l2k2*m2*M2 -
			           2*l1l2*l2k2*m2*M2 + 2*l2k1*l2k2*m2*M2 -
			           2*k1p2*l1k2*m2*p1k1 + 2*l1k2*l1p2*m2*p1k1 +
			           2*k1p2*l2k2*m2*p1k1 + 2*l1p2*l2k2*m2*p1k1 +
			           2*l1k2*l2p2*m2*p1k1 + 2*l2k2*l2p2*m2*p1k1 +
			           2*l1k1*l1l2*l1p2*p1k2 - 2*l1l2*l1p2*l2k1*p1k2 -
			           2*l1k1*l1l2*l2p2*p1k2 + 2*l1l2*l2k1*l2p2*p1k2 +
			           k1p2*l1k1*m2*p1k2 - 2*k1p2*l1l2*m2*p1k2 +
			           2*l1k1*l1p2*m2*p1k2 - k1p2*l2k1*m2*p1k2 -
			           2*l1p2*l2k1*m2*p1k2 - 2*l1k1*l2p2*m2*p1k2 +
			           2*l2k1*l2p2*m2*p1k2 - 2*k1p2*m4*p1k2 -
			           4*l1k1*l1p2*l2k2*p1l1 - 4*l1l2*l1p2*l2k2*p1l1 +
			           4*l1p2*l2k1*l2k2*p1l1 + 2*l1k1*l1k2*l2p2*p1l1 +
			           2*l1k2*l1l2*l2p2*p1l1 - 2*l1k2*l2k1*l2p2*p1l1 -
			           2*l1k1*l2k2*l2p2*p1l1 - 2*l1l2*l2k2*l2p2*p1l1 +
			           2*l2k1*l2k2*l2p2*p1l1 + 2*k1p2*l1k2*m2*p1l1 +
			           2*k1p2*l2k2*m2*p1l1 + 2*l1k1*l1k2*l1p2*p1l2 +
			           2*l1k2*l1l2*l1p2*p1l2 - 2*l1k2*l1p2*l2k1*p1l2 -
			           2*l1k1*l1p2*l2k2*p1l2 - 2*l1l2*l1p2*l2k2*p1l2 +
			           2*l1p2*l2k1*l2k2*p1l2 + 4*l1k1*l1k2*l2p2*p1l2 +
			           4*l1k2*l1l2*l2p2*p1l2 - 4*l1k2*l2k1*l2p2*p1l2 +
			           2*k1p2*l1k2*m2*p1l2 + 2*k1p2*l2k2*m2*p1l2 -
			           4*l1k1*l1k2*m2*p1p2 - 2*l1k2*l1l2*m2*p1p2 +
			           2*l1l2*l2k2*m2*p1p2 - 4*l2k1*l2k2*m2*p1p2 +
			           l1k1*m2*p1k1*p2k2 - 2*l1l2*m2*p1k1*p2k2 -
			           l2k1*m2*p1k1*p2k2 - 2*m4*p1k1*p2k2 + 2*l1k1*l1l2*p1l1*p2k2 -
			           2*l1l2*l2k1*p1l1*p2k2 + 2*l1k1*m2*p1l1*p2k2 -
			           2*l2k1*m2*p1l1*p2k2 - 2*l1k1*l1l2*p1l2*p2k2 +
			           2*l1l2*l2k1*p1l2*p2k2 - 2*l1k1*m2*p1l2*p2k2 +
			           2*l2k1*m2*p1l2*p2k2 - 4*l1k2*M2*pow(l1l2,2) +
			           4*l2k2*M2*pow(l1l2,2) + 2*l1p2*p1k2*pow(l1l2,2) -
			           2*l2p2*p1k2*pow(l1l2,2) + 2*p1l1*p2k2*pow(l1l2,2) -
			           2*p1l2*p2k2*pow(l1l2,2) +
			           k1k2*(-(m2*(-2*l2k1*M2 - l1p2*p1k1 + l2p2*p1k1 -
			                   k1p2*p1l1 + 4*l1p2*p1l1 + 4*l2p2*p1l1 + k1p2*p1l2 +
			                   4*l1p2*p1l2 + 4*l2p2*p1l2 + 6*m2*(M2 - p1p2))) +
			              l1k1*(8*l1l2*M2 - 2*m2*M2 - 4*l1l2*p1p2) +
			              l1l2*(6*m2*(-M2 + p1p2) + l2k1*(-8*M2 + 4*p1p2)) +
			              (8*M2 - 4*p1p2)*pow(l1l2,2))) +
			        2*(l2k1*(l2k1*(2*l1p2*p1l1 + l2p2*p1l1 + l1p2*p1l2 +
			                 2*m2*(M2 - p1p2)) +
			              m2*(l1p2*p1k1 + l2p2*p1k1 + k1p2*p1l1 - l2p2*p1l1 +
			                 k1p2*p1l2 - l1p2*p1l2 + m2*(-2*M2 + p1p2))) +
			           (-2*l1l2*M2 + l2p2*p1l1 + l1p2*p1l2 + 2*l2p2*p1l2 +
			              2*m2*(M2 - p1p2))*pow(l1k1,2) +
			           (4*l2k1*M2 + l1p2*p1k1 - l2p2*p1k1 + k1p2*p1l1 +
			              2*l2p2*p1l1 - k1p2*p1l2 + 2*l1p2*p1l2 +
			              m2*(4*M2 - 2*p1p2))*pow(l1l2,2) +
			           l1k1*(-2*l2k1*(l1p2 + l2p2)*(p1l1 + p1l2) +
			              m2*(l2p2*(p1k1 + p1l1) + l1p2*(p1k1 + p1l2) +
			                 k1p2*(p1l1 + p1l2)) +
			              l1l2*(4*l2k1*M2 + l1p2*p1k1 - l2p2*p1k1 + k1p2*p1l1 +
			                 3*l2p2*p1l1 - k1p2*p1l2 + 3*l1p2*p1l2 + 2*l2p2*p1l2 +
			                 m2*(4*M2 - 3*p1p2)) + m4*(2*M2 - p1p2) -
			              4*M2*pow(l1l2,2)) - 2*M2*pow(l1l2,3) -
			           l1l2*(2*k1p2*m2*p1k1 +
			              l2k1*(-(l2p2*p1k1) + k1p2*p1l1 + 3*l2p2*p1l1 - k1p2*p1l2 +
			                 l1p2*(p1k1 + 2*p1l1 + 3*p1l2) + m2*(4*M2 - 3*p1p2)) +
			              2*M2*pow(l2k1,2)))) -
			     4*f1*f2*M2*(-2*eg1*l1k1*l1k2*l1p2*l2p2 - 6*l1k1*l1l2*l1p2*l2p2 -
			        2*eg1*l1k2*l1l2*l1p2*l2p2 + 4*l1k1*l1p2*l2k1*l2p2 +
			        2*eg1*l1k2*l1p2*l2k1*l2p2 + 6*l1l2*l1p2*l2k1*l2p2 +
			        2*eg1*l1k1*l1p2*l2k2*l2p2 + 2*eg1*l1l2*l1p2*l2k2*l2p2 -
			        2*eg1*l1p2*l2k1*l2k2*l2p2 + 4*eg1*k1k2*l1p2*l2p2*m2 -
			        2*l1k1*l1p2*l2p2*m2 + 2*l1p2*l2k1*l2p2*m2 +
			        12*eg1*k1k2*l1k1*l1l2*M2 - 4*eg1*l1k1*l1k2*l1l2*M2 -
			        12*eg1*k1k2*l1l2*l2k1*M2 + 8*l1k1*l1l2*l2k1*M2 +
			        4*eg1*l1k2*l1l2*l2k1*M2 + 4*eg1*l1k1*l1l2*l2k2*M2 -
			        4*eg1*l1l2*l2k1*l2k2*M2 - 2*eg1*k1k2*l1k1*m2*M2 +
			        6*eg1*l1k1*l1k2*m2*M2 - 12*eg1*k1k2*l1l2*m2*M2 +
			        14*l1k1*l1l2*m2*M2 + 4*eg1*l1k2*l1l2*m2*M2 +
			        2*eg1*k1k2*l2k1*m2*M2 + 2*eg1*l1k2*l2k1*m2*M2 -
			        14*l1l2*l2k1*m2*M2 + 2*eg1*l1k1*l2k2*m2*M2 -
			        4*eg1*l1l2*l2k2*m2*M2 + 6*eg1*l2k1*l2k2*m2*M2 -
			        12*eg1*k1k2*m4*M2 + 6*l1k1*m4*M2 - 6*l2k1*m4*M2 +
			        2*l1k1*l1l2*l1p2*p1k1 - 2*l1l2*l1p2*l2k1*p1k1 -
			        2*l1k1*l1l2*l2p2*p1k1 + 2*l1l2*l2k1*l2p2*p1k1 +
			        eg1*k1k2*l1p2*m2*p1k1 + 2*l1k1*l1p2*m2*p1k1 +
			        2*eg1*l1k2*l1p2*m2*p1k1 + 2*l1p2*l2k1*m2*p1k1 +
			        2*eg1*l1p2*l2k2*m2*p1k1 - eg1*k1k2*l2p2*m2*p1k1 +
			        2*l1k1*l2p2*m2*p1k1 + 2*eg1*l1k2*l2p2*m2*p1k1 +
			        2*l2k1*l2p2*m2*p1k1 + 2*eg1*l2k2*l2p2*m2*p1k1 +
			        2*eg1*l1k1*l1l2*l1p2*p1k2 - 2*eg1*l1l2*l1p2*l2k1*p1k2 -
			        2*eg1*l1k1*l1l2*l2p2*p1k2 + 2*eg1*l1l2*l2k1*l2p2*p1k2 +
			        2*eg1*l1k1*l1p2*m2*p1k2 - 2*eg1*l1p2*l2k1*m2*p1k2 -
			        2*eg1*l1k1*l2p2*m2*p1k2 + 2*eg1*l2k1*l2p2*m2*p1k2 -
			        eg1*l1k1*m2*p1k1*p1k2 + 2*eg1*l1l2*m2*p1k1*p1k2 +
			        eg1*l2k1*m2*p1k1*p1k2 + 2*eg1*m4*p1k1*p1k2 -
			        4*l1k1*l1p2*l2k1*p1l1 - 4*l1l2*l1p2*l2k1*p1l1 -
			        4*eg1*l1k1*l1p2*l2k2*p1l1 - 4*eg1*l1l2*l1p2*l2k2*p1l1 +
			        4*eg1*l1p2*l2k1*l2k2*p1l1 + 2*eg1*l1k1*l1k2*l2p2*p1l1 +
			        6*l1k1*l1l2*l2p2*p1l1 + 2*eg1*l1k2*l1l2*l2p2*p1l1 -
			        4*l1k1*l2k1*l2p2*p1l1 - 2*eg1*l1k2*l2k1*l2p2*p1l1 -
			        6*l1l2*l2k1*l2p2*p1l1 - 2*eg1*l1k1*l2k2*l2p2*p1l1 -
			        2*eg1*l1l2*l2k2*l2p2*p1l1 + 2*eg1*l2k1*l2k2*l2p2*p1l1 -
			        4*eg1*k1k2*l1p2*m2*p1l1 - 4*eg1*k1k2*l2p2*m2*p1l1 +
			        2*l1k1*l2p2*m2*p1l1 - 2*l2k1*l2p2*m2*p1l1 -
			        2*l1k1*l1l2*p1k1*p1l1 + 2*l1l2*l2k1*p1k1*p1l1 -
			        eg1*k1k2*m2*p1k1*p1l1 - 2*l1k1*m2*p1k1*p1l1 -
			        2*eg1*l1k2*m2*p1k1*p1l1 - 2*l2k1*m2*p1k1*p1l1 -
			        2*eg1*l2k2*m2*p1k1*p1l1 - 2*eg1*l1k1*l1l2*p1k2*p1l1 +
			        2*eg1*l1l2*l2k1*p1k2*p1l1 - 2*eg1*l1k1*m2*p1k2*p1l1 +
			        2*eg1*l2k1*m2*p1k2*p1l1 + 2*eg1*l1k1*l1k2*l1p2*p1l2 +
			        6*l1k1*l1l2*l1p2*p1l2 + 2*eg1*l1k2*l1l2*l1p2*p1l2 -
			        4*l1k1*l1p2*l2k1*p1l2 - 2*eg1*l1k2*l1p2*l2k1*p1l2 -
			        6*l1l2*l1p2*l2k1*p1l2 - 2*eg1*l1k1*l1p2*l2k2*p1l2 -
			        2*eg1*l1l2*l1p2*l2k2*p1l2 + 2*eg1*l1p2*l2k1*l2k2*p1l2 +
			        4*eg1*l1k1*l1k2*l2p2*p1l2 + 4*l1k1*l1l2*l2p2*p1l2 +
			        4*eg1*l1k2*l1l2*l2p2*p1l2 - 4*l1k1*l2k1*l2p2*p1l2 -
			        4*eg1*l1k2*l2k1*l2p2*p1l2 - 4*eg1*k1k2*l1p2*m2*p1l2 +
			        2*l1k1*l1p2*m2*p1l2 - 2*l1p2*l2k1*m2*p1l2 -
			        4*eg1*k1k2*l2p2*m2*p1l2 + 2*l1k1*l1l2*p1k1*p1l2 -
			        2*l1l2*l2k1*p1k1*p1l2 + eg1*k1k2*m2*p1k1*p1l2 -
			        2*l1k1*m2*p1k1*p1l2 - 2*eg1*l1k2*m2*p1k1*p1l2 -
			        2*l2k1*m2*p1k1*p1l2 - 2*eg1*l2k2*m2*p1k1*p1l2 +
			        2*eg1*l1k1*l1l2*p1k2*p1l2 - 2*eg1*l1l2*l2k1*p1k2*p1l2 +
			        2*eg1*l1k1*m2*p1k2*p1l2 - 2*eg1*l2k1*m2*p1k2*p1l2 -
			        2*eg1*l1k1*l1k2*p1l1*p1l2 - 6*l1k1*l1l2*p1l1*p1l2 -
			        2*eg1*l1k2*l1l2*p1l1*p1l2 + 4*l1k1*l2k1*p1l1*p1l2 +
			        2*eg1*l1k2*l2k1*p1l1*p1l2 + 6*l1l2*l2k1*p1l1*p1l2 +
			        2*eg1*l1k1*l2k2*p1l1*p1l2 + 2*eg1*l1l2*l2k2*p1l1*p1l2 -
			        2*eg1*l2k1*l2k2*p1l1*p1l2 + 4*eg1*k1k2*m2*p1l1*p1l2 -
			        2*l1k1*m2*p1l1*p1l2 + 2*l2k1*m2*p1l1*p1l2 -
			        12*eg1*k1k2*l1k1*l1l2*p1p2 + 4*eg1*l1k1*l1k2*l1l2*p1p2 +
			        12*eg1*k1k2*l1l2*l2k1*p1p2 - 8*l1k1*l1l2*l2k1*p1p2 -
			        4*eg1*l1k2*l1l2*l2k1*p1p2 - 4*eg1*l1k1*l1l2*l2k2*p1p2 +
			        4*eg1*l1l2*l2k1*l2k2*p1p2 + 2*eg1*k1k2*l1k1*m2*p1p2 -
			        6*eg1*l1k1*l1k2*m2*p1p2 + 12*eg1*k1k2*l1l2*m2*p1p2 -
			        14*l1k1*l1l2*m2*p1p2 - 4*eg1*l1k2*l1l2*m2*p1p2 -
			        2*eg1*k1k2*l2k1*m2*p1p2 - 2*eg1*l1k2*l2k1*m2*p1p2 +
			        14*l1l2*l2k1*m2*p1p2 - 2*eg1*l1k1*l2k2*m2*p1p2 +
			        4*eg1*l1l2*l2k2*m2*p1p2 - 6*eg1*l2k1*l2k2*m2*p1p2 +
			        12*eg1*k1k2*m4*p1p2 - 6*l1k1*m4*p1p2 + 6*l2k1*m4*p1p2 -
			        2*eg1*l1k1*l1l2*l1p2*p2k2 + 2*eg1*l1l2*l1p2*l2k1*p2k2 +
			        2*eg1*l1k1*l1l2*l2p2*p2k2 - 2*eg1*l1l2*l2k1*l2p2*p2k2 -
			        2*eg1*l1k1*l1p2*m2*p2k2 + 2*eg1*l1p2*l2k1*m2*p2k2 +
			        2*eg1*l1k1*l2p2*m2*p2k2 - 2*eg1*l2k1*l2p2*m2*p2k2 +
			        eg1*l1k1*m2*p1k1*p2k2 - 2*eg1*l1l2*m2*p1k1*p2k2 -
			        eg1*l2k1*m2*p1k1*p2k2 - 2*eg1*m4*p1k1*p2k2 +
			        2*eg1*l1k1*l1l2*p1l1*p2k2 - 2*eg1*l1l2*l2k1*p1l1*p2k2 +
			        2*eg1*l1k1*m2*p1l1*p2k2 - 2*eg1*l2k1*m2*p1l1*p2k2 -
			        2*eg1*l1k1*l1l2*p1l2*p2k2 + 2*eg1*l1l2*l2k1*p1l2*p2k2 -
			        2*eg1*l1k1*m2*p1l2*p2k2 + 2*eg1*l2k1*m2*p1l2*p2k2 -
			        k1k2*l1p2*m2*p1k2*pow(eg1,2) + k1k2*l2p2*m2*p1k2*pow(eg1,2) +
			        2*k1k2*m2*p1k1*p1k2*pow(eg1,2) - l1k2*m2*p1k1*p1k2*pow(eg1,2) +
			        l2k2*m2*p1k1*p1k2*pow(eg1,2) + k1k2*m2*p1k2*p1l1*pow(eg1,2) -
			        k1k2*m2*p1k2*p1l2*pow(eg1,2) + k1k2*l1p2*m2*p2k2*pow(eg1,2) -
			        k1k2*l2p2*m2*p2k2*pow(eg1,2) - 2*k1k2*m2*p1k1*p2k2*pow(eg1,2) +
			        l1k2*m2*p1k1*p2k2*pow(eg1,2) - l2k2*m2*p1k1*p2k2*pow(eg1,2) -
			        2*l1k1*m2*p1k2*p2k2*pow(eg1,2) - 4*l1l2*m2*p1k2*p2k2*pow(eg1,2) +
			        2*l2k1*m2*p1k2*p2k2*pow(eg1,2) - k1k2*m2*p1l1*p2k2*pow(eg1,2) +
			        k1k2*m2*p1l2*p2k2*pow(eg1,2) + 2*m2*M2*pow(eg1,2)*pow(k1k2,2) -
			        2*m2*p1p2*pow(eg1,2)*pow(k1k2,2) +
			        (2*l1l2 + eg1*(l1k2 - l2k2))*m2*pow(k1p2,2) -
			        2*l1p2*l2p2*pow(l1k1,2) - 4*l1l2*M2*pow(l1k1,2) +
			        8*m2*M2*pow(l1k1,2) + 2*l2p2*p1l1*pow(l1k1,2) +
			        2*l1p2*p1l2*pow(l1k1,2) + 4*l2p2*p1l2*pow(l1k1,2) -
			        2*p1l1*p1l2*pow(l1k1,2) + 4*l1l2*p1p2*pow(l1k1,2) -
			        8*m2*p1p2*pow(l1k1,2) - 4*l1p2*l2p2*pow(l1l2,2) +
			        12*eg1*k1k2*M2*pow(l1l2,2) - 8*l1k1*M2*pow(l1l2,2) -
			        4*eg1*l1k2*M2*pow(l1l2,2) + 8*l2k1*M2*pow(l1l2,2) +
			        4*eg1*l2k2*M2*pow(l1l2,2) + 12*m2*M2*pow(l1l2,2) +
			        2*l1p2*p1k1*pow(l1l2,2) - 2*l2p2*p1k1*pow(l1l2,2) +
			        2*eg1*l1p2*p1k2*pow(l1l2,2) - 2*eg1*l2p2*p1k2*pow(l1l2,2) +
			        4*l2p2*p1l1*pow(l1l2,2) - 2*p1k1*p1l1*pow(l1l2,2) -
			        2*eg1*p1k2*p1l1*pow(l1l2,2) + 4*l1p2*p1l2*pow(l1l2,2) +
			        2*p1k1*p1l2*pow(l1l2,2) + 2*eg1*p1k2*p1l2*pow(l1l2,2) -
			        4*p1l1*p1l2*pow(l1l2,2) - 12*eg1*k1k2*p1p2*pow(l1l2,2) +
			        8*l1k1*p1p2*pow(l1l2,2) + 4*eg1*l1k2*p1p2*pow(l1l2,2) -
			        8*l2k1*p1p2*pow(l1l2,2) - 4*eg1*l2k2*p1p2*pow(l1l2,2) -
			        12*m2*p1p2*pow(l1l2,2) - 2*eg1*l1p2*p2k2*pow(l1l2,2) +
			        2*eg1*l2p2*p2k2*pow(l1l2,2) + 2*eg1*p1l1*p2k2*pow(l1l2,2) -
			        2*eg1*p1l2*p2k2*pow(l1l2,2) +
			        k1p2*(l1k1*(-2*l1l2*(l1p2 - l2p2 - p1l1 + p1l2) +
			              m2*(-2*l1p2 - 2*l2p2 + eg1*p1k2 + 2*p1l1 + 2*p1l2 -
			                 eg1*p2k2)) + 2*l1l2*
			            (l1p2*l2k1 - l2k1*(l2p2 + p1l1 - p1l2) +
			              m2*(-2*p1k1 - eg1*p1k2 + eg1*p2k2)) +
			           m2*(2*l2k1*(-l1p2 - l2p2 + p1l1 + p1l2) +
			              eg1*(-2*l1p2*l2k2 - 2*l2k2*l2p2 + 2*l2k2*p1k1 -
			                 l2k1*p1k2 - 2*m2*p1k2 + 2*l2k2*p1l1 -
			                 2*l1k2*(l1p2 + l2p2 + p1k1 - p1l1 - p1l2) +
			                 k1k2*(-l1p2 + l2p2 + p1l1 - p1l2) + 2*l2k2*p1l2 +
			                 l2k1*p2k2 + 2*m2*p2k2) -
			              (2*k1k2 - l1k2 + l2k2)*(p1k2 - p2k2)*pow(eg1,2)) -
			           2*(l1p2 - l2p2 - p1l1 + p1l2)*pow(l1l2,2)) -
			        4*M2*pow(l1l2,3) + 4*p1p2*pow(l1l2,3) + 2*l1k1*l2k1*pow(l1p2,2) +
			        2*l1l2*l2k1*pow(l1p2,2) + 2*eg1*l1k1*l2k2*pow(l1p2,2) +
			        2*eg1*l1l2*l2k2*pow(l1p2,2) - 2*eg1*l2k1*l2k2*pow(l1p2,2) +
			        2*eg1*k1k2*m2*pow(l1p2,2) - 2*l1p2*l2p2*pow(l2k1,2) -
			        4*l1l2*M2*pow(l2k1,2) + 8*m2*M2*pow(l2k1,2) +
			        4*l1p2*p1l1*pow(l2k1,2) + 2*l2p2*p1l1*pow(l2k1,2) +
			        2*l1p2*p1l2*pow(l2k1,2) - 2*p1l1*p1l2*pow(l2k1,2) +
			        4*l1l2*p1p2*pow(l2k1,2) - 8*m2*p1p2*pow(l2k1,2) -
			        2*pow(l1p2,2)*pow(l2k1,2) - 2*eg1*l1k1*l1k2*pow(l2p2,2) -
			        2*l1k1*l1l2*pow(l2p2,2) - 2*eg1*l1k2*l1l2*pow(l2p2,2) +
			        2*l1k1*l2k1*pow(l2p2,2) + 2*eg1*l1k2*l2k1*pow(l2p2,2) +
			        2*eg1*k1k2*m2*pow(l2p2,2) - 2*pow(l1k1,2)*pow(l2p2,2) +
			        eg1*l1k2*m2*pow(p1k1,2) + 2*l1l2*m2*pow(p1k1,2) -
			        eg1*l2k2*m2*pow(p1k1,2) + l1k1*m2*pow(eg1,2)*pow(p1k2,2) +
			        2*l1l2*m2*pow(eg1,2)*pow(p1k2,2) -
			        l2k1*m2*pow(eg1,2)*pow(p1k2,2) + 2*l1k1*l2k1*pow(p1l1,2) +
			        2*l1l2*l2k1*pow(p1l1,2) + 2*eg1*l1k1*l2k2*pow(p1l1,2) +
			        2*eg1*l1l2*l2k2*pow(p1l1,2) - 2*eg1*l2k1*l2k2*pow(p1l1,2) +
			        2*eg1*k1k2*m2*pow(p1l1,2) - 2*pow(l2k1,2)*pow(p1l1,2) -
			        2*eg1*l1k1*l1k2*pow(p1l2,2) - 2*l1k1*l1l2*pow(p1l2,2) -
			        2*eg1*l1k2*l1l2*pow(p1l2,2) + 2*l1k1*l2k1*pow(p1l2,2) +
			        2*eg1*l1k2*l2k1*pow(p1l2,2) + 2*eg1*k1k2*m2*pow(p1l2,2) -
			        2*pow(l1k1,2)*pow(p1l2,2) + l1k1*m2*pow(eg1,2)*pow(p2k2,2) +
			        2*l1l2*m2*pow(eg1,2)*pow(p2k2,2) - l2k1*m2*pow(eg1,2)*pow(p2k2,2)) \
			- pow(f2,2)*(-6*eg1*l1k1*l1k2*l1p2*l2p2*M2 - 18*l1k1*l1l2*l1p2*l2p2*M2 -
			        6*eg1*l1k2*l1l2*l1p2*l2p2*M2 + 12*l1k1*l1p2*l2k1*l2p2*M2 +
			        6*eg1*l1k2*l1p2*l2k1*l2p2*M2 + 18*l1l2*l1p2*l2k1*l2p2*M2 +
			        6*eg1*l1k1*l1p2*l2k2*l2p2*M2 + 6*eg1*l1l2*l1p2*l2k2*l2p2*M2 -
			        6*eg1*l1p2*l2k1*l2k2*l2p2*M2 + 12*eg1*k1k2*l1p2*l2p2*m2*M2 -
			        6*l1k1*l1p2*l2p2*m2*M2 + 6*l1p2*l2k1*l2p2*m2*M2 +
			        2*l1k1*l1l2*l1p2*M2*p1k1 - 2*l1l2*l1p2*l2k1*M2*p1k1 -
			        2*l1k1*l1l2*l2p2*M2*p1k1 + 2*l1l2*l2k1*l2p2*M2*p1k1 +
			        eg1*k1k2*l1p2*m2*M2*p1k1 + 2*l1k1*l1p2*m2*M2*p1k1 +
			        2*eg1*l1k2*l1p2*m2*M2*p1k1 + 2*l1p2*l2k1*m2*M2*p1k1 +
			        2*eg1*l1p2*l2k2*m2*M2*p1k1 - eg1*k1k2*l2p2*m2*M2*p1k1 +
			        2*l1k1*l2p2*m2*M2*p1k1 + 2*eg1*l1k2*l2p2*m2*M2*p1k1 +
			        2*l2k1*l2p2*m2*M2*p1k1 + 2*eg1*l2k2*l2p2*m2*M2*p1k1 +
			        2*eg1*l1k1*l1l2*l1p2*M2*p1k2 - 2*eg1*l1l2*l1p2*l2k1*M2*p1k2 -
			        2*eg1*l1k1*l1l2*l2p2*M2*p1k2 + 2*eg1*l1l2*l2k1*l2p2*M2*p1k2 +
			        2*eg1*l1k1*l1p2*m2*M2*p1k2 - 2*eg1*l1p2*l2k1*m2*M2*p1k2 -
			        2*eg1*l1k1*l2p2*m2*M2*p1k2 + 2*eg1*l2k1*l2p2*m2*M2*p1k2 -
			        3*eg1*l1k1*m2*M2*p1k1*p1k2 + 6*eg1*l1l2*m2*M2*p1k1*p1k2 +
			        3*eg1*l2k1*m2*M2*p1k1*p1k2 + 6*eg1*m4*M2*p1k1*p1k2 -
			        4*l1k1*l1p2*l2k1*M2*p1l1 - 4*l1l2*l1p2*l2k1*M2*p1l1 -
			        4*eg1*l1k1*l1p2*l2k2*M2*p1l1 - 4*eg1*l1l2*l1p2*l2k2*M2*p1l1 +
			        4*eg1*l1p2*l2k1*l2k2*M2*p1l1 + 2*eg1*l1k1*l1k2*l2p2*M2*p1l1 +
			        6*l1k1*l1l2*l2p2*M2*p1l1 + 2*eg1*l1k2*l1l2*l2p2*M2*p1l1 -
			        4*l1k1*l2k1*l2p2*M2*p1l1 - 2*eg1*l1k2*l2k1*l2p2*M2*p1l1 -
			        6*l1l2*l2k1*l2p2*M2*p1l1 - 2*eg1*l1k1*l2k2*l2p2*M2*p1l1 -
			        2*eg1*l1l2*l2k2*l2p2*M2*p1l1 + 2*eg1*l2k1*l2k2*l2p2*M2*p1l1 -
			        4*eg1*k1k2*l1p2*m2*M2*p1l1 - 4*eg1*k1k2*l2p2*m2*M2*p1l1 +
			        2*l1k1*l2p2*m2*M2*p1l1 - 2*l2k1*l2p2*m2*M2*p1l1 -
			        6*l1k1*l1l2*M2*p1k1*p1l1 + 6*l1l2*l2k1*M2*p1k1*p1l1 -
			        3*eg1*k1k2*m2*M2*p1k1*p1l1 - 6*l1k1*m2*M2*p1k1*p1l1 -
			        6*eg1*l1k2*m2*M2*p1k1*p1l1 - 6*l2k1*m2*M2*p1k1*p1l1 -
			        6*eg1*l2k2*m2*M2*p1k1*p1l1 - 6*eg1*l1k1*l1l2*M2*p1k2*p1l1 +
			        6*eg1*l1l2*l2k1*M2*p1k2*p1l1 - 6*eg1*l1k1*m2*M2*p1k2*p1l1 +
			        6*eg1*l2k1*m2*M2*p1k2*p1l1 + 2*eg1*l1k1*l1k2*l1p2*M2*p1l2 +
			        6*l1k1*l1l2*l1p2*M2*p1l2 + 2*eg1*l1k2*l1l2*l1p2*M2*p1l2 -
			        4*l1k1*l1p2*l2k1*M2*p1l2 - 2*eg1*l1k2*l1p2*l2k1*M2*p1l2 -
			        6*l1l2*l1p2*l2k1*M2*p1l2 - 2*eg1*l1k1*l1p2*l2k2*M2*p1l2 -
			        2*eg1*l1l2*l1p2*l2k2*M2*p1l2 + 2*eg1*l1p2*l2k1*l2k2*M2*p1l2 +
			        4*eg1*l1k1*l1k2*l2p2*M2*p1l2 + 4*l1k1*l1l2*l2p2*M2*p1l2 +
			        4*eg1*l1k2*l1l2*l2p2*M2*p1l2 - 4*l1k1*l2k1*l2p2*M2*p1l2 -
			        4*eg1*l1k2*l2k1*l2p2*M2*p1l2 - 4*eg1*k1k2*l1p2*m2*M2*p1l2 +
			        2*l1k1*l1p2*m2*M2*p1l2 - 2*l1p2*l2k1*m2*M2*p1l2 -
			        4*eg1*k1k2*l2p2*m2*M2*p1l2 + 6*l1k1*l1l2*M2*p1k1*p1l2 -
			        6*l1l2*l2k1*M2*p1k1*p1l2 + 3*eg1*k1k2*m2*M2*p1k1*p1l2 -
			        6*l1k1*m2*M2*p1k1*p1l2 - 6*eg1*l1k2*m2*M2*p1k1*p1l2 -
			        6*l2k1*m2*M2*p1k1*p1l2 - 6*eg1*l2k2*m2*M2*p1k1*p1l2 +
			        6*eg1*l1k1*l1l2*M2*p1k2*p1l2 - 6*eg1*l1l2*l2k1*M2*p1k2*p1l2 +
			        6*eg1*l1k1*m2*M2*p1k2*p1l2 - 6*eg1*l2k1*m2*M2*p1k2*p1l2 -
			        6*eg1*l1k1*l1k2*M2*p1l1*p1l2 - 18*l1k1*l1l2*M2*p1l1*p1l2 -
			        6*eg1*l1k2*l1l2*M2*p1l1*p1l2 + 12*l1k1*l2k1*M2*p1l1*p1l2 +
			        6*eg1*l1k2*l2k1*M2*p1l1*p1l2 + 18*l1l2*l2k1*M2*p1l1*p1l2 +
			        6*eg1*l1k1*l2k2*M2*p1l1*p1l2 + 6*eg1*l1l2*l2k2*M2*p1l1*p1l2 -
			        6*eg1*l2k1*l2k2*M2*p1l1*p1l2 + 12*eg1*k1k2*m2*M2*p1l1*p1l2 -
			        6*l1k1*m2*M2*p1l1*p1l2 + 6*l2k1*m2*M2*p1l1*p1l2 +
			        2*eg1*l1k1*l1k2*l1p2*l2p2*p1p2 + 6*l1k1*l1l2*l1p2*l2p2*p1p2 +
			        2*eg1*l1k2*l1l2*l1p2*l2p2*p1p2 - 4*l1k1*l1p2*l2k1*l2p2*p1p2 -
			        2*eg1*l1k2*l1p2*l2k1*l2p2*p1p2 - 6*l1l2*l1p2*l2k1*l2p2*p1p2 -
			        2*eg1*l1k1*l1p2*l2k2*l2p2*p1p2 - 2*eg1*l1l2*l1p2*l2k2*l2p2*p1p2 +
			        2*eg1*l1p2*l2k1*l2k2*l2p2*p1p2 - 4*eg1*k1k2*l1p2*l2p2*m2*p1p2 +
			        2*l1k1*l1p2*l2p2*m2*p1p2 - 2*l1p2*l2k1*l2p2*m2*p1p2 -
			        24*eg1*k1k2*l1k1*l1l2*M2*p1p2 + 8*eg1*l1k1*l1k2*l1l2*M2*p1p2 +
			        24*eg1*k1k2*l1l2*l2k1*M2*p1p2 - 16*l1k1*l1l2*l2k1*M2*p1p2 -
			        8*eg1*l1k2*l1l2*l2k1*M2*p1p2 - 8*eg1*l1k1*l1l2*l2k2*M2*p1p2 +
			        8*eg1*l1l2*l2k1*l2k2*M2*p1p2 + 4*eg1*k1k2*l1k1*m2*M2*p1p2 -
			        12*eg1*l1k1*l1k2*m2*M2*p1p2 + 24*eg1*k1k2*l1l2*m2*M2*p1p2 -
			        28*l1k1*l1l2*m2*M2*p1p2 - 8*eg1*l1k2*l1l2*m2*M2*p1p2 -
			        4*eg1*k1k2*l2k1*m2*M2*p1p2 - 4*eg1*l1k2*l2k1*m2*M2*p1p2 +
			        28*l1l2*l2k1*m2*M2*p1p2 - 4*eg1*l1k1*l2k2*m2*M2*p1p2 +
			        8*eg1*l1l2*l2k2*m2*M2*p1p2 - 12*eg1*l2k1*l2k2*m2*M2*p1p2 +
			        24*eg1*k1k2*m4*M2*p1p2 - 12*l1k1*m4*M2*p1p2 +
			        12*l2k1*m4*M2*p1p2 + 2*l1k1*l1l2*l1p2*p1k1*p1p2 -
			        2*l1l2*l1p2*l2k1*p1k1*p1p2 - 2*l1k1*l1l2*l2p2*p1k1*p1p2 +
			        2*l1l2*l2k1*l2p2*p1k1*p1p2 + eg1*k1k2*l1p2*m2*p1k1*p1p2 +
			        2*l1k1*l1p2*m2*p1k1*p1p2 + 2*eg1*l1k2*l1p2*m2*p1k1*p1p2 +
			        2*l1p2*l2k1*m2*p1k1*p1p2 + 2*eg1*l1p2*l2k2*m2*p1k1*p1p2 -
			        eg1*k1k2*l2p2*m2*p1k1*p1p2 + 2*l1k1*l2p2*m2*p1k1*p1p2 +
			        2*eg1*l1k2*l2p2*m2*p1k1*p1p2 + 2*l2k1*l2p2*m2*p1k1*p1p2 +
			        2*eg1*l2k2*l2p2*m2*p1k1*p1p2 + 2*eg1*l1k1*l1l2*l1p2*p1k2*p1p2 -
			        2*eg1*l1l2*l1p2*l2k1*p1k2*p1p2 - 2*eg1*l1k1*l1l2*l2p2*p1k2*p1p2 +
			        2*eg1*l1l2*l2k1*l2p2*p1k2*p1p2 + 2*eg1*l1k1*l1p2*m2*p1k2*p1p2 -
			        2*eg1*l1p2*l2k1*m2*p1k2*p1p2 - 2*eg1*l1k1*l2p2*m2*p1k2*p1p2 +
			        2*eg1*l2k1*l2p2*m2*p1k2*p1p2 + eg1*l1k1*m2*p1k1*p1k2*p1p2 -
			        2*eg1*l1l2*m2*p1k1*p1k2*p1p2 - eg1*l2k1*m2*p1k1*p1k2*p1p2 -
			        2*eg1*m4*p1k1*p1k2*p1p2 - 4*l1k1*l1p2*l2k1*p1l1*p1p2 -
			        4*l1l2*l1p2*l2k1*p1l1*p1p2 - 4*eg1*l1k1*l1p2*l2k2*p1l1*p1p2 -
			        4*eg1*l1l2*l1p2*l2k2*p1l1*p1p2 + 4*eg1*l1p2*l2k1*l2k2*p1l1*p1p2 +
			        2*eg1*l1k1*l1k2*l2p2*p1l1*p1p2 + 6*l1k1*l1l2*l2p2*p1l1*p1p2 +
			        2*eg1*l1k2*l1l2*l2p2*p1l1*p1p2 - 4*l1k1*l2k1*l2p2*p1l1*p1p2 -
			        2*eg1*l1k2*l2k1*l2p2*p1l1*p1p2 - 6*l1l2*l2k1*l2p2*p1l1*p1p2 -
			        2*eg1*l1k1*l2k2*l2p2*p1l1*p1p2 - 2*eg1*l1l2*l2k2*l2p2*p1l1*p1p2 +
			        2*eg1*l2k1*l2k2*l2p2*p1l1*p1p2 - 4*eg1*k1k2*l1p2*m2*p1l1*p1p2 -
			        4*eg1*k1k2*l2p2*m2*p1l1*p1p2 + 2*l1k1*l2p2*m2*p1l1*p1p2 -
			        2*l2k1*l2p2*m2*p1l1*p1p2 + 2*l1k1*l1l2*p1k1*p1l1*p1p2 -
			        2*l1l2*l2k1*p1k1*p1l1*p1p2 + eg1*k1k2*m2*p1k1*p1l1*p1p2 +
			        2*l1k1*m2*p1k1*p1l1*p1p2 + 2*eg1*l1k2*m2*p1k1*p1l1*p1p2 +
			        2*l2k1*m2*p1k1*p1l1*p1p2 + 2*eg1*l2k2*m2*p1k1*p1l1*p1p2 +
			        2*eg1*l1k1*l1l2*p1k2*p1l1*p1p2 - 2*eg1*l1l2*l2k1*p1k2*p1l1*p1p2 +
			        2*eg1*l1k1*m2*p1k2*p1l1*p1p2 - 2*eg1*l2k1*m2*p1k2*p1l1*p1p2 +
			        2*eg1*l1k1*l1k2*l1p2*p1l2*p1p2 + 6*l1k1*l1l2*l1p2*p1l2*p1p2 +
			        2*eg1*l1k2*l1l2*l1p2*p1l2*p1p2 - 4*l1k1*l1p2*l2k1*p1l2*p1p2 -
			        2*eg1*l1k2*l1p2*l2k1*p1l2*p1p2 - 6*l1l2*l1p2*l2k1*p1l2*p1p2 -
			        2*eg1*l1k1*l1p2*l2k2*p1l2*p1p2 - 2*eg1*l1l2*l1p2*l2k2*p1l2*p1p2 +
			        2*eg1*l1p2*l2k1*l2k2*p1l2*p1p2 + 4*eg1*l1k1*l1k2*l2p2*p1l2*p1p2 +
			        4*l1k1*l1l2*l2p2*p1l2*p1p2 + 4*eg1*l1k2*l1l2*l2p2*p1l2*p1p2 -
			        4*l1k1*l2k1*l2p2*p1l2*p1p2 - 4*eg1*l1k2*l2k1*l2p2*p1l2*p1p2 -
			        4*eg1*k1k2*l1p2*m2*p1l2*p1p2 + 2*l1k1*l1p2*m2*p1l2*p1p2 -
			        2*l1p2*l2k1*m2*p1l2*p1p2 - 4*eg1*k1k2*l2p2*m2*p1l2*p1p2 -
			        2*l1k1*l1l2*p1k1*p1l2*p1p2 + 2*l1l2*l2k1*p1k1*p1l2*p1p2 -
			        eg1*k1k2*m2*p1k1*p1l2*p1p2 + 2*l1k1*m2*p1k1*p1l2*p1p2 +
			        2*eg1*l1k2*m2*p1k1*p1l2*p1p2 + 2*l2k1*m2*p1k1*p1l2*p1p2 +
			        2*eg1*l2k2*m2*p1k1*p1l2*p1p2 - 2*eg1*l1k1*l1l2*p1k2*p1l2*p1p2 +
			        2*eg1*l1l2*l2k1*p1k2*p1l2*p1p2 - 2*eg1*l1k1*m2*p1k2*p1l2*p1p2 +
			        2*eg1*l2k1*m2*p1k2*p1l2*p1p2 + 2*eg1*l1k1*l1k2*p1l1*p1l2*p1p2 +
			        6*l1k1*l1l2*p1l1*p1l2*p1p2 + 2*eg1*l1k2*l1l2*p1l1*p1l2*p1p2 -
			        4*l1k1*l2k1*p1l1*p1l2*p1p2 - 2*eg1*l1k2*l2k1*p1l1*p1l2*p1p2 -
			        6*l1l2*l2k1*p1l1*p1l2*p1p2 - 2*eg1*l1k1*l2k2*p1l1*p1l2*p1p2 -
			        2*eg1*l1l2*l2k2*p1l1*p1l2*p1p2 + 2*eg1*l2k1*l2k2*p1l1*p1l2*p1p2 -
			        4*eg1*k1k2*m2*p1l1*p1l2*p1p2 + 2*l1k1*m2*p1l1*p1l2*p1p2 -
			        2*l2k1*m2*p1l1*p1l2*p1p2 - 6*eg1*l1k1*l1l2*l1p2*M2*p2k2 +
			        6*eg1*l1l2*l1p2*l2k1*M2*p2k2 + 6*eg1*l1k1*l1l2*l2p2*M2*p2k2 -
			        6*eg1*l1l2*l2k1*l2p2*M2*p2k2 - 6*eg1*l1k1*l1p2*m2*M2*p2k2 +
			        6*eg1*l1p2*l2k1*m2*M2*p2k2 + 6*eg1*l1k1*l2p2*m2*M2*p2k2 -
			        6*eg1*l2k1*l2p2*m2*M2*p2k2 + eg1*l1k1*m2*M2*p1k1*p2k2 -
			        2*eg1*l1l2*m2*M2*p1k1*p2k2 - eg1*l2k1*m2*M2*p1k1*p2k2 -
			        2*eg1*m4*M2*p1k1*p2k2 + 2*eg1*l1k1*l1l2*M2*p1l1*p2k2 -
			        2*eg1*l1l2*l2k1*M2*p1l1*p2k2 + 2*eg1*l1k1*m2*M2*p1l1*p2k2 -
			        2*eg1*l2k1*m2*M2*p1l1*p2k2 - 2*eg1*l1k1*l1l2*M2*p1l2*p2k2 +
			        2*eg1*l1l2*l2k1*M2*p1l2*p2k2 - 2*eg1*l1k1*m2*M2*p1l2*p2k2 +
			        2*eg1*l2k1*m2*M2*p1l2*p2k2 + 2*eg1*l1k1*l1l2*l1p2*p1p2*p2k2 -
			        2*eg1*l1l2*l1p2*l2k1*p1p2*p2k2 - 2*eg1*l1k1*l1l2*l2p2*p1p2*p2k2 +
			        2*eg1*l1l2*l2k1*l2p2*p1p2*p2k2 + 2*eg1*l1k1*l1p2*m2*p1p2*p2k2 -
			        2*eg1*l1p2*l2k1*m2*p1p2*p2k2 - 2*eg1*l1k1*l2p2*m2*p1p2*p2k2 +
			        2*eg1*l2k1*l2p2*m2*p1p2*p2k2 + eg1*l1k1*m2*p1k1*p1p2*p2k2 -
			        2*eg1*l1l2*m2*p1k1*p1p2*p2k2 - eg1*l2k1*m2*p1k1*p1p2*p2k2 -
			        2*eg1*m4*p1k1*p1p2*p2k2 + 2*eg1*l1k1*l1l2*p1l1*p1p2*p2k2 -
			        2*eg1*l1l2*l2k1*p1l1*p1p2*p2k2 + 2*eg1*l1k1*m2*p1l1*p1p2*p2k2 -
			        2*eg1*l2k1*m2*p1l1*p1p2*p2k2 - 2*eg1*l1k1*l1l2*p1l2*p1p2*p2k2 +
			        2*eg1*l1l2*l2k1*p1l2*p1p2*p2k2 - 2*eg1*l1k1*m2*p1l2*p1p2*p2k2 +
			        2*eg1*l2k1*m2*p1l2*p1p2*p2k2 - k1k2*l1p2*m2*M2*p1k2*pow(eg1,2) +
			        k1k2*l2p2*m2*M2*p1k2*pow(eg1,2) +
			        6*k1k2*m2*M2*p1k1*p1k2*pow(eg1,2) -
			        3*l1k2*m2*M2*p1k1*p1k2*pow(eg1,2) +
			        3*l2k2*m2*M2*p1k1*p1k2*pow(eg1,2) +
			        3*k1k2*m2*M2*p1k2*p1l1*pow(eg1,2) -
			        3*k1k2*m2*M2*p1k2*p1l2*pow(eg1,2) -
			        k1k2*l1p2*m2*p1k2*p1p2*pow(eg1,2) +
			        k1k2*l2p2*m2*p1k2*p1p2*pow(eg1,2) -
			        2*k1k2*m2*p1k1*p1k2*p1p2*pow(eg1,2) +
			        l1k2*m2*p1k1*p1k2*p1p2*pow(eg1,2) -
			        l2k2*m2*p1k1*p1k2*p1p2*pow(eg1,2) -
			        k1k2*m2*p1k2*p1l1*p1p2*pow(eg1,2) +
			        k1k2*m2*p1k2*p1l2*p1p2*pow(eg1,2) +
			        3*k1k2*l1p2*m2*M2*p2k2*pow(eg1,2) -
			        3*k1k2*l2p2*m2*M2*p2k2*pow(eg1,2) -
			        2*k1k2*m2*M2*p1k1*p2k2*pow(eg1,2) +
			        l1k2*m2*M2*p1k1*p2k2*pow(eg1,2) -
			        l2k2*m2*M2*p1k1*p2k2*pow(eg1,2) -
			        2*l1k1*m2*M2*p1k2*p2k2*pow(eg1,2) -
			        4*l1l2*m2*M2*p1k2*p2k2*pow(eg1,2) +
			        2*l2k1*m2*M2*p1k2*p2k2*pow(eg1,2) -
			        k1k2*m2*M2*p1l1*p2k2*pow(eg1,2) +
			        k1k2*m2*M2*p1l2*p2k2*pow(eg1,2) -
			        k1k2*l1p2*m2*p1p2*p2k2*pow(eg1,2) +
			        k1k2*l2p2*m2*p1p2*p2k2*pow(eg1,2) -
			        2*k1k2*m2*p1k1*p1p2*p2k2*pow(eg1,2) +
			        l1k2*m2*p1k1*p1p2*p2k2*pow(eg1,2) -
			        l2k2*m2*p1k1*p1p2*p2k2*pow(eg1,2) -
			        2*l1k1*m2*p1k2*p1p2*p2k2*pow(eg1,2) -
			        4*l1l2*m2*p1k2*p1p2*p2k2*pow(eg1,2) +
			        2*l2k1*m2*p1k2*p1p2*p2k2*pow(eg1,2) -
			        k1k2*m2*p1l1*p1p2*p2k2*pow(eg1,2) +
			        k1k2*m2*p1l2*p1p2*p2k2*pow(eg1,2) -
			        4*m2*M2*p1p2*pow(eg1,2)*pow(k1k2,2) +
			        (2*l1l2 + eg1*(l1k2 - l2k2))*m2*(3*M2 - p1p2)*pow(k1p2,2) -
			        6*l1p2*l2p2*M2*pow(l1k1,2) + 2*l2p2*M2*p1l1*pow(l1k1,2) +
			        2*l1p2*M2*p1l2*pow(l1k1,2) + 4*l2p2*M2*p1l2*pow(l1k1,2) -
			        6*M2*p1l1*p1l2*pow(l1k1,2) + 2*l1p2*l2p2*p1p2*pow(l1k1,2) +
			        8*l1l2*M2*p1p2*pow(l1k1,2) - 16*m2*M2*p1p2*pow(l1k1,2) +
			        2*l2p2*p1l1*p1p2*pow(l1k1,2) + 2*l1p2*p1l2*p1p2*pow(l1k1,2) +
			        4*l2p2*p1l2*p1p2*pow(l1k1,2) + 2*p1l1*p1l2*p1p2*pow(l1k1,2) -
			        12*l1p2*l2p2*M2*pow(l1l2,2) + 2*l1p2*M2*p1k1*pow(l1l2,2) -
			        2*l2p2*M2*p1k1*pow(l1l2,2) + 2*eg1*l1p2*M2*p1k2*pow(l1l2,2) -
			        2*eg1*l2p2*M2*p1k2*pow(l1l2,2) + 4*l2p2*M2*p1l1*pow(l1l2,2) -
			        6*M2*p1k1*p1l1*pow(l1l2,2) - 6*eg1*M2*p1k2*p1l1*pow(l1l2,2) +
			        4*l1p2*M2*p1l2*pow(l1l2,2) + 6*M2*p1k1*p1l2*pow(l1l2,2) +
			        6*eg1*M2*p1k2*p1l2*pow(l1l2,2) - 12*M2*p1l1*p1l2*pow(l1l2,2) +
			        4*l1p2*l2p2*p1p2*pow(l1l2,2) - 24*eg1*k1k2*M2*p1p2*pow(l1l2,2) +
			        16*l1k1*M2*p1p2*pow(l1l2,2) + 8*eg1*l1k2*M2*p1p2*pow(l1l2,2) -
			        16*l2k1*M2*p1p2*pow(l1l2,2) - 8*eg1*l2k2*M2*p1p2*pow(l1l2,2) -
			        24*m2*M2*p1p2*pow(l1l2,2) + 2*l1p2*p1k1*p1p2*pow(l1l2,2) -
			        2*l2p2*p1k1*p1p2*pow(l1l2,2) + 2*eg1*l1p2*p1k2*p1p2*pow(l1l2,2) -
			        2*eg1*l2p2*p1k2*p1p2*pow(l1l2,2) + 4*l2p2*p1l1*p1p2*pow(l1l2,2) +
			        2*p1k1*p1l1*p1p2*pow(l1l2,2) + 2*eg1*p1k2*p1l1*p1p2*pow(l1l2,2) +
			        4*l1p2*p1l2*p1p2*pow(l1l2,2) - 2*p1k1*p1l2*p1p2*pow(l1l2,2) -
			        2*eg1*p1k2*p1l2*p1p2*pow(l1l2,2) + 4*p1l1*p1l2*p1p2*pow(l1l2,2) -
			        6*eg1*l1p2*M2*p2k2*pow(l1l2,2) +
			        6*eg1*l2p2*M2*p2k2*pow(l1l2,2) +
			        2*eg1*M2*p1l1*p2k2*pow(l1l2,2) -
			        2*eg1*M2*p1l2*p2k2*pow(l1l2,2) +
			        2*eg1*l1p2*p1p2*p2k2*pow(l1l2,2) -
			        2*eg1*l2p2*p1p2*p2k2*pow(l1l2,2) +
			        2*eg1*p1l1*p1p2*p2k2*pow(l1l2,2) -
			        2*eg1*p1l2*p1p2*p2k2*pow(l1l2,2) +
			        k1p2*(l1k1*(-2*l1l2*(l1p2*(3*M2 - p1p2) +
			                 l2p2*(-3*M2 + p1p2) - (p1l1 - p1l2)*(M2 + p1p2)) +
			              m2*(-6*l1p2*M2 - 6*l2p2*M2 + eg1*M2*p1k2 +
			                 2*M2*p1l1 + 2*M2*p1l2 + 2*l1p2*p1p2 + 2*l2p2*p1p2 +
			                 eg1*p1k2*p1p2 + 2*p1l1*p1p2 + 2*p1l2*p1p2 -
			                 3*eg1*M2*p2k2 + eg1*p1p2*p2k2)) -
			           2*l1l2*(l1p2*l2k1*(-3*M2 + p1p2) +
			              l2k1*(l2p2*(3*M2 - p1p2) + (p1l1 - p1l2)*(M2 + p1p2)) +
			              m2*(M2*(2*p1k1 + eg1*(p1k2 - 3*p2k2)) +
			                 p1p2*(2*p1k1 + eg1*(p1k2 + p2k2)))) +
			           m2*(2*l2k1*(l1p2*(-3*M2 + p1p2) + l2p2*(-3*M2 + p1p2) +
			                 (p1l1 + p1l2)*(M2 + p1p2)) +
			              eg1*(-6*l1p2*l2k2*M2 - 6*l2k2*l2p2*M2 +
			                 2*l2k2*M2*p1k1 - l2k1*M2*p1k2 - 2*m2*M2*p1k2 +
			                 2*l2k2*M2*p1l1 + 2*l2k2*M2*p1l2 + 2*l1p2*l2k2*p1p2 +
			                 2*l2k2*l2p2*p1p2 + 2*l2k2*p1k1*p1p2 - l2k1*p1k2*p1p2 -
			                 2*m2*p1k2*p1p2 + 2*l2k2*p1l1*p1p2 + 2*l2k2*p1l2*p1p2 +
			                 2*l1k2*(l1p2*(-3*M2 + p1p2) + l2p2*(-3*M2 + p1p2) -
			                    (p1k1 - p1l1 - p1l2)*(M2 + p1p2)) +
			                 k1k2*(l2p2*(3*M2 - p1p2) + l1p2*(-3*M2 + p1p2) +
			                    (p1l1 - p1l2)*(M2 + p1p2)) + 3*l2k1*M2*p2k2 +
			                 6*m2*M2*p2k2 - l2k1*p1p2*p2k2 - 2*m2*p1p2*p2k2) -
			              (2*k1k2 - l1k2 + l2k2)*
			               (M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2))*pow(eg1,2)) -
			           2*(l1p2*(3*M2 - p1p2) + l2p2*(-3*M2 + p1p2) -
			              (p1l1 - p1l2)*(M2 + p1p2))*pow(l1l2,2)) +
			        8*M2*p1p2*pow(l1l2,3) + 6*l1k1*l2k1*M2*pow(l1p2,2) +
			        6*l1l2*l2k1*M2*pow(l1p2,2) + 6*eg1*l1k1*l2k2*M2*pow(l1p2,2) +
			        6*eg1*l1l2*l2k2*M2*pow(l1p2,2) -
			        6*eg1*l2k1*l2k2*M2*pow(l1p2,2) + 6*eg1*k1k2*m2*M2*pow(l1p2,2) -
			        2*l1k1*l2k1*p1p2*pow(l1p2,2) - 2*l1l2*l2k1*p1p2*pow(l1p2,2) -
			        2*eg1*l1k1*l2k2*p1p2*pow(l1p2,2) -
			        2*eg1*l1l2*l2k2*p1p2*pow(l1p2,2) +
			        2*eg1*l2k1*l2k2*p1p2*pow(l1p2,2) - 2*eg1*k1k2*m2*p1p2*pow(l1p2,2) -
			        6*l1p2*l2p2*M2*pow(l2k1,2) + 4*l1p2*M2*p1l1*pow(l2k1,2) +
			        2*l2p2*M2*p1l1*pow(l2k1,2) + 2*l1p2*M2*p1l2*pow(l2k1,2) -
			        6*M2*p1l1*p1l2*pow(l2k1,2) + 2*l1p2*l2p2*p1p2*pow(l2k1,2) +
			        8*l1l2*M2*p1p2*pow(l2k1,2) - 16*m2*M2*p1p2*pow(l2k1,2) +
			        4*l1p2*p1l1*p1p2*pow(l2k1,2) + 2*l2p2*p1l1*p1p2*pow(l2k1,2) +
			        2*l1p2*p1l2*p1p2*pow(l2k1,2) + 2*p1l1*p1l2*p1p2*pow(l2k1,2) -
			        6*M2*pow(l1p2,2)*pow(l2k1,2) + 2*p1p2*pow(l1p2,2)*pow(l2k1,2) -
			        6*eg1*l1k1*l1k2*M2*pow(l2p2,2) - 6*l1k1*l1l2*M2*pow(l2p2,2) -
			        6*eg1*l1k2*l1l2*M2*pow(l2p2,2) + 6*l1k1*l2k1*M2*pow(l2p2,2) +
			        6*eg1*l1k2*l2k1*M2*pow(l2p2,2) + 6*eg1*k1k2*m2*M2*pow(l2p2,2) +
			        2*eg1*l1k1*l1k2*p1p2*pow(l2p2,2) + 2*l1k1*l1l2*p1p2*pow(l2p2,2) +
			        2*eg1*l1k2*l1l2*p1p2*pow(l2p2,2) - 2*l1k1*l2k1*p1p2*pow(l2p2,2) -
			        2*eg1*l1k2*l2k1*p1p2*pow(l2p2,2) - 2*eg1*k1k2*m2*p1p2*pow(l2p2,2) -
			        6*M2*pow(l1k1,2)*pow(l2p2,2) + 2*p1p2*pow(l1k1,2)*pow(l2p2,2) +
			        20*eg1*k1k2*l1k1*l1l2*M4 - 4*eg1*l1k1*l1k2*l1l2*M4 -
			        20*eg1*k1k2*l1l2*l2k1*M4 + 8*l1k1*l1l2*l2k1*M4 +
			        4*eg1*l1k2*l1l2*l2k1*M4 + 4*eg1*l1k1*l1l2*l2k2*M4 -
			        4*eg1*l1l2*l2k1*l2k2*M4 - 2*eg1*k1k2*l1k1*m2*M4 +
			        14*eg1*l1k1*l1k2*m2*M4 - 24*eg1*k1k2*l1l2*m2*M4 +
			        26*l1k1*l1l2*m2*M4 + 8*eg1*l1k2*l1l2*m2*M4 +
			        2*eg1*k1k2*l2k1*m2*M4 + 2*eg1*l1k2*l2k1*m2*M4 -
			        26*l1l2*l2k1*m2*M4 + 2*eg1*l1k1*l2k2*m2*M4 -
			        8*eg1*l1l2*l2k2*m2*M4 + 14*eg1*l2k1*l2k2*m2*M4 -
			        24*eg1*k1k2*m4*M4 + 10*l1k1*m4*M4 -
			        10*l2k1*m4*M4 + 2*m2*pow(eg1,2)*pow(k1k2,2)*M4 -
			        4*l1l2*pow(l1k1,2)*M4 + 16*m2*pow(l1k1,2)*M4 +
			        20*eg1*k1k2*pow(l1l2,2)*M4 - 8*l1k1*pow(l1l2,2)*M4 -
			        4*eg1*l1k2*pow(l1l2,2)*M4 + 8*l2k1*pow(l1l2,2)*M4 +
			        4*eg1*l2k2*pow(l1l2,2)*M4 + 20*m2*pow(l1l2,2)*M4 -
			        4*pow(l1l2,3)*M4 - 4*l1l2*pow(l2k1,2)*M4 +
			        16*m2*pow(l2k1,2)*M4 + 3*eg1*l1k2*m2*M2*pow(p1k1,2) +
			        6*l1l2*m2*M2*pow(p1k1,2) - 3*eg1*l2k2*m2*M2*pow(p1k1,2) -
			        eg1*l1k2*m2*p1p2*pow(p1k1,2) - 2*l1l2*m2*p1p2*pow(p1k1,2) +
			        eg1*l2k2*m2*p1p2*pow(p1k1,2) +
			        3*l1k1*m2*M2*pow(eg1,2)*pow(p1k2,2) +
			        6*l1l2*m2*M2*pow(eg1,2)*pow(p1k2,2) -
			        3*l2k1*m2*M2*pow(eg1,2)*pow(p1k2,2) -
			        l1k1*m2*p1p2*pow(eg1,2)*pow(p1k2,2) -
			        2*l1l2*m2*p1p2*pow(eg1,2)*pow(p1k2,2) +
			        l2k1*m2*p1p2*pow(eg1,2)*pow(p1k2,2) +
			        6*l1k1*l2k1*M2*pow(p1l1,2) + 6*l1l2*l2k1*M2*pow(p1l1,2) +
			        6*eg1*l1k1*l2k2*M2*pow(p1l1,2) +
			        6*eg1*l1l2*l2k2*M2*pow(p1l1,2) -
			        6*eg1*l2k1*l2k2*M2*pow(p1l1,2) + 6*eg1*k1k2*m2*M2*pow(p1l1,2) -
			        2*l1k1*l2k1*p1p2*pow(p1l1,2) - 2*l1l2*l2k1*p1p2*pow(p1l1,2) -
			        2*eg1*l1k1*l2k2*p1p2*pow(p1l1,2) -
			        2*eg1*l1l2*l2k2*p1p2*pow(p1l1,2) +
			        2*eg1*l2k1*l2k2*p1p2*pow(p1l1,2) - 2*eg1*k1k2*m2*p1p2*pow(p1l1,2) -
			        6*M2*pow(l2k1,2)*pow(p1l1,2) + 2*p1p2*pow(l2k1,2)*pow(p1l1,2) -
			        6*eg1*l1k1*l1k2*M2*pow(p1l2,2) - 6*l1k1*l1l2*M2*pow(p1l2,2) -
			        6*eg1*l1k2*l1l2*M2*pow(p1l2,2) + 6*l1k1*l2k1*M2*pow(p1l2,2) +
			        6*eg1*l1k2*l2k1*M2*pow(p1l2,2) + 6*eg1*k1k2*m2*M2*pow(p1l2,2) +
			        2*eg1*l1k1*l1k2*p1p2*pow(p1l2,2) + 2*l1k1*l1l2*p1p2*pow(p1l2,2) +
			        2*eg1*l1k2*l1l2*p1p2*pow(p1l2,2) - 2*l1k1*l2k1*p1p2*pow(p1l2,2) -
			        2*eg1*l1k2*l2k1*p1p2*pow(p1l2,2) - 2*eg1*k1k2*m2*p1p2*pow(p1l2,2) -
			        6*M2*pow(l1k1,2)*pow(p1l2,2) + 2*p1p2*pow(l1k1,2)*pow(p1l2,2) +
			        4*eg1*k1k2*l1k1*l1l2*pow(p1p2,2) -
			        4*eg1*l1k1*l1k2*l1l2*pow(p1p2,2) -
			        4*eg1*k1k2*l1l2*l2k1*pow(p1p2,2) + 8*l1k1*l1l2*l2k1*pow(p1p2,2) +
			        4*eg1*l1k2*l1l2*l2k1*pow(p1p2,2) +
			        4*eg1*l1k1*l1l2*l2k2*pow(p1p2,2) -
			        4*eg1*l1l2*l2k1*l2k2*pow(p1p2,2) - 2*eg1*k1k2*l1k1*m2*pow(p1p2,2) -
			        2*eg1*l1k1*l1k2*m2*pow(p1p2,2) + 2*l1k1*l1l2*m2*pow(p1p2,2) +
			        2*eg1*k1k2*l2k1*m2*pow(p1p2,2) + 2*eg1*l1k2*l2k1*m2*pow(p1p2,2) -
			        2*l1l2*l2k1*m2*pow(p1p2,2) + 2*eg1*l1k1*l2k2*m2*pow(p1p2,2) -
			        2*eg1*l2k1*l2k2*m2*pow(p1p2,2) + 2*l1k1*m4*pow(p1p2,2) -
			        2*l2k1*m4*pow(p1p2,2) + 2*m2*pow(eg1,2)*pow(k1k2,2)*pow(p1p2,2) -
			        4*l1l2*pow(l1k1,2)*pow(p1p2,2) +
			        4*eg1*k1k2*pow(l1l2,2)*pow(p1p2,2) -
			        8*l1k1*pow(l1l2,2)*pow(p1p2,2) -
			        4*eg1*l1k2*pow(l1l2,2)*pow(p1p2,2) +
			        8*l2k1*pow(l1l2,2)*pow(p1p2,2) +
			        4*eg1*l2k2*pow(l1l2,2)*pow(p1p2,2) + 4*m2*pow(l1l2,2)*pow(p1p2,2) -
			        4*pow(l1l2,3)*pow(p1p2,2) - 4*l1l2*pow(l2k1,2)*pow(p1p2,2) +
			        3*l1k1*m2*M2*pow(eg1,2)*pow(p2k2,2) +
			        6*l1l2*m2*M2*pow(eg1,2)*pow(p2k2,2) -
			        3*l2k1*m2*M2*pow(eg1,2)*pow(p2k2,2) -
			        l1k1*m2*p1p2*pow(eg1,2)*pow(p2k2,2) -
			        2*l1l2*m2*p1p2*pow(eg1,2)*pow(p2k2,2) +
			        l2k1*m2*p1p2*pow(eg1,2)*pow(p2k2,2))) +
			  4*eg1*m2*pow(l1k1,-1)*pow(l1k2,-1)*pow(l1k1 + eg1*(-k1k2 + l1k2),-1)*
			   pow(l2k2,-1)/M2*(4*f1*f2*M2*
			      ((p1k2 - p2k2)*(l1p2*l2k2 + l1k1*p1k2 + l1l2*p1k2 - m2*p1k2 +
			           k1k2*(l1p2 - p1l1) - l2k2*p1l1 - l1k1*p2k2 - l1l2*p2k2 +
			           m2*p2k2) + l1k2*(4*l2k2*(M2 - p1p2) +
			           (p1k2 - p2k2)*(-k1p2 + l2p2 + p1k1 + 2*eg1*p1k2 - p1l2 -
			              2*eg1*p2k2))) + 4*M2*
			      (k1k2*l1p2*p1k2 + l1p2*l2k2*p1k2 - 2*l1k1*p1k2*p2k2 -
			        2*l1l2*p1k2*p2k2 + 2*m2*p1k2*p2k2 + k1k2*p1l1*p2k2 +
			        l2k2*p1l1*p2k2 + l1k2*(-(k1p2*p1k2) + l2p2*p1k2 +
			           2*l2k2*(M2 - p1p2) - p1k1*p2k2 - 4*eg1*p1k2*p2k2 + p1l2*p2k2))*
			      pow(f1,2) + pow(f2,2)*(l1p2*l2k2*M2*p1k2 - 3*l2k2*M2*p1k2*p1l1 +
			        l1p2*l2k2*p1k2*p1p2 + l2k2*p1k2*p1l1*p1p2 - 3*l1p2*l2k2*M2*p2k2 -
			        2*l1k1*M2*p1k2*p2k2 - 2*l1l2*M2*p1k2*p2k2 +
			        2*m2*M2*p1k2*p2k2 + l2k2*M2*p1l1*p2k2 + l1p2*l2k2*p1p2*p2k2 -
			        2*l1k1*p1k2*p1p2*p2k2 - 2*l1l2*p1k2*p1p2*p2k2 +
			        2*m2*p1k2*p1p2*p2k2 + l2k2*p1l1*p1p2*p2k2 +
			        k1k2*(l1p2*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)) +
			           p1l1*(M2*(-3*p1k2 + p2k2) + p1p2*(p1k2 + p2k2))) +
			        3*l1k1*M2*pow(p1k2,2) + 3*l1l2*M2*pow(p1k2,2) -
			        3*m2*M2*pow(p1k2,2) - l1k1*p1p2*pow(p1k2,2) -
			        l1l2*p1p2*pow(p1k2,2) + m2*p1p2*pow(p1k2,2) +
			        3*l1k1*M2*pow(p2k2,2) + 3*l1l2*M2*pow(p2k2,2) -
			        3*m2*M2*pow(p2k2,2) - l1k1*p1p2*pow(p2k2,2) -
			        l1l2*p1p2*pow(p2k2,2) + m2*p1p2*pow(p2k2,2) -
			        l1k2*(-(l2p2*M2*p1k2) - 3*M2*p1k1*p1k2 + 3*M2*p1k2*p1l2 -
			           l2p2*p1k2*p1p2 + p1k1*p1k2*p1p2 - p1k2*p1l2*p1p2 +
			           3*l2p2*M2*p2k2 + M2*p1k1*p2k2 + 4*eg1*M2*p1k2*p2k2 -
			           M2*p1l2*p2k2 - l2p2*p1p2*p2k2 + p1k1*p1p2*p2k2 +
			           4*eg1*p1k2*p1p2*p2k2 - p1l2*p1p2*p2k2 +
			           k1p2*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)) -
			           8*l2k2*(-(M2*p1p2) + M4) - 6*eg1*M2*pow(p1k2,2) +
			           2*eg1*p1p2*pow(p1k2,2) - 6*eg1*M2*pow(p2k2,2) +
			           2*eg1*p1p2*pow(p2k2,2)))) +
			  4*eg1*pow(l1k1,-1)*pow(l1k1 + eg1*(-k1k2 + l1k2),-1)*pow(l2k2,-1)*
			   pow(l2k1 + eg1*(k1k2 + l2k2),-1)/M2*
			   (4*M2*pow(f1,2)*(-4*l1k1*l1k2*l1l2*M2 - 2*k1k2*l1l2*l2k1*M2 +
			        6*l1k2*l1l2*l2k1*M2 - 2*l1k1*l1l2*l2k2*M2 +
			        2*l1k1*l2k1*l2k2*M2 - 5*k1k2*l1l2*m2*M2 +
			        8*l1k2*l1l2*m2*M2 + 2*k1k2*l2k1*m2*M2 -
			        5*l1k2*l2k1*m2*M2 + 3*l1k1*l2k2*m2*M2 +
			        2*l1l2*l2k2*m2*M2 - 2*l2k1*l2k2*m2*M2 - 3*k1k2*m4*M2 -
			        4*l2k2*m4*M2 - 2*k1k2*l1l2*l1p2*p1k1 + 2*l1k2*l1l2*l1p2*p1k1 +
			        2*l1k1*l1p2*l2k2*p1k1 - 2*l1k1*l1k2*l2p2*p1k1 -
			        2*l1k2*l1l2*l2p2*p1k1 + 2*l1k2*l2k1*l2p2*p1k1 +
			        4*k1k2*k1p2*m2*p1k1 - 4*k1p2*l1k2*m2*p1k1 + l1k2*l1p2*m2*p1k1 +
			        2*k1p2*l2k2*m2*p1k1 - 2*l1p2*l2k2*m2*p1k1 + k1k2*l2p2*m2*p1k1 +
			        l2k2*l2p2*m2*p1k1 - 2*l1k1*l1p2*l2k1*p1k2 +
			        2*l1k1*l1l2*l2p2*p1k2 - 2*l1k1*l2k1*l2p2*p1k2 +
			        k1p2*l1l2*m2*p1k2 - l1k1*l1p2*m2*p1k2 - k1p2*l2k1*m2*p1k2 +
			        l1p2*l2k1*m2*p1k2 - 3*l1k1*l2p2*m2*p1k2 + l2k1*l2p2*m2*p1k2 -
			        k1p2*m4*p1k2 - 2*k1k2*k1p2*l1l2*p1l1 + 2*k1p2*l1k2*l1l2*p1l1 +
			        4*k1k2*l1p2*l2k1*p1l1 - 4*l1k2*l1p2*l2k1*p1l1 +
			        2*k1p2*l1k1*l2k2*p1l1 - 2*k1k2*l1k1*l2p2*p1l1 +
			        2*l1k1*l1k2*l2p2*p1l1 - 2*k1k2*l1l2*l2p2*p1l1 +
			        4*l1k2*l1l2*l2p2*p1l1 + 2*k1k2*l2k1*l2p2*p1l1 -
			        4*l1k2*l2k1*l2p2*p1l1 + k1p2*l1k2*m2*p1l1 - 2*k1k2*l1p2*m2*p1l1 -
			        2*k1p2*l2k2*m2*p1l1 - 2*k1k2*l2p2*m2*p1l1 -
			        2*l2k2*l2p2*m2*p1l1 - 2*k1p2*l1k1*l1k2*p1l2 -
			        2*k1p2*l1k2*l1l2*p1l2 - 2*k1k2*l1k1*l1p2*p1l2 +
			        2*l1k1*l1k2*l1p2*p1l2 - 2*k1k2*l1l2*l1p2*p1l2 +
			        4*l1k2*l1l2*l1p2*p1l2 + 2*k1p2*l1k2*l2k1*p1l2 +
			        2*k1k2*l1p2*l2k1*p1l2 - 4*l1k2*l1p2*l2k1*p1l2 +
			        k1k2*k1p2*m2*p1l2 - 2*k1k2*l1p2*m2*p1l2 + k1p2*l2k2*m2*p1l2 -
			        2*l1p2*l2k2*m2*p1l2 - 2*k1k2*l2p2*m2*p1l2 +
			        2*k1k2*l1k1*l1l2*p1p2 + 2*l1k1*l1k2*l2k1*p1p2 +
			        3*k1k2*l1l2*m2*p1p2 - 4*l1k2*l1l2*m2*p1p2 -
			        2*k1k2*l2k1*m2*p1p2 + 3*l1k2*l2k1*m2*p1p2 + l1k1*l2k2*m2*p1p2 +
			        3*k1k2*m4*p1p2 + 2*l2k2*m4*p1p2 + l1l2*m2*p1k1*p2k2 -
			        l2k1*m2*p1k1*p2k2 - m4*p1k1*p2k2 - 2*l1k1*l2k1*p1l1*p2k2 -
			        l1k1*m2*p1l1*p2k2 + l2k1*m2*p1l1*p2k2 + 2*l1k1*l1l2*p1l2*p2k2 -
			        2*l1k1*l2k1*p1l2*p2k2 - 3*l1k1*m2*p1l2*p2k2 + l2k1*m2*p1l2*p2k2 +
			        4*(k1k2 - l1k2)*m2*p1k2*p2k2*pow(eg1,2) +
			        2*l2p2*p1k2*pow(l1k1,2) - 2*l2k2*p1p2*pow(l1k1,2) +
			        2*p1l2*p2k2*pow(l1k1,2) +
			        eg1*(2*k1p2*l1k1*l2k2*p1k2 - k1p2*l2k2*m2*p1k2 +
			           l1p2*l2k2*m2*p1k2 + 2*l1k1*l2k2*p1k1*p2k2 -
			           l2k2*m2*p1k1*p2k2 - 4*l1k1*l2k1*p1k2*p2k2 -
			           4*l1k1*m2*p1k2*p2k2 - 2*l1l2*m2*p1k2*p2k2 +
			           2*l2k1*m2*p1k2*p2k2 + 2*m4*p1k2*p2k2 + l2k2*m2*p1l1*p2k2 +
			           k1k2*(2*l1k1*l2k2*M2 - 2*k1p2*l1l2*p1k2 + 2*l1p2*l2k1*p1k2 +
			              4*k1p2*m2*p1k2 + 2*l1p2*m2*p1k2 + l2p2*m2*p1k2 -
			              2*l1k1*l2k2*p1p2 -
			              2*l1k2*(2*m2*M2 - l2p2*p1k1 + l2p2*p1l1 - k1p2*p1l2 +
			                 l1p2*p1l2 - m2*p1p2 + l1l2*(-3*M2 + p1p2) +
			                 l2k1*(M2 + p1p2)) - 2*l1l2*p1k1*p2k2 +
			              4*m2*p1k1*p2k2 + 2*l2k1*p1l1*p2k2 + 2*m2*p1l1*p2k2 +
			              m2*p1l2*p2k2) -
			           l1k2*(-2*l1l2*l1p2*p1k2 - 2*k1p2*l2k1*p1k2 +
			              4*l1p2*l2k1*p1k2 + 2*k1p2*m2*p1k2 + 2*l1p2*m2*p1k2 +
			              l2p2*m2*p1k2 -
			              2*l2k2*(m2*M2 - 2*(k1p2 - l1p2)*(p1k1 - p1l1)) -
			              2*l2k1*p1k1*p2k2 + 2*m2*p1k1*p2k2 - 2*l1l2*p1l1*p2k2 +
			              4*l2k1*p1l1*p2k2 + 2*m2*p1l1*p2k2 + m2*p1l2*p2k2 +
			              2*l1k1*(-(l2p2*p1k2) + l2k2*(M2 + p1p2) - p1l2*p2k2)) -
			           2*(l1l2 - 2*m2)*(M2 - p1p2)*pow(k1k2,2) +
			           2*(-2*l1l2*M2 + m2*M2 - l2p2*p1k1 + l2p2*p1l1 -
			              k1p2*p1l2 + l1p2*p1l2 + l2k1*(M2 + p1p2))*pow(l1k2,2)) +
			        2*k1k2*M2*pow(l1l2,2) - 4*l1k2*M2*pow(l1l2,2) -
			        2*l1k2*M2*pow(l2k1,2)) -
			     4*f1*f2*M2*(2*l1k1*l1k2*l1p2*l2p2 + 4*l1k2*l1l2*l1p2*l2p2 -
			        4*l1k2*l1p2*l2k1*l2p2 - 2*l1p2*l2k2*l2p2*m2 +
			        4*l1k1*l1k2*l1l2*M2 + 2*l1k1*l1k2*l2k1*M2 -
			        6*l1k2*l1l2*l2k1*M2 + 2*l1k1*l1l2*l2k2*M2 -
			        2*l1k1*l2k1*l2k2*M2 - 12*l1k2*l1l2*m2*M2 +
			        8*l1k2*l2k1*m2*M2 - 2*l1k1*l2k2*m2*M2 -
			        2*eg1*l1k2*l2k2*m2*M2 - 2*l1l2*l2k2*m2*M2 +
			        2*l2k1*l2k2*m2*M2 + 6*l2k2*m4*M2 - 2*l1k2*l1l2*l1p2*p1k1 -
			        2*l1k1*l1p2*l2k2*p1k1 - 4*eg1*l1k2*l1p2*l2k2*p1k1 +
			        2*l1k1*l1k2*l2p2*p1k1 + 2*l1k2*l1l2*l2p2*p1k1 -
			        2*l1k2*l2k1*l2p2*p1k1 - l1k2*l1p2*m2*p1k1 + 2*l1p2*l2k2*m2*p1k1 -
			        l2k2*l2p2*m2*p1k1 - 2*eg1*l1k2*l1l2*l1p2*p1k2 +
			        2*l1k1*l1p2*l2k1*p1k2 + 4*eg1*l1k2*l1p2*l2k1*p1k2 -
			        2*eg1*l1k1*l1k2*l2p2*p1k2 - 2*l1k1*l1l2*l2p2*p1k2 +
			        2*l1k1*l2k1*l2p2*p1k2 + l1k1*l1p2*m2*p1k2 +
			        2*eg1*l1k2*l1p2*m2*p1k2 - l1p2*l2k1*m2*p1k2 -
			        eg1*l1p2*l2k2*m2*p1k2 + 3*l1k1*l2p2*m2*p1k2 +
			        eg1*l1k2*l2p2*m2*p1k2 - l2k1*l2p2*m2*p1k2 +
			        2*eg1*l1k2*l2k1*p1k1*p1k2 + 2*eg1*l1k1*l2k2*p1k1*p1k2 -
			        2*eg1*l1k2*m2*p1k1*p1k2 + l1l2*m2*p1k1*p1k2 -
			        l2k1*m2*p1k1*p1k2 - eg1*l2k2*m2*p1k1*p1k2 - m4*p1k1*p1k2 +
			        4*l1k2*l1p2*l2k1*p1l1 + 4*eg1*l1k2*l1p2*l2k2*p1l1 -
			        2*l1k1*l1k2*l2p2*p1l1 - 4*l1k2*l1l2*l2p2*p1l1 +
			        4*l1k2*l2k1*l2p2*p1l1 + 2*l2k2*l2p2*m2*p1l1 +
			        2*l1k2*l1l2*p1k1*p1l1 + 2*l1k1*l2k2*p1k1*p1l1 +
			        4*eg1*l1k2*l2k2*p1k1*p1l1 + l1k2*m2*p1k1*p1l1 -
			        2*l2k2*m2*p1k1*p1l1 + 2*eg1*l1k2*l1l2*p1k2*p1l1 -
			        2*l1k1*l2k1*p1k2*p1l1 - 4*eg1*l1k2*l2k1*p1k2*p1l1 -
			        l1k1*m2*p1k2*p1l1 - 2*eg1*l1k2*m2*p1k2*p1l1 +
			        l2k1*m2*p1k2*p1l1 + eg1*l2k2*m2*p1k2*p1l1 -
			        2*l1k1*l1k2*l1p2*p1l2 - 4*l1k2*l1l2*l1p2*p1l2 +
			        4*l1k2*l1p2*l2k1*p1l2 + 2*l1p2*l2k2*m2*p1l2 -
			        2*l1k1*l1k2*p1k1*p1l2 - 2*l1k2*l1l2*p1k1*p1l2 +
			        2*l1k2*l2k1*p1k1*p1l2 + l2k2*m2*p1k1*p1l2 +
			        2*eg1*l1k1*l1k2*p1k2*p1l2 + 2*l1k1*l1l2*p1k2*p1l2 -
			        2*l1k1*l2k1*p1k2*p1l2 - 3*l1k1*m2*p1k2*p1l2 -
			        eg1*l1k2*m2*p1k2*p1l2 + l2k1*m2*p1k2*p1l2 +
			        2*l1k1*l1k2*p1l1*p1l2 + 4*l1k2*l1l2*p1l1*p1l2 -
			        4*l1k2*l2k1*p1l1*p1l2 - 2*l2k2*m2*p1l1*p1l2 -
			        4*l1k1*l1k2*l1l2*p1p2 - 2*l1k1*l1k2*l2k1*p1p2 +
			        6*l1k2*l1l2*l2k1*p1p2 - 2*l1k1*l1l2*l2k2*p1p2 +
			        2*l1k1*l2k1*l2k2*p1p2 + 12*l1k2*l1l2*m2*p1p2 -
			        8*l1k2*l2k1*m2*p1p2 + 2*l1k1*l2k2*m2*p1p2 +
			        2*eg1*l1k2*l2k2*m2*p1p2 + 2*l1l2*l2k2*m2*p1p2 -
			        2*l2k1*l2k2*m2*p1p2 - 6*l2k2*m4*p1p2 +
			        2*eg1*l1k2*l1l2*l1p2*p2k2 - 2*l1k1*l1p2*l2k1*p2k2 -
			        4*eg1*l1k2*l1p2*l2k1*p2k2 + 2*eg1*l1k1*l1k2*l2p2*p2k2 +
			        2*l1k1*l1l2*l2p2*p2k2 - 2*l1k1*l2k1*l2p2*p2k2 -
			        l1k1*l1p2*m2*p2k2 - 2*eg1*l1k2*l1p2*m2*p2k2 +
			        l1p2*l2k1*m2*p2k2 + eg1*l1p2*l2k2*m2*p2k2 -
			        3*l1k1*l2p2*m2*p2k2 - eg1*l1k2*l2p2*m2*p2k2 +
			        l2k1*l2p2*m2*p2k2 - 2*eg1*l1k2*l2k1*p1k1*p2k2 -
			        2*eg1*l1k1*l2k2*p1k1*p2k2 + 2*eg1*l1k2*m2*p1k1*p2k2 -
			        l1l2*m2*p1k1*p2k2 + l2k1*m2*p1k1*p2k2 + eg1*l2k2*m2*p1k1*p2k2 +
			        m4*p1k1*p2k2 + 4*eg1*l1k1*l2k1*p1k2*p2k2 +
			        4*eg1*l1k1*m2*p1k2*p2k2 + 2*eg1*l1l2*m2*p1k2*p2k2 -
			        2*eg1*l2k1*m2*p1k2*p2k2 - 2*eg1*m4*p1k2*p2k2 -
			        2*eg1*l1k2*l1l2*p1l1*p2k2 + 2*l1k1*l2k1*p1l1*p2k2 +
			        4*eg1*l1k2*l2k1*p1l1*p2k2 + l1k1*m2*p1l1*p2k2 +
			        2*eg1*l1k2*m2*p1l1*p2k2 - l2k1*m2*p1l1*p2k2 -
			        eg1*l2k2*m2*p1l1*p2k2 - 2*eg1*l1k1*l1k2*p1l2*p2k2 -
			        2*l1k1*l1l2*p1l2*p2k2 + 2*l1k1*l2k1*p1l2*p2k2 +
			        3*l1k1*m2*p1l2*p2k2 + eg1*l1k2*m2*p1l2*p2k2 -
			        l2k1*m2*p1l2*p2k2 + 4*l1k2*m2*p1k2*p2k2*pow(eg1,2) +
			        4*eg1*(l1l2 - 2*m2)*(M2 - p1p2)*pow(k1k2,2) +
			        (-2*eg1*l1k2*l2k2 + (-2*l1k2 + l2k2)*m2)*pow(k1p2,2) -
			        2*l2k2*M2*pow(l1k1,2) - 2*l2p2*p1k2*pow(l1k1,2) +
			        2*p1k2*p1l2*pow(l1k1,2) + 2*l2k2*p1p2*pow(l1k1,2) +
			        2*l2p2*p2k2*pow(l1k1,2) - 2*p1l2*p2k2*pow(l1k1,2) +
			        2*eg1*l1p2*l2p2*pow(l1k2,2) + 4*eg1*l1l2*M2*pow(l1k2,2) -
			        2*eg1*m2*M2*pow(l1k2,2) + 2*eg1*l2p2*p1k1*pow(l1k2,2) -
			        2*eg1*l2p2*p1l1*pow(l1k2,2) - 2*eg1*l1p2*p1l2*pow(l1k2,2) -
			        2*eg1*p1k1*p1l2*pow(l1k2,2) + 2*eg1*p1l1*p1l2*pow(l1k2,2) -
			        4*eg1*l1l2*p1p2*pow(l1k2,2) + 2*eg1*m2*p1p2*pow(l1k2,2) +
			        k1p2*(l1k2*(-2*l1k1*l2p2 + 2*l2k1*l2p2 + l1p2*m2 + 4*m2*p1k1 -
			              m2*p1l1 + 2*l1k1*p1l2 - 2*l2k1*p1l2 +
			              2*l1l2*(l1p2 - l2p2 - p1l1 + p1l2) +
			              2*eg1*(2*l1p2*l2k2 + 2*l2k2*(p1k1 - p1l1) -
			                 (l2k1 - m2)*(p1k2 - p2k2))) +
			           2*l1k1*l2k2*(l1p2 - p1l1 + eg1*(-p1k2 + p2k2)) +
			           m2*(-2*l1p2*l2k2 - (l1l2 - l2k1 - m2)*(p1k2 - p2k2) +
			              l2k2*(l2p2 - 2*p1k1 + eg1*p1k2 + 2*p1l1 - p1l2 - eg1*p2k2)) \
			+ 2*eg1*(-l2p2 + p1l2)*pow(l1k2,2)) + 4*l1k2*M2*pow(l1l2,2) -
			        4*l1k2*p1p2*pow(l1l2,2) - 2*l1k2*l2k1*pow(l1p2,2) -
			        2*eg1*l1k2*l2k2*pow(l1p2,2) + 2*l1k2*M2*pow(l2k1,2) -
			        2*l1k2*p1p2*pow(l2k1,2) - 2*eg1*l1k2*l2k2*pow(p1k1,2) -
			        2*l1k2*m2*pow(p1k1,2) + l2k2*m2*pow(p1k1,2) -
			        2*eg1*l1k1*l2k1*pow(p1k2,2) - 2*eg1*l1k1*m2*pow(p1k2,2) -
			        eg1*l1l2*m2*pow(p1k2,2) + eg1*l2k1*m2*pow(p1k2,2) +
			        eg1*m4*pow(p1k2,2) - 2*l1k2*m2*pow(eg1,2)*pow(p1k2,2) -
			        2*l1k2*l2k1*pow(p1l1,2) - 2*eg1*l1k2*l2k2*pow(p1l1,2) -
			        2*eg1*l1k1*l2k1*pow(p2k2,2) - 2*eg1*l1k1*m2*pow(p2k2,2) -
			        eg1*l1l2*m2*pow(p2k2,2) + eg1*l2k1*m2*pow(p2k2,2) +
			        eg1*m4*pow(p2k2,2) - 2*l1k2*m2*pow(eg1,2)*pow(p2k2,2) +
			        k1k2*(2*l1k1*l1l2*M2 - 8*eg1*l1k2*l1l2*M2 + 2*l1l2*l2k1*M2 -
			           4*eg1*l1k1*l2k2*M2 + 6*eg1*l1k2*m2*M2 + 8*l1l2*m2*M2 -
			           4*l2k1*m2*M2 + 6*m4*M2 - 2*eg1*l1k2*l2p2*p1k1 -
			           l2p2*m2*p1k1 - eg1*l2p2*m2*p1k2 - 2*eg1*l1l2*p1k1*p1k2 +
			           4*eg1*m2*p1k1*p1k2 + 2*l1k1*l2p2*p1l1 + 2*eg1*l1k2*l2p2*p1l1 +
			           2*l1l2*l2p2*p1l1 - 2*l2k1*l2p2*p1l1 + 2*l2p2*m2*p1l1 -
			           2*l1l2*p1k1*p1l1 + 2*eg1*l2k1*p1k2*p1l1 + 2*eg1*m2*p1k2*p1l1 +
			           2*l2p2*m2*p1l2 + 2*eg1*l1k2*p1k1*p1l2 + m2*p1k1*p1l2 +
			           eg1*m2*p1k2*p1l2 - 2*l1k1*p1l1*p1l2 - 2*eg1*l1k2*p1l1*p1l2 -
			           2*l1l2*p1l1*p1l2 + 2*l2k1*p1l1*p1l2 - 2*m2*p1l1*p1l2 -
			           2*l1k1*l1l2*p1p2 + 8*eg1*l1k2*l1l2*p1p2 - 2*l1l2*l2k1*p1p2 +
			           4*eg1*l1k1*l2k2*p1p2 - 6*eg1*l1k2*m2*p1p2 - 8*l1l2*m2*p1p2 +
			           4*l2k1*m2*p1p2 - 6*m4*p1p2 -
			           2*l1p2*(l1l2*l2p2 - l2k1*l2p2 + l2p2*m2 - l1l2*p1k1 +
			              2*l2k1*p1l1 - m2*p1l1 + l1k1*(l2p2 - p1l2) - l1l2*p1l2 +
			              l2k1*p1l2 - m2*p1l2 +
			              eg1*(l1k2*(l2p2 - p1l2) + (l2k1 + m2)*(p1k2 - p2k2))) +
			           eg1*l2p2*m2*p2k2 + 2*eg1*l1l2*p1k1*p2k2 - 4*eg1*m2*p1k1*p2k2 -
			           2*eg1*l2k1*p1l1*p2k2 - 2*eg1*m2*p1l1*p2k2 - eg1*m2*p1l2*p2k2 +
			           k1p2*(m2*(l2p2 - 4*p1k1 - p1l2) -
			              2*l1l2*(l1p2 - eg1*p1k2 - p1l1 + eg1*p2k2) +
			              2*eg1*(l1k2*(l2p2 - p1l2) + 2*m2*(-p1k2 + p2k2))) -
			           4*m2*p1k2*p2k2*pow(eg1,2) + 2*m2*pow(k1p2,2) -
			           2*M2*pow(l1l2,2) + 2*p1p2*pow(l1l2,2) +
			           (2*l2k1 - m2)*pow(l1p2,2) - m2*pow(l2p2,2) +
			           2*m2*pow(p1k1,2) + 2*m2*pow(eg1,2)*pow(p1k2,2) +
			           2*l2k1*pow(p1l1,2) - m2*pow(p1l1,2) - m2*pow(p1l2,2) +
			           2*m2*pow(eg1,2)*pow(p2k2,2))) +
			     pow(f2,2)*(-6*l1k1*l1k2*l1p2*l2p2*M2 - 12*l1k2*l1l2*l1p2*l2p2*M2 +
			        12*l1k2*l1p2*l2k1*l2p2*M2 + 6*l1p2*l2k2*l2p2*m2*M2 +
			        2*l1k2*l1l2*l1p2*M2*p1k1 + 2*l1k1*l1p2*l2k2*M2*p1k1 +
			        4*eg1*l1k2*l1p2*l2k2*M2*p1k1 - 2*l1k1*l1k2*l2p2*M2*p1k1 -
			        2*l1k2*l1l2*l2p2*M2*p1k1 + 2*l1k2*l2k1*l2p2*M2*p1k1 +
			        l1k2*l1p2*m2*M2*p1k1 - 2*l1p2*l2k2*m2*M2*p1k1 +
			        l2k2*l2p2*m2*M2*p1k1 + 2*eg1*l1k2*l1l2*l1p2*M2*p1k2 -
			        2*l1k1*l1p2*l2k1*M2*p1k2 - 4*eg1*l1k2*l1p2*l2k1*M2*p1k2 +
			        2*eg1*l1k1*l1k2*l2p2*M2*p1k2 + 2*l1k1*l1l2*l2p2*M2*p1k2 -
			        2*l1k1*l2k1*l2p2*M2*p1k2 - l1k1*l1p2*m2*M2*p1k2 -
			        2*eg1*l1k2*l1p2*m2*M2*p1k2 + l1p2*l2k1*m2*M2*p1k2 +
			        eg1*l1p2*l2k2*m2*M2*p1k2 - 3*l1k1*l2p2*m2*M2*p1k2 -
			        eg1*l1k2*l2p2*m2*M2*p1k2 + l2k1*l2p2*m2*M2*p1k2 -
			        6*eg1*l1k2*l2k1*M2*p1k1*p1k2 - 6*eg1*l1k1*l2k2*M2*p1k1*p1k2 +
			        6*eg1*l1k2*m2*M2*p1k1*p1k2 - 3*l1l2*m2*M2*p1k1*p1k2 +
			        3*l2k1*m2*M2*p1k1*p1k2 + 3*eg1*l2k2*m2*M2*p1k1*p1k2 +
			        3*m4*M2*p1k1*p1k2 - 4*l1k2*l1p2*l2k1*M2*p1l1 -
			        4*eg1*l1k2*l1p2*l2k2*M2*p1l1 + 2*l1k1*l1k2*l2p2*M2*p1l1 +
			        4*l1k2*l1l2*l2p2*M2*p1l1 - 4*l1k2*l2k1*l2p2*M2*p1l1 -
			        2*l2k2*l2p2*m2*M2*p1l1 - 6*l1k2*l1l2*M2*p1k1*p1l1 -
			        6*l1k1*l2k2*M2*p1k1*p1l1 - 12*eg1*l1k2*l2k2*M2*p1k1*p1l1 -
			        3*l1k2*m2*M2*p1k1*p1l1 + 6*l2k2*m2*M2*p1k1*p1l1 -
			        6*eg1*l1k2*l1l2*M2*p1k2*p1l1 + 6*l1k1*l2k1*M2*p1k2*p1l1 +
			        12*eg1*l1k2*l2k1*M2*p1k2*p1l1 + 3*l1k1*m2*M2*p1k2*p1l1 +
			        6*eg1*l1k2*m2*M2*p1k2*p1l1 - 3*l2k1*m2*M2*p1k2*p1l1 -
			        3*eg1*l2k2*m2*M2*p1k2*p1l1 + 2*l1k1*l1k2*l1p2*M2*p1l2 +
			        4*l1k2*l1l2*l1p2*M2*p1l2 - 4*l1k2*l1p2*l2k1*M2*p1l2 -
			        2*l1p2*l2k2*m2*M2*p1l2 + 6*l1k1*l1k2*M2*p1k1*p1l2 +
			        6*l1k2*l1l2*M2*p1k1*p1l2 - 6*l1k2*l2k1*M2*p1k1*p1l2 -
			        3*l2k2*m2*M2*p1k1*p1l2 - 6*eg1*l1k1*l1k2*M2*p1k2*p1l2 -
			        6*l1k1*l1l2*M2*p1k2*p1l2 + 6*l1k1*l2k1*M2*p1k2*p1l2 +
			        9*l1k1*m2*M2*p1k2*p1l2 + 3*eg1*l1k2*m2*M2*p1k2*p1l2 -
			        3*l2k1*m2*M2*p1k2*p1l2 - 6*l1k1*l1k2*M2*p1l1*p1l2 -
			        12*l1k2*l1l2*M2*p1l1*p1l2 + 12*l1k2*l2k1*M2*p1l1*p1l2 +
			        6*l2k2*m2*M2*p1l1*p1l2 + 2*l1k1*l1k2*l1p2*l2p2*p1p2 +
			        4*l1k2*l1l2*l1p2*l2p2*p1p2 - 4*l1k2*l1p2*l2k1*l2p2*p1p2 -
			        2*l1p2*l2k2*l2p2*m2*p1p2 + 8*l1k1*l1k2*l1l2*M2*p1p2 +
			        4*l1k1*l1k2*l2k1*M2*p1p2 - 12*l1k2*l1l2*l2k1*M2*p1p2 +
			        4*l1k1*l1l2*l2k2*M2*p1p2 - 4*l1k1*l2k1*l2k2*M2*p1p2 -
			        24*l1k2*l1l2*m2*M2*p1p2 + 16*l1k2*l2k1*m2*M2*p1p2 -
			        4*l1k1*l2k2*m2*M2*p1p2 - 4*eg1*l1k2*l2k2*m2*M2*p1p2 -
			        4*l1l2*l2k2*m2*M2*p1p2 + 4*l2k1*l2k2*m2*M2*p1p2 +
			        12*l2k2*m4*M2*p1p2 + 2*l1k2*l1l2*l1p2*p1k1*p1p2 +
			        2*l1k1*l1p2*l2k2*p1k1*p1p2 + 4*eg1*l1k2*l1p2*l2k2*p1k1*p1p2 -
			        2*l1k1*l1k2*l2p2*p1k1*p1p2 - 2*l1k2*l1l2*l2p2*p1k1*p1p2 +
			        2*l1k2*l2k1*l2p2*p1k1*p1p2 + l1k2*l1p2*m2*p1k1*p1p2 -
			        2*l1p2*l2k2*m2*p1k1*p1p2 + l2k2*l2p2*m2*p1k1*p1p2 +
			        2*eg1*l1k2*l1l2*l1p2*p1k2*p1p2 - 2*l1k1*l1p2*l2k1*p1k2*p1p2 -
			        4*eg1*l1k2*l1p2*l2k1*p1k2*p1p2 + 2*eg1*l1k1*l1k2*l2p2*p1k2*p1p2 +
			        2*l1k1*l1l2*l2p2*p1k2*p1p2 - 2*l1k1*l2k1*l2p2*p1k2*p1p2 -
			        l1k1*l1p2*m2*p1k2*p1p2 - 2*eg1*l1k2*l1p2*m2*p1k2*p1p2 +
			        l1p2*l2k1*m2*p1k2*p1p2 + eg1*l1p2*l2k2*m2*p1k2*p1p2 -
			        3*l1k1*l2p2*m2*p1k2*p1p2 - eg1*l1k2*l2p2*m2*p1k2*p1p2 +
			        l2k1*l2p2*m2*p1k2*p1p2 + 2*eg1*l1k2*l2k1*p1k1*p1k2*p1p2 +
			        2*eg1*l1k1*l2k2*p1k1*p1k2*p1p2 - 2*eg1*l1k2*m2*p1k1*p1k2*p1p2 +
			        l1l2*m2*p1k1*p1k2*p1p2 - l2k1*m2*p1k1*p1k2*p1p2 -
			        eg1*l2k2*m2*p1k1*p1k2*p1p2 - m4*p1k1*p1k2*p1p2 -
			        4*l1k2*l1p2*l2k1*p1l1*p1p2 - 4*eg1*l1k2*l1p2*l2k2*p1l1*p1p2 +
			        2*l1k1*l1k2*l2p2*p1l1*p1p2 + 4*l1k2*l1l2*l2p2*p1l1*p1p2 -
			        4*l1k2*l2k1*l2p2*p1l1*p1p2 - 2*l2k2*l2p2*m2*p1l1*p1p2 +
			        2*l1k2*l1l2*p1k1*p1l1*p1p2 + 2*l1k1*l2k2*p1k1*p1l1*p1p2 +
			        4*eg1*l1k2*l2k2*p1k1*p1l1*p1p2 + l1k2*m2*p1k1*p1l1*p1p2 -
			        2*l2k2*m2*p1k1*p1l1*p1p2 + 2*eg1*l1k2*l1l2*p1k2*p1l1*p1p2 -
			        2*l1k1*l2k1*p1k2*p1l1*p1p2 - 4*eg1*l1k2*l2k1*p1k2*p1l1*p1p2 -
			        l1k1*m2*p1k2*p1l1*p1p2 - 2*eg1*l1k2*m2*p1k2*p1l1*p1p2 +
			        l2k1*m2*p1k2*p1l1*p1p2 + eg1*l2k2*m2*p1k2*p1l1*p1p2 +
			        2*l1k1*l1k2*l1p2*p1l2*p1p2 + 4*l1k2*l1l2*l1p2*p1l2*p1p2 -
			        4*l1k2*l1p2*l2k1*p1l2*p1p2 - 2*l1p2*l2k2*m2*p1l2*p1p2 -
			        2*l1k1*l1k2*p1k1*p1l2*p1p2 - 2*l1k2*l1l2*p1k1*p1l2*p1p2 +
			        2*l1k2*l2k1*p1k1*p1l2*p1p2 + l2k2*m2*p1k1*p1l2*p1p2 +
			        2*eg1*l1k1*l1k2*p1k2*p1l2*p1p2 + 2*l1k1*l1l2*p1k2*p1l2*p1p2 -
			        2*l1k1*l2k1*p1k2*p1l2*p1p2 - 3*l1k1*m2*p1k2*p1l2*p1p2 -
			        eg1*l1k2*m2*p1k2*p1l2*p1p2 + l2k1*m2*p1k2*p1l2*p1p2 +
			        2*l1k1*l1k2*p1l1*p1l2*p1p2 + 4*l1k2*l1l2*p1l1*p1l2*p1p2 -
			        4*l1k2*l2k1*p1l1*p1l2*p1p2 - 2*l2k2*m2*p1l1*p1l2*p1p2 -
			        6*eg1*l1k2*l1l2*l1p2*M2*p2k2 + 6*l1k1*l1p2*l2k1*M2*p2k2 +
			        12*eg1*l1k2*l1p2*l2k1*M2*p2k2 - 6*eg1*l1k1*l1k2*l2p2*M2*p2k2 -
			        6*l1k1*l1l2*l2p2*M2*p2k2 + 6*l1k1*l2k1*l2p2*M2*p2k2 +
			        3*l1k1*l1p2*m2*M2*p2k2 + 6*eg1*l1k2*l1p2*m2*M2*p2k2 -
			        3*l1p2*l2k1*m2*M2*p2k2 - 3*eg1*l1p2*l2k2*m2*M2*p2k2 +
			        9*l1k1*l2p2*m2*M2*p2k2 + 3*eg1*l1k2*l2p2*m2*M2*p2k2 -
			        3*l2k1*l2p2*m2*M2*p2k2 + 2*eg1*l1k2*l2k1*M2*p1k1*p2k2 +
			        2*eg1*l1k1*l2k2*M2*p1k1*p2k2 - 2*eg1*l1k2*m2*M2*p1k1*p2k2 +
			        l1l2*m2*M2*p1k1*p2k2 - l2k1*m2*M2*p1k1*p2k2 -
			        eg1*l2k2*m2*M2*p1k1*p2k2 - m4*M2*p1k1*p2k2 -
			        4*eg1*l1k1*l2k1*M2*p1k2*p2k2 - 4*eg1*l1k1*m2*M2*p1k2*p2k2 -
			        2*eg1*l1l2*m2*M2*p1k2*p2k2 + 2*eg1*l2k1*m2*M2*p1k2*p2k2 +
			        2*eg1*m4*M2*p1k2*p2k2 + 2*eg1*l1k2*l1l2*M2*p1l1*p2k2 -
			        2*l1k1*l2k1*M2*p1l1*p2k2 - 4*eg1*l1k2*l2k1*M2*p1l1*p2k2 -
			        l1k1*m2*M2*p1l1*p2k2 - 2*eg1*l1k2*m2*M2*p1l1*p2k2 +
			        l2k1*m2*M2*p1l1*p2k2 + eg1*l2k2*m2*M2*p1l1*p2k2 +
			        2*eg1*l1k1*l1k2*M2*p1l2*p2k2 + 2*l1k1*l1l2*M2*p1l2*p2k2 -
			        2*l1k1*l2k1*M2*p1l2*p2k2 - 3*l1k1*m2*M2*p1l2*p2k2 -
			        eg1*l1k2*m2*M2*p1l2*p2k2 + l2k1*m2*M2*p1l2*p2k2 +
			        2*eg1*l1k2*l1l2*l1p2*p1p2*p2k2 - 2*l1k1*l1p2*l2k1*p1p2*p2k2 -
			        4*eg1*l1k2*l1p2*l2k1*p1p2*p2k2 + 2*eg1*l1k1*l1k2*l2p2*p1p2*p2k2 +
			        2*l1k1*l1l2*l2p2*p1p2*p2k2 - 2*l1k1*l2k1*l2p2*p1p2*p2k2 -
			        l1k1*l1p2*m2*p1p2*p2k2 - 2*eg1*l1k2*l1p2*m2*p1p2*p2k2 +
			        l1p2*l2k1*m2*p1p2*p2k2 + eg1*l1p2*l2k2*m2*p1p2*p2k2 -
			        3*l1k1*l2p2*m2*p1p2*p2k2 - eg1*l1k2*l2p2*m2*p1p2*p2k2 +
			        l2k1*l2p2*m2*p1p2*p2k2 + 2*eg1*l1k2*l2k1*p1k1*p1p2*p2k2 +
			        2*eg1*l1k1*l2k2*p1k1*p1p2*p2k2 - 2*eg1*l1k2*m2*p1k1*p1p2*p2k2 +
			        l1l2*m2*p1k1*p1p2*p2k2 - l2k1*m2*p1k1*p1p2*p2k2 -
			        eg1*l2k2*m2*p1k1*p1p2*p2k2 - m4*p1k1*p1p2*p2k2 -
			        4*eg1*l1k1*l2k1*p1k2*p1p2*p2k2 - 4*eg1*l1k1*m2*p1k2*p1p2*p2k2 -
			        2*eg1*l1l2*m2*p1k2*p1p2*p2k2 + 2*eg1*l2k1*m2*p1k2*p1p2*p2k2 +
			        2*eg1*m4*p1k2*p1p2*p2k2 + 2*eg1*l1k2*l1l2*p1l1*p1p2*p2k2 -
			        2*l1k1*l2k1*p1l1*p1p2*p2k2 - 4*eg1*l1k2*l2k1*p1l1*p1p2*p2k2 -
			        l1k1*m2*p1l1*p1p2*p2k2 - 2*eg1*l1k2*m2*p1l1*p1p2*p2k2 +
			        l2k1*m2*p1l1*p1p2*p2k2 + eg1*l2k2*m2*p1l1*p1p2*p2k2 +
			        2*eg1*l1k1*l1k2*p1l2*p1p2*p2k2 + 2*l1k1*l1l2*p1l2*p1p2*p2k2 -
			        2*l1k1*l2k1*p1l2*p1p2*p2k2 - 3*l1k1*m2*p1l2*p1p2*p2k2 -
			        eg1*l1k2*m2*p1l2*p1p2*p2k2 + l2k1*m2*p1l2*p1p2*p2k2 -
			        4*l1k2*m2*M2*p1k2*p2k2*pow(eg1,2) -
			        4*l1k2*m2*p1k2*p1p2*p2k2*pow(eg1,2) -
			        8*eg1*(l1l2 - 2*m2)*M2*(M2 - p1p2)*pow(k1k2,2) +
			        (2*eg1*l1k2*l2k2 + (2*l1k2 - l2k2)*m2)*(3*M2 - p1p2)*
			         pow(k1p2,2) + 2*l2p2*M2*p1k2*pow(l1k1,2) -
			        6*M2*p1k2*p1l2*pow(l1k1,2) - 4*l2k2*M2*p1p2*pow(l1k1,2) +
			        2*l2p2*p1k2*p1p2*pow(l1k1,2) + 2*p1k2*p1l2*p1p2*pow(l1k1,2) -
			        6*l2p2*M2*p2k2*pow(l1k1,2) + 2*M2*p1l2*p2k2*pow(l1k1,2) +
			        2*l2p2*p1p2*p2k2*pow(l1k1,2) + 2*p1l2*p1p2*p2k2*pow(l1k1,2) -
			        6*eg1*l1p2*l2p2*M2*pow(l1k2,2) -
			        2*eg1*l2p2*M2*p1k1*pow(l1k2,2) +
			        2*eg1*l2p2*M2*p1l1*pow(l1k2,2) +
			        2*eg1*l1p2*M2*p1l2*pow(l1k2,2) +
			        6*eg1*M2*p1k1*p1l2*pow(l1k2,2) -
			        6*eg1*M2*p1l1*p1l2*pow(l1k2,2) +
			        2*eg1*l1p2*l2p2*p1p2*pow(l1k2,2) +
			        8*eg1*l1l2*M2*p1p2*pow(l1k2,2) - 4*eg1*m2*M2*p1p2*pow(l1k2,2) -
			        2*eg1*l2p2*p1k1*p1p2*pow(l1k2,2) +
			        2*eg1*l2p2*p1l1*p1p2*pow(l1k2,2) +
			        2*eg1*l1p2*p1l2*p1p2*pow(l1k2,2) -
			        2*eg1*p1k1*p1l2*p1p2*pow(l1k2,2) +
			        2*eg1*p1l1*p1l2*p1p2*pow(l1k2,2) +
			        k1p2*(2*l1k1*l2k2*(l1p2*(-3*M2 + p1p2) + p1l1*(M2 + p1p2) +
			              eg1*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2))) -
			           m2*(2*l1p2*l2k2*(-3*M2 + p1p2) -
			              (l1l2 - l2k1 - m2)*
			               (M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)) +
			              l2k2*(l2p2*(3*M2 - p1p2) -
			                 M2*(2*p1k1 - eg1*p1k2 - 2*p1l1 + p1l2 + 3*eg1*p2k2) +
			                 p1p2*(-2*p1k1 + 2*p1l1 - p1l2 + eg1*(p1k2 + p2k2)))) +
			           l1k2*(6*l1k1*l2p2*M2 - 6*l2k1*l2p2*M2 - 3*l1p2*m2*M2 -
			              4*m2*M2*p1k1 + m2*M2*p1l1 - 2*l1k1*M2*p1l2 +
			              2*l2k1*M2*p1l2 - 2*l1k1*l2p2*p1p2 + 2*l2k1*l2p2*p1p2 +
			              l1p2*m2*p1p2 - 4*m2*p1k1*p1p2 + m2*p1l1*p1p2 -
			              2*l1k1*p1l2*p1p2 + 2*l2k1*p1l2*p1p2 -
			              2*l1l2*(l1p2*(3*M2 - p1p2) + l2p2*(-3*M2 + p1p2) -
			                 (p1l1 - p1l2)*(M2 + p1p2)) -
			              2*eg1*(l1p2*l2k2*(6*M2 - 2*p1p2) +
			                 2*l2k2*(p1k1 - p1l1)*(M2 + p1p2) -
			                 (l2k1 - m2)*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)))) \
			+ 2*eg1*(l2p2*(3*M2 - p1p2) - p1l2*(M2 + p1p2))*pow(l1k2,2)) +
			        8*l1k2*M2*p1p2*pow(l1l2,2) + 6*l1k2*l2k1*M2*pow(l1p2,2) +
			        6*eg1*l1k2*l2k2*M2*pow(l1p2,2) - 2*l1k2*l2k1*p1p2*pow(l1p2,2) -
			        2*eg1*l1k2*l2k2*p1p2*pow(l1p2,2) + 4*l1k2*M2*p1p2*pow(l2k1,2) -
			        4*l1k1*l1k2*l1l2*M4 - 6*l1k1*l1k2*l2k1*M4 +
			        6*l1k2*l1l2*l2k1*M4 + 4*eg1*l1k1*l1k2*l2k2*M4 -
			        2*l1k1*l1l2*l2k2*M4 + 2*l1k1*l2k1*l2k2*M4 +
			        20*l1k2*l1l2*m2*M4 - 14*l1k2*l2k1*m2*M4 +
			        2*eg1*l1k2*l2k2*m2*M4 + 2*l1l2*l2k2*m2*M4 -
			        2*l2k1*l2k2*m2*M4 - 10*l2k2*m4*M4 +
			        6*l2k2*pow(l1k1,2)*M4 - 4*eg1*l1l2*pow(l1k2,2)*M4 -
			        4*eg1*l2k1*pow(l1k2,2)*M4 +
			        2*eg1*m2*pow(l1k2,2)*M4 - 4*l1k2*pow(l1l2,2)*M4 -
			        2*l1k2*pow(l2k1,2)*M4 + 6*eg1*l1k2*l2k2*M2*pow(p1k1,2) +
			        6*l1k2*m2*M2*pow(p1k1,2) - 3*l2k2*m2*M2*pow(p1k1,2) -
			        2*eg1*l1k2*l2k2*p1p2*pow(p1k1,2) - 2*l1k2*m2*p1p2*pow(p1k1,2) +
			        l2k2*m2*p1p2*pow(p1k1,2) + 6*eg1*l1k1*l2k1*M2*pow(p1k2,2) +
			        6*eg1*l1k1*m2*M2*pow(p1k2,2) + 3*eg1*l1l2*m2*M2*pow(p1k2,2) -
			        3*eg1*l2k1*m2*M2*pow(p1k2,2) - 3*eg1*m4*M2*pow(p1k2,2) -
			        2*eg1*l1k1*l2k1*p1p2*pow(p1k2,2) - 2*eg1*l1k1*m2*p1p2*pow(p1k2,2) -
			        eg1*l1l2*m2*p1p2*pow(p1k2,2) + eg1*l2k1*m2*p1p2*pow(p1k2,2) +
			        eg1*m4*p1p2*pow(p1k2,2) + 6*l1k2*m2*M2*pow(eg1,2)*pow(p1k2,2) -
			        2*l1k2*m2*p1p2*pow(eg1,2)*pow(p1k2,2) +
			        6*l1k2*l2k1*M2*pow(p1l1,2) + 6*eg1*l1k2*l2k2*M2*pow(p1l1,2) -
			        2*l1k2*l2k1*p1p2*pow(p1l1,2) - 2*eg1*l1k2*l2k2*p1p2*pow(p1l1,2) -
			        4*l1k1*l1k2*l1l2*pow(p1p2,2) + 2*l1k1*l1k2*l2k1*pow(p1p2,2) +
			        6*l1k2*l1l2*l2k1*pow(p1p2,2) - 4*eg1*l1k1*l1k2*l2k2*pow(p1p2,2) -
			        2*l1k1*l1l2*l2k2*pow(p1p2,2) + 2*l1k1*l2k1*l2k2*pow(p1p2,2) +
			        4*l1k2*l1l2*m2*pow(p1p2,2) - 2*l1k2*l2k1*m2*pow(p1p2,2) +
			        4*l1k1*l2k2*m2*pow(p1p2,2) + 2*eg1*l1k2*l2k2*m2*pow(p1p2,2) +
			        2*l1l2*l2k2*m2*pow(p1p2,2) - 2*l2k1*l2k2*m2*pow(p1p2,2) -
			        2*l2k2*m4*pow(p1p2,2) - 2*l2k2*pow(l1k1,2)*pow(p1p2,2) -
			        4*eg1*l1l2*pow(l1k2,2)*pow(p1p2,2) +
			        4*eg1*l2k1*pow(l1k2,2)*pow(p1p2,2) +
			        2*eg1*m2*pow(l1k2,2)*pow(p1p2,2) - 4*l1k2*pow(l1l2,2)*pow(p1p2,2) -
			        2*l1k2*pow(l2k1,2)*pow(p1p2,2) + 6*eg1*l1k1*l2k1*M2*pow(p2k2,2) +
			        6*eg1*l1k1*m2*M2*pow(p2k2,2) + 3*eg1*l1l2*m2*M2*pow(p2k2,2) -
			        3*eg1*l2k1*m2*M2*pow(p2k2,2) - 3*eg1*m4*M2*pow(p2k2,2) -
			        2*eg1*l1k1*l2k1*p1p2*pow(p2k2,2) - 2*eg1*l1k1*m2*p1p2*pow(p2k2,2) -
			        eg1*l1l2*m2*p1p2*pow(p2k2,2) + eg1*l2k1*m2*p1p2*pow(p2k2,2) +
			        eg1*m4*p1p2*pow(p2k2,2) + 6*l1k2*m2*M2*pow(eg1,2)*pow(p2k2,2) -
			        2*l1k2*m2*p1p2*pow(eg1,2)*pow(p2k2,2) +
			        k1k2*(2*eg1*l1k2*l2p2*M2*p1k1 + l2p2*m2*M2*p1k1 +
			           eg1*l2p2*m2*M2*p1k2 + 6*eg1*l1l2*M2*p1k1*p1k2 -
			           12*eg1*m2*M2*p1k1*p1k2 - 2*l1k1*l2p2*M2*p1l1 -
			           2*eg1*l1k2*l2p2*M2*p1l1 - 2*l1l2*l2p2*M2*p1l1 +
			           2*l2k1*l2p2*M2*p1l1 - 2*l2p2*m2*M2*p1l1 +
			           6*l1l2*M2*p1k1*p1l1 - 6*eg1*l2k1*M2*p1k2*p1l1 -
			           6*eg1*m2*M2*p1k2*p1l1 - 2*l2p2*m2*M2*p1l2 -
			           6*eg1*l1k2*M2*p1k1*p1l2 - 3*m2*M2*p1k1*p1l2 -
			           3*eg1*m2*M2*p1k2*p1l2 + 6*l1k1*M2*p1l1*p1l2 +
			           6*eg1*l1k2*M2*p1l1*p1l2 + 6*l1l2*M2*p1l1*p1l2 -
			           6*l2k1*M2*p1l1*p1l2 + 6*m2*M2*p1l1*p1l2 +
			           4*l1k1*l1l2*M2*p1p2 - 16*eg1*l1k2*l1l2*M2*p1p2 +
			           4*l1l2*l2k1*M2*p1p2 - 8*eg1*l1k1*l2k2*M2*p1p2 +
			           12*eg1*l1k2*m2*M2*p1p2 + 16*l1l2*m2*M2*p1p2 -
			           8*l2k1*m2*M2*p1p2 + 12*m4*M2*p1p2 +
			           2*eg1*l1k2*l2p2*p1k1*p1p2 + l2p2*m2*p1k1*p1p2 +
			           eg1*l2p2*m2*p1k2*p1p2 - 2*eg1*l1l2*p1k1*p1k2*p1p2 +
			           4*eg1*m2*p1k1*p1k2*p1p2 - 2*l1k1*l2p2*p1l1*p1p2 -
			           2*eg1*l1k2*l2p2*p1l1*p1p2 - 2*l1l2*l2p2*p1l1*p1p2 +
			           2*l2k1*l2p2*p1l1*p1p2 - 2*l2p2*m2*p1l1*p1p2 -
			           2*l1l2*p1k1*p1l1*p1p2 + 2*eg1*l2k1*p1k2*p1l1*p1p2 +
			           2*eg1*m2*p1k2*p1l1*p1p2 - 2*l2p2*m2*p1l2*p1p2 +
			           2*eg1*l1k2*p1k1*p1l2*p1p2 + m2*p1k1*p1l2*p1p2 +
			           eg1*m2*p1k2*p1l2*p1p2 - 2*l1k1*p1l1*p1l2*p1p2 -
			           2*eg1*l1k2*p1l1*p1l2*p1p2 - 2*l1l2*p1l1*p1l2*p1p2 +
			           2*l2k1*p1l1*p1l2*p1p2 - 2*m2*p1l1*p1l2*p1p2 -
			           3*eg1*l2p2*m2*M2*p2k2 - 2*eg1*l1l2*M2*p1k1*p2k2 +
			           4*eg1*m2*M2*p1k1*p2k2 + 2*eg1*l2k1*M2*p1l1*p2k2 +
			           2*eg1*m2*M2*p1l1*p2k2 + eg1*m2*M2*p1l2*p2k2 +
			           eg1*l2p2*m2*p1p2*p2k2 - 2*eg1*l1l2*p1k1*p1p2*p2k2 +
			           4*eg1*m2*p1k1*p1p2*p2k2 + 2*eg1*l2k1*p1l1*p1p2*p2k2 +
			           2*eg1*m2*p1l1*p1p2*p2k2 + eg1*m2*p1l2*p1p2*p2k2 +
			           k1p2*(m2*(l2p2*(-3*M2 + p1p2) +
			                 (4*p1k1 + p1l2)*(M2 + p1p2)) -
			              2*l1l2*(l1p2*(-3*M2 + p1p2) + p1l1*(M2 + p1p2) +
			                 eg1*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2))) +
			              2*eg1*(l1k2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) +
			                 2*m2*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)))) +
			           2*l1p2*(3*l1l2*l2p2*M2 - 3*l2k1*l2p2*M2 + 3*l2p2*m2*M2 -
			              l1l2*M2*p1k1 + 2*l2k1*M2*p1l1 - m2*M2*p1l1 -
			              l1l2*M2*p1l2 + l2k1*M2*p1l2 - m2*M2*p1l2 -
			              l1l2*l2p2*p1p2 + l2k1*l2p2*p1p2 - l2p2*m2*p1p2 -
			              l1l2*p1k1*p1p2 + 2*l2k1*p1l1*p1p2 - m2*p1l1*p1p2 -
			              l1l2*p1l2*p1p2 + l2k1*p1l2*p1p2 - m2*p1l2*p1p2 +
			              l1k1*(l2p2*(3*M2 - p1p2) - p1l2*(M2 + p1p2)) +
			              eg1*(l1k2*(l2p2*(3*M2 - p1p2) - p1l2*(M2 + p1p2)) +
			                 (l2k1 + m2)*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)))) \
			+ 4*m2*M2*p1k2*p2k2*pow(eg1,2) + 4*m2*p1k2*p1p2*p2k2*pow(eg1,2) +
			           2*m2*(-3*M2 + p1p2)*pow(k1p2,2) - 4*M2*p1p2*pow(l1l2,2) -
			           (2*l2k1 - m2)*(3*M2 - p1p2)*pow(l1p2,2) +
			           3*m2*M2*pow(l2p2,2) - m2*p1p2*pow(l2p2,2) -
			           6*l1k1*l1l2*M4 + 12*eg1*l1k2*l1l2*M4 +
			           4*eg1*l1k2*l2k1*M4 - 2*l1l2*l2k1*M4 +
			           8*eg1*l1k1*l2k2*M4 - 10*eg1*l1k2*m2*M4 -
			           14*l1l2*m2*M4 + 8*l2k1*m2*M4 -
			           12*m4*M4 + 2*pow(l1l2,2)*M4 -
			           6*m2*M2*pow(p1k1,2) + 2*m2*p1p2*pow(p1k1,2) -
			           6*m2*M2*pow(eg1,2)*pow(p1k2,2) +
			           2*m2*p1p2*pow(eg1,2)*pow(p1k2,2) - 6*l2k1*M2*pow(p1l1,2) +
			           3*m2*M2*pow(p1l1,2) + 2*l2k1*p1p2*pow(p1l1,2) -
			           m2*p1p2*pow(p1l1,2) + 3*m2*M2*pow(p1l2,2) -
			           m2*p1p2*pow(p1l2,2) + 2*l1k1*l1l2*pow(p1p2,2) +
			           4*eg1*l1k2*l1l2*pow(p1p2,2) - 4*eg1*l1k2*l2k1*pow(p1p2,2) -
			           2*l1l2*l2k1*pow(p1p2,2) - 2*eg1*l1k2*m2*pow(p1p2,2) -
			           2*l1l2*m2*pow(p1p2,2) + 2*pow(l1l2,2)*pow(p1p2,2) -
			           6*m2*M2*pow(eg1,2)*pow(p2k2,2) +
			           2*m2*p1p2*pow(eg1,2)*pow(p2k2,2)))) +
			  4*eg1*pow(l1k2,-1)*pow(l1k1 + eg1*(-k1k2 + l1k2),-1)*pow(l2k1,-1)*
			   pow(l2k1 + eg1*(k1k2 + l2k2),-1)/M2*
			   (4*M2*pow(f1,2)*(-2*l1k2*l1l2*l2k1*M2 - 4*l1l2*l2k1*l2k2*M2 -
			        2*l1k2*l1l2*m2*M2 + 3*l1k2*l2k1*m2*M2 - 8*l1l2*l2k2*m2*M2 +
			        4*l1k2*m4*M2 - 2*l1l2*l1p2*l2k2*p1k1 + 2*l1p2*l2k1*l2k2*p1k1 -
			        2*l1k2*l2k1*l2p2*p1k1 + 2*l1l2*l2k2*l2p2*p1k1 -
			        2*k1p2*l1k2*m2*p1k1 + l1k2*l1p2*m2*p1k1 + 4*k1p2*l2k2*m2*p1k1 -
			        2*l1k2*l2p2*m2*p1k1 + l2k2*l2p2*m2*p1k1 + 2*l1l2*l1p2*l2k1*p1k2 +
			        k1p2*l1l2*m2*p1k2 - 3*l1p2*l2k1*m2*p1k2 - l2k1*l2p2*m2*p1k2 -
			        k1p2*m4*p1k2 - 2*k1p2*l1l2*l2k2*p1l1 + 2*k1p2*l2k1*l2k2*p1l1 -
			        4*l1l2*l2k2*l2p2*p1l1 + 2*l2k1*l2k2*l2p2*p1l1 + k1p2*l1k2*m2*p1l1 +
			        2*l1k2*l2p2*m2*p1l1 - 2*k1p2*l1k2*l2k1*p1l2 +
			        2*k1p2*l1l2*l2k2*p1l2 - 4*l1l2*l1p2*l2k2*p1l2 +
			        2*l1p2*l2k1*l2k2*p1l2 - 2*k1p2*l1k2*m2*p1l2 +
			        2*l1k2*l1p2*m2*p1l2 + k1p2*l2k2*m2*p1l2 + l1k2*l2k1*m2*p1p2 +
			        4*l1l2*l2k2*m2*p1p2 - 2*l1k2*m4*p1p2 + l1l2*m2*p1k1*p2k2 -
			        m4*p1k1*p2k2 + 2*l1l2*l2k1*p1l1*p2k2 - 3*l2k1*m2*p1l1*p2k2 -
			        l2k1*m2*p1l2*p2k2 + l1k1*
			         (6*l1l2*l2k2*M2 - 5*l2k2*m2*M2 - 2*l1k2*(l2k1 + m2)*M2 -
			           2*l1p2*l2k2*p1k1 + 2*l1p2*l2k1*p1k2 + 2*l2k1*l2p2*p1k2 +
			           k1p2*m2*p1k2 + l1p2*m2*p1k2 + l2p2*m2*p1k2 -
			           2*k1p2*l2k2*p1l1 - 4*l2k2*l2p2*p1l1 - 4*l1p2*l2k2*p1l2 -
			           4*l2k2*l2p2*p1l2 - 2*l2k1*l2k2*p1p2 + 3*l2k2*m2*p1p2 +
			           m2*p1k1*p2k2 + 2*l2k1*p1l1*p2k2 + m2*p1l1*p2k2 +
			           2*l2k1*p1l2*p2k2 + m2*p1l2*p2k2) + 2*l2k2*M2*pow(l1k1,2) +
			        4*l2k2*M2*pow(l1l2,2) +
			        k1k2*(-3*m4*M2 + 4*k1p2*m2*p1k1 - l1p2*m2*p1k1 +
			           2*l2k1*l2p2*p1l1 - k1p2*m2*p1l1 - 2*l1p2*m2*p1l1 -
			           2*l2p2*m2*p1l1 + 2*l1p2*l2k1*p1l2 - 2*l1p2*m2*p1l2 -
			           2*l2p2*m2*p1l2 + 3*m4*p1p2 +
			           2*l1k1*(l1l2*M2 - l2p2*p1l1 - l1p2*p1l2 - 2*l2p2*p1l2 +
			              m2*(-M2 + p1p2)) +
			           l1l2*(m2*(-5*M2 + 3*p1p2) +
			              2*(l2p2*(p1k1 - p1l1) + k1p2*p1l2 - l1p2*p1l2 - l2k1*p1p2)) +
			           2*M2*pow(l1l2,2)) - 2*l1p2*p1k2*pow(l2k1,2) +
			        2*l1k2*p1p2*pow(l2k1,2) - 2*p1l1*p2k2*pow(l2k1,2)) +
			     4*f1*f2*M2*(4*l1k1*l1p2*l2k2*l2p2 + 4*l1l2*l1p2*l2k2*l2p2 -
			        2*l1p2*l2k1*l2k2*l2p2 - 2*l1k2*l1p2*l2p2*m2 -
			        2*l1k1*l1k2*l2k1*M2 - 2*l1k2*l1l2*l2k1*M2 +
			        6*l1k1*l1l2*l2k2*M2 + 2*l1k1*l2k1*l2k2*M2 -
			        4*l1l2*l2k1*l2k2*M2 - 2*l1k1*l1k2*m2*M2 -
			        2*l1k2*l1l2*m2*M2 + 2*l1k2*l2k1*m2*M2 - 8*l1k1*l2k2*m2*M2 -
			        12*l1l2*l2k2*m2*M2 + 6*l1k2*m4*M2 - 2*l1k1*l1p2*l2k2*p1k1 -
			        2*l1l2*l1p2*l2k2*p1k1 + 2*l1p2*l2k1*l2k2*p1k1 -
			        2*l1k2*l2k1*l2p2*p1k1 + 2*l1l2*l2k2*l2p2*p1k1 + l1k2*l1p2*m2*p1k1 -
			        2*l1k2*l2p2*m2*p1k1 + l2k2*l2p2*m2*p1k1 + 2*l1k1*l1p2*l2k1*p1k2 +
			        2*l1l2*l1p2*l2k1*p1k2 + 2*l1k1*l2k1*l2p2*p1k2 + l1k1*l1p2*m2*p1k2 -
			        3*l1p2*l2k1*m2*p1k2 + l1k1*l2p2*m2*p1k2 - l2k1*l2p2*m2*p1k2 -
			        l1k1*m2*p1k1*p1k2 - l1l2*m2*p1k1*p1k2 + m4*p1k1*p1k2 -
			        4*l1k1*l2k2*l2p2*p1l1 - 4*l1l2*l2k2*l2p2*p1l1 +
			        2*l2k1*l2k2*l2p2*p1l1 + 2*l1k2*l2p2*m2*p1l1 +
			        2*l1k1*l2k2*p1k1*p1l1 + 2*l1l2*l2k2*p1k1*p1l1 -
			        2*l2k1*l2k2*p1k1*p1l1 - l1k2*m2*p1k1*p1l1 - 2*l1k1*l2k1*p1k2*p1l1 -
			        2*l1l2*l2k1*p1k2*p1l1 - l1k1*m2*p1k2*p1l1 + 3*l2k1*m2*p1k2*p1l1 -
			        4*l1k1*l1p2*l2k2*p1l2 - 4*l1l2*l1p2*l2k2*p1l2 +
			        2*l1p2*l2k1*l2k2*p1l2 - 4*l1k1*l2k2*l2p2*p1l2 +
			        2*l1k2*l1p2*m2*p1l2 + 2*l1k2*l2k1*p1k1*p1l2 -
			        2*l1l2*l2k2*p1k1*p1l2 + 2*l1k2*m2*p1k1*p1l2 - l2k2*m2*p1k1*p1l2 -
			        2*l1k1*l2k1*p1k2*p1l2 - l1k1*m2*p1k2*p1l2 + l2k1*m2*p1k2*p1l2 +
			        4*l1k1*l2k2*p1l1*p1l2 + 4*l1l2*l2k2*p1l1*p1l2 -
			        2*l2k1*l2k2*p1l1*p1l2 - 2*l1k2*m2*p1l1*p1l2 +
			        2*l1k1*l1k2*l2k1*p1p2 + 2*l1k2*l1l2*l2k1*p1p2 -
			        6*l1k1*l1l2*l2k2*p1p2 - 2*l1k1*l2k1*l2k2*p1p2 +
			        4*l1l2*l2k1*l2k2*p1p2 + 2*l1k1*l1k2*m2*p1p2 +
			        2*l1k2*l1l2*m2*p1p2 - 2*l1k2*l2k1*m2*p1p2 + 8*l1k1*l2k2*m2*p1p2 +
			        12*l1l2*l2k2*m2*p1p2 - 6*l1k2*m4*p1p2 - 2*l1k1*l1p2*l2k1*p2k2 -
			        2*l1l2*l1p2*l2k1*p2k2 - 2*l1k1*l2k1*l2p2*p2k2 - l1k1*l1p2*m2*p2k2 +
			        3*l1p2*l2k1*m2*p2k2 - l1k1*l2p2*m2*p2k2 + l2k1*l2p2*m2*p2k2 +
			        l1k1*m2*p1k1*p2k2 + l1l2*m2*p1k1*p2k2 - m4*p1k1*p2k2 +
			        2*l1k1*l2k1*p1l1*p2k2 + 2*l1l2*l2k1*p1l1*p2k2 + l1k1*m2*p1l1*p2k2 -
			        3*l2k1*m2*p1l1*p2k2 + 2*l1k1*l2k1*p1l2*p2k2 + l1k1*m2*p1l2*p2k2 -
			        l2k1*m2*p1l2*p2k2 + k1p2*
			         (-2*l1p2*l2k1*l2k2 + 2*l1k2*l2k1*l2p2 + k1k2*l1p2*m2 -
			           l1k2*l1p2*m2 + 2*l1k2*l2p2*m2 - l2k2*l2p2*m2 +
			           4*k1k2*m2*p1k1 - 2*l1k2*m2*p1k1 + 4*l2k2*m2*p1k1 - m4*p1k2 +
			           2*l2k1*l2k2*p1l1 - k1k2*m2*p1l1 + l1k2*m2*p1l1 -
			           2*l1k2*l2k1*p1l2 - 2*l1k2*m2*p1l2 + l2k2*m2*p1l2 +
			           l1k1*(2*l1p2*l2k2 - 2*l2k2*p1l1 + m2*(p1k2 - p2k2)) +
			           m4*p2k2 + l1l2*(2*l1p2*l2k2 - 2*k1k2*l2p2 - 2*l2k2*l2p2 +
			              m2*p1k2 - 2*l2k2*p1l1 + 2*k1k2*p1l2 + 2*l2k2*p1l2 - m2*p2k2)\
			) + (-2*k1k2 + l1k2 - 2*l2k2)*m2*pow(k1p2,2) + 2*l2k2*M2*pow(l1k1,2) -
			        2*l2k2*p1p2*pow(l1k1,2) + 4*l2k2*M2*pow(l1l2,2) -
			        4*l2k2*p1p2*pow(l1l2,2) - 2*l1k2*M2*pow(l2k1,2) -
			        2*l1p2*p1k2*pow(l2k1,2) + 2*p1k2*p1l1*pow(l2k1,2) +
			        2*l1k2*p1p2*pow(l2k1,2) + 2*l1p2*p2k2*pow(l2k1,2) -
			        2*p1l1*p2k2*pow(l2k1,2) + 2*l1k1*l2k2*pow(l2p2,2) +
			        l1k2*m2*pow(p1k1,2) - 2*l2k2*m2*pow(p1k1,2) +
			        2*l1k1*l2k2*pow(p1l2,2) +
			        k1k2*(-2*l1p2*l2k1*l2p2 + 2*l1p2*l2p2*m2 - 6*m4*M2 -
			           l1p2*m2*p1k1 + 2*l2k1*l2p2*p1l1 - 2*l1p2*m2*p1l1 -
			           2*l2p2*m2*p1l1 + m2*p1k1*p1l1 + 2*l1p2*l2k1*p1l2 -
			           2*l1p2*m2*p1l2 - 2*l2p2*m2*p1l2 - 2*l2k1*p1l1*p1l2 +
			           2*m2*p1l1*p1l2 + 6*m4*p1p2 +
			           2*l1l2*(-4*m2*M2 + l2p2*p1k1 - l2p2*p1l1 +
			              l1p2*(l2p2 - p1l2) - p1k1*p1l2 + p1l1*p1l2 +
			              l2k1*(M2 - p1p2) + 4*m2*p1p2) +
			           2*(M2 - p1p2)*pow(l1l2,2) + m2*pow(l1p2,2) + m2*pow(l2p2,2) -
			           2*m2*pow(p1k1,2) + m2*pow(p1l1,2) + m2*pow(p1l2,2) +
			           2*l1k1*(l1l2*M2 - 2*m2*M2 + l1p2*(l2p2 - p1l2) + p1l1*p1l2 -
			              l2p2*(p1l1 + 2*p1l2) - l1l2*p1p2 + 2*m2*p1p2 + pow(l2p2,2) +
			              pow(p1l2,2)))) + pow(f2,2)*
			      (12*l1k1*l1p2*l2k2*l2p2*M2 + 12*l1l2*l1p2*l2k2*l2p2*M2 -
			        6*l1p2*l2k1*l2k2*l2p2*M2 - 6*l1k2*l1p2*l2p2*m2*M2 -
			        2*l1k1*l1p2*l2k2*M2*p1k1 - 2*l1l2*l1p2*l2k2*M2*p1k1 +
			        2*l1p2*l2k1*l2k2*M2*p1k1 - 2*l1k2*l2k1*l2p2*M2*p1k1 +
			        2*l1l2*l2k2*l2p2*M2*p1k1 + l1k2*l1p2*m2*M2*p1k1 -
			        2*l1k2*l2p2*m2*M2*p1k1 + l2k2*l2p2*m2*M2*p1k1 +
			        2*l1k1*l1p2*l2k1*M2*p1k2 + 2*l1l2*l1p2*l2k1*M2*p1k2 +
			        2*l1k1*l2k1*l2p2*M2*p1k2 + l1k1*l1p2*m2*M2*p1k2 -
			        3*l1p2*l2k1*m2*M2*p1k2 + l1k1*l2p2*m2*M2*p1k2 -
			        l2k1*l2p2*m2*M2*p1k2 - 3*l1k1*m2*M2*p1k1*p1k2 -
			        3*l1l2*m2*M2*p1k1*p1k2 + 3*m4*M2*p1k1*p1k2 -
			        4*l1k1*l2k2*l2p2*M2*p1l1 - 4*l1l2*l2k2*l2p2*M2*p1l1 +
			        2*l2k1*l2k2*l2p2*M2*p1l1 + 2*l1k2*l2p2*m2*M2*p1l1 +
			        6*l1k1*l2k2*M2*p1k1*p1l1 + 6*l1l2*l2k2*M2*p1k1*p1l1 -
			        6*l2k1*l2k2*M2*p1k1*p1l1 - 3*l1k2*m2*M2*p1k1*p1l1 -
			        6*l1k1*l2k1*M2*p1k2*p1l1 - 6*l1l2*l2k1*M2*p1k2*p1l1 -
			        3*l1k1*m2*M2*p1k2*p1l1 + 9*l2k1*m2*M2*p1k2*p1l1 -
			        4*l1k1*l1p2*l2k2*M2*p1l2 - 4*l1l2*l1p2*l2k2*M2*p1l2 +
			        2*l1p2*l2k1*l2k2*M2*p1l2 - 4*l1k1*l2k2*l2p2*M2*p1l2 +
			        2*l1k2*l1p2*m2*M2*p1l2 + 6*l1k2*l2k1*M2*p1k1*p1l2 -
			        6*l1l2*l2k2*M2*p1k1*p1l2 + 6*l1k2*m2*M2*p1k1*p1l2 -
			        3*l2k2*m2*M2*p1k1*p1l2 - 6*l1k1*l2k1*M2*p1k2*p1l2 -
			        3*l1k1*m2*M2*p1k2*p1l2 + 3*l2k1*m2*M2*p1k2*p1l2 +
			        12*l1k1*l2k2*M2*p1l1*p1l2 + 12*l1l2*l2k2*M2*p1l1*p1l2 -
			        6*l2k1*l2k2*M2*p1l1*p1l2 - 6*l1k2*m2*M2*p1l1*p1l2 -
			        4*l1k1*l1p2*l2k2*l2p2*p1p2 - 4*l1l2*l1p2*l2k2*l2p2*p1p2 +
			        2*l1p2*l2k1*l2k2*l2p2*p1p2 + 2*l1k2*l1p2*l2p2*m2*p1p2 +
			        4*l1k1*l1k2*l2k1*M2*p1p2 + 4*l1k2*l1l2*l2k1*M2*p1p2 -
			        12*l1k1*l1l2*l2k2*M2*p1p2 - 4*l1k1*l2k1*l2k2*M2*p1p2 +
			        8*l1l2*l2k1*l2k2*M2*p1p2 + 4*l1k1*l1k2*m2*M2*p1p2 +
			        4*l1k2*l1l2*m2*M2*p1p2 - 4*l1k2*l2k1*m2*M2*p1p2 +
			        16*l1k1*l2k2*m2*M2*p1p2 + 24*l1l2*l2k2*m2*M2*p1p2 -
			        12*l1k2*m4*M2*p1p2 - 2*l1k1*l1p2*l2k2*p1k1*p1p2 -
			        2*l1l2*l1p2*l2k2*p1k1*p1p2 + 2*l1p2*l2k1*l2k2*p1k1*p1p2 -
			        2*l1k2*l2k1*l2p2*p1k1*p1p2 + 2*l1l2*l2k2*l2p2*p1k1*p1p2 +
			        l1k2*l1p2*m2*p1k1*p1p2 - 2*l1k2*l2p2*m2*p1k1*p1p2 +
			        l2k2*l2p2*m2*p1k1*p1p2 + 2*l1k1*l1p2*l2k1*p1k2*p1p2 +
			        2*l1l2*l1p2*l2k1*p1k2*p1p2 + 2*l1k1*l2k1*l2p2*p1k2*p1p2 +
			        l1k1*l1p2*m2*p1k2*p1p2 - 3*l1p2*l2k1*m2*p1k2*p1p2 +
			        l1k1*l2p2*m2*p1k2*p1p2 - l2k1*l2p2*m2*p1k2*p1p2 +
			        l1k1*m2*p1k1*p1k2*p1p2 + l1l2*m2*p1k1*p1k2*p1p2 -
			        m4*p1k1*p1k2*p1p2 - 4*l1k1*l2k2*l2p2*p1l1*p1p2 -
			        4*l1l2*l2k2*l2p2*p1l1*p1p2 + 2*l2k1*l2k2*l2p2*p1l1*p1p2 +
			        2*l1k2*l2p2*m2*p1l1*p1p2 - 2*l1k1*l2k2*p1k1*p1l1*p1p2 -
			        2*l1l2*l2k2*p1k1*p1l1*p1p2 + 2*l2k1*l2k2*p1k1*p1l1*p1p2 +
			        l1k2*m2*p1k1*p1l1*p1p2 + 2*l1k1*l2k1*p1k2*p1l1*p1p2 +
			        2*l1l2*l2k1*p1k2*p1l1*p1p2 + l1k1*m2*p1k2*p1l1*p1p2 -
			        3*l2k1*m2*p1k2*p1l1*p1p2 - 4*l1k1*l1p2*l2k2*p1l2*p1p2 -
			        4*l1l2*l1p2*l2k2*p1l2*p1p2 + 2*l1p2*l2k1*l2k2*p1l2*p1p2 -
			        4*l1k1*l2k2*l2p2*p1l2*p1p2 + 2*l1k2*l1p2*m2*p1l2*p1p2 -
			        2*l1k2*l2k1*p1k1*p1l2*p1p2 + 2*l1l2*l2k2*p1k1*p1l2*p1p2 -
			        2*l1k2*m2*p1k1*p1l2*p1p2 + l2k2*m2*p1k1*p1l2*p1p2 +
			        2*l1k1*l2k1*p1k2*p1l2*p1p2 + l1k1*m2*p1k2*p1l2*p1p2 -
			        l2k1*m2*p1k2*p1l2*p1p2 - 4*l1k1*l2k2*p1l1*p1l2*p1p2 -
			        4*l1l2*l2k2*p1l1*p1l2*p1p2 + 2*l2k1*l2k2*p1l1*p1l2*p1p2 +
			        2*l1k2*m2*p1l1*p1l2*p1p2 - 6*l1k1*l1p2*l2k1*M2*p2k2 -
			        6*l1l2*l1p2*l2k1*M2*p2k2 - 6*l1k1*l2k1*l2p2*M2*p2k2 -
			        3*l1k1*l1p2*m2*M2*p2k2 + 9*l1p2*l2k1*m2*M2*p2k2 -
			        3*l1k1*l2p2*m2*M2*p2k2 + 3*l2k1*l2p2*m2*M2*p2k2 +
			        l1k1*m2*M2*p1k1*p2k2 + l1l2*m2*M2*p1k1*p2k2 -
			        m4*M2*p1k1*p2k2 + 2*l1k1*l2k1*M2*p1l1*p2k2 +
			        2*l1l2*l2k1*M2*p1l1*p2k2 + l1k1*m2*M2*p1l1*p2k2 -
			        3*l2k1*m2*M2*p1l1*p2k2 + 2*l1k1*l2k1*M2*p1l2*p2k2 +
			        l1k1*m2*M2*p1l2*p2k2 - l2k1*m2*M2*p1l2*p2k2 +
			        2*l1k1*l1p2*l2k1*p1p2*p2k2 + 2*l1l2*l1p2*l2k1*p1p2*p2k2 +
			        2*l1k1*l2k1*l2p2*p1p2*p2k2 + l1k1*l1p2*m2*p1p2*p2k2 -
			        3*l1p2*l2k1*m2*p1p2*p2k2 + l1k1*l2p2*m2*p1p2*p2k2 -
			        l2k1*l2p2*m2*p1p2*p2k2 + l1k1*m2*p1k1*p1p2*p2k2 +
			        l1l2*m2*p1k1*p1p2*p2k2 - m4*p1k1*p1p2*p2k2 +
			        2*l1k1*l2k1*p1l1*p1p2*p2k2 + 2*l1l2*l2k1*p1l1*p1p2*p2k2 +
			        l1k1*m2*p1l1*p1p2*p2k2 - 3*l2k1*m2*p1l1*p1p2*p2k2 +
			        2*l1k1*l2k1*p1l2*p1p2*p2k2 + l1k1*m2*p1l2*p1p2*p2k2 -
			        l2k1*m2*p1l2*p1p2*p2k2 +
			        k1p2*(6*l1l2*l1p2*l2k2*M2 - 6*l1p2*l2k1*l2k2*M2 -
			           6*k1k2*l1l2*l2p2*M2 + 6*l1k2*l2k1*l2p2*M2 -
			           6*l1l2*l2k2*l2p2*M2 + 3*k1k2*l1p2*m2*M2 -
			           3*l1k2*l1p2*m2*M2 + 6*l1k2*l2p2*m2*M2 -
			           3*l2k2*l2p2*m2*M2 + 4*k1k2*m2*M2*p1k1 -
			           2*l1k2*m2*M2*p1k1 + 4*l2k2*m2*M2*p1k1 +
			           l1l2*m2*M2*p1k2 - m4*M2*p1k2 - 2*l1l2*l2k2*M2*p1l1 +
			           2*l2k1*l2k2*M2*p1l1 - k1k2*m2*M2*p1l1 +
			           l1k2*m2*M2*p1l1 + 2*k1k2*l1l2*M2*p1l2 -
			           2*l1k2*l2k1*M2*p1l2 + 2*l1l2*l2k2*M2*p1l2 -
			           2*l1k2*m2*M2*p1l2 + l2k2*m2*M2*p1l2 -
			           2*l1l2*l1p2*l2k2*p1p2 + 2*l1p2*l2k1*l2k2*p1p2 +
			           2*k1k2*l1l2*l2p2*p1p2 - 2*l1k2*l2k1*l2p2*p1p2 +
			           2*l1l2*l2k2*l2p2*p1p2 - k1k2*l1p2*m2*p1p2 +
			           l1k2*l1p2*m2*p1p2 - 2*l1k2*l2p2*m2*p1p2 + l2k2*l2p2*m2*p1p2 +
			           4*k1k2*m2*p1k1*p1p2 - 2*l1k2*m2*p1k1*p1p2 +
			           4*l2k2*m2*p1k1*p1p2 + l1l2*m2*p1k2*p1p2 - m4*p1k2*p1p2 -
			           2*l1l2*l2k2*p1l1*p1p2 + 2*l2k1*l2k2*p1l1*p1p2 -
			           k1k2*m2*p1l1*p1p2 + l1k2*m2*p1l1*p1p2 +
			           2*k1k2*l1l2*p1l2*p1p2 - 2*l1k2*l2k1*p1l2*p1p2 +
			           2*l1l2*l2k2*p1l2*p1p2 - 2*l1k2*m2*p1l2*p1p2 +
			           l2k2*m2*p1l2*p1p2 - 3*l1l2*m2*M2*p2k2 + 3*m4*M2*p2k2 +
			           l1l2*m2*p1p2*p2k2 - m4*p1p2*p2k2 +
			           l1k1*(l1p2*l2k2*(6*M2 - 2*p1p2) - 2*l2k2*p1l1*(M2 + p1p2) +
			              m2*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2)))) -
			        (2*k1k2 - l1k2 + 2*l2k2)*m2*(3*M2 - p1p2)*pow(k1p2,2) -
			        4*l2k2*M2*p1p2*pow(l1k1,2) - 8*l2k2*M2*p1p2*pow(l1l2,2) -
			        2*l1p2*M2*p1k2*pow(l2k1,2) + 6*M2*p1k2*p1l1*pow(l2k1,2) +
			        4*l1k2*M2*p1p2*pow(l2k1,2) - 2*l1p2*p1k2*p1p2*pow(l2k1,2) -
			        2*p1k2*p1l1*p1p2*pow(l2k1,2) + 6*l1p2*M2*p2k2*pow(l2k1,2) -
			        2*M2*p1l1*p2k2*pow(l2k1,2) - 2*l1p2*p1p2*p2k2*pow(l2k1,2) -
			        2*p1l1*p1p2*p2k2*pow(l2k1,2) + 6*l1k1*l2k2*M2*pow(l2p2,2) -
			        2*l1k1*l2k2*p1p2*pow(l2p2,2) - 2*l1k1*l1k2*l2k1*M4 -
			        2*l1k2*l1l2*l2k1*M4 + 6*l1k1*l1l2*l2k2*M4 +
			        6*l1k1*l2k1*l2k2*M4 - 4*l1l2*l2k1*l2k2*M4 -
			        2*l1k1*l1k2*m2*M4 - 2*l1k2*l1l2*m2*M4 -
			        14*l1k1*l2k2*m2*M4 - 20*l1l2*l2k2*m2*M4 +
			        10*l1k2*m4*M4 + 2*l2k2*pow(l1k1,2)*M4 +
			        4*l2k2*pow(l1l2,2)*M4 - 6*l1k2*pow(l2k1,2)*M4 +
			        3*l1k2*m2*M2*pow(p1k1,2) - 6*l2k2*m2*M2*pow(p1k1,2) -
			        l1k2*m2*p1p2*pow(p1k1,2) + 2*l2k2*m2*p1p2*pow(p1k1,2) +
			        6*l1k1*l2k2*M2*pow(p1l2,2) - 2*l1k1*l2k2*p1p2*pow(p1l2,2) -
			        2*l1k1*l1k2*l2k1*pow(p1p2,2) - 2*l1k2*l1l2*l2k1*pow(p1p2,2) +
			        6*l1k1*l1l2*l2k2*pow(p1p2,2) - 2*l1k1*l2k1*l2k2*pow(p1p2,2) -
			        4*l1l2*l2k1*l2k2*pow(p1p2,2) - 2*l1k1*l1k2*m2*pow(p1p2,2) -
			        2*l1k2*l1l2*m2*pow(p1p2,2) + 4*l1k2*l2k1*m2*pow(p1p2,2) -
			        2*l1k1*l2k2*m2*pow(p1p2,2) - 4*l1l2*l2k2*m2*pow(p1p2,2) +
			        2*l1k2*m4*pow(p1p2,2) + 2*l2k2*pow(l1k1,2)*pow(p1p2,2) +
			        4*l2k2*pow(l1l2,2)*pow(p1p2,2) + 2*l1k2*pow(l2k1,2)*pow(p1p2,2) +
			        k1k2*(-6*l1p2*l2k1*l2p2*M2 + 6*l1p2*l2p2*m2*M2 -
			           l1p2*m2*M2*p1k1 + 2*l2k1*l2p2*M2*p1l1 -
			           2*l1p2*m2*M2*p1l1 - 2*l2p2*m2*M2*p1l1 +
			           3*m2*M2*p1k1*p1l1 + 2*l1p2*l2k1*M2*p1l2 -
			           2*l1p2*m2*M2*p1l2 - 2*l2p2*m2*M2*p1l2 -
			           6*l2k1*M2*p1l1*p1l2 + 6*m2*M2*p1l1*p1l2 +
			           2*l1p2*l2k1*l2p2*p1p2 - 2*l1p2*l2p2*m2*p1p2 + 12*m4*M2*p1p2 -
			           l1p2*m2*p1k1*p1p2 + 2*l2k1*l2p2*p1l1*p1p2 -
			           2*l1p2*m2*p1l1*p1p2 - 2*l2p2*m2*p1l1*p1p2 -
			           m2*p1k1*p1l1*p1p2 + 2*l1p2*l2k1*p1l2*p1p2 -
			           2*l1p2*m2*p1l2*p1p2 - 2*l2p2*m2*p1l2*p1p2 +
			           2*l2k1*p1l1*p1l2*p1p2 - 2*m2*p1l1*p1l2*p1p2 +
			           3*m2*M2*pow(l1p2,2) - m2*p1p2*pow(l1p2,2) +
			           3*m2*M2*pow(l2p2,2) - m2*p1p2*pow(l2p2,2) -
			           12*m4*M4 - 6*m2*M2*pow(p1k1,2) +
			           2*m2*p1p2*pow(p1k1,2) + 3*m2*M2*pow(p1l1,2) -
			           m2*p1p2*pow(p1l1,2) + 3*m2*M2*pow(p1l2,2) -
			           m2*p1p2*pow(p1l2,2) + 2*pow(l1l2,2)*pow(M2 - p1p2,2) +
			           2*l1k1*(3*M2*p1l1*p1l2 - 2*l1l2*M2*p1p2 + 4*m2*M2*p1p2 -
			              p1l1*p1l2*p1p2 - l2p2*(p1l1 + 2*p1l2)*(M2 + p1p2) +
			              l1p2*(l2p2*(3*M2 - p1p2) - p1l2*(M2 + p1p2)) +
			              (3*M2 - p1p2)*pow(l2p2,2) + l1l2*M4 -
			              4*m2*M4 + 3*M2*pow(p1l2,2) - p1p2*pow(p1l2,2) +
			              l1l2*pow(p1p2,2)) -
			           2*l1l2*(-(l2p2*M2*p1k1) + l2p2*M2*p1l1 + 3*M2*p1k1*p1l2 -
			              3*M2*p1l1*p1l2 - 8*m2*M2*p1p2 - l2p2*p1k1*p1p2 +
			              l2p2*p1l1*p1p2 - p1k1*p1l2*p1p2 + p1l1*p1l2*p1p2 +
			              l1p2*(l2p2*(-3*M2 + p1p2) + p1l2*(M2 + p1p2)) +
			              7*m2*M4 + m2*pow(p1p2,2) +
			              l2k1*(2*M2*p1p2 - 3*M4 + pow(p1p2,2))))) +
			     eg1*(-4*M2*pow(f1,2)*(-2*k1p2*l1k1*l2k2*p1k2 +
			           2*l1p2*l2k1*l2k2*p1k2 - 4*l1k1*l2k2*l2p2*p1k2 -
			           2*l1l2*l2k2*l2p2*p1k2 - 2*k1p2*l2k2*m2*p1k2 +
			           l1p2*l2k2*m2*p1k2 + 2*l2k2*l2p2*m2*p1k2 -
			           2*l1k1*l2k2*p1k1*p2k2 - 2*l2k2*m2*p1k1*p2k2 +
			           4*l1k1*l2k1*p1k2*p2k2 + 2*l1k1*m2*p1k2*p2k2 +
			           2*l1l2*m2*p1k2*p2k2 - 4*l2k1*m2*p1k2*p2k2 - 2*m4*p1k2*p2k2 +
			           2*l2k1*l2k2*p1l1*p2k2 + l2k2*m2*p1l1*p2k2 -
			           4*l1k1*l2k2*p1l2*p2k2 - 2*l1l2*l2k2*p1l2*p2k2 +
			           2*l2k2*m2*p1l2*p2k2 +
			           k1k2*(6*l1l2*l2k2*M2 - 4*l2k2*m2*M2 - 2*l1p2*l2k2*p1k1 +
			              2*k1p2*l1l2*p1k2 - 4*k1p2*m2*p1k2 + l1p2*m2*p1k2 +
			              2*l2p2*m2*p1k2 - 2*k1p2*l2k2*p1l1 - 2*l2k2*l2p2*p1l1 -
			              2*l1p2*l2k2*p1l2 - 2*l1l2*l2k2*p1p2 + 2*l2k2*m2*p1p2 +
			              2*l1k2*l2k1*(-M2 + p1p2) + 2*l1l2*p1k1*p2k2 -
			              4*m2*p1k1*p2k2 + m2*p1l1*p2k2 + 2*m2*p1l2*p2k2 +
			              2*l1k1*(-(l2p2*p1k2) + l2k2*(M2 + p1p2) - p1l2*p2k2)) -
			           l1k2*(2*l2k2*(m2*M2 - 2*(k1p2 + l2p2)*(p1k1 + p1l2)) +
			              2*l2k1*(k1p2*p1k2 + l2k2*(M2 + p1p2) + p1k1*p2k2) +
			              m2*(k1p2*p1k2 + l2p2*p1k2 + (p1k1 + p1l2)*p2k2)) +
			           2*(l1l2 - 2*m2)*(M2 - p1p2)*pow(k1k2,2) +
			           2*l1k1*M2*pow(l2k2,2) + 4*l1l2*M2*pow(l2k2,2) -
			           2*m2*M2*pow(l2k2,2) - 2*l1p2*p1k1*pow(l2k2,2) -
			           2*k1p2*p1l1*pow(l2k2,2) - 2*l2p2*p1l1*pow(l2k2,2) -
			           2*l1p2*p1l2*pow(l2k2,2) + 2*l1k1*p1p2*pow(l2k2,2)) +
			        4*f1*f2*M2*(2*l1k2*l2k2*m2*M2 - 4*l1k2*l2k2*l2p2*p1k1 -
			           2*l1p2*l2k1*l2k2*p1k2 + 4*l1k1*l2k2*l2p2*p1k2 +
			           2*l1l2*l2k2*l2p2*p1k2 - l1p2*l2k2*m2*p1k2 +
			           l1k2*l2p2*m2*p1k2 - 2*l2k2*l2p2*m2*p1k2 -
			           2*l1k2*l2k1*p1k1*p1k2 - 2*l1k1*l2k2*p1k1*p1k2 -
			           l1k2*m2*p1k1*p1k2 - 2*l2k2*m2*p1k1*p1k2 +
			           2*l2k1*l2k2*p1k2*p1l1 + l2k2*m2*p1k2*p1l1 -
			           4*l1k2*l2k2*l2p2*p1l2 + 4*l1k2*l2k2*p1k1*p1l2 -
			           4*l1k1*l2k2*p1k2*p1l2 - 2*l1l2*l2k2*p1k2*p1l2 -
			           l1k2*m2*p1k2*p1l2 + 2*l2k2*m2*p1k2*p1l2 -
			           2*l1k2*l2k2*m2*p1p2 + 2*l1p2*l2k1*l2k2*p2k2 -
			           4*l1k1*l2k2*l2p2*p2k2 - 2*l1l2*l2k2*l2p2*p2k2 +
			           l1p2*l2k2*m2*p2k2 - l1k2*l2p2*m2*p2k2 + 2*l2k2*l2p2*m2*p2k2 +
			           2*l1k2*l2k1*p1k1*p2k2 + 2*l1k1*l2k2*p1k1*p2k2 +
			           l1k2*m2*p1k1*p2k2 + 2*l2k2*m2*p1k1*p2k2 -
			           4*l1k1*l2k1*p1k2*p2k2 - 2*l1k1*m2*p1k2*p2k2 -
			           2*l1l2*m2*p1k2*p2k2 + 4*l2k1*m2*p1k2*p2k2 + 2*m4*p1k2*p2k2 -
			           2*l2k1*l2k2*p1l1*p2k2 - l2k2*m2*p1l1*p2k2 +
			           4*l1k1*l2k2*p1l2*p2k2 + 2*l1l2*l2k2*p1l2*p2k2 +
			           l1k2*m2*p1l2*p2k2 - 2*l2k2*m2*p1l2*p2k2 +
			           k1k2*(-8*l1l2*l2k2*M2 + 6*l2k2*m2*M2 + 2*l1k1*l2p2*p1k2 -
			              2*l2p2*m2*p1k2 + 2*l1l2*p1k1*p1k2 - 4*m2*p1k1*p1k2 +
			              2*l2k2*l2p2*p1l1 - 2*l2k2*p1k1*p1l1 + m2*p1k2*p1l1 -
			              2*l1k1*p1k2*p1l2 + 2*m2*p1k2*p1l2 - 2*l2k2*p1l1*p1l2 +
			              4*l1k2*l2k1*(M2 - p1p2) + 8*l1l2*l2k2*p1p2 -
			              6*l2k2*m2*p1p2 - 2*l1k1*l2p2*p2k2 + 2*l2p2*m2*p2k2 -
			              2*l1l2*p1k1*p2k2 + 4*m2*p1k1*p2k2 - m2*p1l1*p2k2 +
			              2*l1k1*p1l2*p2k2 - 2*m2*p1l2*p2k2 +
			              l1p2*(2*l2k2*(-l2p2 + p1k1 + p1l2) + m2*(-p1k2 + p2k2))) -
			           4*(l1l2 - 2*m2)*(M2 - p1p2)*pow(k1k2,2) +
			           2*l1k2*l2k2*pow(k1p2,2) - 2*l1p2*l2p2*pow(l2k2,2) -
			           4*l1l2*M2*pow(l2k2,2) + 2*m2*M2*pow(l2k2,2) +
			           2*l1p2*p1k1*pow(l2k2,2) + 2*l2p2*p1l1*pow(l2k2,2) -
			           2*p1k1*p1l1*pow(l2k2,2) + 2*l1p2*p1l2*pow(l2k2,2) -
			           2*p1l1*p1l2*pow(l2k2,2) + 4*l1l2*p1p2*pow(l2k2,2) -
			           2*m2*p1p2*pow(l2k2,2) +
			           k1p2*(4*l1k2*l2k2*l2p2 - 4*l1k2*l2k2*p1k1 + 2*l1k2*l2k1*p1k2 +
			              2*l1k1*l2k2*p1k2 + l1k2*m2*p1k2 + 2*l2k2*m2*p1k2 -
			              4*l1k2*l2k2*p1l2 - 2*l1k2*l2k1*p2k2 - 2*l1k1*l2k2*p2k2 -
			              l1k2*m2*p2k2 - 2*l2k2*m2*p2k2 -
			              2*k1k2*(l1p2*l2k2 + l1l2*p1k2 - 2*m2*p1k2 - l2k2*p1l1 -
			                 l1l2*p2k2 + 2*m2*p2k2) - 2*l1p2*pow(l2k2,2) +
			              2*p1l1*pow(l2k2,2)) + 2*l1k2*l2k2*pow(l2p2,2) +
			           2*l1k2*l2k2*pow(p1k1,2) + 2*l1k1*l2k1*pow(p1k2,2) +
			           l1k1*m2*pow(p1k2,2) + l1l2*m2*pow(p1k2,2) -
			           2*l2k1*m2*pow(p1k2,2) - m4*pow(p1k2,2) +
			           2*l1k2*l2k2*pow(p1l2,2) + 2*l1k1*l2k1*pow(p2k2,2) +
			           l1k1*m2*pow(p2k2,2) + l1l2*m2*pow(p2k2,2) -
			           2*l2k1*m2*pow(p2k2,2) - m4*pow(p2k2,2)) +
			        pow(f2,2)*(-4*l1k2*l2k2*l2p2*M2*p1k1 -
			           2*l1p2*l2k1*l2k2*M2*p1k2 + 4*l1k1*l2k2*l2p2*M2*p1k2 +
			           2*l1l2*l2k2*l2p2*M2*p1k2 - l1p2*l2k2*m2*M2*p1k2 +
			           l1k2*l2p2*m2*M2*p1k2 - 2*l2k2*l2p2*m2*M2*p1k2 -
			           6*l1k2*l2k1*M2*p1k1*p1k2 - 6*l1k1*l2k2*M2*p1k1*p1k2 -
			           3*l1k2*m2*M2*p1k1*p1k2 - 6*l2k2*m2*M2*p1k1*p1k2 +
			           6*l2k1*l2k2*M2*p1k2*p1l1 + 3*l2k2*m2*M2*p1k2*p1l1 -
			           4*l1k2*l2k2*l2p2*M2*p1l2 + 12*l1k2*l2k2*M2*p1k1*p1l2 -
			           12*l1k1*l2k2*M2*p1k2*p1l2 - 6*l1l2*l2k2*M2*p1k2*p1l2 -
			           3*l1k2*m2*M2*p1k2*p1l2 + 6*l2k2*m2*M2*p1k2*p1l2 -
			           4*l1k2*l2k2*m2*M2*p1p2 - 4*l1k2*l2k2*l2p2*p1k1*p1p2 -
			           2*l1p2*l2k1*l2k2*p1k2*p1p2 + 4*l1k1*l2k2*l2p2*p1k2*p1p2 +
			           2*l1l2*l2k2*l2p2*p1k2*p1p2 - l1p2*l2k2*m2*p1k2*p1p2 +
			           l1k2*l2p2*m2*p1k2*p1p2 - 2*l2k2*l2p2*m2*p1k2*p1p2 +
			           2*l1k2*l2k1*p1k1*p1k2*p1p2 + 2*l1k1*l2k2*p1k1*p1k2*p1p2 +
			           l1k2*m2*p1k1*p1k2*p1p2 + 2*l2k2*m2*p1k1*p1k2*p1p2 -
			           2*l2k1*l2k2*p1k2*p1l1*p1p2 - l2k2*m2*p1k2*p1l1*p1p2 -
			           4*l1k2*l2k2*l2p2*p1l2*p1p2 - 4*l1k2*l2k2*p1k1*p1l2*p1p2 +
			           4*l1k1*l2k2*p1k2*p1l2*p1p2 + 2*l1l2*l2k2*p1k2*p1l2*p1p2 +
			           l1k2*m2*p1k2*p1l2*p1p2 - 2*l2k2*m2*p1k2*p1l2*p1p2 +
			           6*l1p2*l2k1*l2k2*M2*p2k2 - 12*l1k1*l2k2*l2p2*M2*p2k2 -
			           6*l1l2*l2k2*l2p2*M2*p2k2 + 3*l1p2*l2k2*m2*M2*p2k2 -
			           3*l1k2*l2p2*m2*M2*p2k2 + 6*l2k2*l2p2*m2*M2*p2k2 +
			           2*l1k2*l2k1*M2*p1k1*p2k2 + 2*l1k1*l2k2*M2*p1k1*p2k2 +
			           l1k2*m2*M2*p1k1*p2k2 + 2*l2k2*m2*M2*p1k1*p2k2 -
			           4*l1k1*l2k1*M2*p1k2*p2k2 - 2*l1k1*m2*M2*p1k2*p2k2 -
			           2*l1l2*m2*M2*p1k2*p2k2 + 4*l2k1*m2*M2*p1k2*p2k2 +
			           2*m4*M2*p1k2*p2k2 - 2*l2k1*l2k2*M2*p1l1*p2k2 -
			           l2k2*m2*M2*p1l1*p2k2 + 4*l1k1*l2k2*M2*p1l2*p2k2 +
			           2*l1l2*l2k2*M2*p1l2*p2k2 + l1k2*m2*M2*p1l2*p2k2 -
			           2*l2k2*m2*M2*p1l2*p2k2 - 2*l1p2*l2k1*l2k2*p1p2*p2k2 +
			           4*l1k1*l2k2*l2p2*p1p2*p2k2 + 2*l1l2*l2k2*l2p2*p1p2*p2k2 -
			           l1p2*l2k2*m2*p1p2*p2k2 + l1k2*l2p2*m2*p1p2*p2k2 -
			           2*l2k2*l2p2*m2*p1p2*p2k2 + 2*l1k2*l2k1*p1k1*p1p2*p2k2 +
			           2*l1k1*l2k2*p1k1*p1p2*p2k2 + l1k2*m2*p1k1*p1p2*p2k2 +
			           2*l2k2*m2*p1k1*p1p2*p2k2 - 4*l1k1*l2k1*p1k2*p1p2*p2k2 -
			           2*l1k1*m2*p1k2*p1p2*p2k2 - 2*l1l2*m2*p1k2*p1p2*p2k2 +
			           4*l2k1*m2*p1k2*p1p2*p2k2 + 2*m4*p1k2*p1p2*p2k2 -
			           2*l2k1*l2k2*p1l1*p1p2*p2k2 - l2k2*m2*p1l1*p1p2*p2k2 +
			           4*l1k1*l2k2*p1l2*p1p2*p2k2 + 2*l1l2*l2k2*p1l2*p1p2*p2k2 +
			           l1k2*m2*p1l2*p1p2*p2k2 - 2*l2k2*m2*p1l2*p1p2*p2k2 -
			           8*(l1l2 - 2*m2)*M2*(M2 - p1p2)*pow(k1k2,2) +
			           2*l1k2*l2k2*(3*M2 - p1p2)*pow(k1p2,2) -
			           6*l1p2*l2p2*M2*pow(l2k2,2) + 2*l1p2*M2*p1k1*pow(l2k2,2) +
			           2*l2p2*M2*p1l1*pow(l2k2,2) - 6*M2*p1k1*p1l1*pow(l2k2,2) +
			           2*l1p2*M2*p1l2*pow(l2k2,2) - 6*M2*p1l1*p1l2*pow(l2k2,2) +
			           2*l1p2*l2p2*p1p2*pow(l2k2,2) + 8*l1l2*M2*p1p2*pow(l2k2,2) -
			           4*m2*M2*p1p2*pow(l2k2,2) + 2*l1p2*p1k1*p1p2*pow(l2k2,2) +
			           2*l2p2*p1l1*p1p2*pow(l2k2,2) + 2*p1k1*p1l1*p1p2*pow(l2k2,2) +
			           2*l1p2*p1l2*p1p2*pow(l2k2,2) + 2*p1l1*p1l2*p1p2*pow(l2k2,2) +
			           k1p2*(12*l1k2*l2k2*l2p2*M2 - 4*l1k2*l2k2*M2*p1k1 +
			              2*l1k2*l2k1*M2*p1k2 + 2*l1k1*l2k2*M2*p1k2 +
			              l1k2*m2*M2*p1k2 + 2*l2k2*m2*M2*p1k2 -
			              4*l1k2*l2k2*M2*p1l2 - 4*l1k2*l2k2*l2p2*p1p2 -
			              4*l1k2*l2k2*p1k1*p1p2 + 2*l1k2*l2k1*p1k2*p1p2 +
			              2*l1k1*l2k2*p1k2*p1p2 + l1k2*m2*p1k2*p1p2 +
			              2*l2k2*m2*p1k2*p1p2 - 4*l1k2*l2k2*p1l2*p1p2 -
			              6*l1k2*l2k1*M2*p2k2 - 6*l1k1*l2k2*M2*p2k2 -
			              3*l1k2*m2*M2*p2k2 - 6*l2k2*m2*M2*p2k2 +
			              2*l1k2*l2k1*p1p2*p2k2 + 2*l1k1*l2k2*p1p2*p2k2 +
			              l1k2*m2*p1p2*p2k2 + 2*l2k2*m2*p1p2*p2k2 -
			              2*k1k2*(-2*m2*M2*p1k2 - l2k2*M2*p1l1 +
			                 l1p2*l2k2*(3*M2 - p1p2) - 2*m2*p1k2*p1p2 -
			                 l2k2*p1l1*p1p2 + 6*m2*M2*p2k2 - 2*m2*p1p2*p2k2 +
			                 l1l2*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2))) -
			              6*l1p2*M2*pow(l2k2,2) + 2*M2*p1l1*pow(l2k2,2) +
			              2*l1p2*p1p2*pow(l2k2,2) + 2*p1l1*p1p2*pow(l2k2,2)) +
			           6*l1k2*l2k2*M2*pow(l2p2,2) - 2*l1k2*l2k2*p1p2*pow(l2p2,2) -
			           4*l1k2*l2k1*l2k2*M4 + 2*l1k2*l2k2*m2*M4 +
			           4*l1k1*pow(l2k2,2)*M4 - 4*l1l2*pow(l2k2,2)*M4 +
			           2*m2*pow(l2k2,2)*M4 + 6*l1k2*l2k2*M2*pow(p1k1,2) -
			           2*l1k2*l2k2*p1p2*pow(p1k1,2) + 6*l1k1*l2k1*M2*pow(p1k2,2) +
			           3*l1k1*m2*M2*pow(p1k2,2) + 3*l1l2*m2*M2*pow(p1k2,2) -
			           6*l2k1*m2*M2*pow(p1k2,2) - 3*m4*M2*pow(p1k2,2) -
			           2*l1k1*l2k1*p1p2*pow(p1k2,2) - l1k1*m2*p1p2*pow(p1k2,2) -
			           l1l2*m2*p1p2*pow(p1k2,2) + 2*l2k1*m2*p1p2*pow(p1k2,2) +
			           m4*p1p2*pow(p1k2,2) + 6*l1k2*l2k2*M2*pow(p1l2,2) -
			           2*l1k2*l2k2*p1p2*pow(p1l2,2) + 4*l1k2*l2k1*l2k2*pow(p1p2,2) +
			           2*l1k2*l2k2*m2*pow(p1p2,2) - 4*l1k1*pow(l2k2,2)*pow(p1p2,2) -
			           4*l1l2*pow(l2k2,2)*pow(p1p2,2) + 2*m2*pow(l2k2,2)*pow(p1p2,2) -
			           k1k2*(-2*l1k1*l2p2*M2*p1k2 + 2*l2p2*m2*M2*p1k2 -
			              6*l1l2*M2*p1k1*p1k2 + 12*m2*M2*p1k1*p1k2 -
			              2*l2k2*l2p2*M2*p1l1 + 6*l2k2*M2*p1k1*p1l1 -
			              3*m2*M2*p1k2*p1l1 + 6*l1k1*M2*p1k2*p1l2 -
			              6*m2*M2*p1k2*p1l2 + 6*l2k2*M2*p1l1*p1l2 -
			              16*l1l2*l2k2*M2*p1p2 + 12*l2k2*m2*M2*p1p2 -
			              2*l1k1*l2p2*p1k2*p1p2 + 2*l2p2*m2*p1k2*p1p2 +
			              2*l1l2*p1k1*p1k2*p1p2 - 4*m2*p1k1*p1k2*p1p2 -
			              2*l2k2*l2p2*p1l1*p1p2 - 2*l2k2*p1k1*p1l1*p1p2 +
			              m2*p1k2*p1l1*p1p2 - 2*l1k1*p1k2*p1l2*p1p2 +
			              2*m2*p1k2*p1l2*p1p2 - 2*l2k2*p1l1*p1l2*p1p2 +
			              8*l1k2*l2k1*M2*(-M2 + p1p2) + 6*l1k1*l2p2*M2*p2k2 -
			              6*l2p2*m2*M2*p2k2 + 2*l1l2*M2*p1k1*p2k2 -
			              4*m2*M2*p1k1*p2k2 + m2*M2*p1l1*p2k2 -
			              2*l1k1*M2*p1l2*p2k2 + 2*m2*M2*p1l2*p2k2 -
			              2*l1k1*l2p2*p1p2*p2k2 + 2*l2p2*m2*p1p2*p2k2 +
			              2*l1l2*p1k1*p1p2*p2k2 - 4*m2*p1k1*p1p2*p2k2 +
			              m2*p1l1*p1p2*p2k2 - 2*l1k1*p1l2*p1p2*p2k2 +
			              2*m2*p1l2*p1p2*p2k2 +
			              l1p2*(l2k2*(l2p2*(6*M2 - 2*p1p2) -
			                    2*(p1k1 + p1l2)*(M2 + p1p2)) +
			                 m2*(M2*(p1k2 - 3*p2k2) + p1p2*(p1k2 + p2k2))) -
			              4*l1k1*l2k2*M4 + 12*l1l2*l2k2*M4 -
			              10*l2k2*m2*M4 + 4*l1k1*l2k2*pow(p1p2,2) +
			              4*l1l2*l2k2*pow(p1p2,2) - 2*l2k2*m2*pow(p1p2,2)) +
			           6*l1k1*l2k1*M2*pow(p2k2,2) + 3*l1k1*m2*M2*pow(p2k2,2) +
			           3*l1l2*m2*M2*pow(p2k2,2) - 6*l2k1*m2*M2*pow(p2k2,2) -
			           3*m4*M2*pow(p2k2,2) - 2*l1k1*l2k1*p1p2*pow(p2k2,2) -
			           l1k1*m2*p1p2*pow(p2k2,2) - l1l2*m2*p1p2*pow(p2k2,2) +
			           2*l2k1*m2*p1p2*pow(p2k2,2) + m4*p1p2*pow(p2k2,2))) +
			     2*(k1k2 + l2k2)*m2*pow(eg1,2)*
			      (8*M2*p1k2*p2k2*pow(f1,2) - 4*f1*f2*M2*pow(p1k2 - p2k2,2) +
			        pow(f2,2)*(M2*(2*p1k2*p2k2 - 3*pow(p1k2,2) - 3*pow(p2k2,2)) +
			           p1p2*pow(p1k2 + p2k2,2))));
}

long double Melem::melem2_Rinterf(const long double l1k1, const long double l1k2, const long double l2k1, const long double l2k2,
			const long double Q2e, const long double Q2h, const long double S, const long double Sk, const long double Sq2,
			const long double f1, const long double f2)const {

	return -32*pow(f1,2) + (pow(l2k1,-2)*(-64*f1*f2*l1k1*m2 -
		       64*f1*f2*Q2e*m2 - 64*l1k1*pow(f1,2)*m2 +
		       32*S*pow(f1,2)*m2 - 64*Sk*pow(f1,2)*m2 +
		       32*Sq2*pow(f1,2)*m2 - 32*l1k1*pow(f2,2)*m2 -
		       32*Q2e*pow(f2,2)*m2 + 64*f1*f2*m4 +
		       32*pow(f2,2)*m4))/2. +
		  (pow(l1k1,-2)*(64*f1*f2*l2k1*m2 - 64*f1*f2*Q2e*m2 +
		       64*l2k1*pow(f1,2)*m2 - 32*S*pow(f1,2)*m2 -
		       64*Sk*pow(f1,2)*m2 + 32*l2k1*pow(f2,2)*m2 -
		       32*Q2e*pow(f2,2)*m2 + 64*f1*f2*m4 +
		       32*pow(f2,2)*m4))/2. - 8*Q2h*pow(f2,2)/M2 +
		  (pow(l1k1,-1)*(16*Q2h*S*pow(f2,2) - 8*Q2h*Sk*pow(f2,2) +
		       24*Q2h*Sq2*pow(f2,2) - 20*Q2h*pow(f2,2)*m2)/M2)/2. +
		  (pow(l1k2,-1)*(8*l1k1*Q2h*pow(f2,2) - 16*Q2e*Q2h*pow(f2,2) -
		       32*Q2h*S*pow(f2,2) - 24*Q2h*Sq2*pow(f2,2) -
		       16*Q2h*pow(f2,2)*m2)/M2)/2. +
		  (pow(l2k2,-1)*(8*l2k1*Q2h*pow(f2,2) + 16*Q2e*Q2h*pow(f2,2) -
		       32*Q2h*S*pow(f2,2) - 8*Q2h*Sq2*pow(f2,2) +
		       16*Q2h*pow(f2,2)*m2)/M2)/2. +
		  (pow(l2k1,-1)*(16*Q2h*S*pow(f2,2) + 8*Q2h*Sk*pow(f2,2) -
		       8*Q2h*Sq2*pow(f2,2) + 20*Q2h*pow(f2,2)*m2)/M2)/2. +
		  (pow(l1k1,-2)*(8*l2k1*Q2h*pow(f2,2)*m2 +
		       8*Q2e*Q2h*pow(f2,2)*m2 - 8*Q2h*S*pow(f2,2)*m2 -
		       16*Q2h*Sk*pow(f2,2)*m2 - 8*Q2h*pow(f2,2)*m4)*
		     pow(M,-2))/2. + (pow(l2k1,-2)*
		     (-8*l1k1*Q2h*pow(f2,2)*m2 + 8*Q2e*Q2h*pow(f2,2)*m2 +
		       8*Q2h*S*pow(f2,2)*m2 - 16*Q2h*Sk*pow(f2,2)*m2 +
		       8*Q2h*Sq2*pow(f2,2)*m2 - 8*Q2h*pow(f2,2)*m4)*
		     pow(M,-2))/2. + (pow(l1k2,-1)*
		     (128*f1*f2*Q2e + 128*f1*f2*Q2h + 32*l1k1*pow(f1,2) +
		       64*Q2h*pow(f1,2) - 128*S*pow(f1,2) - 96*Sq2*pow(f1,2) +
		       64*Q2e*pow(f2,2) + 32*Q2h*pow(f2,2) - 64*pow(f1,2)*m2 -
		       128*pow(f1,2)*M2))/2. +
		  (pow(l2k1,-1)*(64*f1*f2*l1k1 - 32*f1*f2*Q2e + 32*f1*f2*Q2h +
		       32*l1k1*pow(f1,2) - 16*Q2e*pow(f1,2) + 16*Q2h*pow(f1,2) +
		       64*S*pow(f1,2) + 32*Sk*pow(f1,2) - 32*Sq2*pow(f1,2) +
		       32*l1k1*pow(f2,2) - 16*Q2e*pow(f2,2) + 8*Q2h*pow(f2,2) +
		       64*f1*f2*m2 + 112*pow(f1,2)*m2 +
		       32*pow(f2,2)*m2 - 32*pow(f1,2)*M2))/2. +
		  (pow(l1k1,-1)*(64*f1*f2*l2k1 + 32*f1*f2*Q2e - 32*f1*f2*Q2h +
		       32*l2k1*pow(f1,2) + 16*Q2e*pow(f1,2) - 16*Q2h*pow(f1,2) +
		       64*S*pow(f1,2) - 32*Sk*pow(f1,2) + 96*Sq2*pow(f1,2) +
		       32*l2k1*pow(f2,2) + 16*Q2e*pow(f2,2) - 8*Q2h*pow(f2,2) -
		       64*f1*f2*m2 - 112*pow(f1,2)*m2 -
		       32*pow(f2,2)*m2 + 32*pow(f1,2)*M2))/2. +
		  (pow(l2k2,-1)*(-128*f1*f2*Q2e - 128*f1*f2*Q2h + 32*l2k1*pow(f1,2) -
		       64*Q2h*pow(f1,2) - 128*S*pow(f1,2) - 32*Sq2*pow(f1,2) -
		       64*Q2e*pow(f2,2) - 32*Q2h*pow(f2,2) + 64*pow(f1,2)*m2 +
		       128*pow(f1,2)*M2))/2. +
		  (pow(l2k1,-1)*pow(l2k2,-1)*(16*f1*f2*Q2e*Q2h - 16*l1k1*Q2e*pow(f1,2) -
		       16*l1k2*Q2e*pow(f1,2) + 8*Q2e*Q2h*pow(f1,2) - 32*l1k1*S*pow(f1,2) -
		       32*l1k2*S*pow(f1,2) - 48*Q2e*S*pow(f1,2) - 16*Q2e*Sk*pow(f1,2) +
		       32*Q2e*Sq2*pow(f1,2) + 64*S*Sq2*pow(f1,2) + 4*Q2e*Q2h*pow(f2,2) -
		       32*f1*f2*Q2e*m2 + 32*f1*f2*Q2h*m2 -
		       16*l1k2*pow(f1,2)*m2 - 40*Q2e*pow(f1,2)*m2 +
		       24*Q2h*pow(f1,2)*m2 + 64*S*pow(f1,2)*m2 +
		       32*Sk*pow(f1,2)*m2 + 32*Sq2*pow(f1,2)*m2 -
		       16*Q2e*pow(f2,2)*m2 + 16*Q2h*pow(f2,2)*m2 +
		       64*pow(f1,2)*m4 - 16*Q2e*pow(f1,2)*M2 -
		       24*pow(f1,2)*pow(Q2e,2)))/2. +
		  (pow(l1k1,-1)*pow(l2k1,-1)*(-32*Q2e*Sk*pow(f1,2) + 128*f1*f2*m4 +
		       64*pow(f2,2)*m4 - 64*f1*f2*pow(Q2e,2) -
		       16*pow(f1,2)*pow(Q2e,2) - 32*pow(f2,2)*pow(Q2e,2)))/2. +
		  (pow(l1k1,-1)*pow(l2k1,-1)/M2*
		     (-8*Q2e*Q2h*Sk*pow(f2,2) - 16*Q2h*pow(f2,2)*m4 +
		       4*Q2h*pow(f2,2)*pow(Q2e,2)))/2. +
		  ((-128*f1*f2*Q2e + 64*f1*f2*Q2h + 192*S*pow(f1,2) + 192*Sq2*pow(f1,2) -
		       64*Q2e*pow(f2,2) + 48*Q2h*pow(f2,2) - 384*pow(f1,2)*m2 +
		       64*pow(f1,2)*M2)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1))/2. +
		  (pow(M,-2)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (16*Q2e*Q2h*pow(f2,2) + 48*Q2h*S*pow(f2,2) + 48*Q2h*Sq2*pow(f2,2) -
		       96*Q2h*pow(f2,2)*m2 - 8*pow(f2,2)*pow(Q2h,2)))/2. +
		  (pow(l2k1,-1)*pow(l2k2,-1)/M2*
		     (-4*l1k1*Q2e*Q2h*pow(f2,2) - 4*l1k2*Q2e*Q2h*pow(f2,2) -
		       8*l1k1*Q2h*S*pow(f2,2) - 8*l1k2*Q2h*S*pow(f2,2) -
		       12*Q2e*Q2h*S*pow(f2,2) - 4*Q2e*Q2h*Sk*pow(f2,2) +
		       8*Q2e*Q2h*Sq2*pow(f2,2) + 16*Q2h*S*Sq2*pow(f2,2) -
		       4*l1k2*Q2h*pow(f2,2)*m2 - 6*Q2e*Q2h*pow(f2,2)*m2 +
		       16*Q2h*S*pow(f2,2)*m2 + 8*Q2h*Sk*pow(f2,2)*m2 +
		       8*Q2h*Sq2*pow(f2,2)*m2 + 16*Q2h*pow(f2,2)*m4 -
		       6*Q2h*pow(f2,2)*pow(Q2e,2) + 2*pow(f2,2)*m2*pow(Q2h,2)))/2. +
		  (pow(l2k1,-2)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (-64*f1*f2*l1k2*Q2e*m2 + 64*f1*f2*Q2e*Q2h*m2 -
		       64*l1k2*Q2e*pow(f1,2)*m2 + 32*Q2e*Q2h*pow(f1,2)*m2 -
		       64*l1k2*S*pow(f1,2)*m2 - 64*Q2e*S*pow(f1,2)*m2 -
		       64*Q2h*S*pow(f1,2)*m2 + 64*Q2e*Sk*pow(f1,2)*m2 +
		       32*Q2h*Sk*pow(f1,2)*m2 + 128*S*Sk*pow(f1,2)*m2 -
		       64*l1k2*Sq2*pow(f1,2)*m2 - 64*Q2e*Sq2*pow(f1,2)*m2 -
		       64*Q2h*Sq2*pow(f1,2)*m2 + 128*Sk*Sq2*pow(f1,2)*m2 -
		       32*l1k2*Q2e*pow(f2,2)*m2 + 32*Q2e*Q2h*pow(f2,2)*m2 -
		       64*f1*f2*Q2e*m4 - 64*f1*f2*Q2h*m4 -
		       64*Q2e*pow(f1,2)*m4 - 128*S*pow(f1,2)*m4 -
		       128*Sq2*pow(f1,2)*m4 - 32*Q2e*pow(f2,2)*m4 -
		       32*Q2h*pow(f2,2)*m4 - 64*Q2h*pow(f1,2)*m2*M2 +
		       64*f1*f2*m2*pow(Q2e,2) + 32*pow(f2,2)*m2*pow(Q2e,2) +
		       64*f1*f2*m2*pow(Q2h,2) + 32*pow(f1,2)*m2*pow(Q2h,2) +
		       16*pow(f2,2)*m2*pow(Q2h,2)))/2. +
		  (pow(l2k1,-2)/M2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (-8*l1k2*Q2e*Q2h*pow(f2,2)*m2 -
		       16*l1k2*Q2h*S*pow(f2,2)*m2 -
		       16*Q2e*Q2h*S*pow(f2,2)*m2 +
		       16*Q2e*Q2h*Sk*pow(f2,2)*m2 +
		       32*Q2h*S*Sk*pow(f2,2)*m2 -
		       16*l1k2*Q2h*Sq2*pow(f2,2)*m2 -
		       16*Q2e*Q2h*Sq2*pow(f2,2)*m2 +
		       32*Q2h*Sk*Sq2*pow(f2,2)*m2 - 8*Q2e*Q2h*pow(f2,2)*m4 -
		       32*Q2h*S*pow(f2,2)*m4 - 32*Q2h*Sq2*pow(f2,2)*m4 -
		       8*Q2h*pow(f2,2)*m2*pow(Q2e,2) -
		       16*S*pow(f2,2)*m2*pow(Q2h,2) +
		       8*Sk*pow(f2,2)*m2*pow(Q2h,2) -
		       16*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       8*pow(f2,2)*m4*pow(Q2h,2)))/2. +
		  ((128*f1*f2*Q2e - 64*f1*f2*Q2h + 192*S*pow(f1,2) + 64*Q2e*pow(f2,2) -
		       48*Q2h*pow(f2,2) + 384*pow(f1,2)*m2 - 64*pow(f1,2)*M2\
		)*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1))/2. +
		  (pow(M,-2)*(-16*Q2e*Q2h*pow(f2,2) + 48*Q2h*S*pow(f2,2) +
		       96*Q2h*pow(f2,2)*m2 + 8*pow(f2,2)*pow(Q2h,2))*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1))/2. +
		  (pow(l1k1,-2)*(-64*f1*f2*l2k2*Q2e*m2 -
		       64*f1*f2*Q2e*Q2h*m2 - 64*l2k2*Q2e*pow(f1,2)*m2 -
		       32*Q2e*Q2h*pow(f1,2)*m2 + 64*l2k2*S*pow(f1,2)*m2 -
		       64*Q2e*S*pow(f1,2)*m2 - 64*Q2h*S*pow(f1,2)*m2 -
		       64*Q2e*Sk*pow(f1,2)*m2 - 32*Q2h*Sk*pow(f1,2)*m2 +
		       128*S*Sk*pow(f1,2)*m2 - 32*l2k2*Q2e*pow(f2,2)*m2 -
		       32*Q2e*Q2h*pow(f2,2)*m2 + 64*f1*f2*Q2e*m4 +
		       64*f1*f2*Q2h*m4 + 64*Q2e*pow(f1,2)*m4 -
		       128*S*pow(f1,2)*m4 + 32*Q2e*pow(f2,2)*m4 +
		       32*Q2h*pow(f2,2)*m4 +
		       64*Q2h*pow(f1,2)*m2*M2 -
		       64*f1*f2*m2*pow(Q2e,2) - 32*pow(f2,2)*m2*pow(Q2e,2) -
		       64*f1*f2*m2*pow(Q2h,2) - 32*pow(f1,2)*m2*pow(Q2h,2) -
		       16*pow(f2,2)*m2*pow(Q2h,2))*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1))/2. +
		  (pow(l1k1,-2)/M2*(-8*l2k2*Q2e*Q2h*pow(f2,2)*m2 +
		       16*l2k2*Q2h*S*pow(f2,2)*m2 -
		       16*Q2e*Q2h*S*pow(f2,2)*m2 -
		       16*Q2e*Q2h*Sk*pow(f2,2)*m2 +
		       32*Q2h*S*Sk*pow(f2,2)*m2 + 8*Q2e*Q2h*pow(f2,2)*m4 -
		       32*Q2h*S*pow(f2,2)*m4 +
		       8*Q2h*pow(f2,2)*m2*pow(Q2e,2) -
		       16*S*pow(f2,2)*m2*pow(Q2h,2) -
		       8*Sk*pow(f2,2)*m2*pow(Q2h,2) -
		       8*pow(f2,2)*m4*pow(Q2h,2))*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)\
		)/2. + (pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (256*Q2e*Sq2*pow(f1,2) - 256*Q2h*Sq2*pow(f1,2) -
		       1024*S*Sq2*pow(f1,2) + 64*Q2e*Q2h*pow(f2,2) -
		       512*f1*f2*Q2e*m2 - 512*Q2e*pow(f1,2)*m2 -
		       256*Q2e*pow(f2,2)*m2 + 128*Q2h*pow(f2,2)*m2 +
		       256*Q2e*pow(f1,2)*M2 + 512*pow(f1,2)*m2*M2 -
		       512*f1*f2*pow(Q2e,2) - 256*pow(f1,2)*pow(Q2e,2) -
		       256*pow(f2,2)*pow(Q2e,2) - 256*f1*f2*pow(Q2h,2) -
		       128*pow(f2,2)*pow(Q2h,2) - 1024*pow(f1,2)*pow(S,2)))/2. +
		  (pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (96*f1*f2*Q2e*Q2h - 32*l2k1*Q2h*pow(f1,2) + 48*Q2e*Q2h*pow(f1,2) +
		       128*l2k1*S*pow(f1,2) + 64*l2k2*S*pow(f1,2) - 192*Q2e*S*pow(f1,2) +
		       32*Q2h*S*pow(f1,2) + 32*Q2h*Sk*pow(f1,2) + 128*l2k1*Sq2*pow(f1,2) +
		       64*l2k2*Sq2*pow(f1,2) - 192*Q2e*Sq2*pow(f1,2) +
		       64*Q2h*Sq2*pow(f1,2) - 128*S*Sq2*pow(f1,2) + 40*Q2e*Q2h*pow(f2,2) +
		       64*f1*f2*Q2e*m2 + 64*f1*f2*Q2h*m2 +
		       288*l2k1*pow(f1,2)*m2 + 96*l2k2*pow(f1,2)*m2 -
		       96*Q2e*pow(f1,2)*m2 + 128*Q2h*pow(f1,2)*m2 -
		       64*S*pow(f1,2)*m2 - 480*Sk*pow(f1,2)*m2 +
		       224*Sq2*pow(f1,2)*m2 + 32*Q2e*pow(f2,2)*m2 -
		       32*Q2e*pow(f1,2)*M2 + 64*Q2h*pow(f1,2)*M2 -
		       128*pow(f1,2)*m2*M2 - 32*f1*f2*pow(Q2e,2) -
		       16*pow(f2,2)*pow(Q2e,2) - 64*f1*f2*pow(Q2h,2) -
		       32*pow(f1,2)*pow(Q2h,2) - 16*pow(f2,2)*pow(Q2h,2) -
		       128*pow(f1,2)*pow(S,2)))/2. +
		  (pow(l2k1,-1)*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (96*f1*f2*Q2e*Q2h + 32*l1k1*Q2h*pow(f1,2) + 48*Q2e*Q2h*pow(f1,2) +
		       128*l1k1*S*pow(f1,2) + 64*l1k2*S*pow(f1,2) + 192*Q2e*S*pow(f1,2) -
		       32*Q2h*S*pow(f1,2) + 32*Q2h*Sk*pow(f1,2) + 32*Q2h*Sq2*pow(f1,2) -
		       128*S*Sq2*pow(f1,2) + 40*Q2e*Q2h*pow(f2,2) +
		       64*f1*f2*Q2e*m2 + 64*f1*f2*Q2h*m2 -
		       288*l1k1*pow(f1,2)*m2 - 96*l1k2*pow(f1,2)*m2 -
		       96*Q2e*pow(f1,2)*m2 + 128*Q2h*pow(f1,2)*m2 +
		       64*S*pow(f1,2)*m2 - 480*Sk*pow(f1,2)*m2 +
		       288*Sq2*pow(f1,2)*m2 + 32*Q2e*pow(f2,2)*m2 -
		       32*Q2e*pow(f1,2)*M2 + 64*Q2h*pow(f1,2)*M2 -
		       128*pow(f1,2)*m2*M2 - 32*f1*f2*pow(Q2e,2) -
		       16*pow(f2,2)*pow(Q2e,2) - 64*f1*f2*pow(Q2h,2) -
		       32*pow(f1,2)*pow(Q2h,2) - 16*pow(f2,2)*pow(Q2h,2) -
		       128*pow(f1,2)*pow(S,2)))/2. +
		  (pow(l1k1,-1)*pow(l2k2,-1)*(-64*f1*f2*l1k2*Q2e - 32*f1*f2*l2k1*Q2e -
		       32*f1*f2*l1k2*Q2h + 16*l1k2*Q2e*pow(f1,2) - 16*l1k2*Q2h*pow(f1,2) -
		       16*l2k1*Q2h*pow(f1,2) - 96*l1k2*S*pow(f1,2) + 128*l2k1*S*pow(f1,2) -
		       112*Q2e*S*pow(f1,2) + 32*Q2h*S*pow(f1,2) + 32*l2k1*Sk*pow(f1,2) +
		       16*Q2e*Sk*pow(f1,2) - 64*S*Sk*pow(f1,2) + 96*l2k1*Sq2*pow(f1,2) -
		       112*Q2e*Sq2*pow(f1,2) + 32*Q2h*Sq2*pow(f1,2) - 128*S*Sq2*pow(f1,2) -
		       64*Sk*Sq2*pow(f1,2) - 32*l1k2*Q2e*pow(f2,2) -
		       16*l2k1*Q2e*pow(f2,2) + 8*l2k1*Q2h*pow(f2,2) -
		       32*pow(f1,2)*pow(l2k1,2) + 32*f1*f2*Q2e*m2 +
		       128*f1*f2*Q2h*m2 + 16*l1k2*pow(f1,2)*m2 -
		       48*l2k1*pow(f1,2)*m2 + 56*Q2e*pow(f1,2)*m2 +
		       88*Q2h*pow(f1,2)*m2 - 192*S*pow(f1,2)*m2 +
		       48*Sk*pow(f1,2)*m2 - 80*Sq2*pow(f1,2)*m2 +
		       16*Q2e*pow(f2,2)*m2 + 64*Q2h*pow(f2,2)*m2 +
		       64*l1k2*pow(f1,2)*M2 + 32*l2k1*pow(f1,2)*M2 +
		       16*Q2h*pow(f1,2)*M2 + 24*pow(f1,2)*pow(Q2e,2) +
		       4*pow(f2,2)*pow(Q2h,2) - 128*pow(f1,2)*pow(S,2)))/2. +
		  (pow(l1k2,-1)*pow(l2k1,-1)*(32*f1*f2*l1k1*Q2e + 64*f1*f2*l2k2*Q2e +
		       32*f1*f2*l2k2*Q2h - 16*l2k2*Q2e*pow(f1,2) + 16*l1k1*Q2h*pow(f1,2) +
		       16*l2k2*Q2h*pow(f1,2) + 128*l1k1*S*pow(f1,2) - 96*l2k2*S*pow(f1,2) +
		       112*Q2e*S*pow(f1,2) - 32*Q2h*S*pow(f1,2) - 32*l1k1*Sk*pow(f1,2) +
		       16*Q2e*Sk*pow(f1,2) + 64*S*Sk*pow(f1,2) + 32*l1k1*Sq2*pow(f1,2) -
		       96*l2k2*Sq2*pow(f1,2) - 128*S*Sq2*pow(f1,2) +
		       16*l1k1*Q2e*pow(f2,2) + 32*l2k2*Q2e*pow(f2,2) -
		       8*l1k1*Q2h*pow(f2,2) - 32*pow(f1,2)*pow(l1k1,2) +
		       32*f1*f2*Q2e*m2 + 128*f1*f2*Q2h*m2 +
		       48*l1k1*pow(f1,2)*m2 - 16*l2k2*pow(f1,2)*m2 +
		       56*Q2e*pow(f1,2)*m2 + 88*Q2h*pow(f1,2)*m2 +
		       192*S*pow(f1,2)*m2 + 48*Sk*pow(f1,2)*m2 +
		       112*Sq2*pow(f1,2)*m2 + 16*Q2e*pow(f2,2)*m2 +
		       64*Q2h*pow(f2,2)*m2 - 32*l1k1*pow(f1,2)*M2 -
		       64*l2k2*pow(f1,2)*M2 + 16*Q2h*pow(f1,2)*M2 +
		       24*pow(f1,2)*pow(Q2e,2) + 4*pow(f2,2)*pow(Q2h,2) -
		       128*pow(f1,2)*pow(S,2)))/2. +
		  (pow(l2k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (-64*f1*f2*l1k2*Q2h + 64*f1*f2*l2k2*Q2h - 64*f1*f2*Q2e*Q2h -
		       32*l1k2*Q2h*pow(f1,2) + 32*l2k2*Q2h*pow(f1,2) -
		       32*Q2e*Q2h*pow(f1,2) + 64*l1k2*S*pow(f1,2) + 32*Q2h*Sk*pow(f1,2) -
		       128*S*Sk*pow(f1,2) + 64*l1k2*Sq2*pow(f1,2) + 128*S*Sq2*pow(f1,2) -
		       128*Sk*Sq2*pow(f1,2) - 16*l1k2*Q2h*pow(f2,2) +
		       32*l2k2*Q2h*pow(f2,2) - 24*Q2e*Q2h*pow(f2,2) -
		       128*f1*f2*Q2e*m2 - 96*l1k2*pow(f1,2)*m2 -
		       96*l2k2*pow(f1,2)*m2 - 112*Q2e*pow(f1,2)*m2 -
		       16*Q2h*pow(f1,2)*m2 + 288*Sk*pow(f1,2)*m2 -
		       288*Sq2*pow(f1,2)*m2 - 64*Q2e*pow(f2,2)*m2 +
		       64*l1k2*pow(f1,2)*M2 + 32*Q2e*pow(f1,2)*M2 +
		       32*f1*f2*pow(Q2e,2) + 32*pow(f1,2)*pow(Q2e,2) +
		       16*pow(f2,2)*pow(Q2e,2) - 32*f1*f2*pow(Q2h,2) -
		       16*pow(f1,2)*pow(Q2h,2) - 16*pow(f2,2)*pow(Q2h,2) +
		       128*pow(f1,2)*pow(S,2)))/2. +
		  (pow(l1k1,-1)*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-64*f1*f2*l1k2*Q2h + 64*f1*f2*l2k2*Q2h - 64*f1*f2*Q2e*Q2h -
		       32*l1k2*Q2h*pow(f1,2) + 32*l2k2*Q2h*pow(f1,2) -
		       32*Q2e*Q2h*pow(f1,2) + 64*l2k2*S*pow(f1,2) + 32*Q2h*Sk*pow(f1,2) +
		       128*S*Sk*pow(f1,2) + 128*S*Sq2*pow(f1,2) - 32*l1k2*Q2h*pow(f2,2) +
		       16*l2k2*Q2h*pow(f2,2) - 24*Q2e*Q2h*pow(f2,2) -
		       128*f1*f2*Q2e*m2 + 96*l1k2*pow(f1,2)*m2 +
		       96*l2k2*pow(f1,2)*m2 - 112*Q2e*pow(f1,2)*m2 -
		       16*Q2h*pow(f1,2)*m2 + 288*Sk*pow(f1,2)*m2 -
		       288*Sq2*pow(f1,2)*m2 - 64*Q2e*pow(f2,2)*m2 -
		       64*l2k2*pow(f1,2)*M2 + 32*Q2e*pow(f1,2)*M2 +
		       32*f1*f2*pow(Q2e,2) + 32*pow(f1,2)*pow(Q2e,2) +
		       16*pow(f2,2)*pow(Q2e,2) - 32*f1*f2*pow(Q2h,2) -
		       16*pow(f1,2)*pow(Q2h,2) - 16*pow(f2,2)*pow(Q2h,2) +
		       128*pow(f1,2)*pow(S,2)))/2. +
		  (pow(l2k2,-1)*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (32*f1*f2*Q2e*Q2h + 16*Q2e*Q2h*pow(f1,2) + 128*l1k1*S*pow(f1,2) +
		       128*l1k2*S*pow(f1,2) - 64*Q2e*S*pow(f1,2) + 192*Q2h*S*pow(f1,2) +
		       128*S*Sk*pow(f1,2) + 64*Q2h*Sq2*pow(f1,2) - 128*S*Sq2*pow(f1,2) +
		       32*Q2e*Q2h*pow(f2,2) - 64*f1*f2*Q2h*m2 -
		       128*Q2e*pow(f1,2)*m2 - 32*Q2h*pow(f1,2)*m2 +
		       128*S*pow(f1,2)*m2 - 32*Q2h*pow(f2,2)*m2 +
		       64*Q2e*pow(f1,2)*M2 - 64*Q2h*pow(f1,2)*M2 -
		       64*f1*f2*pow(Q2e,2) - 32*pow(f1,2)*pow(Q2e,2) -
		       32*pow(f2,2)*pow(Q2e,2) + 64*f1*f2*pow(Q2h,2) +
		       32*pow(f1,2)*pow(Q2h,2) + 16*pow(f2,2)*pow(Q2h,2) +
		       128*pow(f1,2)*pow(S,2)))/2. +
		  (pow(l2k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (64*f1*f2*l1k2*Q2e + 64*f1*f2*l2k1*Q2e + 64*f1*f2*l1k2*Q2h +
		       64*l1k2*Q2e*pow(f1,2) + 64*l2k1*Q2e*pow(f1,2) +
		       32*l1k2*Q2h*pow(f1,2) + 64*l2k1*Q2h*pow(f1,2) -
		       192*l2k1*S*pow(f1,2) - 32*Q2h*S*pow(f1,2) - 64*Q2e*Sk*pow(f1,2) -
		       64*Q2h*Sk*pow(f1,2) - 192*l2k1*Sq2*pow(f1,2) + 256*S*Sq2*pow(f1,2) +
		       32*l1k2*Q2e*pow(f2,2) + 32*l2k1*Q2e*pow(f2,2) +
		       16*l1k2*Q2h*pow(f2,2) - 16*l2k1*Q2h*pow(f2,2) -
		       64*f1*f2*Q2e*m2 - 256*f1*f2*Q2h*m2 -
		       96*l1k2*pow(f1,2)*m2 - 288*l2k1*pow(f1,2)*m2 +
		       112*Q2e*pow(f1,2)*m2 - 80*Q2h*pow(f1,2)*m2 +
		       256*S*pow(f1,2)*m2 + 288*Sk*pow(f1,2)*m2 -
		       224*Sq2*pow(f1,2)*m2 - 32*Q2e*pow(f2,2)*m2 -
		       128*Q2h*pow(f2,2)*m2 - 64*l1k2*pow(f1,2)*M2 -
		       64*l2k1*pow(f1,2)*M2 - 32*Q2h*pow(f1,2)*M2 +
		       32*pow(f1,2)*pow(Q2e,2) - 32*pow(f1,2)*pow(Q2h,2) -
		       8*pow(f2,2)*pow(Q2h,2) + 256*pow(f1,2)*pow(S,2)))/2. +
		  (pow(l1k2,-1)*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-64*f1*f2*l1k1*Q2e - 64*f1*f2*l2k2*Q2e - 64*f1*f2*l2k2*Q2h -
		       64*l1k1*Q2e*pow(f1,2) - 64*l2k2*Q2e*pow(f1,2) -
		       64*l1k1*Q2h*pow(f1,2) - 32*l2k2*Q2h*pow(f1,2) -
		       192*l1k1*S*pow(f1,2) + 32*Q2h*S*pow(f1,2) - 64*Q2e*Sk*pow(f1,2) -
		       64*Q2h*Sk*pow(f1,2) + 32*Q2h*Sq2*pow(f1,2) + 256*S*Sq2*pow(f1,2) -
		       32*l1k1*Q2e*pow(f2,2) - 32*l2k2*Q2e*pow(f2,2) +
		       16*l1k1*Q2h*pow(f2,2) - 16*l2k2*Q2h*pow(f2,2) -
		       64*f1*f2*Q2e*m2 - 256*f1*f2*Q2h*m2 +
		       288*l1k1*pow(f1,2)*m2 + 96*l2k2*pow(f1,2)*m2 +
		       112*Q2e*pow(f1,2)*m2 - 80*Q2h*pow(f1,2)*m2 -
		       256*S*pow(f1,2)*m2 + 288*Sk*pow(f1,2)*m2 -
		       480*Sq2*pow(f1,2)*m2 - 32*Q2e*pow(f2,2)*m2 -
		       128*Q2h*pow(f2,2)*m2 + 64*l1k1*pow(f1,2)*M2 +
		       64*l2k2*pow(f1,2)*M2 - 32*Q2h*pow(f1,2)*M2 +
		       32*pow(f1,2)*pow(Q2e,2) - 32*pow(f1,2)*pow(Q2h,2) -
		       8*pow(f2,2)*pow(Q2h,2) + 256*pow(f1,2)*pow(S,2)))/2. +
		  (pow(M,-2)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (64*Q2e*Q2h*Sq2*pow(f2,2) - 256*Q2h*S*Sq2*pow(f2,2) -
		       64*Q2e*Q2h*pow(f2,2)*m2 - 64*Sq2*pow(f2,2)*pow(Q2h,2) +
		       32*pow(f2,2)*pow(Q2h,3) - 256*Q2h*pow(f2,2)*pow(S,2)))/2. +
		  (pow(l1k2,-1)*pow(l2k1,-1)/M2*
		     (-4*l1k1*Q2e*Q2h*pow(f2,2) - 12*l2k2*Q2e*Q2h*pow(f2,2) +
		       32*l1k1*Q2h*S*pow(f2,2) - 24*l2k2*Q2h*S*pow(f2,2) +
		       28*Q2e*Q2h*S*pow(f2,2) - 8*l1k1*Q2h*Sk*pow(f2,2) +
		       4*Q2e*Q2h*Sk*pow(f2,2) + 16*Q2h*S*Sk*pow(f2,2) +
		       8*l1k1*Q2h*Sq2*pow(f2,2) - 24*l2k2*Q2h*Sq2*pow(f2,2) -
		       32*Q2h*S*Sq2*pow(f2,2) - 8*Q2h*pow(f2,2)*pow(l1k1,2) +
		       12*l1k1*Q2h*pow(f2,2)*m2 - 4*l2k2*Q2h*pow(f2,2)*m2 +
		       10*Q2e*Q2h*pow(f2,2)*m2 + 48*Q2h*S*pow(f2,2)*m2 +
		       12*Q2h*Sk*pow(f2,2)*m2 + 28*Q2h*Sq2*pow(f2,2)*m2 +
		       6*Q2h*pow(f2,2)*pow(Q2e,2) + 4*l1k1*pow(f2,2)*pow(Q2h,2) -
		       8*S*pow(f2,2)*pow(Q2h,2) + 6*pow(f2,2)*m2*pow(Q2h,2) -
		       32*Q2h*pow(f2,2)*pow(S,2)))/2. +
		  (pow(l1k1,-1)*pow(l2k2,-1)/M2*
		     (12*l1k2*Q2e*Q2h*pow(f2,2) + 4*l2k1*Q2e*Q2h*pow(f2,2) -
		       24*l1k2*Q2h*S*pow(f2,2) + 32*l2k1*Q2h*S*pow(f2,2) -
		       28*Q2e*Q2h*S*pow(f2,2) + 8*l2k1*Q2h*Sk*pow(f2,2) +
		       4*Q2e*Q2h*Sk*pow(f2,2) - 16*Q2h*S*Sk*pow(f2,2) +
		       24*l2k1*Q2h*Sq2*pow(f2,2) - 28*Q2e*Q2h*Sq2*pow(f2,2) -
		       32*Q2h*S*Sq2*pow(f2,2) - 16*Q2h*Sk*Sq2*pow(f2,2) -
		       8*Q2h*pow(f2,2)*pow(l2k1,2) + 4*l1k2*Q2h*pow(f2,2)*m2 -
		       12*l2k1*Q2h*pow(f2,2)*m2 + 10*Q2e*Q2h*pow(f2,2)*m2 -
		       48*Q2h*S*pow(f2,2)*m2 + 12*Q2h*Sk*pow(f2,2)*m2 -
		       20*Q2h*Sq2*pow(f2,2)*m2 + 6*Q2h*pow(f2,2)*pow(Q2e,2) -
		       4*l2k1*pow(f2,2)*pow(Q2h,2) + 8*S*pow(f2,2)*pow(Q2h,2) +
		       8*Sq2*pow(f2,2)*pow(Q2h,2) + 6*pow(f2,2)*m2*pow(Q2h,2) -
		       32*Q2h*pow(f2,2)*pow(S,2)))/2. +
		  (pow(l2k1,-1)/M2*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (32*l1k1*Q2h*S*pow(f2,2) + 16*l1k2*Q2h*S*pow(f2,2) +
		       48*Q2e*Q2h*S*pow(f2,2) - 32*Q2h*S*Sq2*pow(f2,2) -
		       72*l1k1*Q2h*pow(f2,2)*m2 - 24*l1k2*Q2h*pow(f2,2)*m2 -
		       32*Q2e*Q2h*pow(f2,2)*m2 + 16*Q2h*S*pow(f2,2)*m2 -
		       120*Q2h*Sk*pow(f2,2)*m2 + 72*Q2h*Sq2*pow(f2,2)*m2 +
		       4*Q2h*pow(f2,2)*pow(Q2e,2) + 8*l1k1*pow(f2,2)*pow(Q2h,2) -
		       8*S*pow(f2,2)*pow(Q2h,2) + 8*Sk*pow(f2,2)*pow(Q2h,2) +
		       8*Sq2*pow(f2,2)*pow(Q2h,2) + 24*pow(f2,2)*m2*pow(Q2h,2) -
		       32*Q2h*pow(f2,2)*pow(S,2)))/2. +
		  (pow(l1k1,-1)/M2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (32*l2k1*Q2h*S*pow(f2,2) + 16*l2k2*Q2h*S*pow(f2,2) -
		       48*Q2e*Q2h*S*pow(f2,2) + 32*l2k1*Q2h*Sq2*pow(f2,2) +
		       16*l2k2*Q2h*Sq2*pow(f2,2) - 48*Q2e*Q2h*Sq2*pow(f2,2) -
		       32*Q2h*S*Sq2*pow(f2,2) + 72*l2k1*Q2h*pow(f2,2)*m2 +
		       24*l2k2*Q2h*pow(f2,2)*m2 - 32*Q2e*Q2h*pow(f2,2)*m2 -
		       16*Q2h*S*pow(f2,2)*m2 - 120*Q2h*Sk*pow(f2,2)*m2 +
		       56*Q2h*Sq2*pow(f2,2)*m2 + 4*Q2h*pow(f2,2)*pow(Q2e,2) -
		       8*l2k1*pow(f2,2)*pow(Q2h,2) + 8*S*pow(f2,2)*pow(Q2h,2) +
		       8*Sk*pow(f2,2)*pow(Q2h,2) + 16*Sq2*pow(f2,2)*pow(Q2h,2) +
		       24*pow(f2,2)*m2*pow(Q2h,2) - 32*Q2h*pow(f2,2)*pow(S,2)))/2. +
		  (pow(l2k2,-1)/M2*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (32*l1k1*Q2h*S*pow(f2,2) + 32*l1k2*Q2h*S*pow(f2,2) -
		       16*Q2e*Q2h*S*pow(f2,2) + 32*Q2h*S*Sk*pow(f2,2) -
		       32*Q2h*S*Sq2*pow(f2,2) - 32*Q2e*Q2h*pow(f2,2)*m2 +
		       32*Q2h*S*pow(f2,2)*m2 + 48*S*pow(f2,2)*pow(Q2h,2) +
		       16*Sq2*pow(f2,2)*pow(Q2h,2) + 32*Q2h*pow(f2,2)*pow(S,2)))/2. +
		  (pow(l2k1,-1)/M2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (16*l1k2*Q2h*S*pow(f2,2) - 32*Q2h*S*Sk*pow(f2,2) +
		       16*l1k2*Q2h*Sq2*pow(f2,2) + 32*Q2h*S*Sq2*pow(f2,2) -
		       32*Q2h*Sk*Sq2*pow(f2,2) - 24*l1k2*Q2h*pow(f2,2)*m2 -
		       24*l2k2*Q2h*pow(f2,2)*m2 - 12*Q2e*Q2h*pow(f2,2)*m2 +
		       72*Q2h*Sk*pow(f2,2)*m2 - 72*Q2h*Sq2*pow(f2,2)*m2 +
		       4*Q2h*pow(f2,2)*pow(Q2e,2) + 8*Sk*pow(f2,2)*pow(Q2h,2) -
		       4*pow(f2,2)*m2*pow(Q2h,2) + 32*Q2h*pow(f2,2)*pow(S,2)))/2. +
		  (pow(l1k1,-1)/M2*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (16*l2k2*Q2h*S*pow(f2,2) + 32*Q2h*S*Sk*pow(f2,2) +
		       32*Q2h*S*Sq2*pow(f2,2) + 24*l1k2*Q2h*pow(f2,2)*m2 +
		       24*l2k2*Q2h*pow(f2,2)*m2 - 12*Q2e*Q2h*pow(f2,2)*m2 +
		       72*Q2h*Sk*pow(f2,2)*m2 - 72*Q2h*Sq2*pow(f2,2)*m2 +
		       4*Q2h*pow(f2,2)*pow(Q2e,2) + 8*Sk*pow(f2,2)*pow(Q2h,2) -
		       4*pow(f2,2)*m2*pow(Q2h,2) + 32*Q2h*pow(f2,2)*pow(S,2)))/2. +
		  (pow(l2k2,-1)/M2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (8*l1k2*Q2e*Q2h*pow(f2,2) + 8*l2k1*Q2e*Q2h*pow(f2,2) -
		       48*l2k1*Q2h*S*pow(f2,2) - 16*Q2e*Q2h*Sk*pow(f2,2) -
		       48*l2k1*Q2h*Sq2*pow(f2,2) + 64*Q2h*S*Sq2*pow(f2,2) -
		       24*l1k2*Q2h*pow(f2,2)*m2 - 72*l2k1*Q2h*pow(f2,2)*m2 +
		       36*Q2e*Q2h*pow(f2,2)*m2 + 64*Q2h*S*pow(f2,2)*m2 +
		       72*Q2h*Sk*pow(f2,2)*m2 - 56*Q2h*Sq2*pow(f2,2)*m2 +
		       8*Q2h*pow(f2,2)*pow(Q2e,2) + 16*l2k1*pow(f2,2)*pow(Q2h,2) -
		       8*S*pow(f2,2)*pow(Q2h,2) - 16*Sk*pow(f2,2)*pow(Q2h,2) +
		       12*pow(f2,2)*m2*pow(Q2h,2) - 8*pow(f2,2)*pow(Q2h,3) +
		       64*Q2h*pow(f2,2)*pow(S,2)))/2. +
		  (pow(l1k2,-1)/M2*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-8*l1k1*Q2e*Q2h*pow(f2,2) - 8*l2k2*Q2e*Q2h*pow(f2,2) -
		       48*l1k1*Q2h*S*pow(f2,2) - 16*Q2e*Q2h*Sk*pow(f2,2) +
		       64*Q2h*S*Sq2*pow(f2,2) + 72*l1k1*Q2h*pow(f2,2)*m2 +
		       24*l2k2*Q2h*pow(f2,2)*m2 + 36*Q2e*Q2h*pow(f2,2)*m2 -
		       64*Q2h*S*pow(f2,2)*m2 + 72*Q2h*Sk*pow(f2,2)*m2 -
		       120*Q2h*Sq2*pow(f2,2)*m2 + 8*Q2h*pow(f2,2)*pow(Q2e,2) -
		       16*l1k1*pow(f2,2)*pow(Q2h,2) + 8*S*pow(f2,2)*pow(Q2h,2) -
		       16*Sk*pow(f2,2)*pow(Q2h,2) + 8*Sq2*pow(f2,2)*pow(Q2h,2) +
		       12*pow(f2,2)*m2*pow(Q2h,2) - 8*pow(f2,2)*pow(Q2h,3) +
		       64*Q2h*pow(f2,2)*pow(S,2)))/2. +
		  (pow(l1k1,-1)*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-2)*
		     (-64*f1*f2*Q2e*Q2h*m2 + 128*Q2e*S*pow(f1,2)*m2 -
		       128*Q2h*S*pow(f1,2)*m2 - 32*Q2e*Q2h*pow(f2,2)*m2 +
		       128*f1*f2*Q2h*m4 + 64*Q2h*pow(f1,2)*m4 +
		       64*Q2h*pow(f2,2)*m4 + 64*Q2h*pow(f1,2)*m2*M2 -
		       64*f1*f2*m2*pow(Q2h,2) - 32*pow(f1,2)*m2*pow(Q2h,2) -
		       16*pow(f2,2)*m2*pow(Q2h,2) - 256*pow(f1,2)*m2*pow(S,2))\
		)/2. + (pow(l1k2,-1)*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-2)*
		     (-64*f1*f2*Q2e*Q2h*m2 + 128*Q2e*S*pow(f1,2)*m2 -
		       128*Q2h*S*pow(f1,2)*m2 - 32*Q2e*Q2h*pow(f2,2)*m2 +
		       128*f1*f2*Q2h*m4 + 64*Q2h*pow(f1,2)*m4 +
		       64*Q2h*pow(f2,2)*m4 + 64*Q2h*pow(f1,2)*m2*M2 -
		       64*f1*f2*m2*pow(Q2h,2) - 32*pow(f1,2)*m2*pow(Q2h,2) -
		       16*pow(f2,2)*m2*pow(Q2h,2) - 256*pow(f1,2)*m2*pow(S,2))\
		)/2. + (pow(l1k1,-1)*pow(l2k2,-1)*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (32*f1*f2*l1k2*Q2e*Q2h - 64*l1k2*Q2e*S*pow(f1,2) +
		       64*l1k2*Q2h*S*pow(f1,2) + 32*Q2e*Q2h*S*pow(f1,2) +
		       128*l1k2*S*Sk*pow(f1,2) + 128*Q2h*S*Sk*pow(f1,2) +
		       128*Q2h*S*Sq2*pow(f1,2) + 64*Q2h*Sk*Sq2*pow(f1,2) +
		       16*l1k2*Q2e*Q2h*pow(f2,2) + 64*f1*f2*l1k2*Q2e*m2 -
		       64*f1*f2*l1k2*Q2h*m2 + 64*f1*f2*Q2e*Q2h*m2 -
		       32*l1k2*Q2h*pow(f1,2)*m2 + 32*Q2e*Q2h*pow(f1,2)*m2 -
		       128*Q2h*S*pow(f1,2)*m2 - 64*Q2e*Sk*pow(f1,2)*m2 -
		       64*Q2h*Sk*pow(f1,2)*m2 + 128*S*Sk*pow(f1,2)*m2 +
		       64*Q2e*Sq2*pow(f1,2)*m2 - 64*Q2h*Sq2*pow(f1,2)*m2 +
		       32*l1k2*Q2e*pow(f2,2)*m2 - 48*l1k2*Q2h*pow(f2,2)*m2 +
		       24*Q2e*Q2h*pow(f2,2)*m2 - 64*f1*f2*Q2h*m4 +
		       64*Q2e*pow(f1,2)*m4 - 96*Q2h*pow(f1,2)*m4 -
		       256*S*pow(f1,2)*m4 - 32*Q2h*pow(f2,2)*m4 -
		       32*l1k2*Q2h*pow(f1,2)*M2 - 32*Q2e*Q2h*pow(f1,2)*M2 -
		       64*l1k2*pow(f1,2)*m2*M2 -
		       32*Q2e*pow(f1,2)*m2*M2 -
		       64*Q2h*pow(f1,2)*m2*M2 + 32*f1*f2*Q2h*pow(Q2e,2) -
		       32*S*pow(f1,2)*pow(Q2e,2) + 16*Q2h*pow(f2,2)*pow(Q2e,2) +
		       32*f1*f2*m2*pow(Q2e,2) + 16*pow(f2,2)*m2*pow(Q2e,2) +
		       16*f1*f2*Q2e*pow(Q2h,2) + 8*Q2e*pow(f1,2)*pow(Q2h,2) -
		       8*l1k2*pow(f2,2)*pow(Q2h,2) - 64*f1*f2*m2*pow(Q2h,2) -
		       64*pow(f1,2)*m2*pow(Q2h,2) -
		       48*pow(f2,2)*m2*pow(Q2h,2) - 16*f1*f2*pow(Q2h,3) -
		       8*pow(f1,2)*pow(Q2h,3) - 8*pow(f2,2)*pow(Q2h,3) +
		       256*l1k2*pow(f1,2)*pow(S,2) + 64*Q2e*pow(f1,2)*pow(S,2) +
		       192*Q2h*pow(f1,2)*pow(S,2) + 128*pow(f1,2)*m2*pow(S,2)))/2. +
		  (pow(l1k1,-1)/M2*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-2)*
		     (32*Q2e*Q2h*S*pow(f2,2)*m2 +
		       8*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       32*S*pow(f2,2)*m2*pow(Q2h,2) -
		       64*Q2h*pow(f2,2)*m2*pow(S,2)))/2. +
		  (pow(l1k2,-1)/M2*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-2)*
		     (32*Q2e*Q2h*S*pow(f2,2)*m2 +
		       8*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       32*S*pow(f2,2)*m2*pow(Q2h,2) -
		       64*Q2h*pow(f2,2)*m2*pow(S,2)))/2. +
		  (pow(l1k1,-1)*pow(l1k2,-1)*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-2)*
		     (128*f1*f2*Q2e*Q2h*m4 - 256*Q2e*S*pow(f1,2)*m4 +
		       256*Q2h*S*pow(f1,2)*m4 + 64*Q2e*Q2h*pow(f2,2)*m4 -
		       256*f1*f2*Q2h*m6 - 128*Q2h*pow(f1,2)*m6 -
		       128*Q2h*pow(f2,2)*m6 -
		       128*Q2h*pow(f1,2)*m4*M2 +
		       128*f1*f2*m4*pow(Q2h,2) + 64*pow(f1,2)*m4*pow(Q2h,2) +
		       32*pow(f2,2)*m4*pow(Q2h,2) + 512*pow(f1,2)*m4*pow(S,2))\
		)/2. + (pow(l1k1,-1)*pow(l1k2,-1)/M2*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-2)*
		     (-64*Q2e*Q2h*S*pow(f2,2)*m4 -
		       16*Q2e*pow(f2,2)*m4*pow(Q2h,2) +
		       64*S*pow(f2,2)*m4*pow(Q2h,2) +
		       128*Q2h*pow(f2,2)*m4*pow(S,2)))/2. +
		  (pow(l1k1,-1)*pow(l2k2,-1)/M2*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-16*l1k2*Q2e*Q2h*S*pow(f2,2) + 32*l1k2*Q2h*S*Sk*pow(f2,2) -
		       8*l1k2*Q2e*Q2h*pow(f2,2)*m2 -
		       16*Q2e*Q2h*Sk*pow(f2,2)*m2 +
		       32*Q2h*S*Sk*pow(f2,2)*m2 +
		       16*Q2e*Q2h*Sq2*pow(f2,2)*m2 +
		       16*Q2e*Q2h*pow(f2,2)*m4 - 64*Q2h*S*pow(f2,2)*m4 -
		       8*Q2h*S*pow(f2,2)*pow(Q2e,2) -
		       4*Q2h*pow(f2,2)*m2*pow(Q2e,2) -
		       4*l1k2*Q2e*pow(f2,2)*pow(Q2h,2) + 16*l1k2*S*pow(f2,2)*pow(Q2h,2) +
		       8*Q2e*S*pow(f2,2)*pow(Q2h,2) + 32*S*Sk*pow(f2,2)*pow(Q2h,2) +
		       32*S*Sq2*pow(f2,2)*pow(Q2h,2) + 16*Sk*Sq2*pow(f2,2)*pow(Q2h,2) -
		       32*S*pow(f2,2)*m2*pow(Q2h,2) -
		       16*Sk*pow(f2,2)*m2*pow(Q2h,2) -
		       16*Sq2*pow(f2,2)*m2*pow(Q2h,2) -
		       16*pow(f2,2)*m4*pow(Q2h,2) -
		       4*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) -
		       8*pow(f2,2)*m2*pow(Q2h,3) + 64*l1k2*Q2h*pow(f2,2)*pow(S,2) +
		       16*Q2e*Q2h*pow(f2,2)*pow(S,2) +
		       32*Q2h*pow(f2,2)*m2*pow(S,2) +
		       48*pow(f2,2)*pow(Q2h,2)*pow(S,2)))/2. +
		  (pow(l1k1,-1)*pow(l2k1,-1)*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (32*f1*f2*l1k2*Q2e*Q2h + 16*l1k2*Q2e*Q2h*pow(f1,2) +
		       32*l1k2*Q2e*S*pow(f1,2) - 16*Q2e*Q2h*S*pow(f1,2) +
		       32*Q2e*Q2h*Sk*pow(f1,2) + 128*Q2e*S*Sk*pow(f1,2) -
		       16*Q2e*Q2h*Sq2*pow(f1,2) - 64*Q2e*S*Sq2*pow(f1,2) +
		       8*l1k2*Q2e*Q2h*pow(f2,2) - 32*f1*f2*l1k2*Q2e*m2 +
		       64*f1*f2*l1k2*Q2h*m2 - 16*f1*f2*Q2e*Q2h*m2 -
		       48*l1k2*Q2e*pow(f1,2)*m2 + 32*l1k2*Q2h*pow(f1,2)*m2 +
		       8*Q2e*Q2h*pow(f1,2)*m2 - 64*l1k2*S*pow(f1,2)*m2 +
		       96*Q2e*S*pow(f1,2)*m2 - 32*Q2h*S*pow(f1,2)*m2 -
		       96*l1k2*Sk*pow(f1,2)*m2 - 112*Q2e*Sk*pow(f1,2)*m2 +
		       80*Q2h*Sk*pow(f1,2)*m2 + 112*Q2e*Sq2*pow(f1,2)*m2 -
		       48*Q2h*Sq2*pow(f1,2)*m2 + 128*S*Sq2*pow(f1,2)*m2 +
		       192*Sk*Sq2*pow(f1,2)*m2 - 16*l1k2*Q2e*pow(f2,2)*m2 +
		       16*l1k2*Q2h*pow(f2,2)*m2 - 16*Q2e*Q2h*pow(f2,2)*m2 +
		       64*f1*f2*Q2e*m4 - 64*f1*f2*Q2h*m4 +
		       128*S*pow(f1,2)*m4 + 32*Q2e*pow(f2,2)*m4 -
		       32*Q2h*pow(f2,2)*m4 - 32*l1k2*Q2e*pow(f1,2)*M2 +
		       16*Q2e*Q2h*pow(f1,2)*M2 -
		       64*l1k2*pow(f1,2)*m2*M2 -
		       32*Q2e*pow(f1,2)*m2*M2 +
		       32*Q2h*pow(f1,2)*m2*M2 + 32*f1*f2*l1k2*pow(Q2e,2) +
		       16*S*pow(f1,2)*pow(Q2e,2) - 32*Sk*pow(f1,2)*pow(Q2e,2) +
		       32*Sq2*pow(f1,2)*pow(Q2e,2) + 16*l1k2*pow(f2,2)*pow(Q2e,2) -
		       4*Q2h*pow(f2,2)*pow(Q2e,2) + 16*f1*f2*m2*pow(Q2e,2) -
		       8*pow(f1,2)*m2*pow(Q2e,2) +
		       8*pow(f2,2)*m2*pow(Q2e,2) -
		       16*pow(f1,2)*M2*pow(Q2e,2) - 32*f1*f2*pow(Q2e,3) -
		       16*pow(f1,2)*pow(Q2e,3) - 16*pow(f2,2)*pow(Q2e,3) -
		       16*f1*f2*Q2e*pow(Q2h,2) - 8*Q2e*pow(f1,2)*pow(Q2h,2) -
		       4*Q2e*pow(f2,2)*pow(Q2h,2) - 32*f1*f2*m2*pow(Q2h,2) -
		       16*pow(f1,2)*m2*pow(Q2h,2) -
		       8*pow(f2,2)*m2*pow(Q2h,2) - 192*pow(f1,2)*m2*pow(Sk,2))\
		)/2. + (pow(l1k2,-1)*pow(l2k1,-2)*
		     (-32*f1*f2*Q2e*Q2h*m2 - 16*l1k1*Q2h*pow(f1,2)*m2 -
		       8*Q2e*Q2h*pow(f1,2)*m2 + 32*l1k1*S*pow(f1,2)*m2 -
		       16*Q2e*S*pow(f1,2)*m2 + 80*Q2h*S*pow(f1,2)*m2 -
		       96*l1k1*Sk*pow(f1,2)*m2 - 16*Q2e*Sk*pow(f1,2)*m2 -
		       32*Q2h*Sk*pow(f1,2)*m2 + 64*S*Sk*pow(f1,2)*m2 +
		       32*l1k1*Sq2*pow(f1,2)*m2 - 16*Q2e*Sq2*pow(f1,2)*m2 +
		       80*Q2h*Sq2*pow(f1,2)*m2 + 64*Sk*Sq2*pow(f1,2)*m2 -
		       16*Q2e*Q2h*pow(f2,2)*m2 -
		       32*pow(f1,2)*pow(l1k1,2)*m2 - 32*f1*f2*Q2e*m4 +
		       32*f1*f2*Q2h*m4 - 64*l1k1*pow(f1,2)*m4 +
		       192*S*pow(f1,2)*m4 - 64*Sk*pow(f1,2)*m4 +
		       192*Sq2*pow(f1,2)*m4 - 16*Q2e*pow(f2,2)*m4 +
		       16*Q2h*pow(f2,2)*m4 + 32*Q2h*pow(f1,2)*m2*M2 +
		       8*pow(f1,2)*m2*pow(Q2e,2) - 32*f1*f2*m2*pow(Q2h,2) -
		       16*pow(f1,2)*m2*pow(Q2h,2) -
		       8*pow(f2,2)*m2*pow(Q2h,2) - 64*pow(f1,2)*m2*pow(Sk,2)))/
		   2. + (pow(l1k1,-1)*pow(l2k1,-1)*pow(l2k2,-1)*
		     (16*f1*f2*l1k2*Q2e*Q2h + 8*l1k2*Q2e*Q2h*pow(f1,2) -
		       16*l1k2*Q2e*S*pow(f1,2) - 16*l1k2*Q2e*Sk*pow(f1,2) -
		       8*Q2e*Q2h*Sq2*pow(f1,2) + 32*Q2e*S*Sq2*pow(f1,2) +
		       32*Q2e*Sk*Sq2*pow(f1,2) + 4*l1k2*Q2e*Q2h*pow(f2,2) +
		       16*f1*f2*l1k2*Q2e*m2 + 32*f1*f2*l1k2*Q2h*m2 +
		       16*f1*f2*Q2e*Q2h*m2 - 8*l1k2*Q2e*pow(f1,2)*m2 +
		       16*l1k2*Q2h*pow(f1,2)*m2 + 8*Q2e*Q2h*pow(f1,2)*m2 +
		       32*l1k2*S*pow(f1,2)*m2 + 64*Q2e*S*pow(f1,2)*m2 -
		       32*Q2h*S*pow(f1,2)*m2 - 16*l1k2*Sk*pow(f1,2)*m2 -
		       16*Q2e*Sk*pow(f1,2)*m2 + 128*S*Sk*pow(f1,2)*m2 +
		       56*Q2e*Sq2*pow(f1,2)*m2 - 24*Q2h*Sq2*pow(f1,2)*m2 -
		       64*S*Sq2*pow(f1,2)*m2 + 32*Sk*Sq2*pow(f1,2)*m2 +
		       8*l1k2*Q2e*pow(f2,2)*m2 + 8*l1k2*Q2h*pow(f2,2)*m2 +
		       8*Q2e*Q2h*pow(f2,2)*m2 + 32*f1*f2*Q2e*m4 +
		       32*f1*f2*Q2h*m4 + 32*Q2e*pow(f1,2)*m4 +
		       64*S*pow(f1,2)*m4 + 64*Sk*pow(f1,2)*m4 +
		       16*Q2e*pow(f2,2)*m4 + 16*Q2h*pow(f2,2)*m4 -
		       16*l1k2*Q2e*pow(f1,2)*M2 -
		       32*l1k2*pow(f1,2)*m2*M2 -
		       8*l1k2*pow(f1,2)*pow(Q2e,2) + 4*Q2h*pow(f1,2)*pow(Q2e,2) -
		       16*Sk*pow(f1,2)*pow(Q2e,2) + 16*Sq2*pow(f1,2)*pow(Q2e,2) -
		       32*f1*f2*m2*pow(Q2e,2) - 16*pow(f1,2)*m2*pow(Q2e,2) -
		       16*pow(f2,2)*m2*pow(Q2e,2) - 8*f1*f2*pow(Q2e,3) -
		       8*pow(f1,2)*pow(Q2e,3) - 4*pow(f2,2)*pow(Q2e,3) +
		       32*pow(f1,2)*m2*pow(Sk,2)))/2. +
		  (pow(l1k1,-2)*pow(l2k2,-1)*(32*f1*f2*Q2e*Q2h*m2 -
		       16*l2k1*Q2h*pow(f1,2)*m2 + 8*Q2e*Q2h*pow(f1,2)*m2 -
		       32*l2k1*S*pow(f1,2)*m2 - 16*Q2e*S*pow(f1,2)*m2 +
		       80*Q2h*S*pow(f1,2)*m2 - 96*l2k1*Sk*pow(f1,2)*m2 +
		       16*Q2e*Sk*pow(f1,2)*m2 + 32*Q2h*Sk*pow(f1,2)*m2 +
		       64*S*Sk*pow(f1,2)*m2 + 16*Q2e*Q2h*pow(f2,2)*m2 +
		       32*pow(f1,2)*pow(l2k1,2)*m2 + 32*f1*f2*Q2e*m4 -
		       32*f1*f2*Q2h*m4 - 64*l2k1*pow(f1,2)*m4 +
		       192*S*pow(f1,2)*m4 + 64*Sk*pow(f1,2)*m4 +
		       16*Q2e*pow(f2,2)*m4 - 16*Q2h*pow(f2,2)*m4 -
		       32*Q2h*pow(f1,2)*m2*M2 -
		       8*pow(f1,2)*m2*pow(Q2e,2) + 32*f1*f2*m2*pow(Q2h,2) +
		       16*pow(f1,2)*m2*pow(Q2h,2) +
		       8*pow(f2,2)*m2*pow(Q2h,2) + 64*pow(f1,2)*m2*pow(Sk,2)))/
		   2. + (pow(l1k1,-2)*pow(l2k2,-1)*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (64*Q2e*Q2h*S*pow(f1,2)*m2 +
		       32*Q2e*Q2h*Sk*pow(f1,2)*m2 -
		       128*Q2h*S*Sk*pow(f1,2)*m2 - 64*f1*f2*Q2e*Q2h*m4 +
		       128*Q2e*S*pow(f1,2)*m4 - 128*Q2h*S*pow(f1,2)*m4 +
		       64*Q2e*Sk*pow(f1,2)*m4 - 64*Q2h*Sk*pow(f1,2)*m4 -
		       256*S*Sk*pow(f1,2)*m4 - 32*Q2e*Q2h*pow(f2,2)*m4 +
		       128*f1*f2*Q2h*m6 + 64*Q2h*pow(f1,2)*m6 +
		       64*Q2h*pow(f2,2)*m6 + 64*Q2h*pow(f1,2)*m4*M2 -
		       32*f1*f2*Q2e*m2*pow(Q2h,2) -
		       64*S*pow(f1,2)*m2*pow(Q2h,2) -
		       32*Sk*pow(f1,2)*m2*pow(Q2h,2) -
		       16*Q2e*pow(f2,2)*m2*pow(Q2h,2) +
		       16*pow(f2,2)*m4*pow(Q2h,2) +
		       32*pow(f1,2)*m2*M2*pow(Q2h,2) -
		       32*f1*f2*m2*pow(Q2h,3) - 16*pow(f1,2)*m2*pow(Q2h,3) -
		       8*pow(f2,2)*m2*pow(Q2h,3) -
		       128*Q2h*pow(f1,2)*m2*pow(S,2) -
		       256*pow(f1,2)*m4*pow(S,2) -
		       64*Q2h*pow(f1,2)*m2*pow(Sk,2)))/2. +
		  (pow(l1k1,-1)*pow(l2k1,-1)/M2*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (8*l1k2*Q2e*Q2h*S*pow(f2,2) + 32*Q2e*Q2h*S*Sk*pow(f2,2) -
		       16*Q2e*Q2h*S*Sq2*pow(f2,2) - 8*l1k2*Q2e*Q2h*pow(f2,2)*m2 -
		       16*l1k2*Q2h*S*pow(f2,2)*m2 +
		       24*Q2e*Q2h*S*pow(f2,2)*m2 -
		       24*l1k2*Q2h*Sk*pow(f2,2)*m2 -
		       28*Q2e*Q2h*Sk*pow(f2,2)*m2 +
		       28*Q2e*Q2h*Sq2*pow(f2,2)*m2 +
		       32*Q2h*S*Sq2*pow(f2,2)*m2 +
		       48*Q2h*Sk*Sq2*pow(f2,2)*m2 - 8*Q2e*Q2h*pow(f2,2)*m4 +
		       32*Q2h*S*pow(f2,2)*m4 - 4*l1k2*Q2h*pow(f2,2)*pow(Q2e,2) +
		       4*Q2h*S*pow(f2,2)*pow(Q2e,2) - 8*Q2h*Sk*pow(f2,2)*pow(Q2e,2) +
		       8*Q2h*Sq2*pow(f2,2)*pow(Q2e,2) -
		       4*Q2h*pow(f2,2)*m2*pow(Q2e,2) -
		       4*Q2e*S*pow(f2,2)*pow(Q2h,2) + 8*Q2e*Sk*pow(f2,2)*pow(Q2h,2) -
		       4*Q2e*Sq2*pow(f2,2)*pow(Q2h,2) +
		       4*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       8*S*pow(f2,2)*m2*pow(Q2h,2) +
		       20*Sk*pow(f2,2)*m2*pow(Q2h,2) -
		       12*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       8*pow(f2,2)*m4*pow(Q2h,2) -
		       48*Q2h*pow(f2,2)*m2*pow(Sk,2)))/2. +
		  (pow(l1k2,-1)*pow(l2k1,-2)/M2*
		     (8*l1k1*Q2h*S*pow(f2,2)*m2 - 4*Q2e*Q2h*S*pow(f2,2)*m2 -
		       24*l1k1*Q2h*Sk*pow(f2,2)*m2 -
		       4*Q2e*Q2h*Sk*pow(f2,2)*m2 + 16*Q2h*S*Sk*pow(f2,2)*m2 +
		       8*l1k1*Q2h*Sq2*pow(f2,2)*m2 -
		       4*Q2e*Q2h*Sq2*pow(f2,2)*m2 +
		       16*Q2h*Sk*Sq2*pow(f2,2)*m2 -
		       8*Q2h*pow(f2,2)*pow(l1k1,2)*m2 -
		       16*l1k1*Q2h*pow(f2,2)*m4 + 4*Q2e*Q2h*pow(f2,2)*m4 +
		       48*Q2h*S*pow(f2,2)*m4 - 16*Q2h*Sk*pow(f2,2)*m4 +
		       48*Q2h*Sq2*pow(f2,2)*m4 +
		       2*Q2h*pow(f2,2)*m2*pow(Q2e,2) -
		       4*l1k1*pow(f2,2)*m2*pow(Q2h,2) +
		       2*Q2e*pow(f2,2)*m2*pow(Q2h,2) +
		       20*S*pow(f2,2)*m2*pow(Q2h,2) -
		       8*Sk*pow(f2,2)*m2*pow(Q2h,2) +
		       20*Sq2*pow(f2,2)*m2*pow(Q2h,2) -
		       4*pow(f2,2)*m4*pow(Q2h,2) -
		       16*Q2h*pow(f2,2)*m2*pow(Sk,2)))/2. +
		  (pow(l1k1,-1)*pow(l2k1,-1)*pow(l2k2,-1)/M2*
		     (-4*l1k2*Q2e*Q2h*S*pow(f2,2) - 4*l1k2*Q2e*Q2h*Sk*pow(f2,2) +
		       8*Q2e*Q2h*S*Sq2*pow(f2,2) + 8*Q2e*Q2h*Sk*Sq2*pow(f2,2) -
		       4*l1k2*Q2e*Q2h*pow(f2,2)*m2 +
		       8*l1k2*Q2h*S*pow(f2,2)*m2 +
		       16*Q2e*Q2h*S*pow(f2,2)*m2 -
		       4*l1k2*Q2h*Sk*pow(f2,2)*m2 -
		       4*Q2e*Q2h*Sk*pow(f2,2)*m2 + 32*Q2h*S*Sk*pow(f2,2)*m2 +
		       14*Q2e*Q2h*Sq2*pow(f2,2)*m2 -
		       16*Q2h*S*Sq2*pow(f2,2)*m2 +
		       8*Q2h*Sk*Sq2*pow(f2,2)*m2 + 4*Q2e*Q2h*pow(f2,2)*m4 +
		       16*Q2h*S*pow(f2,2)*m4 + 16*Q2h*Sk*pow(f2,2)*m4 -
		       2*l1k2*Q2h*pow(f2,2)*pow(Q2e,2) - 4*Q2h*Sk*pow(f2,2)*pow(Q2e,2) +
		       4*Q2h*Sq2*pow(f2,2)*pow(Q2e,2) - Q2h*pow(f2,2)*pow(Q2e,3) -
		       2*Q2e*Sq2*pow(f2,2)*pow(Q2h,2) -
		       8*S*pow(f2,2)*m2*pow(Q2h,2) -
		       6*Sq2*pow(f2,2)*m2*pow(Q2h,2) -
		       4*pow(f2,2)*m4*pow(Q2h,2) + pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) +
		       8*Q2h*pow(f2,2)*m2*pow(Sk,2)))/2. +
		  (pow(l1k1,-2)*pow(l2k2,-1)/M2*
		     (-8*l2k1*Q2h*S*pow(f2,2)*m2 - 4*Q2e*Q2h*S*pow(f2,2)*m2 -
		       24*l2k1*Q2h*Sk*pow(f2,2)*m2 +
		       4*Q2e*Q2h*Sk*pow(f2,2)*m2 + 16*Q2h*S*Sk*pow(f2,2)*m2 +
		       8*Q2h*pow(f2,2)*pow(l2k1,2)*m2 -
		       16*l2k1*Q2h*pow(f2,2)*m4 - 4*Q2e*Q2h*pow(f2,2)*m4 +
		       48*Q2h*S*pow(f2,2)*m4 + 16*Q2h*Sk*pow(f2,2)*m4 -
		       2*Q2h*pow(f2,2)*m2*pow(Q2e,2) -
		       4*l2k1*pow(f2,2)*m2*pow(Q2h,2) -
		       2*Q2e*pow(f2,2)*m2*pow(Q2h,2) +
		       20*S*pow(f2,2)*m2*pow(Q2h,2) +
		       8*Sk*pow(f2,2)*m2*pow(Q2h,2) +
		       4*pow(f2,2)*m4*pow(Q2h,2) +
		       16*Q2h*pow(f2,2)*m2*pow(Sk,2)))/2. +
		  (pow(l1k1,-2)*pow(l2k2,-1)/M2*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (32*Q2e*Q2h*S*pow(f2,2)*m4 +
		       16*Q2e*Q2h*Sk*pow(f2,2)*m4 -
		       64*Q2h*S*Sk*pow(f2,2)*m4 +
		       16*Q2e*S*pow(f2,2)*m2*pow(Q2h,2) +
		       8*Q2e*Sk*pow(f2,2)*m2*pow(Q2h,2) -
		       32*S*Sk*pow(f2,2)*m2*pow(Q2h,2) +
		       8*Q2e*pow(f2,2)*m4*pow(Q2h,2) -
		       32*S*pow(f2,2)*m4*pow(Q2h,2) -
		       16*Sk*pow(f2,2)*m4*pow(Q2h,2) +
		       4*Q2e*pow(f2,2)*m2*pow(Q2h,3) -
		       16*S*pow(f2,2)*m2*pow(Q2h,3) -
		       8*Sk*pow(f2,2)*m2*pow(Q2h,3) -
		       64*Q2h*pow(f2,2)*m4*pow(S,2) -
		       32*pow(f2,2)*m2*pow(Q2h,2)*pow(S,2) -
		       16*pow(f2,2)*m2*pow(Q2h,2)*pow(Sk,2)))/2. +
		  (pow(l1k1,-1)*pow(l1k2,-1)*(16*f1*f2*Q2e*Q2h + 16*l2k1*Q2e*pow(f1,2) +
		       16*l2k2*Q2e*pow(f1,2) + 8*Q2e*Q2h*pow(f1,2) - 32*l2k1*S*pow(f1,2) -
		       32*l2k2*S*pow(f1,2) + 48*Q2e*S*pow(f1,2) - 16*Q2e*Sk*pow(f1,2) -
		       32*l2k1*Sq2*pow(f1,2) - 32*l2k2*Sq2*pow(f1,2) +
		       80*Q2e*Sq2*pow(f1,2) - 64*S*Sq2*pow(f1,2) + 4*Q2e*Q2h*pow(f2,2) -
		       32*f1*f2*Q2e*m2 + 32*f1*f2*Q2h*m2 +
		       16*l2k2*pow(f1,2)*m2 - 40*Q2e*pow(f1,2)*m2 +
		       24*Q2h*pow(f1,2)*m2 - 64*S*pow(f1,2)*m2 +
		       32*Sk*pow(f1,2)*m2 - 32*Sq2*pow(f1,2)*m2 -
		       16*Q2e*pow(f2,2)*m2 + 16*Q2h*pow(f2,2)*m2 +
		       64*pow(f1,2)*m4 - 16*Q2e*pow(f1,2)*M2 -
		       24*pow(f1,2)*pow(Q2e,2) - 64*pow(f1,2)*pow(Sq2,2)))/2. +
		  (pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (32*f1*f2*Q2e*Q2h + 16*Q2e*Q2h*pow(f1,2) + 128*l2k1*S*pow(f1,2) +
		       128*l2k2*S*pow(f1,2) + 64*Q2e*S*pow(f1,2) - 192*Q2h*S*pow(f1,2) -
		       128*S*Sk*pow(f1,2) + 128*l2k1*Sq2*pow(f1,2) +
		       128*l2k2*Sq2*pow(f1,2) + 64*Q2e*Sq2*pow(f1,2) -
		       128*Q2h*Sq2*pow(f1,2) + 384*S*Sq2*pow(f1,2) - 128*Sk*Sq2*pow(f1,2) +
		       32*Q2e*Q2h*pow(f2,2) - 64*f1*f2*Q2h*m2 -
		       128*Q2e*pow(f1,2)*m2 - 32*Q2h*pow(f1,2)*m2 -
		       128*S*pow(f1,2)*m2 - 128*Sq2*pow(f1,2)*m2 -
		       32*Q2h*pow(f2,2)*m2 + 64*Q2e*pow(f1,2)*M2 -
		       64*Q2h*pow(f1,2)*M2 - 64*f1*f2*pow(Q2e,2) -
		       32*pow(f1,2)*pow(Q2e,2) - 32*pow(f2,2)*pow(Q2e,2) +
		       64*f1*f2*pow(Q2h,2) + 32*pow(f1,2)*pow(Q2h,2) +
		       16*pow(f2,2)*pow(Q2h,2) + 128*pow(f1,2)*pow(S,2) +
		       256*pow(f1,2)*pow(Sq2,2)))/2. +
		  (pow(l2k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-32*l1k2*Q2e*Q2h*pow(f1,2) + 128*Q2e*S*Sq2*pow(f1,2) -
		       384*Q2h*S*Sq2*pow(f1,2) - 16*l1k2*Q2e*Q2h*pow(f2,2) +
		       128*f1*f2*l1k2*Q2e*m2 + 128*f1*f2*l1k2*Q2h*m2 -
		       64*f1*f2*Q2e*Q2h*m2 + 128*l1k2*Q2e*pow(f1,2)*m2 +
		       64*l1k2*Q2h*pow(f1,2)*m2 - 96*Q2e*Q2h*pow(f1,2)*m2 -
		       128*Q2e*Sk*pow(f1,2)*m2 - 128*Q2h*Sk*pow(f1,2)*m2 +
		       128*Q2e*Sq2*pow(f1,2)*m2 - 256*S*Sq2*pow(f1,2)*m2 +
		       64*l1k2*Q2e*pow(f2,2)*m2 + 32*l1k2*Q2h*pow(f2,2)*m2 -
		       48*Q2e*Q2h*pow(f2,2)*m2 + 128*f1*f2*Q2h*m4 +
		       64*Q2h*pow(f1,2)*m4 + 64*Q2h*pow(f2,2)*m4 -
		       64*l1k2*Q2e*pow(f1,2)*M2 +
		       64*l1k2*Q2h*pow(f1,2)*M2 + 32*Q2e*Q2h*pow(f1,2)*M2 -
		       128*l1k2*pow(f1,2)*m2*M2 -
		       64*Q2e*pow(f1,2)*m2*M2 +
		       128*Q2h*pow(f1,2)*m2*M2 + 64*f1*f2*l1k2*pow(Q2e,2) -
		       32*f1*f2*Q2h*pow(Q2e,2) + 64*l1k2*pow(f1,2)*pow(Q2e,2) -
		       32*Q2h*pow(f1,2)*pow(Q2e,2) - 128*S*pow(f1,2)*pow(Q2e,2) -
		       64*Sk*pow(f1,2)*pow(Q2e,2) - 64*Sq2*pow(f1,2)*pow(Q2e,2) +
		       32*l1k2*pow(f2,2)*pow(Q2e,2) - 24*Q2h*pow(f2,2)*pow(Q2e,2) -
		       64*f1*f2*m2*pow(Q2e,2) + 64*pow(f1,2)*m2*pow(Q2e,2) -
		       32*pow(f2,2)*m2*pow(Q2e,2) -
		       32*pow(f1,2)*M2*pow(Q2e,2) + 32*f1*f2*pow(Q2e,3) +
		       32*pow(f1,2)*pow(Q2e,3) + 16*pow(f2,2)*pow(Q2e,3) -
		       64*f1*f2*l1k2*pow(Q2h,2) - 96*f1*f2*Q2e*pow(Q2h,2) -
		       32*l1k2*pow(f1,2)*pow(Q2h,2) - 48*Q2e*pow(f1,2)*pow(Q2h,2) +
		       128*S*pow(f1,2)*pow(Q2h,2) + 64*Sk*pow(f1,2)*pow(Q2h,2) +
		       64*Sq2*pow(f1,2)*pow(Q2h,2) - 16*l1k2*pow(f2,2)*pow(Q2h,2) -
		       40*Q2e*pow(f2,2)*pow(Q2h,2) + 32*pow(f2,2)*m2*pow(Q2h,2) +
		       64*pow(f1,2)*M2*pow(Q2h,2) - 32*f1*f2*pow(Q2h,3) +
		       16*pow(f1,2)*pow(Q2h,3) + 128*Q2e*pow(f1,2)*pow(S,2) -
		       384*Q2h*pow(f1,2)*pow(S,2) - 256*pow(f1,2)*m2*pow(S,2) -
		       128*Q2h*pow(f1,2)*pow(Sq2,2)))/2. +
		  (pow(l2k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (32*l1k2*Q2e*Q2h*pow(f1,2) + 128*Q2e*S*Sq2*pow(f1,2) -
		       384*Q2h*S*Sq2*pow(f1,2) + 16*l1k2*Q2e*Q2h*pow(f2,2) -
		       128*f1*f2*l1k2*Q2e*m2 - 128*f1*f2*l1k2*Q2h*m2 -
		       64*f1*f2*Q2e*Q2h*m2 - 128*l1k2*Q2e*pow(f1,2)*m2 -
		       64*l1k2*Q2h*pow(f1,2)*m2 - 64*Q2e*Q2h*pow(f1,2)*m2 +
		       128*Q2e*Sk*pow(f1,2)*m2 + 128*Q2h*Sk*pow(f1,2)*m2 -
		       128*Q2h*Sq2*pow(f1,2)*m2 - 256*S*Sq2*pow(f1,2)*m2 -
		       64*l1k2*Q2e*pow(f2,2)*m2 - 32*l1k2*Q2h*pow(f2,2)*m2 -
		       32*Q2e*Q2h*pow(f2,2)*m2 + 128*f1*f2*Q2h*m4 +
		       64*Q2h*pow(f1,2)*m4 + 64*Q2h*pow(f2,2)*m4 +
		       64*l1k2*Q2e*pow(f1,2)*M2 -
		       64*l1k2*Q2h*pow(f1,2)*M2 - 32*Q2e*Q2h*pow(f1,2)*M2 +
		       128*l1k2*pow(f1,2)*m2*M2 +
		       64*Q2h*pow(f1,2)*m2*M2 - 64*f1*f2*l1k2*pow(Q2e,2) -
		       64*f1*f2*Q2h*pow(Q2e,2) - 64*l1k2*pow(f1,2)*pow(Q2e,2) +
		       16*Q2h*pow(f1,2)*pow(Q2e,2) - 128*S*pow(f1,2)*pow(Q2e,2) +
		       64*Sk*pow(f1,2)*pow(Q2e,2) - 128*Sq2*pow(f1,2)*pow(Q2e,2) -
		       32*l1k2*pow(f2,2)*pow(Q2e,2) - 32*Q2h*pow(f2,2)*pow(Q2e,2) +
		       64*f1*f2*pow(Q2e,3) + 32*pow(f2,2)*pow(Q2e,3) +
		       64*f1*f2*l1k2*pow(Q2h,2) - 64*f1*f2*Q2e*pow(Q2h,2) +
		       32*l1k2*pow(f1,2)*pow(Q2h,2) - 16*Q2e*pow(f1,2)*pow(Q2h,2) +
		       128*S*pow(f1,2)*pow(Q2h,2) - 64*Sk*pow(f1,2)*pow(Q2h,2) +
		       128*Sq2*pow(f1,2)*pow(Q2h,2) + 16*l1k2*pow(f2,2)*pow(Q2h,2) -
		       40*Q2e*pow(f2,2)*pow(Q2h,2) + 64*f1*f2*m2*pow(Q2h,2) +
		       96*pow(f1,2)*m2*pow(Q2h,2) +
		       48*pow(f2,2)*m2*pow(Q2h,2) +
		       96*pow(f1,2)*M2*pow(Q2h,2) - 64*f1*f2*pow(Q2h,3) -
		       32*pow(f1,2)*pow(Q2h,3) - 8*pow(f2,2)*pow(Q2h,3) +
		       128*Q2e*pow(f1,2)*pow(S,2) - 384*Q2h*pow(f1,2)*pow(S,2) -
		       256*pow(f1,2)*m2*pow(S,2) - 128*Q2h*pow(f1,2)*pow(Sq2,2)))/2. +
		  (pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-32*l2k2*Q2e*Q2h*pow(f1,2) - 128*Q2e*S*Sq2*pow(f1,2) +
		       384*Q2h*S*Sq2*pow(f1,2) - 16*l2k2*Q2e*Q2h*pow(f2,2) +
		       128*f1*f2*l2k2*Q2e*m2 + 128*f1*f2*l2k2*Q2h*m2 +
		       64*f1*f2*Q2e*Q2h*m2 + 128*l2k2*Q2e*pow(f1,2)*m2 +
		       64*l2k2*Q2h*pow(f1,2)*m2 + 96*Q2e*Q2h*pow(f1,2)*m2 +
		       128*Q2e*Sk*pow(f1,2)*m2 + 128*Q2h*Sk*pow(f1,2)*m2 -
		       128*Q2e*Sq2*pow(f1,2)*m2 + 256*S*Sq2*pow(f1,2)*m2 +
		       64*l2k2*Q2e*pow(f2,2)*m2 + 32*l2k2*Q2h*pow(f2,2)*m2 +
		       48*Q2e*Q2h*pow(f2,2)*m2 - 128*f1*f2*Q2h*m4 -
		       64*Q2h*pow(f1,2)*m4 - 64*Q2h*pow(f2,2)*m4 -
		       64*l2k2*Q2e*pow(f1,2)*M2 +
		       64*l2k2*Q2h*pow(f1,2)*M2 - 32*Q2e*Q2h*pow(f1,2)*M2 -
		       128*l2k2*pow(f1,2)*m2*M2 +
		       64*Q2e*pow(f1,2)*m2*M2 -
		       128*Q2h*pow(f1,2)*m2*M2 + 64*f1*f2*l2k2*pow(Q2e,2) +
		       32*f1*f2*Q2h*pow(Q2e,2) + 64*l2k2*pow(f1,2)*pow(Q2e,2) +
		       32*Q2h*pow(f1,2)*pow(Q2e,2) - 128*S*pow(f1,2)*pow(Q2e,2) +
		       64*Sk*pow(f1,2)*pow(Q2e,2) - 64*Sq2*pow(f1,2)*pow(Q2e,2) +
		       32*l2k2*pow(f2,2)*pow(Q2e,2) + 24*Q2h*pow(f2,2)*pow(Q2e,2) +
		       64*f1*f2*m2*pow(Q2e,2) - 64*pow(f1,2)*m2*pow(Q2e,2) +
		       32*pow(f2,2)*m2*pow(Q2e,2) +
		       32*pow(f1,2)*M2*pow(Q2e,2) - 32*f1*f2*pow(Q2e,3) -
		       32*pow(f1,2)*pow(Q2e,3) - 16*pow(f2,2)*pow(Q2e,3) -
		       64*f1*f2*l2k2*pow(Q2h,2) + 96*f1*f2*Q2e*pow(Q2h,2) -
		       32*l2k2*pow(f1,2)*pow(Q2h,2) + 48*Q2e*pow(f1,2)*pow(Q2h,2) +
		       128*S*pow(f1,2)*pow(Q2h,2) - 64*Sk*pow(f1,2)*pow(Q2h,2) +
		       64*Sq2*pow(f1,2)*pow(Q2h,2) - 16*l2k2*pow(f2,2)*pow(Q2h,2) +
		       40*Q2e*pow(f2,2)*pow(Q2h,2) - 32*pow(f2,2)*m2*pow(Q2h,2) -
		       64*pow(f1,2)*M2*pow(Q2h,2) + 32*f1*f2*pow(Q2h,3) -
		       16*pow(f1,2)*pow(Q2h,3) - 128*Q2e*pow(f1,2)*pow(S,2) +
		       384*Q2h*pow(f1,2)*pow(S,2) + 256*pow(f1,2)*m2*pow(S,2) +
		       128*Q2h*pow(f1,2)*pow(Sq2,2)))/2. +
		  (pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (32*l2k2*Q2e*Q2h*pow(f1,2) - 128*Q2e*S*Sq2*pow(f1,2) +
		       384*Q2h*S*Sq2*pow(f1,2) + 16*l2k2*Q2e*Q2h*pow(f2,2) -
		       128*f1*f2*l2k2*Q2e*m2 - 128*f1*f2*l2k2*Q2h*m2 +
		       64*f1*f2*Q2e*Q2h*m2 - 128*l2k2*Q2e*pow(f1,2)*m2 -
		       64*l2k2*Q2h*pow(f1,2)*m2 + 64*Q2e*Q2h*pow(f1,2)*m2 -
		       128*Q2e*Sk*pow(f1,2)*m2 - 128*Q2h*Sk*pow(f1,2)*m2 +
		       128*Q2h*Sq2*pow(f1,2)*m2 + 256*S*Sq2*pow(f1,2)*m2 -
		       64*l2k2*Q2e*pow(f2,2)*m2 - 32*l2k2*Q2h*pow(f2,2)*m2 +
		       32*Q2e*Q2h*pow(f2,2)*m2 - 128*f1*f2*Q2h*m4 -
		       64*Q2h*pow(f1,2)*m4 - 64*Q2h*pow(f2,2)*m4 +
		       64*l2k2*Q2e*pow(f1,2)*M2 -
		       64*l2k2*Q2h*pow(f1,2)*M2 + 32*Q2e*Q2h*pow(f1,2)*M2 +
		       128*l2k2*pow(f1,2)*m2*M2 -
		       64*Q2h*pow(f1,2)*m2*M2 - 64*f1*f2*l2k2*pow(Q2e,2) +
		       64*f1*f2*Q2h*pow(Q2e,2) - 64*l2k2*pow(f1,2)*pow(Q2e,2) -
		       16*Q2h*pow(f1,2)*pow(Q2e,2) - 128*S*pow(f1,2)*pow(Q2e,2) -
		       64*Sk*pow(f1,2)*pow(Q2e,2) - 32*l2k2*pow(f2,2)*pow(Q2e,2) +
		       32*Q2h*pow(f2,2)*pow(Q2e,2) - 64*f1*f2*pow(Q2e,3) -
		       32*pow(f2,2)*pow(Q2e,3) + 64*f1*f2*l2k2*pow(Q2h,2) +
		       64*f1*f2*Q2e*pow(Q2h,2) + 32*l2k2*pow(f1,2)*pow(Q2h,2) +
		       16*Q2e*pow(f1,2)*pow(Q2h,2) + 128*S*pow(f1,2)*pow(Q2h,2) +
		       64*Sk*pow(f1,2)*pow(Q2h,2) + 16*l2k2*pow(f2,2)*pow(Q2h,2) +
		       40*Q2e*pow(f2,2)*pow(Q2h,2) - 64*f1*f2*m2*pow(Q2h,2) -
		       96*pow(f1,2)*m2*pow(Q2h,2) -
		       48*pow(f2,2)*m2*pow(Q2h,2) -
		       96*pow(f1,2)*M2*pow(Q2h,2) + 64*f1*f2*pow(Q2h,3) +
		       32*pow(f1,2)*pow(Q2h,3) + 8*pow(f2,2)*pow(Q2h,3) -
		       128*Q2e*pow(f1,2)*pow(S,2) + 384*Q2h*pow(f1,2)*pow(S,2) +
		       256*pow(f1,2)*m2*pow(S,2) + 128*Q2h*pow(f1,2)*pow(Sq2,2)))/2. +
		  (pow(l1k1,-1)*pow(l1k2,-1)/M2*
		     (4*l2k1*Q2e*Q2h*pow(f2,2) + 4*l2k2*Q2e*Q2h*pow(f2,2) -
		       8*l2k1*Q2h*S*pow(f2,2) - 8*l2k2*Q2h*S*pow(f2,2) +
		       12*Q2e*Q2h*S*pow(f2,2) - 4*Q2e*Q2h*Sk*pow(f2,2) -
		       8*l2k1*Q2h*Sq2*pow(f2,2) - 8*l2k2*Q2h*Sq2*pow(f2,2) +
		       20*Q2e*Q2h*Sq2*pow(f2,2) - 16*Q2h*S*Sq2*pow(f2,2) +
		       4*l2k2*Q2h*pow(f2,2)*m2 - 6*Q2e*Q2h*pow(f2,2)*m2 -
		       16*Q2h*S*pow(f2,2)*m2 + 8*Q2h*Sk*pow(f2,2)*m2 -
		       8*Q2h*Sq2*pow(f2,2)*m2 + 16*Q2h*pow(f2,2)*m4 -
		       6*Q2h*pow(f2,2)*pow(Q2e,2) + 2*pow(f2,2)*m2*pow(Q2h,2) -
		       16*Q2h*pow(f2,2)*pow(Sq2,2)))/2. +
		  (pow(l1k2,-1)/M2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (32*l2k1*Q2h*S*pow(f2,2) + 32*l2k2*Q2h*S*pow(f2,2) +
		       16*Q2e*Q2h*S*pow(f2,2) - 32*Q2h*S*Sk*pow(f2,2) +
		       32*l2k1*Q2h*Sq2*pow(f2,2) + 32*l2k2*Q2h*Sq2*pow(f2,2) +
		       16*Q2e*Q2h*Sq2*pow(f2,2) + 96*Q2h*S*Sq2*pow(f2,2) -
		       32*Q2h*Sk*Sq2*pow(f2,2) - 32*Q2e*Q2h*pow(f2,2)*m2 -
		       32*Q2h*S*pow(f2,2)*m2 - 32*Q2h*Sq2*pow(f2,2)*m2 -
		       48*S*pow(f2,2)*pow(Q2h,2) - 32*Sq2*pow(f2,2)*pow(Q2h,2) +
		       32*Q2h*pow(f2,2)*pow(S,2) + 64*Q2h*pow(f2,2)*pow(Sq2,2)))/2. +
		  (pow(l1k2,-1)*pow(l2k1,-1)*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-32*f1*f2*l1k1*Q2e*Q2h - 16*l1k1*Q2e*Q2h*pow(f1,2) +
		       128*l1k1*Q2e*S*pow(f1,2) - 224*l1k1*Q2h*S*pow(f1,2) +
		       16*Q2e*Q2h*S*pow(f1,2) + 128*l1k1*Q2h*Sk*pow(f1,2) -
		       48*Q2e*Q2h*Sk*pow(f1,2) + 128*l1k1*S*Sk*pow(f1,2) +
		       64*Q2e*S*Sk*pow(f1,2) - 128*Q2h*S*Sk*pow(f1,2) -
		       64*l1k1*Q2h*Sq2*pow(f1,2) + 80*Q2e*Q2h*Sq2*pow(f1,2) -
		       128*l1k1*S*Sq2*pow(f1,2) - 64*Q2e*S*Sq2*pow(f1,2) +
		       192*Q2h*S*Sq2*pow(f1,2) - 64*Q2h*Sk*Sq2*pow(f1,2) -
		       16*l1k1*Q2e*Q2h*pow(f2,2) + 64*Q2h*pow(f1,2)*pow(l1k1,2) +
		       128*S*pow(f1,2)*pow(l1k1,2) + 64*f1*f2*l1k1*Q2e*m2 +
		       32*f1*f2*Q2e*Q2h*m2 - 128*l1k1*Q2e*pow(f1,2)*m2 +
		       32*l1k1*Q2h*pow(f1,2)*m2 + 64*Q2e*Q2h*pow(f1,2)*m2 +
		       64*l1k1*S*pow(f1,2)*m2 + 224*Q2e*S*pow(f1,2)*m2 -
		       160*Q2h*S*pow(f1,2)*m2 - 384*l1k1*Sk*pow(f1,2)*m2 -
		       224*Q2e*Sk*pow(f1,2)*m2 - 32*Q2h*Sk*pow(f1,2)*m2 +
		       384*l1k1*Sq2*pow(f1,2)*m2 + 288*Q2e*Sq2*pow(f1,2)*m2 -
		       128*Q2h*Sq2*pow(f1,2)*m2 - 128*S*Sq2*pow(f1,2)*m2 +
		       384*Sk*Sq2*pow(f1,2)*m2 + 32*l1k1*Q2e*pow(f2,2)*m2 +
		       16*Q2e*Q2h*pow(f2,2)*m2 -
		       192*pow(f1,2)*pow(l1k1,2)*m2 - 64*f1*f2*Q2e*m4 +
		       128*f1*f2*Q2h*m4 - 64*Q2e*pow(f1,2)*m4 +
		       96*Q2h*pow(f1,2)*m4 + 128*S*pow(f1,2)*m4 -
		       32*Q2e*pow(f2,2)*m4 + 64*Q2h*pow(f2,2)*m4 +
		       32*l1k1*Q2h*pow(f1,2)*M2 + 16*Q2e*Q2h*pow(f1,2)*M2 -
		       32*Q2h*pow(f1,2)*m2*M2 - 16*f1*f2*Q2h*pow(Q2e,2) -
		       8*Q2h*pow(f1,2)*pow(Q2e,2) + 32*S*pow(f1,2)*pow(Q2e,2) -
		       8*Q2h*pow(f2,2)*pow(Q2e,2) + 96*f1*f2*m2*pow(Q2e,2) +
		       48*pow(f2,2)*m2*pow(Q2e,2) + 16*f1*f2*Q2e*pow(Q2h,2) -
		       8*Q2e*pow(f1,2)*pow(Q2h,2) - 48*S*pow(f1,2)*pow(Q2h,2) +
		       32*Sk*pow(f1,2)*pow(Q2h,2) - 80*Sq2*pow(f1,2)*pow(Q2h,2) +
		       8*l1k1*pow(f2,2)*pow(Q2h,2) + 12*Q2e*pow(f2,2)*pow(Q2h,2) -
		       128*f1*f2*m2*pow(Q2h,2) - 80*pow(f1,2)*m2*pow(Q2h,2) -
		       72*pow(f2,2)*m2*pow(Q2h,2) -
		       48*pow(f1,2)*M2*pow(Q2h,2) + 32*f1*f2*pow(Q2h,3) +
		       16*pow(f1,2)*pow(Q2h,3) + 4*pow(f2,2)*pow(Q2h,3) -
		       128*l1k1*pow(f1,2)*pow(S,2) - 64*Q2e*pow(f1,2)*pow(S,2) +
		       320*Q2h*pow(f1,2)*pow(S,2) - 128*pow(f1,2)*m2*pow(S,2) +
		       64*Q2h*pow(f1,2)*pow(Sk,2) - 192*pow(f1,2)*m2*pow(Sk,2) -
		       192*pow(f1,2)*m2*pow(Sq2,2)))/2. +
		  (pow(l1k2,-1)*pow(l2k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (32*f1*f2*l2k2*Q2e*Q2h + 64*l2k2*Q2e*S*pow(f1,2) -
		       64*l2k2*Q2h*S*pow(f1,2) + 32*Q2e*Q2h*S*pow(f1,2) -
		       128*l2k2*S*Sk*pow(f1,2) + 128*Q2h*S*Sk*pow(f1,2) +
		       64*l2k2*Q2e*Sq2*pow(f1,2) - 64*l2k2*Q2h*Sq2*pow(f1,2) +
		       32*Q2e*Q2h*Sq2*pow(f1,2) + 512*l2k2*S*Sq2*pow(f1,2) -
		       128*Q2e*S*Sq2*pow(f1,2) - 256*Q2h*S*Sq2*pow(f1,2) -
		       128*l2k2*Sk*Sq2*pow(f1,2) + 64*Q2h*Sk*Sq2*pow(f1,2) +
		       16*l2k2*Q2e*Q2h*pow(f2,2) + 64*f1*f2*l2k2*Q2e*m2 -
		       64*f1*f2*l2k2*Q2h*m2 - 64*f1*f2*Q2e*Q2h*m2 -
		       32*l2k2*Q2h*pow(f1,2)*m2 - 32*Q2e*Q2h*pow(f1,2)*m2 -
		       128*Q2h*S*pow(f1,2)*m2 + 64*Q2e*Sk*pow(f1,2)*m2 +
		       64*Q2h*Sk*pow(f1,2)*m2 + 128*S*Sk*pow(f1,2)*m2 -
		       64*Q2e*Sq2*pow(f1,2)*m2 - 64*Q2h*Sq2*pow(f1,2)*m2 -
		       256*S*Sq2*pow(f1,2)*m2 + 128*Sk*Sq2*pow(f1,2)*m2 +
		       32*l2k2*Q2e*pow(f2,2)*m2 - 48*l2k2*Q2h*pow(f2,2)*m2 -
		       24*Q2e*Q2h*pow(f2,2)*m2 + 64*f1*f2*Q2h*m4 -
		       64*Q2e*pow(f1,2)*m4 + 96*Q2h*pow(f1,2)*m4 -
		       256*S*pow(f1,2)*m4 - 256*Sq2*pow(f1,2)*m4 +
		       32*Q2h*pow(f2,2)*m4 - 32*l2k2*Q2h*pow(f1,2)*M2 +
		       32*Q2e*Q2h*pow(f1,2)*M2 -
		       64*l2k2*pow(f1,2)*m2*M2 +
		       32*Q2e*pow(f1,2)*m2*M2 +
		       64*Q2h*pow(f1,2)*m2*M2 - 32*f1*f2*Q2h*pow(Q2e,2) -
		       32*S*pow(f1,2)*pow(Q2e,2) - 32*Sq2*pow(f1,2)*pow(Q2e,2) -
		       16*Q2h*pow(f2,2)*pow(Q2e,2) - 32*f1*f2*m2*pow(Q2e,2) -
		       16*pow(f2,2)*m2*pow(Q2e,2) - 16*f1*f2*Q2e*pow(Q2h,2) -
		       8*Q2e*pow(f1,2)*pow(Q2h,2) - 8*l2k2*pow(f2,2)*pow(Q2h,2) +
		       64*f1*f2*m2*pow(Q2h,2) + 64*pow(f1,2)*m2*pow(Q2h,2) +
		       48*pow(f2,2)*m2*pow(Q2h,2) + 16*f1*f2*pow(Q2h,3) +
		       8*pow(f1,2)*pow(Q2h,3) + 8*pow(f2,2)*pow(Q2h,3) +
		       256*l2k2*pow(f1,2)*pow(S,2) - 64*Q2e*pow(f1,2)*pow(S,2) -
		       192*Q2h*pow(f1,2)*pow(S,2) - 128*pow(f1,2)*m2*pow(S,2) +
		       256*l2k2*pow(f1,2)*pow(Sq2,2) - 64*Q2e*pow(f1,2)*pow(Sq2,2) -
		       64*Q2h*pow(f1,2)*pow(Sq2,2) - 128*pow(f1,2)*m2*pow(Sq2,2)))/2. \
		+ (pow(l2k1,-1)*pow(l2k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (-32*f1*f2*l1k2*Q2e*m2 - 64*f1*f2*l1k2*Q2h*m2 +
		       16*l1k2*Q2e*pow(f1,2)*m2 - 32*l1k2*Q2h*pow(f1,2)*m2 +
		       16*Q2e*Q2h*pow(f1,2)*m2 + 128*l1k2*S*pow(f1,2)*m2 -
		       160*Q2e*S*pow(f1,2)*m2 + 32*Q2h*S*pow(f1,2)*m2 -
		       112*Q2e*Sk*pow(f1,2)*m2 + 48*Q2h*Sk*pow(f1,2)*m2 -
		       256*S*Sk*pow(f1,2)*m2 + 32*l1k2*Sq2*pow(f1,2)*m2 -
		       32*Q2e*Sq2*pow(f1,2)*m2 - 128*S*Sq2*pow(f1,2)*m2 -
		       64*Sk*Sq2*pow(f1,2)*m2 - 16*l1k2*Q2e*pow(f2,2)*m2 -
		       16*l1k2*Q2h*pow(f2,2)*m2 + 64*f1*f2*Q2e*m4 -
		       64*f1*f2*Q2h*m4 - 64*Q2e*pow(f1,2)*m4 +
		       64*Q2h*pow(f1,2)*m4 - 384*S*pow(f1,2)*m4 -
		       384*Sq2*pow(f1,2)*m4 + 32*Q2e*pow(f2,2)*m4 -
		       32*Q2h*pow(f2,2)*m4 +
		       64*l1k2*pow(f1,2)*m2*M2 +
		       64*f1*f2*m2*pow(Q2e,2) + 16*pow(f1,2)*m2*pow(Q2e,2) +
		       32*pow(f2,2)*m2*pow(Q2e,2) -
		       256*pow(f1,2)*m2*pow(S,2) - 64*pow(f1,2)*m2*pow(Sq2,2))\
		)/2. + (pow(l1k1,-1)*pow(l1k2,-1)*pow(l2k1,-1)*
		     (16*f1*f2*l2k2*Q2e*Q2h + 8*l2k2*Q2e*Q2h*pow(f1,2) +
		       16*l2k2*Q2e*S*pow(f1,2) - 16*l2k2*Q2e*Sk*pow(f1,2) +
		       16*l2k2*Q2e*Sq2*pow(f1,2) + 8*Q2e*Q2h*Sq2*pow(f1,2) +
		       32*Q2e*S*Sq2*pow(f1,2) - 32*Q2e*Sk*Sq2*pow(f1,2) +
		       4*l2k2*Q2e*Q2h*pow(f2,2) + 16*f1*f2*l2k2*Q2e*m2 +
		       32*f1*f2*l2k2*Q2h*m2 - 16*f1*f2*Q2e*Q2h*m2 -
		       8*l2k2*Q2e*pow(f1,2)*m2 + 16*l2k2*Q2h*pow(f1,2)*m2 -
		       8*Q2e*Q2h*pow(f1,2)*m2 - 32*l2k2*S*pow(f1,2)*m2 +
		       64*Q2e*S*pow(f1,2)*m2 - 32*Q2h*S*pow(f1,2)*m2 -
		       16*l2k2*Sk*pow(f1,2)*m2 + 16*Q2e*Sk*pow(f1,2)*m2 +
		       128*S*Sk*pow(f1,2)*m2 - 32*l2k2*Sq2*pow(f1,2)*m2 +
		       8*Q2e*Sq2*pow(f1,2)*m2 - 8*Q2h*Sq2*pow(f1,2)*m2 -
		       64*S*Sq2*pow(f1,2)*m2 + 96*Sk*Sq2*pow(f1,2)*m2 +
		       8*l2k2*Q2e*pow(f2,2)*m2 + 8*l2k2*Q2h*pow(f2,2)*m2 -
		       8*Q2e*Q2h*pow(f2,2)*m2 - 32*f1*f2*Q2e*m4 -
		       32*f1*f2*Q2h*m4 - 32*Q2e*pow(f1,2)*m4 +
		       64*S*pow(f1,2)*m4 - 64*Sk*pow(f1,2)*m4 +
		       64*Sq2*pow(f1,2)*m4 - 16*Q2e*pow(f2,2)*m4 -
		       16*Q2h*pow(f2,2)*m4 - 16*l2k2*Q2e*pow(f1,2)*M2 -
		       32*l2k2*pow(f1,2)*m2*M2 -
		       8*l2k2*pow(f1,2)*pow(Q2e,2) - 4*Q2h*pow(f1,2)*pow(Q2e,2) +
		       16*Sk*pow(f1,2)*pow(Q2e,2) - 16*Sq2*pow(f1,2)*pow(Q2e,2) +
		       32*f1*f2*m2*pow(Q2e,2) + 16*pow(f1,2)*m2*pow(Q2e,2) +
		       16*pow(f2,2)*m2*pow(Q2e,2) + 8*f1*f2*pow(Q2e,3) +
		       8*pow(f1,2)*pow(Q2e,3) + 4*pow(f2,2)*pow(Q2e,3) -
		       32*pow(f1,2)*m2*pow(Sk,2) + 32*Q2e*pow(f1,2)*pow(Sq2,2) -
		       64*pow(f1,2)*m2*pow(Sq2,2)))/2. +
		  (pow(l1k1,-1)*pow(l2k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (32*f1*f2*l2k2*Q2e*Q2h + 16*l2k2*Q2e*Q2h*pow(f1,2) -
		       32*l2k2*Q2e*S*pow(f1,2) - 16*Q2e*Q2h*S*pow(f1,2) -
		       32*Q2e*Q2h*Sk*pow(f1,2) + 128*Q2e*S*Sk*pow(f1,2) -
		       32*l2k2*Q2e*Sq2*pow(f1,2) - 64*Q2e*S*Sq2*pow(f1,2) +
		       128*Q2e*Sk*Sq2*pow(f1,2) + 8*l2k2*Q2e*Q2h*pow(f2,2) -
		       32*f1*f2*l2k2*Q2e*m2 + 64*f1*f2*l2k2*Q2h*m2 +
		       16*f1*f2*Q2e*Q2h*m2 - 48*l2k2*Q2e*pow(f1,2)*m2 +
		       32*l2k2*Q2h*pow(f1,2)*m2 - 8*Q2e*Q2h*pow(f1,2)*m2 +
		       64*l2k2*S*pow(f1,2)*m2 + 96*Q2e*S*pow(f1,2)*m2 -
		       32*Q2h*S*pow(f1,2)*m2 - 96*l2k2*Sk*pow(f1,2)*m2 +
		       112*Q2e*Sk*pow(f1,2)*m2 - 80*Q2h*Sk*pow(f1,2)*m2 +
		       64*l2k2*Sq2*pow(f1,2)*m2 - 16*Q2e*Sq2*pow(f1,2)*m2 +
		       16*Q2h*Sq2*pow(f1,2)*m2 + 128*S*Sq2*pow(f1,2)*m2 -
		       192*Sk*Sq2*pow(f1,2)*m2 - 16*l2k2*Q2e*pow(f2,2)*m2 +
		       16*l2k2*Q2h*pow(f2,2)*m2 + 16*Q2e*Q2h*pow(f2,2)*m2 -
		       64*f1*f2*Q2e*m4 + 64*f1*f2*Q2h*m4 +
		       128*S*pow(f1,2)*m4 + 128*Sq2*pow(f1,2)*m4 -
		       32*Q2e*pow(f2,2)*m4 + 32*Q2h*pow(f2,2)*m4 -
		       32*l2k2*Q2e*pow(f1,2)*M2 - 16*Q2e*Q2h*pow(f1,2)*M2 -
		       64*l2k2*pow(f1,2)*m2*M2 +
		       32*Q2e*pow(f1,2)*m2*M2 -
		       32*Q2h*pow(f1,2)*m2*M2 + 32*f1*f2*l2k2*pow(Q2e,2) +
		       16*S*pow(f1,2)*pow(Q2e,2) + 32*Sk*pow(f1,2)*pow(Q2e,2) -
		       16*Sq2*pow(f1,2)*pow(Q2e,2) + 16*l2k2*pow(f2,2)*pow(Q2e,2) +
		       4*Q2h*pow(f2,2)*pow(Q2e,2) - 16*f1*f2*m2*pow(Q2e,2) +
		       8*pow(f1,2)*m2*pow(Q2e,2) -
		       8*pow(f2,2)*m2*pow(Q2e,2) +
		       16*pow(f1,2)*M2*pow(Q2e,2) + 32*f1*f2*pow(Q2e,3) +
		       16*pow(f1,2)*pow(Q2e,3) + 16*pow(f2,2)*pow(Q2e,3) +
		       16*f1*f2*Q2e*pow(Q2h,2) + 8*Q2e*pow(f1,2)*pow(Q2h,2) +
		       4*Q2e*pow(f2,2)*pow(Q2h,2) + 32*f1*f2*m2*pow(Q2h,2) +
		       16*pow(f1,2)*m2*pow(Q2h,2) +
		       8*pow(f2,2)*m2*pow(Q2h,2) +
		       192*pow(f1,2)*m2*pow(Sk,2) - 64*Q2e*pow(f1,2)*pow(Sq2,2) +
		       128*pow(f1,2)*m2*pow(Sq2,2)))/2. +
		  (pow(l1k1,-1)*pow(l1k2,-1)*pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-32*f1*f2*l2k2*Q2e*m2 - 64*f1*f2*l2k2*Q2h*m2 +
		       16*l2k2*Q2e*pow(f1,2)*m2 - 32*l2k2*Q2h*pow(f1,2)*m2 -
		       16*Q2e*Q2h*pow(f1,2)*m2 - 128*l2k2*S*pow(f1,2)*m2 -
		       160*Q2e*S*pow(f1,2)*m2 + 32*Q2h*S*pow(f1,2)*m2 +
		       112*Q2e*Sk*pow(f1,2)*m2 - 48*Q2h*Sk*pow(f1,2)*m2 -
		       256*S*Sk*pow(f1,2)*m2 - 96*l2k2*Sq2*pow(f1,2)*m2 -
		       128*Q2e*Sq2*pow(f1,2)*m2 + 32*Q2h*Sq2*pow(f1,2)*m2 +
		       384*S*Sq2*pow(f1,2)*m2 - 192*Sk*Sq2*pow(f1,2)*m2 -
		       16*l2k2*Q2e*pow(f2,2)*m2 - 16*l2k2*Q2h*pow(f2,2)*m2 -
		       64*f1*f2*Q2e*m4 + 64*f1*f2*Q2h*m4 +
		       64*Q2e*pow(f1,2)*m4 - 64*Q2h*pow(f1,2)*m4 -
		       384*S*pow(f1,2)*m4 - 32*Q2e*pow(f2,2)*m4 +
		       32*Q2h*pow(f2,2)*m4 +
		       64*l2k2*pow(f1,2)*m2*M2 -
		       64*f1*f2*m2*pow(Q2e,2) - 16*pow(f1,2)*m2*pow(Q2e,2) -
		       32*pow(f2,2)*m2*pow(Q2e,2) +
		       256*pow(f1,2)*m2*pow(S,2) + 192*pow(f1,2)*m2*pow(Sq2,2)\
		))/2. + (pow(l1k1,-1)*pow(l2k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (-32*f1*f2*l2k1*Q2e*Q2h - 16*l2k1*Q2e*Q2h*pow(f1,2) -
		       128*l2k1*Q2e*S*pow(f1,2) + 224*l2k1*Q2h*S*pow(f1,2) +
		       16*Q2e*Q2h*S*pow(f1,2) + 128*l2k1*Q2h*Sk*pow(f1,2) +
		       48*Q2e*Q2h*Sk*pow(f1,2) - 128*l2k1*S*Sk*pow(f1,2) +
		       64*Q2e*S*Sk*pow(f1,2) - 128*Q2h*S*Sk*pow(f1,2) -
		       128*l2k1*Q2e*Sq2*pow(f1,2) + 160*l2k1*Q2h*Sq2*pow(f1,2) -
		       64*Q2e*Q2h*Sq2*pow(f1,2) - 128*l2k1*S*Sq2*pow(f1,2) +
		       64*Q2e*S*Sq2*pow(f1,2) - 448*Q2h*S*Sq2*pow(f1,2) -
		       128*l2k1*Sk*Sq2*pow(f1,2) + 64*Q2e*Sk*Sq2*pow(f1,2) -
		       64*Q2h*Sk*Sq2*pow(f1,2) - 16*l2k1*Q2e*Q2h*pow(f2,2) -
		       64*Q2h*pow(f1,2)*pow(l2k1,2) + 128*S*pow(f1,2)*pow(l2k1,2) +
		       128*Sq2*pow(f1,2)*pow(l2k1,2) + 64*f1*f2*l2k1*Q2e*m2 -
		       32*f1*f2*Q2e*Q2h*m2 - 128*l2k1*Q2e*pow(f1,2)*m2 +
		       32*l2k1*Q2h*pow(f1,2)*m2 - 64*Q2e*Q2h*pow(f1,2)*m2 -
		       64*l2k1*S*pow(f1,2)*m2 + 224*Q2e*S*pow(f1,2)*m2 -
		       160*Q2h*S*pow(f1,2)*m2 - 384*l2k1*Sk*pow(f1,2)*m2 +
		       224*Q2e*Sk*pow(f1,2)*m2 + 32*Q2h*Sk*pow(f1,2)*m2 +
		       320*l2k1*Sq2*pow(f1,2)*m2 - 64*Q2e*Sq2*pow(f1,2)*m2 -
		       32*Q2h*Sq2*pow(f1,2)*m2 + 128*S*Sq2*pow(f1,2)*m2 -
		       384*Sk*Sq2*pow(f1,2)*m2 + 32*l2k1*Q2e*pow(f2,2)*m2 -
		       16*Q2e*Q2h*pow(f2,2)*m2 +
		       192*pow(f1,2)*pow(l2k1,2)*m2 + 64*f1*f2*Q2e*m4 -
		       128*f1*f2*Q2h*m4 + 64*Q2e*pow(f1,2)*m4 -
		       96*Q2h*pow(f1,2)*m4 + 128*S*pow(f1,2)*m4 +
		       128*Sq2*pow(f1,2)*m4 + 32*Q2e*pow(f2,2)*m4 -
		       64*Q2h*pow(f2,2)*m4 + 32*l2k1*Q2h*pow(f1,2)*M2 -
		       16*Q2e*Q2h*pow(f1,2)*M2 +
		       32*Q2h*pow(f1,2)*m2*M2 + 16*f1*f2*Q2h*pow(Q2e,2) +
		       8*Q2h*pow(f1,2)*pow(Q2e,2) + 32*S*pow(f1,2)*pow(Q2e,2) +
		       32*Sq2*pow(f1,2)*pow(Q2e,2) + 8*Q2h*pow(f2,2)*pow(Q2e,2) -
		       96*f1*f2*m2*pow(Q2e,2) - 48*pow(f2,2)*m2*pow(Q2e,2) -
		       16*f1*f2*Q2e*pow(Q2h,2) + 8*Q2e*pow(f1,2)*pow(Q2h,2) -
		       48*S*pow(f1,2)*pow(Q2h,2) - 32*Sk*pow(f1,2)*pow(Q2h,2) +
		       32*Sq2*pow(f1,2)*pow(Q2h,2) + 8*l2k1*pow(f2,2)*pow(Q2h,2) -
		       12*Q2e*pow(f2,2)*pow(Q2h,2) + 128*f1*f2*m2*pow(Q2h,2) +
		       80*pow(f1,2)*m2*pow(Q2h,2) +
		       72*pow(f2,2)*m2*pow(Q2h,2) +
		       48*pow(f1,2)*M2*pow(Q2h,2) - 32*f1*f2*pow(Q2h,3) -
		       16*pow(f1,2)*pow(Q2h,3) - 4*pow(f2,2)*pow(Q2h,3) -
		       128*l2k1*pow(f1,2)*pow(S,2) + 64*Q2e*pow(f1,2)*pow(S,2) -
		       320*Q2h*pow(f1,2)*pow(S,2) + 128*pow(f1,2)*m2*pow(S,2) -
		       64*Q2h*pow(f1,2)*pow(Sk,2) + 192*pow(f1,2)*m2*pow(Sk,2) -
		       128*Q2h*pow(f1,2)*pow(Sq2,2) + 192*pow(f1,2)*m2*pow(Sq2,2)))/2. \
		+ (pow(l2k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		     (64*f1*f2*Q2e*Q2h*m2 + 128*Q2e*S*pow(f1,2)*m2 -
		       128*Q2h*S*pow(f1,2)*m2 + 128*Q2e*Sq2*pow(f1,2)*m2 -
		       128*Q2h*Sq2*pow(f1,2)*m2 + 512*S*Sq2*pow(f1,2)*m2 +
		       32*Q2e*Q2h*pow(f2,2)*m2 - 128*f1*f2*Q2h*m4 -
		       64*Q2h*pow(f1,2)*m4 - 64*Q2h*pow(f2,2)*m4 -
		       64*Q2h*pow(f1,2)*m2*M2 +
		       64*f1*f2*m2*pow(Q2h,2) + 32*pow(f1,2)*m2*pow(Q2h,2) +
		       16*pow(f2,2)*m2*pow(Q2h,2) +
		       256*pow(f1,2)*m2*pow(S,2) + 256*pow(f1,2)*m2*pow(Sq2,2)\
		))/2. + (pow(l2k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		     (64*f1*f2*Q2e*Q2h*m2 + 128*Q2e*S*pow(f1,2)*m2 -
		       128*Q2h*S*pow(f1,2)*m2 + 128*Q2e*Sq2*pow(f1,2)*m2 -
		       128*Q2h*Sq2*pow(f1,2)*m2 + 512*S*Sq2*pow(f1,2)*m2 +
		       32*Q2e*Q2h*pow(f2,2)*m2 - 128*f1*f2*Q2h*m4 -
		       64*Q2h*pow(f1,2)*m4 - 64*Q2h*pow(f2,2)*m4 -
		       64*Q2h*pow(f1,2)*m2*M2 +
		       64*f1*f2*m2*pow(Q2h,2) + 32*pow(f1,2)*m2*pow(Q2h,2) +
		       16*pow(f2,2)*m2*pow(Q2h,2) +
		       256*pow(f1,2)*m2*pow(S,2) + 256*pow(f1,2)*m2*pow(Sq2,2)\
		))/2. + (pow(l1k1,-1)*pow(l2k1,-1)*pow(l2k2,-1)*
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (64*Q2e*Q2h*S*pow(f1,2)*m2 +
		       32*Q2e*Q2h*Sk*pow(f1,2)*m2 -
		       128*Q2e*S*Sk*pow(f1,2)*m2 +
		       32*Q2e*Q2h*Sq2*pow(f1,2)*m2 -
		       128*Q2e*S*Sq2*pow(f1,2)*m2 +
		       64*Q2e*Sk*Sq2*pow(f1,2)*m2 -
		       64*Q2h*Sk*Sq2*pow(f1,2)*m2 + 32*Q2e*Q2h*pow(f1,2)*m4 -
		       128*Q2e*S*pow(f1,2)*m4 + 128*Q2h*S*pow(f1,2)*m4 -
		       64*Q2e*Sk*pow(f1,2)*m4 + 64*Q2h*Sk*pow(f1,2)*m4 -
		       256*S*Sk*pow(f1,2)*m4 - 64*Q2e*Sq2*pow(f1,2)*m4 +
		       64*Q2h*Sq2*pow(f1,2)*m4 - 256*S*Sq2*pow(f1,2)*m4 -
		       256*Sk*Sq2*pow(f1,2)*m4 + 128*f1*f2*Q2h*m6 +
		       64*Q2h*pow(f1,2)*m6 + 64*Q2h*pow(f2,2)*m6 +
		       32*Q2e*Q2h*pow(f1,2)*m2*M2 +
		       64*Q2h*pow(f1,2)*m4*M2 -
		       48*f1*f2*Q2h*m2*pow(Q2e,2) -
		       8*Q2h*pow(f1,2)*m2*pow(Q2e,2) -
		       64*S*pow(f1,2)*m2*pow(Q2e,2) -
		       32*Sk*pow(f1,2)*m2*pow(Q2e,2) -
		       32*Sq2*pow(f1,2)*m2*pow(Q2e,2) -
		       24*Q2h*pow(f2,2)*m2*pow(Q2e,2) -
		       64*f1*f2*m4*pow(Q2e,2) - 32*pow(f1,2)*m4*pow(Q2e,2) -
		       32*pow(f2,2)*m4*pow(Q2e,2) + 16*f1*f2*m2*pow(Q2e,3) +
		       8*pow(f1,2)*m2*pow(Q2e,3) +
		       8*pow(f2,2)*m2*pow(Q2e,3) -
		       32*f1*f2*Q2e*m2*pow(Q2h,2) -
		       16*Q2e*pow(f1,2)*m2*pow(Q2h,2) -
		       8*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       64*f1*f2*m4*pow(Q2h,2) - 32*pow(f1,2)*m4*pow(Q2h,2) -
		       16*pow(f2,2)*m4*pow(Q2h,2) -
		       128*Q2e*pow(f1,2)*m2*pow(S,2) -
		       256*pow(f1,2)*m4*pow(S,2) -
		       96*Q2e*pow(f1,2)*m2*pow(Sk,2) +
		       32*Q2h*pow(f1,2)*m2*pow(Sk,2) -
		       96*Q2e*pow(f1,2)*m2*pow(Sq2,2) +
		       32*Q2h*pow(f1,2)*m2*pow(Sq2,2)))/2. +
		  (pow(l1k1,-1)*pow(l1k2,-1)*pow(l2k1,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-64*Q2e*Q2h*S*pow(f1,2)*m2 +
		       32*Q2e*Q2h*Sk*pow(f1,2)*m2 +
		       128*Q2e*S*Sk*pow(f1,2)*m2 -
		       32*Q2e*Q2h*Sq2*pow(f1,2)*m2 -
		       128*Q2e*S*Sq2*pow(f1,2)*m2 +
		       192*Q2e*Sk*Sq2*pow(f1,2)*m2 -
		       64*Q2h*Sk*Sq2*pow(f1,2)*m2 + 32*Q2e*Q2h*pow(f1,2)*m4 +
		       128*Q2e*S*pow(f1,2)*m4 - 128*Q2h*S*pow(f1,2)*m4 -
		       64*Q2e*Sk*pow(f1,2)*m4 + 64*Q2h*Sk*pow(f1,2)*m4 +
		       256*S*Sk*pow(f1,2)*m4 + 64*Q2e*Sq2*pow(f1,2)*m4 -
		       64*Q2h*Sq2*pow(f1,2)*m4 - 256*S*Sq2*pow(f1,2)*m4 +
		       128*f1*f2*Q2h*m6 + 64*Q2h*pow(f1,2)*m6 +
		       64*Q2h*pow(f2,2)*m6 +
		       32*Q2e*Q2h*pow(f1,2)*m2*M2 +
		       64*Q2h*pow(f1,2)*m4*M2 -
		       48*f1*f2*Q2h*m2*pow(Q2e,2) -
		       8*Q2h*pow(f1,2)*m2*pow(Q2e,2) +
		       64*S*pow(f1,2)*m2*pow(Q2e,2) -
		       32*Sk*pow(f1,2)*m2*pow(Q2e,2) +
		       32*Sq2*pow(f1,2)*m2*pow(Q2e,2) -
		       24*Q2h*pow(f2,2)*m2*pow(Q2e,2) -
		       64*f1*f2*m4*pow(Q2e,2) - 32*pow(f1,2)*m4*pow(Q2e,2) -
		       32*pow(f2,2)*m4*pow(Q2e,2) + 16*f1*f2*m2*pow(Q2e,3) +
		       8*pow(f1,2)*m2*pow(Q2e,3) +
		       8*pow(f2,2)*m2*pow(Q2e,3) -
		       32*f1*f2*Q2e*m2*pow(Q2h,2) -
		       16*Q2e*pow(f1,2)*m2*pow(Q2h,2) -
		       8*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       64*f1*f2*m4*pow(Q2h,2) - 32*pow(f1,2)*m4*pow(Q2h,2) -
		       16*pow(f2,2)*m4*pow(Q2h,2) -
		       128*Q2e*pow(f1,2)*m2*pow(S,2) -
		       256*pow(f1,2)*m4*pow(S,2) -
		       96*Q2e*pow(f1,2)*m2*pow(Sk,2) +
		       32*Q2h*pow(f1,2)*m2*pow(Sk,2) -
		       96*Q2e*pow(f1,2)*m2*pow(Sq2,2) +
		       32*Q2h*pow(f1,2)*m2*pow(Sq2,2)))/2. +
		  (pow(l1k2,-1)*pow(l2k1,-1)/M2*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (32*l1k1*Q2e*Q2h*S*pow(f2,2) + 32*l1k1*Q2h*S*Sk*pow(f2,2) +
		       16*Q2e*Q2h*S*Sk*pow(f2,2) - 32*l1k1*Q2h*S*Sq2*pow(f2,2) -
		       16*Q2e*Q2h*S*Sq2*pow(f2,2) + 32*Q2h*S*pow(f2,2)*pow(l1k1,2) -
		       40*l1k1*Q2e*Q2h*pow(f2,2)*m2 +
		       16*l1k1*Q2h*S*pow(f2,2)*m2 +
		       56*Q2e*Q2h*S*pow(f2,2)*m2 -
		       96*l1k1*Q2h*Sk*pow(f2,2)*m2 -
		       56*Q2e*Q2h*Sk*pow(f2,2)*m2 +
		       96*l1k1*Q2h*Sq2*pow(f2,2)*m2 +
		       72*Q2e*Q2h*Sq2*pow(f2,2)*m2 -
		       32*Q2h*S*Sq2*pow(f2,2)*m2 +
		       96*Q2h*Sk*Sq2*pow(f2,2)*m2 -
		       48*Q2h*pow(f2,2)*pow(l1k1,2)*m2 -
		       8*Q2e*Q2h*pow(f2,2)*m4 + 32*Q2h*S*pow(f2,2)*m4 +
		       8*Q2h*S*pow(f2,2)*pow(Q2e,2) -
		       12*Q2h*pow(f2,2)*m2*pow(Q2e,2) -
		       56*l1k1*S*pow(f2,2)*pow(Q2h,2) + 4*Q2e*S*pow(f2,2)*pow(Q2h,2) +
		       32*l1k1*Sk*pow(f2,2)*pow(Q2h,2) - 12*Q2e*Sk*pow(f2,2)*pow(Q2h,2) -
		       32*S*Sk*pow(f2,2)*pow(Q2h,2) - 16*l1k1*Sq2*pow(f2,2)*pow(Q2h,2) +
		       20*Q2e*Sq2*pow(f2,2)*pow(Q2h,2) + 48*S*Sq2*pow(f2,2)*pow(Q2h,2) -
		       16*Sk*Sq2*pow(f2,2)*pow(Q2h,2) +
		       16*pow(f2,2)*pow(l1k1,2)*pow(Q2h,2) +
		       8*l1k1*pow(f2,2)*m2*pow(Q2h,2) +
		       12*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       40*S*pow(f2,2)*m2*pow(Q2h,2) -
		       8*Sk*pow(f2,2)*m2*pow(Q2h,2) -
		       32*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       8*pow(f2,2)*m4*pow(Q2h,2) - 4*Q2e*pow(f2,2)*pow(Q2h,3) -
		       12*S*pow(f2,2)*pow(Q2h,3) + 8*Sk*pow(f2,2)*pow(Q2h,3) -
		       20*Sq2*pow(f2,2)*pow(Q2h,3) - 4*pow(f2,2)*m2*pow(Q2h,3) -
		       32*l1k1*Q2h*pow(f2,2)*pow(S,2) - 16*Q2e*Q2h*pow(f2,2)*pow(S,2) -
		       32*Q2h*pow(f2,2)*m2*pow(S,2) +
		       80*pow(f2,2)*pow(Q2h,2)*pow(S,2) -
		       48*Q2h*pow(f2,2)*m2*pow(Sk,2) +
		       16*pow(f2,2)*pow(Q2h,2)*pow(Sk,2) -
		       48*Q2h*pow(f2,2)*m2*pow(Sq2,2)))/2. +
		  (pow(l2k1,-1)*pow(l2k2,-1)/M2*
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (8*l1k2*Q2e*Q2h*pow(f2,2)*m2 +
		       32*l1k2*Q2h*S*pow(f2,2)*m2 -
		       40*Q2e*Q2h*S*pow(f2,2)*m2 -
		       28*Q2e*Q2h*Sk*pow(f2,2)*m2 -
		       64*Q2h*S*Sk*pow(f2,2)*m2 +
		       8*l1k2*Q2h*Sq2*pow(f2,2)*m2 -
		       8*Q2e*Q2h*Sq2*pow(f2,2)*m2 -
		       32*Q2h*S*Sq2*pow(f2,2)*m2 -
		       16*Q2h*Sk*Sq2*pow(f2,2)*m2 - 24*Q2e*Q2h*pow(f2,2)*m4 -
		       96*Q2h*S*pow(f2,2)*m4 - 96*Q2h*Sq2*pow(f2,2)*m4 -
		       4*Q2h*pow(f2,2)*m2*pow(Q2e,2) +
		       4*Q2e*pow(f2,2)*m2*pow(Q2h,2) +
		       8*S*pow(f2,2)*m2*pow(Q2h,2) +
		       12*Sk*pow(f2,2)*m2*pow(Q2h,2) +
		       24*pow(f2,2)*m4*pow(Q2h,2) -
		       64*Q2h*pow(f2,2)*m2*pow(S,2) -
		       16*Q2h*pow(f2,2)*m2*pow(Sq2,2)))/2. +
		  (pow(l1k1,-1)*pow(l1k2,-1)*pow(l2k1,-1)/M2*
		     (4*l2k2*Q2e*Q2h*S*pow(f2,2) - 4*l2k2*Q2e*Q2h*Sk*pow(f2,2) +
		       4*l2k2*Q2e*Q2h*Sq2*pow(f2,2) + 8*Q2e*Q2h*S*Sq2*pow(f2,2) -
		       8*Q2e*Q2h*Sk*Sq2*pow(f2,2) - 4*l2k2*Q2e*Q2h*pow(f2,2)*m2 -
		       8*l2k2*Q2h*S*pow(f2,2)*m2 +
		       16*Q2e*Q2h*S*pow(f2,2)*m2 -
		       4*l2k2*Q2h*Sk*pow(f2,2)*m2 +
		       4*Q2e*Q2h*Sk*pow(f2,2)*m2 + 32*Q2h*S*Sk*pow(f2,2)*m2 -
		       8*l2k2*Q2h*Sq2*pow(f2,2)*m2 +
		       2*Q2e*Q2h*Sq2*pow(f2,2)*m2 -
		       16*Q2h*S*Sq2*pow(f2,2)*m2 +
		       24*Q2h*Sk*Sq2*pow(f2,2)*m2 - 4*Q2e*Q2h*pow(f2,2)*m4 +
		       16*Q2h*S*pow(f2,2)*m4 - 16*Q2h*Sk*pow(f2,2)*m4 +
		       16*Q2h*Sq2*pow(f2,2)*m4 - 2*l2k2*Q2h*pow(f2,2)*pow(Q2e,2) +
		       4*Q2h*Sk*pow(f2,2)*pow(Q2e,2) - 4*Q2h*Sq2*pow(f2,2)*pow(Q2e,2) +
		       Q2h*pow(f2,2)*pow(Q2e,3) + 2*Q2e*Sq2*pow(f2,2)*pow(Q2h,2) -
		       8*S*pow(f2,2)*m2*pow(Q2h,2) -
		       2*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       4*pow(f2,2)*m4*pow(Q2h,2) - pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) -
		       8*Q2h*pow(f2,2)*m2*pow(Sk,2) +
		       8*Q2e*Q2h*pow(f2,2)*pow(Sq2,2) -
		       16*Q2h*pow(f2,2)*m2*pow(Sq2,2)))/2. +
		  (pow(l1k1,-1)*pow(l2k1,-1)/M2*
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (-8*l2k2*Q2e*Q2h*S*pow(f2,2) + 32*Q2e*Q2h*S*Sk*pow(f2,2) -
		       8*l2k2*Q2e*Q2h*Sq2*pow(f2,2) - 16*Q2e*Q2h*S*Sq2*pow(f2,2) +
		       32*Q2e*Q2h*Sk*Sq2*pow(f2,2) - 8*l2k2*Q2e*Q2h*pow(f2,2)*m2 +
		       16*l2k2*Q2h*S*pow(f2,2)*m2 +
		       24*Q2e*Q2h*S*pow(f2,2)*m2 -
		       24*l2k2*Q2h*Sk*pow(f2,2)*m2 +
		       28*Q2e*Q2h*Sk*pow(f2,2)*m2 +
		       16*l2k2*Q2h*Sq2*pow(f2,2)*m2 -
		       4*Q2e*Q2h*Sq2*pow(f2,2)*m2 +
		       32*Q2h*S*Sq2*pow(f2,2)*m2 -
		       48*Q2h*Sk*Sq2*pow(f2,2)*m2 + 8*Q2e*Q2h*pow(f2,2)*m4 +
		       32*Q2h*S*pow(f2,2)*m4 + 32*Q2h*Sq2*pow(f2,2)*m4 -
		       4*l2k2*Q2h*pow(f2,2)*pow(Q2e,2) + 4*Q2h*S*pow(f2,2)*pow(Q2e,2) +
		       8*Q2h*Sk*pow(f2,2)*pow(Q2e,2) - 4*Q2h*Sq2*pow(f2,2)*pow(Q2e,2) +
		       4*Q2h*pow(f2,2)*m2*pow(Q2e,2) -
		       4*Q2e*S*pow(f2,2)*pow(Q2h,2) - 8*Q2e*Sk*pow(f2,2)*pow(Q2h,2) -
		       4*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       8*S*pow(f2,2)*m2*pow(Q2h,2) -
		       20*Sk*pow(f2,2)*m2*pow(Q2h,2) +
		       4*Sq2*pow(f2,2)*m2*pow(Q2h,2) -
		       8*pow(f2,2)*m4*pow(Q2h,2) +
		       48*Q2h*pow(f2,2)*m2*pow(Sk,2) -
		       16*Q2e*Q2h*pow(f2,2)*pow(Sq2,2) +
		       32*Q2h*pow(f2,2)*m2*pow(Sq2,2)))/2. +
		  (pow(l1k1,-1)*pow(l1k2,-1)/M2*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (8*l2k2*Q2e*Q2h*pow(f2,2)*m2 -
		       32*l2k2*Q2h*S*pow(f2,2)*m2 -
		       40*Q2e*Q2h*S*pow(f2,2)*m2 +
		       28*Q2e*Q2h*Sk*pow(f2,2)*m2 -
		       64*Q2h*S*Sk*pow(f2,2)*m2 -
		       24*l2k2*Q2h*Sq2*pow(f2,2)*m2 -
		       32*Q2e*Q2h*Sq2*pow(f2,2)*m2 +
		       96*Q2h*S*Sq2*pow(f2,2)*m2 -
		       48*Q2h*Sk*Sq2*pow(f2,2)*m2 + 24*Q2e*Q2h*pow(f2,2)*m4 -
		       96*Q2h*S*pow(f2,2)*m4 +
		       4*Q2h*pow(f2,2)*m2*pow(Q2e,2) -
		       4*Q2e*pow(f2,2)*m2*pow(Q2h,2) +
		       8*S*pow(f2,2)*m2*pow(Q2h,2) -
		       12*Sk*pow(f2,2)*m2*pow(Q2h,2) +
		       8*Sq2*pow(f2,2)*m2*pow(Q2h,2) -
		       24*pow(f2,2)*m4*pow(Q2h,2) +
		       64*Q2h*pow(f2,2)*m2*pow(S,2) +
		       48*Q2h*pow(f2,2)*m2*pow(Sq2,2)))/2. +
		  (pow(l2k1,-1)/M2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		     (32*Q2e*Q2h*S*pow(f2,2)*m2 +
		       32*Q2e*Q2h*Sq2*pow(f2,2)*m2 +
		       128*Q2h*S*Sq2*pow(f2,2)*m2 -
		       8*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       32*S*pow(f2,2)*m2*pow(Q2h,2) -
		       32*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       64*Q2h*pow(f2,2)*m2*pow(S,2) +
		       64*Q2h*pow(f2,2)*m2*pow(Sq2,2)))/2. +
		  (pow(l2k2,-1)/M2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		     (32*Q2e*Q2h*S*pow(f2,2)*m2 +
		       32*Q2e*Q2h*Sq2*pow(f2,2)*m2 +
		       128*Q2h*S*Sq2*pow(f2,2)*m2 -
		       8*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       32*S*pow(f2,2)*m2*pow(Q2h,2) -
		       32*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       64*Q2h*pow(f2,2)*m2*pow(S,2) +
		       64*Q2h*pow(f2,2)*m2*pow(Sq2,2)))/2. +
		  (pow(l1k2,-1)*pow(l2k1,-2)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (-64*Q2e*Q2h*S*pow(f1,2)*m2 +
		       32*Q2e*Q2h*Sk*pow(f1,2)*m2 +
		       128*Q2h*S*Sk*pow(f1,2)*m2 -
		       64*Q2e*Q2h*Sq2*pow(f1,2)*m2 -
		       256*Q2h*S*Sq2*pow(f1,2)*m2 +
		       128*Q2h*Sk*Sq2*pow(f1,2)*m2 - 64*f1*f2*Q2e*Q2h*m4 -
		       128*Q2e*S*pow(f1,2)*m4 + 128*Q2h*S*pow(f1,2)*m4 +
		       64*Q2e*Sk*pow(f1,2)*m4 - 64*Q2h*Sk*pow(f1,2)*m4 +
		       256*S*Sk*pow(f1,2)*m4 - 128*Q2e*Sq2*pow(f1,2)*m4 +
		       128*Q2h*Sq2*pow(f1,2)*m4 - 512*S*Sq2*pow(f1,2)*m4 +
		       256*Sk*Sq2*pow(f1,2)*m4 - 32*Q2e*Q2h*pow(f2,2)*m4 +
		       128*f1*f2*Q2h*m6 + 64*Q2h*pow(f1,2)*m6 +
		       64*Q2h*pow(f2,2)*m6 + 64*Q2h*pow(f1,2)*m4*M2 -
		       32*f1*f2*Q2e*m2*pow(Q2h,2) +
		       64*S*pow(f1,2)*m2*pow(Q2h,2) -
		       32*Sk*pow(f1,2)*m2*pow(Q2h,2) +
		       64*Sq2*pow(f1,2)*m2*pow(Q2h,2) -
		       16*Q2e*pow(f2,2)*m2*pow(Q2h,2) +
		       16*pow(f2,2)*m4*pow(Q2h,2) +
		       32*pow(f1,2)*m2*M2*pow(Q2h,2) -
		       32*f1*f2*m2*pow(Q2h,3) - 16*pow(f1,2)*m2*pow(Q2h,3) -
		       8*pow(f2,2)*m2*pow(Q2h,3) -
		       128*Q2h*pow(f1,2)*m2*pow(S,2) -
		       256*pow(f1,2)*m4*pow(S,2) -
		       64*Q2h*pow(f1,2)*m2*pow(Sk,2) -
		       128*Q2h*pow(f1,2)*m2*pow(Sq2,2) -
		       256*pow(f1,2)*m4*pow(Sq2,2)))/2. +
		  (pow(l2k1,-1)*pow(l2k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		     (128*f1*f2*Q2e*Q2h*m4 + 256*Q2e*S*pow(f1,2)*m4 -
		       256*Q2h*S*pow(f1,2)*m4 + 256*Q2e*Sq2*pow(f1,2)*m4 -
		       256*Q2h*Sq2*pow(f1,2)*m4 + 1024*S*Sq2*pow(f1,2)*m4 +
		       64*Q2e*Q2h*pow(f2,2)*m4 - 256*f1*f2*Q2h*m6 -
		       128*Q2h*pow(f1,2)*m6 - 128*Q2h*pow(f2,2)*m6 -
		       128*Q2h*pow(f1,2)*m4*M2 +
		       128*f1*f2*m4*pow(Q2h,2) + 64*pow(f1,2)*m4*pow(Q2h,2) +
		       32*pow(f2,2)*m4*pow(Q2h,2) +
		       512*pow(f1,2)*m4*pow(S,2) + 512*pow(f1,2)*m4*pow(Sq2,2)\
		))/2. + (pow(l2k1,-1)*pow(l2k2,-1)/M2*
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		     (64*Q2e*Q2h*S*pow(f2,2)*m4 +
		       64*Q2e*Q2h*Sq2*pow(f2,2)*m4 +
		       256*Q2h*S*Sq2*pow(f2,2)*m4 -
		       16*Q2e*pow(f2,2)*m4*pow(Q2h,2) -
		       64*S*pow(f2,2)*m4*pow(Q2h,2) -
		       64*Sq2*pow(f2,2)*m4*pow(Q2h,2) +
		       128*Q2h*pow(f2,2)*m4*pow(S,2) +
		       128*Q2h*pow(f2,2)*m4*pow(Sq2,2)))/2. +
		  (pow(l1k1,-1)*pow(l2k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-256*Q2e*Q2h*S*Sk*pow(f1,2) - 128*Q2e*Q2h*S*Sq2*pow(f1,2) -
		       64*Q2e*Q2h*Sk*Sq2*pow(f1,2) + 64*Q2e*Q2h*Sk*pow(f1,2)*m2 +
		       512*Q2h*S*Sq2*pow(f1,2)*m2 +
		       128*Q2h*Sk*Sq2*pow(f1,2)*m2 +
		       64*Q2e*Q2h*pow(f1,2)*m4 - 128*Q2e*Sq2*pow(f1,2)*m4 +
		       128*Q2h*Sq2*pow(f1,2)*m4 + 512*S*Sq2*pow(f1,2)*m4 -
		       256*f1*f2*Q2h*m6 - 128*Q2h*pow(f1,2)*m6 -
		       128*Q2h*pow(f2,2)*m6 -
		       128*Q2h*pow(f1,2)*m4*M2 +
		       64*Q2h*S*pow(f1,2)*pow(Q2e,2) + 32*Q2h*Sk*pow(f1,2)*pow(Q2e,2) +
		       32*Q2h*Sq2*pow(f1,2)*pow(Q2e,2) +
		       64*f1*f2*Q2h*m2*pow(Q2e,2) +
		       32*Q2h*pow(f2,2)*m2*pow(Q2e,2) +
		       128*f1*f2*m4*pow(Q2e,2) + 64*pow(f2,2)*m4*pow(Q2e,2) -
		       128*Q2e*S*pow(f1,2)*pow(Q2h,2) - 64*Q2e*Sk*pow(f1,2)*pow(Q2h,2) +
		       256*S*Sk*pow(f1,2)*pow(Q2h,2) - 32*Q2e*Sq2*pow(f1,2)*pow(Q2h,2) +
		       256*S*Sq2*pow(f1,2)*pow(Q2h,2) + 64*Sk*Sq2*pow(f1,2)*pow(Q2h,2) +
		       224*f1*f2*Q2e*m2*pow(Q2h,2) +
		       112*Q2e*pow(f1,2)*m2*pow(Q2h,2) -
		       64*Sk*pow(f1,2)*m2*pow(Q2h,2) +
		       112*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       64*pow(f1,2)*m4*pow(Q2h,2) -
		       32*pow(f2,2)*m4*pow(Q2h,2) -
		       128*pow(f1,2)*m2*M2*pow(Q2h,2) +
		       16*f1*f2*pow(Q2e,2)*pow(Q2h,2) + 8*pow(f1,2)*pow(Q2e,2)*pow(Q2h,2) +
		       8*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) + 32*f1*f2*Q2e*pow(Q2h,3) +
		       64*S*pow(f1,2)*pow(Q2h,3) + 32*Sk*pow(f1,2)*pow(Q2h,3) +
		       16*Q2e*pow(f2,2)*pow(Q2h,3) - 32*f1*f2*m2*pow(Q2h,3) -
		       48*pow(f1,2)*m2*pow(Q2h,3) -
		       48*pow(f2,2)*m2*pow(Q2h,3) -
		       32*pow(f1,2)*M2*pow(Q2h,3) + 16*f1*f2*pow(Q2h,4) +
		       8*pow(f1,2)*pow(Q2h,4) - 256*Q2e*Q2h*pow(f1,2)*pow(S,2) +
		       512*Q2h*pow(f1,2)*m2*pow(S,2) +
		       512*pow(f1,2)*m4*pow(S,2) +
		       384*pow(f1,2)*pow(Q2h,2)*pow(S,2) - 64*Q2e*Q2h*pow(f1,2)*pow(Sk,2) -
		       128*Q2h*pow(f1,2)*m2*pow(Sk,2) +
		       64*pow(f1,2)*pow(Q2h,2)*pow(Sk,2) +
		       128*Q2h*pow(f1,2)*m2*pow(Sq2,2) +
		       64*pow(f1,2)*pow(Q2h,2)*pow(Sq2,2)))/2. +
		  (pow(l1k2,-1)*pow(l2k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (256*Q2e*Q2h*S*Sk*pow(f1,2) - 384*Q2e*Q2h*S*Sq2*pow(f1,2) +
		       192*Q2e*Q2h*Sk*Sq2*pow(f1,2) + 64*Q2e*Q2h*Sk*pow(f1,2)*m2 +
		       512*Q2h*S*Sq2*pow(f1,2)*m2 +
		       128*Q2h*Sk*Sq2*pow(f1,2)*m2 +
		       64*Q2e*Q2h*pow(f1,2)*m4 - 128*Q2e*Sq2*pow(f1,2)*m4 +
		       128*Q2h*Sq2*pow(f1,2)*m4 + 512*S*Sq2*pow(f1,2)*m4 -
		       256*f1*f2*Q2h*m6 - 128*Q2h*pow(f1,2)*m6 -
		       128*Q2h*pow(f2,2)*m6 -
		       128*Q2h*pow(f1,2)*m4*M2 -
		       64*Q2h*S*pow(f1,2)*pow(Q2e,2) + 32*Q2h*Sk*pow(f1,2)*pow(Q2e,2) -
		       32*Q2h*Sq2*pow(f1,2)*pow(Q2e,2) +
		       64*f1*f2*Q2h*m2*pow(Q2e,2) +
		       32*Q2h*pow(f2,2)*m2*pow(Q2e,2) +
		       128*f1*f2*m4*pow(Q2e,2) + 64*pow(f2,2)*m4*pow(Q2e,2) +
		       128*Q2e*S*pow(f1,2)*pow(Q2h,2) - 64*Q2e*Sk*pow(f1,2)*pow(Q2h,2) -
		       256*S*Sk*pow(f1,2)*pow(Q2h,2) + 96*Q2e*Sq2*pow(f1,2)*pow(Q2h,2) +
		       512*S*Sq2*pow(f1,2)*pow(Q2h,2) - 192*Sk*Sq2*pow(f1,2)*pow(Q2h,2) +
		       224*f1*f2*Q2e*m2*pow(Q2h,2) +
		       112*Q2e*pow(f1,2)*m2*pow(Q2h,2) -
		       64*Sk*pow(f1,2)*m2*pow(Q2h,2) +
		       112*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       64*pow(f1,2)*m4*pow(Q2h,2) -
		       32*pow(f2,2)*m4*pow(Q2h,2) -
		       128*pow(f1,2)*m2*M2*pow(Q2h,2) +
		       16*f1*f2*pow(Q2e,2)*pow(Q2h,2) + 8*pow(f1,2)*pow(Q2e,2)*pow(Q2h,2) +
		       8*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) + 32*f1*f2*Q2e*pow(Q2h,3) -
		       64*S*pow(f1,2)*pow(Q2h,3) + 32*Sk*pow(f1,2)*pow(Q2h,3) -
		       64*Sq2*pow(f1,2)*pow(Q2h,3) + 16*Q2e*pow(f2,2)*pow(Q2h,3) -
		       32*f1*f2*m2*pow(Q2h,3) - 48*pow(f1,2)*m2*pow(Q2h,3) -
		       48*pow(f2,2)*m2*pow(Q2h,3) -
		       32*pow(f1,2)*M2*pow(Q2h,3) + 16*f1*f2*pow(Q2h,4) +
		       8*pow(f1,2)*pow(Q2h,4) - 256*Q2e*Q2h*pow(f1,2)*pow(S,2) +
		       512*Q2h*pow(f1,2)*m2*pow(S,2) +
		       512*pow(f1,2)*m4*pow(S,2) +
		       384*pow(f1,2)*pow(Q2h,2)*pow(S,2) - 64*Q2e*Q2h*pow(f1,2)*pow(Sk,2) -
		       128*Q2h*pow(f1,2)*m2*pow(Sk,2) +
		       64*pow(f1,2)*pow(Q2h,2)*pow(Sk,2) -
		       128*Q2e*Q2h*pow(f1,2)*pow(Sq2,2) +
		       128*Q2h*pow(f1,2)*m2*pow(Sq2,2) +
		       192*pow(f1,2)*pow(Q2h,2)*pow(Sq2,2)))/2. +
		  (pow(l2k1,-1)/M2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (32*Q2e*Q2h*S*Sq2*pow(f2,2) - 16*l1k2*Q2e*Q2h*pow(f2,2)*m2 +
		       32*Q2e*Q2h*Sk*pow(f2,2)*m2 -
		       64*Q2h*S*Sq2*pow(f2,2)*m2 - 8*l1k2*Q2h*pow(f2,2)*pow(Q2e,2) -
		       32*Q2h*S*pow(f2,2)*pow(Q2e,2) + 16*Q2h*Sk*pow(f2,2)*pow(Q2e,2) -
		       32*Q2h*Sq2*pow(f2,2)*pow(Q2e,2) - 8*Q2h*pow(f2,2)*pow(Q2e,3) +
		       8*l1k2*Q2e*pow(f2,2)*pow(Q2h,2) - 96*S*Sq2*pow(f2,2)*pow(Q2h,2) -
		       8*Q2e*pow(f2,2)*m2*pow(Q2h,2) +
		       32*Sk*pow(f2,2)*m2*pow(Q2h,2) -
		       32*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       12*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) + 4*Q2e*pow(f2,2)*pow(Q2h,3) +
		       32*S*pow(f2,2)*pow(Q2h,3) - 16*Sk*pow(f2,2)*pow(Q2h,3) +
		       32*Sq2*pow(f2,2)*pow(Q2h,3) + 16*pow(f2,2)*m2*pow(Q2h,3) +
		       32*Q2e*Q2h*pow(f2,2)*pow(S,2) -
		       64*Q2h*pow(f2,2)*m2*pow(S,2) -
		       96*pow(f2,2)*pow(Q2h,2)*pow(S,2) - 32*pow(f2,2)*pow(Q2h,2)*pow(Sq2,2)\
		))/2. + (pow(l2k2,-1)/M2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (32*Q2e*Q2h*S*Sq2*pow(f2,2) + 16*l1k2*Q2e*Q2h*pow(f2,2)*m2 -
		       32*Q2e*Q2h*Sk*pow(f2,2)*m2 +
		       32*Q2e*Q2h*Sq2*pow(f2,2)*m2 -
		       64*Q2h*S*Sq2*pow(f2,2)*m2 + 8*l1k2*Q2h*pow(f2,2)*pow(Q2e,2) -
		       32*Q2h*S*pow(f2,2)*pow(Q2e,2) - 16*Q2h*Sk*pow(f2,2)*pow(Q2e,2) -
		       16*Q2h*Sq2*pow(f2,2)*pow(Q2e,2) +
		       24*Q2h*pow(f2,2)*m2*pow(Q2e,2) + 4*Q2h*pow(f2,2)*pow(Q2e,3) -
		       8*l1k2*Q2e*pow(f2,2)*pow(Q2h,2) - 96*S*Sq2*pow(f2,2)*pow(Q2h,2) -
		       16*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       32*Sk*pow(f2,2)*m2*pow(Q2h,2) -
		       4*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) + 32*S*pow(f2,2)*pow(Q2h,3) +
		       16*Sk*pow(f2,2)*pow(Q2h,3) + 16*Sq2*pow(f2,2)*pow(Q2h,3) +
		       8*pow(f2,2)*pow(Q2h,4) + 32*Q2e*Q2h*pow(f2,2)*pow(S,2) -
		       64*Q2h*pow(f2,2)*m2*pow(S,2) -
		       96*pow(f2,2)*pow(Q2h,2)*pow(S,2) - 32*pow(f2,2)*pow(Q2h,2)*pow(Sq2,2)\
		))/2. + (pow(l1k1,-1)*pow(l2k2,-1)/M2*
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (-32*l2k1*Q2e*Q2h*S*pow(f2,2) - 32*l2k1*Q2h*S*Sk*pow(f2,2) +
		       16*Q2e*Q2h*S*Sk*pow(f2,2) - 32*l2k1*Q2e*Q2h*Sq2*pow(f2,2) -
		       32*l2k1*Q2h*S*Sq2*pow(f2,2) + 16*Q2e*Q2h*S*Sq2*pow(f2,2) -
		       32*l2k1*Q2h*Sk*Sq2*pow(f2,2) + 16*Q2e*Q2h*Sk*Sq2*pow(f2,2) +
		       32*Q2h*S*pow(f2,2)*pow(l2k1,2) + 32*Q2h*Sq2*pow(f2,2)*pow(l2k1,2) -
		       40*l2k1*Q2e*Q2h*pow(f2,2)*m2 -
		       16*l2k1*Q2h*S*pow(f2,2)*m2 +
		       56*Q2e*Q2h*S*pow(f2,2)*m2 -
		       96*l2k1*Q2h*Sk*pow(f2,2)*m2 +
		       56*Q2e*Q2h*Sk*pow(f2,2)*m2 +
		       80*l2k1*Q2h*Sq2*pow(f2,2)*m2 -
		       16*Q2e*Q2h*Sq2*pow(f2,2)*m2 +
		       32*Q2h*S*Sq2*pow(f2,2)*m2 -
		       96*Q2h*Sk*Sq2*pow(f2,2)*m2 +
		       48*Q2h*pow(f2,2)*pow(l2k1,2)*m2 +
		       8*Q2e*Q2h*pow(f2,2)*m4 + 32*Q2h*S*pow(f2,2)*m4 +
		       32*Q2h*Sq2*pow(f2,2)*m4 + 8*Q2h*S*pow(f2,2)*pow(Q2e,2) +
		       8*Q2h*Sq2*pow(f2,2)*pow(Q2e,2) +
		       12*Q2h*pow(f2,2)*m2*pow(Q2e,2) +
		       56*l2k1*S*pow(f2,2)*pow(Q2h,2) + 4*Q2e*S*pow(f2,2)*pow(Q2h,2) +
		       32*l2k1*Sk*pow(f2,2)*pow(Q2h,2) + 12*Q2e*Sk*pow(f2,2)*pow(Q2h,2) -
		       32*S*Sk*pow(f2,2)*pow(Q2h,2) + 40*l2k1*Sq2*pow(f2,2)*pow(Q2h,2) -
		       16*Q2e*Sq2*pow(f2,2)*pow(Q2h,2) - 112*S*Sq2*pow(f2,2)*pow(Q2h,2) -
		       16*Sk*Sq2*pow(f2,2)*pow(Q2h,2) -
		       16*pow(f2,2)*pow(l2k1,2)*pow(Q2h,2) +
		       8*l2k1*pow(f2,2)*m2*pow(Q2h,2) -
		       12*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       40*S*pow(f2,2)*m2*pow(Q2h,2) +
		       8*Sk*pow(f2,2)*m2*pow(Q2h,2) -
		       8*Sq2*pow(f2,2)*m2*pow(Q2h,2) -
		       8*pow(f2,2)*m4*pow(Q2h,2) + 4*Q2e*pow(f2,2)*pow(Q2h,3) -
		       12*S*pow(f2,2)*pow(Q2h,3) - 8*Sk*pow(f2,2)*pow(Q2h,3) +
		       8*Sq2*pow(f2,2)*pow(Q2h,3) + 4*pow(f2,2)*m2*pow(Q2h,3) -
		       32*l2k1*Q2h*pow(f2,2)*pow(S,2) + 16*Q2e*Q2h*pow(f2,2)*pow(S,2) +
		       32*Q2h*pow(f2,2)*m2*pow(S,2) -
		       80*pow(f2,2)*pow(Q2h,2)*pow(S,2) +
		       48*Q2h*pow(f2,2)*m2*pow(Sk,2) -
		       16*pow(f2,2)*pow(Q2h,2)*pow(Sk,2) +
		       48*Q2h*pow(f2,2)*m2*pow(Sq2,2) -
		       32*pow(f2,2)*pow(Q2h,2)*pow(Sq2,2)))/2. +
		  (pow(l1k2,-1)*pow(l2k1,-1)/M2*
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (16*l2k2*Q2e*Q2h*S*pow(f2,2) - 32*l2k2*Q2h*S*Sk*pow(f2,2) +
		       16*l2k2*Q2e*Q2h*Sq2*pow(f2,2) + 128*l2k2*Q2h*S*Sq2*pow(f2,2) -
		       32*Q2e*Q2h*S*Sq2*pow(f2,2) - 32*l2k2*Q2h*Sk*Sq2*pow(f2,2) -
		       8*l2k2*Q2e*Q2h*pow(f2,2)*m2 +
		       16*Q2e*Q2h*Sk*pow(f2,2)*m2 +
		       32*Q2h*S*Sk*pow(f2,2)*m2 -
		       16*Q2e*Q2h*Sq2*pow(f2,2)*m2 -
		       64*Q2h*S*Sq2*pow(f2,2)*m2 +
		       32*Q2h*Sk*Sq2*pow(f2,2)*m2 - 16*Q2e*Q2h*pow(f2,2)*m4 -
		       64*Q2h*S*pow(f2,2)*m4 - 64*Q2h*Sq2*pow(f2,2)*m4 -
		       8*Q2h*S*pow(f2,2)*pow(Q2e,2) - 8*Q2h*Sq2*pow(f2,2)*pow(Q2e,2) +
		       4*Q2h*pow(f2,2)*m2*pow(Q2e,2) -
		       4*l2k2*Q2e*pow(f2,2)*pow(Q2h,2) - 16*l2k2*S*pow(f2,2)*pow(Q2h,2) +
		       8*Q2e*S*pow(f2,2)*pow(Q2h,2) + 32*S*Sk*pow(f2,2)*pow(Q2h,2) -
		       16*l2k2*Sq2*pow(f2,2)*pow(Q2h,2) + 8*Q2e*Sq2*pow(f2,2)*pow(Q2h,2) -
		       64*S*Sq2*pow(f2,2)*pow(Q2h,2) + 16*Sk*Sq2*pow(f2,2)*pow(Q2h,2) -
		       32*S*pow(f2,2)*m2*pow(Q2h,2) +
		       16*Sk*pow(f2,2)*m2*pow(Q2h,2) -
		       16*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       16*pow(f2,2)*m4*pow(Q2h,2) +
		       4*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) +
		       8*pow(f2,2)*m2*pow(Q2h,3) + 64*l2k2*Q2h*pow(f2,2)*pow(S,2) -
		       16*Q2e*Q2h*pow(f2,2)*pow(S,2) -
		       32*Q2h*pow(f2,2)*m2*pow(S,2) -
		       48*pow(f2,2)*pow(Q2h,2)*pow(S,2) +
		       64*l2k2*Q2h*pow(f2,2)*pow(Sq2,2) - 16*Q2e*Q2h*pow(f2,2)*pow(Sq2,2) -
		       32*Q2h*pow(f2,2)*m2*pow(Sq2,2) -
		       16*pow(f2,2)*pow(Q2h,2)*pow(Sq2,2)))/2. +
		  (pow(l1k1,-1)/M2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-32*Q2e*Q2h*S*Sq2*pow(f2,2) - 16*l2k2*Q2e*Q2h*pow(f2,2)*m2 -
		       32*Q2e*Q2h*Sk*pow(f2,2)*m2 +
		       64*Q2h*S*Sq2*pow(f2,2)*m2 - 8*l2k2*Q2h*pow(f2,2)*pow(Q2e,2) -
		       32*Q2h*S*pow(f2,2)*pow(Q2e,2) - 16*Q2h*Sk*pow(f2,2)*pow(Q2e,2) +
		       8*Q2h*pow(f2,2)*pow(Q2e,3) + 8*l2k2*Q2e*pow(f2,2)*pow(Q2h,2) +
		       96*S*Sq2*pow(f2,2)*pow(Q2h,2) +
		       8*Q2e*pow(f2,2)*m2*pow(Q2h,2) -
		       32*Sk*pow(f2,2)*m2*pow(Q2h,2) +
		       32*Sq2*pow(f2,2)*m2*pow(Q2h,2) -
		       12*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) - 4*Q2e*pow(f2,2)*pow(Q2h,3) +
		       32*S*pow(f2,2)*pow(Q2h,3) + 16*Sk*pow(f2,2)*pow(Q2h,3) -
		       16*pow(f2,2)*m2*pow(Q2h,3) - 32*Q2e*Q2h*pow(f2,2)*pow(S,2) +
		       64*Q2h*pow(f2,2)*m2*pow(S,2) +
		       96*pow(f2,2)*pow(Q2h,2)*pow(S,2) + 32*pow(f2,2)*pow(Q2h,2)*pow(Sq2,2)\
		))/2. + (pow(l1k2,-1)/M2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-32*Q2e*Q2h*S*Sq2*pow(f2,2) + 16*l2k2*Q2e*Q2h*pow(f2,2)*m2 +
		       32*Q2e*Q2h*Sk*pow(f2,2)*m2 -
		       32*Q2e*Q2h*Sq2*pow(f2,2)*m2 +
		       64*Q2h*S*Sq2*pow(f2,2)*m2 + 8*l2k2*Q2h*pow(f2,2)*pow(Q2e,2) -
		       32*Q2h*S*pow(f2,2)*pow(Q2e,2) + 16*Q2h*Sk*pow(f2,2)*pow(Q2e,2) -
		       16*Q2h*Sq2*pow(f2,2)*pow(Q2e,2) -
		       24*Q2h*pow(f2,2)*m2*pow(Q2e,2) - 4*Q2h*pow(f2,2)*pow(Q2e,3) -
		       8*l2k2*Q2e*pow(f2,2)*pow(Q2h,2) + 96*S*Sq2*pow(f2,2)*pow(Q2h,2) +
		       16*Q2e*pow(f2,2)*m2*pow(Q2h,2) +
		       32*Sk*pow(f2,2)*m2*pow(Q2h,2) +
		       4*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) + 32*S*pow(f2,2)*pow(Q2h,3) -
		       16*Sk*pow(f2,2)*pow(Q2h,3) + 16*Sq2*pow(f2,2)*pow(Q2h,3) -
		       8*pow(f2,2)*pow(Q2h,4) - 32*Q2e*Q2h*pow(f2,2)*pow(S,2) +
		       64*Q2h*pow(f2,2)*m2*pow(S,2) +
		       96*pow(f2,2)*pow(Q2h,2)*pow(S,2) + 32*pow(f2,2)*pow(Q2h,2)*pow(Sq2,2)\
		))/2. + (pow(l1k2,-1)*pow(l2k1,-2)/M2*
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (-32*Q2e*Q2h*S*pow(f2,2)*m4 +
		       16*Q2e*Q2h*Sk*pow(f2,2)*m4 +
		       64*Q2h*S*Sk*pow(f2,2)*m4 -
		       32*Q2e*Q2h*Sq2*pow(f2,2)*m4 -
		       128*Q2h*S*Sq2*pow(f2,2)*m4 +
		       64*Q2h*Sk*Sq2*pow(f2,2)*m4 -
		       16*Q2e*S*pow(f2,2)*m2*pow(Q2h,2) +
		       8*Q2e*Sk*pow(f2,2)*m2*pow(Q2h,2) +
		       32*S*Sk*pow(f2,2)*m2*pow(Q2h,2) -
		       16*Q2e*Sq2*pow(f2,2)*m2*pow(Q2h,2) -
		       64*S*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       32*Sk*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       8*Q2e*pow(f2,2)*m4*pow(Q2h,2) +
		       32*S*pow(f2,2)*m4*pow(Q2h,2) -
		       16*Sk*pow(f2,2)*m4*pow(Q2h,2) +
		       32*Sq2*pow(f2,2)*m4*pow(Q2h,2) +
		       4*Q2e*pow(f2,2)*m2*pow(Q2h,3) +
		       16*S*pow(f2,2)*m2*pow(Q2h,3) -
		       8*Sk*pow(f2,2)*m2*pow(Q2h,3) +
		       16*Sq2*pow(f2,2)*m2*pow(Q2h,3) -
		       64*Q2h*pow(f2,2)*m4*pow(S,2) -
		       32*pow(f2,2)*m2*pow(Q2h,2)*pow(S,2) -
		       16*pow(f2,2)*m2*pow(Q2h,2)*pow(Sk,2) -
		       64*Q2h*pow(f2,2)*m4*pow(Sq2,2) -
		       32*pow(f2,2)*m2*pow(Q2h,2)*pow(Sq2,2)))/2. +
		  (pow(l1k1,-1)*pow(l1k2,-1)*pow(l2k1,-1)/M2*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (32*Q2e*Q2h*S*Sk*pow(f2,2)*m2 -
		       32*Q2e*Q2h*S*Sq2*pow(f2,2)*m2 +
		       48*Q2e*Q2h*Sk*Sq2*pow(f2,2)*m2 +
		       32*Q2e*Q2h*S*pow(f2,2)*m4 -
		       16*Q2e*Q2h*Sk*pow(f2,2)*m4 +
		       64*Q2h*S*Sk*pow(f2,2)*m4 +
		       16*Q2e*Q2h*Sq2*pow(f2,2)*m4 -
		       64*Q2h*S*Sq2*pow(f2,2)*m4 +
		       16*Q2h*S*pow(f2,2)*m2*pow(Q2e,2) -
		       8*Q2h*Sk*pow(f2,2)*m2*pow(Q2e,2) +
		       8*Q2h*Sq2*pow(f2,2)*m2*pow(Q2e,2) -
		       16*Q2e*S*pow(f2,2)*m2*pow(Q2h,2) +
		       8*Q2e*Sk*pow(f2,2)*m2*pow(Q2h,2) -
		       8*Q2e*Sq2*pow(f2,2)*m2*pow(Q2h,2) -
		       16*Sk*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       8*Q2e*pow(f2,2)*m4*pow(Q2h,2) -
		       32*S*pow(f2,2)*m4*pow(Q2h,2) +
		       16*Sk*pow(f2,2)*m4*pow(Q2h,2) -
		       16*Sq2*pow(f2,2)*m4*pow(Q2h,2) +
		       4*pow(f2,2)*m2*pow(Q2e,2)*pow(Q2h,2) -
		       32*Q2e*Q2h*pow(f2,2)*m2*pow(S,2) -
		       64*Q2h*pow(f2,2)*m4*pow(S,2) -
		       24*Q2e*Q2h*pow(f2,2)*m2*pow(Sk,2) +
		       8*pow(f2,2)*m2*pow(Q2h,2)*pow(Sk,2) -
		       24*Q2e*Q2h*pow(f2,2)*m2*pow(Sq2,2) +
		       8*pow(f2,2)*m2*pow(Q2h,2)*pow(Sq2,2)))/2. +
		  (pow(l1k1,-1)*pow(l2k1,-1)*pow(l2k2,-1)/M2*
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (-32*Q2e*Q2h*S*Sk*pow(f2,2)*m2 -
		       32*Q2e*Q2h*S*Sq2*pow(f2,2)*m2 +
		       16*Q2e*Q2h*Sk*Sq2*pow(f2,2)*m2 -
		       32*Q2e*Q2h*S*pow(f2,2)*m4 -
		       16*Q2e*Q2h*Sk*pow(f2,2)*m4 -
		       64*Q2h*S*Sk*pow(f2,2)*m4 -
		       16*Q2e*Q2h*Sq2*pow(f2,2)*m4 -
		       64*Q2h*S*Sq2*pow(f2,2)*m4 -
		       64*Q2h*Sk*Sq2*pow(f2,2)*m4 -
		       16*Q2h*S*pow(f2,2)*m2*pow(Q2e,2) -
		       8*Q2h*Sk*pow(f2,2)*m2*pow(Q2e,2) -
		       8*Q2h*Sq2*pow(f2,2)*m2*pow(Q2e,2) +
		       16*Q2e*S*pow(f2,2)*m2*pow(Q2h,2) +
		       8*Q2e*Sk*pow(f2,2)*m2*pow(Q2h,2) +
		       8*Q2e*Sq2*pow(f2,2)*m2*pow(Q2h,2) -
		       16*Sk*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       8*Q2e*pow(f2,2)*m4*pow(Q2h,2) +
		       32*S*pow(f2,2)*m4*pow(Q2h,2) +
		       16*Sk*pow(f2,2)*m4*pow(Q2h,2) +
		       16*Sq2*pow(f2,2)*m4*pow(Q2h,2) +
		       4*pow(f2,2)*m2*pow(Q2e,2)*pow(Q2h,2) -
		       32*Q2e*Q2h*pow(f2,2)*m2*pow(S,2) -
		       64*Q2h*pow(f2,2)*m4*pow(S,2) -
		       24*Q2e*Q2h*pow(f2,2)*m2*pow(Sk,2) +
		       8*pow(f2,2)*m2*pow(Q2h,2)*pow(Sk,2) -
		       24*Q2e*Q2h*pow(f2,2)*m2*pow(Sq2,2) +
		       8*pow(f2,2)*m2*pow(Q2h,2)*pow(Sq2,2)))/2. +
		  (pow(l1k1,-1)*pow(l2k2,-1)/M2*
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-32*Q2e*Q2h*Sq2*pow(f2,2)*m4 +
		       128*Q2h*S*Sq2*pow(f2,2)*m4 -
		       16*Q2h*pow(f2,2)*m4*pow(Q2e,2) -
		       64*Q2e*S*Sk*pow(f2,2)*pow(Q2h,2) -
		       32*Q2e*S*Sq2*pow(f2,2)*pow(Q2h,2) -
		       16*Q2e*Sk*Sq2*pow(f2,2)*pow(Q2h,2) +
		       16*Q2e*Sk*pow(f2,2)*m2*pow(Q2h,2) +
		       128*S*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       32*Sk*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       16*Q2e*pow(f2,2)*m4*pow(Q2h,2) +
		       32*Sq2*pow(f2,2)*m4*pow(Q2h,2) +
		       16*S*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) +
		       8*Sk*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) +
		       8*Sq2*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) -
		       8*pow(f2,2)*m2*pow(Q2e,2)*pow(Q2h,2) -
		       32*Q2e*S*pow(f2,2)*pow(Q2h,3) - 16*Q2e*Sk*pow(f2,2)*pow(Q2h,3) +
		       64*S*Sk*pow(f2,2)*pow(Q2h,3) - 8*Q2e*Sq2*pow(f2,2)*pow(Q2h,3) +
		       64*S*Sq2*pow(f2,2)*pow(Q2h,3) + 16*Sk*Sq2*pow(f2,2)*pow(Q2h,3) -
		       16*Sk*pow(f2,2)*m2*pow(Q2h,3) -
		       16*pow(f2,2)*m4*pow(Q2h,3) - 4*Q2e*pow(f2,2)*pow(Q2h,4) +
		       16*S*pow(f2,2)*pow(Q2h,4) + 8*Sk*pow(f2,2)*pow(Q2h,4) -
		       8*pow(f2,2)*m2*pow(Q2h,4) +
		       128*Q2h*pow(f2,2)*m4*pow(S,2) -
		       64*Q2e*pow(f2,2)*pow(Q2h,2)*pow(S,2) +
		       128*pow(f2,2)*m2*pow(Q2h,2)*pow(S,2) +
		       96*pow(f2,2)*pow(Q2h,3)*pow(S,2) -
		       16*Q2e*pow(f2,2)*pow(Q2h,2)*pow(Sk,2) -
		       32*pow(f2,2)*m2*pow(Q2h,2)*pow(Sk,2) +
		       16*pow(f2,2)*pow(Q2h,3)*pow(Sk,2) +
		       32*pow(f2,2)*m2*pow(Q2h,2)*pow(Sq2,2) +
		       16*pow(f2,2)*pow(Q2h,3)*pow(Sq2,2)))/2. +
		  (pow(l1k2,-1)*pow(l2k1,-1)/M2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     pow(2*l2k1 + 2*l2k2 - Q2e + Q2h,-1)*
		     (-32*Q2e*Q2h*Sq2*pow(f2,2)*m4 +
		       128*Q2h*S*Sq2*pow(f2,2)*m4 -
		       16*Q2h*pow(f2,2)*m4*pow(Q2e,2) +
		       64*Q2e*S*Sk*pow(f2,2)*pow(Q2h,2) -
		       96*Q2e*S*Sq2*pow(f2,2)*pow(Q2h,2) +
		       48*Q2e*Sk*Sq2*pow(f2,2)*pow(Q2h,2) +
		       16*Q2e*Sk*pow(f2,2)*m2*pow(Q2h,2) +
		       128*S*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       32*Sk*Sq2*pow(f2,2)*m2*pow(Q2h,2) +
		       16*Q2e*pow(f2,2)*m4*pow(Q2h,2) +
		       32*Sq2*pow(f2,2)*m4*pow(Q2h,2) -
		       16*S*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) +
		       8*Sk*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) -
		       8*Sq2*pow(f2,2)*pow(Q2e,2)*pow(Q2h,2) -
		       8*pow(f2,2)*m2*pow(Q2e,2)*pow(Q2h,2) +
		       32*Q2e*S*pow(f2,2)*pow(Q2h,3) - 16*Q2e*Sk*pow(f2,2)*pow(Q2h,3) -
		       64*S*Sk*pow(f2,2)*pow(Q2h,3) + 24*Q2e*Sq2*pow(f2,2)*pow(Q2h,3) +
		       128*S*Sq2*pow(f2,2)*pow(Q2h,3) - 48*Sk*Sq2*pow(f2,2)*pow(Q2h,3) -
		       16*Sk*pow(f2,2)*m2*pow(Q2h,3) -
		       16*pow(f2,2)*m4*pow(Q2h,3) - 4*Q2e*pow(f2,2)*pow(Q2h,4) -
		       16*S*pow(f2,2)*pow(Q2h,4) + 8*Sk*pow(f2,2)*pow(Q2h,4) -
		       16*Sq2*pow(f2,2)*pow(Q2h,4) - 8*pow(f2,2)*m2*pow(Q2h,4) +
		       128*Q2h*pow(f2,2)*m4*pow(S,2) -
		       64*Q2e*pow(f2,2)*pow(Q2h,2)*pow(S,2) +
		       128*pow(f2,2)*m2*pow(Q2h,2)*pow(S,2) +
		       96*pow(f2,2)*pow(Q2h,3)*pow(S,2) -
		       16*Q2e*pow(f2,2)*pow(Q2h,2)*pow(Sk,2) -
		       32*pow(f2,2)*m2*pow(Q2h,2)*pow(Sk,2) +
		       16*pow(f2,2)*pow(Q2h,3)*pow(Sk,2) -
		       32*Q2e*pow(f2,2)*pow(Q2h,2)*pow(Sq2,2) +
		       32*pow(f2,2)*m2*pow(Q2h,2)*pow(Sq2,2) +
		       48*pow(f2,2)*pow(Q2h,3)*pow(Sq2,2)))/2.;
}

long double Melem::melem2_l1k1_1(const long double l1k1,const long double l1k2,const long double k1k2,const long double Q2e,
		const long double Q2h,const long double Q2k,const long double S,const long double Sk,const long double Sq2,
		const long double f1, const long double f2)const {

	long double p1p2 = M2 + Q2h/2.;
	long double G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
	long double G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
	long double G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;

	return G1*(48*k1k2*m2 + 32*l1k1*m2 - 24*l1k2*m2 + 8*k1k2*Q2e - 8*l1k1*Q2e -
		     20*l1k2*Q2e + 16*m2*Q2e - 4*l1k1*Q2h + 28*l1k2*Q2h + 32*m2*Q2h +
		     2*Q2e*Q2h - 40*m2*Q2k + 4*Q2h*Q2k - 4*pow(Q2e,2) - 6*pow(Q2h,2) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (-112*k1k2*l1k2*m2 + 32*k1k2*m4 - 24*k1k2*l1k2*Q2e +
		        8*k1k2*m2*Q2e - 8*l1k2*m2*Q2e + 8*k1k2*l1k2*Q2h -
		        24*k1k2*m2*Q2h - 48*l1k2*m2*Q2h + 32*m4*Q2h - 8*k1k2*Q2e*Q2h -
		        8*l1k2*Q2e*Q2h - 16*k1k2*m2*Q2k + 56*l1k2*m2*Q2k - 32*m4*Q2k -
		        12*k1k2*Q2e*Q2k + 24*l1k2*Q2e*Q2k - 8*m2*Q2e*Q2k +
		        4*k1k2*Q2h*Q2k - 48*l1k2*Q2h*Q2k + 8*m2*Q2h*Q2k + 8*Q2e*Q2h*Q2k +
		        8*Q2e*pow(k1k2,2) + 80*m2*pow(l1k2,2) + 32*Q2e*pow(l1k2,2) -
		        80*Q2h*pow(l1k2,2) + 8*k1k2*pow(Q2e,2) - 4*l1k2*pow(Q2e,2) -
		        6*Q2h*pow(Q2e,2) - 6*Q2k*pow(Q2e,2) + 4*pow(Q2e,3) +
		        52*l1k2*pow(Q2h,2) - 8*m2*pow(Q2h,2) + 2*Q2e*pow(Q2h,2) +
		        6*Q2k*pow(Q2h,2) - 8*pow(Q2h,3) + 8*m2*pow(Q2k,2) +
		        6*Q2e*pow(Q2k,2) - 6*Q2h*pow(Q2k,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*
		      (24*l1k1*l1k2*m2 + 16*l1k1*m4 - 16*m6 - 12*l1k1*l1k2*Q2e +
		        36*l1k1*m2*Q2e + 24*l1k2*m2*Q2e + 28*l1k1*l1k2*Q2h +
		        32*l1k1*m2*Q2h + 8*m4*Q2h - 12*l1k2*Q2e*Q2h + 10*m2*Q2e*Q2h -
		        68*l1k1*m2*Q2k - 16*l1k2*m2*Q2k - 8*l1k1*Q2e*Q2k +
		        14*l1k2*Q2e*Q2k - 20*m2*Q2e*Q2k + 8*l1k1*Q2h*Q2k -
		        30*l1k2*Q2h*Q2k - 16*m2*Q2h*Q2k - 5*Q2e*Q2h*Q2k +
		        72*m2*pow(l1k1,2) + 8*Q2e*pow(l1k1,2) - 8*Q2h*pow(l1k1,2) +
		        40*m2*pow(l1k2,2) + 12*Q2e*pow(l1k2,2) - 36*Q2h*pow(l1k2,2) +
		        4*l1k1*pow(Q2e,2) + 2*l1k2*pow(Q2e,2) + Q2h*pow(Q2e,2) -
		        2*Q2k*pow(Q2e,2) + pow(Q2e,3) - 8*l1k1*pow(Q2h,2) +
		        30*l1k2*pow(Q2h,2) + 10*m2*pow(Q2h,2) + 3*Q2e*pow(Q2h,2) +
		        9*Q2k*pow(Q2h,2) - 7*pow(Q2h,3) + 16*m2*pow(Q2k,2) +
		        3*Q2e*pow(Q2k,2) - 3*Q2h*pow(Q2k,2)) +
		     pow(l1k1,-1)*(32*k1k2*l1k2*m2 + 16*k1k2*l1k2*Q2e + 8*m4*Q2e +
		        16*l1k2*m2*Q2h - 8*m4*Q2h + 4*l1k2*Q2e*Q2h - 16*l1k2*m2*Q2k +
		        16*m4*Q2k + 6*k1k2*Q2e*Q2k - 6*l1k2*Q2e*Q2k - 2*m2*Q2e*Q2k -
		        2*k1k2*Q2h*Q2k + 2*l1k2*Q2h*Q2k - 2*m2*Q2h*Q2k - 4*Q2e*Q2h*Q2k -
		        8*Q2e*pow(k1k2,2) - 32*m2*pow(l1k2,2) - 16*Q2e*pow(l1k2,2) -
		        4*k1k2*pow(Q2e,2) + 4*m2*pow(Q2e,2) + 2*Q2h*pow(Q2e,2) +
		        3*Q2k*pow(Q2e,2) - 3*pow(Q2e,3) - 2*Q2e*pow(Q2h,2) +
		        Q2k*pow(Q2h,2) - 3*Q2e*pow(Q2k,2) + 2*Q2h*pow(Q2k,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		      (-64*k1k2*l1k2*m4 + 64*l1k2*m6 - 32*l1k2*m4*Q2h +
		        32*k1k2*l1k2*m2*Q2k - 16*k1k2*l1k2*Q2h*Q2k + 32*l1k2*m2*Q2h*Q2k -
		        32*l1k2*m2*pow(k1k2,2) + 16*l1k2*Q2h*pow(k1k2,2) +
		        64*k1k2*m2*pow(l1k2,2) - 32*k1k2*Q2h*pow(l1k2,2) +
		        64*m2*Q2h*pow(l1k2,2) - 64*m2*Q2k*pow(l1k2,2) +
		        32*Q2h*Q2k*pow(l1k2,2) - 64*m2*pow(l1k2,3) + 32*Q2h*pow(l1k2,3) +
		        16*k1k2*l1k2*pow(Q2h,2) - 16*l1k2*m2*pow(Q2h,2) -
		        16*l1k2*Q2k*pow(Q2h,2) - 32*pow(l1k2,2)*pow(Q2h,2) +
		        8*l1k2*pow(Q2h,3) - 16*l1k2*m2*pow(Q2k,2) + 8*l1k2*Q2h*pow(Q2k,2)) \
		+ pow(k1k2 - l1k1 - l1k2,-2)*(16*l1k1*l1k2*m4 - 16*l1k1*m6 +
		        8*l1k1*m4*Q2e + 8*l1k2*m4*Q2e - 8*m6*Q2e - 8*l1k1*l1k2*m2*Q2h +
		        16*l1k1*m4*Q2h + 8*l1k2*m4*Q2h - 8*m6*Q2h - 2*l1k1*m2*Q2e*Q2h -
		        4*l1k2*m2*Q2e*Q2h + 8*m4*Q2e*Q2h - 16*l1k1*m4*Q2k -
		        16*l1k2*m4*Q2k + 16*m6*Q2k - 2*l1k1*m2*Q2e*Q2k - 4*m4*Q2e*Q2k +
		        6*l1k1*m2*Q2h*Q2k + 8*l1k2*m2*Q2h*Q2k - 12*m4*Q2h*Q2k +
		        l1k1*Q2e*Q2h*Q2k + 2*m2*Q2e*Q2h*Q2k + 16*m4*pow(l1k1,2) +
		        4*m2*Q2e*pow(l1k1,2) - 4*m2*Q2h*pow(l1k1,2) -
		        2*Q2e*Q2h*pow(l1k1,2) - 8*m2*Q2k*pow(l1k1,2) +
		        4*Q2h*Q2k*pow(l1k1,2) + 8*m2*pow(l1k1,3) - 4*Q2h*pow(l1k1,3) +
		        8*l1k1*m2*pow(l1k2,2) + 4*m2*Q2e*pow(l1k2,2) -
		        4*l1k1*Q2h*pow(l1k2,2) + 4*m2*Q2h*pow(l1k2,2) -
		        2*Q2e*Q2h*pow(l1k2,2) - 8*m2*Q2k*pow(l1k2,2) +
		        4*Q2h*Q2k*pow(l1k2,2) - 4*l1k1*m2*pow(Q2h,2) -
		        4*l1k2*m2*pow(Q2h,2) + 4*m4*pow(Q2h,2) - l1k1*Q2e*pow(Q2h,2) -
		        2*m2*Q2e*pow(Q2h,2) + l1k1*Q2k*pow(Q2h,2) + 2*m2*Q2k*pow(Q2h,2) -
		        2*pow(l1k1,2)*pow(Q2h,2) - 2*pow(l1k2,2)*pow(Q2h,2) +
		        2*l1k1*m2*pow(Q2k,2) + 4*m4*pow(Q2k,2) - l1k1*Q2h*pow(Q2k,2) -
		        2*m2*Q2h*pow(Q2k,2)) +
		     pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-2)*
		      (4*l1k2*m4*Q2e*Q2h - 4*m6*Q2e*Q2h - 4*l1k2*m4*Q2e*Q2k +
		        4*m6*Q2e*Q2k - 4*l1k2*m4*Q2h*Q2k + 4*m6*Q2h*Q2k +
		        2*l1k2*m2*Q2e*Q2h*Q2k - 2*m4*Q2e*Q2h*Q2k +
		        2*m2*Q2e*Q2h*pow(l1k2,2) - 2*m2*Q2e*Q2k*pow(l1k2,2) -
		        2*m2*Q2h*Q2k*pow(l1k2,2) + Q2e*Q2h*Q2k*pow(l1k2,2) -
		        2*l1k2*m2*Q2e*pow(Q2h,2) + 2*m4*Q2e*pow(Q2h,2) +
		        2*l1k2*m2*Q2k*pow(Q2h,2) - 2*m4*Q2k*pow(Q2h,2) -
		        Q2e*pow(l1k2,2)*pow(Q2h,2) + Q2k*pow(l1k2,2)*pow(Q2h,2) +
		        4*l1k2*m4*pow(Q2k,2) - 4*m6*pow(Q2k,2) -
		        2*l1k2*m2*Q2h*pow(Q2k,2) + 2*m4*Q2h*pow(Q2k,2) +
		        2*m2*pow(l1k2,2)*pow(Q2k,2) - Q2h*pow(l1k2,2)*pow(Q2k,2)) +
		     pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (16*k1k2*m4*Q2e - 16*k1k2*m4*Q2h - 16*k1k2*m2*Q2e*Q2h +
		        8*m4*Q2e*Q2h - 12*k1k2*m2*Q2e*Q2k - 16*m4*Q2e*Q2k +
		        12*k1k2*m2*Q2h*Q2k + 16*m4*Q2h*Q2k + 8*k1k2*Q2e*Q2h*Q2k +
		        4*m2*Q2e*Q2h*Q2k - 4*Q2e*Q2h*pow(k1k2,2) + 8*k1k2*m2*pow(Q2e,2) -
		        8*k1k2*Q2h*pow(Q2e,2) + 4*m2*Q2h*pow(Q2e,2) -
		        6*k1k2*Q2k*pow(Q2e,2) - 6*m2*Q2k*pow(Q2e,2) +
		        7*Q2h*Q2k*pow(Q2e,2) + 4*pow(k1k2,2)*pow(Q2e,2) +
		        4*k1k2*pow(Q2e,3) - 4*Q2h*pow(Q2e,3) - 3*Q2k*pow(Q2e,3) +
		        2*pow(Q2e,4) + 8*k1k2*m2*pow(Q2h,2) - 8*m4*pow(Q2h,2) +
		        4*k1k2*Q2e*pow(Q2h,2) - 6*m2*Q2e*pow(Q2h,2) -
		        2*k1k2*Q2k*pow(Q2h,2) + 2*m2*Q2k*pow(Q2h,2) -
		        5*Q2e*Q2k*pow(Q2h,2) + 3*pow(Q2e,2)*pow(Q2h,2) + 2*m2*pow(Q2h,3) -
		        Q2e*pow(Q2h,3) + Q2k*pow(Q2h,3) + 6*m2*Q2e*pow(Q2k,2) -
		        6*m2*Q2h*pow(Q2k,2) - 5*Q2e*Q2h*pow(Q2k,2) +
		        3*pow(Q2e,2)*pow(Q2k,2) + 2*pow(Q2h,2)*pow(Q2k,2)) +
		     pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*
		      (8*l1k2*m4*Q2e - 8*m6*Q2e - 8*l1k2*m4*Q2h + 8*l1k2*m2*Q2e*Q2h +
		        4*m4*Q2e*Q2h + 8*m6*Q2k - 14*l1k2*m2*Q2e*Q2k - 4*m4*Q2e*Q2k -
		        10*l1k2*m2*Q2h*Q2k + 8*m4*Q2h*Q2k - l1k2*Q2e*Q2h*Q2k -
		        2*m2*Q2e*Q2h*Q2k + 4*m2*Q2e*pow(l1k2,2) - 4*Q2e*Q2h*pow(l1k2,2) -
		        4*m2*Q2k*pow(l1k2,2) + 2*Q2e*Q2k*pow(l1k2,2) +
		        2*Q2h*Q2k*pow(l1k2,2) + 2*l1k2*m2*pow(Q2e,2) +
		        l1k2*Q2h*pow(Q2e,2) + 3*m2*Q2h*pow(Q2e,2) - l1k2*Q2k*pow(Q2e,2) -
		        2*m2*Q2k*pow(Q2e,2) - (Q2h*Q2k*pow(Q2e,2))/2. +
		        (Q2k*pow(Q2e,3))/2. + 4*l1k2*m2*pow(Q2h,2) - 4*m4*pow(Q2h,2) +
		        l1k2*Q2e*pow(Q2h,2) - 2*m2*Q2e*pow(Q2h,2) + m2*pow(Q2h,3) +
		        10*l1k2*m2*pow(Q2k,2) - 4*m4*pow(Q2k,2) + 2*m2*Q2e*pow(Q2k,2) +
		        (Q2e*pow(Q2k,3))/2. - (Q2h*pow(Q2k,3))/2.) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (32*l1k2*m4*Q2e - 32*l1k2*m4*Q2h + 16*l1k2*m2*Q2e*Q2h -
		        16*m4*Q2e*Q2h - 16*l1k2*m2*Q2e*Q2k + 16*m4*Q2e*Q2k +
		        48*l1k2*m2*Q2h*Q2k - 16*m4*Q2h*Q2k + 48*l1k2*Q2e*Q2h*Q2k -
		        32*m2*Q2e*pow(l1k2,2) + 96*m2*Q2h*pow(l1k2,2) +
		        80*Q2e*Q2h*pow(l1k2,2) - 64*m2*Q2k*pow(l1k2,2) -
		        32*Q2e*Q2k*pow(l1k2,2) + 64*Q2h*Q2k*pow(l1k2,2) -
		        64*m2*pow(l1k2,3) - 32*Q2e*pow(l1k2,3) + 64*Q2h*pow(l1k2,3) +
		        8*l1k2*m2*pow(Q2e,2) + 12*l1k2*Q2h*pow(Q2e,2) -
		        4*m2*Q2h*pow(Q2e,2) - 8*l1k2*Q2k*pow(Q2e,2) +
		        4*m2*Q2k*pow(Q2e,2) + 2*Q2h*Q2k*pow(Q2e,2) -
		        16*pow(l1k2,2)*pow(Q2e,2) - 4*l1k2*pow(Q2e,3) + 2*Q2h*pow(Q2e,3) -
		        2*Q2k*pow(Q2e,3) - 40*l1k2*m2*pow(Q2h,2) + 16*m4*pow(Q2h,2) -
		        48*l1k2*Q2e*pow(Q2h,2) - 56*l1k2*Q2k*pow(Q2h,2) -
		        4*m2*Q2k*pow(Q2h,2) - 12*Q2e*Q2k*pow(Q2h,2) -
		        96*pow(l1k2,2)*pow(Q2h,2) - 2*pow(Q2e,2)*pow(Q2h,2) +
		        48*l1k2*pow(Q2h,3) + 4*m2*pow(Q2h,3) + 8*Q2e*pow(Q2h,3) +
		        12*Q2k*pow(Q2h,3) - 8*pow(Q2h,4) - 16*l1k2*m2*pow(Q2k,2) -
		        12*l1k2*Q2e*pow(Q2k,2) + 20*l1k2*Q2h*pow(Q2k,2) +
		        6*Q2e*Q2h*pow(Q2k,2) - 6*pow(Q2h,2)*pow(Q2k,2) - 2*Q2e*pow(Q2k,3) +
		        2*Q2h*pow(Q2k,3)) + pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*
		      pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (-16*m6*Q2e*Q2h - 16*m4*Q2e*Q2h*Q2k + 8*m6*pow(Q2e,2) -
		        12*m4*Q2h*pow(Q2e,2) + 4*m4*Q2k*pow(Q2e,2) +
		        m2*Q2h*Q2k*pow(Q2e,2) + 4*m4*pow(Q2e,3) - 4*m2*Q2h*pow(Q2e,3) +
		        m2*Q2k*pow(Q2e,3) + Q2h*Q2k*pow(Q2e,3) - (Q2k*pow(Q2e,4))/2. +
		        8*m6*pow(Q2h,2) + 16*m4*Q2e*pow(Q2h,2) + 12*m4*Q2k*pow(Q2h,2) -
		        m2*Q2e*Q2k*pow(Q2h,2) + 7*m2*pow(Q2e,2)*pow(Q2h,2) -
		        (Q2k*pow(Q2e,2)*pow(Q2h,2))/2. - 8*m4*pow(Q2h,3) -
		        4*m2*Q2e*pow(Q2h,3) - m2*Q2k*pow(Q2h,3) + m2*pow(Q2h,4) +
		        4*m4*Q2e*pow(Q2k,2) - 4*m4*Q2h*pow(Q2k,2) -
		        m2*pow(Q2e,2)*pow(Q2k,2) + m2*pow(Q2h,2)*pow(Q2k,2) +
		        Q2e*Q2h*pow(Q2k,3) - (pow(Q2e,2)*pow(Q2k,3))/2. -
		        (pow(Q2h,2)*pow(Q2k,3))/2.)) +
		  (G2 + G3)*((32*k1k2*m2*M2 - 16*l1k1*m2*M2 + 16*l1k2*m2*M2 +
		        24*k1k2*M2*Q2e - 8*l1k1*M2*Q2e - 8*l1k2*M2*Q2e +
		        8*m2*M2*Q2e - 8*k1k2*M2*Q2h - 16*l1k1*M2*Q2h +
		        24*l1k2*M2*Q2h + 8*m2*M2*Q2h + 8*M2*Q2e*Q2h -
		        16*m2*M2*Q2k - 8*M2*Q2e*Q2k - 24*m2*Q2h*Q2k +
		        16*M2*Q2h*Q2k + 16*l1k1*Q2h*S + 16*l1k2*Q2h*S - 16*Q2e*Q2h*S -
		        32*m2*Q2k*S + 24*Q2e*Q2k*S - 16*Q2h*Q2k*S - 48*m2*Q2h*Sk -
		        64*m2*S*Sk + 48*Q2e*S*Sk - 16*Q2h*S*Sk + 16*l1k1*Q2h*Sq2 +
		        32*m2*Q2h*Sq2 + 8*Q2e*Q2h*Sq2 - 8*Q2h*Q2k*Sq2 + 64*k1k2*S*Sq2 +
		        64*l1k1*S*Sq2 - 128*l1k2*S*Sq2 + 128*m2*S*Sq2 + 64*Q2h*S*Sq2 -
		        64*Q2k*S*Sq2 + 16*l1k2*pow(Q2h,2) + 8*m2*pow(Q2h,2) -
		        24*M2*pow(Q2h,2) + 4*Q2k*pow(Q2h,2) - 8*S*pow(Q2h,2) +
		        8*Sk*pow(Q2h,2) - 16*Sq2*pow(Q2h,2) - 4*pow(Q2h,3) +
		        64*k1k2*pow(S,2) + 96*l1k1*pow(S,2) - 128*l1k2*pow(S,2) +
		        128*m2*pow(S,2) + 16*Q2e*pow(S,2) + 80*Q2h*pow(S,2) -
		        64*Q2k*pow(S,2))/2. + (pow(k1k2 - l1k1 - l1k2,-2)*
		        (-16*l1k1*l1k2*m2*M2*Q2h + 16*l1k1*m4*M2*Q2h -
		          8*l1k1*m2*M2*Q2e*Q2h - 8*l1k2*m2*M2*Q2e*Q2h +
		          8*m4*M2*Q2e*Q2h + 16*l1k1*m2*M2*Q2h*Q2k +
		          16*l1k2*m2*M2*Q2h*Q2k - 16*m4*M2*Q2h*Q2k +
		          2*l1k1*M2*Q2e*Q2h*Q2k + 4*m2*M2*Q2e*Q2h*Q2k +
		          32*l1k1*l1k2*m2*Q2h*S - 32*l1k1*m4*Q2h*S +
		          16*l1k1*m2*Q2e*Q2h*S + 16*l1k2*m2*Q2e*Q2h*S -
		          16*m4*Q2e*Q2h*S - 32*l1k1*m2*Q2h*Q2k*S -
		          32*l1k2*m2*Q2h*Q2k*S + 32*m4*Q2h*Q2k*S - 4*l1k1*Q2e*Q2h*Q2k*S -
		          8*m2*Q2e*Q2h*Q2k*S - 16*m2*M2*Q2h*pow(l1k1,2) -
		          4*M2*Q2e*Q2h*pow(l1k1,2) + 8*M2*Q2h*Q2k*pow(l1k1,2) +
		          32*m2*Q2h*S*pow(l1k1,2) + 8*Q2e*Q2h*S*pow(l1k1,2) -
		          16*Q2h*Q2k*S*pow(l1k1,2) - 8*M2*Q2h*pow(l1k1,3) +
		          16*Q2h*S*pow(l1k1,3) - 8*l1k1*M2*Q2h*pow(l1k2,2) -
		          4*M2*Q2e*Q2h*pow(l1k2,2) + 8*M2*Q2h*Q2k*pow(l1k2,2) +
		          16*l1k1*Q2h*S*pow(l1k2,2) + 8*Q2e*Q2h*S*pow(l1k2,2) -
		          16*Q2h*Q2k*S*pow(l1k2,2) - 8*l1k1*m2*M2*pow(Q2h,2) -
		          8*l1k2*m2*M2*pow(Q2h,2) + 8*m4*M2*pow(Q2h,2) -
		          2*l1k1*M2*Q2e*pow(Q2h,2) - 4*m2*M2*Q2e*pow(Q2h,2) +
		          2*l1k1*M2*Q2k*pow(Q2h,2) + 4*m2*M2*Q2k*pow(Q2h,2) +
		          16*l1k1*m2*S*pow(Q2h,2) + 16*l1k2*m2*S*pow(Q2h,2) -
		          16*m4*S*pow(Q2h,2) + 4*l1k1*Q2e*S*pow(Q2h,2) +
		          8*m2*Q2e*S*pow(Q2h,2) - 4*l1k1*Q2k*S*pow(Q2h,2) -
		          8*m2*Q2k*S*pow(Q2h,2) - 4*M2*pow(l1k1,2)*pow(Q2h,2) +
		          8*S*pow(l1k1,2)*pow(Q2h,2) - 4*M2*pow(l1k2,2)*pow(Q2h,2) +
		          8*S*pow(l1k2,2)*pow(Q2h,2) - 2*l1k1*M2*Q2h*pow(Q2k,2) -
		          4*m2*M2*Q2h*pow(Q2k,2) + 4*l1k1*Q2h*S*pow(Q2k,2) +
		          8*m2*Q2h*S*pow(Q2k,2) + 64*l1k1*l1k2*m2*pow(S,2) -
		          64*l1k1*m4*pow(S,2) + 32*l1k1*m2*Q2e*pow(S,2) +
		          32*l1k2*m2*Q2e*pow(S,2) - 32*m4*Q2e*pow(S,2) +
		          32*l1k1*m2*Q2h*pow(S,2) + 32*l1k2*m2*Q2h*pow(S,2) -
		          32*m4*Q2h*pow(S,2) + 8*l1k1*Q2e*Q2h*pow(S,2) +
		          16*m2*Q2e*Q2h*pow(S,2) - 64*l1k1*m2*Q2k*pow(S,2) -
		          64*l1k2*m2*Q2k*pow(S,2) + 64*m4*Q2k*pow(S,2) -
		          8*l1k1*Q2e*Q2k*pow(S,2) - 16*m2*Q2e*Q2k*pow(S,2) -
		          8*l1k1*Q2h*Q2k*pow(S,2) - 16*m2*Q2h*Q2k*pow(S,2) +
		          64*m2*pow(l1k1,2)*pow(S,2) + 16*Q2e*pow(l1k1,2)*pow(S,2) +
		          16*Q2h*pow(l1k1,2)*pow(S,2) - 32*Q2k*pow(l1k1,2)*pow(S,2) +
		          32*pow(l1k1,3)*pow(S,2) + 32*l1k1*pow(l1k2,2)*pow(S,2) +
		          16*Q2e*pow(l1k2,2)*pow(S,2) + 16*Q2h*pow(l1k2,2)*pow(S,2) -
		          32*Q2k*pow(l1k2,2)*pow(S,2) + 8*l1k1*pow(Q2k,2)*pow(S,2) +
		          16*m2*pow(Q2k,2)*pow(S,2)))/2. +
		     (pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-2)*
		        (4*l1k2*m2*M2*Q2e*Q2h*Q2k - 4*m4*M2*Q2e*Q2h*Q2k -
		          8*l1k2*m2*Q2e*Q2h*Q2k*S + 8*m4*Q2e*Q2h*Q2k*S +
		          2*M2*Q2e*Q2h*Q2k*pow(l1k2,2) - 4*Q2e*Q2h*Q2k*S*pow(l1k2,2) -
		          4*l1k2*m2*M2*Q2e*pow(Q2h,2) + 4*m4*M2*Q2e*pow(Q2h,2) +
		          4*l1k2*m2*M2*Q2k*pow(Q2h,2) - 4*m4*M2*Q2k*pow(Q2h,2) +
		          8*l1k2*m2*Q2e*S*pow(Q2h,2) - 8*m4*Q2e*S*pow(Q2h,2) -
		          8*l1k2*m2*Q2k*S*pow(Q2h,2) + 8*m4*Q2k*S*pow(Q2h,2) -
		          2*M2*Q2e*pow(l1k2,2)*pow(Q2h,2) +
		          2*M2*Q2k*pow(l1k2,2)*pow(Q2h,2) +
		          4*Q2e*S*pow(l1k2,2)*pow(Q2h,2) - 4*Q2k*S*pow(l1k2,2)*pow(Q2h,2) -
		          4*l1k2*m2*M2*Q2h*pow(Q2k,2) + 4*m4*M2*Q2h*pow(Q2k,2) +
		          8*l1k2*m2*Q2h*S*pow(Q2k,2) - 8*m4*Q2h*S*pow(Q2k,2) -
		          2*M2*Q2h*pow(l1k2,2)*pow(Q2k,2) +
		          4*Q2h*S*pow(l1k2,2)*pow(Q2k,2) + 16*l1k2*m2*Q2e*Q2h*pow(S,2) -
		          16*m4*Q2e*Q2h*pow(S,2) - 16*l1k2*m2*Q2e*Q2k*pow(S,2) +
		          16*m4*Q2e*Q2k*pow(S,2) - 16*l1k2*m2*Q2h*Q2k*pow(S,2) +
		          16*m4*Q2h*Q2k*pow(S,2) + 8*Q2e*Q2h*pow(l1k2,2)*pow(S,2) -
		          8*Q2e*Q2k*pow(l1k2,2)*pow(S,2) - 8*Q2h*Q2k*pow(l1k2,2)*pow(S,2) +
		          16*l1k2*m2*pow(Q2k,2)*pow(S,2) - 16*m4*pow(Q2k,2)*pow(S,2) +
		          8*pow(l1k2,2)*pow(Q2k,2)*pow(S,2)))/2. +
		     (pow(k1k2 - l1k1 - l1k2,-1)*
		        (48*l1k1*l1k2*m2*M2 + 32*l1k2*m4*M2 + 16*l1k1*l1k2*M2*Q2e +
		          40*l1k1*m2*M2*Q2e + 32*l1k2*m2*M2*Q2e + 16*l1k2*m4*Q2h +
		          16*l1k1*l1k2*M2*Q2h + 32*l1k2*m2*M2*Q2h + 16*m4*M2*Q2h +
		          4*l1k1*M2*Q2e*Q2h + 12*l1k2*M2*Q2e*Q2h +
		          12*m2*M2*Q2e*Q2h - 56*l1k1*m2*M2*Q2k -
		          32*l1k2*m2*M2*Q2k - 24*l1k1*M2*Q2e*Q2k -
		          12*l1k2*M2*Q2e*Q2k - 20*m2*M2*Q2e*Q2k -
		          24*l1k1*m2*Q2h*Q2k - 16*l1k2*m2*Q2h*Q2k + 16*m4*Q2h*Q2k +
		          24*l1k1*M2*Q2h*Q2k - 20*l1k2*M2*Q2h*Q2k -
		          12*m2*M2*Q2h*Q2k + 4*l1k2*Q2e*Q2h*Q2k - 10*m2*Q2e*Q2h*Q2k -
		          16*M2*Q2e*Q2h*Q2k + 16*l1k1*l1k2*Q2h*S + 96*l1k1*m2*Q2h*S +
		          32*l1k2*m2*Q2h*S - 64*m4*Q2h*S + 24*l1k1*Q2e*Q2h*S -
		          8*l1k2*Q2e*Q2h*S + 40*m2*Q2e*Q2h*S - 32*l1k1*m2*Q2k*S -
		          64*l1k2*m2*Q2k*S + 64*m4*Q2k*S + 24*l1k1*Q2e*Q2k*S +
		          24*l1k2*Q2e*Q2k*S + 16*m2*Q2e*Q2k*S - 32*l1k1*Q2h*Q2k*S -
		          24*l1k2*Q2h*Q2k*S - 64*m2*Q2h*Q2k*S + 4*Q2e*Q2h*Q2k*S -
		          32*l1k2*m4*Sk - 48*l1k1*m2*Q2h*Sk - 32*l1k2*m2*Q2h*Sk +
		          32*m4*Q2h*Sk + 8*l1k2*Q2e*Q2h*Sk - 20*m2*Q2e*Q2h*Sk -
		          16*l1k2*m2*Q2k*Sk + 24*m2*Q2e*Q2k*Sk - 4*Q2e*Q2h*Q2k*Sk -
		          64*l1k1*m2*S*Sk - 128*l1k2*m2*S*Sk + 128*m4*S*Sk +
		          48*l1k1*Q2e*S*Sk + 48*l1k2*Q2e*S*Sk + 32*m2*Q2e*S*Sk -
		          16*l1k1*Q2h*S*Sk - 16*l1k2*Q2h*S*Sk + 40*Q2e*Q2h*S*Sk +
		          32*m2*Q2k*S*Sk - 24*Q2e*Q2k*S*Sk + 8*Q2h*Q2k*S*Sk -
		          16*l1k1*l1k2*Q2h*Sq2 + 64*l1k1*m2*Q2h*Sq2 +
		          16*l1k2*m2*Q2h*Sq2 - 16*m4*Q2h*Sq2 + 24*l1k1*Q2e*Q2h*Sq2 +
		          16*m2*Q2e*Q2h*Sq2 - 16*l1k1*Q2h*Q2k*Sq2 + 16*l1k2*Q2h*Q2k*Sq2 -
		          32*m2*Q2h*Q2k*Sq2 - 8*Q2e*Q2h*Q2k*Sq2 - 64*l1k1*l1k2*S*Sq2 +
		          128*l1k1*m2*S*Sq2 - 64*m4*S*Sq2 + 32*l1k2*Q2e*S*Sq2 +
		          64*l1k1*Q2h*S*Sq2 - 128*l1k2*Q2h*S*Sq2 + 64*m2*Q2h*S*Sq2 -
		          56*Q2e*Q2h*S*Sq2 - 64*l1k1*Q2k*S*Sq2 + 64*l1k2*Q2k*S*Sq2 -
		          96*m2*Q2k*S*Sq2 - 8*Q2e*Q2k*S*Sq2 - 24*Q2h*Q2k*S*Sq2 +
		          32*m2*Q2h*Sk*Sq2 + 24*Q2e*Q2h*Sk*Sq2 + 48*m2*M2*pow(l1k1,2) +
		          24*M2*Q2e*pow(l1k1,2) - 24*M2*Q2h*pow(l1k1,2) +
		          32*Q2h*S*pow(l1k1,2) + 16*Q2h*Sq2*pow(l1k1,2) +
		          64*S*Sq2*pow(l1k1,2) + 16*m2*M2*pow(l1k2,2) -
		          8*M2*Q2e*pow(l1k2,2) - 40*M2*Q2h*pow(l1k2,2) +
		          16*Q2h*S*pow(l1k2,2) + 32*Q2h*Sq2*pow(l1k2,2) +
		          128*S*Sq2*pow(l1k2,2) + 16*l1k1*M2*pow(Q2e,2) +
		          8*m2*M2*pow(Q2e,2) + 6*M2*Q2h*pow(Q2e,2) -
		          6*M2*Q2k*pow(Q2e,2) + 2*Q2h*Q2k*pow(Q2e,2) -
		          4*Q2h*S*pow(Q2e,2) + 8*Q2k*S*pow(Q2e,2) + 4*Q2h*Sk*pow(Q2e,2) +
		          16*S*Sk*pow(Q2e,2) + 4*Q2h*Sq2*pow(Q2e,2) + 2*M2*pow(Q2e,3) +
		          16*l1k1*l1k2*pow(Q2h,2) + 8*l1k1*m2*pow(Q2h,2) +
		          16*l1k2*m2*pow(Q2h,2) - 8*m4*pow(Q2h,2) -
		          28*l1k1*M2*pow(Q2h,2) + 28*l1k2*M2*pow(Q2h,2) -
		          4*m2*M2*pow(Q2h,2) + 8*m2*Q2e*pow(Q2h,2) +
		          10*M2*Q2e*pow(Q2h,2) + 4*l1k1*Q2k*pow(Q2h,2) -
		          16*l1k2*Q2k*pow(Q2h,2) - 6*m2*Q2k*pow(Q2h,2) +
		          26*M2*Q2k*pow(Q2h,2) - 4*Q2e*Q2k*pow(Q2h,2) +
		          8*l1k2*S*pow(Q2h,2) + 8*m2*S*pow(Q2h,2) + 12*Q2e*S*pow(Q2h,2) -
		          12*Q2k*S*pow(Q2h,2) + 8*l1k1*Sk*pow(Q2h,2) -
		          4*m2*Sk*pow(Q2h,2) - 4*Q2e*Sk*pow(Q2h,2) - 40*S*Sk*pow(Q2h,2) -
		          8*l1k1*Sq2*pow(Q2h,2) - 24*l1k2*Sq2*pow(Q2h,2) +
		          24*m2*Sq2*pow(Q2h,2) + 12*Q2e*Sq2*pow(Q2h,2) +
		          72*S*Sq2*pow(Q2h,2) - 24*Sk*Sq2*pow(Q2h,2) -
		          16*pow(l1k2,2)*pow(Q2h,2) - 2*pow(Q2e,2)*pow(Q2h,2) -
		          4*l1k1*pow(Q2h,3) + 16*l1k2*pow(Q2h,3) - 4*m2*pow(Q2h,3) -
		          22*M2*pow(Q2h,3) + 4*Q2e*pow(Q2h,3) + 6*Q2k*pow(Q2h,3) -
		          4*S*pow(Q2h,3) + 4*Sk*pow(Q2h,3) - 12*Sq2*pow(Q2h,3) -
		          4*pow(Q2h,4) + 16*m2*M2*pow(Q2k,2) + 6*m2*Q2e*pow(Q2k,2) +
		          8*M2*Q2e*pow(Q2k,2) + 6*m2*Q2h*pow(Q2k,2) -
		          8*M2*Q2h*pow(Q2k,2) + 16*m2*S*pow(Q2k,2) -
		          12*Q2e*S*pow(Q2k,2) + 8*Q2h*S*pow(Q2k,2) + 4*Q2h*Sq2*pow(Q2k,2) +
		          16*S*Sq2*pow(Q2k,2) - 2*pow(Q2h,2)*pow(Q2k,2) -
		          64*l1k1*l1k2*pow(S,2) + 192*l1k1*m2*pow(S,2) -
		          64*l1k2*m2*pow(S,2) - 64*m4*pow(S,2) + 32*l1k1*Q2e*pow(S,2) +
		          32*l1k2*Q2e*pow(S,2) + 64*m2*Q2e*pow(S,2) +
		          96*l1k1*Q2h*pow(S,2) - 128*l1k2*Q2h*pow(S,2) +
		          64*m2*Q2h*pow(S,2) - 32*Q2e*Q2h*pow(S,2) -
		          96*l1k1*Q2k*pow(S,2) + 64*l1k2*Q2k*pow(S,2) -
		          128*m2*Q2k*pow(S,2) - 16*Q2e*Q2k*pow(S,2) -
		          32*Q2h*Q2k*pow(S,2) + 128*pow(l1k1,2)*pow(S,2) +
		          128*pow(l1k2,2)*pow(S,2) + 64*pow(Q2h,2)*pow(S,2) +
		          16*pow(Q2k,2)*pow(S,2) - 32*l1k2*m2*pow(Sk,2) +
		          24*m2*Q2e*pow(Sk,2) - 24*m2*Q2h*pow(Sk,2) -
		          8*Q2e*Q2h*pow(Sk,2) + 8*pow(Q2h,2)*pow(Sk,2) -
		          16*m2*Q2h*pow(Sq2,2) - 24*Q2e*Q2h*pow(Sq2,2) +
		          24*pow(Q2h,2)*pow(Sq2,2)))/2. +
		     (pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		        (-32*k1k2*l1k2*m2*M2 - 64*k1k2*m4*M2 + 64*l1k2*m4*M2 -
		          64*k1k2*l1k2*M2*Q2e - 16*k1k2*m2*M2*Q2e -
		          16*l1k2*m2*M2*Q2e - 32*k1k2*m4*Q2h + 32*l1k2*m4*Q2h +
		          32*k1k2*l1k2*M2*Q2h - 80*k1k2*m2*M2*Q2h +
		          64*l1k2*m2*M2*Q2h - 32*m4*M2*Q2h - 24*k1k2*M2*Q2e*Q2h +
		          16*l1k2*M2*Q2e*Q2h + 16*l1k2*m2*M2*Q2k + 32*m4*M2*Q2k -
		          24*k1k2*M2*Q2e*Q2k + 16*l1k2*M2*Q2e*Q2k +
		          32*k1k2*m2*Q2h*Q2k + 16*l1k2*m2*Q2h*Q2k - 16*m4*Q2h*Q2k +
		          8*k1k2*M2*Q2h*Q2k - 64*l1k2*M2*Q2h*Q2k +
		          32*m2*M2*Q2h*Q2k - 8*k1k2*Q2e*Q2h*Q2k + 8*l1k2*Q2e*Q2h*Q2k -
		          8*m2*Q2e*Q2h*Q2k + 24*M2*Q2e*Q2h*Q2k + 32*k1k2*l1k2*Q2h*S +
		          64*l1k2*m2*Q2h*S - 64*m4*Q2h*S - 32*k1k2*Q2e*Q2h*S +
		          48*l1k2*Q2e*Q2h*S - 32*m2*Q2e*Q2h*S - 128*k1k2*m2*Q2k*S +
		          64*l1k2*m2*Q2k*S + 128*m4*Q2k*S + 48*k1k2*Q2e*Q2k*S +
		          32*m2*Q2e*Q2k*S - 16*k1k2*Q2h*Q2k*S - 32*l1k2*Q2h*Q2k*S +
		          32*m2*Q2h*Q2k*S + 40*Q2e*Q2h*Q2k*S + 64*k1k2*m4*Sk -
		          64*l1k2*m4*Sk + 64*k1k2*m2*Q2h*Sk + 32*l1k2*m2*Q2h*Sk -
		          32*m4*Q2h*Sk - 16*k1k2*Q2e*Q2h*Sk + 16*l1k2*Q2e*Q2h*Sk -
		          16*m2*Q2e*Q2h*Sk + 32*k1k2*m2*Q2k*Sk - 32*l1k2*m2*Q2k*Sk -
		          32*m4*Q2k*Sk + 32*m2*Q2h*Q2k*Sk - 256*k1k2*m2*S*Sk +
		          128*l1k2*m2*S*Sk + 256*m4*S*Sk + 96*k1k2*Q2e*S*Sk +
		          64*m2*Q2e*S*Sk - 32*k1k2*Q2h*S*Sk + 64*Q2e*Q2h*S*Sk +
		          64*m2*Q2k*S*Sk - 32*k1k2*m2*Q2h*Sq2 - 32*l1k2*m2*Q2h*Sq2 -
		          32*m4*Q2h*Sq2 - 16*k1k2*Q2e*Q2h*Sq2 - 32*m2*Q2e*Q2h*Sq2 -
		          128*k1k2*m2*Q2k*Sq2 + 64*l1k2*m2*Q2k*Sq2 + 128*m4*Q2k*Sq2 +
		          48*k1k2*Q2e*Q2k*Sq2 + 32*m2*Q2e*Q2k*Sq2 - 16*k1k2*Q2h*Q2k*Sq2 -
		          16*m2*Q2h*Q2k*Sq2 + 16*Q2e*Q2h*Q2k*Sq2 - 384*k1k2*l1k2*S*Sq2 +
		          512*k1k2*m2*S*Sq2 - 512*l1k2*m2*S*Sq2 - 128*m4*S*Sq2 +
		          64*l1k2*Q2e*S*Sq2 + 192*k1k2*Q2h*S*Sq2 - 512*l1k2*Q2h*S*Sq2 +
		          256*m2*Q2h*S*Sq2 - 80*Q2e*Q2h*S*Sq2 - 128*k1k2*Q2k*S*Sq2 +
		          256*l1k2*Q2k*S*Sq2 - 320*m2*Q2k*S*Sq2 + 16*Q2e*Q2k*S*Sq2 -
		          112*Q2h*Q2k*S*Sq2 - 256*k1k2*m2*Sk*Sq2 + 128*l1k2*m2*Sk*Sq2 +
		          256*m4*Sk*Sq2 + 96*k1k2*Q2e*Sk*Sq2 + 64*m2*Q2e*Sk*Sq2 -
		          32*k1k2*Q2h*Sk*Sq2 - 64*m2*Q2h*Sk*Sq2 + 48*Q2e*Q2h*Sk*Sq2 +
		          64*m2*Q2k*Sk*Sq2 - 32*m2*M2*pow(k1k2,2) +
		          16*M2*Q2e*pow(k1k2,2) + 128*S*Sq2*pow(k1k2,2) -
		          32*m2*M2*pow(l1k2,2) - 96*M2*Q2h*pow(l1k2,2) -
		          64*Q2h*S*pow(l1k2,2) + 512*S*Sq2*pow(l1k2,2) +
		          32*k1k2*M2*pow(Q2e,2) - 32*l1k2*M2*pow(Q2e,2) -
		          8*M2*Q2h*pow(Q2e,2) - 28*M2*Q2k*pow(Q2e,2) -
		          8*Q2h*Q2k*pow(Q2e,2) - 24*Q2h*S*pow(Q2e,2) + 8*Q2k*S*pow(Q2e,2) -
		          16*Q2h*Sk*pow(Q2e,2) + 16*S*Sk*pow(Q2e,2) +
		          8*Q2k*Sq2*pow(Q2e,2) + 16*Sk*Sq2*pow(Q2e,2) +
		          12*M2*pow(Q2e,3) - 32*k1k2*m2*pow(Q2h,2) +
		          16*l1k2*m2*pow(Q2h,2) - 8*k1k2*M2*pow(Q2h,2) +
		          96*l1k2*M2*pow(Q2h,2) - 48*m2*M2*pow(Q2h,2) +
		          8*m2*Q2e*pow(Q2h,2) + 8*M2*Q2e*pow(Q2h,2) -
		          24*l1k2*Q2k*pow(Q2h,2) + 8*m2*Q2k*pow(Q2h,2) +
		          20*M2*Q2k*pow(Q2h,2) + 8*Q2e*Q2k*pow(Q2h,2) +
		          16*l1k2*S*pow(Q2h,2) - 32*m2*S*pow(Q2h,2) + 8*Q2e*S*pow(Q2h,2) -
		          32*Q2k*S*pow(Q2h,2) - 16*l1k2*Sk*pow(Q2h,2) +
		          16*Q2e*Sk*pow(Q2h,2) - 80*S*Sk*pow(Q2h,2) +
		          16*l1k2*Sq2*pow(Q2h,2) - 16*Q2k*Sq2*pow(Q2h,2) +
		          176*S*Sq2*pow(Q2h,2) - 64*Sk*Sq2*pow(Q2h,2) -
		          32*pow(l1k2,2)*pow(Q2h,2) + 24*l1k2*pow(Q2h,3) -
		          16*m2*pow(Q2h,3) - 28*M2*pow(Q2h,3) + 4*Q2k*pow(Q2h,3) -
		          8*Sq2*pow(Q2h,3) - 4*pow(Q2h,4) + 16*m2*M2*pow(Q2k,2) +
		          12*M2*Q2e*pow(Q2k,2) + 8*m2*Q2h*pow(Q2k,2) -
		          12*M2*Q2h*pow(Q2k,2) + 4*Q2e*Q2h*pow(Q2k,2) +
		          32*m2*S*pow(Q2k,2) - 16*m2*Sk*pow(Q2k,2) +
		          32*m2*Sq2*pow(Q2k,2) - 4*pow(Q2h,2)*pow(Q2k,2) -
		          256*k1k2*l1k2*pow(S,2) + 256*k1k2*m2*pow(S,2) -
		          384*l1k2*m2*pow(S,2) + 64*l1k2*Q2e*pow(S,2) +
		          128*k1k2*Q2h*pow(S,2) - 384*l1k2*Q2h*pow(S,2) +
		          192*m2*Q2h*pow(S,2) - 64*Q2e*Q2h*pow(S,2) -
		          64*k1k2*Q2k*pow(S,2) + 192*l1k2*Q2k*pow(S,2) -
		          192*m2*Q2k*pow(S,2) + 32*Q2e*Q2k*pow(S,2) -
		          96*Q2h*Q2k*pow(S,2) + 64*pow(k1k2,2)*pow(S,2) +
		          384*pow(l1k2,2)*pow(S,2) + 128*pow(Q2h,2)*pow(S,2) +
		          64*k1k2*m2*pow(Sk,2) - 64*l1k2*m2*pow(Sk,2) +
		          32*m2*Q2h*pow(Sk,2) - 16*Q2e*Q2h*pow(Sk,2) -
		          32*m2*Q2k*pow(Sk,2) + 16*pow(Q2h,2)*pow(Sk,2) -
		          128*k1k2*l1k2*pow(Sq2,2) + 256*k1k2*m2*pow(Sq2,2) -
		          128*l1k2*m2*pow(Sq2,2) - 128*m4*pow(Sq2,2) +
		          64*k1k2*Q2h*pow(Sq2,2) - 128*l1k2*Q2h*pow(Sq2,2) +
		          96*m2*Q2h*pow(Sq2,2) - 32*Q2e*Q2h*pow(Sq2,2) -
		          64*k1k2*Q2k*pow(Sq2,2) + 64*l1k2*Q2k*pow(Sq2,2) -
		          128*m2*Q2k*pow(Sq2,2) - 16*Q2e*Q2k*pow(Sq2,2) -
		          16*Q2h*Q2k*pow(Sq2,2) + 64*pow(k1k2,2)*pow(Sq2,2) +
		          128*pow(l1k2,2)*pow(Sq2,2) + 64*pow(Q2h,2)*pow(Sq2,2)))/2. +
		     (pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*
		        (16*l1k2*m4*M2*Q2e + 8*l1k2*m4*Q2e*Q2h +
		          8*l1k2*m2*M2*Q2e*Q2h + 8*m4*M2*Q2e*Q2h -
		          16*l1k2*m4*M2*Q2k - 8*l1k2*m4*Q2h*Q2k -
		          24*l1k2*m2*M2*Q2h*Q2k - 8*m4*M2*Q2h*Q2k -
		          10*l1k2*m2*Q2e*Q2h*Q2k + 4*m4*Q2e*Q2h*Q2k -
		          4*l1k2*M2*Q2e*Q2h*Q2k + 8*m2*M2*Q2e*Q2h*Q2k +
		          24*l1k2*m2*Q2e*Q2h*S - 16*m4*Q2e*Q2h*S -
		          16*l1k2*m2*Q2e*Q2k*S + 16*m4*Q2e*Q2k*S -
		          32*l1k2*m2*Q2h*Q2k*S + 48*m4*Q2h*Q2k*S + 8*l1k2*Q2e*Q2h*Q2k*S +
		          8*m2*Q2e*Q2h*Q2k*S - 16*l1k2*m4*Q2e*Sk -
		          20*l1k2*m2*Q2e*Q2h*Sk + 8*m4*Q2e*Q2h*Sk + 16*l1k2*m4*Q2k*Sk +
		          16*l1k2*m2*Q2e*Q2k*Sk - 8*l1k2*m2*Q2h*Q2k*Sk -
		          16*m4*Q2h*Q2k*Sk - 4*l1k2*Q2e*Q2h*Q2k*Sk -
		          8*m2*Q2e*Q2h*Q2k*Sk - 32*l1k2*m2*Q2e*S*Sk + 32*m4*Q2e*S*Sk -
		          32*l1k2*m2*Q2h*S*Sk + 32*m4*Q2h*S*Sk + 8*l1k2*Q2e*Q2h*S*Sk +
		          16*m2*Q2e*Q2h*S*Sk + 64*l1k2*m2*Q2k*S*Sk - 64*m4*Q2k*S*Sk -
		          24*l1k2*Q2e*Q2k*S*Sk - 32*m2*Q2e*Q2k*S*Sk +
		          8*l1k2*Q2h*Q2k*S*Sk - 16*Q2e*Q2h*Q2k*S*Sk +
		          8*l1k2*m2*Q2e*Q2h*Sq2 - 8*m4*Q2e*Q2h*Sq2 -
		          8*l1k2*m2*Q2h*Q2k*Sq2 + 8*m4*Q2h*Q2k*Sq2 +
		          10*m2*Q2e*Q2h*Q2k*Sq2 - 32*m4*Q2e*S*Sq2 + 32*m4*Q2k*S*Sq2 -
		          16*l1k2*Q2e*Q2k*S*Sq2 - 32*m2*Q2e*Q2k*S*Sq2 +
		          16*m2*Q2h*Q2k*S*Sq2 - 24*Q2e*Q2h*Q2k*S*Sq2 +
		          16*m2*Q2e*Q2h*Sk*Sq2 - 24*m2*Q2e*Q2k*Sk*Sq2 +
		          8*m2*Q2h*Q2k*Sk*Sq2 - 4*Q2e*Q2h*Q2k*Sk*Sq2 +
		          8*m2*M2*Q2e*pow(l1k2,2) - 8*M2*Q2e*Q2h*pow(l1k2,2) -
		          8*m2*M2*Q2k*pow(l1k2,2) + 4*M2*Q2e*Q2k*pow(l1k2,2) +
		          4*M2*Q2h*Q2k*pow(l1k2,2) + 8*Q2e*Q2h*S*pow(l1k2,2) -
		          8*Q2h*Q2k*S*pow(l1k2,2) + 4*l1k2*m2*M2*pow(Q2e,2) +
		          2*l1k2*M2*Q2h*pow(Q2e,2) - 2*m2*M2*Q2h*pow(Q2e,2) -
		          2*m2*M2*Q2k*pow(Q2e,2) - 4*m2*Q2h*Q2k*pow(Q2e,2) +
		          M2*Q2h*Q2k*pow(Q2e,2) - 4*l1k2*Q2h*S*pow(Q2e,2) +
		          4*l1k2*Q2k*S*pow(Q2e,2) + 8*m2*Q2k*S*pow(Q2e,2) +
		          4*Q2h*Q2k*S*pow(Q2e,2) - 8*m2*Q2h*Sk*pow(Q2e,2) +
		          12*m2*Q2k*Sk*pow(Q2e,2) + 2*Q2h*Q2k*Sk*pow(Q2e,2) +
		          8*l1k2*S*Sk*pow(Q2e,2) + 16*m2*S*Sk*pow(Q2e,2) +
		          16*Q2h*S*Sk*pow(Q2e,2) + 2*Q2h*Q2k*Sq2*pow(Q2e,2) +
		          4*Q2h*S*Sq2*pow(Q2e,2) + 4*Q2k*S*Sq2*pow(Q2e,2) +
		          4*Q2h*Sk*Sq2*pow(Q2e,2) + 8*l1k2*m2*M2*pow(Q2h,2) +
		          8*l1k2*m2*Q2e*pow(Q2h,2) + 2*l1k2*M2*Q2e*pow(Q2h,2) +
		          4*m2*M2*Q2e*pow(Q2h,2) - 6*l1k2*m2*Q2k*pow(Q2h,2) +
		          8*m4*Q2k*pow(Q2h,2) - 14*m2*M2*Q2k*pow(Q2h,2) +
		          2*l1k2*Q2e*Q2k*pow(Q2h,2) + 5*m2*Q2e*Q2k*pow(Q2h,2) -
		          M2*Q2e*Q2k*pow(Q2h,2) - 8*l1k2*m2*S*pow(Q2h,2) -
		          16*m4*S*pow(Q2h,2) - 4*l1k2*Q2k*S*pow(Q2h,2) -
		          8*m2*Q2k*S*pow(Q2h,2) - 4*Q2e*Q2k*S*pow(Q2h,2) +
		          4*l1k2*m2*Sk*pow(Q2h,2) + 8*m4*Sk*pow(Q2h,2) +
		          4*l1k2*Q2e*Sk*pow(Q2h,2) + 8*m2*Q2e*Sk*pow(Q2h,2) -
		          4*m2*Q2k*Sk*pow(Q2h,2) - 2*Q2e*Q2k*Sk*pow(Q2h,2) -
		          16*Q2e*S*Sk*pow(Q2h,2) + 16*Q2k*S*Sk*pow(Q2h,2) +
		          2*m2*Q2e*Sq2*pow(Q2h,2) - 4*l1k2*Q2k*Sq2*pow(Q2h,2) -
		          14*m2*Q2k*Sq2*pow(Q2h,2) - 2*Q2e*Q2k*Sq2*pow(Q2h,2) -
		          4*Q2e*S*Sq2*pow(Q2h,2) + 20*Q2k*S*Sq2*pow(Q2h,2) -
		          4*Q2e*Sk*Sq2*pow(Q2h,2) + 4*Q2k*Sk*Sq2*pow(Q2h,2) +
		          m2*pow(Q2e,2)*pow(Q2h,2) - 4*m4*pow(Q2h,3) +
		          2*m2*M2*pow(Q2h,3) - m2*Q2e*pow(Q2h,3) - m2*Q2k*pow(Q2h,3) +
		          2*m2*Sq2*pow(Q2h,3) + 4*l1k2*m2*M2*pow(Q2k,2) +
		          6*l1k2*m2*Q2e*pow(Q2k,2) + 2*l1k2*M2*Q2e*pow(Q2k,2) +
		          2*l1k2*m2*Q2h*pow(Q2k,2) - 8*m4*Q2h*pow(Q2k,2) -
		          2*l1k2*M2*Q2h*pow(Q2k,2) + 4*m2*M2*Q2h*pow(Q2k,2) -
		          2*l1k2*Q2e*Q2h*pow(Q2k,2) - 2*m2*Q2e*Q2h*pow(Q2k,2) +
		          2*M2*Q2e*Q2h*pow(Q2k,2) + 32*l1k2*m2*S*pow(Q2k,2) -
		          32*m4*S*pow(Q2k,2) - 12*l1k2*Q2e*S*pow(Q2k,2) -
		          16*m2*Q2e*S*pow(Q2k,2) + 8*l1k2*Q2h*S*pow(Q2k,2) +
		          8*m2*Q2h*S*pow(Q2k,2) - 4*Q2e*Q2h*S*pow(Q2k,2) +
		          8*l1k2*m2*Sk*pow(Q2k,2) - 2*Q2e*Q2h*Sk*pow(Q2k,2) -
		          12*m2*Q2e*Sq2*pow(Q2k,2) + 4*l1k2*Q2h*Sq2*pow(Q2k,2) +
		          12*m2*Q2h*Sq2*pow(Q2k,2) + 16*l1k2*S*Sq2*pow(Q2k,2) +
		          16*m2*S*Sq2*pow(Q2k,2) + 4*Q2e*S*Sq2*pow(Q2k,2) -
		          4*Q2h*S*Sq2*pow(Q2k,2) + 3*m2*pow(Q2e,2)*pow(Q2k,2) -
		          M2*pow(Q2e,2)*pow(Q2k,2) - m2*pow(Q2h,2)*pow(Q2k,2) -
		          M2*pow(Q2h,2)*pow(Q2k,2) + 4*S*pow(Q2h,2)*pow(Q2k,2) +
		          2*Sk*pow(Q2h,2)*pow(Q2k,2) + M2*Q2e*pow(Q2k,3) -
		          M2*Q2h*pow(Q2k,3) - 32*m4*Q2e*pow(S,2) -
		          32*l1k2*m2*Q2h*pow(S,2) + 32*l1k2*m2*Q2k*pow(S,2) +
		          32*m4*Q2k*pow(S,2) - 32*l1k2*Q2e*Q2k*pow(S,2) -
		          48*m2*Q2e*Q2k*pow(S,2) + 16*l1k2*Q2h*Q2k*pow(S,2) +
		          16*m2*Q2h*Q2k*pow(S,2) - 32*Q2e*Q2h*Q2k*pow(S,2) +
		          16*m2*pow(Q2e,2)*pow(S,2) + 12*Q2h*pow(Q2e,2)*pow(S,2) +
		          4*Q2k*pow(Q2e,2)*pow(S,2) - 12*Q2e*pow(Q2h,2)*pow(S,2) +
		          28*Q2k*pow(Q2h,2)*pow(S,2) + 16*l1k2*pow(Q2k,2)*pow(S,2) +
		          16*m2*pow(Q2k,2)*pow(S,2) + 4*Q2e*pow(Q2k,2)*pow(S,2) -
		          4*Q2h*pow(Q2k,2)*pow(S,2) + 8*l1k2*m2*Q2e*pow(Sk,2) -
		          24*l1k2*m2*Q2h*pow(Sk,2) - 8*m2*Q2e*Q2h*pow(Sk,2) +
		          16*l1k2*m2*Q2k*pow(Sk,2) - 4*Q2e*Q2h*Q2k*pow(Sk,2) +
		          12*m2*pow(Q2e,2)*pow(Sk,2) + 4*Q2h*pow(Q2e,2)*pow(Sk,2) -
		          4*m2*pow(Q2h,2)*pow(Sk,2) - 4*Q2e*pow(Q2h,2)*pow(Sk,2) +
		          4*Q2k*pow(Q2h,2)*pow(Sk,2) - 8*m2*Q2e*Q2h*pow(Sq2,2) +
		          8*m2*Q2h*Q2k*pow(Sq2,2) - 4*Q2e*Q2h*Q2k*pow(Sq2,2) +
		          4*Q2k*pow(Q2h,2)*pow(Sq2,2)))/2. +
		     (pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		        (-32*l1k2*m2*M2*Q2e*Q2h + 32*l1k2*M2*Q2e*Q2h*Q2k -
		          16*m2*M2*Q2e*Q2h*Q2k + 16*l1k2*Q2e*Q2h*Q2k*Sk +
		          32*l1k2*m2*Q2e*Q2h*Sq2 - 24*l1k2*Q2e*Q2h*Q2k*Sq2 +
		          16*m2*Q2e*Q2h*Q2k*Sq2 + 128*l1k2*m2*Q2e*S*Sq2 -
		          128*l1k2*m2*Q2h*S*Sq2 + 256*l1k2*Q2e*Q2h*S*Sq2 -
		          64*m2*Q2e*Q2h*S*Sq2 - 64*l1k2*Q2e*Q2k*S*Sq2 +
		          64*m2*Q2e*Q2k*S*Sq2 + 192*l1k2*Q2h*Q2k*S*Sq2 -
		          64*m2*Q2h*Q2k*S*Sq2 + 96*Q2e*Q2h*Q2k*S*Sq2 -
		          32*l1k2*Q2e*Q2h*Sk*Sq2 - 16*Q2e*Q2h*Q2k*Sk*Sq2 +
		          32*M2*Q2e*Q2h*pow(l1k2,2) + 64*M2*Q2h*Q2k*pow(l1k2,2) -
		          32*Q2e*Q2h*Sq2*pow(l1k2,2) - 64*Q2h*Q2k*Sq2*pow(l1k2,2) -
		          128*Q2e*S*Sq2*pow(l1k2,2) + 384*Q2h*S*Sq2*pow(l1k2,2) -
		          256*Q2k*S*Sq2*pow(l1k2,2) + 64*M2*Q2h*pow(l1k2,3) -
		          64*Q2h*Sq2*pow(l1k2,3) - 256*S*Sq2*pow(l1k2,3) +
		          8*l1k2*M2*Q2h*pow(Q2e,2) + 4*M2*Q2h*Q2k*pow(Q2e,2) -
		          16*Q2h*Q2k*S*pow(Q2e,2) - 64*Q2h*S*Sk*pow(Q2e,2) -
		          8*l1k2*Q2h*Sq2*pow(Q2e,2) - 12*Q2h*Q2k*Sq2*pow(Q2e,2) -
		          32*l1k2*S*Sq2*pow(Q2e,2) + 48*Q2h*S*Sq2*pow(Q2e,2) -
		          16*Q2k*S*Sq2*pow(Q2e,2) - 32*Q2h*Sk*Sq2*pow(Q2e,2) +
		          32*l1k2*m2*M2*pow(Q2h,2) - 16*l1k2*m2*Q2e*pow(Q2h,2) -
		          72*l1k2*M2*Q2e*pow(Q2h,2) + 16*m2*M2*Q2e*pow(Q2h,2) -
		          64*l1k2*M2*Q2k*pow(Q2h,2) + 16*m2*M2*Q2k*pow(Q2h,2) +
		          12*l1k2*Q2e*Q2k*pow(Q2h,2) - 8*m2*Q2e*Q2k*pow(Q2h,2) -
		          36*M2*Q2e*Q2k*pow(Q2h,2) + 32*Q2e*Q2k*S*pow(Q2h,2) -
		          16*l1k2*Q2k*Sk*pow(Q2h,2) - 8*Q2e*Q2k*Sk*pow(Q2h,2) +
		          128*Q2e*S*Sk*pow(Q2h,2) - 32*l1k2*m2*Sq2*pow(Q2h,2) +
		          8*l1k2*Q2e*Sq2*pow(Q2h,2) - 16*m2*Q2e*Sq2*pow(Q2h,2) +
		          56*l1k2*Q2k*Sq2*pow(Q2h,2) - 16*m2*Q2k*Sq2*pow(Q2h,2) +
		          16*Q2e*Q2k*Sq2*pow(Q2h,2) - 288*l1k2*S*Sq2*pow(Q2h,2) +
		          64*m2*S*Sq2*pow(Q2h,2) - 160*Q2e*S*Sq2*pow(Q2h,2) -
		          80*Q2k*S*Sq2*pow(Q2h,2) + 32*l1k2*Sk*Sq2*pow(Q2h,2) +
		          80*Q2e*Sk*Sq2*pow(Q2h,2) + 16*Q2k*Sk*Sq2*pow(Q2h,2) -
		          96*M2*pow(l1k2,2)*pow(Q2h,2) + 16*Q2e*pow(l1k2,2)*pow(Q2h,2) +
		          32*Q2k*pow(l1k2,2)*pow(Q2h,2) + 96*Sq2*pow(l1k2,2)*pow(Q2h,2) +
		          32*pow(l1k2,3)*pow(Q2h,2) + 4*l1k2*pow(Q2e,2)*pow(Q2h,2) -
		          4*M2*pow(Q2e,2)*pow(Q2h,2) + 2*Q2k*pow(Q2e,2)*pow(Q2h,2) +
		          4*Sq2*pow(Q2e,2)*pow(Q2h,2) + 16*l1k2*m2*pow(Q2h,3) +
		          80*l1k2*M2*pow(Q2h,3) - 16*m2*M2*pow(Q2h,3) -
		          20*l1k2*Q2e*pow(Q2h,3) + 8*m2*Q2e*pow(Q2h,3) +
		          28*M2*Q2e*pow(Q2h,3) - 28*l1k2*Q2k*pow(Q2h,3) +
		          8*m2*Q2k*pow(Q2h,3) + 32*M2*Q2k*pow(Q2h,3) -
		          8*Q2e*Q2k*pow(Q2h,3) - 16*Q2k*S*pow(Q2h,3) +
		          8*Q2k*Sk*pow(Q2h,3) - 64*S*Sk*pow(Q2h,3) -
		          16*l1k2*Sq2*pow(Q2h,3) + 16*m2*Sq2*pow(Q2h,3) +
		          4*Q2e*Sq2*pow(Q2h,3) - 4*Q2k*Sq2*pow(Q2h,3) +
		          112*S*Sq2*pow(Q2h,3) - 48*Sk*Sq2*pow(Q2h,3) -
		          48*pow(l1k2,2)*pow(Q2h,3) - 2*pow(Q2e,2)*pow(Q2h,3) +
		          24*l1k2*pow(Q2h,4) - 8*m2*pow(Q2h,4) - 24*M2*pow(Q2h,4) +
		          6*Q2e*pow(Q2h,4) + 6*Q2k*pow(Q2h,4) - 8*Sq2*pow(Q2h,4) -
		          4*pow(Q2h,5) - 8*l1k2*M2*Q2e*pow(Q2k,2) +
		          24*l1k2*M2*Q2h*pow(Q2k,2) + 12*M2*Q2e*Q2h*pow(Q2k,2) +
		          8*Q2e*Q2h*Sk*pow(Q2k,2) - 16*l1k2*Q2h*Sq2*pow(Q2k,2) -
		          4*Q2e*Q2h*Sq2*pow(Q2k,2) - 64*l1k2*S*Sq2*pow(Q2k,2) +
		          8*l1k2*pow(Q2h,2)*pow(Q2k,2) - 12*M2*pow(Q2h,2)*pow(Q2k,2) +
		          2*Q2e*pow(Q2h,2)*pow(Q2k,2) - 8*Sk*pow(Q2h,2)*pow(Q2k,2) +
		          4*Sq2*pow(Q2h,2)*pow(Q2k,2) - 2*pow(Q2h,3)*pow(Q2k,2) -
		          4*M2*Q2e*pow(Q2k,3) + 4*M2*Q2h*pow(Q2k,3) +
		          128*l1k2*m2*Q2e*pow(S,2) - 128*l1k2*m2*Q2h*pow(S,2) +
		          256*l1k2*Q2e*Q2h*pow(S,2) - 64*m2*Q2e*Q2h*pow(S,2) -
		          64*l1k2*Q2e*Q2k*pow(S,2) + 64*m2*Q2e*Q2k*pow(S,2) +
		          192*l1k2*Q2h*Q2k*pow(S,2) - 64*m2*Q2h*Q2k*pow(S,2) +
		          96*Q2e*Q2h*Q2k*pow(S,2) - 128*Q2e*pow(l1k2,2)*pow(S,2) +
		          384*Q2h*pow(l1k2,2)*pow(S,2) - 256*Q2k*pow(l1k2,2)*pow(S,2) -
		          256*pow(l1k2,3)*pow(S,2) - 32*l1k2*pow(Q2e,2)*pow(S,2) +
		          16*Q2h*pow(Q2e,2)*pow(S,2) - 16*Q2k*pow(Q2e,2)*pow(S,2) -
		          288*l1k2*pow(Q2h,2)*pow(S,2) + 64*m2*pow(Q2h,2)*pow(S,2) -
		          96*Q2e*pow(Q2h,2)*pow(S,2) - 80*Q2k*pow(Q2h,2)*pow(S,2) +
		          80*pow(Q2h,3)*pow(S,2) - 64*l1k2*pow(Q2k,2)*pow(S,2) +
		          32*l1k2*Q2e*Q2h*pow(Sk,2) + 16*Q2e*Q2h*Q2k*pow(Sk,2) -
		          32*l1k2*pow(Q2h,2)*pow(Sk,2) - 16*Q2e*pow(Q2h,2)*pow(Sk,2) -
		          16*Q2k*pow(Q2h,2)*pow(Sk,2) + 16*pow(Q2h,3)*pow(Sk,2) +
		          64*l1k2*Q2e*Q2h*pow(Sq2,2) + 32*Q2e*Q2h*Q2k*pow(Sq2,2) +
		          16*Q2h*pow(Q2e,2)*pow(Sq2,2) - 64*l1k2*pow(Q2h,2)*pow(Sq2,2) -
		          64*Q2e*pow(Q2h,2)*pow(Sq2,2) - 32*Q2k*pow(Q2h,2)*pow(Sq2,2) +
		          48*pow(Q2h,3)*pow(Sq2,2)))/2. +
		     (pow(l1k1,-1)*(64*k1k2*l1k2*m2*M2 + 32*k1k2*l1k2*M2*Q2e +
		          32*l1k2*m2*M2*Q2h + 8*l1k2*M2*Q2e*Q2h - 4*m2*M2*Q2e*Q2h -
		          32*l1k2*m2*M2*Q2k + 12*k1k2*M2*Q2e*Q2k -
		          12*l1k2*M2*Q2e*Q2k + 4*m2*M2*Q2e*Q2k - 8*k1k2*m2*Q2h*Q2k +
		          32*m4*Q2h*Q2k - 4*k1k2*M2*Q2h*Q2k + 4*l1k2*M2*Q2h*Q2k -
		          12*m2*M2*Q2h*Q2k + 4*k1k2*Q2e*Q2h*Q2k + 14*m2*Q2e*Q2h*Q2k -
		          14*M2*Q2e*Q2h*Q2k + 32*k1k2*m2*Q2h*S - 64*l1k2*m2*Q2h*S +
		          16*k1k2*Q2e*Q2h*S - 24*l1k2*Q2e*Q2h*S + 24*m2*Q2e*Q2h*S -
		          8*l1k2*Q2h*Q2k*S - 48*m2*Q2h*Q2k*S - 16*k1k2*m2*Q2h*Sk +
		          64*m4*Q2h*Sk + 8*k1k2*Q2e*Q2h*Sk + 28*m2*Q2e*Q2h*Sk -
		          64*m4*Q2k*Sk - 40*m2*Q2e*Q2k*Sk - 8*Q2e*Q2h*Q2k*Sk -
		          32*m2*Q2h*S*Sk + 8*Q2e*Q2h*S*Sk + 32*m2*Q2k*S*Sk -
		          24*Q2e*Q2k*S*Sk + 8*Q2h*Q2k*S*Sk + 32*k1k2*m2*Q2h*Sq2 -
		          32*l1k2*m2*Q2h*Sq2 - 8*l1k2*Q2e*Q2h*Sq2 + 16*m2*Q2e*Q2h*Sq2 +
		          16*m2*Q2e*Q2k*Sq2 - 8*l1k2*Q2h*Q2k*Sq2 - 40*m2*Q2h*Q2k*Sq2 +
		          8*Q2e*Q2h*Q2k*Sq2 + 64*m2*Q2h*S*Sq2 + 16*Q2e*Q2h*S*Sq2 -
		          16*Q2e*Q2k*S*Sq2 + 16*Q2h*Q2k*S*Sq2 + 32*m2*Q2e*Sk*Sq2 -
		          32*m2*Q2h*Sk*Sq2 + 8*Q2e*Q2h*Sk*Sq2 + 32*m2*Q2k*Sk*Sq2 -
		          24*Q2e*Q2k*Sk*Sq2 + 8*Q2h*Q2k*Sk*Sq2 - 16*M2*Q2e*pow(k1k2,2) -
		          64*m2*M2*pow(l1k2,2) - 32*M2*Q2e*pow(l1k2,2) -
		          8*k1k2*M2*pow(Q2e,2) + 6*M2*Q2h*pow(Q2e,2) +
		          12*M2*Q2k*pow(Q2e,2) + 4*Q2k*S*pow(Q2e,2) +
		          16*Q2k*Sk*pow(Q2e,2) + 8*S*Sk*pow(Q2e,2) - 4*Q2h*Sq2*pow(Q2e,2) -
		          4*Q2k*Sq2*pow(Q2e,2) - 8*Sk*Sq2*pow(Q2e,2) - 8*M2*pow(Q2e,3) -
		          16*m4*pow(Q2h,2) - 8*m2*Q2e*pow(Q2h,2) -
		          4*M2*Q2e*pow(Q2h,2) - 2*m2*Q2k*pow(Q2h,2) +
		          2*M2*Q2k*pow(Q2h,2) - 2*Q2e*Q2k*pow(Q2h,2) +
		          8*m2*S*pow(Q2h,2) + 4*Q2k*S*pow(Q2h,2) - 4*m2*Sk*pow(Q2h,2) -
		          4*Q2e*Sk*pow(Q2h,2) + 4*Q2k*Sk*pow(Q2h,2) +
		          8*m2*Sq2*pow(Q2h,2) - 16*m4*pow(Q2k,2) -
		          4*m2*M2*pow(Q2k,2) - 10*m2*Q2e*pow(Q2k,2) -
		          6*M2*Q2e*pow(Q2k,2) - 2*m2*Q2h*pow(Q2k,2) +
		          4*M2*Q2h*pow(Q2k,2) - 4*Q2e*Q2h*pow(Q2k,2) +
		          16*m2*S*pow(Q2k,2) - 12*Q2e*S*pow(Q2k,2) + 4*Q2h*S*pow(Q2k,2) +
		          16*m2*Sq2*pow(Q2k,2) - 12*Q2e*Sq2*pow(Q2k,2) +
		          4*Q2h*Sq2*pow(Q2k,2) + 16*S*Sq2*pow(Q2k,2) +
		          4*pow(Q2e,2)*pow(Q2k,2) + 2*pow(Q2h,2)*pow(Q2k,2) +
		          32*m2*Q2e*pow(S,2) + 32*m2*Q2h*pow(S,2) + 24*Q2e*Q2h*pow(S,2) -
		          24*Q2e*Q2k*pow(S,2) + 24*Q2h*Q2k*pow(S,2) +
		          8*pow(Q2k,2)*pow(S,2) - 64*m4*pow(Sk,2) - 40*m2*Q2e*pow(Sk,2) +
		          8*m2*Q2h*pow(Sk,2) + 16*pow(Q2e,2)*pow(Sk,2) -
		          32*m2*Q2e*pow(Sq2,2) + 32*m2*Q2h*pow(Sq2,2) +
		          8*Q2e*Q2h*pow(Sq2,2) + 8*Q2e*Q2k*pow(Sq2,2) -
		          8*Q2h*Q2k*pow(Sq2,2) + 8*pow(Q2k,2)*pow(Sq2,2)))/2. +
		     (pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		        (64*k1k2*l1k2*m2*M2*Q2h - 64*l1k2*m4*M2*Q2h -
		          32*k1k2*l1k2*M2*Q2h*Q2k + 128*k1k2*l1k2*m2*Q2h*S -
		          128*l1k2*m4*Q2h*S - 64*k1k2*l1k2*Q2h*Q2k*S +
		          128*k1k2*l1k2*m2*Q2h*Sq2 - 128*l1k2*m4*Q2h*Sq2 -
		          64*k1k2*l1k2*Q2h*Q2k*Sq2 - 512*k1k2*l1k2*m2*S*Sq2 +
		          512*l1k2*m4*S*Sq2 - 256*k1k2*l1k2*Q2h*S*Sq2 +
		          256*k1k2*l1k2*Q2k*S*Sq2 + 256*l1k2*Q2h*Q2k*S*Sq2 +
		          32*l1k2*M2*Q2h*pow(k1k2,2) + 64*l1k2*Q2h*S*pow(k1k2,2) +
		          64*l1k2*Q2h*Sq2*pow(k1k2,2) - 256*l1k2*S*Sq2*pow(k1k2,2) -
		          64*k1k2*M2*Q2h*pow(l1k2,2) + 64*M2*Q2h*Q2k*pow(l1k2,2) -
		          128*k1k2*Q2h*S*pow(l1k2,2) + 128*Q2h*Q2k*S*pow(l1k2,2) -
		          128*k1k2*Q2h*Sq2*pow(l1k2,2) + 128*Q2h*Q2k*Sq2*pow(l1k2,2) +
		          512*k1k2*S*Sq2*pow(l1k2,2) + 512*Q2h*S*Sq2*pow(l1k2,2) -
		          512*Q2k*S*Sq2*pow(l1k2,2) + 64*M2*Q2h*pow(l1k2,3) +
		          128*Q2h*S*pow(l1k2,3) + 128*Q2h*Sq2*pow(l1k2,3) -
		          512*S*Sq2*pow(l1k2,3) + 32*k1k2*l1k2*M2*pow(Q2h,2) -
		          32*l1k2*M2*Q2k*pow(Q2h,2) + 64*k1k2*l1k2*S*pow(Q2h,2) -
		          64*l1k2*Q2k*S*pow(Q2h,2) + 64*k1k2*l1k2*Sq2*pow(Q2h,2) -
		          64*l1k2*Q2k*Sq2*pow(Q2h,2) - 128*l1k2*S*Sq2*pow(Q2h,2) -
		          64*M2*pow(l1k2,2)*pow(Q2h,2) - 128*S*pow(l1k2,2)*pow(Q2h,2) -
		          128*Sq2*pow(l1k2,2)*pow(Q2h,2) + 16*l1k2*M2*pow(Q2h,3) +
		          32*l1k2*S*pow(Q2h,3) + 32*l1k2*Sq2*pow(Q2h,3) +
		          16*l1k2*M2*Q2h*pow(Q2k,2) + 32*l1k2*Q2h*S*pow(Q2k,2) +
		          32*l1k2*Q2h*Sq2*pow(Q2k,2) - 128*l1k2*S*Sq2*pow(Q2k,2) -
		          256*k1k2*l1k2*m2*pow(S,2) + 256*l1k2*m4*pow(S,2) -
		          128*k1k2*l1k2*Q2h*pow(S,2) + 128*k1k2*l1k2*Q2k*pow(S,2) +
		          128*l1k2*Q2h*Q2k*pow(S,2) - 128*l1k2*pow(k1k2,2)*pow(S,2) +
		          256*k1k2*pow(l1k2,2)*pow(S,2) + 256*Q2h*pow(l1k2,2)*pow(S,2) -
		          256*Q2k*pow(l1k2,2)*pow(S,2) - 256*pow(l1k2,3)*pow(S,2) -
		          64*l1k2*pow(Q2h,2)*pow(S,2) - 64*l1k2*pow(Q2k,2)*pow(S,2) -
		          256*k1k2*l1k2*m2*pow(Sq2,2) + 256*l1k2*m4*pow(Sq2,2) -
		          128*k1k2*l1k2*Q2h*pow(Sq2,2) + 128*k1k2*l1k2*Q2k*pow(Sq2,2) +
		          128*l1k2*Q2h*Q2k*pow(Sq2,2) - 128*l1k2*pow(k1k2,2)*pow(Sq2,2) +
		          256*k1k2*pow(l1k2,2)*pow(Sq2,2) +
		          256*Q2h*pow(l1k2,2)*pow(Sq2,2) - 256*Q2k*pow(l1k2,2)*pow(Sq2,2) -
		          256*pow(l1k2,3)*pow(Sq2,2) - 64*l1k2*pow(Q2h,2)*pow(Sq2,2) -
		          64*l1k2*pow(Q2k,2)*pow(Sq2,2)))/2. +
		     (pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		        (-16*k1k2*m2*M2*Q2e*Q2h + 16*k1k2*m2*M2*Q2e*Q2k -
		          16*k1k2*m2*M2*Q2h*Q2k - 4*k1k2*m2*Q2e*Q2h*Q2k -
		          8*m4*Q2e*Q2h*Q2k + 16*k1k2*M2*Q2e*Q2h*Q2k +
		          24*m2*M2*Q2e*Q2h*Q2k - 16*k1k2*m2*Q2e*Q2h*S -
		          32*m4*Q2e*Q2h*S - 32*k1k2*m2*Q2e*Q2k*S + 32*m4*Q2e*Q2k*S +
		          32*k1k2*m2*Q2h*Q2k*S - 32*m4*Q2h*Q2k*S -
		          16*k1k2*Q2e*Q2h*Q2k*S - 8*k1k2*m2*Q2e*Q2h*Sk -
		          16*m4*Q2e*Q2h*Sk + 48*k1k2*m2*Q2e*Q2k*Sk -
		          48*k1k2*m2*Q2h*Q2k*Sk - 8*m2*Q2e*Q2h*Q2k*Sk -
		          64*k1k2*m2*Q2e*S*Sk + 64*m4*Q2e*S*Sk + 64*k1k2*m2*Q2h*S*Sk -
		          64*m4*Q2h*S*Sk - 32*k1k2*Q2e*Q2h*S*Sk + 32*m2*Q2e*Q2k*S*Sk -
		          32*m2*Q2h*Q2k*S*Sk - 16*k1k2*m2*Q2e*Q2h*Sq2 -
		          32*m4*Q2e*Q2h*Sq2 - 32*k1k2*m2*Q2e*Q2k*Sq2 +
		          32*m4*Q2e*Q2k*Sq2 + 32*k1k2*m2*Q2h*Q2k*Sq2 -
		          32*m4*Q2h*Q2k*Sq2 - 16*k1k2*Q2e*Q2h*Q2k*Sq2 +
		          20*m2*Q2e*Q2h*Q2k*Sq2 + 128*k1k2*m2*Q2e*S*Sq2 -
		          128*k1k2*m2*Q2h*S*Sq2 - 64*m2*Q2e*Q2h*S*Sq2 -
		          64*k1k2*Q2e*Q2k*S*Sq2 - 96*m2*Q2e*Q2k*S*Sq2 +
		          64*k1k2*Q2h*Q2k*S*Sq2 + 96*m2*Q2h*Q2k*S*Sq2 -
		          48*Q2e*Q2h*Q2k*S*Sq2 - 64*k1k2*m2*Q2e*Sk*Sq2 +
		          64*m4*Q2e*Sk*Sq2 + 64*k1k2*m2*Q2h*Sk*Sq2 - 64*m4*Q2h*Sk*Sq2 -
		          32*k1k2*Q2e*Q2h*Sk*Sq2 - 16*m2*Q2e*Q2k*Sk*Sq2 +
		          16*m2*Q2h*Q2k*Sk*Sq2 - 8*M2*Q2e*Q2h*pow(k1k2,2) -
		          24*k1k2*M2*Q2h*pow(Q2e,2) - 8*m2*M2*Q2h*pow(Q2e,2) -
		          12*k1k2*M2*Q2k*pow(Q2e,2) - 4*k1k2*Q2h*Q2k*pow(Q2e,2) -
		          12*m2*Q2h*Q2k*pow(Q2e,2) + 26*M2*Q2h*Q2k*pow(Q2e,2) -
		          8*k1k2*Q2h*S*pow(Q2e,2) - 16*m2*Q2h*S*pow(Q2e,2) +
		          16*k1k2*Q2k*S*pow(Q2e,2) + 12*Q2h*Q2k*S*pow(Q2e,2) -
		          8*k1k2*Q2h*Sk*pow(Q2e,2) - 24*m2*Q2h*Sk*pow(Q2e,2) +
		          24*m2*Q2k*Sk*pow(Q2e,2) + 8*Q2h*Q2k*Sk*pow(Q2e,2) +
		          32*k1k2*S*Sk*pow(Q2e,2) - 8*Q2k*S*Sk*pow(Q2e,2) -
		          16*m2*Q2h*Sq2*pow(Q2e,2) + 16*k1k2*Q2k*Sq2*pow(Q2e,2) +
		          4*Q2h*Q2k*Sq2*pow(Q2e,2) + 64*m2*S*Sq2*pow(Q2e,2) +
		          16*Q2h*S*Sq2*pow(Q2e,2) + 32*k1k2*Sk*Sq2*pow(Q2e,2) -
		          8*Q2k*Sk*Sq2*pow(Q2e,2) + 8*M2*pow(k1k2,2)*pow(Q2e,2) +
		          16*k1k2*M2*pow(Q2e,3) - 12*M2*Q2h*pow(Q2e,3) -
		          12*M2*Q2k*pow(Q2e,3) - 4*Q2h*Q2k*pow(Q2e,3) -
		          12*Q2h*S*pow(Q2e,3) + 4*Q2k*S*pow(Q2e,3) - 8*Q2h*Sk*pow(Q2e,3) +
		          8*S*Sk*pow(Q2e,3) + 4*Q2k*Sq2*pow(Q2e,3) + 8*Sk*Sq2*pow(Q2e,3) +
		          6*M2*pow(Q2e,4) + 16*k1k2*m2*M2*pow(Q2h,2) +
		          8*m4*Q2e*pow(Q2h,2) + 8*k1k2*M2*Q2e*pow(Q2h,2) +
		          4*m2*M2*Q2e*pow(Q2h,2) + 4*k1k2*m2*Q2k*pow(Q2h,2) +
		          8*m4*Q2k*pow(Q2h,2) - 4*k1k2*M2*Q2k*pow(Q2h,2) -
		          24*m2*M2*Q2k*pow(Q2h,2) + 4*k1k2*Q2e*Q2k*pow(Q2h,2) +
		          10*m2*Q2e*Q2k*pow(Q2h,2) - 16*M2*Q2e*Q2k*pow(Q2h,2) +
		          16*k1k2*m2*S*pow(Q2h,2) + 32*m4*S*pow(Q2h,2) +
		          8*k1k2*Q2e*S*pow(Q2h,2) + 16*m2*Q2e*S*pow(Q2h,2) -
		          20*Q2e*Q2k*S*pow(Q2h,2) + 8*k1k2*m2*Sk*pow(Q2h,2) +
		          16*m4*Sk*pow(Q2h,2) + 8*k1k2*Q2e*Sk*pow(Q2h,2) +
		          24*m2*Q2e*Sk*pow(Q2h,2) - 16*m2*Q2k*Sk*pow(Q2h,2) -
		          12*Q2e*Q2k*Sk*pow(Q2h,2) - 8*Q2e*S*Sk*pow(Q2h,2) +
		          8*Q2k*S*Sk*pow(Q2h,2) + 16*k1k2*m2*Sq2*pow(Q2h,2) +
		          32*m4*Sq2*pow(Q2h,2) + 12*m2*Q2e*Sq2*pow(Q2h,2) -
		          20*m2*Q2k*Sq2*pow(Q2h,2) - 8*Q2e*Q2k*Sq2*pow(Q2h,2) -
		          16*Q2e*S*Sq2*pow(Q2h,2) + 48*Q2k*S*Sq2*pow(Q2h,2) -
		          8*Q2e*Sk*Sq2*pow(Q2h,2) + 8*Q2k*Sk*Sq2*pow(Q2h,2) +
		          6*m2*pow(Q2e,2)*pow(Q2h,2) + 8*M2*pow(Q2e,2)*pow(Q2h,2) +
		          6*Q2k*pow(Q2e,2)*pow(Q2h,2) + 16*S*pow(Q2e,2)*pow(Q2h,2) +
		          12*Sk*pow(Q2e,2)*pow(Q2h,2) - 8*m4*pow(Q2h,3) +
		          4*m2*M2*pow(Q2h,3) - 6*m2*Q2e*pow(Q2h,3) -
		          2*M2*Q2e*pow(Q2h,3) + 2*m2*Q2k*pow(Q2h,3) +
		          2*M2*Q2k*pow(Q2h,3) - 2*Q2e*Q2k*pow(Q2h,3) -
		          4*Q2e*S*pow(Q2h,3) + 4*Q2k*S*pow(Q2h,3) - 4*Q2e*Sk*pow(Q2h,3) +
		          4*Q2k*Sk*pow(Q2h,3) + 4*m2*Sq2*pow(Q2h,3) +
		          12*k1k2*m2*Q2e*pow(Q2k,2) - 4*m2*M2*Q2e*pow(Q2k,2) -
		          12*k1k2*m2*Q2h*pow(Q2k,2) + 4*m2*M2*Q2h*pow(Q2k,2) -
		          10*M2*Q2e*Q2h*pow(Q2k,2) + 16*m2*Q2e*S*pow(Q2k,2) -
		          16*m2*Q2h*S*pow(Q2k,2) - 8*m2*Q2e*Sq2*pow(Q2k,2) +
		          8*m2*Q2h*Sq2*pow(Q2k,2) + 16*Q2e*S*Sq2*pow(Q2k,2) -
		          16*Q2h*S*Sq2*pow(Q2k,2) + 6*m2*pow(Q2e,2)*pow(Q2k,2) +
		          6*M2*pow(Q2e,2)*pow(Q2k,2) + 4*Q2h*pow(Q2e,2)*pow(Q2k,2) -
		          4*S*pow(Q2e,2)*pow(Q2k,2) - 4*Sq2*pow(Q2e,2)*pow(Q2k,2) -
		          6*m2*pow(Q2h,2)*pow(Q2k,2) + 4*M2*pow(Q2h,2)*pow(Q2k,2) -
		          6*Q2e*pow(Q2h,2)*pow(Q2k,2) + 4*S*pow(Q2h,2)*pow(Q2k,2) +
		          4*Sq2*pow(Q2h,2)*pow(Q2k,2) + 2*pow(Q2h,3)*pow(Q2k,2) +
		          64*k1k2*m2*Q2e*pow(S,2) - 64*k1k2*m2*Q2h*pow(S,2) -
		          32*m2*Q2e*Q2h*pow(S,2) - 32*k1k2*Q2e*Q2k*pow(S,2) -
		          64*m2*Q2e*Q2k*pow(S,2) + 32*k1k2*Q2h*Q2k*pow(S,2) +
		          64*m2*Q2h*Q2k*pow(S,2) - 48*Q2e*Q2h*Q2k*pow(S,2) +
		          32*m2*pow(Q2e,2)*pow(S,2) + 8*Q2h*pow(Q2e,2)*pow(S,2) +
		          8*Q2k*pow(Q2e,2)*pow(S,2) - 8*Q2e*pow(Q2h,2)*pow(S,2) +
		          40*Q2k*pow(Q2h,2)*pow(S,2) + 8*Q2e*pow(Q2k,2)*pow(S,2) -
		          8*Q2h*pow(Q2k,2)*pow(S,2) + 48*k1k2*m2*Q2e*pow(Sk,2) -
		          48*k1k2*m2*Q2h*pow(Sk,2) - 16*m2*Q2e*Q2h*pow(Sk,2) +
		          24*m2*pow(Q2e,2)*pow(Sk,2) - 8*m2*pow(Q2h,2)*pow(Sk,2) +
		          64*k1k2*m2*Q2e*pow(Sq2,2) - 64*k1k2*m2*Q2h*pow(Sq2,2) -
		          32*m2*Q2e*Q2h*pow(Sq2,2) - 32*k1k2*Q2e*Q2k*pow(Sq2,2) -
		          32*m2*Q2e*Q2k*pow(Sq2,2) + 32*k1k2*Q2h*Q2k*pow(Sq2,2) +
		          32*m2*Q2h*Q2k*pow(Sq2,2) + 32*m2*pow(Q2e,2)*pow(Sq2,2) -
		          8*Q2k*pow(Q2e,2)*pow(Sq2,2) + 8*Q2k*pow(Q2h,2)*pow(Sq2,2) +
		          8*Q2e*pow(Q2k,2)*pow(Sq2,2) - 8*Q2h*pow(Q2k,2)*pow(Sq2,2)))/2. +
		     (pow(l1k1,-1)*pow(k1k2 - l1k1 - l1k2,-1)*
		        pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		        (-64*m4*Q2e*Q2h*S*Sq2 - 64*m2*Q2e*Q2h*Q2k*S*Sq2 -
		          8*m4*M2*Q2h*pow(Q2e,2) - 8*m2*M2*Q2h*Q2k*pow(Q2e,2) -
		          4*m2*Q2h*Q2k*Sk*pow(Q2e,2) + 16*Q2h*Q2k*S*Sk*pow(Q2e,2) +
		          8*m4*Q2h*Sq2*pow(Q2e,2) + 6*m2*Q2h*Q2k*Sq2*pow(Q2e,2) +
		          32*m4*S*Sq2*pow(Q2e,2) + 16*m2*Q2h*S*Sq2*pow(Q2e,2) +
		          16*m2*Q2k*S*Sq2*pow(Q2e,2) + 28*Q2h*Q2k*S*Sq2*pow(Q2e,2) +
		          8*m2*Q2h*Sk*Sq2*pow(Q2e,2) + 4*Q2h*Q2k*Sk*Sq2*pow(Q2e,2) -
		          M2*Q2h*Q2k*pow(Q2e,3) - 4*Q2h*Q2k*S*pow(Q2e,3) -
		          2*Q2h*Q2k*Sk*pow(Q2e,3) - 16*Q2h*S*Sk*pow(Q2e,3) -
		          2*Q2h*Q2k*Sq2*pow(Q2e,3) - 4*Q2h*S*Sq2*pow(Q2e,3) -
		          4*Q2k*S*Sq2*pow(Q2e,3) - 4*Q2h*Sk*Sq2*pow(Q2e,3) +
		          16*m4*M2*Q2e*pow(Q2h,2) + 24*m2*M2*Q2e*Q2k*pow(Q2h,2) +
		          8*m2*Q2e*Q2k*Sk*pow(Q2h,2) - 32*Q2e*Q2k*S*Sk*pow(Q2h,2) -
		          16*m4*Q2e*Sq2*pow(Q2h,2) + 4*m2*Q2e*Q2k*Sq2*pow(Q2h,2) +
		          32*m4*S*Sq2*pow(Q2h,2) - 16*m2*Q2e*S*Sq2*pow(Q2h,2) +
		          48*m2*Q2k*S*Sq2*pow(Q2h,2) - 44*Q2e*Q2k*S*Sq2*pow(Q2h,2) -
		          16*m2*Q2e*Sk*Sq2*pow(Q2h,2) - 8*Q2e*Q2k*Sk*Sq2*pow(Q2h,2) -
		          4*m4*pow(Q2e,2)*pow(Q2h,2) - 2*m2*M2*pow(Q2e,2)*pow(Q2h,2) -
		          3*m2*Q2k*pow(Q2e,2)*pow(Q2h,2) +
		          2*M2*Q2k*pow(Q2e,2)*pow(Q2h,2) + 8*Q2k*S*pow(Q2e,2)*pow(Q2h,2) +
		          4*Q2k*Sk*pow(Q2e,2)*pow(Q2h,2) + 32*S*Sk*pow(Q2e,2)*pow(Q2h,2) -
		          6*m2*Sq2*pow(Q2e,2)*pow(Q2h,2) +
		          4*Q2k*Sq2*pow(Q2e,2)*pow(Q2h,2) + 8*S*Sq2*pow(Q2e,2)*pow(Q2h,2) +
		          8*Sk*Sq2*pow(Q2e,2)*pow(Q2h,2) - 8*m4*M2*pow(Q2h,3) +
		          8*m4*Q2e*pow(Q2h,3) - 16*m2*M2*Q2k*pow(Q2h,3) +
		          4*m2*Q2e*Q2k*pow(Q2h,3) - M2*Q2e*Q2k*pow(Q2h,3) -
		          4*Q2e*Q2k*S*pow(Q2h,3) - 4*m2*Q2k*Sk*pow(Q2h,3) -
		          2*Q2e*Q2k*Sk*pow(Q2h,3) - 16*Q2e*S*Sk*pow(Q2h,3) +
		          16*Q2k*S*Sk*pow(Q2h,3) + 8*m4*Sq2*pow(Q2h,3) +
		          4*m2*Q2e*Sq2*pow(Q2h,3) - 10*m2*Q2k*Sq2*pow(Q2h,3) -
		          2*Q2e*Q2k*Sq2*pow(Q2h,3) - 4*Q2e*S*Sq2*pow(Q2h,3) +
		          20*Q2k*S*Sq2*pow(Q2h,3) + 8*m2*Sk*Sq2*pow(Q2h,3) -
		          4*Q2e*Sk*Sq2*pow(Q2h,3) + 4*Q2k*Sk*Sq2*pow(Q2h,3) +
		          m2*pow(Q2e,2)*pow(Q2h,3) - 4*m4*pow(Q2h,4) +
		          2*m2*M2*pow(Q2h,4) - m2*Q2e*pow(Q2h,4) - m2*Q2k*pow(Q2h,4) +
		          2*m2*Sq2*pow(Q2h,4) - 8*m2*M2*Q2e*Q2h*pow(Q2k,2) -
		          8*m2*Q2e*Q2h*Sq2*pow(Q2k,2) + 16*m2*Q2e*S*Sq2*pow(Q2k,2) -
		          16*m2*Q2h*S*Sq2*pow(Q2k,2) + 8*Q2e*Q2h*S*Sq2*pow(Q2k,2) +
		          2*m2*M2*pow(Q2e,2)*pow(Q2k,2) -
		          3*M2*Q2h*pow(Q2e,2)*pow(Q2k,2) + 4*Q2h*S*pow(Q2e,2)*pow(Q2k,2) +
		          2*Q2h*Sk*pow(Q2e,2)*pow(Q2k,2) - 4*S*Sq2*pow(Q2e,2)*pow(Q2k,2) +
		          M2*pow(Q2e,3)*pow(Q2k,2) + 6*m2*M2*pow(Q2h,2)*pow(Q2k,2) +
		          m2*Q2e*pow(Q2h,2)*pow(Q2k,2) + 3*M2*Q2e*pow(Q2h,2)*pow(Q2k,2) -
		          8*Q2e*S*pow(Q2h,2)*pow(Q2k,2) - 4*Q2e*Sk*pow(Q2h,2)*pow(Q2k,2) +
		          8*m2*Sq2*pow(Q2h,2)*pow(Q2k,2) - 4*S*Sq2*pow(Q2h,2)*pow(Q2k,2) -
		          m2*pow(Q2h,3)*pow(Q2k,2) - M2*pow(Q2h,3)*pow(Q2k,2) +
		          4*S*pow(Q2h,3)*pow(Q2k,2) + 2*Sk*pow(Q2h,3)*pow(Q2k,2) +
		          2*M2*Q2e*Q2h*pow(Q2k,3) - M2*pow(Q2e,2)*pow(Q2k,3) -
		          M2*pow(Q2h,2)*pow(Q2k,3) - 64*m4*Q2e*Q2h*pow(S,2) -
		          64*m2*Q2e*Q2h*Q2k*pow(S,2) + 32*m4*pow(Q2e,2)*pow(S,2) +
		          16*m2*Q2h*pow(Q2e,2)*pow(S,2) + 16*m2*Q2k*pow(Q2e,2)*pow(S,2) +
		          36*Q2h*Q2k*pow(Q2e,2)*pow(S,2) - 12*Q2h*pow(Q2e,3)*pow(S,2) -
		          4*Q2k*pow(Q2e,3)*pow(S,2) + 32*m4*pow(Q2h,2)*pow(S,2) -
		          16*m2*Q2e*pow(Q2h,2)*pow(S,2) + 48*m2*Q2k*pow(Q2h,2)*pow(S,2) -
		          60*Q2e*Q2k*pow(Q2h,2)*pow(S,2) +
		          24*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) - 12*Q2e*pow(Q2h,3)*pow(S,2) +
		          28*Q2k*pow(Q2h,3)*pow(S,2) + 16*m2*Q2e*pow(Q2k,2)*pow(S,2) -
		          16*m2*Q2h*pow(Q2k,2)*pow(S,2) + 8*Q2e*Q2h*pow(Q2k,2)*pow(S,2) -
		          4*pow(Q2e,2)*pow(Q2k,2)*pow(S,2) -
		          4*pow(Q2h,2)*pow(Q2k,2)*pow(S,2) -
		          8*m2*Q2h*pow(Q2e,2)*pow(Sk,2) + 4*Q2h*Q2k*pow(Q2e,2)*pow(Sk,2) -
		          4*Q2h*pow(Q2e,3)*pow(Sk,2) + 16*m2*Q2e*pow(Q2h,2)*pow(Sk,2) -
		          8*Q2e*Q2k*pow(Q2h,2)*pow(Sk,2) +
		          8*pow(Q2e,2)*pow(Q2h,2)*pow(Sk,2) - 8*m2*pow(Q2h,3)*pow(Sk,2) -
		          4*Q2e*pow(Q2h,3)*pow(Sk,2) + 4*Q2k*pow(Q2h,3)*pow(Sk,2) -
		          24*m2*Q2e*Q2h*Q2k*pow(Sq2,2) + 8*m2*Q2h*pow(Q2e,2)*pow(Sq2,2) +
		          4*Q2h*Q2k*pow(Q2e,2)*pow(Sq2,2) -
		          4*m2*Q2e*pow(Q2h,2)*pow(Sq2,2) +
		          24*m2*Q2k*pow(Q2h,2)*pow(Sq2,2) -
		          8*Q2e*Q2k*pow(Q2h,2)*pow(Sq2,2) - 4*m2*pow(Q2h,3)*pow(Sq2,2) +
		          4*Q2k*pow(Q2h,3)*pow(Sq2,2) + 12*m2*Q2e*pow(Q2k,2)*pow(Sq2,2) -
		          12*m2*Q2h*pow(Q2k,2)*pow(Sq2,2)))/2.);
}

long double Melem::melem2_l1k2_1(const long double l1k1,const long double l1k2,const long double k1k2,const long double Q2e,
		const long double Q2h,const long double Q2k,const long double S,const long double Sk,const long double Sq2,
		const long double f1, const long double f2)const {

	long double p1p2 = M2 + Q2h/2.;
	long double G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
	long double G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
	long double G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;

	return G1*(-24*k1k2*m2 - 24*l1k1*m2 + 8*l1k2*m2 + 4*k1k2*Q2e - 20*l1k1*Q2e -
		     4*m2*Q2e + 4*k1k2*Q2h + 28*l1k1*Q2h - 40*l1k2*Q2h - 8*m2*Q2h +
		     12*m2*Q2k + 6*Q2e*Q2k - 22*Q2h*Q2k - 2*pow(Q2e,2) + 26*pow(Q2h,2) +
		     pow(l1k2,-1)*(48*k1k2*l1k1*m2 + 16*k1k2*m4 + 8*k1k2*l1k1*Q2e +
		        4*k1k2*m2*Q2e + 16*l1k1*m2*Q2e + 8*m4*Q2e - 12*k1k2*m2*Q2h +
		        32*l1k1*m2*Q2h + 8*m4*Q2h - 4*k1k2*Q2e*Q2h + 2*l1k1*Q2e*Q2h -
		        8*k1k2*m2*Q2k - 40*l1k1*m2*Q2k - 6*m2*Q2e*Q2k + 4*l1k1*Q2h*Q2k +
		        2*m2*Q2h*Q2k - 4*Q2e*pow(k1k2,2) + 32*m2*pow(l1k1,2) -
		        8*Q2e*pow(l1k1,2) - 4*Q2h*pow(l1k1,2) - 4*l1k1*pow(Q2e,2) +
		        4*m2*pow(Q2e,2) - Q2h*pow(Q2e,2) - pow(Q2e,3) -
		        6*l1k1*pow(Q2h,2) - 4*m2*pow(Q2h,2) - Q2e*pow(Q2h,2) +
		        4*Q2k*pow(Q2h,2) - 4*pow(Q2h,3) + 4*m2*pow(Q2k,2) - Q2h*pow(Q2k,2)) \
		+ pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (144*k1k2*l1k2*m2 - 64*k1k2*m4 + 32*m6 + 24*k1k2*l1k2*Q2e +
		        48*k1k2*m2*Q2e - 32*l1k2*m2*Q2e - 24*k1k2*l1k2*Q2h -
		        32*k1k2*m2*Q2h + 120*l1k2*m2*Q2h - 48*m4*Q2h - 8*k1k2*Q2e*Q2h +
		        64*l1k2*Q2e*Q2h + 20*m2*Q2e*Q2h + 32*k1k2*m2*Q2k -
		        88*l1k2*m2*Q2k + 32*m4*Q2k + 12*k1k2*Q2e*Q2k - 24*l1k2*Q2e*Q2k -
		        20*m2*Q2e*Q2k - 12*k1k2*Q2h*Q2k + 64*l1k2*Q2h*Q2k +
		        36*m2*Q2h*Q2k + 28*Q2e*Q2h*Q2k - 16*m2*pow(k1k2,2) -
		        8*Q2e*pow(k1k2,2) + 8*Q2h*pow(k1k2,2) - 112*m2*pow(l1k2,2) -
		        32*Q2e*pow(l1k2,2) + 96*Q2h*pow(l1k2,2) + 4*k1k2*pow(Q2e,2) -
		        12*l1k2*pow(Q2e,2) + 4*m2*pow(Q2e,2) + 8*Q2h*pow(Q2e,2) -
		        6*Q2k*pow(Q2e,2) - 2*pow(Q2e,3) + 12*k1k2*pow(Q2h,2) -
		        108*l1k2*pow(Q2h,2) - 24*m2*pow(Q2h,2) - 32*Q2e*pow(Q2h,2) -
		        38*Q2k*pow(Q2h,2) + 38*pow(Q2h,3) - 16*m2*pow(Q2k,2) -
		        6*Q2e*pow(Q2k,2) + 10*Q2h*pow(Q2k,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*
		      (40*l1k1*l1k2*m2 + 12*l1k1*l1k2*Q2e + 24*l1k1*m2*Q2e -
		        12*l1k2*m2*Q2e + 24*m4*Q2e - 36*l1k1*l1k2*Q2h + 48*l1k2*m2*Q2h -
		        24*m4*Q2h - 12*l1k1*Q2e*Q2h + 36*l1k2*Q2e*Q2h + 16*m2*Q2e*Q2h -
		        16*l1k1*m2*Q2k - 36*l1k2*m2*Q2k + 14*l1k1*Q2e*Q2k -
		        14*l1k2*Q2e*Q2k - 22*m2*Q2e*Q2k - 30*l1k1*Q2h*Q2k +
		        34*l1k2*Q2h*Q2k + 14*m2*Q2h*Q2k + 23*Q2e*Q2h*Q2k +
		        24*m2*pow(l1k1,2) - 12*Q2e*pow(l1k1,2) + 28*Q2h*pow(l1k1,2) -
		        32*m2*pow(l1k2,2) - 16*Q2e*pow(l1k2,2) + 32*Q2h*pow(l1k2,2) +
		        2*l1k1*pow(Q2e,2) - 8*l1k2*pow(Q2e,2) + 6*m2*pow(Q2e,2) +
		        7*Q2h*pow(Q2e,2) - 5*Q2k*pow(Q2e,2) - 2*pow(Q2e,3) +
		        30*l1k1*pow(Q2h,2) - 48*l1k2*pow(Q2h,2) - 16*m2*pow(Q2h,2) -
		        23*Q2e*pow(Q2h,2) - 28*Q2k*pow(Q2h,2) + 24*pow(Q2h,3) +
		        2*m2*pow(Q2k,2) - 6*Q2e*pow(Q2k,2) + 10*Q2h*pow(Q2k,2)) +
		     pow(k1k2 - l1k1 - l1k2,-2)*
		      (4*l1k1*l1k2*m2*Q2e + 8*l1k1*m4*Q2e + 4*l1k1*l1k2*m2*Q2h +
		        8*l1k1*m4*Q2h - 2*l1k1*l1k2*Q2e*Q2h - 4*l1k1*m2*Q2e*Q2h +
		        2*l1k2*m2*Q2e*Q2h + 4*m4*Q2e*Q2h - 8*l1k1*l1k2*m2*Q2k -
		        16*l1k1*m4*Q2k - 2*l1k2*m2*Q2e*Q2k - 4*m4*Q2e*Q2k +
		        4*l1k1*l1k2*Q2h*Q2k + 8*l1k1*m2*Q2h*Q2k - 2*l1k2*m2*Q2h*Q2k -
		        4*m4*Q2h*Q2k + l1k2*Q2e*Q2h*Q2k + 2*m2*Q2e*Q2h*Q2k +
		        8*l1k2*m2*pow(l1k1,2) + 16*m4*pow(l1k1,2) -
		        4*l1k2*Q2h*pow(l1k1,2) - 8*m2*Q2h*pow(l1k1,2) -
		        2*l1k1*l1k2*pow(Q2h,2) - 4*l1k1*m2*pow(Q2h,2) -
		        l1k2*Q2e*pow(Q2h,2) - 2*m2*Q2e*pow(Q2h,2) + l1k2*Q2k*pow(Q2h,2) +
		        2*m2*Q2k*pow(Q2h,2) + 2*l1k2*m2*pow(Q2k,2) + 4*m4*pow(Q2k,2) -
		        l1k2*Q2h*pow(Q2k,2) - 2*m2*Q2h*pow(Q2k,2)) +
		     pow(k1k2 - l1k1 - l1k2,-2)*pow(l1k2,-1)*
		      (-8*l1k1*m6*Q2e - 8*l1k1*m6*Q2h + 8*l1k1*m4*Q2e*Q2h -
		        4*m6*Q2e*Q2h + 16*l1k1*m6*Q2k - 4*l1k1*m4*Q2e*Q2k +
		        4*m6*Q2e*Q2k - 12*l1k1*m4*Q2h*Q2k + 4*m6*Q2h*Q2k +
		        2*l1k1*m2*Q2e*Q2h*Q2k - 2*m4*Q2e*Q2h*Q2k - 16*m6*pow(l1k1,2) +
		        8*m4*Q2e*pow(l1k1,2) + 16*m4*Q2h*pow(l1k1,2) -
		        2*m2*Q2e*Q2h*pow(l1k1,2) - 16*m4*Q2k*pow(l1k1,2) -
		        2*m2*Q2e*Q2k*pow(l1k1,2) + 6*m2*Q2h*Q2k*pow(l1k1,2) +
		        Q2e*Q2h*Q2k*pow(l1k1,2) + 16*m4*pow(l1k1,3) +
		        4*m2*Q2e*pow(l1k1,3) - 4*m2*Q2h*pow(l1k1,3) -
		        2*Q2e*Q2h*pow(l1k1,3) - 8*m2*Q2k*pow(l1k1,3) +
		        4*Q2h*Q2k*pow(l1k1,3) + 8*m2*pow(l1k1,4) - 4*Q2h*pow(l1k1,4) +
		        4*l1k1*m4*pow(Q2h,2) - 2*l1k1*m2*Q2e*pow(Q2h,2) +
		        2*m4*Q2e*pow(Q2h,2) + 2*l1k1*m2*Q2k*pow(Q2h,2) -
		        2*m4*Q2k*pow(Q2h,2) - 4*m2*pow(l1k1,2)*pow(Q2h,2) -
		        Q2e*pow(l1k1,2)*pow(Q2h,2) + Q2k*pow(l1k1,2)*pow(Q2h,2) -
		        2*pow(l1k1,3)*pow(Q2h,2) + 4*l1k1*m4*pow(Q2k,2) -
		        4*m6*pow(Q2k,2) - 2*l1k1*m2*Q2h*pow(Q2k,2) +
		        2*m4*Q2h*pow(Q2k,2) + 2*m2*pow(l1k1,2)*pow(Q2k,2) -
		        Q2h*pow(l1k1,2)*pow(Q2k,2)) +
		     pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (-8*m4*Q2e*Q2h - 4*k1k2*m2*Q2e*Q2k + 4*k1k2*m2*Q2h*Q2k -
		        4*m2*Q2e*Q2h*Q2k + 4*k1k2*m2*pow(Q2e,2) + 4*m2*Q2h*pow(Q2e,2) -
		        2*m2*Q2k*pow(Q2e,2) + Q2h*pow(Q2e,3) - 4*k1k2*m2*pow(Q2h,2) +
		        8*m4*pow(Q2h,2) - 2*m2*Q2e*pow(Q2h,2) + 6*m2*Q2k*pow(Q2h,2) -
		        4*Q2e*Q2k*pow(Q2h,2) - pow(Q2e,2)*pow(Q2h,2) - 2*m2*pow(Q2h,3) +
		        4*Q2e*pow(Q2h,3) + 4*Q2k*pow(Q2h,3) - 4*pow(Q2h,4) +
		        2*m2*Q2e*pow(Q2k,2) - 2*m2*Q2h*pow(Q2k,2) + Q2e*Q2h*pow(Q2k,2) -
		        pow(Q2h,2)*pow(Q2k,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		      (64*k1k2*l1k2*m4 - 64*l1k2*m6 - 32*k1k2*l1k2*m2*Q2e +
		        32*k1k2*m4*Q2e - 32*m6*Q2e + 32*k1k2*l1k2*m2*Q2h -
		        32*k1k2*m4*Q2h + 32*l1k2*m4*Q2h + 32*m6*Q2h +
		        16*k1k2*l1k2*Q2e*Q2h - 32*l1k2*m2*Q2e*Q2h + 16*m4*Q2e*Q2h -
		        32*k1k2*l1k2*m2*Q2k - 16*k1k2*m2*Q2e*Q2k + 32*l1k2*m2*Q2e*Q2k +
		        16*k1k2*l1k2*Q2h*Q2k + 16*k1k2*m2*Q2h*Q2k - 64*l1k2*m2*Q2h*Q2k +
		        8*k1k2*Q2e*Q2h*Q2k - 16*l1k2*Q2e*Q2h*Q2k - 16*m2*Q2e*Q2h*Q2k +
		        32*l1k2*m2*pow(k1k2,2) + 16*m2*Q2e*pow(k1k2,2) -
		        16*l1k2*Q2h*pow(k1k2,2) - 16*m2*Q2h*pow(k1k2,2) -
		        8*Q2e*Q2h*pow(k1k2,2) - 64*k1k2*m2*pow(l1k2,2) +
		        32*m2*Q2e*pow(l1k2,2) + 32*k1k2*Q2h*pow(l1k2,2) -
		        96*m2*Q2h*pow(l1k2,2) - 16*Q2e*Q2h*pow(l1k2,2) +
		        64*m2*Q2k*pow(l1k2,2) - 32*Q2h*Q2k*pow(l1k2,2) +
		        64*m2*pow(l1k2,3) - 32*Q2h*pow(l1k2,3) - 32*k1k2*l1k2*pow(Q2h,2) +
		        48*l1k2*m2*pow(Q2h,2) - 16*m4*pow(Q2h,2) -
		        8*k1k2*Q2e*pow(Q2h,2) + 16*l1k2*Q2e*pow(Q2h,2) +
		        8*m2*Q2e*pow(Q2h,2) - 8*k1k2*Q2k*pow(Q2h,2) +
		        32*l1k2*Q2k*pow(Q2h,2) + 16*m2*Q2k*pow(Q2h,2) +
		        8*Q2e*Q2k*pow(Q2h,2) + 8*pow(k1k2,2)*pow(Q2h,2) +
		        48*pow(l1k2,2)*pow(Q2h,2) + 8*k1k2*pow(Q2h,3) -
		        24*l1k2*pow(Q2h,3) - 8*m2*pow(Q2h,3) - 4*Q2e*pow(Q2h,3) -
		        8*Q2k*pow(Q2h,3) + 4*pow(Q2h,4) + 16*l1k2*m2*pow(Q2k,2) +
		        8*m2*Q2e*pow(Q2k,2) - 8*l1k2*Q2h*pow(Q2k,2) -
		        8*m2*Q2h*pow(Q2k,2) - 4*Q2e*Q2h*pow(Q2k,2) +
		        4*pow(Q2h,2)*pow(Q2k,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (-32*l1k2*m4*Q2e + 32*l1k2*m4*Q2h - 80*l1k2*m2*Q2e*Q2h +
		        48*m4*Q2e*Q2h + 48*l1k2*m2*Q2e*Q2k - 16*m4*Q2e*Q2k -
		        80*l1k2*m2*Q2h*Q2k + 16*m4*Q2h*Q2k - 96*l1k2*Q2e*Q2h*Q2k -
		        32*m2*Q2e*Q2h*Q2k + 64*m2*Q2e*pow(l1k2,2) -
		        128*m2*Q2h*pow(l1k2,2) - 128*Q2e*Q2h*pow(l1k2,2) +
		        64*m2*Q2k*pow(l1k2,2) + 32*Q2e*Q2k*pow(l1k2,2) -
		        64*Q2h*Q2k*pow(l1k2,2) + 64*m2*pow(l1k2,3) + 32*Q2e*pow(l1k2,3) -
		        64*Q2h*pow(l1k2,3) + 8*l1k2*m2*pow(Q2e,2) - 16*m4*pow(Q2e,2) -
		        60*l1k2*Q2h*pow(Q2e,2) + 24*l1k2*Q2k*pow(Q2e,2) +
		        4*m2*Q2k*pow(Q2e,2) - 30*Q2h*Q2k*pow(Q2e,2) +
		        32*pow(l1k2,2)*pow(Q2e,2) + 12*l1k2*pow(Q2e,3) - 4*m2*pow(Q2e,3) -
		        10*Q2h*pow(Q2e,3) + 6*Q2k*pow(Q2e,3) + 2*pow(Q2e,4) +
		        88*l1k2*m2*pow(Q2h,2) - 32*m4*pow(Q2h,2) +
		        136*l1k2*Q2e*pow(Q2h,2) + 28*m2*Q2e*pow(Q2h,2) +
		        88*l1k2*Q2k*pow(Q2h,2) + 28*m2*Q2k*pow(Q2h,2) +
		        64*Q2e*Q2k*pow(Q2h,2) + 128*pow(l1k2,2)*pow(Q2h,2) +
		        32*pow(Q2e,2)*pow(Q2h,2) - 96*l1k2*pow(Q2h,3) - 24*m2*pow(Q2h,3) -
		        56*Q2e*pow(Q2h,3) - 40*Q2k*pow(Q2h,3) + 32*pow(Q2h,4) +
		        16*l1k2*m2*pow(Q2k,2) + 12*l1k2*Q2e*pow(Q2k,2) +
		        8*m2*Q2e*pow(Q2k,2) - 20*l1k2*Q2h*pow(Q2k,2) -
		        8*m2*Q2h*pow(Q2k,2) - 22*Q2e*Q2h*pow(Q2k,2) +
		        6*pow(Q2e,2)*pow(Q2k,2) + 16*pow(Q2h,2)*pow(Q2k,2) +
		        2*Q2e*pow(Q2k,3) - 2*Q2h*pow(Q2k,3)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
		      (-16*l1k1*m6 - 8*m6*Q2e + 8*l1k1*m4*Q2h + 10*l1k1*m2*Q2e*Q2h -
		        4*m4*Q2e*Q2h + 8*m6*Q2k - 20*l1k1*m2*Q2e*Q2k + 4*m4*Q2e*Q2k -
		        16*l1k1*m2*Q2h*Q2k - 5*l1k1*Q2e*Q2h*Q2k - 2*m2*Q2e*Q2h*Q2k +
		        16*m4*pow(l1k1,2) + 36*m2*Q2e*pow(l1k1,2) +
		        32*m2*Q2h*pow(l1k1,2) - 68*m2*Q2k*pow(l1k1,2) -
		        8*Q2e*Q2k*pow(l1k1,2) + 8*Q2h*Q2k*pow(l1k1,2) +
		        72*m2*pow(l1k1,3) + 8*Q2e*pow(l1k1,3) - 8*Q2h*pow(l1k1,3) +
		        l1k1*Q2h*pow(Q2e,2) + m2*Q2h*pow(Q2e,2) - 2*l1k1*Q2k*pow(Q2e,2) +
		        (Q2h*Q2k*pow(Q2e,2))/2. + 4*pow(l1k1,2)*pow(Q2e,2) +
		        l1k1*pow(Q2e,3) + Q2h*pow(Q2e,3) - (Q2k*pow(Q2e,3))/2. +
		        10*l1k1*m2*pow(Q2h,2) + 4*m4*pow(Q2h,2) + 3*l1k1*Q2e*pow(Q2h,2) -
		        2*m2*Q2e*pow(Q2h,2) + 9*l1k1*Q2k*pow(Q2h,2) -
		        2*m2*Q2k*pow(Q2h,2) - 6*Q2e*Q2k*pow(Q2h,2) -
		        8*pow(l1k1,2)*pow(Q2h,2) - pow(Q2e,2)*pow(Q2h,2) -
		        7*l1k1*pow(Q2h,3) + 3*m2*pow(Q2h,3) + 4*Q2e*pow(Q2h,3) +
		        6*Q2k*pow(Q2h,3) - 4*pow(Q2h,4) + 16*l1k1*m2*pow(Q2k,2) -
		        4*m4*pow(Q2k,2) + 3*l1k1*Q2e*pow(Q2k,2) + 2*m2*Q2e*pow(Q2k,2) -
		        3*l1k1*Q2h*pow(Q2k,2) + 3*Q2e*Q2h*pow(Q2k,2) -
		        3*pow(Q2h,2)*pow(Q2k,2) - (Q2e*pow(Q2k,3))/2. + (Q2h*pow(Q2k,3))/2.) \
		+ pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
		      pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (-16*m6*Q2e*Q2h + 8*m6*pow(Q2e,2) - 4*m4*Q2h*pow(Q2e,2) -
		        4*m4*Q2k*pow(Q2e,2) + 3*m2*Q2h*Q2k*pow(Q2e,2) + 4*m4*pow(Q2e,3) -
		        2*m2*Q2h*pow(Q2e,3) - m2*Q2k*pow(Q2e,3) - Q2h*Q2k*pow(Q2e,3) -
		        Q2h*pow(Q2e,4) + (Q2k*pow(Q2e,4))/2. + 8*m6*pow(Q2h,2) +
		        4*m4*Q2k*pow(Q2h,2) + m2*Q2e*Q2k*pow(Q2h,2) +
		        5*m2*pow(Q2e,2)*pow(Q2h,2) + (13*Q2k*pow(Q2e,2)*pow(Q2h,2))/2. +
		        2*pow(Q2e,3)*pow(Q2h,2) - 6*m2*Q2e*pow(Q2h,3) -
		        3*m2*Q2k*pow(Q2h,3) - 12*Q2e*Q2k*pow(Q2h,3) -
		        5*pow(Q2e,2)*pow(Q2h,3) + 3*m2*pow(Q2h,4) + 8*Q2e*pow(Q2h,4) +
		        6*Q2k*pow(Q2h,4) - 4*pow(Q2h,5) + 4*m4*Q2e*pow(Q2k,2) -
		        4*m4*Q2h*pow(Q2k,2) - m2*pow(Q2e,2)*pow(Q2k,2) -
		        3*Q2h*pow(Q2e,2)*pow(Q2k,2) + m2*pow(Q2h,2)*pow(Q2k,2) +
		        6*Q2e*pow(Q2h,2)*pow(Q2k,2) - 3*pow(Q2h,3)*pow(Q2k,2) -
		        Q2e*Q2h*pow(Q2k,3) + (pow(Q2e,2)*pow(Q2k,3))/2. +
		        (pow(Q2h,2)*pow(Q2k,3))/2.)) +
		  (G2 + G3)*((pow(k1k2 - l1k1 - l1k2,-2)*
		        (-4*l1k1*l1k2*M2*Q2e*Q2h - 8*l1k1*m2*M2*Q2e*Q2h +
		          8*l1k1*l1k2*M2*Q2h*Q2k + 16*l1k1*m2*M2*Q2h*Q2k +
		          2*l1k2*M2*Q2e*Q2h*Q2k + 4*m2*M2*Q2e*Q2h*Q2k +
		          8*l1k1*l1k2*Q2e*Q2h*S + 16*l1k1*m2*Q2e*Q2h*S -
		          16*l1k1*l1k2*Q2h*Q2k*S - 32*l1k1*m2*Q2h*Q2k*S -
		          4*l1k2*Q2e*Q2h*Q2k*S - 8*m2*Q2e*Q2h*Q2k*S -
		          8*l1k2*M2*Q2h*pow(l1k1,2) - 16*m2*M2*Q2h*pow(l1k1,2) +
		          16*l1k2*Q2h*S*pow(l1k1,2) + 32*m2*Q2h*S*pow(l1k1,2) -
		          4*l1k1*l1k2*M2*pow(Q2h,2) - 8*l1k1*m2*M2*pow(Q2h,2) -
		          2*l1k2*M2*Q2e*pow(Q2h,2) - 4*m2*M2*Q2e*pow(Q2h,2) +
		          2*l1k2*M2*Q2k*pow(Q2h,2) + 4*m2*M2*Q2k*pow(Q2h,2) +
		          8*l1k1*l1k2*S*pow(Q2h,2) + 16*l1k1*m2*S*pow(Q2h,2) +
		          4*l1k2*Q2e*S*pow(Q2h,2) + 8*m2*Q2e*S*pow(Q2h,2) -
		          4*l1k2*Q2k*S*pow(Q2h,2) - 8*m2*Q2k*S*pow(Q2h,2) -
		          2*l1k2*M2*Q2h*pow(Q2k,2) - 4*m2*M2*Q2h*pow(Q2k,2) +
		          4*l1k2*Q2h*S*pow(Q2k,2) + 8*m2*Q2h*S*pow(Q2k,2) +
		          16*l1k1*l1k2*Q2e*pow(S,2) + 32*l1k1*m2*Q2e*pow(S,2) +
		          16*l1k1*l1k2*Q2h*pow(S,2) + 32*l1k1*m2*Q2h*pow(S,2) +
		          8*l1k2*Q2e*Q2h*pow(S,2) + 16*m2*Q2e*Q2h*pow(S,2) -
		          32*l1k1*l1k2*Q2k*pow(S,2) - 64*l1k1*m2*Q2k*pow(S,2) -
		          8*l1k2*Q2e*Q2k*pow(S,2) - 16*m2*Q2e*Q2k*pow(S,2) -
		          8*l1k2*Q2h*Q2k*pow(S,2) - 16*m2*Q2h*Q2k*pow(S,2) +
		          32*l1k2*pow(l1k1,2)*pow(S,2) + 64*m2*pow(l1k1,2)*pow(S,2) +
		          8*l1k2*pow(Q2k,2)*pow(S,2) + 16*m2*pow(Q2k,2)*pow(S,2)))/2. +
		     (pow(k1k2 - l1k1 - l1k2,-2)*pow(l1k2,-1)*
		        (8*l1k1*m4*M2*Q2e*Q2h - 16*l1k1*m4*M2*Q2h*Q2k +
		          4*l1k1*m2*M2*Q2e*Q2h*Q2k - 4*m4*M2*Q2e*Q2h*Q2k -
		          16*l1k1*m4*Q2e*Q2h*S + 32*l1k1*m4*Q2h*Q2k*S -
		          8*l1k1*m2*Q2e*Q2h*Q2k*S + 8*m4*Q2e*Q2h*Q2k*S +
		          16*m4*M2*Q2h*pow(l1k1,2) - 8*m2*M2*Q2e*Q2h*pow(l1k1,2) +
		          16*m2*M2*Q2h*Q2k*pow(l1k1,2) +
		          2*M2*Q2e*Q2h*Q2k*pow(l1k1,2) - 32*m4*Q2h*S*pow(l1k1,2) +
		          16*m2*Q2e*Q2h*S*pow(l1k1,2) - 32*m2*Q2h*Q2k*S*pow(l1k1,2) -
		          4*Q2e*Q2h*Q2k*S*pow(l1k1,2) - 16*m2*M2*Q2h*pow(l1k1,3) -
		          4*M2*Q2e*Q2h*pow(l1k1,3) + 8*M2*Q2h*Q2k*pow(l1k1,3) +
		          32*m2*Q2h*S*pow(l1k1,3) + 8*Q2e*Q2h*S*pow(l1k1,3) -
		          16*Q2h*Q2k*S*pow(l1k1,3) - 8*M2*Q2h*pow(l1k1,4) +
		          16*Q2h*S*pow(l1k1,4) + 8*l1k1*m4*M2*pow(Q2h,2) -
		          4*l1k1*m2*M2*Q2e*pow(Q2h,2) + 4*m4*M2*Q2e*pow(Q2h,2) +
		          4*l1k1*m2*M2*Q2k*pow(Q2h,2) - 4*m4*M2*Q2k*pow(Q2h,2) -
		          16*l1k1*m4*S*pow(Q2h,2) + 8*l1k1*m2*Q2e*S*pow(Q2h,2) -
		          8*m4*Q2e*S*pow(Q2h,2) - 8*l1k1*m2*Q2k*S*pow(Q2h,2) +
		          8*m4*Q2k*S*pow(Q2h,2) - 8*m2*M2*pow(l1k1,2)*pow(Q2h,2) -
		          2*M2*Q2e*pow(l1k1,2)*pow(Q2h,2) +
		          2*M2*Q2k*pow(l1k1,2)*pow(Q2h,2) +
		          16*m2*S*pow(l1k1,2)*pow(Q2h,2) +
		          4*Q2e*S*pow(l1k1,2)*pow(Q2h,2) - 4*Q2k*S*pow(l1k1,2)*pow(Q2h,2) -
		          4*M2*pow(l1k1,3)*pow(Q2h,2) + 8*S*pow(l1k1,3)*pow(Q2h,2) -
		          4*l1k1*m2*M2*Q2h*pow(Q2k,2) + 4*m4*M2*Q2h*pow(Q2k,2) +
		          8*l1k1*m2*Q2h*S*pow(Q2k,2) - 8*m4*Q2h*S*pow(Q2k,2) -
		          2*M2*Q2h*pow(l1k1,2)*pow(Q2k,2) +
		          4*Q2h*S*pow(l1k1,2)*pow(Q2k,2) - 32*l1k1*m4*Q2e*pow(S,2) -
		          32*l1k1*m4*Q2h*pow(S,2) + 16*l1k1*m2*Q2e*Q2h*pow(S,2) -
		          16*m4*Q2e*Q2h*pow(S,2) + 64*l1k1*m4*Q2k*pow(S,2) -
		          16*l1k1*m2*Q2e*Q2k*pow(S,2) + 16*m4*Q2e*Q2k*pow(S,2) -
		          16*l1k1*m2*Q2h*Q2k*pow(S,2) + 16*m4*Q2h*Q2k*pow(S,2) -
		          64*m4*pow(l1k1,2)*pow(S,2) + 32*m2*Q2e*pow(l1k1,2)*pow(S,2) +
		          32*m2*Q2h*pow(l1k1,2)*pow(S,2) +
		          8*Q2e*Q2h*pow(l1k1,2)*pow(S,2) -
		          64*m2*Q2k*pow(l1k1,2)*pow(S,2) -
		          8*Q2e*Q2k*pow(l1k1,2)*pow(S,2) - 8*Q2h*Q2k*pow(l1k1,2)*pow(S,2) +
		          64*m2*pow(l1k1,3)*pow(S,2) + 16*Q2e*pow(l1k1,3)*pow(S,2) +
		          16*Q2h*pow(l1k1,3)*pow(S,2) - 32*Q2k*pow(l1k1,3)*pow(S,2) +
		          32*pow(l1k1,4)*pow(S,2) + 16*l1k1*m2*pow(Q2k,2)*pow(S,2) -
		          16*m4*pow(Q2k,2)*pow(S,2) + 8*pow(l1k1,2)*pow(Q2k,2)*pow(S,2)))/2. \
		+ (48*k1k2*m2*M2 + 16*l1k1*m2*M2 - 80*l1k2*m2*M2 + 32*m4*M2 -
		        8*l1k1*M2*Q2e - 32*l1k2*M2*Q2e - 8*m2*M2*Q2e + 16*m4*Q2h +
		        16*k1k2*M2*Q2h + 24*l1k1*M2*Q2h - 48*l1k2*M2*Q2h +
		        64*m2*M2*Q2h + 16*M2*Q2e*Q2h - 24*m2*M2*Q2k -
		        4*M2*Q2e*Q2k + 8*m2*Q2h*Q2k - 28*M2*Q2h*Q2k + 4*Q2e*Q2h*Q2k +
		        16*k1k2*Q2h*S + 16*l1k1*Q2h*S - 32*l1k2*Q2h*S - 32*m2*Q2h*S +
		        32*m2*Q2k*S - 24*Q2h*Q2k*S - 32*m4*Sk + 16*m2*Q2h*Sk +
		        8*Q2e*Q2h*Sk - 16*m2*Q2k*Sk + 64*m2*S*Sk - 48*m2*Q2h*Sq2 -
		        8*Q2e*Q2h*Sq2 + 32*m2*Q2k*Sq2 - 8*Q2h*Q2k*Sq2 - 192*k1k2*S*Sq2 -
		        128*l1k1*S*Sq2 + 256*l1k2*S*Sq2 - 256*m2*S*Sq2 + 32*Q2e*S*Sq2 -
		        256*Q2h*S*Sq2 + 128*Q2k*S*Sq2 + 64*m2*Sk*Sq2 -
		        16*M2*pow(Q2e,2) + 16*l1k1*pow(Q2h,2) - 16*l1k2*pow(Q2h,2) +
		        8*m2*pow(Q2h,2) + 48*M2*pow(Q2h,2) - 12*Q2k*pow(Q2h,2) +
		        8*S*pow(Q2h,2) - 8*Sk*pow(Q2h,2) + 8*Sq2*pow(Q2h,2) +
		        12*pow(Q2h,3) - 128*k1k2*pow(S,2) - 128*l1k1*pow(S,2) +
		        192*l1k2*pow(S,2) - 192*m2*pow(S,2) + 32*Q2e*pow(S,2) -
		        192*Q2h*pow(S,2) + 96*Q2k*pow(S,2) - 32*m2*pow(Sk,2) -
		        64*k1k2*pow(Sq2,2) + 64*l1k2*pow(Sq2,2) - 64*m2*pow(Sq2,2) -
		        64*Q2h*pow(Sq2,2) + 32*Q2k*pow(Sq2,2))/2. +
		     (pow(k1k2 - l1k1 - l1k2,-1)*
		        (16*l1k1*l1k2*m2*M2 + 32*l1k1*m4*M2 - 8*l1k1*l1k2*M2*Q2e +
		          32*l1k1*m2*M2*Q2e + 8*l1k2*m2*M2*Q2e + 16*m4*M2*Q2e +
		          16*l1k1*m4*Q2h - 40*l1k1*l1k2*M2*Q2h + 32*l1k1*m2*M2*Q2h +
		          8*m4*Q2e*Q2h + 12*l1k1*M2*Q2e*Q2h + 8*l1k2*M2*Q2e*Q2h -
		          8*m2*M2*Q2e*Q2h - 32*l1k1*m2*M2*Q2k - 8*l1k2*m2*M2*Q2k -
		          16*m4*M2*Q2k - 12*l1k1*M2*Q2e*Q2k + 4*l1k2*M2*Q2e*Q2k -
		          16*l1k1*m2*Q2h*Q2k - 8*m4*Q2h*Q2k - 20*l1k1*M2*Q2h*Q2k +
		          36*l1k2*M2*Q2h*Q2k - 24*m2*M2*Q2h*Q2k + 4*l1k1*Q2e*Q2h*Q2k -
		          10*m2*Q2e*Q2h*Q2k + 12*M2*Q2e*Q2h*Q2k + 16*l1k1*l1k2*Q2h*S +
		          32*l1k1*m2*Q2h*S - 8*l1k1*Q2e*Q2h*S + 8*l1k2*Q2e*Q2h*S +
		          24*m2*Q2e*Q2h*S - 64*l1k1*m2*Q2k*S + 24*l1k1*Q2e*Q2k*S -
		          16*m2*Q2e*Q2k*S - 24*l1k1*Q2h*Q2k*S - 8*l1k2*Q2h*Q2k*S -
		          32*m2*Q2h*Q2k*S + 8*Q2e*Q2h*Q2k*S - 32*l1k1*m4*Sk -
		          16*m4*Q2e*Sk - 32*l1k1*m2*Q2h*Sk + 8*l1k1*Q2e*Q2h*Sk -
		          20*m2*Q2e*Q2h*Sk - 16*l1k1*m2*Q2k*Sk + 16*m4*Q2k*Sk +
		          16*m2*Q2e*Q2k*Sk - 8*m2*Q2h*Q2k*Sk + 4*Q2e*Q2h*Q2k*Sk -
		          128*l1k1*m2*S*Sk + 48*l1k1*Q2e*S*Sk - 32*m2*Q2e*S*Sk -
		          16*l1k1*Q2h*S*Sk - 32*m2*Q2h*S*Sk + 8*Q2e*Q2h*S*Sk +
		          64*m2*Q2k*S*Sk - 24*Q2e*Q2k*S*Sk + 8*Q2h*Q2k*S*Sk +
		          32*l1k1*l1k2*Q2h*Sq2 + 16*l1k1*m2*Q2h*Sq2 -
		          16*l1k2*Q2e*Q2h*Sq2 + 24*m2*Q2e*Q2h*Sq2 + 16*l1k1*Q2h*Q2k*Sq2 -
		          32*l1k2*Q2h*Q2k*Sq2 - 8*m2*Q2h*Q2k*Sq2 - 12*Q2e*Q2h*Q2k*Sq2 +
		          128*l1k1*l1k2*S*Sq2 + 32*l1k1*Q2e*S*Sq2 - 64*l1k2*Q2e*S*Sq2 +
		          64*m2*Q2e*S*Sq2 - 128*l1k1*Q2h*S*Sq2 + 192*l1k2*Q2h*S*Sq2 -
		          64*m2*Q2h*S*Sq2 + 128*Q2e*Q2h*S*Sq2 + 64*l1k1*Q2k*S*Sq2 -
		          128*l1k2*Q2k*S*Sq2 - 48*Q2e*Q2k*S*Sq2 + 96*Q2h*Q2k*S*Sq2 -
		          16*Q2e*Q2h*Sk*Sq2 + 48*m2*M2*pow(l1k1,2) +
		          16*M2*Q2e*pow(l1k1,2) + 16*M2*Q2h*pow(l1k1,2) +
		          16*Q2h*S*pow(l1k1,2) - 16*Q2h*Sq2*pow(l1k1,2) -
		          64*S*Sq2*pow(l1k1,2) + 32*M2*Q2h*pow(l1k2,2) -
		          32*Q2h*Sq2*pow(l1k2,2) - 128*S*Sq2*pow(l1k2,2) +
		          4*m2*M2*pow(Q2e,2) + 6*M2*Q2h*pow(Q2e,2) -
		          4*Q2h*S*pow(Q2e,2) + 4*Q2k*S*pow(Q2e,2) + 8*S*Sk*pow(Q2e,2) -
		          4*Q2h*Sq2*pow(Q2e,2) - 16*S*Sq2*pow(Q2e,2) -
		          16*l1k1*l1k2*pow(Q2h,2) + 16*l1k1*m2*pow(Q2h,2) +
		          28*l1k1*M2*pow(Q2h,2) - 48*l1k2*M2*pow(Q2h,2) +
		          24*m2*M2*pow(Q2h,2) + 8*l1k2*Q2e*pow(Q2h,2) -
		          34*M2*Q2e*pow(Q2h,2) - 16*l1k1*Q2k*pow(Q2h,2) +
		          16*l1k2*Q2k*pow(Q2h,2) - 6*m2*Q2k*pow(Q2h,2) -
		          32*M2*Q2k*pow(Q2h,2) + 8*Q2e*Q2k*pow(Q2h,2) +
		          8*l1k1*S*pow(Q2h,2) - 8*m2*S*pow(Q2h,2) - 4*Q2k*S*pow(Q2h,2) +
		          4*m2*Sk*pow(Q2h,2) + 4*Q2e*Sk*pow(Q2h,2) - 8*Q2k*Sk*pow(Q2h,2) -
		          24*l1k1*Sq2*pow(Q2h,2) + 48*l1k2*Sq2*pow(Q2h,2) -
		          16*m2*Sq2*pow(Q2h,2) + 4*Q2e*Sq2*pow(Q2h,2) +
		          24*Q2k*Sq2*pow(Q2h,2) - 144*S*Sq2*pow(Q2h,2) +
		          16*Sk*Sq2*pow(Q2h,2) + 16*pow(l1k1,2)*pow(Q2h,2) +
		          16*pow(l1k2,2)*pow(Q2h,2) + 2*pow(Q2e,2)*pow(Q2h,2) +
		          16*l1k1*pow(Q2h,3) - 24*l1k2*pow(Q2h,3) + 8*m2*pow(Q2h,3) +
		          40*M2*pow(Q2h,3) - 10*Q2e*pow(Q2h,3) - 14*Q2k*pow(Q2h,3) -
		          8*Sq2*pow(Q2h,3) + 12*pow(Q2h,4) + 4*m2*M2*pow(Q2k,2) +
		          6*m2*Q2e*pow(Q2k,2) - 2*M2*Q2e*pow(Q2k,2) +
		          2*m2*Q2h*pow(Q2k,2) + 10*M2*Q2h*pow(Q2k,2) -
		          2*Q2e*Q2h*pow(Q2k,2) + 32*m2*S*pow(Q2k,2) -
		          12*Q2e*S*pow(Q2k,2) + 8*Q2h*S*pow(Q2k,2) + 8*m2*Sk*pow(Q2k,2) -
		          4*Q2h*Sq2*pow(Q2k,2) - 16*S*Sq2*pow(Q2k,2) +
		          4*pow(Q2h,2)*pow(Q2k,2) + 128*l1k1*l1k2*pow(S,2) -
		          64*l1k1*m2*pow(S,2) + 32*l1k1*Q2e*pow(S,2) -
		          64*l1k2*Q2e*pow(S,2) + 64*m2*Q2e*pow(S,2) -
		          128*l1k1*Q2h*pow(S,2) + 192*l1k2*Q2h*pow(S,2) -
		          96*m2*Q2h*pow(S,2) + 128*Q2e*Q2h*pow(S,2) +
		          64*l1k1*Q2k*pow(S,2) - 128*l1k2*Q2k*pow(S,2) +
		          32*m2*Q2k*pow(S,2) - 64*Q2e*Q2k*pow(S,2) +
		          112*Q2h*Q2k*pow(S,2) - 64*pow(l1k1,2)*pow(S,2) -
		          128*pow(l1k2,2)*pow(S,2) - 16*pow(Q2e,2)*pow(S,2) -
		          144*pow(Q2h,2)*pow(S,2) - 16*pow(Q2k,2)*pow(S,2) -
		          32*l1k1*m2*pow(Sk,2) + 8*m2*Q2e*pow(Sk,2) -
		          24*m2*Q2h*pow(Sk,2) + 16*Q2e*Q2h*pow(Sk,2) +
		          16*m2*Q2k*pow(Sk,2) - 16*pow(Q2h,2)*pow(Sk,2) +
		          32*Q2e*Q2h*pow(Sq2,2) - 32*pow(Q2h,2)*pow(Sq2,2)))/2. +
		     (pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		        (32*l1k2*m2*M2*Q2e*Q2h - 64*l1k2*M2*Q2e*Q2h*Q2k +
		          16*m2*M2*Q2e*Q2h*Q2k - 16*l1k2*Q2e*Q2h*Q2k*Sk -
		          32*l1k2*m2*Q2e*Q2h*Sq2 + 56*l1k2*Q2e*Q2h*Q2k*Sq2 -
		          16*m2*Q2e*Q2h*Q2k*Sq2 - 128*l1k2*m2*Q2e*S*Sq2 +
		          128*l1k2*m2*Q2h*S*Sq2 - 512*l1k2*Q2e*Q2h*S*Sq2 +
		          192*m2*Q2e*Q2h*S*Sq2 + 192*l1k2*Q2e*Q2k*S*Sq2 -
		          64*m2*Q2e*Q2k*S*Sq2 - 320*l1k2*Q2h*Q2k*S*Sq2 +
		          64*m2*Q2h*Q2k*S*Sq2 - 224*Q2e*Q2h*Q2k*S*Sq2 +
		          32*l1k2*Q2e*Q2h*Sk*Sq2 + 16*Q2e*Q2h*Q2k*Sk*Sq2 -
		          64*M2*Q2e*Q2h*pow(l1k2,2) - 64*M2*Q2h*Q2k*pow(l1k2,2) +
		          64*Q2e*Q2h*Sq2*pow(l1k2,2) + 64*Q2h*Q2k*Sq2*pow(l1k2,2) +
		          256*Q2e*S*Sq2*pow(l1k2,2) - 512*Q2h*S*Sq2*pow(l1k2,2) +
		          256*Q2k*S*Sq2*pow(l1k2,2) - 64*M2*Q2h*pow(l1k2,3) +
		          64*Q2h*Sq2*pow(l1k2,3) + 256*S*Sq2*pow(l1k2,3) -
		          24*l1k2*M2*Q2h*pow(Q2e,2) + 16*m2*M2*Q2h*pow(Q2e,2) -
		          20*M2*Q2h*Q2k*pow(Q2e,2) + 16*Q2h*Q2k*S*pow(Q2e,2) -
		          8*Q2h*Q2k*Sk*pow(Q2e,2) + 64*Q2h*S*Sk*pow(Q2e,2) +
		          24*l1k2*Q2h*Sq2*pow(Q2e,2) - 16*m2*Q2h*Sq2*pow(Q2e,2) +
		          24*Q2h*Q2k*Sq2*pow(Q2e,2) + 96*l1k2*S*Sq2*pow(Q2e,2) -
		          64*m2*S*Sq2*pow(Q2e,2) - 192*Q2h*S*Sq2*pow(Q2e,2) +
		          48*Q2k*S*Sq2*pow(Q2e,2) + 48*Q2h*Sk*Sq2*pow(Q2e,2) -
		          4*M2*Q2h*pow(Q2e,3) + 4*Q2h*Sq2*pow(Q2e,3) +
		          16*S*Sq2*pow(Q2e,3) - 32*l1k2*m2*M2*pow(Q2h,2) +
		          16*l1k2*m2*Q2e*pow(Q2h,2) + 136*l1k2*M2*Q2e*pow(Q2h,2) -
		          48*m2*M2*Q2e*pow(Q2h,2) + 96*l1k2*M2*Q2k*pow(Q2h,2) -
		          16*m2*M2*Q2k*pow(Q2h,2) - 28*l1k2*Q2e*Q2k*pow(Q2h,2) +
		          8*m2*Q2e*Q2k*pow(Q2h,2) + 84*M2*Q2e*Q2k*pow(Q2h,2) -
		          32*Q2e*Q2k*S*pow(Q2h,2) + 16*l1k2*Q2k*Sk*pow(Q2h,2) +
		          24*Q2e*Q2k*Sk*pow(Q2h,2) - 128*Q2e*S*Sk*pow(Q2h,2) +
		          32*l1k2*m2*Sq2*pow(Q2h,2) - 72*l1k2*Q2e*Sq2*pow(Q2h,2) +
		          48*m2*Q2e*Sq2*pow(Q2h,2) - 88*l1k2*Q2k*Sq2*pow(Q2h,2) +
		          16*m2*Q2k*Sq2*pow(Q2h,2) - 56*Q2e*Q2k*Sq2*pow(Q2h,2) +
		          480*l1k2*S*Sq2*pow(Q2h,2) - 128*m2*S*Sq2*pow(Q2h,2) +
		          432*Q2e*S*Sq2*pow(Q2h,2) + 176*Q2k*S*Sq2*pow(Q2h,2) -
		          32*l1k2*Sk*Sq2*pow(Q2h,2) - 112*Q2e*Sk*Sq2*pow(Q2h,2) -
		          16*Q2k*Sk*Sq2*pow(Q2h,2) + 128*M2*pow(l1k2,2)*pow(Q2h,2) -
		          32*Q2e*pow(l1k2,2)*pow(Q2h,2) - 32*Q2k*pow(l1k2,2)*pow(Q2h,2) -
		          128*Sq2*pow(l1k2,2)*pow(Q2h,2) - 32*pow(l1k2,3)*pow(Q2h,2) -
		          12*l1k2*pow(Q2e,2)*pow(Q2h,2) + 8*m2*pow(Q2e,2)*pow(Q2h,2) +
		          44*M2*pow(Q2e,2)*pow(Q2h,2) - 8*Q2k*pow(Q2e,2)*pow(Q2h,2) -
		          12*Sq2*pow(Q2e,2)*pow(Q2h,2) - 2*pow(Q2e,3)*pow(Q2h,2) -
		          16*l1k2*m2*pow(Q2h,3) - 128*l1k2*M2*pow(Q2h,3) +
		          32*m2*M2*pow(Q2h,3) + 52*l1k2*Q2e*pow(Q2h,3) -
		          24*m2*Q2e*pow(Q2h,3) - 104*M2*Q2e*pow(Q2h,3) +
		          44*l1k2*Q2k*pow(Q2h,3) - 8*m2*Q2k*pow(Q2h,3) -
		          64*M2*Q2k*pow(Q2h,3) + 28*Q2e*Q2k*pow(Q2h,3) +
		          16*Q2k*S*pow(Q2h,3) - 16*Q2k*Sk*pow(Q2h,3) + 64*S*Sk*pow(Q2h,3) +
		          64*l1k2*Sq2*pow(Q2h,3) - 32*m2*Sq2*pow(Q2h,3) +
		          8*Q2e*Sq2*pow(Q2h,3) + 32*Q2k*Sq2*pow(Q2h,3) -
		          256*S*Sq2*pow(Q2h,3) + 64*Sk*Sq2*pow(Q2h,3) +
		          64*pow(l1k2,2)*pow(Q2h,3) + 14*pow(Q2e,2)*pow(Q2h,3) -
		          48*l1k2*pow(Q2h,4) + 16*m2*pow(Q2h,4) + 64*M2*pow(Q2h,4) -
		          28*Q2e*pow(Q2h,4) - 20*Q2k*pow(Q2h,4) + 16*pow(Q2h,5) +
		          8*l1k2*M2*Q2e*pow(Q2k,2) - 24*l1k2*M2*Q2h*pow(Q2k,2) -
		          28*M2*Q2e*Q2h*pow(Q2k,2) - 8*Q2e*Q2h*Sk*pow(Q2k,2) +
		          16*l1k2*Q2h*Sq2*pow(Q2k,2) + 12*Q2e*Q2h*Sq2*pow(Q2k,2) +
		          64*l1k2*S*Sq2*pow(Q2k,2) + 32*Q2e*S*Sq2*pow(Q2k,2) -
		          32*Q2h*S*Sq2*pow(Q2k,2) + 4*M2*pow(Q2e,2)*pow(Q2k,2) -
		          8*l1k2*pow(Q2h,2)*pow(Q2k,2) + 24*M2*pow(Q2h,2)*pow(Q2k,2) -
		          6*Q2e*pow(Q2h,2)*pow(Q2k,2) + 8*Sk*pow(Q2h,2)*pow(Q2k,2) -
		          12*Sq2*pow(Q2h,2)*pow(Q2k,2) + 6*pow(Q2h,3)*pow(Q2k,2) +
		          4*M2*Q2e*pow(Q2k,3) - 4*M2*Q2h*pow(Q2k,3) -
		          128*l1k2*m2*Q2e*pow(S,2) + 128*l1k2*m2*Q2h*pow(S,2) -
		          512*l1k2*Q2e*Q2h*pow(S,2) + 192*m2*Q2e*Q2h*pow(S,2) +
		          192*l1k2*Q2e*Q2k*pow(S,2) - 64*m2*Q2e*Q2k*pow(S,2) -
		          320*l1k2*Q2h*Q2k*pow(S,2) + 64*m2*Q2h*Q2k*pow(S,2) -
		          224*Q2e*Q2h*Q2k*pow(S,2) + 256*Q2e*pow(l1k2,2)*pow(S,2) -
		          512*Q2h*pow(l1k2,2)*pow(S,2) + 256*Q2k*pow(l1k2,2)*pow(S,2) +
		          256*pow(l1k2,3)*pow(S,2) + 96*l1k2*pow(Q2e,2)*pow(S,2) -
		          64*m2*pow(Q2e,2)*pow(S,2) - 160*Q2h*pow(Q2e,2)*pow(S,2) +
		          48*Q2k*pow(Q2e,2)*pow(S,2) + 16*pow(Q2e,3)*pow(S,2) +
		          480*l1k2*pow(Q2h,2)*pow(S,2) - 128*m2*pow(Q2h,2)*pow(S,2) +
		          368*Q2e*pow(Q2h,2)*pow(S,2) + 176*Q2k*pow(Q2h,2)*pow(S,2) -
		          224*pow(Q2h,3)*pow(S,2) + 64*l1k2*pow(Q2k,2)*pow(S,2) +
		          32*Q2e*pow(Q2k,2)*pow(S,2) - 32*Q2h*pow(Q2k,2)*pow(S,2) -
		          32*l1k2*Q2e*Q2h*pow(Sk,2) - 16*Q2e*Q2h*Q2k*pow(Sk,2) -
		          16*Q2h*pow(Q2e,2)*pow(Sk,2) + 32*l1k2*pow(Q2h,2)*pow(Sk,2) +
		          48*Q2e*pow(Q2h,2)*pow(Sk,2) + 16*Q2k*pow(Q2h,2)*pow(Sk,2) -
		          32*pow(Q2h,3)*pow(Sk,2) - 64*l1k2*Q2e*Q2h*pow(Sq2,2) -
		          32*Q2e*Q2h*Q2k*pow(Sq2,2) - 48*Q2h*pow(Q2e,2)*pow(Sq2,2) +
		          64*l1k2*pow(Q2h,2)*pow(Sq2,2) + 128*Q2e*pow(Q2h,2)*pow(Sq2,2) +
		          32*Q2k*pow(Q2h,2)*pow(Sq2,2) - 80*pow(Q2h,3)*pow(Sq2,2)))/2. +
		     (pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
		        (16*l1k1*m4*M2*Q2h + 12*l1k1*m2*M2*Q2e*Q2h +
		          8*m4*M2*Q2e*Q2h - 20*l1k1*m2*M2*Q2e*Q2k +
		          16*l1k1*m4*Q2h*Q2k - 12*l1k1*m2*M2*Q2h*Q2k -
		          8*m4*M2*Q2h*Q2k - 10*l1k1*m2*Q2e*Q2h*Q2k +
		          4*m4*Q2e*Q2h*Q2k - 16*l1k1*M2*Q2e*Q2h*Q2k -
		          64*l1k1*m4*Q2h*S + 40*l1k1*m2*Q2e*Q2h*S - 16*m4*Q2e*Q2h*S +
		          64*l1k1*m4*Q2k*S + 16*l1k1*m2*Q2e*Q2k*S + 16*m4*Q2e*Q2k*S -
		          64*l1k1*m2*Q2h*Q2k*S + 48*m4*Q2h*Q2k*S + 4*l1k1*Q2e*Q2h*Q2k*S +
		          8*m2*Q2e*Q2h*Q2k*S + 32*l1k1*m4*Q2h*Sk -
		          20*l1k1*m2*Q2e*Q2h*Sk + 8*m4*Q2e*Q2h*Sk +
		          24*l1k1*m2*Q2e*Q2k*Sk - 16*m4*Q2h*Q2k*Sk -
		          4*l1k1*Q2e*Q2h*Q2k*Sk - 8*m2*Q2e*Q2h*Q2k*Sk +
		          128*l1k1*m4*S*Sk + 32*l1k1*m2*Q2e*S*Sk + 32*m4*Q2e*S*Sk +
		          32*m4*Q2h*S*Sk + 40*l1k1*Q2e*Q2h*S*Sk + 16*m2*Q2e*Q2h*S*Sk +
		          32*l1k1*m2*Q2k*S*Sk - 64*m4*Q2k*S*Sk - 24*l1k1*Q2e*Q2k*S*Sk -
		          32*m2*Q2e*Q2k*S*Sk + 8*l1k1*Q2h*Q2k*S*Sk - 16*Q2e*Q2h*Q2k*S*Sk -
		          16*l1k1*m4*Q2h*Sq2 + 16*l1k1*m2*Q2e*Q2h*Sq2 -
		          8*m4*Q2e*Q2h*Sq2 - 32*l1k1*m2*Q2h*Q2k*Sq2 + 8*m4*Q2h*Q2k*Sq2 -
		          8*l1k1*Q2e*Q2h*Q2k*Sq2 + 18*m2*Q2e*Q2h*Q2k*Sq2 -
		          64*l1k1*m4*S*Sq2 - 32*m4*Q2e*S*Sq2 + 64*l1k1*m2*Q2h*S*Sq2 -
		          56*l1k1*Q2e*Q2h*S*Sq2 - 32*m2*Q2e*Q2h*S*Sq2 -
		          96*l1k1*m2*Q2k*S*Sq2 + 32*m4*Q2k*S*Sq2 - 8*l1k1*Q2e*Q2k*S*Sq2 -
		          24*l1k1*Q2h*Q2k*S*Sq2 - 16*m2*Q2h*Q2k*S*Sq2 +
		          24*Q2e*Q2h*Q2k*S*Sq2 + 32*l1k1*m2*Q2h*Sk*Sq2 +
		          24*l1k1*Q2e*Q2h*Sk*Sq2 + 16*m2*Q2e*Q2h*Sk*Sq2 -
		          24*m2*Q2e*Q2k*Sk*Sq2 + 8*m2*Q2h*Q2k*Sk*Sq2 -
		          12*Q2e*Q2h*Q2k*Sk*Sq2 + 40*m2*M2*Q2e*pow(l1k1,2) +
		          4*M2*Q2e*Q2h*pow(l1k1,2) - 56*m2*M2*Q2k*pow(l1k1,2) -
		          24*M2*Q2e*Q2k*pow(l1k1,2) - 24*m2*Q2h*Q2k*pow(l1k1,2) +
		          24*M2*Q2h*Q2k*pow(l1k1,2) + 96*m2*Q2h*S*pow(l1k1,2) +
		          24*Q2e*Q2h*S*pow(l1k1,2) - 32*m2*Q2k*S*pow(l1k1,2) +
		          24*Q2e*Q2k*S*pow(l1k1,2) - 32*Q2h*Q2k*S*pow(l1k1,2) -
		          48*m2*Q2h*Sk*pow(l1k1,2) - 64*m2*S*Sk*pow(l1k1,2) +
		          48*Q2e*S*Sk*pow(l1k1,2) - 16*Q2h*S*Sk*pow(l1k1,2) +
		          64*m2*Q2h*Sq2*pow(l1k1,2) + 24*Q2e*Q2h*Sq2*pow(l1k1,2) -
		          16*Q2h*Q2k*Sq2*pow(l1k1,2) + 128*m2*S*Sq2*pow(l1k1,2) +
		          64*Q2h*S*Sq2*pow(l1k1,2) - 64*Q2k*S*Sq2*pow(l1k1,2) +
		          48*m2*M2*pow(l1k1,3) + 24*M2*Q2e*pow(l1k1,3) -
		          24*M2*Q2h*pow(l1k1,3) + 32*Q2h*S*pow(l1k1,3) +
		          16*Q2h*Sq2*pow(l1k1,3) + 64*S*Sq2*pow(l1k1,3) +
		          8*l1k1*m2*M2*pow(Q2e,2) + 6*l1k1*M2*Q2h*pow(Q2e,2) -
		          2*m2*M2*Q2h*pow(Q2e,2) - 6*l1k1*M2*Q2k*pow(Q2e,2) -
		          2*m2*M2*Q2k*pow(Q2e,2) + 2*l1k1*Q2h*Q2k*pow(Q2e,2) -
		          4*m2*Q2h*Q2k*pow(Q2e,2) + 3*M2*Q2h*Q2k*pow(Q2e,2) -
		          4*l1k1*Q2h*S*pow(Q2e,2) + 8*l1k1*Q2k*S*pow(Q2e,2) +
		          8*m2*Q2k*S*pow(Q2e,2) - 4*Q2h*Q2k*S*pow(Q2e,2) +
		          4*l1k1*Q2h*Sk*pow(Q2e,2) - 8*m2*Q2h*Sk*pow(Q2e,2) +
		          12*m2*Q2k*Sk*pow(Q2e,2) + 2*Q2h*Q2k*Sk*pow(Q2e,2) +
		          16*l1k1*S*Sk*pow(Q2e,2) + 16*m2*S*Sk*pow(Q2e,2) -
		          16*Q2h*S*Sk*pow(Q2e,2) + 4*l1k1*Q2h*Sq2*pow(Q2e,2) -
		          4*Q2h*Q2k*Sq2*pow(Q2e,2) + 28*Q2h*S*Sq2*pow(Q2e,2) -
		          4*Q2k*S*Sq2*pow(Q2e,2) - 12*Q2h*Sk*Sq2*pow(Q2e,2) +
		          16*M2*pow(l1k1,2)*pow(Q2e,2) + 2*l1k1*M2*pow(Q2e,3) -
		          8*l1k1*m4*pow(Q2h,2) - 4*l1k1*m2*M2*pow(Q2h,2) +
		          8*l1k1*m2*Q2e*pow(Q2h,2) + 10*l1k1*M2*Q2e*pow(Q2h,2) +
		          12*m2*M2*Q2e*pow(Q2h,2) - 6*l1k1*m2*Q2k*pow(Q2h,2) +
		          8*m4*Q2k*pow(Q2h,2) + 26*l1k1*M2*Q2k*pow(Q2h,2) -
		          6*m2*M2*Q2k*pow(Q2h,2) - 4*l1k1*Q2e*Q2k*pow(Q2h,2) +
		          m2*Q2e*Q2k*pow(Q2h,2) - 19*M2*Q2e*Q2k*pow(Q2h,2) +
		          8*l1k1*m2*S*pow(Q2h,2) - 16*m4*S*pow(Q2h,2) +
		          12*l1k1*Q2e*S*pow(Q2h,2) - 12*l1k1*Q2k*S*pow(Q2h,2) -
		          8*m2*Q2k*S*pow(Q2h,2) + 12*Q2e*Q2k*S*pow(Q2h,2) -
		          4*l1k1*m2*Sk*pow(Q2h,2) + 8*m4*Sk*pow(Q2h,2) -
		          4*l1k1*Q2e*Sk*pow(Q2h,2) + 8*m2*Q2e*Sk*pow(Q2h,2) -
		          4*m2*Q2k*Sk*pow(Q2h,2) - 6*Q2e*Q2k*Sk*pow(Q2h,2) -
		          40*l1k1*S*Sk*pow(Q2h,2) + 48*Q2e*S*Sk*pow(Q2h,2) +
		          16*Q2k*S*Sk*pow(Q2h,2) + 24*l1k1*m2*Sq2*pow(Q2h,2) +
		          12*l1k1*Q2e*Sq2*pow(Q2h,2) - 6*m2*Q2e*Sq2*pow(Q2h,2) -
		          22*m2*Q2k*Sq2*pow(Q2h,2) + 6*Q2e*Q2k*Sq2*pow(Q2h,2) +
		          72*l1k1*S*Sq2*pow(Q2h,2) + 32*m2*S*Sq2*pow(Q2h,2) -
		          84*Q2e*S*Sq2*pow(Q2h,2) - 20*Q2k*S*Sq2*pow(Q2h,2) -
		          24*l1k1*Sk*Sq2*pow(Q2h,2) + 36*Q2e*Sk*Sq2*pow(Q2h,2) +
		          12*Q2k*Sk*Sq2*pow(Q2h,2) + 8*m2*pow(l1k1,2)*pow(Q2h,2) -
		          28*M2*pow(l1k1,2)*pow(Q2h,2) + 4*Q2k*pow(l1k1,2)*pow(Q2h,2) +
		          8*Sk*pow(l1k1,2)*pow(Q2h,2) - 8*Sq2*pow(l1k1,2)*pow(Q2h,2) -
		          2*l1k1*pow(Q2e,2)*pow(Q2h,2) + m2*pow(Q2e,2)*pow(Q2h,2) -
		          2*M2*pow(Q2e,2)*pow(Q2h,2) + Q2k*pow(Q2e,2)*pow(Q2h,2) +
		          2*Sq2*pow(Q2e,2)*pow(Q2h,2) - 4*l1k1*m2*pow(Q2h,3) -
		          4*m4*pow(Q2h,3) - 22*l1k1*M2*pow(Q2h,3) -
		          6*m2*M2*pow(Q2h,3) + 4*l1k1*Q2e*pow(Q2h,3) +
		          3*m2*Q2e*pow(Q2h,3) + 14*M2*Q2e*pow(Q2h,3) +
		          6*l1k1*Q2k*pow(Q2h,3) + 3*m2*Q2k*pow(Q2h,3) +
		          16*M2*Q2k*pow(Q2h,3) - 4*Q2e*Q2k*pow(Q2h,3) -
		          4*l1k1*S*pow(Q2h,3) - 8*Q2k*S*pow(Q2h,3) + 4*l1k1*Sk*pow(Q2h,3) +
		          4*Q2k*Sk*pow(Q2h,3) - 32*S*Sk*pow(Q2h,3) -
		          12*l1k1*Sq2*pow(Q2h,3) + 10*m2*Sq2*pow(Q2h,3) +
		          2*Q2e*Sq2*pow(Q2h,3) - 2*Q2k*Sq2*pow(Q2h,3) +
		          56*S*Sq2*pow(Q2h,3) - 24*Sk*Sq2*pow(Q2h,3) -
		          4*pow(l1k1,2)*pow(Q2h,3) - pow(Q2e,2)*pow(Q2h,3) -
		          4*l1k1*pow(Q2h,4) - 4*m2*pow(Q2h,4) - 12*M2*pow(Q2h,4) +
		          3*Q2e*pow(Q2h,4) + 3*Q2k*pow(Q2h,4) - 4*Sq2*pow(Q2h,4) -
		          2*pow(Q2h,5) + 16*l1k1*m2*M2*pow(Q2k,2) +
		          6*l1k1*m2*Q2e*pow(Q2k,2) + 8*l1k1*M2*Q2e*pow(Q2k,2) +
		          6*l1k1*m2*Q2h*pow(Q2k,2) - 8*m4*Q2h*pow(Q2k,2) -
		          8*l1k1*M2*Q2h*pow(Q2k,2) + 4*m2*M2*Q2h*pow(Q2k,2) -
		          2*m2*Q2e*Q2h*pow(Q2k,2) + 8*M2*Q2e*Q2h*pow(Q2k,2) +
		          16*l1k1*m2*S*pow(Q2k,2) - 32*m4*S*pow(Q2k,2) -
		          12*l1k1*Q2e*S*pow(Q2k,2) - 16*m2*Q2e*S*pow(Q2k,2) +
		          8*l1k1*Q2h*S*pow(Q2k,2) + 8*m2*Q2h*S*pow(Q2k,2) -
		          4*Q2e*Q2h*S*pow(Q2k,2) + 2*Q2e*Q2h*Sk*pow(Q2k,2) -
		          12*m2*Q2e*Sq2*pow(Q2k,2) + 4*l1k1*Q2h*Sq2*pow(Q2k,2) +
		          12*m2*Q2h*Sq2*pow(Q2k,2) - 2*Q2e*Q2h*Sq2*pow(Q2k,2) +
		          16*l1k1*S*Sq2*pow(Q2k,2) + 16*m2*S*Sq2*pow(Q2k,2) +
		          4*Q2e*S*Sq2*pow(Q2k,2) - 4*Q2h*S*Sq2*pow(Q2k,2) +
		          3*m2*pow(Q2e,2)*pow(Q2k,2) - M2*pow(Q2e,2)*pow(Q2k,2) -
		          2*l1k1*pow(Q2h,2)*pow(Q2k,2) - m2*pow(Q2h,2)*pow(Q2k,2) -
		          7*M2*pow(Q2h,2)*pow(Q2k,2) + Q2e*pow(Q2h,2)*pow(Q2k,2) +
		          4*S*pow(Q2h,2)*pow(Q2k,2) - 2*Sk*pow(Q2h,2)*pow(Q2k,2) +
		          2*Sq2*pow(Q2h,2)*pow(Q2k,2) - pow(Q2h,3)*pow(Q2k,2) -
		          M2*Q2e*pow(Q2k,3) + M2*Q2h*pow(Q2k,3) -
		          64*l1k1*m4*pow(S,2) + 64*l1k1*m2*Q2e*pow(S,2) -
		          32*m4*Q2e*pow(S,2) + 64*l1k1*m2*Q2h*pow(S,2) -
		          32*l1k1*Q2e*Q2h*pow(S,2) - 32*m2*Q2e*Q2h*pow(S,2) -
		          128*l1k1*m2*Q2k*pow(S,2) + 32*m4*Q2k*pow(S,2) -
		          16*l1k1*Q2e*Q2k*pow(S,2) - 16*m2*Q2e*Q2k*pow(S,2) -
		          32*l1k1*Q2h*Q2k*pow(S,2) - 16*m2*Q2h*Q2k*pow(S,2) +
		          16*Q2e*Q2h*Q2k*pow(S,2) + 192*m2*pow(l1k1,2)*pow(S,2) +
		          32*Q2e*pow(l1k1,2)*pow(S,2) + 96*Q2h*pow(l1k1,2)*pow(S,2) -
		          96*Q2k*pow(l1k1,2)*pow(S,2) + 128*pow(l1k1,3)*pow(S,2) +
		          16*m2*pow(Q2e,2)*pow(S,2) + 20*Q2h*pow(Q2e,2)*pow(S,2) -
		          4*Q2k*pow(Q2e,2)*pow(S,2) + 64*l1k1*pow(Q2h,2)*pow(S,2) +
		          32*m2*pow(Q2h,2)*pow(S,2) - 60*Q2e*pow(Q2h,2)*pow(S,2) -
		          12*Q2k*pow(Q2h,2)*pow(S,2) + 40*pow(Q2h,3)*pow(S,2) +
		          16*l1k1*pow(Q2k,2)*pow(S,2) + 16*m2*pow(Q2k,2)*pow(S,2) +
		          4*Q2e*pow(Q2k,2)*pow(S,2) - 4*Q2h*pow(Q2k,2)*pow(S,2) +
		          24*l1k1*m2*Q2e*pow(Sk,2) - 24*l1k1*m2*Q2h*pow(Sk,2) -
		          8*l1k1*Q2e*Q2h*pow(Sk,2) - 8*m2*Q2e*Q2h*pow(Sk,2) +
		          4*Q2e*Q2h*Q2k*pow(Sk,2) + 12*m2*pow(Q2e,2)*pow(Sk,2) +
		          4*Q2h*pow(Q2e,2)*pow(Sk,2) + 8*l1k1*pow(Q2h,2)*pow(Sk,2) -
		          4*m2*pow(Q2h,2)*pow(Sk,2) - 12*Q2e*pow(Q2h,2)*pow(Sk,2) -
		          4*Q2k*pow(Q2h,2)*pow(Sk,2) + 8*pow(Q2h,3)*pow(Sk,2) -
		          16*l1k1*m2*Q2h*pow(Sq2,2) - 24*l1k1*Q2e*Q2h*pow(Sq2,2) -
		          8*m2*Q2e*Q2h*pow(Sq2,2) + 8*m2*Q2h*Q2k*pow(Sq2,2) +
		          12*Q2e*Q2h*Q2k*pow(Sq2,2) + 8*Q2h*pow(Q2e,2)*pow(Sq2,2) +
		          24*l1k1*pow(Q2h,2)*pow(Sq2,2) - 32*Q2e*pow(Q2h,2)*pow(Sq2,2) -
		          12*Q2k*pow(Q2h,2)*pow(Sq2,2) + 24*pow(Q2h,3)*pow(Sq2,2)))/2. +
		     (pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		        (32*k1k2*l1k2*m2*M2 + 64*k1k2*m4*M2 - 64*l1k2*m4*M2 +
		          64*k1k2*l1k2*M2*Q2e + 32*k1k2*m2*M2*Q2e +
		          32*l1k2*m2*M2*Q2e - 32*m4*M2*Q2e + 32*k1k2*m4*Q2h -
		          32*l1k2*m4*Q2h - 64*k1k2*l1k2*M2*Q2h + 96*k1k2*m2*M2*Q2h -
		          80*l1k2*m2*M2*Q2h + 32*m4*M2*Q2h - 16*m4*Q2e*Q2h -
		          24*k1k2*M2*Q2e*Q2h + 32*l1k2*M2*Q2e*Q2h -
		          40*m2*M2*Q2e*Q2h - 16*l1k2*m2*M2*Q2k - 32*m4*M2*Q2k +
		          24*k1k2*M2*Q2e*Q2k - 16*l1k2*M2*Q2e*Q2k -
		          8*m2*M2*Q2e*Q2k - 32*k1k2*m2*Q2h*Q2k - 16*l1k2*m2*Q2h*Q2k +
		          16*m4*Q2h*Q2k - 24*k1k2*M2*Q2h*Q2k + 96*l1k2*M2*Q2h*Q2k -
		          24*m2*M2*Q2h*Q2k + 8*k1k2*Q2e*Q2h*Q2k - 8*l1k2*Q2e*Q2h*Q2k +
		          16*M2*Q2e*Q2h*Q2k - 96*k1k2*l1k2*Q2h*S + 64*k1k2*m2*Q2h*S -
		          64*l1k2*m2*Q2h*S + 16*k1k2*Q2e*Q2h*S - 16*l1k2*Q2e*Q2h*S +
		          128*k1k2*m2*Q2k*S - 64*l1k2*m2*Q2k*S - 128*m4*Q2k*S -
		          48*k1k2*Q2e*Q2k*S - 64*m2*Q2e*Q2k*S - 16*k1k2*Q2h*Q2k*S +
		          96*l1k2*Q2h*Q2k*S - 24*Q2e*Q2h*Q2k*S - 64*k1k2*m4*Sk +
		          64*l1k2*m4*Sk + 32*m4*Q2e*Sk - 64*k1k2*m2*Q2h*Sk -
		          32*l1k2*m2*Q2h*Sk + 16*k1k2*Q2e*Q2h*Sk - 16*l1k2*Q2e*Q2h*Sk -
		          32*k1k2*m2*Q2k*Sk + 32*l1k2*m2*Q2k*Sk + 32*m4*Q2k*Sk +
		          16*m2*Q2e*Q2k*Sk - 48*m2*Q2h*Q2k*Sk + 256*k1k2*m2*S*Sk -
		          128*l1k2*m2*S*Sk - 256*m4*S*Sk - 96*k1k2*Q2e*S*Sk -
		          128*m2*Q2e*S*Sk + 32*k1k2*Q2h*S*Sk + 64*m2*Q2h*S*Sk -
		          64*Q2e*Q2h*S*Sk - 64*m2*Q2k*S*Sk - 64*k1k2*l1k2*Q2h*Sq2 +
		          96*k1k2*m2*Q2h*Sq2 + 32*l1k2*m2*Q2h*Sq2 - 32*m4*Q2h*Sq2 +
		          16*k1k2*Q2e*Q2h*Sq2 + 48*m2*Q2e*Q2h*Sq2 + 128*k1k2*m2*Q2k*Sq2 -
		          64*l1k2*m2*Q2k*Sq2 - 128*m4*Q2k*Sq2 - 48*k1k2*Q2e*Q2k*Sq2 -
		          64*m2*Q2e*Q2k*Sq2 - 16*k1k2*Q2h*Q2k*Sq2 + 64*l1k2*Q2h*Q2k*Sq2 +
		          48*m2*Q2h*Q2k*Sq2 - 16*Q2e*Q2h*Q2k*Sq2 + 640*k1k2*l1k2*S*Sq2 -
		          768*k1k2*m2*S*Sq2 + 512*l1k2*m2*S*Sq2 + 384*m4*S*Sq2 +
		          192*k1k2*Q2e*S*Sq2 - 320*l1k2*Q2e*S*Sq2 + 256*m2*Q2e*S*Sq2 -
		          512*k1k2*Q2h*S*Sq2 + 1024*l1k2*Q2h*S*Sq2 - 512*m2*Q2h*S*Sq2 +
		          368*Q2e*Q2h*S*Sq2 + 256*k1k2*Q2k*S*Sq2 - 512*l1k2*Q2k*S*Sq2 +
		          320*m2*Q2k*S*Sq2 - 144*Q2e*Q2k*S*Sq2 + 368*Q2h*Q2k*S*Sq2 +
		          256*k1k2*m2*Sk*Sq2 - 128*l1k2*m2*Sk*Sq2 - 256*m4*Sk*Sq2 -
		          96*k1k2*Q2e*Sk*Sq2 - 128*m2*Q2e*Sk*Sq2 + 32*k1k2*Q2h*Sk*Sq2 +
		          128*m2*Q2h*Sk*Sq2 - 48*Q2e*Q2h*Sk*Sq2 - 64*m2*Q2k*Sk*Sq2 +
		          32*m2*M2*pow(k1k2,2) - 16*M2*Q2e*pow(k1k2,2) +
		          16*M2*Q2h*pow(k1k2,2) + 32*Q2h*S*pow(k1k2,2) +
		          32*Q2h*Sq2*pow(k1k2,2) - 256*S*Sq2*pow(k1k2,2) +
		          32*m2*M2*pow(l1k2,2) + 128*M2*Q2h*pow(l1k2,2) +
		          128*Q2h*S*pow(l1k2,2) + 64*Q2h*Sq2*pow(l1k2,2) -
		          768*S*Sq2*pow(l1k2,2) + 32*l1k2*M2*pow(Q2e,2) +
		          8*m2*M2*pow(Q2e,2) - 16*M2*Q2h*pow(Q2e,2) +
		          20*M2*Q2k*pow(Q2e,2) + 4*Q2h*Q2k*pow(Q2e,2) -
		          8*Q2k*S*pow(Q2e,2) + 8*Q2h*Sk*pow(Q2e,2) - 16*S*Sk*pow(Q2e,2) -
		          8*Q2k*Sq2*pow(Q2e,2) - 32*S*Sq2*pow(Q2e,2) -
		          16*Sk*Sq2*pow(Q2e,2) + 4*M2*pow(Q2e,3) +
		          32*k1k2*m2*pow(Q2h,2) - 16*l1k2*m2*pow(Q2h,2) +
		          16*m4*pow(Q2h,2) + 40*k1k2*M2*pow(Q2h,2) -
		          176*l1k2*M2*pow(Q2h,2) + 80*m2*M2*pow(Q2h,2) +
		          16*l1k2*Q2e*pow(Q2h,2) - 16*m2*Q2e*pow(Q2h,2) -
		          48*M2*Q2e*pow(Q2h,2) + 24*l1k2*Q2k*pow(Q2h,2) -
		          68*M2*Q2k*pow(Q2h,2) + 8*Q2e*Q2k*pow(Q2h,2) +
		          48*k1k2*S*pow(Q2h,2) - 112*l1k2*S*pow(Q2h,2) +
		          64*m2*S*pow(Q2h,2) + 8*Q2e*S*pow(Q2h,2) - 16*Q2k*S*pow(Q2h,2) +
		          16*l1k2*Sk*pow(Q2h,2) + 16*m2*Sk*pow(Q2h,2) +
		          80*S*Sk*pow(Q2h,2) + 32*k1k2*Sq2*pow(Q2h,2) -
		          80*l1k2*Sq2*pow(Q2h,2) - 16*m2*Sq2*pow(Q2h,2) -
		          8*Q2e*Sq2*pow(Q2h,2) - 16*Q2k*Sq2*pow(Q2h,2) -
		          496*S*Sq2*pow(Q2h,2) + 64*Sk*Sq2*pow(Q2h,2) +
		          32*pow(l1k2,2)*pow(Q2h,2) - 40*l1k2*pow(Q2h,3) +
		          24*m2*pow(Q2h,3) + 84*M2*pow(Q2h,3) - 12*Q2e*pow(Q2h,3) -
		          16*Q2k*pow(Q2h,3) + 24*S*pow(Q2h,3) - 8*Sk*pow(Q2h,3) +
		          32*Sq2*pow(Q2h,3) + 16*pow(Q2h,4) - 16*m2*M2*pow(Q2k,2) -
		          12*M2*Q2e*pow(Q2k,2) - 8*m2*Q2h*pow(Q2k,2) +
		          20*M2*Q2h*pow(Q2k,2) - 4*Q2e*Q2h*pow(Q2k,2) -
		          32*m2*S*pow(Q2k,2) + 16*Q2h*S*pow(Q2k,2) +
		          16*m2*Sk*pow(Q2k,2) - 32*m2*Sq2*pow(Q2k,2) +
		          16*Q2h*Sq2*pow(Q2k,2) - 64*S*Sq2*pow(Q2k,2) +
		          4*pow(Q2h,2)*pow(Q2k,2) + 384*k1k2*l1k2*pow(S,2) -
		          384*k1k2*m2*pow(S,2) + 384*l1k2*m2*pow(S,2) +
		          128*m4*pow(S,2) + 128*k1k2*Q2e*pow(S,2) -
		          256*l1k2*Q2e*pow(S,2) + 192*m2*Q2e*pow(S,2) -
		          320*k1k2*Q2h*pow(S,2) + 704*l1k2*Q2h*pow(S,2) -
		          384*m2*Q2h*pow(S,2) + 288*Q2e*Q2h*pow(S,2) +
		          128*k1k2*Q2k*pow(S,2) - 320*l1k2*Q2k*pow(S,2) +
		          192*m2*Q2k*pow(S,2) - 128*Q2e*Q2k*pow(S,2) +
		          256*Q2h*Q2k*pow(S,2) - 128*pow(k1k2,2)*pow(S,2) -
		          512*pow(l1k2,2)*pow(S,2) - 32*pow(Q2e,2)*pow(S,2) -
		          352*pow(Q2h,2)*pow(S,2) - 32*pow(Q2k,2)*pow(S,2) -
		          64*k1k2*m2*pow(Sk,2) + 64*l1k2*m2*pow(Sk,2) +
		          32*m2*Q2e*pow(Sk,2) - 64*m2*Q2h*pow(Sk,2) +
		          16*Q2e*Q2h*pow(Sk,2) + 32*m2*Q2k*pow(Sk,2) -
		          16*pow(Q2h,2)*pow(Sk,2) + 256*k1k2*l1k2*pow(Sq2,2) -
		          384*k1k2*m2*pow(Sq2,2) + 128*l1k2*m2*pow(Sq2,2) +
		          256*m4*pow(Sq2,2) + 64*k1k2*Q2e*pow(Sq2,2) -
		          64*l1k2*Q2e*pow(Sq2,2) + 64*m2*Q2e*pow(Sq2,2) -
		          192*k1k2*Q2h*pow(Sq2,2) + 320*l1k2*Q2h*pow(Sq2,2) -
		          160*m2*Q2h*pow(Sq2,2) + 96*Q2e*Q2h*pow(Sq2,2) +
		          128*k1k2*Q2k*pow(Sq2,2) - 192*l1k2*Q2k*pow(Sq2,2) +
		          128*m2*Q2k*pow(Sq2,2) - 16*Q2e*Q2k*pow(Sq2,2) +
		          112*Q2h*Q2k*pow(Sq2,2) - 128*pow(k1k2,2)*pow(Sq2,2) -
		          256*pow(l1k2,2)*pow(Sq2,2) - 160*pow(Q2h,2)*pow(Sq2,2) -
		          32*pow(Q2k,2)*pow(Sq2,2)))/2. +
		     (pow(l1k2,-1)*(32*k1k2*l1k1*m2*M2 - 32*k1k2*m4*M2 +
		          24*k1k2*l1k1*M2*Q2e - 8*k1k2*m2*M2*Q2e +
		          8*l1k1*m2*M2*Q2e - 16*k1k2*m4*Q2h - 8*k1k2*l1k1*M2*Q2h -
		          40*k1k2*m2*M2*Q2h + 8*l1k1*m2*M2*Q2h - 16*m4*M2*Q2h -
		          12*k1k2*M2*Q2e*Q2h + 8*l1k1*M2*Q2e*Q2h - 4*m2*M2*Q2e*Q2h -
		          16*l1k1*m2*M2*Q2k + 16*m4*M2*Q2k - 8*l1k1*M2*Q2e*Q2k +
		          4*m2*M2*Q2e*Q2k + 8*k1k2*m2*Q2h*Q2k - 24*l1k1*m2*Q2h*Q2k +
		          24*m4*Q2h*Q2k + 16*l1k1*M2*Q2h*Q2k + 4*m2*M2*Q2h*Q2k +
		          10*m2*Q2e*Q2h*Q2k - 2*M2*Q2e*Q2h*Q2k + 32*k1k2*m2*Q2h*S -
		          32*m4*Q2h*S - 16*l1k1*Q2e*Q2h*S + 8*m2*Q2e*Q2h*S -
		          64*k1k2*m2*Q2k*S - 32*l1k1*m2*Q2k*S + 64*m4*Q2k*S +
		          24*k1k2*Q2e*Q2k*S + 24*l1k1*Q2e*Q2k*S + 16*m2*Q2e*Q2k*S -
		          8*k1k2*Q2h*Q2k*S - 16*l1k1*Q2h*Q2k*S - 32*m2*Q2h*Q2k*S +
		          20*Q2e*Q2h*Q2k*S + 32*k1k2*m4*Sk + 16*k1k2*m2*Q2h*Sk -
		          48*l1k1*m2*Q2h*Sk + 48*m4*Q2h*Sk + 20*m2*Q2e*Q2h*Sk +
		          16*k1k2*m2*Q2k*Sk - 80*m4*Q2k*Sk - 40*m2*Q2e*Q2k*Sk +
		          16*m2*Q2h*Q2k*Sk - 8*Q2e*Q2h*Q2k*Sk - 128*k1k2*m2*S*Sk -
		          64*l1k1*m2*S*Sk + 128*m4*S*Sk + 48*k1k2*Q2e*S*Sk +
		          48*l1k1*Q2e*S*Sk + 32*m2*Q2e*S*Sk - 16*k1k2*Q2h*S*Sk -
		          16*l1k1*Q2h*S*Sk - 32*m2*Q2h*S*Sk + 40*Q2e*Q2h*S*Sk +
		          64*m2*Q2k*S*Sk - 24*Q2e*Q2k*S*Sk + 8*Q2h*Q2k*S*Sk +
		          16*k1k2*m2*Q2h*Sq2 + 32*l1k1*m2*Q2h*Sq2 - 16*m4*Q2h*Sq2 -
		          8*k1k2*Q2e*Q2h*Sq2 + 8*l1k1*Q2e*Q2h*Sq2 - 64*k1k2*m2*Q2k*Sq2 +
		          64*m4*Q2k*Sq2 + 24*k1k2*Q2e*Q2k*Sq2 + 32*m2*Q2e*Q2k*Sq2 -
		          8*k1k2*Q2h*Q2k*Sq2 - 8*l1k1*Q2h*Q2k*Sq2 - 48*m2*Q2h*Q2k*Sq2 +
		          16*Q2e*Q2h*Q2k*Sq2 + 64*k1k2*l1k1*S*Sq2 + 256*k1k2*m2*S*Sq2 +
		          128*l1k1*m2*S*Sq2 - 64*m4*S*Sq2 + 96*k1k2*Q2h*S*Sq2 +
		          64*l1k1*Q2h*S*Sq2 + 192*m2*Q2h*S*Sq2 - 24*Q2e*Q2h*S*Sq2 -
		          64*k1k2*Q2k*S*Sq2 - 64*l1k1*Q2k*S*Sq2 - 160*m2*Q2k*S*Sq2 -
		          8*Q2e*Q2k*S*Sq2 - 40*Q2h*Q2k*S*Sq2 - 128*k1k2*m2*Sk*Sq2 +
		          128*m4*Sk*Sq2 + 48*k1k2*Q2e*Sk*Sq2 + 64*m2*Q2e*Sk*Sq2 -
		          16*k1k2*Q2h*Sk*Sq2 - 64*m2*Q2h*Sk*Sq2 + 32*Q2e*Q2h*Sk*Sq2 +
		          64*m2*Q2k*Sk*Sq2 - 24*Q2e*Q2k*Sk*Sq2 + 8*Q2h*Q2k*Sk*Sq2 -
		          16*m2*M2*pow(k1k2,2) - 8*M2*Q2e*pow(k1k2,2) +
		          64*S*Sq2*pow(k1k2,2) - 16*m2*M2*pow(l1k1,2) -
		          8*M2*Q2e*pow(l1k1,2) - 16*M2*Q2h*pow(l1k1,2) +
		          16*Q2h*S*pow(l1k1,2) + 16*Q2h*Sq2*pow(l1k1,2) +
		          64*S*Sq2*pow(l1k1,2) + 8*k1k2*M2*pow(Q2e,2) +
		          2*M2*Q2h*pow(Q2e,2) - 2*M2*Q2k*pow(Q2e,2) -
		          4*Q2h*Q2k*pow(Q2e,2) - 12*Q2h*S*pow(Q2e,2) + 8*Q2k*S*pow(Q2e,2) -
		          8*Q2h*Sk*pow(Q2e,2) + 16*Q2k*Sk*pow(Q2e,2) + 16*S*Sk*pow(Q2e,2) -
		          4*Q2h*Sq2*pow(Q2e,2) - 2*M2*pow(Q2e,3) -
		          16*k1k2*m2*pow(Q2h,2) + 8*l1k1*m2*pow(Q2h,2) -
		          16*m4*pow(Q2h,2) - 4*k1k2*M2*pow(Q2h,2) -
		          24*l1k1*M2*pow(Q2h,2) - 24*m2*M2*pow(Q2h,2) -
		          4*m2*Q2e*pow(Q2h,2) + 4*l1k1*Q2k*pow(Q2h,2) +
		          2*m2*Q2k*pow(Q2h,2) + 12*M2*Q2k*pow(Q2h,2) +
		          2*Q2e*Q2k*pow(Q2h,2) - 8*l1k1*S*pow(Q2h,2) - 8*m2*S*pow(Q2h,2) +
		          4*Q2e*S*pow(Q2h,2) - 12*Q2k*S*pow(Q2h,2) + 8*l1k1*Sk*pow(Q2h,2) -
		          4*m2*Sk*pow(Q2h,2) + 4*Q2e*Sk*pow(Q2h,2) + 4*Q2k*Sk*pow(Q2h,2) -
		          40*S*Sk*pow(Q2h,2) - 16*l1k1*Sq2*pow(Q2h,2) +
		          8*m2*Sq2*pow(Q2h,2) - 8*Q2k*Sq2*pow(Q2h,2) +
		          88*S*Sq2*pow(Q2h,2) - 32*Sk*Sq2*pow(Q2h,2) - 4*l1k1*pow(Q2h,3) -
		          8*m2*pow(Q2h,3) - 14*M2*pow(Q2h,3) + 2*Q2k*pow(Q2h,3) -
		          4*Sq2*pow(Q2h,3) - 2*pow(Q2h,4) - 16*m4*pow(Q2k,2) +
		          4*m2*M2*pow(Q2k,2) - 10*m2*Q2e*pow(Q2k,2) +
		          2*m2*Q2h*pow(Q2k,2) - 2*M2*Q2h*pow(Q2k,2) -
		          2*Q2e*Q2h*pow(Q2k,2) + 32*m2*S*pow(Q2k,2) -
		          12*Q2e*S*pow(Q2k,2) + 4*Q2h*S*pow(Q2k,2) - 8*m2*Sk*pow(Q2k,2) +
		          32*m2*Sq2*pow(Q2k,2) - 12*Q2e*Sq2*pow(Q2k,2) +
		          4*Q2h*Sq2*pow(Q2k,2) + 16*S*Sq2*pow(Q2k,2) +
		          4*pow(Q2e,2)*pow(Q2k,2) + 64*k1k2*l1k1*pow(S,2) +
		          128*k1k2*m2*pow(S,2) + 128*l1k1*m2*pow(S,2) +
		          16*l1k1*Q2e*pow(S,2) + 32*m2*Q2e*pow(S,2) +
		          64*k1k2*Q2h*pow(S,2) + 80*l1k1*Q2h*pow(S,2) +
		          128*m2*Q2h*pow(S,2) - 8*Q2e*Q2h*pow(S,2) -
		          32*k1k2*Q2k*pow(S,2) - 64*l1k1*Q2k*pow(S,2) -
		          96*m2*Q2k*pow(S,2) - 8*Q2e*Q2k*pow(S,2) - 24*Q2h*Q2k*pow(S,2) +
		          32*pow(k1k2,2)*pow(S,2) + 96*pow(l1k1,2)*pow(S,2) +
		          64*pow(Q2h,2)*pow(S,2) + 8*pow(Q2k,2)*pow(S,2) +
		          32*k1k2*m2*pow(Sk,2) - 64*m4*pow(Sk,2) - 40*m2*Q2e*pow(Sk,2) +
		          24*m2*Q2h*pow(Sk,2) - 8*Q2e*Q2h*pow(Sk,2) -
		          16*m2*Q2k*pow(Sk,2) + 16*pow(Q2e,2)*pow(Sk,2) +
		          8*pow(Q2h,2)*pow(Sk,2) + 128*k1k2*m2*pow(Sq2,2) -
		          64*m4*pow(Sq2,2) - 32*m2*Q2e*pow(Sq2,2) +
		          32*k1k2*Q2h*pow(Sq2,2) + 80*m2*Q2h*pow(Sq2,2) -
		          8*Q2e*Q2h*pow(Sq2,2) - 32*k1k2*Q2k*pow(Sq2,2) -
		          64*m2*Q2k*pow(Sq2,2) - 16*Q2h*Q2k*pow(Sq2,2) +
		          32*pow(k1k2,2)*pow(Sq2,2) + 32*pow(Q2h,2)*pow(Sq2,2) +
		          8*pow(Q2k,2)*pow(Sq2,2)))/2. +
		     (pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		        (-64*k1k2*l1k2*m2*M2*Q2h + 64*l1k2*m4*M2*Q2h +
		          32*k1k2*l1k2*M2*Q2e*Q2h - 32*k1k2*m2*M2*Q2e*Q2h +
		          32*m4*M2*Q2e*Q2h + 32*k1k2*l1k2*M2*Q2h*Q2k +
		          16*k1k2*M2*Q2e*Q2h*Q2k - 32*l1k2*M2*Q2e*Q2h*Q2k -
		          128*k1k2*l1k2*m2*Q2h*S + 128*l1k2*m4*Q2h*S +
		          64*k1k2*l1k2*Q2e*Q2h*S - 64*k1k2*m2*Q2e*Q2h*S +
		          64*m4*Q2e*Q2h*S + 64*k1k2*l1k2*Q2h*Q2k*S +
		          32*k1k2*Q2e*Q2h*Q2k*S - 64*l1k2*Q2e*Q2h*Q2k*S -
		          128*k1k2*l1k2*m2*Q2h*Sq2 + 128*l1k2*m4*Q2h*Sq2 +
		          64*k1k2*l1k2*Q2e*Q2h*Sq2 - 64*k1k2*m2*Q2e*Q2h*Sq2 +
		          64*m4*Q2e*Q2h*Sq2 + 64*k1k2*l1k2*Q2h*Q2k*Sq2 +
		          32*k1k2*Q2e*Q2h*Q2k*Sq2 - 64*l1k2*Q2e*Q2h*Q2k*Sq2 +
		          512*k1k2*l1k2*m2*S*Sq2 - 512*l1k2*m4*S*Sq2 -
		          256*k1k2*l1k2*Q2e*S*Sq2 + 256*k1k2*m2*Q2e*S*Sq2 -
		          256*m4*Q2e*S*Sq2 + 512*k1k2*l1k2*Q2h*S*Sq2 -
		          256*k1k2*m2*Q2h*S*Sq2 + 256*m4*Q2h*S*Sq2 +
		          128*k1k2*Q2e*Q2h*S*Sq2 - 256*l1k2*Q2e*Q2h*S*Sq2 -
		          256*k1k2*l1k2*Q2k*S*Sq2 - 128*k1k2*Q2e*Q2k*S*Sq2 +
		          256*l1k2*Q2e*Q2k*S*Sq2 + 128*k1k2*Q2h*Q2k*S*Sq2 -
		          512*l1k2*Q2h*Q2k*S*Sq2 - 128*Q2e*Q2h*Q2k*S*Sq2 -
		          32*l1k2*M2*Q2h*pow(k1k2,2) - 16*M2*Q2e*Q2h*pow(k1k2,2) -
		          64*l1k2*Q2h*S*pow(k1k2,2) - 32*Q2e*Q2h*S*pow(k1k2,2) -
		          64*l1k2*Q2h*Sq2*pow(k1k2,2) - 32*Q2e*Q2h*Sq2*pow(k1k2,2) +
		          256*l1k2*S*Sq2*pow(k1k2,2) + 128*Q2e*S*Sq2*pow(k1k2,2) -
		          128*Q2h*S*Sq2*pow(k1k2,2) + 64*k1k2*M2*Q2h*pow(l1k2,2) -
		          32*M2*Q2e*Q2h*pow(l1k2,2) - 64*M2*Q2h*Q2k*pow(l1k2,2) +
		          128*k1k2*Q2h*S*pow(l1k2,2) - 64*Q2e*Q2h*S*pow(l1k2,2) -
		          128*Q2h*Q2k*S*pow(l1k2,2) + 128*k1k2*Q2h*Sq2*pow(l1k2,2) -
		          64*Q2e*Q2h*Sq2*pow(l1k2,2) - 128*Q2h*Q2k*Sq2*pow(l1k2,2) -
		          512*k1k2*S*Sq2*pow(l1k2,2) + 256*Q2e*S*Sq2*pow(l1k2,2) -
		          768*Q2h*S*Sq2*pow(l1k2,2) + 512*Q2k*S*Sq2*pow(l1k2,2) -
		          64*M2*Q2h*pow(l1k2,3) - 128*Q2h*S*pow(l1k2,3) -
		          128*Q2h*Sq2*pow(l1k2,3) + 512*S*Sq2*pow(l1k2,3) -
		          64*k1k2*l1k2*M2*pow(Q2h,2) + 32*k1k2*m2*M2*pow(Q2h,2) -
		          32*m4*M2*pow(Q2h,2) - 16*k1k2*M2*Q2e*pow(Q2h,2) +
		          32*l1k2*M2*Q2e*pow(Q2h,2) - 16*k1k2*M2*Q2k*pow(Q2h,2) +
		          64*l1k2*M2*Q2k*pow(Q2h,2) + 16*M2*Q2e*Q2k*pow(Q2h,2) -
		          128*k1k2*l1k2*S*pow(Q2h,2) + 64*k1k2*m2*S*pow(Q2h,2) -
		          64*m4*S*pow(Q2h,2) - 32*k1k2*Q2e*S*pow(Q2h,2) +
		          64*l1k2*Q2e*S*pow(Q2h,2) - 32*k1k2*Q2k*S*pow(Q2h,2) +
		          128*l1k2*Q2k*S*pow(Q2h,2) + 32*Q2e*Q2k*S*pow(Q2h,2) -
		          128*k1k2*l1k2*Sq2*pow(Q2h,2) + 64*k1k2*m2*Sq2*pow(Q2h,2) -
		          64*m4*Sq2*pow(Q2h,2) - 32*k1k2*Q2e*Sq2*pow(Q2h,2) +
		          64*l1k2*Q2e*Sq2*pow(Q2h,2) - 32*k1k2*Q2k*Sq2*pow(Q2h,2) +
		          128*l1k2*Q2k*Sq2*pow(Q2h,2) + 32*Q2e*Q2k*Sq2*pow(Q2h,2) -
		          128*k1k2*S*Sq2*pow(Q2h,2) + 384*l1k2*S*Sq2*pow(Q2h,2) +
		          64*Q2e*S*Sq2*pow(Q2h,2) + 128*Q2k*S*Sq2*pow(Q2h,2) +
		          16*M2*pow(k1k2,2)*pow(Q2h,2) + 32*S*pow(k1k2,2)*pow(Q2h,2) +
		          32*Sq2*pow(k1k2,2)*pow(Q2h,2) + 96*M2*pow(l1k2,2)*pow(Q2h,2) +
		          192*S*pow(l1k2,2)*pow(Q2h,2) + 192*Sq2*pow(l1k2,2)*pow(Q2h,2) +
		          16*k1k2*M2*pow(Q2h,3) - 48*l1k2*M2*pow(Q2h,3) -
		          8*M2*Q2e*pow(Q2h,3) - 16*M2*Q2k*pow(Q2h,3) +
		          32*k1k2*S*pow(Q2h,3) - 96*l1k2*S*pow(Q2h,3) -
		          16*Q2e*S*pow(Q2h,3) - 32*Q2k*S*pow(Q2h,3) +
		          32*k1k2*Sq2*pow(Q2h,3) - 96*l1k2*Sq2*pow(Q2h,3) -
		          16*Q2e*Sq2*pow(Q2h,3) - 32*Q2k*Sq2*pow(Q2h,3) -
		          64*S*Sq2*pow(Q2h,3) + 8*M2*pow(Q2h,4) + 16*S*pow(Q2h,4) +
		          16*Sq2*pow(Q2h,4) - 16*l1k2*M2*Q2h*pow(Q2k,2) -
		          8*M2*Q2e*Q2h*pow(Q2k,2) - 32*l1k2*Q2h*S*pow(Q2k,2) -
		          16*Q2e*Q2h*S*pow(Q2k,2) - 32*l1k2*Q2h*Sq2*pow(Q2k,2) -
		          16*Q2e*Q2h*Sq2*pow(Q2k,2) + 128*l1k2*S*Sq2*pow(Q2k,2) +
		          64*Q2e*S*Sq2*pow(Q2k,2) - 64*Q2h*S*Sq2*pow(Q2k,2) +
		          8*M2*pow(Q2h,2)*pow(Q2k,2) + 16*S*pow(Q2h,2)*pow(Q2k,2) +
		          16*Sq2*pow(Q2h,2)*pow(Q2k,2) + 256*k1k2*l1k2*m2*pow(S,2) -
		          256*l1k2*m4*pow(S,2) - 128*k1k2*l1k2*Q2e*pow(S,2) +
		          128*k1k2*m2*Q2e*pow(S,2) - 128*m4*Q2e*pow(S,2) +
		          256*k1k2*l1k2*Q2h*pow(S,2) - 128*k1k2*m2*Q2h*pow(S,2) +
		          128*m4*Q2h*pow(S,2) + 64*k1k2*Q2e*Q2h*pow(S,2) -
		          128*l1k2*Q2e*Q2h*pow(S,2) - 128*k1k2*l1k2*Q2k*pow(S,2) -
		          64*k1k2*Q2e*Q2k*pow(S,2) + 128*l1k2*Q2e*Q2k*pow(S,2) +
		          64*k1k2*Q2h*Q2k*pow(S,2) - 256*l1k2*Q2h*Q2k*pow(S,2) -
		          64*Q2e*Q2h*Q2k*pow(S,2) + 128*l1k2*pow(k1k2,2)*pow(S,2) +
		          64*Q2e*pow(k1k2,2)*pow(S,2) - 64*Q2h*pow(k1k2,2)*pow(S,2) -
		          256*k1k2*pow(l1k2,2)*pow(S,2) + 128*Q2e*pow(l1k2,2)*pow(S,2) -
		          384*Q2h*pow(l1k2,2)*pow(S,2) + 256*Q2k*pow(l1k2,2)*pow(S,2) +
		          256*pow(l1k2,3)*pow(S,2) - 64*k1k2*pow(Q2h,2)*pow(S,2) +
		          192*l1k2*pow(Q2h,2)*pow(S,2) + 32*Q2e*pow(Q2h,2)*pow(S,2) +
		          64*Q2k*pow(Q2h,2)*pow(S,2) - 32*pow(Q2h,3)*pow(S,2) +
		          64*l1k2*pow(Q2k,2)*pow(S,2) + 32*Q2e*pow(Q2k,2)*pow(S,2) -
		          32*Q2h*pow(Q2k,2)*pow(S,2) + 256*k1k2*l1k2*m2*pow(Sq2,2) -
		          256*l1k2*m4*pow(Sq2,2) - 128*k1k2*l1k2*Q2e*pow(Sq2,2) +
		          128*k1k2*m2*Q2e*pow(Sq2,2) - 128*m4*Q2e*pow(Sq2,2) +
		          256*k1k2*l1k2*Q2h*pow(Sq2,2) - 128*k1k2*m2*Q2h*pow(Sq2,2) +
		          128*m4*Q2h*pow(Sq2,2) + 64*k1k2*Q2e*Q2h*pow(Sq2,2) -
		          128*l1k2*Q2e*Q2h*pow(Sq2,2) - 128*k1k2*l1k2*Q2k*pow(Sq2,2) -
		          64*k1k2*Q2e*Q2k*pow(Sq2,2) + 128*l1k2*Q2e*Q2k*pow(Sq2,2) +
		          64*k1k2*Q2h*Q2k*pow(Sq2,2) - 256*l1k2*Q2h*Q2k*pow(Sq2,2) -
		          64*Q2e*Q2h*Q2k*pow(Sq2,2) + 128*l1k2*pow(k1k2,2)*pow(Sq2,2) +
		          64*Q2e*pow(k1k2,2)*pow(Sq2,2) - 64*Q2h*pow(k1k2,2)*pow(Sq2,2) -
		          256*k1k2*pow(l1k2,2)*pow(Sq2,2) +
		          128*Q2e*pow(l1k2,2)*pow(Sq2,2) - 384*Q2h*pow(l1k2,2)*pow(Sq2,2) +
		          256*Q2k*pow(l1k2,2)*pow(Sq2,2) + 256*pow(l1k2,3)*pow(Sq2,2) -
		          64*k1k2*pow(Q2h,2)*pow(Sq2,2) + 192*l1k2*pow(Q2h,2)*pow(Sq2,2) +
		          32*Q2e*pow(Q2h,2)*pow(Sq2,2) + 64*Q2k*pow(Q2h,2)*pow(Sq2,2) -
		          32*pow(Q2h,3)*pow(Sq2,2) + 64*l1k2*pow(Q2k,2)*pow(Sq2,2) +
		          32*Q2e*pow(Q2k,2)*pow(Sq2,2) - 32*Q2h*pow(Q2k,2)*pow(Sq2,2)))/2. +
		     (pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		        (32*k1k2*m4*M2*Q2e - 32*k1k2*m4*M2*Q2h +
		          16*k1k2*m4*Q2e*Q2h + 16*k1k2*m2*M2*Q2e*Q2h +
		          16*m4*M2*Q2e*Q2h + 16*k1k2*m2*M2*Q2e*Q2k -
		          16*m4*M2*Q2e*Q2k - 16*k1k2*m2*M2*Q2h*Q2k +
		          16*m4*M2*Q2h*Q2k - 20*k1k2*m2*Q2e*Q2h*Q2k +
		          8*m2*M2*Q2e*Q2h*Q2k - 16*k1k2*m2*Q2e*Q2h*S +
		          32*k1k2*m2*Q2e*Q2k*S - 32*m4*Q2e*Q2k*S -
		          32*k1k2*m2*Q2h*Q2k*S + 32*m4*Q2h*Q2k*S +
		          16*k1k2*Q2e*Q2h*Q2k*S - 32*k1k2*m4*Q2e*Sk + 32*k1k2*m4*Q2h*Sk -
		          40*k1k2*m2*Q2e*Q2h*Sk + 32*k1k2*m2*Q2e*Q2k*Sk +
		          16*m4*Q2e*Q2k*Sk - 32*k1k2*m2*Q2h*Q2k*Sk - 16*m4*Q2h*Q2k*Sk -
		          24*m2*Q2e*Q2h*Q2k*Sk + 64*k1k2*m2*Q2e*S*Sk - 64*m4*Q2e*S*Sk -
		          64*k1k2*m2*Q2h*S*Sk + 64*m4*Q2h*S*Sk + 32*k1k2*Q2e*Q2h*S*Sk +
		          32*m2*Q2e*Q2h*S*Sk - 16*m4*Q2e*Q2h*Sq2 +
		          32*k1k2*m2*Q2e*Q2k*Sq2 - 32*m4*Q2e*Q2k*Sq2 -
		          32*k1k2*m2*Q2h*Q2k*Sq2 + 32*m4*Q2h*Q2k*Sq2 +
		          16*k1k2*Q2e*Q2h*Q2k*Sq2 + 44*m2*Q2e*Q2h*Q2k*Sq2 -
		          128*k1k2*m2*Q2e*S*Sq2 + 64*m4*Q2e*S*Sq2 +
		          128*k1k2*m2*Q2h*S*Sq2 - 64*m4*Q2h*S*Sq2 -
		          96*k1k2*Q2e*Q2h*S*Sq2 - 192*m2*Q2e*Q2h*S*Sq2 +
		          64*m2*Q2e*Q2k*S*Sq2 - 64*m2*Q2h*Q2k*S*Sq2 +
		          16*Q2e*Q2h*Q2k*S*Sq2 + 64*k1k2*m2*Q2e*Sk*Sq2 -
		          64*m4*Q2e*Sk*Sq2 - 64*k1k2*m2*Q2h*Sk*Sq2 + 64*m4*Q2h*Sk*Sq2 +
		          32*k1k2*Q2e*Q2h*Sk*Sq2 + 64*m2*Q2e*Q2h*Sk*Sq2 -
		          48*m2*Q2e*Q2k*Sk*Sq2 + 48*m2*Q2h*Q2k*Sk*Sq2 +
		          16*m2*M2*Q2e*pow(k1k2,2) - 16*m2*M2*Q2h*pow(k1k2,2) -
		          64*Q2e*S*Sq2*pow(k1k2,2) + 64*Q2h*S*Sq2*pow(k1k2,2) +
		          8*k1k2*m2*M2*pow(Q2e,2) + 4*k1k2*M2*Q2h*pow(Q2e,2) -
		          8*m2*M2*Q2h*pow(Q2e,2) - 8*m2*Q2h*Q2k*pow(Q2e,2) +
		          8*k1k2*Q2h*S*pow(Q2e,2) - 8*k1k2*Q2k*S*pow(Q2e,2) -
		          16*m2*Q2k*S*pow(Q2e,2) - 4*Q2h*Q2k*S*pow(Q2e,2) -
		          16*m2*Q2h*Sk*pow(Q2e,2) + 24*m2*Q2k*Sk*pow(Q2e,2) +
		          8*Q2h*Q2k*Sk*pow(Q2e,2) - 16*k1k2*S*Sk*pow(Q2e,2) -
		          32*m2*S*Sk*pow(Q2e,2) - 24*Q2h*S*Sk*pow(Q2e,2) -
		          8*Q2k*S*Sk*pow(Q2e,2) + 8*k1k2*Q2h*Sq2*pow(Q2e,2) -
		          8*k1k2*Q2k*Sq2*pow(Q2e,2) - 16*m2*Q2k*Sq2*pow(Q2e,2) +
		          64*m2*S*Sq2*pow(Q2e,2) + 56*Q2h*S*Sq2*pow(Q2e,2) -
		          8*Q2k*S*Sq2*pow(Q2e,2) - 16*k1k2*Sk*Sq2*pow(Q2e,2) -
		          32*m2*Sk*Sq2*pow(Q2e,2) - 16*Q2h*Sk*Sq2*pow(Q2e,2) -
		          8*Q2k*Sk*Sq2*pow(Q2e,2) - 2*M2*Q2h*pow(Q2e,3) +
		          2*M2*Q2k*pow(Q2e,3) - 16*k1k2*m4*pow(Q2h,2) -
		          24*k1k2*m2*M2*pow(Q2h,2) - 16*m4*M2*pow(Q2h,2) +
		          16*k1k2*m2*Q2e*pow(Q2h,2) + 8*m4*Q2e*pow(Q2h,2) +
		          28*m2*M2*Q2e*pow(Q2h,2) + 20*k1k2*m2*Q2k*pow(Q2h,2) -
		          8*m2*M2*Q2k*pow(Q2h,2) + 2*m2*Q2e*Q2k*pow(Q2h,2) -
		          14*M2*Q2e*Q2k*pow(Q2h,2) + 16*k1k2*m2*S*pow(Q2h,2) -
		          8*k1k2*Q2e*S*pow(Q2h,2) + 16*m2*Q2e*S*pow(Q2h,2) -
		          8*k1k2*Q2k*S*pow(Q2h,2) + 16*m2*Q2k*S*pow(Q2h,2) +
		          16*Q2e*Q2k*S*pow(Q2h,2) + 40*k1k2*m2*Sk*pow(Q2h,2) +
		          16*m2*Q2e*Sk*pow(Q2h,2) - 12*Q2e*Q2k*Sk*pow(Q2h,2) -
		          16*k1k2*S*Sk*pow(Q2h,2) + 64*Q2e*S*Sk*pow(Q2h,2) +
		          8*Q2k*S*Sk*pow(Q2h,2) + 16*m4*Sq2*pow(Q2h,2) -
		          8*k1k2*Q2e*Sq2*pow(Q2h,2) - 4*m2*Q2e*Sq2*pow(Q2h,2) -
		          8*k1k2*Q2k*Sq2*pow(Q2h,2) - 28*m2*Q2k*Sq2*pow(Q2h,2) +
		          8*Q2e*Q2k*Sq2*pow(Q2h,2) + 96*k1k2*S*Sq2*pow(Q2h,2) +
		          128*m2*S*Sq2*pow(Q2h,2) - 144*Q2e*S*Sq2*pow(Q2h,2) -
		          8*Q2k*S*Sq2*pow(Q2h,2) - 16*k1k2*Sk*Sq2*pow(Q2h,2) -
		          32*m2*Sk*Sq2*pow(Q2h,2) + 48*Q2e*Sk*Sq2*pow(Q2h,2) +
		          8*Q2k*Sk*Sq2*pow(Q2h,2) + 2*m2*pow(Q2e,2)*pow(Q2h,2) -
		          2*Q2k*pow(Q2e,2)*pow(Q2h,2) - 4*Sk*pow(Q2e,2)*pow(Q2h,2) -
		          16*k1k2*m2*pow(Q2h,3) - 8*m4*pow(Q2h,3) -
		          4*k1k2*M2*pow(Q2h,3) - 20*m2*M2*pow(Q2h,3) +
		          6*m2*Q2e*pow(Q2h,3) + 16*M2*Q2e*pow(Q2h,3) +
		          6*m2*Q2k*pow(Q2h,3) + 12*M2*Q2k*pow(Q2h,3) -
		          16*m2*S*pow(Q2h,3) - 12*Q2k*S*pow(Q2h,3) + 4*Q2e*Sk*pow(Q2h,3) +
		          4*Q2k*Sk*pow(Q2h,3) - 40*S*Sk*pow(Q2h,3) + 4*m2*Sq2*pow(Q2h,3) +
		          4*Q2e*Sq2*pow(Q2h,3) - 8*Q2k*Sq2*pow(Q2h,3) +
		          88*S*Sq2*pow(Q2h,3) - 32*Sk*Sq2*pow(Q2h,3) - 8*m2*pow(Q2h,4) -
		          14*M2*pow(Q2h,4) + 2*Q2e*pow(Q2h,4) + 2*Q2k*pow(Q2h,4) -
		          4*Sq2*pow(Q2h,4) - 2*pow(Q2h,5) + 12*k1k2*m2*Q2e*pow(Q2k,2) -
		          12*m2*M2*Q2e*pow(Q2k,2) - 12*k1k2*m2*Q2h*pow(Q2k,2) +
		          12*m2*M2*Q2h*pow(Q2k,2) - 4*m2*Q2e*Q2h*pow(Q2k,2) +
		          2*M2*Q2e*Q2h*pow(Q2k,2) + 8*m2*Q2e*Sk*pow(Q2k,2) -
		          8*m2*Q2h*Sk*pow(Q2k,2) - 24*m2*Q2e*Sq2*pow(Q2k,2) +
		          24*m2*Q2h*Sq2*pow(Q2k,2) + 16*Q2e*S*Sq2*pow(Q2k,2) -
		          16*Q2h*S*Sq2*pow(Q2k,2) + 6*m2*pow(Q2e,2)*pow(Q2k,2) +
		          2*Q2h*pow(Q2e,2)*pow(Q2k,2) - 4*S*pow(Q2e,2)*pow(Q2k,2) -
		          4*Sq2*pow(Q2e,2)*pow(Q2k,2) - 2*m2*pow(Q2h,2)*pow(Q2k,2) -
		          2*M2*pow(Q2h,2)*pow(Q2k,2) - 2*Q2e*pow(Q2h,2)*pow(Q2k,2) +
		          4*S*pow(Q2h,2)*pow(Q2k,2) + 4*Sq2*pow(Q2h,2)*pow(Q2k,2) -
		          64*k1k2*m2*Q2e*pow(S,2) + 64*k1k2*m2*Q2h*pow(S,2) -
		          64*k1k2*Q2e*Q2h*pow(S,2) - 128*m2*Q2e*Q2h*pow(S,2) +
		          32*m2*Q2e*Q2k*pow(S,2) - 32*m2*Q2h*Q2k*pow(S,2) +
		          16*Q2e*Q2h*Q2k*pow(S,2) - 32*Q2e*pow(k1k2,2)*pow(S,2) +
		          32*Q2h*pow(k1k2,2)*pow(S,2) + 32*m2*pow(Q2e,2)*pow(S,2) +
		          40*Q2h*pow(Q2e,2)*pow(S,2) - 8*Q2k*pow(Q2e,2)*pow(S,2) +
		          64*k1k2*pow(Q2h,2)*pow(S,2) + 96*m2*pow(Q2h,2)*pow(S,2) -
		          104*Q2e*pow(Q2h,2)*pow(S,2) - 8*Q2k*pow(Q2h,2)*pow(S,2) +
		          64*pow(Q2h,3)*pow(S,2) + 8*Q2e*pow(Q2k,2)*pow(S,2) -
		          8*Q2h*pow(Q2k,2)*pow(S,2) + 16*k1k2*m2*Q2e*pow(Sk,2) -
		          16*k1k2*m2*Q2h*pow(Sk,2) - 32*m2*Q2e*Q2h*pow(Sk,2) +
		          16*m2*Q2e*Q2k*pow(Sk,2) - 16*m2*Q2h*Q2k*pow(Sk,2) +
		          24*m2*pow(Q2e,2)*pow(Sk,2) + 8*Q2h*pow(Q2e,2)*pow(Sk,2) +
		          8*m2*pow(Q2h,2)*pow(Sk,2) - 16*Q2e*pow(Q2h,2)*pow(Sk,2) +
		          8*pow(Q2h,3)*pow(Sk,2) - 64*k1k2*m2*Q2e*pow(Sq2,2) +
		          64*m4*Q2e*pow(Sq2,2) + 64*k1k2*m2*Q2h*pow(Sq2,2) -
		          64*m4*Q2h*pow(Sq2,2) - 32*k1k2*Q2e*Q2h*pow(Sq2,2) -
		          80*m2*Q2e*Q2h*pow(Sq2,2) + 32*m2*Q2e*Q2k*pow(Sq2,2) -
		          32*m2*Q2h*Q2k*pow(Sq2,2) - 32*Q2e*pow(k1k2,2)*pow(Sq2,2) +
		          32*Q2h*pow(k1k2,2)*pow(Sq2,2) + 32*m2*pow(Q2e,2)*pow(Sq2,2) +
		          16*Q2h*pow(Q2e,2)*pow(Sq2,2) + 32*k1k2*pow(Q2h,2)*pow(Sq2,2) +
		          48*m2*pow(Q2h,2)*pow(Sq2,2) - 48*Q2e*pow(Q2h,2)*pow(Sq2,2) +
		          32*pow(Q2h,3)*pow(Sq2,2) + 8*Q2e*pow(Q2k,2)*pow(Sq2,2) -
		          8*Q2h*pow(Q2k,2)*pow(Sq2,2)))/2. +
		     (pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-1)*
		        pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		        (-64*m4*Q2e*Q2h*S*Sq2 - 8*m4*M2*Q2h*pow(Q2e,2) -
		          4*m2*Q2h*Q2k*Sk*pow(Q2e,2) + 16*Q2h*Q2k*S*Sk*pow(Q2e,2) +
		          8*m4*Q2h*Sq2*pow(Q2e,2) - 2*m2*Q2h*Q2k*Sq2*pow(Q2e,2) +
		          32*m4*S*Sq2*pow(Q2e,2) + 48*m2*Q2h*S*Sq2*pow(Q2e,2) -
		          16*m2*Q2k*S*Sq2*pow(Q2e,2) - 28*Q2h*Q2k*S*Sq2*pow(Q2e,2) +
		          8*m2*Q2h*Sk*Sq2*pow(Q2e,2) + 12*Q2h*Q2k*Sk*Sq2*pow(Q2e,2) -
		          3*M2*Q2h*Q2k*pow(Q2e,3) + 4*Q2h*Q2k*S*pow(Q2e,3) -
		          2*Q2h*Q2k*Sk*pow(Q2e,3) + 16*Q2h*S*Sk*pow(Q2e,3) +
		          4*Q2h*Q2k*Sq2*pow(Q2e,3) - 28*Q2h*S*Sq2*pow(Q2e,3) +
		          4*Q2k*S*Sq2*pow(Q2e,3) + 12*Q2h*Sk*Sq2*pow(Q2e,3) +
		          16*m4*M2*Q2e*pow(Q2h,2) + 8*m2*M2*Q2e*Q2k*pow(Q2h,2) +
		          8*m2*Q2e*Q2k*Sk*pow(Q2h,2) - 32*Q2e*Q2k*S*Sk*pow(Q2h,2) -
		          16*m4*Q2e*Sq2*pow(Q2h,2) + 20*m2*Q2e*Q2k*Sq2*pow(Q2h,2) +
		          32*m4*S*Sq2*pow(Q2h,2) - 80*m2*Q2e*S*Sq2*pow(Q2h,2) +
		          16*m2*Q2k*S*Sq2*pow(Q2h,2) + 44*Q2e*Q2k*S*Sq2*pow(Q2h,2) -
		          16*m2*Q2e*Sk*Sq2*pow(Q2h,2) - 24*Q2e*Q2k*Sk*Sq2*pow(Q2h,2) -
		          4*m4*pow(Q2e,2)*pow(Q2h,2) - 10*m2*M2*pow(Q2e,2)*pow(Q2h,2) +
		          m2*Q2k*pow(Q2e,2)*pow(Q2h,2) +
		          22*M2*Q2k*pow(Q2e,2)*pow(Q2h,2) -
		          16*Q2k*S*pow(Q2e,2)*pow(Q2h,2) + 8*Q2k*Sk*pow(Q2e,2)*pow(Q2h,2) -
		          64*S*Sk*pow(Q2e,2)*pow(Q2h,2) + 2*m2*Sq2*pow(Q2e,2)*pow(Q2h,2) -
		          10*Q2k*Sq2*pow(Q2e,2)*pow(Q2h,2) +
		          112*S*Sq2*pow(Q2e,2)*pow(Q2h,2) -
		          48*Sk*Sq2*pow(Q2e,2)*pow(Q2h,2) + 2*M2*pow(Q2e,3)*pow(Q2h,2) -
		          Q2k*pow(Q2e,3)*pow(Q2h,2) - 2*Sq2*pow(Q2e,3)*pow(Q2h,2) -
		          8*m4*M2*pow(Q2h,3) + 8*m4*Q2e*pow(Q2h,3) +
		          16*m2*M2*Q2e*pow(Q2h,3) - 8*m2*M2*Q2k*pow(Q2h,3) -
		          4*m2*Q2e*Q2k*pow(Q2h,3) - 35*M2*Q2e*Q2k*pow(Q2h,3) +
		          20*Q2e*Q2k*S*pow(Q2h,3) - 4*m2*Q2k*Sk*pow(Q2h,3) -
		          10*Q2e*Q2k*Sk*pow(Q2h,3) + 80*Q2e*S*Sk*pow(Q2h,3) +
		          16*Q2k*S*Sk*pow(Q2h,3) + 8*m4*Sq2*pow(Q2h,3) -
		          12*m2*Q2e*Sq2*pow(Q2h,3) - 18*m2*Q2k*Sq2*pow(Q2h,3) +
		          8*Q2e*Q2k*Sq2*pow(Q2h,3) + 32*m2*S*Sq2*pow(Q2h,3) -
		          140*Q2e*S*Sq2*pow(Q2h,3) - 20*Q2k*S*Sq2*pow(Q2h,3) +
		          8*m2*Sk*Sq2*pow(Q2h,3) + 60*Q2e*Sk*Sq2*pow(Q2h,3) +
		          12*Q2k*Sk*Sq2*pow(Q2h,3) - 3*m2*pow(Q2e,2)*pow(Q2h,3) -
		          16*M2*pow(Q2e,2)*pow(Q2h,3) + 5*Q2k*pow(Q2e,2)*pow(Q2h,3) +
		          pow(Q2e,3)*pow(Q2h,3) - 4*m4*pow(Q2h,4) - 6*m2*M2*pow(Q2h,4) +
		          7*m2*Q2e*pow(Q2h,4) + 26*M2*Q2e*pow(Q2h,4) +
		          3*m2*Q2k*pow(Q2h,4) + 16*M2*Q2k*pow(Q2h,4) -
		          7*Q2e*Q2k*pow(Q2h,4) - 8*Q2k*S*pow(Q2h,4) + 4*Q2k*Sk*pow(Q2h,4) -
		          32*S*Sk*pow(Q2h,4) + 10*m2*Sq2*pow(Q2h,4) +
		          6*Q2e*Sq2*pow(Q2h,4) - 2*Q2k*Sq2*pow(Q2h,4) +
		          56*S*Sq2*pow(Q2h,4) - 24*Sk*Sq2*pow(Q2h,4) -
		          4*pow(Q2e,2)*pow(Q2h,4) - 4*m2*pow(Q2h,5) - 12*M2*pow(Q2h,5) +
		          5*Q2e*pow(Q2h,5) + 3*Q2k*pow(Q2h,5) - 4*Sq2*pow(Q2h,5) -
		          2*pow(Q2h,6) - 8*m2*M2*Q2e*Q2h*pow(Q2k,2) -
		          8*m2*Q2e*Q2h*Sq2*pow(Q2k,2) + 16*m2*Q2e*S*Sq2*pow(Q2k,2) -
		          16*m2*Q2h*S*Sq2*pow(Q2k,2) + 8*Q2e*Q2h*S*Sq2*pow(Q2k,2) +
		          2*m2*M2*pow(Q2e,2)*pow(Q2k,2) -
		          9*M2*Q2h*pow(Q2e,2)*pow(Q2k,2) + 4*Q2h*S*pow(Q2e,2)*pow(Q2k,2) -
		          2*Q2h*Sk*pow(Q2e,2)*pow(Q2k,2) + 2*Q2h*Sq2*pow(Q2e,2)*pow(Q2k,2) -
		          4*S*Sq2*pow(Q2e,2)*pow(Q2k,2) + M2*pow(Q2e,3)*pow(Q2k,2) +
		          6*m2*M2*pow(Q2h,2)*pow(Q2k,2) + m2*Q2e*pow(Q2h,2)*pow(Q2k,2) +
		          15*M2*Q2e*pow(Q2h,2)*pow(Q2k,2) -
		          8*Q2e*S*pow(Q2h,2)*pow(Q2k,2) + 4*Q2e*Sk*pow(Q2h,2)*pow(Q2k,2) +
		          8*m2*Sq2*pow(Q2h,2)*pow(Q2k,2) -
		          4*Q2e*Sq2*pow(Q2h,2)*pow(Q2k,2) - 4*S*Sq2*pow(Q2h,2)*pow(Q2k,2) -
		          pow(Q2e,2)*pow(Q2h,2)*pow(Q2k,2) - m2*pow(Q2h,3)*pow(Q2k,2) -
		          7*M2*pow(Q2h,3)*pow(Q2k,2) + 2*Q2e*pow(Q2h,3)*pow(Q2k,2) +
		          4*S*pow(Q2h,3)*pow(Q2k,2) - 2*Sk*pow(Q2h,3)*pow(Q2k,2) +
		          2*Sq2*pow(Q2h,3)*pow(Q2k,2) - pow(Q2h,4)*pow(Q2k,2) -
		          2*M2*Q2e*Q2h*pow(Q2k,3) + M2*pow(Q2e,2)*pow(Q2k,3) +
		          M2*pow(Q2h,2)*pow(Q2k,3) - 64*m4*Q2e*Q2h*pow(S,2) +
		          32*m4*pow(Q2e,2)*pow(S,2) + 48*m2*Q2h*pow(Q2e,2)*pow(S,2) -
		          16*m2*Q2k*pow(Q2e,2)*pow(S,2) - 20*Q2h*Q2k*pow(Q2e,2)*pow(S,2) -
		          20*Q2h*pow(Q2e,3)*pow(S,2) + 4*Q2k*pow(Q2e,3)*pow(S,2) +
		          32*m4*pow(Q2h,2)*pow(S,2) - 80*m2*Q2e*pow(Q2h,2)*pow(S,2) +
		          16*m2*Q2k*pow(Q2h,2)*pow(S,2) + 28*Q2e*Q2k*pow(Q2h,2)*pow(S,2) +
		          80*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) + 32*m2*pow(Q2h,3)*pow(S,2) -
		          100*Q2e*pow(Q2h,3)*pow(S,2) - 12*Q2k*pow(Q2h,3)*pow(S,2) +
		          40*pow(Q2h,4)*pow(S,2) + 16*m2*Q2e*pow(Q2k,2)*pow(S,2) -
		          16*m2*Q2h*pow(Q2k,2)*pow(S,2) + 8*Q2e*Q2h*pow(Q2k,2)*pow(S,2) -
		          4*pow(Q2e,2)*pow(Q2k,2)*pow(S,2) -
		          4*pow(Q2h,2)*pow(Q2k,2)*pow(S,2) -
		          8*m2*Q2h*pow(Q2e,2)*pow(Sk,2) - 4*Q2h*Q2k*pow(Q2e,2)*pow(Sk,2) -
		          4*Q2h*pow(Q2e,3)*pow(Sk,2) + 16*m2*Q2e*pow(Q2h,2)*pow(Sk,2) +
		          8*Q2e*Q2k*pow(Q2h,2)*pow(Sk,2) +
		          16*pow(Q2e,2)*pow(Q2h,2)*pow(Sk,2) - 8*m2*pow(Q2h,3)*pow(Sk,2) -
		          20*Q2e*pow(Q2h,3)*pow(Sk,2) - 4*Q2k*pow(Q2h,3)*pow(Sk,2) +
		          8*pow(Q2h,4)*pow(Sk,2) - 24*m2*Q2e*Q2h*Q2k*pow(Sq2,2) +
		          8*m2*Q2h*pow(Q2e,2)*pow(Sq2,2) -
		          12*Q2h*Q2k*pow(Q2e,2)*pow(Sq2,2) - 8*Q2h*pow(Q2e,3)*pow(Sq2,2) -
		          4*m2*Q2e*pow(Q2h,2)*pow(Sq2,2) +
		          24*m2*Q2k*pow(Q2h,2)*pow(Sq2,2) +
		          24*Q2e*Q2k*pow(Q2h,2)*pow(Sq2,2) +
		          40*pow(Q2e,2)*pow(Q2h,2)*pow(Sq2,2) -
		          4*m2*pow(Q2h,3)*pow(Sq2,2) - 56*Q2e*pow(Q2h,3)*pow(Sq2,2) -
		          12*Q2k*pow(Q2h,3)*pow(Sq2,2) + 24*pow(Q2h,4)*pow(Sq2,2) +
		          12*m2*Q2e*pow(Q2k,2)*pow(Sq2,2) - 12*m2*Q2h*pow(Q2k,2)*pow(Sq2,2)))/2.);
}

long double Melem::melem2_l2k1_1(const long double l1k1,const long double l1k2,const long double k1k2,const long double Q2e,
		const long double Q2h,const long double Q2k,const long double S,const long double Sk,const long double Sq2,
		const long double f1, const long double f2)const {

	long double p1p2 = M2 + Q2h/2.;
	long double G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
	long double G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
	long double G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;

	return G1*(48*k1k2*m2 + 32*l1k1*m2 - 24*l1k2*m2 + 8*k1k2*Q2e - 8*l1k1*Q2e -
		     20*l1k2*Q2e - 4*l1k1*Q2h + 28*l1k2*Q2h + 32*m2*Q2h + 4*Q2e*Q2h -
		     24*m2*Q2k - 4*Q2e*Q2k + 2*Q2h*Q2k - 6*pow(Q2h,2) +
		     pow(k1k2 - l1k1 - l1k2,-2)*
		      (16*l1k1*l1k2*m4 - 16*l1k1*m6 - 8*l1k1*l1k2*m2*Q2h +
		        16*l1k1*m4*Q2h + 8*l1k2*m4*Q2h - 8*m6*Q2h - 8*l1k1*m4*Q2k -
		        8*l1k2*m4*Q2k + 8*m6*Q2k + 4*l1k1*m2*Q2h*Q2k +
		        4*l1k2*m2*Q2h*Q2k - 4*m4*Q2h*Q2k + 16*m4*pow(l1k1,2) -
		        4*m2*Q2h*pow(l1k1,2) - 4*m2*Q2k*pow(l1k1,2) +
		        2*Q2h*Q2k*pow(l1k1,2) + 8*m2*pow(l1k1,3) - 4*Q2h*pow(l1k1,3) +
		        8*l1k1*m2*pow(l1k2,2) - 4*l1k1*Q2h*pow(l1k2,2) +
		        4*m2*Q2h*pow(l1k2,2) - 4*m2*Q2k*pow(l1k2,2) +
		        2*Q2h*Q2k*pow(l1k2,2) - 4*l1k1*m2*pow(Q2h,2) -
		        4*l1k2*m2*pow(Q2h,2) + 4*m4*pow(Q2h,2) -
		        2*pow(l1k1,2)*pow(Q2h,2) - 2*pow(l1k2,2)*pow(Q2h,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		      (-64*k1k2*l1k2*m4 + 64*l1k2*m6 + 32*k1k2*l1k2*m2*Q2e -
		        32*k1k2*m4*Q2e + 32*m6*Q2e - 32*l1k2*m4*Q2h -
		        16*k1k2*l1k2*Q2e*Q2h + 16*k1k2*m2*Q2e*Q2h + 16*l1k2*m2*Q2e*Q2h -
		        16*m4*Q2e*Q2h + 32*k1k2*m4*Q2k - 32*m6*Q2k -
		        16*l1k2*m2*Q2e*Q2k - 16*k1k2*m2*Q2h*Q2k + 16*l1k2*m2*Q2h*Q2k +
		        16*m4*Q2h*Q2k + 8*l1k2*Q2e*Q2h*Q2k - 32*l1k2*m2*pow(k1k2,2) -
		        16*m2*Q2e*pow(k1k2,2) + 16*l1k2*Q2h*pow(k1k2,2) +
		        8*Q2e*Q2h*pow(k1k2,2) + 16*m2*Q2k*pow(k1k2,2) -
		        8*Q2h*Q2k*pow(k1k2,2) + 64*k1k2*m2*pow(l1k2,2) -
		        32*m2*Q2e*pow(l1k2,2) - 32*k1k2*Q2h*pow(l1k2,2) +
		        64*m2*Q2h*pow(l1k2,2) + 16*Q2e*Q2h*pow(l1k2,2) -
		        32*m2*Q2k*pow(l1k2,2) + 16*Q2h*Q2k*pow(l1k2,2) -
		        64*m2*pow(l1k2,3) + 32*Q2h*pow(l1k2,3) + 16*k1k2*l1k2*pow(Q2h,2) -
		        16*l1k2*m2*pow(Q2h,2) - 8*l1k2*Q2e*pow(Q2h,2) -
		        8*l1k2*Q2k*pow(Q2h,2) - 32*pow(l1k2,2)*pow(Q2h,2) +
		        8*l1k2*pow(Q2h,3)) + pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (-112*k1k2*l1k2*m2 + 32*k1k2*m4 - 24*k1k2*l1k2*Q2e -
		        32*k1k2*m2*Q2e + 16*l1k2*m2*Q2e + 8*k1k2*l1k2*Q2h -
		        24*k1k2*m2*Q2h - 48*l1k2*m2*Q2h + 32*m4*Q2h - 12*k1k2*Q2e*Q2h -
		        40*l1k2*Q2e*Q2h - 4*m2*Q2e*Q2h + 24*k1k2*m2*Q2k +
		        32*l1k2*m2*Q2k - 32*m4*Q2k + 8*l1k2*Q2e*Q2k + 4*m2*Q2e*Q2k +
		        8*k1k2*Q2h*Q2k - 16*l1k2*Q2h*Q2k + 12*m2*Q2h*Q2k +
		        8*Q2e*pow(k1k2,2) + 80*m2*pow(l1k2,2) + 32*Q2e*pow(l1k2,2) -
		        80*Q2h*pow(l1k2,2) - 4*k1k2*pow(Q2e,2) + 12*l1k2*pow(Q2e,2) -
		        4*m2*pow(Q2e,2) - 2*Q2h*pow(Q2e,2) + 2*pow(Q2e,3) +
		        52*l1k2*pow(Q2h,2) - 8*m2*pow(Q2h,2) + 8*Q2e*pow(Q2h,2) -
		        8*pow(Q2h,3) + 2*Q2e*pow(Q2k,2) - 2*Q2h*pow(Q2k,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*
		      (24*l1k1*l1k2*m2 + 16*l1k1*m4 - 16*m6 - 12*l1k1*l1k2*Q2e +
		        12*l1k2*m2*Q2e - 8*m4*Q2e + 28*l1k1*l1k2*Q2h + 32*l1k1*m2*Q2h +
		        8*m4*Q2h + 4*l1k1*Q2e*Q2h - 26*l1k2*Q2e*Q2h - 6*m2*Q2e*Q2h -
		        32*l1k1*m2*Q2k - 4*l1k2*m2*Q2k + 8*m4*Q2k - 4*l1k1*Q2e*Q2k +
		        8*l1k2*Q2e*Q2k - 4*m2*Q2e*Q2k + 4*l1k1*Q2h*Q2k - 16*l1k2*Q2h*Q2k -
		        5*Q2e*Q2h*Q2k + 72*m2*pow(l1k1,2) + 8*Q2e*pow(l1k1,2) -
		        8*Q2h*pow(l1k1,2) + 40*m2*pow(l1k2,2) + 12*Q2e*pow(l1k2,2) -
		        36*Q2h*pow(l1k2,2) + 8*l1k2*pow(Q2e,2) - Q2h*pow(Q2e,2) +
		        pow(Q2e,3) - 8*l1k1*pow(Q2h,2) + 30*l1k2*pow(Q2h,2) +
		        10*m2*pow(Q2h,2) + 7*Q2e*pow(Q2h,2) + 5*Q2k*pow(Q2h,2) -
		        7*pow(Q2h,3) + Q2e*pow(Q2k,2) - Q2h*pow(Q2k,2)) +
		     pow(2*l1k1 + Q2e - Q2k,-1)*
		      (64*k1k2*l1k2*m2 + 32*k1k2*l1k2*Q2e - 8*k1k2*m2*Q2e + 16*m4*Q2e +
		        32*l1k2*m2*Q2h - 16*m4*Q2h + 4*k1k2*Q2e*Q2h + 12*l1k2*Q2e*Q2h -
		        28*m2*Q2e*Q2h + 8*k1k2*m2*Q2k - 32*l1k2*m2*Q2k + 32*m4*Q2k +
		        8*k1k2*Q2e*Q2k - 16*l1k2*Q2e*Q2k + 8*m2*Q2e*Q2k - 8*k1k2*Q2h*Q2k +
		        24*m2*Q2h*Q2k + 2*Q2e*Q2h*Q2k - 16*Q2e*pow(k1k2,2) -
		        64*m2*pow(l1k2,2) - 32*Q2e*pow(l1k2,2) - 4*k1k2*pow(Q2e,2) +
		        4*l1k2*pow(Q2e,2) + 12*m2*pow(Q2e,2) - 4*Q2h*pow(Q2e,2) +
		        4*Q2k*pow(Q2e,2) - 4*pow(Q2e,3) - 4*Q2e*pow(Q2h,2) +
		        2*Q2k*pow(Q2h,2) - 16*m2*pow(Q2k,2) - 6*Q2e*pow(Q2k,2) +
		        2*Q2h*pow(Q2k,2)) + pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		      pow(2*l1k1 + Q2e - Q2k,-1)*
		      (32*k1k2*m4*Q2e*Q2h - 32*m6*Q2e*Q2h - 32*k1k2*m4*Q2e*Q2k +
		        32*m6*Q2e*Q2k - 32*k1k2*m4*Q2h*Q2k + 32*m6*Q2h*Q2k +
		        16*k1k2*m2*Q2e*Q2h*Q2k - 16*m4*Q2e*Q2h*Q2k +
		        16*m2*Q2e*Q2h*pow(k1k2,2) - 16*m2*Q2e*Q2k*pow(k1k2,2) -
		        16*m2*Q2h*Q2k*pow(k1k2,2) + 8*Q2e*Q2h*Q2k*pow(k1k2,2) -
		        16*k1k2*m2*Q2e*pow(Q2h,2) + 16*m4*Q2e*pow(Q2h,2) +
		        16*k1k2*m2*Q2k*pow(Q2h,2) - 16*m4*Q2k*pow(Q2h,2) -
		        8*Q2e*pow(k1k2,2)*pow(Q2h,2) + 8*Q2k*pow(k1k2,2)*pow(Q2h,2) +
		        32*k1k2*m4*pow(Q2k,2) - 32*m6*pow(Q2k,2) -
		        16*k1k2*m2*Q2h*pow(Q2k,2) + 16*m4*Q2h*pow(Q2k,2) +
		        16*m2*pow(k1k2,2)*pow(Q2k,2) - 8*Q2h*pow(k1k2,2)*pow(Q2k,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (32*l1k2*m4*Q2e - 32*l1k2*m4*Q2h + 48*l1k2*m2*Q2e*Q2h -
		        32*m4*Q2e*Q2h - 16*l1k2*m2*Q2e*Q2k + 16*l1k2*m2*Q2h*Q2k +
		        32*l1k2*Q2e*Q2h*Q2k - 64*m2*Q2e*pow(l1k2,2) +
		        96*m2*Q2h*pow(l1k2,2) + 112*Q2e*Q2h*pow(l1k2,2) -
		        32*m2*Q2k*pow(l1k2,2) - 16*Q2e*Q2k*pow(l1k2,2) +
		        32*Q2h*Q2k*pow(l1k2,2) - 64*m2*pow(l1k2,3) - 32*Q2e*pow(l1k2,3) +
		        64*Q2h*pow(l1k2,3) - 8*l1k2*m2*pow(Q2e,2) + 16*m4*pow(Q2e,2) +
		        44*l1k2*Q2h*pow(Q2e,2) - 4*m2*Q2h*pow(Q2e,2) -
		        8*l1k2*Q2k*pow(Q2e,2) + 4*Q2h*Q2k*pow(Q2e,2) -
		        32*pow(l1k2,2)*pow(Q2e,2) - 12*l1k2*pow(Q2e,3) + 4*m2*pow(Q2e,3) +
		        4*Q2h*pow(Q2e,3) - 2*pow(Q2e,4) - 40*l1k2*m2*pow(Q2h,2) +
		        16*m4*pow(Q2h,2) - 80*l1k2*Q2e*pow(Q2h,2) - 4*m2*Q2e*pow(Q2h,2) -
		        24*l1k2*Q2k*pow(Q2h,2) - 8*Q2e*Q2k*pow(Q2h,2) -
		        96*pow(l1k2,2)*pow(Q2h,2) - 10*pow(Q2e,2)*pow(Q2h,2) +
		        48*l1k2*pow(Q2h,3) + 4*m2*pow(Q2h,3) + 16*Q2e*pow(Q2h,3) +
		        4*Q2k*pow(Q2h,3) - 8*pow(Q2h,4) - 4*l1k2*Q2e*pow(Q2k,2) +
		        4*l1k2*Q2h*pow(Q2k,2) + 4*Q2e*Q2h*pow(Q2k,2) -
		        2*pow(Q2e,2)*pow(Q2k,2) - 2*pow(Q2h,2)*pow(Q2k,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
		      (16*l1k2*m4*Q2e - 16*l1k2*m4*Q2h - 16*l1k2*m2*Q2e*Q2h +
		        16*m4*Q2e*Q2h - 12*l1k2*m2*Q2e*Q2k - 8*m4*Q2e*Q2k +
		        12*l1k2*m2*Q2h*Q2k + 8*m4*Q2h*Q2k + 4*l1k2*Q2e*Q2h*Q2k -
		        10*m2*Q2e*Q2h*Q2k - 4*Q2e*Q2h*pow(l1k2,2) +
		        8*l1k2*m2*pow(Q2e,2) - 8*m4*pow(Q2e,2) - 4*l1k2*Q2h*pow(Q2e,2) +
		        12*m2*Q2h*pow(Q2e,2) - 2*l1k2*Q2k*pow(Q2e,2) +
		        4*m2*Q2k*pow(Q2e,2) + Q2h*Q2k*pow(Q2e,2) +
		        4*pow(l1k2,2)*pow(Q2e,2) - 4*m2*pow(Q2e,3) - Q2h*pow(Q2e,3) +
		        pow(Q2e,4) + 8*l1k2*m2*pow(Q2h,2) - 8*m4*pow(Q2h,2) +
		        4*l1k2*Q2e*pow(Q2h,2) - 10*m2*Q2e*pow(Q2h,2) -
		        2*l1k2*Q2k*pow(Q2h,2) + 6*m2*Q2k*pow(Q2h,2) -
		        2*Q2e*Q2k*pow(Q2h,2) + pow(Q2e,2)*pow(Q2h,2) + 2*m2*pow(Q2h,3) -
		        Q2e*pow(Q2h,3) + Q2k*pow(Q2h,3) - 2*Q2e*Q2h*pow(Q2k,2) +
		        pow(Q2e,2)*pow(Q2k,2) + pow(Q2h,2)*pow(Q2k,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
		      (32*k1k2*m4*Q2e - 32*m6*Q2e - 32*k1k2*m4*Q2h +
		        16*k1k2*m2*Q2e*Q2h + 32*m6*Q2k - 40*k1k2*m2*Q2e*Q2k -
		        24*k1k2*m2*Q2h*Q2k + 48*m4*Q2h*Q2k - 12*k1k2*Q2e*Q2h*Q2k -
		        16*m2*Q2e*Q2h*Q2k + 16*m2*Q2e*pow(k1k2,2) -
		        16*Q2e*Q2h*pow(k1k2,2) - 16*m2*Q2k*pow(k1k2,2) +
		        8*Q2e*Q2k*pow(k1k2,2) + 8*Q2h*Q2k*pow(k1k2,2) +
		        8*k1k2*m2*pow(Q2e,2) + 4*k1k2*Q2h*pow(Q2e,2) +
		        12*m2*Q2h*pow(Q2e,2) - 4*k1k2*Q2k*pow(Q2e,2) -
		        8*m2*Q2k*pow(Q2e,2) - 2*Q2h*Q2k*pow(Q2e,2) + 2*Q2k*pow(Q2e,3) +
		        16*k1k2*m2*pow(Q2h,2) - 16*m4*pow(Q2h,2) +
		        12*k1k2*Q2e*pow(Q2h,2) - 8*k1k2*Q2k*pow(Q2h,2) -
		        8*m2*Q2k*pow(Q2h,2) + 4*m2*pow(Q2h,3) + 24*k1k2*m2*pow(Q2k,2) -
		        32*m4*pow(Q2k,2) + 8*m2*Q2e*pow(Q2k,2) + 8*k1k2*Q2h*pow(Q2k,2) +
		        8*m2*Q2h*pow(Q2k,2) + 2*Q2e*pow(Q2k,3) - 2*Q2h*pow(Q2k,3)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      pow(2*l1k1 + Q2e - Q2k,-1)*
		      (-32*m6*Q2e*Q2h - 32*m4*Q2e*Q2h*Q2k + 16*m6*pow(Q2e,2) -
		        24*m4*Q2h*pow(Q2e,2) + 8*m4*Q2k*pow(Q2e,2) +
		        2*m2*Q2h*Q2k*pow(Q2e,2) + 8*m4*pow(Q2e,3) - 8*m2*Q2h*pow(Q2e,3) +
		        2*m2*Q2k*pow(Q2e,3) + 2*Q2h*Q2k*pow(Q2e,3) - Q2k*pow(Q2e,4) +
		        16*m6*pow(Q2h,2) + 32*m4*Q2e*pow(Q2h,2) + 24*m4*Q2k*pow(Q2h,2) -
		        2*m2*Q2e*Q2k*pow(Q2h,2) + 14*m2*pow(Q2e,2)*pow(Q2h,2) -
		        Q2k*pow(Q2e,2)*pow(Q2h,2) - 16*m4*pow(Q2h,3) -
		        8*m2*Q2e*pow(Q2h,3) - 2*m2*Q2k*pow(Q2h,3) + 2*m2*pow(Q2h,4) +
		        8*m4*Q2e*pow(Q2k,2) - 8*m4*Q2h*pow(Q2k,2) -
		        2*m2*pow(Q2e,2)*pow(Q2k,2) + 2*m2*pow(Q2h,2)*pow(Q2k,2) +
		        2*Q2e*Q2h*pow(Q2k,3) - pow(Q2e,2)*pow(Q2k,3) - pow(Q2h,2)*pow(Q2k,3))\
		) + (G2 + G3)*(16*k1k2*m2*M2 - 8*l1k1*m2*M2 + 8*l1k2*m2*M2 +
		     12*k1k2*M2*Q2e - 4*l1k1*M2*Q2e - 4*l1k2*M2*Q2e + 8*m2*M2*Q2e -
		     4*k1k2*M2*Q2h - 8*l1k1*M2*Q2h + 12*l1k2*M2*Q2h + 4*m2*M2*Q2h +
		     8*M2*Q2e*Q2h - 12*m2*M2*Q2k - 6*M2*Q2e*Q2k - 12*m2*Q2h*Q2k +
		     4*M2*Q2h*Q2k + 8*l1k1*Q2h*S + 8*l1k2*Q2h*S - 12*Q2e*Q2h*S -
		     16*m2*Q2k*S + 12*Q2e*Q2k*S - 4*Q2h*Q2k*S - 24*m2*Q2h*Sk -
		     32*m2*S*Sk + 24*Q2e*S*Sk - 8*Q2h*S*Sk + 8*l1k1*Q2h*Sq2 +
		     16*m2*Q2h*Sq2 + 32*k1k2*S*Sq2 + 32*l1k1*S*Sq2 - 64*l1k2*S*Sq2 +
		     64*m2*S*Sq2 - 16*Q2e*S*Sq2 + 32*Q2h*S*Sq2 - 16*Q2k*S*Sq2 +
		     2*M2*pow(Q2e,2) + 8*l1k2*pow(Q2h,2) + 4*m2*pow(Q2h,2) -
		     12*M2*pow(Q2h,2) + 2*Q2k*pow(Q2h,2) - 4*S*pow(Q2h,2) +
		     4*Sk*pow(Q2h,2) - 8*Sq2*pow(Q2h,2) - 2*pow(Q2h,3) + 32*k1k2*pow(S,2) +
		     48*l1k1*pow(S,2) - 64*l1k2*pow(S,2) + 64*m2*pow(S,2) -
		     16*Q2e*pow(S,2) + 40*Q2h*pow(S,2) - 8*Q2k*pow(S,2) +
		     pow(k1k2 - l1k1 - l1k2,-2)*
		      (-8*l1k1*l1k2*m2*M2*Q2h + 8*l1k1*m4*M2*Q2h +
		        4*l1k1*m2*M2*Q2h*Q2k + 4*l1k2*m2*M2*Q2h*Q2k -
		        4*m4*M2*Q2h*Q2k + 16*l1k1*l1k2*m2*Q2h*S - 16*l1k1*m4*Q2h*S -
		        8*l1k1*m2*Q2h*Q2k*S - 8*l1k2*m2*Q2h*Q2k*S + 8*m4*Q2h*Q2k*S -
		        8*m2*M2*Q2h*pow(l1k1,2) + 2*M2*Q2h*Q2k*pow(l1k1,2) +
		        16*m2*Q2h*S*pow(l1k1,2) - 4*Q2h*Q2k*S*pow(l1k1,2) -
		        4*M2*Q2h*pow(l1k1,3) + 8*Q2h*S*pow(l1k1,3) -
		        4*l1k1*M2*Q2h*pow(l1k2,2) + 2*M2*Q2h*Q2k*pow(l1k2,2) +
		        8*l1k1*Q2h*S*pow(l1k2,2) - 4*Q2h*Q2k*S*pow(l1k2,2) -
		        4*l1k1*m2*M2*pow(Q2h,2) - 4*l1k2*m2*M2*pow(Q2h,2) +
		        4*m4*M2*pow(Q2h,2) + 8*l1k1*m2*S*pow(Q2h,2) +
		        8*l1k2*m2*S*pow(Q2h,2) - 8*m4*S*pow(Q2h,2) -
		        2*M2*pow(l1k1,2)*pow(Q2h,2) + 4*S*pow(l1k1,2)*pow(Q2h,2) -
		        2*M2*pow(l1k2,2)*pow(Q2h,2) + 4*S*pow(l1k2,2)*pow(Q2h,2) +
		        32*l1k1*l1k2*m2*pow(S,2) - 32*l1k1*m4*pow(S,2) +
		        16*l1k1*m2*Q2h*pow(S,2) + 16*l1k2*m2*Q2h*pow(S,2) -
		        16*m4*Q2h*pow(S,2) - 16*l1k1*m2*Q2k*pow(S,2) -
		        16*l1k2*m2*Q2k*pow(S,2) + 16*m4*Q2k*pow(S,2) +
		        32*m2*pow(l1k1,2)*pow(S,2) + 8*Q2h*pow(l1k1,2)*pow(S,2) -
		        8*Q2k*pow(l1k1,2)*pow(S,2) + 16*pow(l1k1,3)*pow(S,2) +
		        16*l1k1*pow(l1k2,2)*pow(S,2) + 8*Q2h*pow(l1k2,2)*pow(S,2) -
		        8*Q2k*pow(l1k2,2)*pow(S,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*
		      (24*l1k1*l1k2*m2*M2 + 16*l1k2*m4*M2 + 8*l1k1*l1k2*M2*Q2e +
		        8*l1k1*m2*M2*Q2e + 4*l1k2*m2*M2*Q2e + 8*l1k2*m4*Q2h +
		        8*l1k1*l1k2*M2*Q2h + 16*l1k2*m2*M2*Q2h + 8*m4*M2*Q2h +
		        8*l1k1*M2*Q2e*Q2h + 2*l1k2*M2*Q2e*Q2h + 6*m2*M2*Q2e*Q2h -
		        16*l1k1*m2*M2*Q2k - 4*l1k2*m2*M2*Q2k - 6*l1k1*M2*Q2e*Q2k -
		        2*l1k2*M2*Q2e*Q2k + 2*m2*M2*Q2e*Q2k - 12*l1k1*m2*Q2h*Q2k -
		        8*l1k2*m2*Q2h*Q2k + 8*m4*Q2h*Q2k + 6*l1k1*M2*Q2h*Q2k -
		        6*l1k2*M2*Q2h*Q2k - 6*m2*M2*Q2h*Q2k + 2*l1k2*Q2e*Q2h*Q2k +
		        m2*Q2e*Q2h*Q2k - 7*M2*Q2e*Q2h*Q2k + 8*l1k1*l1k2*Q2h*S +
		        48*l1k1*m2*Q2h*S + 16*l1k2*m2*Q2h*S - 32*m4*Q2h*S +
		        4*l1k1*Q2e*Q2h*S - 8*l1k2*Q2e*Q2h*S - 4*m2*Q2e*Q2h*S -
		        16*l1k1*m2*Q2k*S - 32*l1k2*m2*Q2k*S + 32*m4*Q2k*S +
		        12*l1k1*Q2e*Q2k*S + 12*l1k2*Q2e*Q2k*S + 16*m2*Q2e*Q2k*S -
		        8*l1k1*Q2h*Q2k*S - 8*l1k2*Q2h*Q2k*S - 8*m2*Q2h*Q2k*S +
		        8*Q2e*Q2h*Q2k*S - 16*l1k2*m4*Sk - 24*l1k1*m2*Q2h*Sk -
		        16*l1k2*m2*Q2h*Sk + 16*m4*Q2h*Sk + 4*l1k2*Q2e*Q2h*Sk +
		        2*m2*Q2e*Q2h*Sk - 8*l1k2*m2*Q2k*Sk + 12*m2*Q2e*Q2k*Sk -
		        12*m2*Q2h*Q2k*Sk - 2*Q2e*Q2h*Q2k*Sk - 32*l1k1*m2*S*Sk -
		        64*l1k2*m2*S*Sk + 64*m4*S*Sk + 24*l1k1*Q2e*S*Sk +
		        24*l1k2*Q2e*S*Sk + 32*m2*Q2e*S*Sk - 8*l1k1*Q2h*S*Sk -
		        8*l1k2*Q2h*S*Sk + 24*Q2e*Q2h*S*Sk - 8*l1k1*l1k2*Q2h*Sq2 +
		        32*l1k1*m2*Q2h*Sq2 + 8*l1k2*m2*Q2h*Sq2 - 8*m4*Q2h*Sq2 +
		        8*l1k1*Q2e*Q2h*Sq2 + 4*l1k2*Q2e*Q2h*Sq2 - 8*m2*Q2e*Q2h*Sq2 -
		        4*l1k1*Q2h*Q2k*Sq2 + 4*l1k2*Q2h*Q2k*Sq2 + 2*Q2e*Q2h*Q2k*Sq2 -
		        32*l1k1*l1k2*S*Sq2 + 64*l1k1*m2*S*Sq2 - 32*m4*S*Sq2 -
		        16*l1k1*Q2e*S*Sq2 + 32*l1k2*Q2e*S*Sq2 - 32*m2*Q2e*S*Sq2 +
		        32*l1k1*Q2h*S*Sq2 - 64*l1k2*Q2h*S*Sq2 + 32*m2*Q2h*S*Sq2 -
		        44*Q2e*Q2h*S*Sq2 - 16*l1k1*Q2k*S*Sq2 + 16*l1k2*Q2k*S*Sq2 -
		        16*m2*Q2k*S*Sq2 - 4*Q2e*Q2k*S*Sq2 + 4*Q2h*Q2k*S*Sq2 +
		        16*m2*Q2h*Sk*Sq2 + 12*Q2e*Q2h*Sk*Sq2 + 24*m2*M2*pow(l1k1,2) +
		        12*M2*Q2e*pow(l1k1,2) - 12*M2*Q2h*pow(l1k1,2) +
		        16*Q2h*S*pow(l1k1,2) + 8*Q2h*Sq2*pow(l1k1,2) +
		        32*S*Sq2*pow(l1k1,2) + 8*m2*M2*pow(l1k2,2) -
		        4*M2*Q2e*pow(l1k2,2) - 20*M2*Q2h*pow(l1k2,2) +
		        8*Q2h*S*pow(l1k2,2) + 16*Q2h*Sq2*pow(l1k2,2) +
		        64*S*Sq2*pow(l1k2,2) + 2*l1k1*M2*pow(Q2e,2) -
		        4*l1k2*M2*pow(Q2e,2) - M2*Q2h*pow(Q2e,2) + M2*Q2k*pow(Q2e,2) +
		        Q2h*Q2k*pow(Q2e,2) - 4*Q2h*S*pow(Q2e,2) - 2*Q2k*S*pow(Q2e,2) +
		        2*Q2h*Sk*pow(Q2e,2) - 4*S*Sk*pow(Q2e,2) - 2*Q2h*Sq2*pow(Q2e,2) +
		        8*S*Sq2*pow(Q2e,2) + 8*l1k1*l1k2*pow(Q2h,2) +
		        4*l1k1*m2*pow(Q2h,2) + 8*l1k2*m2*pow(Q2h,2) - 4*m4*pow(Q2h,2) -
		        14*l1k1*M2*pow(Q2h,2) + 14*l1k2*M2*pow(Q2h,2) -
		        2*m2*M2*pow(Q2h,2) - 4*l1k2*Q2e*pow(Q2h,2) +
		        2*m2*Q2e*pow(Q2h,2) + 12*M2*Q2e*pow(Q2h,2) +
		        2*l1k1*Q2k*pow(Q2h,2) - 4*l1k2*Q2k*pow(Q2h,2) - m2*Q2k*pow(Q2h,2) +
		        6*M2*Q2k*pow(Q2h,2) - 3*Q2e*Q2k*pow(Q2h,2) + 4*l1k2*S*pow(Q2h,2) +
		        4*m2*S*pow(Q2h,2) + 6*Q2e*S*pow(Q2h,2) - 6*Q2k*S*pow(Q2h,2) +
		        4*l1k1*Sk*pow(Q2h,2) - 2*m2*Sk*pow(Q2h,2) - 4*Q2e*Sk*pow(Q2h,2) +
		        2*Q2k*Sk*pow(Q2h,2) - 20*S*Sk*pow(Q2h,2) - 4*l1k1*Sq2*pow(Q2h,2) -
		        12*l1k2*Sq2*pow(Q2h,2) + 12*m2*Sq2*pow(Q2h,2) +
		        8*Q2e*Sq2*pow(Q2h,2) - 2*Q2k*Sq2*pow(Q2h,2) + 36*S*Sq2*pow(Q2h,2) -
		        12*Sk*Sq2*pow(Q2h,2) - 8*pow(l1k2,2)*pow(Q2h,2) -
		        pow(Q2e,2)*pow(Q2h,2) - 2*l1k1*pow(Q2h,3) + 8*l1k2*pow(Q2h,3) -
		        2*m2*pow(Q2h,3) - 11*M2*pow(Q2h,3) + 3*Q2e*pow(Q2h,3) +
		        2*Q2k*pow(Q2h,3) - 2*S*pow(Q2h,3) + 2*Sk*pow(Q2h,3) -
		        6*Sq2*pow(Q2h,3) - 2*pow(Q2h,4) + 3*m2*Q2e*pow(Q2k,2) +
		        M2*Q2e*pow(Q2k,2) - 3*m2*Q2h*pow(Q2k,2) - M2*Q2h*pow(Q2k,2) -
		        32*l1k1*l1k2*pow(S,2) + 96*l1k1*m2*pow(S,2) -
		        32*l1k2*m2*pow(S,2) - 32*m4*pow(S,2) - 16*l1k1*Q2e*pow(S,2) +
		        32*l1k2*Q2e*pow(S,2) - 16*m2*Q2e*pow(S,2) + 48*l1k1*Q2h*pow(S,2) -
		        64*l1k2*Q2h*pow(S,2) + 32*m2*Q2h*pow(S,2) - 40*Q2e*Q2h*pow(S,2) -
		        16*l1k1*Q2k*pow(S,2) + 16*l1k2*Q2k*pow(S,2) - 16*m2*Q2k*pow(S,2) -
		        8*Q2e*Q2k*pow(S,2) + 8*Q2h*Q2k*pow(S,2) + 64*pow(l1k1,2)*pow(S,2) +
		        64*pow(l1k2,2)*pow(S,2) + 8*pow(Q2e,2)*pow(S,2) +
		        32*pow(Q2h,2)*pow(S,2) - 16*l1k2*m2*pow(Sk,2) +
		        12*m2*Q2e*pow(Sk,2) - 12*m2*Q2h*pow(Sk,2) - 4*Q2e*Q2h*pow(Sk,2) +
		        4*pow(Q2h,2)*pow(Sk,2) - 8*m2*Q2h*pow(Sq2,2) -
		        12*Q2e*Q2h*pow(Sq2,2) + 12*pow(Q2h,2)*pow(Sq2,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (-16*k1k2*l1k2*m2*M2 - 32*k1k2*m4*M2 + 32*l1k2*m4*M2 -
		        32*k1k2*l1k2*M2*Q2e - 16*k1k2*m2*M2*Q2e -
		        16*l1k2*m2*M2*Q2e + 16*m4*M2*Q2e - 16*k1k2*m4*Q2h +
		        16*l1k2*m4*Q2h + 16*k1k2*l1k2*M2*Q2h - 40*k1k2*m2*M2*Q2h +
		        32*l1k2*m2*M2*Q2h - 16*m4*M2*Q2h + 8*m4*Q2e*Q2h -
		        12*k1k2*M2*Q2e*Q2h - 8*l1k2*M2*Q2e*Q2h + 12*m2*M2*Q2e*Q2h +
		        8*k1k2*m2*M2*Q2k + 16*l1k2*m2*M2*Q2k + 4*k1k2*M2*Q2e*Q2k +
		        8*l1k2*M2*Q2e*Q2k + 12*m2*M2*Q2e*Q2k + 16*k1k2*m2*Q2h*Q2k +
		        8*l1k2*m2*Q2h*Q2k - 16*m4*Q2h*Q2k + 4*k1k2*M2*Q2h*Q2k -
		        16*l1k2*M2*Q2h*Q2k + 4*m2*M2*Q2h*Q2k - 4*k1k2*Q2e*Q2h*Q2k +
		        4*l1k2*Q2e*Q2h*Q2k + 4*M2*Q2e*Q2h*Q2k + 16*k1k2*l1k2*Q2h*S +
		        32*l1k2*m2*Q2h*S - 32*m4*Q2h*S - 24*k1k2*Q2e*Q2h*S +
		        24*l1k2*Q2e*Q2h*S - 64*k1k2*m2*Q2k*S + 32*l1k2*m2*Q2k*S +
		        64*m4*Q2k*S + 24*k1k2*Q2e*Q2k*S + 32*m2*Q2e*Q2k*S -
		        16*l1k2*Q2h*Q2k*S + 8*Q2e*Q2h*Q2k*S + 32*k1k2*m4*Sk -
		        32*l1k2*m4*Sk - 16*m4*Q2e*Sk + 32*k1k2*m2*Q2h*Sk +
		        16*l1k2*m2*Q2h*Sk - 16*m4*Q2h*Sk - 8*k1k2*Q2e*Q2h*Sk +
		        8*l1k2*Q2e*Q2h*Sk + 16*k1k2*m2*Q2k*Sk - 16*l1k2*m2*Q2k*Sk -
		        8*m2*Q2e*Q2k*Sk + 8*m2*Q2h*Q2k*Sk - 4*Q2e*Q2h*Q2k*Sk -
		        128*k1k2*m2*S*Sk + 64*l1k2*m2*S*Sk + 128*m4*S*Sk +
		        48*k1k2*Q2e*S*Sk + 64*m2*Q2e*S*Sk - 16*k1k2*Q2h*S*Sk +
		        32*Q2e*Q2h*S*Sk - 16*k1k2*m2*Q2h*Sq2 - 16*l1k2*m2*Q2h*Sq2 -
		        16*m4*Q2h*Sq2 - 24*k1k2*Q2e*Q2h*Sq2 + 16*l1k2*Q2e*Q2h*Sq2 -
		        24*m2*Q2e*Q2h*Sq2 - 64*k1k2*m2*Q2k*Sq2 + 32*l1k2*m2*Q2k*Sq2 +
		        64*m4*Q2k*Sq2 + 24*k1k2*Q2e*Q2k*Sq2 + 32*m2*Q2e*Q2k*Sq2 +
		        8*k1k2*Q2h*Q2k*Sq2 - 16*l1k2*Q2h*Q2k*Sq2 + 8*Q2e*Q2h*Q2k*Sq2 -
		        192*k1k2*l1k2*S*Sq2 + 256*k1k2*m2*S*Sq2 - 256*l1k2*m2*S*Sq2 -
		        64*m4*S*Sq2 - 32*k1k2*Q2e*S*Sq2 + 96*l1k2*Q2e*S*Sq2 -
		        128*m2*Q2e*S*Sq2 + 96*k1k2*Q2h*S*Sq2 - 256*l1k2*Q2h*S*Sq2 +
		        128*m2*Q2h*S*Sq2 - 104*Q2e*Q2h*S*Sq2 - 32*k1k2*Q2k*S*Sq2 +
		        64*l1k2*Q2k*S*Sq2 - 32*m2*Q2k*S*Sq2 - 8*Q2e*Q2k*S*Sq2 +
		        8*Q2h*Q2k*S*Sq2 - 128*k1k2*m2*Sk*Sq2 + 64*l1k2*m2*Sk*Sq2 +
		        128*m4*Sk*Sq2 + 48*k1k2*Q2e*Sk*Sq2 + 64*m2*Q2e*Sk*Sq2 -
		        16*k1k2*Q2h*Sk*Sq2 - 32*m2*Q2h*Sk*Sq2 + 24*Q2e*Q2h*Sk*Sq2 -
		        16*m2*M2*pow(k1k2,2) + 8*M2*Q2e*pow(k1k2,2) +
		        64*S*Sq2*pow(k1k2,2) - 16*m2*M2*pow(l1k2,2) -
		        48*M2*Q2h*pow(l1k2,2) - 32*Q2h*S*pow(l1k2,2) +
		        256*S*Sq2*pow(l1k2,2) - 16*l1k2*M2*pow(Q2e,2) -
		        4*m2*M2*pow(Q2e,2) - 2*M2*Q2k*pow(Q2e,2) -
		        2*Q2h*Q2k*pow(Q2e,2) + 4*Q2k*S*pow(Q2e,2) - 4*Q2h*Sk*pow(Q2e,2) +
		        8*S*Sk*pow(Q2e,2) + 4*Q2k*Sq2*pow(Q2e,2) + 16*S*Sq2*pow(Q2e,2) +
		        8*Sk*Sq2*pow(Q2e,2) - 2*M2*pow(Q2e,3) - 16*k1k2*m2*pow(Q2h,2) +
		        8*l1k2*m2*pow(Q2h,2) - 4*k1k2*M2*pow(Q2h,2) +
		        48*l1k2*M2*pow(Q2h,2) - 24*m2*M2*pow(Q2h,2) -
		        8*l1k2*Q2e*pow(Q2h,2) + 8*m2*Q2e*pow(Q2h,2) +
		        16*M2*Q2e*pow(Q2h,2) - 4*l1k2*Q2k*pow(Q2h,2) -
		        2*M2*Q2k*pow(Q2h,2) + 2*Q2e*Q2k*pow(Q2h,2) + 8*l1k2*S*pow(Q2h,2) -
		        16*m2*S*pow(Q2h,2) - 12*Q2k*S*pow(Q2h,2) - 8*l1k2*Sk*pow(Q2h,2) +
		        4*Q2e*Sk*pow(Q2h,2) + 4*Q2k*Sk*pow(Q2h,2) - 40*S*Sk*pow(Q2h,2) +
		        8*l1k2*Sq2*pow(Q2h,2) + 4*Q2e*Sq2*pow(Q2h,2) -
		        12*Q2k*Sq2*pow(Q2h,2) + 88*S*Sq2*pow(Q2h,2) - 32*Sk*Sq2*pow(Q2h,2) -
		        16*pow(l1k2,2)*pow(Q2h,2) + 12*l1k2*pow(Q2h,3) - 8*m2*pow(Q2h,3) -
		        14*M2*pow(Q2h,3) + 2*Q2e*pow(Q2h,3) - 4*Sq2*pow(Q2h,3) -
		        2*pow(Q2h,4) + 2*M2*Q2e*pow(Q2k,2) - 2*M2*Q2h*pow(Q2k,2) -
		        128*k1k2*l1k2*pow(S,2) + 128*k1k2*m2*pow(S,2) -
		        192*l1k2*m2*pow(S,2) - 32*k1k2*Q2e*pow(S,2) +
		        96*l1k2*Q2e*pow(S,2) - 96*m2*Q2e*pow(S,2) + 64*k1k2*Q2h*pow(S,2) -
		        192*l1k2*Q2h*pow(S,2) + 96*m2*Q2h*pow(S,2) - 80*Q2e*Q2h*pow(S,2) +
		        32*l1k2*Q2k*pow(S,2) + 32*pow(k1k2,2)*pow(S,2) +
		        192*pow(l1k2,2)*pow(S,2) + 16*pow(Q2e,2)*pow(S,2) +
		        64*pow(Q2h,2)*pow(S,2) + 32*k1k2*m2*pow(Sk,2) -
		        32*l1k2*m2*pow(Sk,2) - 16*m2*Q2e*pow(Sk,2) +
		        16*m2*Q2h*pow(Sk,2) - 8*Q2e*Q2h*pow(Sk,2) +
		        8*pow(Q2h,2)*pow(Sk,2) - 64*k1k2*l1k2*pow(Sq2,2) +
		        128*k1k2*m2*pow(Sq2,2) - 64*l1k2*m2*pow(Sq2,2) -
		        64*m4*pow(Sq2,2) - 32*m2*Q2e*pow(Sq2,2) + 32*k1k2*Q2h*pow(Sq2,2) -
		        64*l1k2*Q2h*pow(Sq2,2) + 48*m2*Q2h*pow(Sq2,2) -
		        32*Q2e*Q2h*pow(Sq2,2) - 32*k1k2*Q2k*pow(Sq2,2) +
		        32*l1k2*Q2k*pow(Sq2,2) - 32*m2*Q2k*pow(Sq2,2) -
		        8*Q2e*Q2k*pow(Sq2,2) + 8*Q2h*Q2k*pow(Sq2,2) +
		        32*pow(k1k2,2)*pow(Sq2,2) + 64*pow(l1k2,2)*pow(Sq2,2) +
		        32*pow(Q2h,2)*pow(Sq2,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		      (32*k1k2*l1k2*m2*M2*Q2h - 32*l1k2*m4*M2*Q2h -
		        16*k1k2*l1k2*M2*Q2e*Q2h + 16*k1k2*m2*M2*Q2e*Q2h -
		        16*m4*M2*Q2e*Q2h - 16*k1k2*m2*M2*Q2h*Q2k +
		        16*m4*M2*Q2h*Q2k + 8*l1k2*M2*Q2e*Q2h*Q2k +
		        64*k1k2*l1k2*m2*Q2h*S - 64*l1k2*m4*Q2h*S -
		        32*k1k2*l1k2*Q2e*Q2h*S + 32*k1k2*m2*Q2e*Q2h*S - 32*m4*Q2e*Q2h*S -
		        32*k1k2*m2*Q2h*Q2k*S + 32*m4*Q2h*Q2k*S + 16*l1k2*Q2e*Q2h*Q2k*S +
		        64*k1k2*l1k2*m2*Q2h*Sq2 - 64*l1k2*m4*Q2h*Sq2 -
		        32*k1k2*l1k2*Q2e*Q2h*Sq2 + 32*k1k2*m2*Q2e*Q2h*Sq2 -
		        32*m4*Q2e*Q2h*Sq2 - 32*k1k2*m2*Q2h*Q2k*Sq2 + 32*m4*Q2h*Q2k*Sq2 +
		        16*l1k2*Q2e*Q2h*Q2k*Sq2 - 256*k1k2*l1k2*m2*S*Sq2 +
		        256*l1k2*m4*S*Sq2 + 128*k1k2*l1k2*Q2e*S*Sq2 -
		        128*k1k2*m2*Q2e*S*Sq2 + 128*m4*Q2e*S*Sq2 -
		        128*k1k2*l1k2*Q2h*S*Sq2 + 64*l1k2*Q2e*Q2h*S*Sq2 +
		        128*k1k2*m2*Q2k*S*Sq2 - 128*m4*Q2k*S*Sq2 - 64*l1k2*Q2e*Q2k*S*Sq2 +
		        64*l1k2*Q2h*Q2k*S*Sq2 + 16*l1k2*M2*Q2h*pow(k1k2,2) +
		        8*M2*Q2e*Q2h*pow(k1k2,2) - 8*M2*Q2h*Q2k*pow(k1k2,2) +
		        32*l1k2*Q2h*S*pow(k1k2,2) + 16*Q2e*Q2h*S*pow(k1k2,2) -
		        16*Q2h*Q2k*S*pow(k1k2,2) + 32*l1k2*Q2h*Sq2*pow(k1k2,2) +
		        16*Q2e*Q2h*Sq2*pow(k1k2,2) - 16*Q2h*Q2k*Sq2*pow(k1k2,2) -
		        128*l1k2*S*Sq2*pow(k1k2,2) - 64*Q2e*S*Sq2*pow(k1k2,2) +
		        64*Q2k*S*Sq2*pow(k1k2,2) - 32*k1k2*M2*Q2h*pow(l1k2,2) +
		        16*M2*Q2e*Q2h*pow(l1k2,2) + 16*M2*Q2h*Q2k*pow(l1k2,2) -
		        64*k1k2*Q2h*S*pow(l1k2,2) + 32*Q2e*Q2h*S*pow(l1k2,2) +
		        32*Q2h*Q2k*S*pow(l1k2,2) - 64*k1k2*Q2h*Sq2*pow(l1k2,2) +
		        32*Q2e*Q2h*Sq2*pow(l1k2,2) + 32*Q2h*Q2k*Sq2*pow(l1k2,2) +
		        256*k1k2*S*Sq2*pow(l1k2,2) - 128*Q2e*S*Sq2*pow(l1k2,2) +
		        256*Q2h*S*Sq2*pow(l1k2,2) - 128*Q2k*S*Sq2*pow(l1k2,2) +
		        32*M2*Q2h*pow(l1k2,3) + 64*Q2h*S*pow(l1k2,3) +
		        64*Q2h*Sq2*pow(l1k2,3) - 256*S*Sq2*pow(l1k2,3) +
		        16*k1k2*l1k2*M2*pow(Q2h,2) - 8*l1k2*M2*Q2e*pow(Q2h,2) -
		        8*l1k2*M2*Q2k*pow(Q2h,2) + 32*k1k2*l1k2*S*pow(Q2h,2) -
		        16*l1k2*Q2e*S*pow(Q2h,2) - 16*l1k2*Q2k*S*pow(Q2h,2) +
		        32*k1k2*l1k2*Sq2*pow(Q2h,2) - 16*l1k2*Q2e*Sq2*pow(Q2h,2) -
		        16*l1k2*Q2k*Sq2*pow(Q2h,2) - 64*l1k2*S*Sq2*pow(Q2h,2) -
		        32*M2*pow(l1k2,2)*pow(Q2h,2) - 64*S*pow(l1k2,2)*pow(Q2h,2) -
		        64*Sq2*pow(l1k2,2)*pow(Q2h,2) + 8*l1k2*M2*pow(Q2h,3) +
		        16*l1k2*S*pow(Q2h,3) + 16*l1k2*Sq2*pow(Q2h,3) -
		        128*k1k2*l1k2*m2*pow(S,2) + 128*l1k2*m4*pow(S,2) +
		        64*k1k2*l1k2*Q2e*pow(S,2) - 64*k1k2*m2*Q2e*pow(S,2) +
		        64*m4*Q2e*pow(S,2) - 64*k1k2*l1k2*Q2h*pow(S,2) +
		        32*l1k2*Q2e*Q2h*pow(S,2) + 64*k1k2*m2*Q2k*pow(S,2) -
		        64*m4*Q2k*pow(S,2) - 32*l1k2*Q2e*Q2k*pow(S,2) +
		        32*l1k2*Q2h*Q2k*pow(S,2) - 64*l1k2*pow(k1k2,2)*pow(S,2) -
		        32*Q2e*pow(k1k2,2)*pow(S,2) + 32*Q2k*pow(k1k2,2)*pow(S,2) +
		        128*k1k2*pow(l1k2,2)*pow(S,2) - 64*Q2e*pow(l1k2,2)*pow(S,2) +
		        128*Q2h*pow(l1k2,2)*pow(S,2) - 64*Q2k*pow(l1k2,2)*pow(S,2) -
		        128*pow(l1k2,3)*pow(S,2) - 32*l1k2*pow(Q2h,2)*pow(S,2) -
		        128*k1k2*l1k2*m2*pow(Sq2,2) + 128*l1k2*m4*pow(Sq2,2) +
		        64*k1k2*l1k2*Q2e*pow(Sq2,2) - 64*k1k2*m2*Q2e*pow(Sq2,2) +
		        64*m4*Q2e*pow(Sq2,2) - 64*k1k2*l1k2*Q2h*pow(Sq2,2) +
		        32*l1k2*Q2e*Q2h*pow(Sq2,2) + 64*k1k2*m2*Q2k*pow(Sq2,2) -
		        64*m4*Q2k*pow(Sq2,2) - 32*l1k2*Q2e*Q2k*pow(Sq2,2) +
		        32*l1k2*Q2h*Q2k*pow(Sq2,2) - 64*l1k2*pow(k1k2,2)*pow(Sq2,2) -
		        32*Q2e*pow(k1k2,2)*pow(Sq2,2) + 32*Q2k*pow(k1k2,2)*pow(Sq2,2) +
		        128*k1k2*pow(l1k2,2)*pow(Sq2,2) - 64*Q2e*pow(l1k2,2)*pow(Sq2,2) +
		        128*Q2h*pow(l1k2,2)*pow(Sq2,2) - 64*Q2k*pow(l1k2,2)*pow(Sq2,2) -
		        128*pow(l1k2,3)*pow(Sq2,2) - 32*l1k2*pow(Q2h,2)*pow(Sq2,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
		      (-8*l1k2*m2*M2*Q2e*Q2h + 8*l1k2*m2*M2*Q2e*Q2k -
		        8*l1k2*m2*M2*Q2h*Q2k - 2*l1k2*m2*Q2e*Q2h*Q2k -
		        4*m4*Q2e*Q2h*Q2k + 4*l1k2*M2*Q2e*Q2h*Q2k +
		        12*m2*M2*Q2e*Q2h*Q2k + 8*l1k2*m2*Q2e*Q2h*S + 16*m4*Q2e*Q2h*S +
		        16*l1k2*m2*Q2e*Q2k*S - 16*m4*Q2e*Q2k*S - 16*l1k2*m2*Q2h*Q2k*S +
		        16*m4*Q2h*Q2k*S + 8*l1k2*Q2e*Q2h*Q2k*S + 12*m2*Q2e*Q2h*Q2k*S -
		        4*l1k2*m2*Q2e*Q2h*Sk - 8*m4*Q2e*Q2h*Sk + 24*l1k2*m2*Q2e*Q2k*Sk -
		        24*l1k2*m2*Q2h*Q2k*Sk + 6*m2*Q2e*Q2h*Q2k*Sk +
		        32*l1k2*m2*Q2e*S*Sk - 32*m4*Q2e*S*Sk - 32*l1k2*m2*Q2h*S*Sk +
		        32*m4*Q2h*S*Sk + 16*l1k2*Q2e*Q2h*S*Sk + 16*m2*Q2e*Q2h*S*Sk +
		        8*Q2e*Q2h*Q2k*S*Sk + 10*m2*Q2e*Q2h*Q2k*Sq2 - 16*m2*Q2e*Q2k*S*Sq2 +
		        16*m2*Q2h*Q2k*S*Sq2 - 24*Q2e*Q2h*Q2k*S*Sq2 -
		        24*m2*Q2e*Q2k*Sk*Sq2 + 24*m2*Q2h*Q2k*Sk*Sq2 -
		        4*M2*Q2e*Q2h*pow(l1k2,2) - 8*l1k2*M2*Q2h*pow(Q2e,2) -
		        2*l1k2*M2*Q2k*pow(Q2e,2) - 4*m2*M2*Q2k*pow(Q2e,2) -
		        2*l1k2*Q2h*Q2k*pow(Q2e,2) - 5*m2*Q2h*Q2k*pow(Q2e,2) +
		        5*M2*Q2h*Q2k*pow(Q2e,2) + 4*l1k2*Q2h*S*pow(Q2e,2) +
		        4*m2*Q2h*S*pow(Q2e,2) - 8*l1k2*Q2k*S*pow(Q2e,2) -
		        8*m2*Q2k*S*pow(Q2e,2) - 8*Q2h*Q2k*S*pow(Q2e,2) -
		        4*l1k2*Q2h*Sk*pow(Q2e,2) - 10*m2*Q2h*Sk*pow(Q2e,2) +
		        2*Q2h*Q2k*Sk*pow(Q2e,2) - 16*l1k2*S*Sk*pow(Q2e,2) -
		        16*m2*S*Sk*pow(Q2e,2) - 8*Q2h*S*Sk*pow(Q2e,2) -
		        4*Q2k*S*Sk*pow(Q2e,2) + 4*l1k2*Q2h*Sq2*pow(Q2e,2) -
		        2*Q2h*Q2k*Sq2*pow(Q2e,2) + 8*Q2k*S*Sq2*pow(Q2e,2) +
		        4*M2*pow(l1k2,2)*pow(Q2e,2) + 4*l1k2*M2*pow(Q2e,3) -
		        M2*Q2h*pow(Q2e,3) - M2*Q2k*pow(Q2e,3) - Q2h*Q2k*pow(Q2e,3) +
		        4*Q2h*S*pow(Q2e,3) + 2*Q2k*S*pow(Q2e,3) - 2*Q2h*Sk*pow(Q2e,3) +
		        4*S*Sk*pow(Q2e,3) + 4*Q2h*Sq2*pow(Q2e,3) +
		        8*l1k2*m2*M2*pow(Q2h,2) + 4*m4*Q2e*pow(Q2h,2) +
		        4*l1k2*M2*Q2e*pow(Q2h,2) - 2*m2*M2*Q2e*pow(Q2h,2) +
		        2*l1k2*m2*Q2k*pow(Q2h,2) + 4*m4*Q2k*pow(Q2h,2) -
		        2*l1k2*M2*Q2k*pow(Q2h,2) - 8*m2*M2*Q2k*pow(Q2h,2) +
		        2*l1k2*Q2e*Q2k*pow(Q2h,2) + 4*m2*Q2e*Q2k*pow(Q2h,2) -
		        5*M2*Q2e*Q2k*pow(Q2h,2) - 8*l1k2*m2*S*pow(Q2h,2) -
		        16*m4*S*pow(Q2h,2) - 4*l1k2*Q2e*S*pow(Q2h,2) -
		        4*m2*Q2e*S*pow(Q2h,2) - 4*m2*Q2k*S*pow(Q2h,2) +
		        8*Q2e*Q2k*S*pow(Q2h,2) + 4*l1k2*m2*Sk*pow(Q2h,2) +
		        8*m4*Sk*pow(Q2h,2) + 4*l1k2*Q2e*Sk*pow(Q2h,2) +
		        10*m2*Q2e*Sk*pow(Q2h,2) - 6*m2*Q2k*Sk*pow(Q2h,2) -
		        4*Q2e*Q2k*Sk*pow(Q2h,2) + 4*Q2e*S*Sk*pow(Q2h,2) -
		        4*Q2k*S*Sk*pow(Q2h,2) - 4*l1k2*Q2e*Sq2*pow(Q2h,2) -
		        2*m2*Q2e*Sq2*pow(Q2h,2) - 10*m2*Q2k*Sq2*pow(Q2h,2) +
		        4*Q2e*Q2k*Sq2*pow(Q2h,2) + 16*Q2k*S*Sq2*pow(Q2h,2) +
		        3*m2*pow(Q2e,2)*pow(Q2h,2) + 2*M2*pow(Q2e,2)*pow(Q2h,2) +
		        2*Q2k*pow(Q2e,2)*pow(Q2h,2) - 6*S*pow(Q2e,2)*pow(Q2h,2) +
		        4*Sk*pow(Q2e,2)*pow(Q2h,2) - 6*Sq2*pow(Q2e,2)*pow(Q2h,2) -
		        4*m4*pow(Q2h,3) + 2*m2*M2*pow(Q2h,3) - 3*m2*Q2e*pow(Q2h,3) -
		        M2*Q2e*pow(Q2h,3) + m2*Q2k*pow(Q2h,3) + M2*Q2k*pow(Q2h,3) -
		        Q2e*Q2k*pow(Q2h,3) + 2*Q2e*S*pow(Q2h,3) - 2*Q2k*S*pow(Q2h,3) -
		        2*Q2e*Sk*pow(Q2h,3) + 2*Q2k*Sk*pow(Q2h,3) + 2*m2*Sq2*pow(Q2h,3) +
		        2*Q2e*Sq2*pow(Q2h,3) - 2*Q2k*Sq2*pow(Q2h,3) +
		        6*l1k2*m2*Q2e*pow(Q2k,2) + 2*m2*M2*Q2e*pow(Q2k,2) -
		        6*l1k2*m2*Q2h*pow(Q2k,2) - 2*m2*M2*Q2h*pow(Q2k,2) +
		        2*m2*Q2e*Q2h*pow(Q2k,2) - 2*M2*Q2e*Q2h*pow(Q2k,2) +
		        4*Q2e*Q2h*S*pow(Q2k,2) + 12*m2*Q2e*Sk*pow(Q2k,2) -
		        12*m2*Q2h*Sk*pow(Q2k,2) - 12*m2*Q2e*Sq2*pow(Q2k,2) +
		        12*m2*Q2h*Sq2*pow(Q2k,2) + M2*pow(Q2e,2)*pow(Q2k,2) +
		        Q2h*pow(Q2e,2)*pow(Q2k,2) - 2*S*pow(Q2e,2)*pow(Q2k,2) -
		        2*m2*pow(Q2h,2)*pow(Q2k,2) + M2*pow(Q2h,2)*pow(Q2k,2) -
		        2*Q2e*pow(Q2h,2)*pow(Q2k,2) - 2*S*pow(Q2h,2)*pow(Q2k,2) +
		        pow(Q2h,3)*pow(Q2k,2) + 3*m2*Q2e*pow(Q2k,3) -
		        3*m2*Q2h*pow(Q2k,3) + 32*l1k2*m2*Q2e*pow(S,2) -
		        32*l1k2*m2*Q2h*pow(S,2) - 16*l1k2*Q2e*Q2k*pow(S,2) -
		        16*m2*Q2e*Q2k*pow(S,2) + 16*l1k2*Q2h*Q2k*pow(S,2) +
		        16*m2*Q2h*Q2k*pow(S,2) - 32*Q2e*Q2h*Q2k*pow(S,2) +
		        4*Q2h*pow(Q2e,2)*pow(S,2) + 12*Q2k*pow(Q2e,2)*pow(S,2) -
		        4*Q2e*pow(Q2h,2)*pow(S,2) + 20*Q2k*pow(Q2h,2)*pow(S,2) -
		        4*Q2e*pow(Q2k,2)*pow(S,2) + 4*Q2h*pow(Q2k,2)*pow(S,2) +
		        24*l1k2*m2*Q2e*pow(Sk,2) - 24*l1k2*m2*Q2h*pow(Sk,2) +
		        4*m2*Q2e*Q2h*pow(Sk,2) + 12*m2*Q2e*Q2k*pow(Sk,2) -
		        12*m2*Q2h*Q2k*pow(Sk,2) - 4*m2*pow(Q2h,2)*pow(Sk,2) -
		        4*Q2h*pow(Q2e,2)*pow(Sq2,2) + 4*Q2e*pow(Q2h,2)*pow(Sq2,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (-16*l1k2*m2*M2*Q2e*Q2h + 16*l1k2*M2*Q2e*Q2h*Q2k +
		        8*l1k2*Q2e*Q2h*Q2k*Sk + 16*l1k2*m2*Q2e*Q2h*Sq2 -
		        12*l1k2*Q2e*Q2h*Q2k*Sq2 + 64*l1k2*m2*Q2e*S*Sq2 -
		        64*l1k2*m2*Q2h*S*Sq2 + 192*l1k2*Q2e*Q2h*S*Sq2 -
		        64*m2*Q2e*Q2h*S*Sq2 - 32*l1k2*Q2e*Q2k*S*Sq2 +
		        32*l1k2*Q2h*Q2k*S*Sq2 - 16*l1k2*Q2e*Q2h*Sk*Sq2 +
		        32*M2*Q2e*Q2h*pow(l1k2,2) + 16*M2*Q2h*Q2k*pow(l1k2,2) -
		        32*Q2e*Q2h*Sq2*pow(l1k2,2) - 16*Q2h*Q2k*Sq2*pow(l1k2,2) -
		        128*Q2e*S*Sq2*pow(l1k2,2) + 192*Q2h*S*Sq2*pow(l1k2,2) -
		        64*Q2k*S*Sq2*pow(l1k2,2) + 32*M2*Q2h*pow(l1k2,3) -
		        32*Q2h*Sq2*pow(l1k2,3) - 128*S*Sq2*pow(l1k2,3) +
		        12*l1k2*M2*Q2h*pow(Q2e,2) - 8*m2*M2*Q2h*pow(Q2e,2) +
		        4*M2*Q2h*Q2k*pow(Q2e,2) - 8*Q2h*Q2k*S*pow(Q2e,2) +
		        4*Q2h*Q2k*Sk*pow(Q2e,2) - 32*Q2h*S*Sk*pow(Q2e,2) -
		        12*l1k2*Q2h*Sq2*pow(Q2e,2) + 8*m2*Q2h*Sq2*pow(Q2e,2) -
		        6*Q2h*Q2k*Sq2*pow(Q2e,2) - 48*l1k2*S*Sq2*pow(Q2e,2) +
		        32*m2*S*Sq2*pow(Q2e,2) + 72*Q2h*S*Sq2*pow(Q2e,2) -
		        24*Q2h*Sk*Sq2*pow(Q2e,2) + 2*M2*Q2h*pow(Q2e,3) -
		        2*Q2h*Sq2*pow(Q2e,3) - 8*S*Sq2*pow(Q2e,3) +
		        16*l1k2*m2*M2*pow(Q2h,2) - 8*l1k2*m2*Q2e*pow(Q2h,2) -
		        52*l1k2*M2*Q2e*pow(Q2h,2) + 16*m2*M2*Q2e*pow(Q2h,2) -
		        16*l1k2*M2*Q2k*pow(Q2h,2) + 6*l1k2*Q2e*Q2k*pow(Q2h,2) -
		        8*M2*Q2e*Q2k*pow(Q2h,2) + 16*Q2e*Q2k*S*pow(Q2h,2) -
		        8*l1k2*Q2k*Sk*pow(Q2h,2) - 8*Q2e*Q2k*Sk*pow(Q2h,2) +
		        64*Q2e*S*Sk*pow(Q2h,2) - 16*l1k2*m2*Sq2*pow(Q2h,2) +
		        20*l1k2*Q2e*Sq2*pow(Q2h,2) - 16*m2*Q2e*Sq2*pow(Q2h,2) +
		        12*l1k2*Q2k*Sq2*pow(Q2h,2) + 12*Q2e*Q2k*Sq2*pow(Q2h,2) -
		        144*l1k2*S*Sq2*pow(Q2h,2) + 32*m2*S*Sq2*pow(Q2h,2) -
		        120*Q2e*S*Sq2*pow(Q2h,2) + 16*l1k2*Sk*Sq2*pow(Q2h,2) +
		        48*Q2e*Sk*Sq2*pow(Q2h,2) - 48*M2*pow(l1k2,2)*pow(Q2h,2) +
		        16*Q2e*pow(l1k2,2)*pow(Q2h,2) + 8*Q2k*pow(l1k2,2)*pow(Q2h,2) +
		        48*Sq2*pow(l1k2,2)*pow(Q2h,2) + 16*pow(l1k2,3)*pow(Q2h,2) +
		        6*l1k2*pow(Q2e,2)*pow(Q2h,2) - 4*m2*pow(Q2e,2)*pow(Q2h,2) -
		        16*M2*pow(Q2e,2)*pow(Q2h,2) + Q2k*pow(Q2e,2)*pow(Q2h,2) +
		        pow(Q2e,3)*pow(Q2h,2) + 8*l1k2*m2*pow(Q2h,3) +
		        40*l1k2*M2*pow(Q2h,3) - 8*m2*M2*pow(Q2h,3) -
		        18*l1k2*Q2e*pow(Q2h,3) + 8*m2*Q2e*pow(Q2h,3) +
		        26*M2*Q2e*pow(Q2h,3) - 6*l1k2*Q2k*pow(Q2h,3) +
		        4*M2*Q2k*pow(Q2h,3) - 2*Q2e*Q2k*pow(Q2h,3) - 8*Q2k*S*pow(Q2h,3) +
		        4*Q2k*Sk*pow(Q2h,3) - 32*S*Sk*pow(Q2h,3) - 8*l1k2*Sq2*pow(Q2h,3) +
		        8*m2*Sq2*pow(Q2h,3) + 6*Q2e*Sq2*pow(Q2h,3) - 6*Q2k*Sq2*pow(Q2h,3) +
		        56*S*Sq2*pow(Q2h,3) - 24*Sk*Sq2*pow(Q2h,3) -
		        24*pow(l1k2,2)*pow(Q2h,3) - 4*pow(Q2e,2)*pow(Q2h,3) +
		        12*l1k2*pow(Q2h,4) - 4*m2*pow(Q2h,4) - 12*M2*pow(Q2h,4) +
		        5*Q2e*pow(Q2h,4) + Q2k*pow(Q2h,4) - 4*Sq2*pow(Q2h,4) -
		        2*pow(Q2h,5) - 4*l1k2*M2*Q2e*pow(Q2k,2) +
		        4*l1k2*M2*Q2h*pow(Q2k,2) + 4*M2*Q2e*Q2h*pow(Q2k,2) -
		        2*M2*pow(Q2e,2)*pow(Q2k,2) - 2*M2*pow(Q2h,2)*pow(Q2k,2) +
		        64*l1k2*m2*Q2e*pow(S,2) - 64*l1k2*m2*Q2h*pow(S,2) +
		        192*l1k2*Q2e*Q2h*pow(S,2) - 64*m2*Q2e*Q2h*pow(S,2) -
		        32*l1k2*Q2e*Q2k*pow(S,2) + 32*l1k2*Q2h*Q2k*pow(S,2) -
		        128*Q2e*pow(l1k2,2)*pow(S,2) + 192*Q2h*pow(l1k2,2)*pow(S,2) -
		        64*Q2k*pow(l1k2,2)*pow(S,2) - 128*pow(l1k2,3)*pow(S,2) -
		        48*l1k2*pow(Q2e,2)*pow(S,2) + 32*m2*pow(Q2e,2)*pow(S,2) +
		        56*Q2h*pow(Q2e,2)*pow(S,2) - 8*pow(Q2e,3)*pow(S,2) -
		        144*l1k2*pow(Q2h,2)*pow(S,2) + 32*m2*pow(Q2h,2)*pow(S,2) -
		        88*Q2e*pow(Q2h,2)*pow(S,2) + 40*pow(Q2h,3)*pow(S,2) +
		        16*l1k2*Q2e*Q2h*pow(Sk,2) + 8*Q2h*pow(Q2e,2)*pow(Sk,2) -
		        16*l1k2*pow(Q2h,2)*pow(Sk,2) - 16*Q2e*pow(Q2h,2)*pow(Sk,2) +
		        8*pow(Q2h,3)*pow(Sk,2) + 32*l1k2*Q2e*Q2h*pow(Sq2,2) +
		        24*Q2h*pow(Q2e,2)*pow(Sq2,2) - 32*l1k2*pow(Q2h,2)*pow(Sq2,2) -
		        48*Q2e*pow(Q2h,2)*pow(Sq2,2) + 24*pow(Q2h,3)*pow(Sq2,2)) +
		     pow(2*l1k1 + Q2e - Q2k,-1)*
		      (64*k1k2*l1k2*m2*M2 + 32*k1k2*l1k2*M2*Q2e - 8*k1k2*m2*M2*Q2e -
		        16*m4*M2*Q2e + 32*l1k2*m2*M2*Q2h - 8*m4*Q2e*Q2h +
		        4*k1k2*M2*Q2e*Q2h + 12*l1k2*M2*Q2e*Q2h - 20*m2*M2*Q2e*Q2h +
		        8*k1k2*m2*M2*Q2k - 32*l1k2*m2*M2*Q2k + 16*m4*M2*Q2k +
		        8*k1k2*M2*Q2e*Q2k - 16*l1k2*M2*Q2e*Q2k + 12*m2*M2*Q2e*Q2k -
		        8*k1k2*m2*Q2h*Q2k + 40*m4*Q2h*Q2k - 8*k1k2*M2*Q2h*Q2k +
		        4*m2*M2*Q2h*Q2k + 4*k1k2*Q2e*Q2h*Q2k + 22*m2*Q2e*Q2h*Q2k -
		        2*M2*Q2e*Q2h*Q2k + 32*k1k2*m2*Q2h*S - 64*l1k2*m2*Q2h*S +
		        24*k1k2*Q2e*Q2h*S - 32*l1k2*Q2e*Q2h*S + 8*m2*Q2e*Q2h*S -
		        8*k1k2*Q2h*Q2k*S - 32*m2*Q2h*Q2k*S + 4*Q2e*Q2h*Q2k*S +
		        16*m4*Q2e*Sk - 16*k1k2*m2*Q2h*Sk + 64*m4*Q2h*Sk +
		        8*k1k2*Q2e*Q2h*Sk + 44*m2*Q2e*Q2h*Sk - 80*m4*Q2k*Sk -
		        32*m2*Q2e*Q2k*Sk - 16*m2*Q2h*Q2k*Sk - 4*Q2e*Q2h*Q2k*Sk -
		        32*m2*Q2h*S*Sk + 16*Q2e*Q2h*S*Sk + 32*m2*Q2k*S*Sk +
		        32*k1k2*m2*Q2h*Sq2 - 32*l1k2*m2*Q2h*Sq2 + 16*k1k2*Q2e*Q2h*Sq2 -
		        24*l1k2*Q2e*Q2h*Sq2 + 8*m2*Q2e*Q2h*Sq2 - 16*k1k2*Q2h*Q2k*Sq2 +
		        8*l1k2*Q2h*Q2k*Sq2 - 32*m2*Q2h*Q2k*Sq2 + 8*Q2e*Q2h*Q2k*Sq2 +
		        64*m2*Q2e*S*Sq2 + 64*m2*Q2h*S*Sq2 + 48*Q2e*Q2h*S*Sq2 -
		        64*m2*Q2k*S*Sq2 - 16*Q2h*Q2k*S*Sq2 - 32*m2*Q2h*Sk*Sq2 +
		        8*Q2e*Q2h*Sk*Sq2 + 64*m2*Q2k*Sk*Sq2 - 24*Q2e*Q2k*Sk*Sq2 +
		        8*Q2h*Q2k*Sk*Sq2 - 16*M2*Q2e*pow(k1k2,2) -
		        64*m2*M2*pow(l1k2,2) - 32*M2*Q2e*pow(l1k2,2) -
		        4*k1k2*M2*pow(Q2e,2) + 4*l1k2*M2*pow(Q2e,2) -
		        4*m2*M2*pow(Q2e,2) - 6*M2*Q2h*pow(Q2e,2) +
		        8*M2*Q2k*pow(Q2e,2) - 2*Q2h*Q2k*pow(Q2e,2) - 8*Q2k*S*pow(Q2e,2) -
		        4*Q2h*Sk*pow(Q2e,2) + 16*Q2k*Sk*pow(Q2e,2) - 16*S*Sk*pow(Q2e,2) -
		        4*Q2h*Sq2*pow(Q2e,2) - 4*Q2k*Sq2*pow(Q2e,2) - 8*Sk*Sq2*pow(Q2e,2) -
		        2*M2*pow(Q2e,3) - 16*m4*pow(Q2h,2) - 16*m2*Q2e*pow(Q2h,2) -
		        4*M2*Q2e*pow(Q2h,2) + 6*m2*Q2k*pow(Q2h,2) +
		        2*M2*Q2k*pow(Q2h,2) - 2*Q2e*Q2k*pow(Q2h,2) + 8*m2*S*pow(Q2h,2) +
		        8*Q2e*S*pow(Q2h,2) - 4*Q2k*S*pow(Q2h,2) - 4*m2*Sk*pow(Q2h,2) -
		        4*Q2e*Sk*pow(Q2h,2) + 4*Q2k*Sk*pow(Q2h,2) + 8*m2*Sq2*pow(Q2h,2) +
		        4*Q2e*Sq2*pow(Q2h,2) - 4*Q2k*Sq2*pow(Q2h,2) - 16*m4*pow(Q2k,2) -
		        8*m2*M2*pow(Q2k,2) - 10*m2*Q2e*pow(Q2k,2) -
		        8*M2*Q2e*pow(Q2k,2) - 10*m2*Q2h*pow(Q2k,2) +
		        4*M2*Q2h*pow(Q2k,2) - 2*Q2e*Q2h*pow(Q2k,2) + 16*m2*S*pow(Q2k,2) -
		        8*m2*Sk*pow(Q2k,2) + 32*m2*Sq2*pow(Q2k,2) -
		        12*Q2e*Sq2*pow(Q2k,2) + 4*Q2h*Sq2*pow(Q2k,2) +
		        4*pow(Q2e,2)*pow(Q2k,2) + 2*pow(Q2h,2)*pow(Q2k,2) +
		        64*m2*Q2e*pow(S,2) + 32*m2*Q2h*pow(S,2) + 32*Q2e*Q2h*pow(S,2) -
		        32*m2*Q2k*pow(S,2) - 16*Q2e*Q2k*pow(S,2) + 16*Q2h*Q2k*pow(S,2) -
		        64*m4*pow(Sk,2) - 24*m2*Q2e*pow(Sk,2) + 8*m2*Q2h*pow(Sk,2) -
		        16*m2*Q2k*pow(Sk,2) + 16*pow(Q2e,2)*pow(Sk,2) +
		        32*m2*Q2h*pow(Sq2,2) + 24*Q2e*Q2h*pow(Sq2,2) -
		        32*m2*Q2k*pow(Sq2,2) + 8*Q2e*Q2k*pow(Sq2,2) -
		        24*Q2h*Q2k*pow(Sq2,2) + 8*pow(Q2k,2)*pow(Sq2,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-1)*
		      (32*k1k2*m4*M2*Q2e + 16*k1k2*m4*Q2e*Q2h +
		        16*k1k2*m2*M2*Q2e*Q2h + 16*m4*M2*Q2e*Q2h -
		        32*k1k2*m4*M2*Q2k - 16*k1k2*m4*Q2h*Q2k -
		        48*k1k2*m2*M2*Q2h*Q2k - 16*m4*M2*Q2h*Q2k -
		        20*k1k2*m2*Q2e*Q2h*Q2k + 8*m4*Q2e*Q2h*Q2k -
		        16*k1k2*M2*Q2e*Q2h*Q2k + 8*m2*M2*Q2e*Q2h*Q2k -
		        48*k1k2*m2*Q2e*Q2h*S + 32*m4*Q2e*Q2h*S + 32*k1k2*m2*Q2e*Q2k*S -
		        32*m4*Q2e*Q2k*S + 64*k1k2*m2*Q2h*Q2k*S - 96*m4*Q2h*Q2k*S -
		        32*k1k2*Q2e*Q2h*Q2k*S - 32*m2*Q2e*Q2h*Q2k*S - 32*k1k2*m4*Q2e*Sk -
		        40*k1k2*m2*Q2e*Q2h*Sk + 16*m4*Q2e*Q2h*Sk + 32*k1k2*m4*Q2k*Sk +
		        32*k1k2*m2*Q2e*Q2k*Sk - 16*k1k2*m2*Q2h*Q2k*Sk -
		        32*m4*Q2h*Q2k*Sk - 8*k1k2*Q2e*Q2h*Q2k*Sk - 16*m2*Q2e*Q2h*Q2k*Sk +
		        64*k1k2*m2*Q2e*S*Sk - 64*m4*Q2e*S*Sk + 64*k1k2*m2*Q2h*S*Sk -
		        64*m4*Q2h*S*Sk - 16*k1k2*Q2e*Q2h*S*Sk - 32*m2*Q2e*Q2h*S*Sk -
		        128*k1k2*m2*Q2k*S*Sk + 128*m4*Q2k*S*Sk + 48*k1k2*Q2e*Q2k*S*Sk +
		        64*m2*Q2e*Q2k*S*Sk - 16*k1k2*Q2h*Q2k*S*Sk + 32*Q2e*Q2h*Q2k*S*Sk -
		        32*k1k2*m2*Q2e*Q2h*Sq2 + 16*m4*Q2e*Q2h*Sq2 +
		        32*k1k2*m2*Q2e*Q2k*Sq2 - 32*m4*Q2e*Q2k*Sq2 +
		        48*k1k2*m2*Q2h*Q2k*Sq2 - 80*m4*Q2h*Q2k*Sq2 -
		        32*k1k2*Q2e*Q2h*Q2k*Sq2 - 12*m2*Q2e*Q2h*Q2k*Sq2 -
		        64*m4*Q2e*S*Sq2 - 128*k1k2*m2*Q2h*S*Sq2 - 64*k1k2*Q2e*Q2h*S*Sq2 -
		        64*m2*Q2e*Q2h*S*Sq2 + 128*k1k2*m2*Q2k*S*Sq2 + 64*m4*Q2k*S*Sq2 -
		        32*k1k2*Q2e*Q2k*S*Sq2 - 64*m2*Q2e*Q2k*S*Sq2 +
		        128*k1k2*Q2h*Q2k*S*Sq2 + 96*m2*Q2h*Q2k*S*Sq2 -
		        80*Q2e*Q2h*Q2k*S*Sq2 + 64*k1k2*m2*Q2e*Sk*Sq2 - 64*m4*Q2e*Sk*Sq2 +
		        64*k1k2*m2*Q2h*Sk*Sq2 - 64*m4*Q2h*Sk*Sq2 -
		        16*k1k2*Q2e*Q2h*Sk*Sq2 - 128*k1k2*m2*Q2k*Sk*Sq2 +
		        128*m4*Q2k*Sk*Sq2 + 48*k1k2*Q2e*Q2k*Sk*Sq2 +
		        16*m2*Q2e*Q2k*Sk*Sq2 - 16*k1k2*Q2h*Q2k*Sk*Sq2 +
		        16*m2*Q2h*Q2k*Sk*Sq2 + 24*Q2e*Q2h*Q2k*Sk*Sq2 +
		        16*m2*M2*Q2e*pow(k1k2,2) - 16*M2*Q2e*Q2h*pow(k1k2,2) -
		        16*m2*M2*Q2k*pow(k1k2,2) + 8*M2*Q2e*Q2k*pow(k1k2,2) +
		        8*M2*Q2h*Q2k*pow(k1k2,2) - 16*Q2e*Q2h*S*pow(k1k2,2) +
		        16*Q2h*Q2k*S*pow(k1k2,2) - 16*Q2e*Q2h*Sq2*pow(k1k2,2) +
		        16*Q2h*Q2k*Sq2*pow(k1k2,2) + 8*k1k2*m2*M2*pow(Q2e,2) +
		        4*k1k2*M2*Q2h*pow(Q2e,2) - 4*m2*M2*Q2h*pow(Q2e,2) -
		        4*m2*M2*Q2k*pow(Q2e,2) - 8*m2*Q2h*Q2k*pow(Q2e,2) +
		        2*M2*Q2h*Q2k*pow(Q2e,2) + 8*k1k2*Q2h*S*pow(Q2e,2) -
		        8*k1k2*Q2k*S*pow(Q2e,2) - 16*m2*Q2k*S*pow(Q2e,2) -
		        8*Q2h*Q2k*S*pow(Q2e,2) - 16*m2*Q2h*Sk*pow(Q2e,2) +
		        24*m2*Q2k*Sk*pow(Q2e,2) + 4*Q2h*Q2k*Sk*pow(Q2e,2) -
		        16*k1k2*S*Sk*pow(Q2e,2) - 32*m2*S*Sk*pow(Q2e,2) -
		        32*Q2h*S*Sk*pow(Q2e,2) + 8*k1k2*Q2h*Sq2*pow(Q2e,2) -
		        8*k1k2*Q2k*Sq2*pow(Q2e,2) - 16*m2*Q2k*Sq2*pow(Q2e,2) -
		        4*Q2h*Q2k*Sq2*pow(Q2e,2) + 64*m2*S*Sq2*pow(Q2e,2) +
		        40*Q2h*S*Sq2*pow(Q2e,2) + 8*Q2k*S*Sq2*pow(Q2e,2) -
		        16*k1k2*Sk*Sq2*pow(Q2e,2) - 32*m2*Sk*Sq2*pow(Q2e,2) -
		        24*Q2h*Sk*Sq2*pow(Q2e,2) + 16*k1k2*m2*M2*pow(Q2h,2) +
		        16*k1k2*m2*Q2e*pow(Q2h,2) + 12*k1k2*M2*Q2e*pow(Q2h,2) +
		        16*m2*M2*Q2e*pow(Q2h,2) - 12*k1k2*m2*Q2k*pow(Q2h,2) +
		        16*m4*Q2k*pow(Q2h,2) - 8*k1k2*M2*Q2k*pow(Q2h,2) -
		        36*m2*M2*Q2k*pow(Q2h,2) + 4*k1k2*Q2e*Q2k*pow(Q2h,2) +
		        10*m2*Q2e*Q2k*pow(Q2h,2) - 2*M2*Q2e*Q2k*pow(Q2h,2) +
		        16*k1k2*m2*S*pow(Q2h,2) + 32*m4*S*pow(Q2h,2) +
		        16*k1k2*Q2e*S*pow(Q2h,2) + 16*m2*Q2e*S*pow(Q2h,2) -
		        8*k1k2*Q2k*S*pow(Q2h,2) + 8*Q2e*Q2k*S*pow(Q2h,2) +
		        8*k1k2*m2*Sk*pow(Q2h,2) + 16*m4*Sk*pow(Q2h,2) +
		        8*k1k2*Q2e*Sk*pow(Q2h,2) + 16*m2*Q2e*Sk*pow(Q2h,2) -
		        8*m2*Q2k*Sk*pow(Q2h,2) - 4*Q2e*Q2k*Sk*pow(Q2h,2) +
		        32*Q2e*S*Sk*pow(Q2h,2) - 32*Q2k*S*Sk*pow(Q2h,2) +
		        16*k1k2*m2*Sq2*pow(Q2h,2) + 32*m4*Sq2*pow(Q2h,2) +
		        16*k1k2*Q2e*Sq2*pow(Q2h,2) + 20*m2*Q2e*Sq2*pow(Q2h,2) -
		        16*k1k2*Q2k*Sq2*pow(Q2h,2) - 28*m2*Q2k*Sq2*pow(Q2h,2) +
		        4*Q2e*Q2k*Sq2*pow(Q2h,2) - 40*Q2e*S*Sq2*pow(Q2h,2) +
		        72*Q2k*S*Sq2*pow(Q2h,2) + 24*Q2e*Sk*Sq2*pow(Q2h,2) -
		        24*Q2k*Sk*Sq2*pow(Q2h,2) + 2*m2*pow(Q2e,2)*pow(Q2h,2) -
		        8*m4*pow(Q2h,3) + 4*m2*M2*pow(Q2h,3) - 2*m2*Q2e*pow(Q2h,3) -
		        2*m2*Q2k*pow(Q2h,3) + 4*m2*Sq2*pow(Q2h,3) +
		        8*k1k2*m2*M2*pow(Q2k,2) + 12*k1k2*m2*Q2e*pow(Q2k,2) +
		        4*k1k2*M2*Q2e*pow(Q2k,2) + 4*k1k2*m2*Q2h*pow(Q2k,2) -
		        16*m4*Q2h*pow(Q2k,2) + 4*k1k2*M2*Q2h*pow(Q2k,2) +
		        16*m2*M2*Q2h*pow(Q2k,2) - 4*k1k2*Q2e*Q2h*pow(Q2k,2) -
		        4*m2*Q2e*Q2h*pow(Q2k,2) + 4*M2*Q2e*Q2h*pow(Q2k,2) -
		        64*k1k2*m2*S*pow(Q2k,2) + 64*m4*S*pow(Q2k,2) +
		        24*k1k2*Q2e*S*pow(Q2k,2) + 32*m2*Q2e*S*pow(Q2k,2) +
		        8*Q2e*Q2h*S*pow(Q2k,2) + 16*k1k2*m2*Sk*pow(Q2k,2) -
		        4*Q2e*Q2h*Sk*pow(Q2k,2) - 64*k1k2*m2*Sq2*pow(Q2k,2) +
		        64*m4*Sq2*pow(Q2k,2) + 24*k1k2*Q2e*Sq2*pow(Q2k,2) +
		        8*m2*Q2e*Sq2*pow(Q2k,2) + 8*k1k2*Q2h*Sq2*pow(Q2k,2) +
		        24*m2*Q2h*Sq2*pow(Q2k,2) + 8*Q2e*Q2h*Sq2*pow(Q2k,2) -
		        32*k1k2*S*Sq2*pow(Q2k,2) - 32*m2*S*Sq2*pow(Q2k,2) +
		        8*Q2e*S*Sq2*pow(Q2k,2) - 8*Q2h*S*Sq2*pow(Q2k,2) +
		        6*m2*pow(Q2e,2)*pow(Q2k,2) - 2*M2*pow(Q2e,2)*pow(Q2k,2) -
		        2*m2*pow(Q2h,2)*pow(Q2k,2) - 2*M2*pow(Q2h,2)*pow(Q2k,2) -
		        8*S*pow(Q2h,2)*pow(Q2k,2) + 4*Sk*pow(Q2h,2)*pow(Q2k,2) -
		        8*Sq2*pow(Q2h,2)*pow(Q2k,2) + 2*M2*Q2e*pow(Q2k,3) -
		        2*M2*Q2h*pow(Q2k,3) - 64*m4*Q2e*pow(S,2) -
		        64*k1k2*m2*Q2h*pow(S,2) - 32*k1k2*Q2e*Q2h*pow(S,2) -
		        32*m2*Q2e*Q2h*pow(S,2) + 64*k1k2*m2*Q2k*pow(S,2) +
		        64*m4*Q2k*pow(S,2) - 32*k1k2*Q2e*Q2k*pow(S,2) -
		        64*m2*Q2e*Q2k*pow(S,2) + 64*k1k2*Q2h*Q2k*pow(S,2) +
		        64*m2*Q2h*Q2k*pow(S,2) - 64*Q2e*Q2h*Q2k*pow(S,2) +
		        32*m2*pow(Q2e,2)*pow(S,2) + 24*Q2h*pow(Q2e,2)*pow(S,2) +
		        8*Q2k*pow(Q2e,2)*pow(S,2) - 24*Q2e*pow(Q2h,2)*pow(S,2) +
		        56*Q2k*pow(Q2h,2)*pow(S,2) + 8*Q2e*pow(Q2k,2)*pow(S,2) -
		        8*Q2h*pow(Q2k,2)*pow(S,2) + 16*k1k2*m2*Q2e*pow(Sk,2) -
		        48*k1k2*m2*Q2h*pow(Sk,2) - 16*m2*Q2e*Q2h*pow(Sk,2) +
		        32*k1k2*m2*Q2k*pow(Sk,2) - 8*Q2e*Q2h*Q2k*pow(Sk,2) +
		        24*m2*pow(Q2e,2)*pow(Sk,2) + 8*Q2h*pow(Q2e,2)*pow(Sk,2) -
		        8*m2*pow(Q2h,2)*pow(Sk,2) - 8*Q2e*pow(Q2h,2)*pow(Sk,2) +
		        8*Q2k*pow(Q2h,2)*pow(Sk,2) - 64*k1k2*m2*Q2h*pow(Sq2,2) -
		        32*k1k2*Q2e*Q2h*pow(Sq2,2) - 48*m2*Q2e*Q2h*pow(Sq2,2) +
		        64*k1k2*m2*Q2k*pow(Sq2,2) + 64*k1k2*Q2h*Q2k*pow(Sq2,2) +
		        48*m2*Q2h*Q2k*pow(Sq2,2) - 24*Q2e*Q2h*Q2k*pow(Sq2,2) +
		        32*m2*pow(Q2e,2)*pow(Sq2,2) + 16*Q2h*pow(Q2e,2)*pow(Sq2,2) -
		        16*Q2e*pow(Q2h,2)*pow(Sq2,2) + 24*Q2k*pow(Q2h,2)*pow(Sq2,2) -
		        32*k1k2*pow(Q2k,2)*pow(Sq2,2) - 32*m2*pow(Q2k,2)*pow(Sq2,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      pow(2*l1k1 + Q2e - Q2k,-1)*
		      (-64*m4*Q2e*Q2h*S*Sq2 - 64*m2*Q2e*Q2h*Q2k*S*Sq2 -
		        8*m4*M2*Q2h*pow(Q2e,2) - 8*m2*M2*Q2h*Q2k*pow(Q2e,2) -
		        4*m2*Q2h*Q2k*Sk*pow(Q2e,2) - 16*Q2h*Q2k*S*Sk*pow(Q2e,2) +
		        8*m4*Q2h*Sq2*pow(Q2e,2) + 6*m2*Q2h*Q2k*Sq2*pow(Q2e,2) +
		        32*m4*S*Sq2*pow(Q2e,2) + 16*m2*Q2h*S*Sq2*pow(Q2e,2) +
		        16*m2*Q2k*S*Sq2*pow(Q2e,2) + 44*Q2h*Q2k*S*Sq2*pow(Q2e,2) +
		        8*m2*Q2h*Sk*Sq2*pow(Q2e,2) - 12*Q2h*Q2k*Sk*Sq2*pow(Q2e,2) -
		        M2*Q2h*Q2k*pow(Q2e,3) + 4*Q2h*Q2k*S*pow(Q2e,3) -
		        2*Q2h*Q2k*Sk*pow(Q2e,3) + 16*Q2h*S*Sk*pow(Q2e,3) +
		        2*Q2h*Q2k*Sq2*pow(Q2e,3) - 20*Q2h*S*Sq2*pow(Q2e,3) -
		        4*Q2k*S*Sq2*pow(Q2e,3) + 12*Q2h*Sk*Sq2*pow(Q2e,3) +
		        16*m4*M2*Q2e*pow(Q2h,2) + 24*m2*M2*Q2e*Q2k*pow(Q2h,2) +
		        8*m2*Q2e*Q2k*Sk*pow(Q2h,2) + 32*Q2e*Q2k*S*Sk*pow(Q2h,2) -
		        16*m4*Q2e*Sq2*pow(Q2h,2) + 4*m2*Q2e*Q2k*Sq2*pow(Q2h,2) +
		        32*m4*S*Sq2*pow(Q2h,2) - 16*m2*Q2e*S*Sq2*pow(Q2h,2) +
		        48*m2*Q2k*S*Sq2*pow(Q2h,2) - 76*Q2e*Q2k*S*Sq2*pow(Q2h,2) -
		        16*m2*Q2e*Sk*Sq2*pow(Q2h,2) + 24*Q2e*Q2k*Sk*Sq2*pow(Q2h,2) -
		        4*m4*pow(Q2e,2)*pow(Q2h,2) - 2*m2*M2*pow(Q2e,2)*pow(Q2h,2) -
		        3*m2*Q2k*pow(Q2e,2)*pow(Q2h,2) + 2*M2*Q2k*pow(Q2e,2)*pow(Q2h,2) -
		        8*Q2k*S*pow(Q2e,2)*pow(Q2h,2) + 4*Q2k*Sk*pow(Q2e,2)*pow(Q2h,2) -
		        32*S*Sk*pow(Q2e,2)*pow(Q2h,2) - 6*m2*Sq2*pow(Q2e,2)*pow(Q2h,2) -
		        4*Q2k*Sq2*pow(Q2e,2)*pow(Q2h,2) + 40*S*Sq2*pow(Q2e,2)*pow(Q2h,2) -
		        24*Sk*Sq2*pow(Q2e,2)*pow(Q2h,2) - 8*m4*M2*pow(Q2h,3) +
		        8*m4*Q2e*pow(Q2h,3) - 16*m2*M2*Q2k*pow(Q2h,3) +
		        4*m2*Q2e*Q2k*pow(Q2h,3) - M2*Q2e*Q2k*pow(Q2h,3) +
		        4*Q2e*Q2k*S*pow(Q2h,3) - 4*m2*Q2k*Sk*pow(Q2h,3) -
		        2*Q2e*Q2k*Sk*pow(Q2h,3) + 16*Q2e*S*Sk*pow(Q2h,3) -
		        16*Q2k*S*Sk*pow(Q2h,3) + 8*m4*Sq2*pow(Q2h,3) +
		        4*m2*Q2e*Sq2*pow(Q2h,3) - 10*m2*Q2k*Sq2*pow(Q2h,3) +
		        2*Q2e*Q2k*Sq2*pow(Q2h,3) - 20*Q2e*S*Sq2*pow(Q2h,3) +
		        36*Q2k*S*Sq2*pow(Q2h,3) + 8*m2*Sk*Sq2*pow(Q2h,3) +
		        12*Q2e*Sk*Sq2*pow(Q2h,3) - 12*Q2k*Sk*Sq2*pow(Q2h,3) +
		        m2*pow(Q2e,2)*pow(Q2h,3) - 4*m4*pow(Q2h,4) +
		        2*m2*M2*pow(Q2h,4) - m2*Q2e*pow(Q2h,4) - m2*Q2k*pow(Q2h,4) +
		        2*m2*Sq2*pow(Q2h,4) - 8*m2*M2*Q2e*Q2h*pow(Q2k,2) -
		        8*m2*Q2e*Q2h*Sq2*pow(Q2k,2) + 16*m2*Q2e*S*Sq2*pow(Q2k,2) -
		        16*m2*Q2h*S*Sq2*pow(Q2k,2) + 8*Q2e*Q2h*S*Sq2*pow(Q2k,2) +
		        2*m2*M2*pow(Q2e,2)*pow(Q2k,2) -
		        3*M2*Q2h*pow(Q2e,2)*pow(Q2k,2) - 4*Q2h*S*pow(Q2e,2)*pow(Q2k,2) +
		        2*Q2h*Sk*pow(Q2e,2)*pow(Q2k,2) - 4*Q2h*Sq2*pow(Q2e,2)*pow(Q2k,2) -
		        4*S*Sq2*pow(Q2e,2)*pow(Q2k,2) + M2*pow(Q2e,3)*pow(Q2k,2) +
		        6*m2*M2*pow(Q2h,2)*pow(Q2k,2) + m2*Q2e*pow(Q2h,2)*pow(Q2k,2) +
		        3*M2*Q2e*pow(Q2h,2)*pow(Q2k,2) + 8*Q2e*S*pow(Q2h,2)*pow(Q2k,2) -
		        4*Q2e*Sk*pow(Q2h,2)*pow(Q2k,2) + 8*m2*Sq2*pow(Q2h,2)*pow(Q2k,2) +
		        8*Q2e*Sq2*pow(Q2h,2)*pow(Q2k,2) - 4*S*Sq2*pow(Q2h,2)*pow(Q2k,2) -
		        m2*pow(Q2h,3)*pow(Q2k,2) - M2*pow(Q2h,3)*pow(Q2k,2) -
		        4*S*pow(Q2h,3)*pow(Q2k,2) + 2*Sk*pow(Q2h,3)*pow(Q2k,2) -
		        4*Sq2*pow(Q2h,3)*pow(Q2k,2) + 2*M2*Q2e*Q2h*pow(Q2k,3) -
		        M2*pow(Q2e,2)*pow(Q2k,3) - M2*pow(Q2h,2)*pow(Q2k,3) -
		        64*m4*Q2e*Q2h*pow(S,2) - 64*m2*Q2e*Q2h*Q2k*pow(S,2) +
		        32*m4*pow(Q2e,2)*pow(S,2) + 16*m2*Q2h*pow(Q2e,2)*pow(S,2) +
		        16*m2*Q2k*pow(Q2e,2)*pow(S,2) + 36*Q2h*Q2k*pow(Q2e,2)*pow(S,2) -
		        12*Q2h*pow(Q2e,3)*pow(S,2) - 4*Q2k*pow(Q2e,3)*pow(S,2) +
		        32*m4*pow(Q2h,2)*pow(S,2) - 16*m2*Q2e*pow(Q2h,2)*pow(S,2) +
		        48*m2*Q2k*pow(Q2h,2)*pow(S,2) - 60*Q2e*Q2k*pow(Q2h,2)*pow(S,2) +
		        24*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) - 12*Q2e*pow(Q2h,3)*pow(S,2) +
		        28*Q2k*pow(Q2h,3)*pow(S,2) + 16*m2*Q2e*pow(Q2k,2)*pow(S,2) -
		        16*m2*Q2h*pow(Q2k,2)*pow(S,2) + 8*Q2e*Q2h*pow(Q2k,2)*pow(S,2) -
		        4*pow(Q2e,2)*pow(Q2k,2)*pow(S,2) -
		        4*pow(Q2h,2)*pow(Q2k,2)*pow(S,2) - 8*m2*Q2h*pow(Q2e,2)*pow(Sk,2) +
		        4*Q2h*Q2k*pow(Q2e,2)*pow(Sk,2) - 4*Q2h*pow(Q2e,3)*pow(Sk,2) +
		        16*m2*Q2e*pow(Q2h,2)*pow(Sk,2) - 8*Q2e*Q2k*pow(Q2h,2)*pow(Sk,2) +
		        8*pow(Q2e,2)*pow(Q2h,2)*pow(Sk,2) - 8*m2*pow(Q2h,3)*pow(Sk,2) -
		        4*Q2e*pow(Q2h,3)*pow(Sk,2) + 4*Q2k*pow(Q2h,3)*pow(Sk,2) -
		        24*m2*Q2e*Q2h*Q2k*pow(Sq2,2) + 8*m2*Q2h*pow(Q2e,2)*pow(Sq2,2) +
		        12*Q2h*Q2k*pow(Q2e,2)*pow(Sq2,2) - 8*Q2h*pow(Q2e,3)*pow(Sq2,2) -
		        4*m2*Q2e*pow(Q2h,2)*pow(Sq2,2) + 24*m2*Q2k*pow(Q2h,2)*pow(Sq2,2) -
		        24*Q2e*Q2k*pow(Q2h,2)*pow(Sq2,2) +
		        16*pow(Q2e,2)*pow(Q2h,2)*pow(Sq2,2) - 4*m2*pow(Q2h,3)*pow(Sq2,2) -
		        8*Q2e*pow(Q2h,3)*pow(Sq2,2) + 12*Q2k*pow(Q2h,3)*pow(Sq2,2) +
		        12*m2*Q2e*pow(Q2k,2)*pow(Sq2,2) - 12*m2*Q2h*pow(Q2k,2)*pow(Sq2,2)) \
		+ pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*pow(2*l1k1 + Q2e - Q2k,-1)*
		      (16*k1k2*m2*M2*Q2e*Q2h*Q2k - 16*m4*M2*Q2e*Q2h*Q2k +
		        32*k1k2*m2*Q2e*Q2h*Q2k*S - 32*m4*Q2e*Q2h*Q2k*S +
		        32*k1k2*m2*Q2e*Q2h*Q2k*Sq2 - 32*m4*Q2e*Q2h*Q2k*Sq2 +
		        128*k1k2*m2*Q2e*Q2h*S*Sq2 - 128*m4*Q2e*Q2h*S*Sq2 -
		        128*k1k2*m2*Q2e*Q2k*S*Sq2 + 128*m4*Q2e*Q2k*S*Sq2 -
		        128*k1k2*m2*Q2h*Q2k*S*Sq2 + 128*m4*Q2h*Q2k*S*Sq2 +
		        8*M2*Q2e*Q2h*Q2k*pow(k1k2,2) + 16*Q2e*Q2h*Q2k*S*pow(k1k2,2) +
		        16*Q2e*Q2h*Q2k*Sq2*pow(k1k2,2) + 64*Q2e*Q2h*S*Sq2*pow(k1k2,2) -
		        64*Q2e*Q2k*S*Sq2*pow(k1k2,2) - 64*Q2h*Q2k*S*Sq2*pow(k1k2,2) -
		        16*k1k2*m2*M2*Q2e*pow(Q2h,2) + 16*m4*M2*Q2e*pow(Q2h,2) +
		        16*k1k2*m2*M2*Q2k*pow(Q2h,2) - 16*m4*M2*Q2k*pow(Q2h,2) -
		        32*k1k2*m2*Q2e*S*pow(Q2h,2) + 32*m4*Q2e*S*pow(Q2h,2) +
		        32*k1k2*m2*Q2k*S*pow(Q2h,2) - 32*m4*Q2k*S*pow(Q2h,2) -
		        32*k1k2*m2*Q2e*Sq2*pow(Q2h,2) + 32*m4*Q2e*Sq2*pow(Q2h,2) +
		        32*k1k2*m2*Q2k*Sq2*pow(Q2h,2) - 32*m4*Q2k*Sq2*pow(Q2h,2) -
		        8*M2*Q2e*pow(k1k2,2)*pow(Q2h,2) +
		        8*M2*Q2k*pow(k1k2,2)*pow(Q2h,2) - 16*Q2e*S*pow(k1k2,2)*pow(Q2h,2) +
		        16*Q2k*S*pow(k1k2,2)*pow(Q2h,2) - 16*Q2e*Sq2*pow(k1k2,2)*pow(Q2h,2) +
		        16*Q2k*Sq2*pow(k1k2,2)*pow(Q2h,2) - 16*k1k2*m2*M2*Q2h*pow(Q2k,2) +
		        16*m4*M2*Q2h*pow(Q2k,2) - 32*k1k2*m2*Q2h*S*pow(Q2k,2) +
		        32*m4*Q2h*S*pow(Q2k,2) - 32*k1k2*m2*Q2h*Sq2*pow(Q2k,2) +
		        32*m4*Q2h*Sq2*pow(Q2k,2) + 128*k1k2*m2*S*Sq2*pow(Q2k,2) -
		        128*m4*S*Sq2*pow(Q2k,2) - 8*M2*Q2h*pow(k1k2,2)*pow(Q2k,2) -
		        16*Q2h*S*pow(k1k2,2)*pow(Q2k,2) - 16*Q2h*Sq2*pow(k1k2,2)*pow(Q2k,2) +
		        64*S*Sq2*pow(k1k2,2)*pow(Q2k,2) + 64*k1k2*m2*Q2e*Q2h*pow(S,2) -
		        64*m4*Q2e*Q2h*pow(S,2) - 64*k1k2*m2*Q2e*Q2k*pow(S,2) +
		        64*m4*Q2e*Q2k*pow(S,2) - 64*k1k2*m2*Q2h*Q2k*pow(S,2) +
		        64*m4*Q2h*Q2k*pow(S,2) + 32*Q2e*Q2h*pow(k1k2,2)*pow(S,2) -
		        32*Q2e*Q2k*pow(k1k2,2)*pow(S,2) - 32*Q2h*Q2k*pow(k1k2,2)*pow(S,2) +
		        64*k1k2*m2*pow(Q2k,2)*pow(S,2) - 64*m4*pow(Q2k,2)*pow(S,2) +
		        32*pow(k1k2,2)*pow(Q2k,2)*pow(S,2) + 64*k1k2*m2*Q2e*Q2h*pow(Sq2,2) -
		        64*m4*Q2e*Q2h*pow(Sq2,2) - 64*k1k2*m2*Q2e*Q2k*pow(Sq2,2) +
		        64*m4*Q2e*Q2k*pow(Sq2,2) - 64*k1k2*m2*Q2h*Q2k*pow(Sq2,2) +
		        64*m4*Q2h*Q2k*pow(Sq2,2) + 32*Q2e*Q2h*pow(k1k2,2)*pow(Sq2,2) -
		        32*Q2e*Q2k*pow(k1k2,2)*pow(Sq2,2) -
		        32*Q2h*Q2k*pow(k1k2,2)*pow(Sq2,2) +
		        64*k1k2*m2*pow(Q2k,2)*pow(Sq2,2) - 64*m4*pow(Q2k,2)*pow(Sq2,2) +
		        32*pow(k1k2,2)*pow(Q2k,2)*pow(Sq2,2)));
}

long double Melem::melem2_l2k2_1(const long double l1k1,const long double l1k2,const long double k1k2,const long double Q2e,
		const long double Q2h,const long double Q2k,const long double S,const long double Sk,const long double Sq2,
		const long double f1, const long double f2)const {

	long double p1p2 = M2 + Q2h/2.;
	long double G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
	long double G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
	long double G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;

	return G1*(-48*l1k1*m2 + 24*l1k2*m2 - 16*m4 + 4*k1k2*Q2e - 8*l1k1*Q2e -
		     4*m2*Q2e - 4*l1k2*Q2h + 12*m2*Q2h + 2*Q2e*Q2h + 8*m2*Q2k +
		     2*Q2e*Q2k + pow(k1k2 - l1k1 - l1k2,-2)*
		      (-16*l1k1*l1k2*m4 + 16*l1k1*m6 - 8*l1k1*m4*Q2e - 8*l1k2*m4*Q2e +
		        8*m6*Q2e + 8*l1k1*l1k2*m2*Q2h - 8*l1k1*m4*Q2h +
		        4*l1k1*m2*Q2e*Q2h + 4*l1k2*m2*Q2e*Q2h - 4*m4*Q2e*Q2h +
		        8*l1k1*m4*Q2k + 8*l1k2*m4*Q2k - 8*m6*Q2k - 4*l1k1*m2*Q2h*Q2k -
		        4*l1k2*m2*Q2h*Q2k + 4*m4*Q2h*Q2k - 16*m4*pow(l1k1,2) -
		        4*m2*Q2e*pow(l1k1,2) + 8*m2*Q2h*pow(l1k1,2) +
		        2*Q2e*Q2h*pow(l1k1,2) + 4*m2*Q2k*pow(l1k1,2) -
		        2*Q2h*Q2k*pow(l1k1,2) - 8*m2*pow(l1k1,3) + 4*Q2h*pow(l1k1,3) -
		        8*l1k1*m2*pow(l1k2,2) - 4*m2*Q2e*pow(l1k2,2) +
		        4*l1k1*Q2h*pow(l1k2,2) + 2*Q2e*Q2h*pow(l1k2,2) +
		        4*m2*Q2k*pow(l1k2,2) - 2*Q2h*Q2k*pow(l1k2,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (16*k1k2*l1k2*m2 + 64*l1k2*m4 + 8*k1k2*l1k2*Q2e - 48*l1k2*m2*Q2e -
		        8*k1k2*l1k2*Q2h + 24*l1k2*m2*Q2h + 4*l1k2*Q2e*Q2h -
		        24*l1k2*m2*Q2k - 8*l1k2*Q2e*Q2k + 4*m2*Q2e*Q2k + 8*l1k2*Q2h*Q2k -
		        4*m2*Q2h*Q2k - 128*m2*pow(l1k2,2) - 16*Q2e*pow(l1k2,2) +
		        16*Q2h*pow(l1k2,2) - 4*l1k2*pow(Q2e,2) - 4*m2*pow(Q2e,2) -
		        8*l1k2*pow(Q2h,2) + 4*m2*pow(Q2h,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		      (-16*k1k2*l1k2*m2*Q2e - 32*l1k2*m4*Q2e + 16*k1k2*l1k2*m2*Q2h +
		        32*l1k2*m4*Q2h + 8*k1k2*l1k2*Q2e*Q2h + 8*l1k2*m2*Q2e*Q2h +
		        8*l1k2*m2*Q2e*Q2k - 8*l1k2*m2*Q2h*Q2k - 4*l1k2*Q2e*Q2h*Q2k -
		        32*k1k2*m2*pow(l1k2,2) - 64*m4*pow(l1k2,2) +
		        16*m2*Q2e*pow(l1k2,2) + 16*k1k2*Q2h*pow(l1k2,2) -
		        8*Q2e*Q2h*pow(l1k2,2) + 16*m2*Q2k*pow(l1k2,2) -
		        8*Q2h*Q2k*pow(l1k2,2) + 32*m2*pow(l1k2,3) - 16*Q2h*pow(l1k2,3) -
		        8*k1k2*l1k2*pow(Q2h,2) - 8*l1k2*m2*pow(Q2h,2) +
		        4*l1k2*Q2e*pow(Q2h,2) + 4*l1k2*Q2k*pow(Q2h,2) +
		        16*pow(l1k2,2)*pow(Q2h,2) - 4*l1k2*pow(Q2h,3)) +
		     pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
		      (-64*l1k2*m4 + 32*l1k1*m2*Q2e + 48*l1k2*m2*Q2e - 32*m4*Q2e +
		        8*l1k2*m2*Q2h - 8*l1k2*Q2e*Q2h - 32*l1k1*m2*Q2k +
		        8*l1k2*m2*Q2k - 16*l1k1*Q2e*Q2k + 8*l1k2*Q2e*Q2k +
		        8*l1k2*Q2h*Q2k - 8*Q2e*Q2h*Q2k + 64*m2*pow(l1k1,2) +
		        32*Q2e*pow(l1k1,2) + 96*m2*pow(l1k2,2) + 24*Q2e*pow(l1k2,2) +
		        8*Q2h*pow(l1k2,2) + 16*l1k1*pow(Q2e,2) + 8*l1k2*pow(Q2e,2) -
		        8*m2*pow(Q2e,2) + 4*pow(Q2e,3) - 16*l1k2*pow(Q2h,2) +
		        8*m2*pow(Q2h,2) + 8*Q2e*pow(Q2h,2) + 4*Q2e*pow(Q2k,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*
		      (-24*l1k1*l1k2*m2 + 16*l1k2*m4 + 12*l1k1*l1k2*Q2e -
		        32*l1k1*m2*Q2e - 24*l1k2*m2*Q2e + 8*m4*Q2e - 28*l1k1*l1k2*Q2h -
		        8*l1k1*m2*Q2h + 4*l1k2*m2*Q2h + 2*l1k1*Q2e*Q2h + 6*l1k2*Q2e*Q2h +
		        2*m2*Q2e*Q2h + 32*l1k1*m2*Q2k + 4*l1k2*m2*Q2k - 8*m4*Q2k +
		        4*l1k1*Q2e*Q2k - 8*l1k2*Q2e*Q2k + 4*m2*Q2e*Q2k - 4*l1k1*Q2h*Q2k +
		        16*l1k2*Q2h*Q2k + 4*Q2e*Q2h*Q2k - 64*m2*pow(l1k1,2) -
		        8*Q2e*pow(l1k1,2) + 4*Q2h*pow(l1k1,2) - 32*m2*pow(l1k2,2) -
		        12*Q2e*pow(l1k2,2) + 32*Q2h*pow(l1k2,2) - 4*l1k1*pow(Q2e,2) -
		        2*l1k2*pow(Q2e,2) + Q2h*pow(Q2e,2) - pow(Q2e,3) +
		        6*l1k1*pow(Q2h,2) - 16*l1k2*pow(Q2h,2) - 6*m2*pow(Q2h,2) -
		        4*Q2e*pow(Q2h,2) - 4*Q2k*pow(Q2h,2) + 4*pow(Q2h,3) -
		        Q2e*pow(Q2k,2) + Q2h*pow(Q2k,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (-32*l1k2*m4*Q2e + 32*l1k2*m4*Q2h - 16*l1k2*m2*Q2e*Q2h +
		        16*l1k2*m2*Q2e*Q2k - 16*l1k2*m2*Q2h*Q2k - 32*l1k2*Q2e*Q2h*Q2k +
		        32*m2*Q2e*pow(l1k2,2) - 64*m2*Q2h*pow(l1k2,2) -
		        64*Q2e*Q2h*pow(l1k2,2) + 32*m2*Q2k*pow(l1k2,2) +
		        16*Q2e*Q2k*pow(l1k2,2) - 32*Q2h*Q2k*pow(l1k2,2) +
		        64*m2*pow(l1k2,3) + 32*Q2e*pow(l1k2,3) - 64*Q2h*pow(l1k2,3) -
		        8*l1k2*m2*pow(Q2e,2) - 12*l1k2*Q2h*pow(Q2e,2) +
		        8*l1k2*Q2k*pow(Q2e,2) + 16*pow(l1k2,2)*pow(Q2e,2) +
		        4*l1k2*pow(Q2e,3) + 24*l1k2*m2*pow(Q2h,2) +
		        40*l1k2*Q2e*pow(Q2h,2) + 24*l1k2*Q2k*pow(Q2h,2) +
		        64*pow(l1k2,2)*pow(Q2h,2) - 32*l1k2*pow(Q2h,3) +
		        4*l1k2*Q2e*pow(Q2k,2) - 4*l1k2*Q2h*pow(Q2k,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
		      (8*m4*Q2e*Q2h + 4*l1k2*m2*Q2e*Q2k - 4*l1k2*m2*Q2h*Q2k +
		        2*m2*Q2e*Q2h*Q2k - 4*l1k2*m2*pow(Q2e,2) - 2*m2*Q2h*pow(Q2e,2) -
		        Q2h*pow(Q2e,3) + 4*l1k2*m2*pow(Q2h,2) - 8*m4*pow(Q2h,2) +
		        2*m2*Q2e*pow(Q2h,2) - 2*m2*Q2k*pow(Q2h,2) +
		        4*Q2e*Q2k*pow(Q2h,2) + pow(Q2e,2)*pow(Q2h,2) - 4*Q2e*pow(Q2h,3) -
		        4*Q2k*pow(Q2h,3) + 4*pow(Q2h,4) - Q2e*Q2h*pow(Q2k,2) +
		        pow(Q2h,2)*pow(Q2k,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
		      (-64*l1k2*m6 + 64*l1k2*m4*Q2e - 32*l1k2*m4*Q2h +
		        40*l1k2*m2*Q2e*Q2h + 16*m4*Q2e*Q2h - 32*l1k2*m2*Q2e*Q2k +
		        20*l1k2*Q2e*Q2h*Q2k + 4*m2*Q2e*Q2h*Q2k + 128*m4*pow(l1k2,2) -
		        96*m2*Q2e*pow(l1k2,2) + 64*m2*Q2h*pow(l1k2,2) +
		        24*Q2e*Q2h*pow(l1k2,2) - 64*m2*Q2k*pow(l1k2,2) -
		        16*Q2e*Q2k*pow(l1k2,2) - 32*Q2h*Q2k*pow(l1k2,2) -
		        160*m2*pow(l1k2,3) - 32*Q2e*pow(l1k2,3) - 32*Q2h*pow(l1k2,3) +
		        12*l1k2*Q2h*pow(Q2e,2) - 4*m2*Q2h*pow(Q2e,2) -
		        8*l1k2*Q2k*pow(Q2e,2) - 16*pow(l1k2,2)*pow(Q2e,2) -
		        4*l1k2*pow(Q2e,3) - 2*Q2h*pow(Q2e,3) - 16*l1k2*m2*pow(Q2h,2) -
		        16*m4*pow(Q2h,2) - 20*l1k2*Q2e*pow(Q2h,2) + 4*m2*Q2e*pow(Q2h,2) +
		        12*l1k2*Q2k*pow(Q2h,2) - 4*m2*Q2k*pow(Q2h,2) +
		        8*Q2e*Q2k*pow(Q2h,2) + 56*pow(l1k2,2)*pow(Q2h,2) +
		        2*pow(Q2e,2)*pow(Q2h,2) - 4*l1k2*pow(Q2h,3) - 8*Q2e*pow(Q2h,3) -
		        8*Q2k*pow(Q2h,3) + 8*pow(Q2h,4) + 8*l1k2*m2*pow(Q2k,2) -
		        4*l1k2*Q2e*pow(Q2k,2) - 4*l1k2*Q2h*pow(Q2k,2) -
		        2*Q2e*Q2h*pow(Q2k,2) + 2*pow(Q2h,2)*pow(Q2k,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		      pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
		      (64*l1k2*m6*Q2e - 64*l1k2*m6*Q2h - 32*l1k2*m4*Q2e*Q2k +
		        32*l1k2*m4*Q2h*Q2k + 32*l1k2*m2*Q2e*Q2h*Q2k +
		        128*m6*pow(l1k2,2) - 64*m4*Q2e*pow(l1k2,2) +
		        64*m4*Q2h*pow(l1k2,2) + 64*m2*Q2e*Q2h*pow(l1k2,2) -
		        64*m4*Q2k*pow(l1k2,2) - 32*m2*Q2e*Q2k*pow(l1k2,2) +
		        96*m2*Q2h*Q2k*pow(l1k2,2) + 16*Q2e*Q2h*Q2k*pow(l1k2,2) -
		        128*m4*pow(l1k2,3) - 32*m2*Q2e*pow(l1k2,3) +
		        160*m2*Q2h*pow(l1k2,3) + 16*Q2e*Q2h*pow(l1k2,3) -
		        64*m2*Q2k*pow(l1k2,3) + 32*Q2h*Q2k*pow(l1k2,3) -
		        64*m2*pow(l1k2,4) + 32*Q2h*pow(l1k2,4) -
		        24*l1k2*m2*Q2e*pow(Q2h,2) - 32*l1k2*m2*Q2k*pow(Q2h,2) -
		        8*l1k2*Q2e*Q2k*pow(Q2h,2) - 112*m2*pow(l1k2,2)*pow(Q2h,2) -
		        16*Q2e*pow(l1k2,2)*pow(Q2h,2) - 32*Q2k*pow(l1k2,2)*pow(Q2h,2) -
		        48*pow(l1k2,3)*pow(Q2h,2) + 24*l1k2*m2*pow(Q2h,3) +
		        4*l1k2*Q2e*pow(Q2h,3) + 8*l1k2*Q2k*pow(Q2h,3) +
		        24*pow(l1k2,2)*pow(Q2h,3) - 4*l1k2*pow(Q2h,4) -
		        8*l1k2*m2*Q2e*pow(Q2k,2) + 8*l1k2*m2*Q2h*pow(Q2k,2) +
		        4*l1k2*Q2e*Q2h*pow(Q2k,2) - 16*m2*pow(l1k2,2)*pow(Q2k,2) +
		        8*Q2h*pow(l1k2,2)*pow(Q2k,2) - 4*l1k2*pow(Q2h,2)*pow(Q2k,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
		      (32*m6*Q2e*Q2h - 16*m6*pow(Q2e,2) + 8*m4*Q2h*pow(Q2e,2) +
		        8*m4*Q2k*pow(Q2e,2) - 6*m2*Q2h*Q2k*pow(Q2e,2) - 8*m4*pow(Q2e,3) +
		        4*m2*Q2h*pow(Q2e,3) + 2*m2*Q2k*pow(Q2e,3) + 2*Q2h*Q2k*pow(Q2e,3) +
		        2*Q2h*pow(Q2e,4) - Q2k*pow(Q2e,4) - 16*m6*pow(Q2h,2) -
		        8*m4*Q2k*pow(Q2h,2) - 2*m2*Q2e*Q2k*pow(Q2h,2) -
		        10*m2*pow(Q2e,2)*pow(Q2h,2) - 13*Q2k*pow(Q2e,2)*pow(Q2h,2) -
		        4*pow(Q2e,3)*pow(Q2h,2) + 12*m2*Q2e*pow(Q2h,3) +
		        6*m2*Q2k*pow(Q2h,3) + 24*Q2e*Q2k*pow(Q2h,3) +
		        10*pow(Q2e,2)*pow(Q2h,3) - 6*m2*pow(Q2h,4) - 16*Q2e*pow(Q2h,4) -
		        12*Q2k*pow(Q2h,4) + 8*pow(Q2h,5) - 8*m4*Q2e*pow(Q2k,2) +
		        8*m4*Q2h*pow(Q2k,2) + 2*m2*pow(Q2e,2)*pow(Q2k,2) +
		        6*Q2h*pow(Q2e,2)*pow(Q2k,2) - 2*m2*pow(Q2h,2)*pow(Q2k,2) -
		        12*Q2e*pow(Q2h,2)*pow(Q2k,2) + 6*pow(Q2h,3)*pow(Q2k,2) +
		        2*Q2e*Q2h*pow(Q2k,3) - pow(Q2e,2)*pow(Q2k,3) - pow(Q2h,2)*pow(Q2k,3))\
		) + (G2 + G3)*(8*k1k2*m2*M2 - 16*l1k1*m2*M2 - 16*l1k2*m2*M2 +
		     16*m4*M2 + 4*k1k2*M2*Q2e - 12*l1k1*M2*Q2e + 4*l1k2*M2*Q2e +
		     4*m2*M2*Q2e + 8*m4*Q2h + 4*l1k1*M2*Q2h - 8*l1k2*M2*Q2h +
		     16*m2*M2*Q2h + 4*M2*Q2e*Q2h + 4*m2*M2*Q2k + 2*M2*Q2e*Q2k -
		     4*m2*Q2h*Q2k - 8*l1k2*Q2h*S - 16*m2*Q2h*S + 32*m2*Q2k*S -
		     12*Q2e*Q2k*S + 4*Q2h*Q2k*S - 16*m4*Sk - 8*m2*Q2h*Sk - 8*m2*Q2k*Sk +
		     64*m2*S*Sk - 24*Q2e*S*Sk + 8*Q2h*S*Sk - 8*m2*Q2h*Sq2 + 4*Q2e*Q2h*Sq2 +
		     32*m2*Q2k*Sq2 - 12*Q2e*Q2k*Sq2 + 4*Q2h*Q2k*Sq2 - 32*k1k2*S*Sq2 -
		     32*l1k1*S*Sq2 + 64*l1k2*S*Sq2 - 128*m2*S*Sq2 - 32*Q2h*S*Sq2 +
		     16*Q2k*S*Sq2 + 64*m2*Sk*Sq2 - 24*Q2e*Sk*Sq2 + 8*Q2h*Sk*Sq2 -
		     4*M2*pow(Q2e,2) + 8*m2*pow(Q2h,2) + 2*M2*pow(Q2h,2) -
		     16*k1k2*pow(S,2) - 32*l1k1*pow(S,2) + 48*l1k2*pow(S,2) -
		     64*m2*pow(S,2) - 24*Q2h*pow(S,2) + 8*Q2k*pow(S,2) +
		     pow(k1k2 - l1k1 - l1k2,-2)*
		      (8*l1k1*l1k2*m2*M2*Q2h - 8*l1k1*m4*M2*Q2h +
		        4*l1k1*m2*M2*Q2e*Q2h + 4*l1k2*m2*M2*Q2e*Q2h -
		        4*m4*M2*Q2e*Q2h - 4*l1k1*m2*M2*Q2h*Q2k -
		        4*l1k2*m2*M2*Q2h*Q2k + 4*m4*M2*Q2h*Q2k -
		        16*l1k1*l1k2*m2*Q2h*S + 16*l1k1*m4*Q2h*S - 8*l1k1*m2*Q2e*Q2h*S -
		        8*l1k2*m2*Q2e*Q2h*S + 8*m4*Q2e*Q2h*S + 8*l1k1*m2*Q2h*Q2k*S +
		        8*l1k2*m2*Q2h*Q2k*S - 8*m4*Q2h*Q2k*S +
		        8*m2*M2*Q2h*pow(l1k1,2) + 2*M2*Q2e*Q2h*pow(l1k1,2) -
		        2*M2*Q2h*Q2k*pow(l1k1,2) - 16*m2*Q2h*S*pow(l1k1,2) -
		        4*Q2e*Q2h*S*pow(l1k1,2) + 4*Q2h*Q2k*S*pow(l1k1,2) +
		        4*M2*Q2h*pow(l1k1,3) - 8*Q2h*S*pow(l1k1,3) +
		        4*l1k1*M2*Q2h*pow(l1k2,2) + 2*M2*Q2e*Q2h*pow(l1k2,2) -
		        2*M2*Q2h*Q2k*pow(l1k2,2) - 8*l1k1*Q2h*S*pow(l1k2,2) -
		        4*Q2e*Q2h*S*pow(l1k2,2) + 4*Q2h*Q2k*S*pow(l1k2,2) -
		        32*l1k1*l1k2*m2*pow(S,2) + 32*l1k1*m4*pow(S,2) -
		        16*l1k1*m2*Q2e*pow(S,2) - 16*l1k2*m2*Q2e*pow(S,2) +
		        16*m4*Q2e*pow(S,2) + 16*l1k1*m2*Q2k*pow(S,2) +
		        16*l1k2*m2*Q2k*pow(S,2) - 16*m4*Q2k*pow(S,2) -
		        32*m2*pow(l1k1,2)*pow(S,2) - 8*Q2e*pow(l1k1,2)*pow(S,2) +
		        8*Q2k*pow(l1k1,2)*pow(S,2) - 16*pow(l1k1,3)*pow(S,2) -
		        16*l1k1*pow(l1k2,2)*pow(S,2) - 8*Q2e*pow(l1k2,2)*pow(S,2) +
		        8*Q2k*pow(l1k2,2)*pow(S,2)) - 16*m2*pow(Sk,2) - 16*k1k2*pow(Sq2,2) +
		     16*l1k2*pow(Sq2,2) - 64*m2*pow(Sq2,2) - 8*Q2h*pow(Sq2,2) +
		     8*Q2k*pow(Sq2,2) + pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
		      (32*l1k2*m4*M2 + 32*l1k1*m2*M2*Q2e + 48*l1k2*m2*M2*Q2e +
		        16*l1k2*m4*Q2h + 8*l1k2*m2*M2*Q2h - 8*l1k2*M2*Q2e*Q2h -
		        32*l1k1*m2*M2*Q2k + 8*l1k2*m2*M2*Q2k - 16*l1k1*M2*Q2e*Q2k +
		        8*l1k2*M2*Q2e*Q2k - 32*l1k2*m2*Q2h*Q2k + 8*l1k2*M2*Q2h*Q2k -
		        20*m2*Q2e*Q2h*Q2k - 8*M2*Q2e*Q2h*Q2k + 64*l1k1*m2*Q2h*S +
		        32*l1k1*Q2e*Q2h*S - 8*l1k2*Q2e*Q2h*S + 16*m2*Q2e*Q2h*S -
		        32*l1k2*m2*Q2k*S + 8*l1k2*Q2h*Q2k*S - 16*m2*Q2h*Q2k*S -
		        8*Q2e*Q2h*Q2k*S - 32*l1k2*m4*Sk - 64*l1k2*m2*Q2h*Sk -
		        40*m2*Q2e*Q2h*Sk - 16*l1k2*m2*Q2k*Sk + 64*m4*Q2k*Sk +
		        64*m2*Q2e*Q2k*Sk - 40*m2*Q2h*Q2k*Sk + 4*Q2e*Q2h*Q2k*Sk -
		        64*l1k2*m2*S*Sk - 16*l1k1*l1k2*Q2h*Sq2 + 32*l1k1*m2*Q2h*Sq2 +
		        48*l1k2*m2*Q2h*Sq2 + 16*l1k1*Q2e*Q2h*Sq2 + 16*l1k2*Q2e*Q2h*Sq2 +
		        16*m2*Q2e*Q2h*Sq2 + 32*l1k2*m2*Q2k*Sq2 - 64*m4*Q2k*Sq2 -
		        24*l1k2*Q2e*Q2k*Sq2 - 32*m2*Q2e*Q2k*Sq2 + 24*l1k2*Q2h*Q2k*Sq2 +
		        8*m2*Q2h*Q2k*Sq2 + 4*Q2e*Q2h*Q2k*Sq2 + 32*l1k2*Q2h*S*Sq2 -
		        64*m2*Q2h*S*Sq2 - 32*Q2e*Q2h*S*Sq2 + 64*l1k2*m2*Sk*Sq2 -
		        128*m4*Sk*Sq2 - 48*l1k2*Q2e*Sk*Sq2 - 64*m2*Q2e*Sk*Sq2 +
		        16*l1k2*Q2h*Sk*Sq2 + 32*m2*Q2h*Sk*Sq2 + 16*Q2e*Q2h*Sk*Sq2 +
		        64*m2*M2*pow(l1k1,2) + 32*M2*Q2e*pow(l1k1,2) +
		        64*m2*M2*pow(l1k2,2) + 32*M2*Q2e*pow(l1k2,2) +
		        16*Q2h*S*pow(l1k2,2) + 32*Q2h*Sq2*pow(l1k2,2) +
		        16*l1k1*M2*pow(Q2e,2) + 8*l1k2*M2*pow(Q2e,2) +
		        8*m2*M2*pow(Q2e,2) + 6*Q2h*Q2k*pow(Q2e,2) + 8*Q2h*S*pow(Q2e,2) +
		        12*Q2h*Sk*pow(Q2e,2) - 16*Q2k*Sk*pow(Q2e,2) + 8*Q2h*Sq2*pow(Q2e,2) +
		        4*M2*pow(Q2e,3) + 24*l1k2*m2*pow(Q2h,2) -
		        16*l1k2*M2*pow(Q2h,2) + 8*m2*M2*pow(Q2h,2) +
		        12*m2*Q2e*pow(Q2h,2) + 8*M2*Q2e*pow(Q2h,2) +
		        4*l1k2*Q2k*pow(Q2h,2) + 20*m2*Q2k*pow(Q2h,2) -
		        6*Q2e*Q2k*pow(Q2h,2) + 8*l1k2*Sk*pow(Q2h,2) + 32*m2*Sk*pow(Q2h,2) -
		        8*Q2e*Sk*pow(Q2h,2) - 24*l1k2*Sq2*pow(Q2h,2) -
		        8*m2*Sq2*pow(Q2h,2) - 4*Q2e*Sq2*pow(Q2h,2) -
		        2*pow(Q2e,2)*pow(Q2h,2) - 4*l1k2*pow(Q2h,3) - 8*m2*pow(Q2h,3) +
		        4*Q2e*pow(Q2h,3) + 16*m4*pow(Q2k,2) + 16*m2*Q2e*pow(Q2k,2) +
		        4*M2*Q2e*pow(Q2k,2) - 12*m2*Q2h*pow(Q2k,2) +
		        2*Q2e*Q2h*pow(Q2k,2) - 4*pow(Q2e,2)*pow(Q2k,2) -
		        64*l1k2*m2*pow(S,2) - 64*m2*Q2h*pow(S,2) - 32*Q2e*Q2h*pow(S,2) -
		        32*l1k2*m2*pow(Sk,2) + 64*m4*pow(Sk,2) + 64*m2*Q2e*pow(Sk,2) -
		        32*m2*Q2h*pow(Sk,2) - 16*pow(Q2e,2)*pow(Sk,2) -
		        64*l1k2*m2*pow(Sq2,2) + 64*m4*pow(Sq2,2) + 32*m2*Q2e*pow(Sq2,2) +
		        32*l1k2*Q2h*pow(Sq2,2) - 32*m2*Q2h*pow(Sq2,2) -
		        16*Q2e*Q2h*pow(Sq2,2) - 32*pow(l1k2,2)*pow(Sq2,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*
		      (-24*l1k1*l1k2*m2*M2 - 16*l1k2*m4*M2 - 8*l1k1*l1k2*M2*Q2e -
		        20*l1k1*m2*M2*Q2e - 16*l1k2*m2*M2*Q2e - 8*l1k2*m4*Q2h -
		        8*l1k1*l1k2*M2*Q2h + 4*l1k1*m2*M2*Q2h - 12*l1k2*m2*M2*Q2h +
		        2*l1k1*M2*Q2e*Q2h - 2*l1k2*M2*Q2e*Q2h + 16*l1k1*m2*M2*Q2k +
		        4*l1k2*m2*M2*Q2k + 6*l1k1*M2*Q2e*Q2k + 2*l1k2*M2*Q2e*Q2k +
		        12*l1k1*m2*Q2h*Q2k + 8*l1k2*m2*Q2h*Q2k - 8*m4*Q2h*Q2k -
		        6*l1k1*M2*Q2h*Q2k + 6*l1k2*M2*Q2h*Q2k + 4*m2*M2*Q2h*Q2k -
		        2*l1k2*Q2e*Q2h*Q2k + 5*m2*Q2e*Q2h*Q2k + 6*M2*Q2e*Q2h*Q2k -
		        8*l1k1*l1k2*Q2h*S - 32*l1k1*m2*Q2h*S + 16*m4*Q2h*S -
		        8*l1k1*Q2e*Q2h*S + 4*l1k2*Q2e*Q2h*S - 12*m2*Q2e*Q2h*S +
		        16*l1k1*m2*Q2k*S + 32*l1k2*m2*Q2k*S - 32*m4*Q2k*S -
		        12*l1k1*Q2e*Q2k*S - 12*l1k2*Q2e*Q2k*S - 8*m2*Q2e*Q2k*S +
		        8*l1k1*Q2h*Q2k*S + 8*l1k2*Q2h*Q2k*S + 16*l1k2*m4*Sk +
		        24*l1k1*m2*Q2h*Sk + 16*l1k2*m2*Q2h*Sk - 16*m4*Q2h*Sk -
		        4*l1k2*Q2e*Q2h*Sk + 10*m2*Q2e*Q2h*Sk + 8*l1k2*m2*Q2k*Sk -
		        12*m2*Q2e*Q2k*Sk + 12*m2*Q2h*Q2k*Sk + 2*Q2e*Q2h*Q2k*Sk +
		        32*l1k1*m2*S*Sk + 64*l1k2*m2*S*Sk - 64*m4*S*Sk -
		        24*l1k1*Q2e*S*Sk - 24*l1k2*Q2e*S*Sk - 16*m2*Q2e*S*Sk +
		        8*l1k1*Q2h*S*Sk + 8*l1k2*Q2h*S*Sk - 16*m2*Q2h*S*Sk -
		        8*Q2e*Q2h*S*Sk + 8*l1k1*l1k2*Q2h*Sq2 - 32*l1k1*m2*Q2h*Sq2 -
		        8*l1k2*m2*Q2h*Sq2 + 8*m4*Q2h*Sq2 - 12*l1k1*Q2e*Q2h*Sq2 -
		        8*m2*Q2e*Q2h*Sq2 + 4*l1k1*Q2h*Q2k*Sq2 - 4*l1k2*Q2h*Q2k*Sq2 -
		        2*Q2e*Q2h*Q2k*Sq2 + 32*l1k1*l1k2*S*Sq2 - 64*l1k1*m2*S*Sq2 +
		        32*m4*S*Sq2 - 16*l1k2*Q2e*S*Sq2 - 16*l1k1*Q2h*S*Sq2 +
		        48*l1k2*Q2h*S*Sq2 + 28*Q2e*Q2h*S*Sq2 + 16*l1k1*Q2k*S*Sq2 -
		        16*l1k2*Q2k*S*Sq2 + 16*m2*Q2k*S*Sq2 + 4*Q2e*Q2k*S*Sq2 -
		        4*Q2h*Q2k*S*Sq2 - 16*m2*Q2h*Sk*Sq2 - 12*Q2e*Q2h*Sk*Sq2 -
		        24*m2*M2*pow(l1k1,2) - 12*M2*Q2e*pow(l1k1,2) +
		        8*M2*Q2h*pow(l1k1,2) - 8*Q2h*S*pow(l1k1,2) -
		        8*Q2h*Sq2*pow(l1k1,2) - 32*S*Sq2*pow(l1k1,2) -
		        8*m2*M2*pow(l1k2,2) + 4*M2*Q2e*pow(l1k2,2) +
		        16*M2*Q2h*pow(l1k2,2) - 16*Q2h*Sq2*pow(l1k2,2) -
		        64*S*Sq2*pow(l1k2,2) - 8*l1k1*M2*pow(Q2e,2) -
		        4*m2*M2*pow(Q2e,2) + M2*Q2h*pow(Q2e,2) - M2*Q2k*pow(Q2e,2) -
		        Q2h*Q2k*pow(Q2e,2) + 2*Q2h*S*pow(Q2e,2) - 4*Q2k*S*pow(Q2e,2) -
		        2*Q2h*Sk*pow(Q2e,2) - 8*S*Sk*pow(Q2e,2) - 2*Q2h*Sq2*pow(Q2e,2) -
		        M2*pow(Q2e,3) - 8*l1k1*l1k2*pow(Q2h,2) - 4*l1k1*m2*pow(Q2h,2) -
		        8*l1k2*m2*pow(Q2h,2) + 4*m4*pow(Q2h,2) + 10*l1k1*M2*pow(Q2h,2) -
		        10*l1k2*M2*pow(Q2h,2) - 4*m2*Q2e*pow(Q2h,2) -
		        6*M2*Q2e*pow(Q2h,2) - 2*l1k1*Q2k*pow(Q2h,2) +
		        4*l1k2*Q2k*pow(Q2h,2) - 5*m2*Q2k*pow(Q2h,2) -
		        5*M2*Q2k*pow(Q2h,2) + 2*Q2e*Q2k*pow(Q2h,2) + 4*l1k1*S*pow(Q2h,2) +
		        12*m2*S*pow(Q2h,2) - 2*Q2e*S*pow(Q2h,2) + 4*Q2k*S*pow(Q2h,2) -
		        4*l1k1*Sk*pow(Q2h,2) - 10*m2*Sk*pow(Q2h,2) + 2*Q2e*Sk*pow(Q2h,2) -
		        2*Q2k*Sk*pow(Q2h,2) + 16*S*Sk*pow(Q2h,2) + 8*l1k1*Sq2*pow(Q2h,2) +
		        8*l1k2*Sq2*pow(Q2h,2) + 4*m2*Sq2*pow(Q2h,2) +
		        2*Q2k*Sq2*pow(Q2h,2) - 28*S*Sq2*pow(Q2h,2) + 12*Sk*Sq2*pow(Q2h,2) +
		        8*pow(l1k2,2)*pow(Q2h,2) + pow(Q2e,2)*pow(Q2h,2) +
		        2*l1k1*pow(Q2h,3) - 4*l1k2*pow(Q2h,3) + 4*m2*pow(Q2h,3) +
		        6*M2*pow(Q2h,3) - 2*Q2e*pow(Q2h,3) - Q2k*pow(Q2h,3) +
		        2*Sq2*pow(Q2h,3) + pow(Q2h,4) - 3*m2*Q2e*pow(Q2k,2) -
		        M2*Q2e*pow(Q2k,2) + 3*m2*Q2h*pow(Q2k,2) + M2*Q2h*pow(Q2k,2) +
		        32*l1k1*l1k2*pow(S,2) - 64*l1k1*m2*pow(S,2) +
		        64*l1k2*m2*pow(S,2) - 8*l1k1*Q2e*pow(S,2) - 16*l1k2*Q2e*pow(S,2) -
		        16*m2*Q2e*pow(S,2) - 24*l1k1*Q2h*pow(S,2) + 48*l1k2*Q2h*pow(S,2) +
		        20*Q2e*Q2h*pow(S,2) + 16*l1k1*Q2k*pow(S,2) - 16*l1k2*Q2k*pow(S,2) +
		        16*m2*Q2k*pow(S,2) + 4*Q2e*Q2k*pow(S,2) - 4*Q2h*Q2k*pow(S,2) -
		        48*pow(l1k1,2)*pow(S,2) - 48*pow(l1k2,2)*pow(S,2) -
		        20*pow(Q2h,2)*pow(S,2) + 16*l1k2*m2*pow(Sk,2) -
		        12*m2*Q2e*pow(Sk,2) + 12*m2*Q2h*pow(Sk,2) + 4*Q2e*Q2h*pow(Sk,2) -
		        4*pow(Q2h,2)*pow(Sk,2) + 8*m2*Q2h*pow(Sq2,2) +
		        12*Q2e*Q2h*pow(Sq2,2) - 12*pow(Q2h,2)*pow(Sq2,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (-16*k1k2*l1k2*m2*M2 - 32*l1k2*m4*M2 + 8*k1k2*l1k2*M2*Q2e -
		        8*k1k2*m2*M2*Q2e - 24*l1k2*m2*M2*Q2e - 16*m4*M2*Q2e -
		        16*l1k2*m4*Q2h - 8*k1k2*l1k2*M2*Q2h + 8*k1k2*m2*M2*Q2h -
		        32*l1k2*m2*M2*Q2h + 16*m4*M2*Q2h - 8*m4*Q2e*Q2h +
		        8*l1k2*M2*Q2e*Q2h - 4*m2*M2*Q2e*Q2h - 8*l1k2*m2*M2*Q2k -
		        8*l1k2*M2*Q2e*Q2k - 12*m2*M2*Q2e*Q2k + 16*l1k2*m2*Q2h*Q2k +
		        8*l1k2*M2*Q2h*Q2k + 12*m2*M2*Q2h*Q2k - 4*l1k2*Q2e*Q2h*Q2k +
		        10*m2*Q2e*Q2h*Q2k - 16*k1k2*l1k2*Q2h*S - 32*l1k2*m2*Q2h*S -
		        8*l1k2*Q2e*Q2h*S + 8*m2*Q2e*Q2h*S - 64*l1k2*m2*Q2k*S +
		        24*l1k2*Q2e*Q2k*S - 16*m2*Q2e*Q2k*S + 16*m2*Q2h*Q2k*S -
		        8*Q2e*Q2h*Q2k*S + 32*l1k2*m4*Sk + 16*m4*Q2e*Sk +
		        32*l1k2*m2*Q2h*Sk - 16*m4*Q2h*Sk - 8*l1k2*Q2e*Q2h*Sk +
		        20*m2*Q2e*Q2h*Sk + 16*l1k2*m2*Q2k*Sk - 16*m2*Q2e*Q2k*Sk +
		        16*m2*Q2h*Q2k*Sk - 128*l1k2*m2*S*Sk + 48*l1k2*Q2e*S*Sk -
		        32*m2*Q2e*S*Sk - 16*l1k2*Q2h*S*Sk + 32*m2*Q2h*S*Sk -
		        16*Q2e*Q2h*S*Sk - 16*k1k2*l1k2*Q2h*Sq2 - 48*l1k2*m2*Q2h*Sq2 -
		        8*l1k2*Q2e*Q2h*Sq2 - 64*l1k2*m2*Q2k*Sq2 + 24*l1k2*Q2e*Q2k*Sq2 -
		        16*m2*Q2e*Q2k*Sq2 + 16*m2*Q2h*Q2k*Sq2 - 8*Q2e*Q2h*Q2k*Sq2 +
		        128*k1k2*l1k2*S*Sq2 + 384*l1k2*m2*S*Sq2 + 32*k1k2*Q2e*S*Sq2 -
		        64*l1k2*Q2e*S*Sq2 + 64*m2*Q2e*S*Sq2 - 32*k1k2*Q2h*S*Sq2 +
		        160*l1k2*Q2h*S*Sq2 - 64*m2*Q2h*S*Sq2 + 32*Q2e*Q2h*S*Sq2 -
		        64*l1k2*Q2k*S*Sq2 + 16*Q2e*Q2k*S*Sq2 - 16*Q2h*Q2k*S*Sq2 -
		        128*l1k2*m2*Sk*Sq2 + 48*l1k2*Q2e*Sk*Sq2 - 32*m2*Q2e*Sk*Sq2 -
		        16*l1k2*Q2h*Sk*Sq2 + 32*m2*Q2h*Sk*Sq2 - 16*Q2e*Q2h*Sk*Sq2 -
		        32*m2*M2*pow(l1k2,2) - 24*M2*Q2e*pow(l1k2,2) +
		        24*M2*Q2h*pow(l1k2,2) + 32*Q2h*S*pow(l1k2,2) +
		        16*Q2h*Sq2*pow(l1k2,2) - 192*S*Sq2*pow(l1k2,2) -
		        4*m2*M2*pow(Q2e,2) - 2*M2*Q2h*pow(Q2e,2) - 4*Q2h*S*pow(Q2e,2) +
		        4*Q2k*S*pow(Q2e,2) + 8*S*Sk*pow(Q2e,2) - 4*Q2h*Sq2*pow(Q2e,2) +
		        4*Q2k*Sq2*pow(Q2e,2) + 8*Sk*Sq2*pow(Q2e,2) -
		        16*l1k2*m2*pow(Q2h,2) + 8*m4*pow(Q2h,2) -
		        16*l1k2*M2*pow(Q2h,2) + 8*m2*M2*pow(Q2h,2) -
		        8*m2*Q2e*pow(Q2h,2) - 10*m2*Q2k*pow(Q2h,2) -
		        16*l1k2*S*pow(Q2h,2) - 8*m2*S*pow(Q2h,2) + 4*Q2e*S*pow(Q2h,2) +
		        4*Q2k*S*pow(Q2h,2) - 20*m2*Sk*pow(Q2h,2) + 8*S*Sk*pow(Q2h,2) -
		        8*l1k2*Sq2*pow(Q2h,2) + 4*Q2e*Sq2*pow(Q2h,2) +
		        4*Q2k*Sq2*pow(Q2h,2) - 32*S*Sq2*pow(Q2h,2) + 8*Sk*Sq2*pow(Q2h,2) +
		        8*m2*pow(Q2h,3) + 2*M2*pow(Q2h,3) - 6*m2*Q2e*pow(Q2k,2) +
		        6*m2*Q2h*pow(Q2k,2) + 64*k1k2*l1k2*pow(S,2) +
		        192*l1k2*m2*pow(S,2) + 16*k1k2*Q2e*pow(S,2) -
		        48*l1k2*Q2e*pow(S,2) + 32*m2*Q2e*pow(S,2) - 16*k1k2*Q2h*pow(S,2) +
		        112*l1k2*Q2h*pow(S,2) - 32*m2*Q2h*pow(S,2) + 24*Q2e*Q2h*pow(S,2) -
		        32*l1k2*Q2k*pow(S,2) + 8*Q2e*Q2k*pow(S,2) - 8*Q2h*Q2k*pow(S,2) -
		        128*pow(l1k2,2)*pow(S,2) - 24*pow(Q2h,2)*pow(S,2) +
		        32*l1k2*m2*pow(Sk,2) - 8*m2*Q2e*pow(Sk,2) + 8*m2*Q2h*pow(Sk,2) +
		        64*k1k2*l1k2*pow(Sq2,2) + 192*l1k2*m2*pow(Sq2,2) +
		        16*k1k2*Q2e*pow(Sq2,2) - 16*l1k2*Q2e*pow(Sq2,2) +
		        32*m2*Q2e*pow(Sq2,2) - 16*k1k2*Q2h*pow(Sq2,2) +
		        48*l1k2*Q2h*pow(Sq2,2) - 32*m2*Q2h*pow(Sq2,2) +
		        8*Q2e*Q2h*pow(Sq2,2) - 32*l1k2*Q2k*pow(Sq2,2) +
		        8*Q2e*Q2k*pow(Sq2,2) - 8*Q2h*Q2k*pow(Sq2,2) -
		        64*pow(l1k2,2)*pow(Sq2,2) - 8*pow(Q2h,2)*pow(Sq2,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		      (8*k1k2*l1k2*M2*Q2e*Q2h + 16*l1k2*m2*M2*Q2e*Q2h -
		        4*l1k2*M2*Q2e*Q2h*Q2k + 16*k1k2*l1k2*Q2e*Q2h*S +
		        32*l1k2*m2*Q2e*Q2h*S - 8*l1k2*Q2e*Q2h*Q2k*S +
		        16*k1k2*l1k2*Q2e*Q2h*Sq2 + 32*l1k2*m2*Q2e*Q2h*Sq2 -
		        8*l1k2*Q2e*Q2h*Q2k*Sq2 - 64*k1k2*l1k2*Q2e*S*Sq2 -
		        128*l1k2*m2*Q2e*S*Sq2 + 64*k1k2*l1k2*Q2h*S*Sq2 +
		        128*l1k2*m2*Q2h*S*Sq2 - 32*l1k2*Q2e*Q2h*S*Sq2 +
		        32*l1k2*Q2e*Q2k*S*Sq2 - 32*l1k2*Q2h*Q2k*S*Sq2 +
		        16*k1k2*M2*Q2h*pow(l1k2,2) + 32*m2*M2*Q2h*pow(l1k2,2) -
		        8*M2*Q2e*Q2h*pow(l1k2,2) - 8*M2*Q2h*Q2k*pow(l1k2,2) +
		        32*k1k2*Q2h*S*pow(l1k2,2) + 64*m2*Q2h*S*pow(l1k2,2) -
		        16*Q2e*Q2h*S*pow(l1k2,2) - 16*Q2h*Q2k*S*pow(l1k2,2) +
		        32*k1k2*Q2h*Sq2*pow(l1k2,2) + 64*m2*Q2h*Sq2*pow(l1k2,2) -
		        16*Q2e*Q2h*Sq2*pow(l1k2,2) - 16*Q2h*Q2k*Sq2*pow(l1k2,2) -
		        128*k1k2*S*Sq2*pow(l1k2,2) - 256*m2*S*Sq2*pow(l1k2,2) +
		        64*Q2e*S*Sq2*pow(l1k2,2) - 128*Q2h*S*Sq2*pow(l1k2,2) +
		        64*Q2k*S*Sq2*pow(l1k2,2) - 16*M2*Q2h*pow(l1k2,3) -
		        32*Q2h*S*pow(l1k2,3) - 32*Q2h*Sq2*pow(l1k2,3) +
		        128*S*Sq2*pow(l1k2,3) - 8*k1k2*l1k2*M2*pow(Q2h,2) -
		        16*l1k2*m2*M2*pow(Q2h,2) + 4*l1k2*M2*Q2e*pow(Q2h,2) +
		        4*l1k2*M2*Q2k*pow(Q2h,2) - 16*k1k2*l1k2*S*pow(Q2h,2) -
		        32*l1k2*m2*S*pow(Q2h,2) + 8*l1k2*Q2e*S*pow(Q2h,2) +
		        8*l1k2*Q2k*S*pow(Q2h,2) - 16*k1k2*l1k2*Sq2*pow(Q2h,2) -
		        32*l1k2*m2*Sq2*pow(Q2h,2) + 8*l1k2*Q2e*Sq2*pow(Q2h,2) +
		        8*l1k2*Q2k*Sq2*pow(Q2h,2) + 32*l1k2*S*Sq2*pow(Q2h,2) +
		        16*M2*pow(l1k2,2)*pow(Q2h,2) + 32*S*pow(l1k2,2)*pow(Q2h,2) +
		        32*Sq2*pow(l1k2,2)*pow(Q2h,2) - 4*l1k2*M2*pow(Q2h,3) -
		        8*l1k2*S*pow(Q2h,3) - 8*l1k2*Sq2*pow(Q2h,3) -
		        32*k1k2*l1k2*Q2e*pow(S,2) - 64*l1k2*m2*Q2e*pow(S,2) +
		        32*k1k2*l1k2*Q2h*pow(S,2) + 64*l1k2*m2*Q2h*pow(S,2) -
		        16*l1k2*Q2e*Q2h*pow(S,2) + 16*l1k2*Q2e*Q2k*pow(S,2) -
		        16*l1k2*Q2h*Q2k*pow(S,2) - 64*k1k2*pow(l1k2,2)*pow(S,2) -
		        128*m2*pow(l1k2,2)*pow(S,2) + 32*Q2e*pow(l1k2,2)*pow(S,2) -
		        64*Q2h*pow(l1k2,2)*pow(S,2) + 32*Q2k*pow(l1k2,2)*pow(S,2) +
		        64*pow(l1k2,3)*pow(S,2) + 16*l1k2*pow(Q2h,2)*pow(S,2) -
		        32*k1k2*l1k2*Q2e*pow(Sq2,2) - 64*l1k2*m2*Q2e*pow(Sq2,2) +
		        32*k1k2*l1k2*Q2h*pow(Sq2,2) + 64*l1k2*m2*Q2h*pow(Sq2,2) -
		        16*l1k2*Q2e*Q2h*pow(Sq2,2) + 16*l1k2*Q2e*Q2k*pow(Sq2,2) -
		        16*l1k2*Q2h*Q2k*pow(Sq2,2) - 64*k1k2*pow(l1k2,2)*pow(Sq2,2) -
		        128*m2*pow(l1k2,2)*pow(Sq2,2) + 32*Q2e*pow(l1k2,2)*pow(Sq2,2) -
		        64*Q2h*pow(l1k2,2)*pow(Sq2,2) + 32*Q2k*pow(l1k2,2)*pow(Sq2,2) +
		        64*pow(l1k2,3)*pow(Sq2,2) + 16*l1k2*pow(Q2h,2)*pow(Sq2,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
		      (32*l1k2*m4*M2*Q2h + 24*l1k2*m2*M2*Q2e*Q2h -
		        40*l1k2*m2*M2*Q2e*Q2k - 32*l1k2*m4*Q2h*Q2k +
		        24*l1k2*m2*M2*Q2h*Q2k + 20*l1k2*m2*Q2e*Q2h*Q2k -
		        8*m4*Q2e*Q2h*Q2k + 32*l1k2*M2*Q2e*Q2h*Q2k +
		        16*l1k2*m2*Q2e*Q2h*S + 128*l1k2*m4*Q2k*S + 32*l1k2*m2*Q2e*Q2k*S +
		        32*m4*Q2e*Q2k*S + 64*l1k2*m2*Q2h*Q2k*S - 32*m4*Q2h*Q2k*S -
		        24*l1k2*Q2e*Q2h*Q2k*S + 24*m2*Q2e*Q2h*Q2k*S - 64*l1k2*m4*Q2h*Sk +
		        40*l1k2*m2*Q2e*Q2h*Sk - 16*m4*Q2e*Q2h*Sk -
		        48*l1k2*m2*Q2e*Q2k*Sk + 96*l1k2*m2*Q2h*Q2k*Sk +
		        8*l1k2*Q2e*Q2h*Q2k*Sk + 60*m2*Q2e*Q2h*Q2k*Sk + 256*l1k2*m4*S*Sk +
		        64*l1k2*m2*Q2e*S*Sk + 64*m4*Q2e*S*Sk + 128*l1k2*m2*Q2h*S*Sk -
		        64*m4*Q2h*S*Sk - 16*l1k2*Q2e*Q2h*S*Sk - 64*l1k2*m2*Q2k*S*Sk +
		        48*l1k2*Q2e*Q2k*S*Sk - 32*m2*Q2e*Q2k*S*Sk - 16*l1k2*Q2h*Q2k*S*Sk +
		        32*m2*Q2h*Q2k*S*Sk - 16*Q2e*Q2h*Q2k*S*Sk + 32*l1k2*m4*Q2h*Sq2 -
		        16*l1k2*m2*Q2e*Q2h*Sq2 + 16*m4*Q2e*Q2h*Sq2 +
		        128*l1k2*m4*Q2k*Sq2 + 32*l1k2*m2*Q2e*Q2k*Sq2 +
		        32*m4*Q2e*Q2k*Sq2 - 32*m4*Q2h*Q2k*Sq2 - 56*l1k2*Q2e*Q2h*Q2k*Sq2 -
		        28*m2*Q2e*Q2h*Q2k*Sq2 - 384*l1k2*m4*S*Sq2 - 64*m4*Q2e*S*Sq2 -
		        128*l1k2*m2*Q2h*S*Sq2 + 64*m4*Q2h*S*Sq2 + 16*l1k2*Q2e*Q2h*S*Sq2 +
		        128*m2*Q2e*Q2h*S*Sq2 + 64*l1k2*m2*Q2k*S*Sq2 +
		        48*l1k2*Q2e*Q2k*S*Sq2 - 112*l1k2*Q2h*Q2k*S*Sq2 +
		        256*l1k2*m4*Sk*Sq2 + 64*l1k2*m2*Q2e*Sk*Sq2 + 64*m4*Q2e*Sk*Sq2 +
		        64*l1k2*m2*Q2h*Sk*Sq2 - 64*m4*Q2h*Sk*Sq2 -
		        64*l1k2*Q2e*Q2h*Sk*Sq2 - 32*m2*Q2e*Q2h*Sk*Sq2 -
		        64*l1k2*m2*Q2k*Sk*Sq2 + 48*l1k2*Q2e*Q2k*Sk*Sq2 +
		        16*m2*Q2e*Q2k*Sk*Sq2 - 16*l1k2*Q2h*Q2k*Sk*Sq2 -
		        16*m2*Q2h*Q2k*Sk*Sq2 - 16*Q2e*Q2h*Q2k*Sk*Sq2 -
		        80*m2*M2*Q2e*pow(l1k2,2) + 48*m2*M2*Q2h*pow(l1k2,2) +
		        40*M2*Q2e*Q2h*pow(l1k2,2) - 32*m2*M2*Q2k*pow(l1k2,2) -
		        24*M2*Q2e*Q2k*pow(l1k2,2) + 48*m2*Q2h*Q2k*pow(l1k2,2) -
		        24*M2*Q2h*Q2k*pow(l1k2,2) - 64*m2*Q2k*S*pow(l1k2,2) +
		        48*Q2e*Q2k*S*pow(l1k2,2) - 64*Q2h*Q2k*S*pow(l1k2,2) +
		        96*m2*Q2h*Sk*pow(l1k2,2) - 128*m2*S*Sk*pow(l1k2,2) +
		        96*Q2e*S*Sk*pow(l1k2,2) - 32*Q2h*S*Sk*pow(l1k2,2) -
		        128*m2*Q2h*Sq2*pow(l1k2,2) - 48*Q2e*Q2h*Sq2*pow(l1k2,2) -
		        64*m2*Q2k*Sq2*pow(l1k2,2) + 48*Q2e*Q2k*Sq2*pow(l1k2,2) -
		        80*Q2h*Q2k*Sq2*pow(l1k2,2) + 256*m2*S*Sq2*pow(l1k2,2) +
		        64*Q2e*S*Sq2*pow(l1k2,2) - 256*Q2h*S*Sq2*pow(l1k2,2) +
		        64*Q2k*S*Sq2*pow(l1k2,2) - 128*m2*Sk*Sq2*pow(l1k2,2) +
		        96*Q2e*Sk*Sq2*pow(l1k2,2) - 32*Q2h*Sk*Sq2*pow(l1k2,2) -
		        96*m2*M2*pow(l1k2,3) - 48*M2*Q2e*pow(l1k2,3) -
		        16*M2*Q2h*pow(l1k2,3) - 64*Q2h*S*pow(l1k2,3) -
		        96*Q2h*Sq2*pow(l1k2,3) + 128*S*Sq2*pow(l1k2,3) -
		        16*l1k2*m2*M2*pow(Q2e,2) + 20*l1k2*M2*Q2h*pow(Q2e,2) +
		        12*m2*M2*Q2h*pow(Q2e,2) - 20*l1k2*M2*Q2k*pow(Q2e,2) -
		        4*m2*M2*Q2k*pow(Q2e,2) - 4*l1k2*Q2h*Q2k*pow(Q2e,2) +
		        8*m2*Q2h*Q2k*pow(Q2e,2) - 2*M2*Q2h*Q2k*pow(Q2e,2) -
		        8*l1k2*Q2h*S*pow(Q2e,2) + 16*l1k2*Q2k*S*pow(Q2e,2) +
		        16*m2*Q2k*S*pow(Q2e,2) - 20*Q2h*Q2k*S*pow(Q2e,2) -
		        8*l1k2*Q2h*Sk*pow(Q2e,2) + 16*m2*Q2h*Sk*pow(Q2e,2) -
		        24*m2*Q2k*Sk*pow(Q2e,2) - 8*Q2h*Q2k*Sk*pow(Q2e,2) +
		        32*l1k2*S*Sk*pow(Q2e,2) + 32*m2*S*Sk*pow(Q2e,2) -
		        48*Q2h*S*Sk*pow(Q2e,2) + 16*Q2k*S*Sk*pow(Q2e,2) -
		        16*l1k2*Q2h*Sq2*pow(Q2e,2) + 16*l1k2*Q2k*Sq2*pow(Q2e,2) +
		        16*m2*Q2k*Sq2*pow(Q2e,2) - 16*Q2h*Q2k*Sq2*pow(Q2e,2) -
		        64*m2*S*Sq2*pow(Q2e,2) - 24*Q2h*S*Sq2*pow(Q2e,2) +
		        8*Q2k*S*Sq2*pow(Q2e,2) + 32*l1k2*Sk*Sq2*pow(Q2e,2) +
		        32*m2*Sk*Sq2*pow(Q2e,2) - 24*Q2h*Sk*Sq2*pow(Q2e,2) +
		        16*Q2k*Sk*Sq2*pow(Q2e,2) - 32*M2*pow(l1k2,2)*pow(Q2e,2) -
		        4*l1k2*M2*pow(Q2e,3) + 2*M2*Q2h*pow(Q2e,3) -
		        2*M2*Q2k*pow(Q2e,3) + 16*l1k2*m4*pow(Q2h,2) -
		        16*l1k2*m2*Q2e*pow(Q2h,2) - 32*l1k2*M2*Q2e*pow(Q2h,2) -
		        24*m2*M2*Q2e*pow(Q2h,2) - 52*l1k2*m2*Q2k*pow(Q2h,2) +
		        8*m4*Q2k*pow(Q2h,2) + 12*l1k2*M2*Q2k*pow(Q2h,2) +
		        4*m2*M2*Q2k*pow(Q2h,2) + 8*l1k2*Q2e*Q2k*pow(Q2h,2) -
		        20*m2*Q2e*Q2k*pow(Q2h,2) + 14*M2*Q2e*Q2k*pow(Q2h,2) -
		        48*l1k2*m2*S*pow(Q2h,2) + 8*l1k2*Q2e*S*pow(Q2h,2) -
		        24*m2*Q2e*S*pow(Q2h,2) + 8*l1k2*Q2k*S*pow(Q2h,2) -
		        40*m2*Q2k*S*pow(Q2h,2) + 28*Q2e*Q2k*S*pow(Q2h,2) -
		        88*l1k2*m2*Sk*pow(Q2h,2) + 16*m4*Sk*pow(Q2h,2) +
		        8*l1k2*Q2e*Sk*pow(Q2h,2) - 36*m2*Q2e*Sk*pow(Q2h,2) -
		        16*l1k2*Q2k*Sk*pow(Q2h,2) - 36*m2*Q2k*Sk*pow(Q2h,2) +
		        12*Q2e*Q2k*Sk*pow(Q2h,2) - 48*l1k2*S*Sk*pow(Q2h,2) -
		        32*m2*S*Sk*pow(Q2h,2) + 80*Q2e*S*Sk*pow(Q2h,2) +
		        32*l1k2*m2*Sq2*pow(Q2h,2) - 16*m4*Sq2*pow(Q2h,2) +
		        32*l1k2*Q2e*Sq2*pow(Q2h,2) + 4*m2*Q2e*Sq2*pow(Q2h,2) +
		        40*l1k2*Q2k*Sq2*pow(Q2h,2) + 12*m2*Q2k*Sq2*pow(Q2h,2) +
		        20*Q2e*Q2k*Sq2*pow(Q2h,2) + 48*l1k2*S*Sq2*pow(Q2h,2) -
		        64*m2*S*Sq2*pow(Q2h,2) + 48*Q2e*S*Sq2*pow(Q2h,2) -
		        8*Q2k*S*Sq2*pow(Q2h,2) + 32*Q2e*Sk*Sq2*pow(Q2h,2) -
		        16*m2*pow(l1k2,2)*pow(Q2h,2) + 56*M2*pow(l1k2,2)*pow(Q2h,2) -
		        8*Q2k*pow(l1k2,2)*pow(Q2h,2) + 48*S*pow(l1k2,2)*pow(Q2h,2) -
		        16*Sk*pow(l1k2,2)*pow(Q2h,2) + 112*Sq2*pow(l1k2,2)*pow(Q2h,2) +
		        4*l1k2*pow(Q2e,2)*pow(Q2h,2) - 2*m2*pow(Q2e,2)*pow(Q2h,2) +
		        2*M2*pow(Q2e,2)*pow(Q2h,2) + 2*Q2k*pow(Q2e,2)*pow(Q2h,2) +
		        4*S*pow(Q2e,2)*pow(Q2h,2) + 4*Sk*pow(Q2e,2)*pow(Q2h,2) +
		        4*Sq2*pow(Q2e,2)*pow(Q2h,2) + 24*l1k2*m2*pow(Q2h,3) +
		        12*m2*M2*pow(Q2h,3) - 8*l1k2*Q2e*pow(Q2h,3) +
		        2*m2*Q2e*pow(Q2h,3) - 16*M2*Q2e*pow(Q2h,3) +
		        4*l1k2*Q2k*pow(Q2h,3) + 12*m2*Q2k*pow(Q2h,3) -
		        10*M2*Q2k*pow(Q2h,3) - 8*l1k2*S*pow(Q2h,3) + 24*m2*S*pow(Q2h,3) -
		        4*Q2e*S*pow(Q2h,3) - 8*Q2k*S*pow(Q2h,3) + 8*l1k2*Sk*pow(Q2h,3) +
		        20*m2*Sk*pow(Q2h,3) - 4*Q2e*Sk*pow(Q2h,3) - 4*Q2k*Sk*pow(Q2h,3) -
		        32*S*Sk*pow(Q2h,3) - 24*l1k2*Sq2*pow(Q2h,3) - 4*m2*Sq2*pow(Q2h,3) -
		        8*Q2e*Sq2*pow(Q2h,3) - 4*Q2k*Sq2*pow(Q2h,3) - 24*S*Sq2*pow(Q2h,3) -
		        8*Sk*Sq2*pow(Q2h,3) + 8*pow(l1k2,2)*pow(Q2h,3) +
		        12*M2*pow(Q2h,4) - 2*Q2e*pow(Q2h,4) - 2*Q2k*pow(Q2h,4) +
		        4*Sq2*pow(Q2h,4) + 2*pow(Q2h,5) + 8*l1k2*m2*M2*pow(Q2k,2) -
		        12*l1k2*m2*Q2e*pow(Q2k,2) - 4*l1k2*M2*Q2e*pow(Q2k,2) +
		        36*l1k2*m2*Q2h*pow(Q2k,2) - 4*l1k2*M2*Q2h*pow(Q2k,2) +
		        20*m2*Q2e*Q2h*pow(Q2k,2) - 2*M2*Q2e*Q2h*pow(Q2k,2) -
		        32*l1k2*m2*S*pow(Q2k,2) + 24*l1k2*Q2e*S*pow(Q2k,2) -
		        16*m2*Q2e*S*pow(Q2k,2) - 16*l1k2*Q2h*S*pow(Q2k,2) +
		        16*m2*Q2h*S*pow(Q2k,2) - 8*Q2e*Q2h*S*pow(Q2k,2) -
		        24*m2*Q2e*Sk*pow(Q2k,2) + 24*m2*Q2h*Sk*pow(Q2k,2) -
		        32*l1k2*m2*Sq2*pow(Q2k,2) + 24*l1k2*Q2e*Sq2*pow(Q2k,2) +
		        8*m2*Q2e*Sq2*pow(Q2k,2) - 16*l1k2*Q2h*Sq2*pow(Q2k,2) -
		        8*m2*Q2h*Sq2*pow(Q2k,2) - 8*Q2e*Q2h*Sq2*pow(Q2k,2) -
		        6*m2*pow(Q2e,2)*pow(Q2k,2) - 2*Q2h*pow(Q2e,2)*pow(Q2k,2) +
		        8*S*pow(Q2e,2)*pow(Q2k,2) + 8*Sq2*pow(Q2e,2)*pow(Q2k,2) -
		        4*l1k2*pow(Q2h,2)*pow(Q2k,2) - 14*m2*pow(Q2h,2)*pow(Q2k,2) +
		        2*M2*pow(Q2h,2)*pow(Q2k,2) + 2*Q2e*pow(Q2h,2)*pow(Q2k,2) -
		        6*m2*Q2e*pow(Q2k,3) + 6*m2*Q2h*pow(Q2k,3) -
		        128*l1k2*m4*pow(S,2) + 64*l1k2*Q2e*Q2h*pow(S,2) +
		        96*m2*Q2e*Q2h*pow(S,2) + 32*l1k2*Q2e*Q2k*pow(S,2) -
		        64*l1k2*Q2h*Q2k*pow(S,2) + 32*Q2e*pow(l1k2,2)*pow(S,2) -
		        96*Q2h*pow(l1k2,2)*pow(S,2) - 32*m2*pow(Q2e,2)*pow(S,2) -
		        40*Q2h*pow(Q2e,2)*pow(S,2) + 8*Q2k*pow(Q2e,2)*pow(S,2) -
		        32*l1k2*pow(Q2h,2)*pow(S,2) - 64*m2*pow(Q2h,2)*pow(S,2) +
		        80*Q2e*pow(Q2h,2)*pow(S,2) - 8*Q2k*pow(Q2h,2)*pow(S,2) -
		        40*pow(Q2h,3)*pow(S,2) - 48*l1k2*m2*Q2e*pow(Sk,2) +
		        48*l1k2*m2*Q2h*pow(Sk,2) + 16*l1k2*Q2e*Q2h*pow(Sk,2) +
		        40*m2*Q2e*Q2h*pow(Sk,2) - 24*m2*Q2e*Q2k*pow(Sk,2) +
		        24*m2*Q2h*Q2k*pow(Sk,2) - 24*m2*pow(Q2e,2)*pow(Sk,2) -
		        8*Q2h*pow(Q2e,2)*pow(Sk,2) - 16*l1k2*pow(Q2h,2)*pow(Sk,2) -
		        16*m2*pow(Q2h,2)*pow(Sk,2) + 16*Q2e*pow(Q2h,2)*pow(Sk,2) -
		        8*pow(Q2h,3)*pow(Sk,2) - 256*l1k2*m4*pow(Sq2,2) -
		        64*m4*Q2e*pow(Sq2,2) - 96*l1k2*m2*Q2h*pow(Sq2,2) +
		        64*m4*Q2h*pow(Sq2,2) + 48*m2*Q2e*Q2h*pow(Sq2,2) +
		        64*l1k2*m2*Q2k*pow(Sq2,2) + 16*l1k2*Q2e*Q2k*pow(Sq2,2) -
		        48*l1k2*Q2h*Q2k*pow(Sq2,2) + 256*m2*pow(l1k2,2)*pow(Sq2,2) +
		        32*Q2e*pow(l1k2,2)*pow(Sq2,2) - 160*Q2h*pow(l1k2,2)*pow(Sq2,2) +
		        64*Q2k*pow(l1k2,2)*pow(Sq2,2) + 128*pow(l1k2,3)*pow(Sq2,2) -
		        32*m2*pow(Q2e,2)*pow(Sq2,2) + 32*l1k2*pow(Q2h,2)*pow(Sq2,2) -
		        16*m2*pow(Q2h,2)*pow(Sq2,2) + 8*Q2e*pow(Q2h,2)*pow(Sq2,2) -
		        8*pow(Q2h,3)*pow(Sq2,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (16*l1k2*m2*M2*Q2e*Q2h - 16*l1k2*M2*Q2e*Q2h*Q2k -
		        8*l1k2*Q2e*Q2h*Q2k*Sk - 16*l1k2*m2*Q2e*Q2h*Sq2 +
		        12*l1k2*Q2e*Q2h*Q2k*Sq2 - 64*l1k2*m2*Q2e*S*Sq2 +
		        64*l1k2*m2*Q2h*S*Sq2 - 128*l1k2*Q2e*Q2h*S*Sq2 +
		        32*l1k2*Q2e*Q2k*S*Sq2 - 32*l1k2*Q2h*Q2k*S*Sq2 +
		        16*l1k2*Q2e*Q2h*Sk*Sq2 - 16*M2*Q2e*Q2h*pow(l1k2,2) -
		        16*M2*Q2h*Q2k*pow(l1k2,2) + 16*Q2e*Q2h*Sq2*pow(l1k2,2) +
		        16*Q2h*Q2k*Sq2*pow(l1k2,2) + 64*Q2e*S*Sq2*pow(l1k2,2) -
		        128*Q2h*S*Sq2*pow(l1k2,2) + 64*Q2k*S*Sq2*pow(l1k2,2) -
		        32*M2*Q2h*pow(l1k2,3) + 32*Q2h*Sq2*pow(l1k2,3) +
		        128*S*Sq2*pow(l1k2,3) - 4*l1k2*M2*Q2h*pow(Q2e,2) +
		        8*Q2h*Q2k*S*pow(Q2e,2) + 32*Q2h*S*Sk*pow(Q2e,2) +
		        4*l1k2*Q2h*Sq2*pow(Q2e,2) + 4*Q2h*Q2k*Sq2*pow(Q2e,2) +
		        16*l1k2*S*Sq2*pow(Q2e,2) - 16*Q2h*S*Sq2*pow(Q2e,2) +
		        16*Q2h*Sk*Sq2*pow(Q2e,2) - 16*l1k2*m2*M2*pow(Q2h,2) +
		        8*l1k2*m2*Q2e*pow(Q2h,2) + 36*l1k2*M2*Q2e*pow(Q2h,2) +
		        16*l1k2*M2*Q2k*pow(Q2h,2) - 6*l1k2*Q2e*Q2k*pow(Q2h,2) -
		        16*Q2e*Q2k*S*pow(Q2h,2) + 8*l1k2*Q2k*Sk*pow(Q2h,2) -
		        64*Q2e*S*Sk*pow(Q2h,2) + 16*l1k2*m2*Sq2*pow(Q2h,2) -
		        4*l1k2*Q2e*Sq2*pow(Q2h,2) - 12*l1k2*Q2k*Sq2*pow(Q2h,2) -
		        8*Q2e*Q2k*Sq2*pow(Q2h,2) + 112*l1k2*S*Sq2*pow(Q2h,2) +
		        32*Q2e*S*Sq2*pow(Q2h,2) - 16*l1k2*Sk*Sq2*pow(Q2h,2) -
		        32*Q2e*Sk*Sq2*pow(Q2h,2) + 32*M2*pow(l1k2,2)*pow(Q2h,2) -
		        8*Q2e*pow(l1k2,2)*pow(Q2h,2) - 8*Q2k*pow(l1k2,2)*pow(Q2h,2) -
		        32*Sq2*pow(l1k2,2)*pow(Q2h,2) - 16*pow(l1k2,3)*pow(Q2h,2) -
		        2*l1k2*pow(Q2e,2)*pow(Q2h,2) - 8*l1k2*m2*pow(Q2h,3) -
		        32*l1k2*M2*pow(Q2h,3) + 10*l1k2*Q2e*pow(Q2h,3) +
		        6*l1k2*Q2k*pow(Q2h,3) + 8*Q2k*S*pow(Q2h,3) + 32*S*Sk*pow(Q2h,3) +
		        4*Q2k*Sq2*pow(Q2h,3) - 16*S*Sq2*pow(Q2h,3) + 16*Sk*Sq2*pow(Q2h,3) +
		        16*pow(l1k2,2)*pow(Q2h,3) - 8*l1k2*pow(Q2h,4) +
		        4*l1k2*M2*Q2e*pow(Q2k,2) - 4*l1k2*M2*Q2h*pow(Q2k,2) -
		        64*l1k2*m2*Q2e*pow(S,2) + 64*l1k2*m2*Q2h*pow(S,2) -
		        128*l1k2*Q2e*Q2h*pow(S,2) + 32*l1k2*Q2e*Q2k*pow(S,2) -
		        32*l1k2*Q2h*Q2k*pow(S,2) + 64*Q2e*pow(l1k2,2)*pow(S,2) -
		        128*Q2h*pow(l1k2,2)*pow(S,2) + 64*Q2k*pow(l1k2,2)*pow(S,2) +
		        128*pow(l1k2,3)*pow(S,2) + 16*l1k2*pow(Q2e,2)*pow(S,2) +
		        112*l1k2*pow(Q2h,2)*pow(S,2) - 16*l1k2*Q2e*Q2h*pow(Sk,2) +
		        16*l1k2*pow(Q2h,2)*pow(Sk,2) - 32*l1k2*Q2e*Q2h*pow(Sq2,2) -
		        8*Q2h*pow(Q2e,2)*pow(Sq2,2) + 32*l1k2*pow(Q2h,2)*pow(Sq2,2) +
		        16*Q2e*pow(Q2h,2)*pow(Sq2,2) - 8*pow(Q2h,3)*pow(Sq2,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
		      (-16*l1k2*m4*M2*Q2e + 16*l1k2*m4*M2*Q2h - 8*l1k2*m4*Q2e*Q2h -
		        16*l1k2*m2*M2*Q2e*Q2k + 16*l1k2*m2*M2*Q2h*Q2k +
		        10*l1k2*m2*Q2e*Q2h*Q2k - 4*m4*Q2e*Q2h*Q2k - 8*l1k2*m2*Q2e*Q2h*S +
		        16*l1k2*m2*Q2e*Q2k*S - 16*m4*Q2e*Q2k*S - 16*l1k2*m2*Q2h*Q2k*S +
		        16*m4*Q2h*Q2k*S + 8*l1k2*Q2e*Q2h*Q2k*S - 12*m2*Q2e*Q2h*Q2k*S +
		        16*l1k2*m4*Q2e*Sk - 16*l1k2*m4*Q2h*Sk + 20*l1k2*m2*Q2e*Q2h*Sk -
		        8*m4*Q2e*Q2h*Sk - 16*l1k2*m2*Q2e*Q2k*Sk + 16*l1k2*m2*Q2h*Q2k*Sk +
		        30*m2*Q2e*Q2h*Q2k*Sk + 32*l1k2*m2*Q2e*S*Sk - 32*m4*Q2e*S*Sk -
		        32*l1k2*m2*Q2h*S*Sk + 32*m4*Q2h*S*Sk + 16*l1k2*Q2e*Q2h*S*Sk +
		        16*m2*Q2e*Q2k*S*Sk - 16*m2*Q2h*Q2k*S*Sk + 8*Q2e*Q2h*Q2k*S*Sk -
		        8*l1k2*m2*Q2e*Q2h*Sq2 + 8*m4*Q2e*Q2h*Sq2 -
		        26*m2*Q2e*Q2h*Q2k*Sq2 + 32*m4*Q2e*S*Sq2 - 32*m4*Q2h*S*Sq2 +
		        16*l1k2*Q2e*Q2h*S*Sq2 + 32*m2*Q2e*Q2h*S*Sq2 -
		        16*m2*Q2e*Q2h*Sk*Sq2 + 24*m2*Q2e*Q2k*Sk*Sq2 -
		        24*m2*Q2h*Q2k*Sk*Sq2 - 8*m2*M2*Q2e*pow(l1k2,2) +
		        8*m2*M2*Q2h*pow(l1k2,2) - 4*l1k2*m2*M2*pow(Q2e,2) -
		        2*l1k2*M2*Q2h*pow(Q2e,2) + 6*m2*M2*Q2h*pow(Q2e,2) -
		        2*m2*M2*Q2k*pow(Q2e,2) + 4*m2*Q2h*Q2k*pow(Q2e,2) -
		        M2*Q2h*Q2k*pow(Q2e,2) + 4*l1k2*Q2h*S*pow(Q2e,2) -
		        4*l1k2*Q2k*S*pow(Q2e,2) - 8*m2*Q2k*S*pow(Q2e,2) +
		        2*Q2h*Q2k*S*pow(Q2e,2) + 8*m2*Q2h*Sk*pow(Q2e,2) -
		        12*m2*Q2k*Sk*pow(Q2e,2) - 4*Q2h*Q2k*Sk*pow(Q2e,2) -
		        8*l1k2*S*Sk*pow(Q2e,2) - 16*m2*S*Sk*pow(Q2e,2) -
		        8*Q2h*S*Sk*pow(Q2e,2) - 8*Q2k*S*Sk*pow(Q2e,2) -
		        2*Q2h*Q2k*Sq2*pow(Q2e,2) - 12*Q2h*S*Sq2*pow(Q2e,2) +
		        4*Q2k*S*Sq2*pow(Q2e,2) - 4*Q2h*Sk*Sq2*pow(Q2e,2) +
		        M2*Q2h*pow(Q2e,3) - M2*Q2k*pow(Q2e,3) + 8*l1k2*m4*pow(Q2h,2) +
		        4*l1k2*m2*M2*pow(Q2h,2) - 8*l1k2*m2*Q2e*pow(Q2h,2) -
		        12*m2*M2*Q2e*pow(Q2h,2) - 10*l1k2*m2*Q2k*pow(Q2h,2) +
		        4*m4*Q2k*pow(Q2h,2) + 2*m2*M2*Q2k*pow(Q2h,2) -
		        10*m2*Q2e*Q2k*pow(Q2h,2) + 7*M2*Q2e*Q2k*pow(Q2h,2) +
		        8*l1k2*m2*S*pow(Q2h,2) - 4*l1k2*Q2e*S*pow(Q2h,2) +
		        12*m2*Q2e*S*pow(Q2h,2) - 4*l1k2*Q2k*S*pow(Q2h,2) +
		        20*m2*Q2k*S*pow(Q2h,2) + 2*Q2e*Q2k*S*pow(Q2h,2) -
		        20*l1k2*m2*Sk*pow(Q2h,2) + 8*m4*Sk*pow(Q2h,2) -
		        18*m2*Q2e*Sk*pow(Q2h,2) - 18*m2*Q2k*Sk*pow(Q2h,2) +
		        6*Q2e*Q2k*Sk*pow(Q2h,2) - 8*l1k2*S*Sk*pow(Q2h,2) +
		        16*m2*S*Sk*pow(Q2h,2) + 24*Q2e*S*Sk*pow(Q2h,2) +
		        8*l1k2*m2*Sq2*pow(Q2h,2) - 8*m4*Sq2*pow(Q2h,2) +
		        14*m2*Q2e*Sq2*pow(Q2h,2) + 26*m2*Q2k*Sq2*pow(Q2h,2) +
		        4*Q2e*Q2k*Sq2*pow(Q2h,2) - 16*l1k2*S*Sq2*pow(Q2h,2) -
		        32*m2*S*Sq2*pow(Q2h,2) + 24*Q2e*S*Sq2*pow(Q2h,2) -
		        4*Q2k*S*Sq2*pow(Q2h,2) + 16*m2*Sk*Sq2*pow(Q2h,2) +
		        8*Q2e*Sk*Sq2*pow(Q2h,2) - m2*pow(Q2e,2)*pow(Q2h,2) +
		        M2*pow(Q2e,2)*pow(Q2h,2) + Q2k*pow(Q2e,2)*pow(Q2h,2) -
		        2*S*pow(Q2e,2)*pow(Q2h,2) + 2*Sk*pow(Q2e,2)*pow(Q2h,2) +
		        8*l1k2*m2*pow(Q2h,3) + 2*l1k2*M2*pow(Q2h,3) +
		        6*m2*M2*pow(Q2h,3) + m2*Q2e*pow(Q2h,3) - 8*M2*Q2e*pow(Q2h,3) +
		        6*m2*Q2k*pow(Q2h,3) - 5*M2*Q2k*pow(Q2h,3) - 12*m2*S*pow(Q2h,3) +
		        2*Q2e*S*pow(Q2h,3) - 4*Q2k*S*pow(Q2h,3) + 10*m2*Sk*pow(Q2h,3) -
		        2*Q2e*Sk*pow(Q2h,3) - 2*Q2k*Sk*pow(Q2h,3) - 16*S*Sk*pow(Q2h,3) -
		        14*m2*Sq2*pow(Q2h,3) - 2*Q2e*Sq2*pow(Q2h,3) -
		        2*Q2k*Sq2*pow(Q2h,3) - 12*S*Sq2*pow(Q2h,3) - 4*Sk*Sq2*pow(Q2h,3) +
		        6*M2*pow(Q2h,4) - Q2e*pow(Q2h,4) - Q2k*pow(Q2h,4) +
		        2*Sq2*pow(Q2h,4) + pow(Q2h,5) - 6*l1k2*m2*Q2e*pow(Q2k,2) +
		        6*l1k2*m2*Q2h*pow(Q2k,2) + 10*m2*Q2e*Q2h*pow(Q2k,2) -
		        M2*Q2e*Q2h*pow(Q2k,2) + 8*m2*Q2e*S*pow(Q2k,2) -
		        8*m2*Q2h*S*pow(Q2k,2) + 4*Q2e*Q2h*S*pow(Q2k,2) -
		        12*m2*Q2e*Sk*pow(Q2k,2) + 12*m2*Q2h*Sk*pow(Q2k,2) +
		        12*m2*Q2e*Sq2*pow(Q2k,2) - 12*m2*Q2h*Sq2*pow(Q2k,2) -
		        3*m2*pow(Q2e,2)*pow(Q2k,2) - Q2h*pow(Q2e,2)*pow(Q2k,2) -
		        4*S*pow(Q2e,2)*pow(Q2k,2) - 7*m2*pow(Q2h,2)*pow(Q2k,2) +
		        M2*pow(Q2h,2)*pow(Q2k,2) + Q2e*pow(Q2h,2)*pow(Q2k,2) -
		        3*m2*Q2e*pow(Q2k,3) + 3*m2*Q2h*pow(Q2k,3) +
		        32*l1k2*m2*Q2e*pow(S,2) - 32*l1k2*m2*Q2h*pow(S,2) +
		        16*l1k2*Q2e*Q2h*pow(S,2) + 48*m2*Q2e*Q2h*pow(S,2) +
		        16*l1k2*Q2e*Q2k*pow(S,2) - 16*l1k2*Q2h*Q2k*pow(S,2) +
		        16*Q2e*pow(l1k2,2)*pow(S,2) - 16*Q2h*pow(l1k2,2)*pow(S,2) -
		        16*m2*pow(Q2e,2)*pow(S,2) - 20*Q2h*pow(Q2e,2)*pow(S,2) +
		        4*Q2k*pow(Q2e,2)*pow(S,2) - 16*l1k2*pow(Q2h,2)*pow(S,2) -
		        32*m2*pow(Q2h,2)*pow(S,2) + 40*Q2e*pow(Q2h,2)*pow(S,2) -
		        4*Q2k*pow(Q2h,2)*pow(S,2) - 20*pow(Q2h,3)*pow(S,2) -
		        8*l1k2*m2*Q2e*pow(Sk,2) + 8*l1k2*m2*Q2h*pow(Sk,2) +
		        20*m2*Q2e*Q2h*pow(Sk,2) - 12*m2*Q2e*Q2k*pow(Sk,2) +
		        12*m2*Q2h*Q2k*pow(Sk,2) - 12*m2*pow(Q2e,2)*pow(Sk,2) -
		        4*Q2h*pow(Q2e,2)*pow(Sk,2) - 8*m2*pow(Q2h,2)*pow(Sk,2) +
		        8*Q2e*pow(Q2h,2)*pow(Sk,2) - 4*pow(Q2h,3)*pow(Sk,2) +
		        8*m2*Q2e*Q2h*pow(Sq2,2) - 8*m2*pow(Q2h,2)*pow(Sq2,2) +
		        4*Q2e*pow(Q2h,2)*pow(Sq2,2) - 4*pow(Q2h,3)*pow(Sq2,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
		      (64*m4*Q2e*Q2h*S*Sq2 + 8*m4*M2*Q2h*pow(Q2e,2) +
		        4*m2*Q2h*Q2k*Sk*pow(Q2e,2) + 16*Q2h*Q2k*S*Sk*pow(Q2e,2) -
		        8*m4*Q2h*Sq2*pow(Q2e,2) + 2*m2*Q2h*Q2k*Sq2*pow(Q2e,2) -
		        32*m4*S*Sq2*pow(Q2e,2) - 48*m2*Q2h*S*Sq2*pow(Q2e,2) +
		        16*m2*Q2k*S*Sq2*pow(Q2e,2) + 12*Q2h*Q2k*S*Sq2*pow(Q2e,2) -
		        8*m2*Q2h*Sk*Sq2*pow(Q2e,2) + 4*Q2h*Q2k*Sk*Sq2*pow(Q2e,2) +
		        3*M2*Q2h*Q2k*pow(Q2e,3) + 4*Q2h*Q2k*S*pow(Q2e,3) +
		        2*Q2h*Q2k*Sk*pow(Q2e,3) + 16*Q2h*S*Sk*pow(Q2e,3) +
		        12*Q2h*S*Sq2*pow(Q2e,3) - 4*Q2k*S*Sq2*pow(Q2e,3) +
		        4*Q2h*Sk*Sq2*pow(Q2e,3) - 16*m4*M2*Q2e*pow(Q2h,2) -
		        8*m2*M2*Q2e*Q2k*pow(Q2h,2) - 8*m2*Q2e*Q2k*Sk*pow(Q2h,2) -
		        32*Q2e*Q2k*S*Sk*pow(Q2h,2) + 16*m4*Q2e*Sq2*pow(Q2h,2) -
		        20*m2*Q2e*Q2k*Sq2*pow(Q2h,2) - 32*m4*S*Sq2*pow(Q2h,2) +
		        80*m2*Q2e*S*Sq2*pow(Q2h,2) - 16*m2*Q2k*S*Sq2*pow(Q2h,2) -
		        12*Q2e*Q2k*S*Sq2*pow(Q2h,2) + 16*m2*Q2e*Sk*Sq2*pow(Q2h,2) -
		        8*Q2e*Q2k*Sk*Sq2*pow(Q2h,2) + 4*m4*pow(Q2e,2)*pow(Q2h,2) +
		        10*m2*M2*pow(Q2e,2)*pow(Q2h,2) - m2*Q2k*pow(Q2e,2)*pow(Q2h,2) -
		        22*M2*Q2k*pow(Q2e,2)*pow(Q2h,2) - 16*Q2k*S*pow(Q2e,2)*pow(Q2h,2) -
		        8*Q2k*Sk*pow(Q2e,2)*pow(Q2h,2) - 64*S*Sk*pow(Q2e,2)*pow(Q2h,2) -
		        2*m2*Sq2*pow(Q2e,2)*pow(Q2h,2) - 6*Q2k*Sq2*pow(Q2e,2)*pow(Q2h,2) -
		        48*S*Sq2*pow(Q2e,2)*pow(Q2h,2) - 16*Sk*Sq2*pow(Q2e,2)*pow(Q2h,2) -
		        2*M2*pow(Q2e,3)*pow(Q2h,2) + Q2k*pow(Q2e,3)*pow(Q2h,2) +
		        2*Sq2*pow(Q2e,3)*pow(Q2h,2) + 8*m4*M2*pow(Q2h,3) -
		        8*m4*Q2e*pow(Q2h,3) - 16*m2*M2*Q2e*pow(Q2h,3) +
		        8*m2*M2*Q2k*pow(Q2h,3) + 4*m2*Q2e*Q2k*pow(Q2h,3) +
		        35*M2*Q2e*Q2k*pow(Q2h,3) + 20*Q2e*Q2k*S*pow(Q2h,3) +
		        4*m2*Q2k*Sk*pow(Q2h,3) + 10*Q2e*Q2k*Sk*pow(Q2h,3) +
		        80*Q2e*S*Sk*pow(Q2h,3) + 16*Q2k*S*Sk*pow(Q2h,3) -
		        8*m4*Sq2*pow(Q2h,3) + 12*m2*Q2e*Sq2*pow(Q2h,3) +
		        18*m2*Q2k*Sq2*pow(Q2h,3) + 12*Q2e*Q2k*Sq2*pow(Q2h,3) -
		        32*m2*S*Sq2*pow(Q2h,3) + 60*Q2e*S*Sq2*pow(Q2h,3) +
		        4*Q2k*S*Sq2*pow(Q2h,3) - 8*m2*Sk*Sq2*pow(Q2h,3) +
		        20*Q2e*Sk*Sq2*pow(Q2h,3) + 4*Q2k*Sk*Sq2*pow(Q2h,3) +
		        3*m2*pow(Q2e,2)*pow(Q2h,3) + 16*M2*pow(Q2e,2)*pow(Q2h,3) -
		        5*Q2k*pow(Q2e,2)*pow(Q2h,3) - pow(Q2e,3)*pow(Q2h,3) +
		        4*m4*pow(Q2h,4) + 6*m2*M2*pow(Q2h,4) - 7*m2*Q2e*pow(Q2h,4) -
		        26*M2*Q2e*pow(Q2h,4) - 3*m2*Q2k*pow(Q2h,4) -
		        16*M2*Q2k*pow(Q2h,4) + 7*Q2e*Q2k*pow(Q2h,4) - 8*Q2k*S*pow(Q2h,4) -
		        4*Q2k*Sk*pow(Q2h,4) - 32*S*Sk*pow(Q2h,4) - 10*m2*Sq2*pow(Q2h,4) -
		        6*Q2e*Sq2*pow(Q2h,4) - 6*Q2k*Sq2*pow(Q2h,4) - 24*S*Sq2*pow(Q2h,4) -
		        8*Sk*Sq2*pow(Q2h,4) + 4*pow(Q2e,2)*pow(Q2h,4) + 4*m2*pow(Q2h,5) +
		        12*M2*pow(Q2h,5) - 5*Q2e*pow(Q2h,5) - 3*Q2k*pow(Q2h,5) +
		        4*Sq2*pow(Q2h,5) + 2*pow(Q2h,6) + 8*m2*M2*Q2e*Q2h*pow(Q2k,2) +
		        8*m2*Q2e*Q2h*Sq2*pow(Q2k,2) - 16*m2*Q2e*S*Sq2*pow(Q2k,2) +
		        16*m2*Q2h*S*Sq2*pow(Q2k,2) - 8*Q2e*Q2h*S*Sq2*pow(Q2k,2) -
		        2*m2*M2*pow(Q2e,2)*pow(Q2k,2) +
		        9*M2*Q2h*pow(Q2e,2)*pow(Q2k,2) + 4*Q2h*S*pow(Q2e,2)*pow(Q2k,2) +
		        2*Q2h*Sk*pow(Q2e,2)*pow(Q2k,2) + 2*Q2h*Sq2*pow(Q2e,2)*pow(Q2k,2) +
		        4*S*Sq2*pow(Q2e,2)*pow(Q2k,2) - M2*pow(Q2e,3)*pow(Q2k,2) -
		        6*m2*M2*pow(Q2h,2)*pow(Q2k,2) - m2*Q2e*pow(Q2h,2)*pow(Q2k,2) -
		        15*M2*Q2e*pow(Q2h,2)*pow(Q2k,2) - 8*Q2e*S*pow(Q2h,2)*pow(Q2k,2) -
		        4*Q2e*Sk*pow(Q2h,2)*pow(Q2k,2) - 8*m2*Sq2*pow(Q2h,2)*pow(Q2k,2) -
		        4*Q2e*Sq2*pow(Q2h,2)*pow(Q2k,2) + 4*S*Sq2*pow(Q2h,2)*pow(Q2k,2) +
		        pow(Q2e,2)*pow(Q2h,2)*pow(Q2k,2) + m2*pow(Q2h,3)*pow(Q2k,2) +
		        7*M2*pow(Q2h,3)*pow(Q2k,2) - 2*Q2e*pow(Q2h,3)*pow(Q2k,2) +
		        4*S*pow(Q2h,3)*pow(Q2k,2) + 2*Sk*pow(Q2h,3)*pow(Q2k,2) +
		        2*Sq2*pow(Q2h,3)*pow(Q2k,2) + pow(Q2h,4)*pow(Q2k,2) +
		        2*M2*Q2e*Q2h*pow(Q2k,3) - M2*pow(Q2e,2)*pow(Q2k,3) -
		        M2*pow(Q2h,2)*pow(Q2k,3) + 64*m4*Q2e*Q2h*pow(S,2) -
		        32*m4*pow(Q2e,2)*pow(S,2) - 48*m2*Q2h*pow(Q2e,2)*pow(S,2) +
		        16*m2*Q2k*pow(Q2e,2)*pow(S,2) + 20*Q2h*Q2k*pow(Q2e,2)*pow(S,2) +
		        20*Q2h*pow(Q2e,3)*pow(S,2) - 4*Q2k*pow(Q2e,3)*pow(S,2) -
		        32*m4*pow(Q2h,2)*pow(S,2) + 80*m2*Q2e*pow(Q2h,2)*pow(S,2) -
		        16*m2*Q2k*pow(Q2h,2)*pow(S,2) - 28*Q2e*Q2k*pow(Q2h,2)*pow(S,2) -
		        80*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) - 32*m2*pow(Q2h,3)*pow(S,2) +
		        100*Q2e*pow(Q2h,3)*pow(S,2) + 12*Q2k*pow(Q2h,3)*pow(S,2) -
		        40*pow(Q2h,4)*pow(S,2) - 16*m2*Q2e*pow(Q2k,2)*pow(S,2) +
		        16*m2*Q2h*pow(Q2k,2)*pow(S,2) - 8*Q2e*Q2h*pow(Q2k,2)*pow(S,2) +
		        4*pow(Q2e,2)*pow(Q2k,2)*pow(S,2) +
		        4*pow(Q2h,2)*pow(Q2k,2)*pow(S,2) + 8*m2*Q2h*pow(Q2e,2)*pow(Sk,2) +
		        4*Q2h*Q2k*pow(Q2e,2)*pow(Sk,2) + 4*Q2h*pow(Q2e,3)*pow(Sk,2) -
		        16*m2*Q2e*pow(Q2h,2)*pow(Sk,2) - 8*Q2e*Q2k*pow(Q2h,2)*pow(Sk,2) -
		        16*pow(Q2e,2)*pow(Q2h,2)*pow(Sk,2) + 8*m2*pow(Q2h,3)*pow(Sk,2) +
		        20*Q2e*pow(Q2h,3)*pow(Sk,2) + 4*Q2k*pow(Q2h,3)*pow(Sk,2) -
		        8*pow(Q2h,4)*pow(Sk,2) + 24*m2*Q2e*Q2h*Q2k*pow(Sq2,2) -
		        8*m2*Q2h*pow(Q2e,2)*pow(Sq2,2) + 4*Q2h*Q2k*pow(Q2e,2)*pow(Sq2,2) +
		        4*m2*Q2e*pow(Q2h,2)*pow(Sq2,2) - 24*m2*Q2k*pow(Q2h,2)*pow(Sq2,2) -
		        8*Q2e*Q2k*pow(Q2h,2)*pow(Sq2,2) -
		        8*pow(Q2e,2)*pow(Q2h,2)*pow(Sq2,2) + 4*m2*pow(Q2h,3)*pow(Sq2,2) +
		        16*Q2e*pow(Q2h,3)*pow(Sq2,2) + 4*Q2k*pow(Q2h,3)*pow(Sq2,2) -
		        8*pow(Q2h,4)*pow(Sq2,2) - 12*m2*Q2e*pow(Q2k,2)*pow(Sq2,2) +
		        12*m2*Q2h*pow(Q2k,2)*pow(Sq2,2)) +
		     pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
		      (-32*l1k2*m4*M2*Q2e*Q2h + 16*l1k2*m2*M2*Q2e*Q2h*Q2k -
		        64*l1k2*m4*Q2e*Q2h*S + 32*l1k2*m2*Q2e*Q2h*Q2k*S -
		        64*l1k2*m4*Q2e*Q2h*Sq2 + 32*l1k2*m2*Q2e*Q2h*Q2k*Sq2 +
		        256*l1k2*m4*Q2e*S*Sq2 - 256*l1k2*m4*Q2h*S*Sq2 +
		        128*l1k2*m2*Q2e*Q2h*S*Sq2 - 128*l1k2*m2*Q2e*Q2k*S*Sq2 +
		        128*l1k2*m2*Q2h*Q2k*S*Sq2 + 64*l1k2*Q2e*Q2h*Q2k*S*Sq2 -
		        64*m4*M2*Q2h*pow(l1k2,2) + 32*m2*M2*Q2e*Q2h*pow(l1k2,2) +
		        32*m2*M2*Q2h*Q2k*pow(l1k2,2) + 16*M2*Q2e*Q2h*Q2k*pow(l1k2,2) -
		        128*m4*Q2h*S*pow(l1k2,2) + 64*m2*Q2e*Q2h*S*pow(l1k2,2) +
		        64*m2*Q2h*Q2k*S*pow(l1k2,2) + 32*Q2e*Q2h*Q2k*S*pow(l1k2,2) -
		        128*m4*Q2h*Sq2*pow(l1k2,2) + 64*m2*Q2e*Q2h*Sq2*pow(l1k2,2) +
		        64*m2*Q2h*Q2k*Sq2*pow(l1k2,2) + 32*Q2e*Q2h*Q2k*Sq2*pow(l1k2,2) +
		        512*m4*S*Sq2*pow(l1k2,2) - 256*m2*Q2e*S*Sq2*pow(l1k2,2) +
		        512*m2*Q2h*S*Sq2*pow(l1k2,2) + 128*Q2e*Q2h*S*Sq2*pow(l1k2,2) -
		        256*m2*Q2k*S*Sq2*pow(l1k2,2) - 128*Q2e*Q2k*S*Sq2*pow(l1k2,2) +
		        256*Q2h*Q2k*S*Sq2*pow(l1k2,2) + 64*m2*M2*Q2h*pow(l1k2,3) +
		        16*M2*Q2e*Q2h*pow(l1k2,3) + 32*M2*Q2h*Q2k*pow(l1k2,3) +
		        128*m2*Q2h*S*pow(l1k2,3) + 32*Q2e*Q2h*S*pow(l1k2,3) +
		        64*Q2h*Q2k*S*pow(l1k2,3) + 128*m2*Q2h*Sq2*pow(l1k2,3) +
		        32*Q2e*Q2h*Sq2*pow(l1k2,3) + 64*Q2h*Q2k*Sq2*pow(l1k2,3) -
		        512*m2*S*Sq2*pow(l1k2,3) - 128*Q2e*S*Sq2*pow(l1k2,3) +
		        384*Q2h*S*Sq2*pow(l1k2,3) - 256*Q2k*S*Sq2*pow(l1k2,3) +
		        32*M2*Q2h*pow(l1k2,4) + 64*Q2h*S*pow(l1k2,4) +
		        64*Q2h*Sq2*pow(l1k2,4) - 256*S*Sq2*pow(l1k2,4) +
		        32*l1k2*m4*M2*pow(Q2h,2) - 16*l1k2*m2*M2*Q2e*pow(Q2h,2) -
		        16*l1k2*m2*M2*Q2k*pow(Q2h,2) - 8*l1k2*M2*Q2e*Q2k*pow(Q2h,2) +
		        64*l1k2*m4*S*pow(Q2h,2) - 32*l1k2*m2*Q2e*S*pow(Q2h,2) -
		        32*l1k2*m2*Q2k*S*pow(Q2h,2) - 16*l1k2*Q2e*Q2k*S*pow(Q2h,2) +
		        64*l1k2*m4*Sq2*pow(Q2h,2) - 32*l1k2*m2*Q2e*Sq2*pow(Q2h,2) -
		        32*l1k2*m2*Q2k*Sq2*pow(Q2h,2) - 16*l1k2*Q2e*Q2k*Sq2*pow(Q2h,2) -
		        128*l1k2*m2*S*Sq2*pow(Q2h,2) - 32*l1k2*Q2e*S*Sq2*pow(Q2h,2) -
		        64*l1k2*Q2k*S*Sq2*pow(Q2h,2) - 64*m2*M2*pow(l1k2,2)*pow(Q2h,2) -
		        16*M2*Q2e*pow(l1k2,2)*pow(Q2h,2) -
		        32*M2*Q2k*pow(l1k2,2)*pow(Q2h,2) -
		        128*m2*S*pow(l1k2,2)*pow(Q2h,2) - 32*Q2e*S*pow(l1k2,2)*pow(Q2h,2) -
		        64*Q2k*S*pow(l1k2,2)*pow(Q2h,2) -
		        128*m2*Sq2*pow(l1k2,2)*pow(Q2h,2) -
		        32*Q2e*Sq2*pow(l1k2,2)*pow(Q2h,2) -
		        64*Q2k*Sq2*pow(l1k2,2)*pow(Q2h,2) -
		        192*S*Sq2*pow(l1k2,2)*pow(Q2h,2) - 48*M2*pow(l1k2,3)*pow(Q2h,2) -
		        96*S*pow(l1k2,3)*pow(Q2h,2) - 96*Sq2*pow(l1k2,3)*pow(Q2h,2) +
		        16*l1k2*m2*M2*pow(Q2h,3) + 4*l1k2*M2*Q2e*pow(Q2h,3) +
		        8*l1k2*M2*Q2k*pow(Q2h,3) + 32*l1k2*m2*S*pow(Q2h,3) +
		        8*l1k2*Q2e*S*pow(Q2h,3) + 16*l1k2*Q2k*S*pow(Q2h,3) +
		        32*l1k2*m2*Sq2*pow(Q2h,3) + 8*l1k2*Q2e*Sq2*pow(Q2h,3) +
		        16*l1k2*Q2k*Sq2*pow(Q2h,3) + 32*l1k2*S*Sq2*pow(Q2h,3) +
		        24*M2*pow(l1k2,2)*pow(Q2h,3) + 48*S*pow(l1k2,2)*pow(Q2h,3) +
		        48*Sq2*pow(l1k2,2)*pow(Q2h,3) - 4*l1k2*M2*pow(Q2h,4) -
		        8*l1k2*S*pow(Q2h,4) - 8*l1k2*Sq2*pow(Q2h,4) +
		        4*l1k2*M2*Q2e*Q2h*pow(Q2k,2) + 8*l1k2*Q2e*Q2h*S*pow(Q2k,2) +
		        8*l1k2*Q2e*Q2h*Sq2*pow(Q2k,2) - 32*l1k2*Q2e*S*Sq2*pow(Q2k,2) +
		        32*l1k2*Q2h*S*Sq2*pow(Q2k,2) + 8*M2*Q2h*pow(l1k2,2)*pow(Q2k,2) +
		        16*Q2h*S*pow(l1k2,2)*pow(Q2k,2) + 16*Q2h*Sq2*pow(l1k2,2)*pow(Q2k,2) -
		        64*S*Sq2*pow(l1k2,2)*pow(Q2k,2) - 4*l1k2*M2*pow(Q2h,2)*pow(Q2k,2) -
		        8*l1k2*S*pow(Q2h,2)*pow(Q2k,2) - 8*l1k2*Sq2*pow(Q2h,2)*pow(Q2k,2) +
		        128*l1k2*m4*Q2e*pow(S,2) - 128*l1k2*m4*Q2h*pow(S,2) +
		        64*l1k2*m2*Q2e*Q2h*pow(S,2) - 64*l1k2*m2*Q2e*Q2k*pow(S,2) +
		        64*l1k2*m2*Q2h*Q2k*pow(S,2) + 32*l1k2*Q2e*Q2h*Q2k*pow(S,2) +
		        256*m4*pow(l1k2,2)*pow(S,2) - 128*m2*Q2e*pow(l1k2,2)*pow(S,2) +
		        256*m2*Q2h*pow(l1k2,2)*pow(S,2) + 64*Q2e*Q2h*pow(l1k2,2)*pow(S,2) -
		        128*m2*Q2k*pow(l1k2,2)*pow(S,2) - 64*Q2e*Q2k*pow(l1k2,2)*pow(S,2) +
		        128*Q2h*Q2k*pow(l1k2,2)*pow(S,2) - 256*m2*pow(l1k2,3)*pow(S,2) -
		        64*Q2e*pow(l1k2,3)*pow(S,2) + 192*Q2h*pow(l1k2,3)*pow(S,2) -
		        128*Q2k*pow(l1k2,3)*pow(S,2) - 128*pow(l1k2,4)*pow(S,2) -
		        64*l1k2*m2*pow(Q2h,2)*pow(S,2) - 16*l1k2*Q2e*pow(Q2h,2)*pow(S,2) -
		        32*l1k2*Q2k*pow(Q2h,2)*pow(S,2) -
		        96*pow(l1k2,2)*pow(Q2h,2)*pow(S,2) + 16*l1k2*pow(Q2h,3)*pow(S,2) -
		        16*l1k2*Q2e*pow(Q2k,2)*pow(S,2) + 16*l1k2*Q2h*pow(Q2k,2)*pow(S,2) -
		        32*pow(l1k2,2)*pow(Q2k,2)*pow(S,2) + 128*l1k2*m4*Q2e*pow(Sq2,2) -
		        128*l1k2*m4*Q2h*pow(Sq2,2) + 64*l1k2*m2*Q2e*Q2h*pow(Sq2,2) -
		        64*l1k2*m2*Q2e*Q2k*pow(Sq2,2) + 64*l1k2*m2*Q2h*Q2k*pow(Sq2,2) +
		        32*l1k2*Q2e*Q2h*Q2k*pow(Sq2,2) + 256*m4*pow(l1k2,2)*pow(Sq2,2) -
		        128*m2*Q2e*pow(l1k2,2)*pow(Sq2,2) +
		        256*m2*Q2h*pow(l1k2,2)*pow(Sq2,2) +
		        64*Q2e*Q2h*pow(l1k2,2)*pow(Sq2,2) -
		        128*m2*Q2k*pow(l1k2,2)*pow(Sq2,2) -
		        64*Q2e*Q2k*pow(l1k2,2)*pow(Sq2,2) +
		        128*Q2h*Q2k*pow(l1k2,2)*pow(Sq2,2) - 256*m2*pow(l1k2,3)*pow(Sq2,2) -
		        64*Q2e*pow(l1k2,3)*pow(Sq2,2) + 192*Q2h*pow(l1k2,3)*pow(Sq2,2) -
		        128*Q2k*pow(l1k2,3)*pow(Sq2,2) - 128*pow(l1k2,4)*pow(Sq2,2) -
		        64*l1k2*m2*pow(Q2h,2)*pow(Sq2,2) -
		        16*l1k2*Q2e*pow(Q2h,2)*pow(Sq2,2) -
		        32*l1k2*Q2k*pow(Q2h,2)*pow(Sq2,2) -
		        96*pow(l1k2,2)*pow(Q2h,2)*pow(Sq2,2) +
		        16*l1k2*pow(Q2h,3)*pow(Sq2,2) - 16*l1k2*Q2e*pow(Q2k,2)*pow(Sq2,2) +
		        16*l1k2*Q2h*pow(Q2k,2)*pow(Sq2,2) -
		        32*pow(l1k2,2)*pow(Q2k,2)*pow(Sq2,2)));
}

long double Melem::melem2_test(const long double l1k1,const long double l1k2,const long double k1k2,const long double Q2e,
		const long double Q2h,const long double Q2k,const long double S,const long double Sk,const long double Sq2,
		const long double f1, const long double f2)const {

	long double p1p2 = M2 + Q2h/2.;
	long double G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
	long double G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
	long double G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;

	return G1*(16*m2*pow(l1k1,-2) + 8*m2*pow(l1k2,-2) -
		     4*m2*(2*l1k1 - 4*m2 + Q2e + 2*Q2h - Q2k)*pow(k1k2 - l1k1 - l1k2,-1)*
		      pow(l1k2,-2) + (8*m6 - 4*m4*Q2h)*pow(l1k1,-2)*
		      pow(-k1k2 + l1k1 + l1k2,-2) +
		     (8*m6 - 4*m4*Q2h)*pow(l1k2,-2)*pow(-k1k2 + l1k1 + l1k2,-2) -
		     4*m2*(4*m2 + Q2h + Q2k)*pow(l1k1,-2)*pow(-k1k2 + l1k1 + l1k2,-1) +
		     64*m2*pow(2*l1k1 + Q2e - Q2k,-2) +
		     16*(2*m6 - m4*Q2h)*pow(l1k2,-2)*pow(2*l1k1 + Q2e - Q2k,-2) -
		     16*m2*(4*m2 + Q2h + Q2k)*pow(l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-2) -
		     64*(-2*m6 + m4*Q2h)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		      pow(2*l1k1 + Q2e - Q2k,-2) +
		     32*m2*(4*m2 + Q2h + Q2k)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      pow(2*l1k1 + Q2e - Q2k,-2) +
		     8*m2*(-2*k1k2 + 4*m2 - Q2e - 2*Q2h + Q2k)*pow(l1k2,-2)*
		      pow(2*l1k1 + Q2e - Q2k,-1) +
		     8*pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-2)*(4*m6 - m2*pow(Q2h,2))*
		      pow(2*l1k1 + Q2e - Q2k,-1) +
		     64*m2*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) +
		     16*(2*m6 - m4*Q2h)*pow(l1k1,-2)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) +
		     16*m2*(2*l1k2 - 4*m2 + Q2e + Q2h)*pow(l1k1,-1)*
		      pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) -
		     64*(-2*m6 + m4*Q2h)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*
		      pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) -
		     32*m2*(2*l1k2 - 4*m2 + Q2e + Q2h)*
		      pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) +
		     32*pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      (-4*m6 + m2*pow(Q2h,2))*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) -
		     8*m2*(4*m2 + Q2h + Q2k)*pow(l1k1,-2)*
		      pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1) +
		     32*pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      pow(2*l1k1 + Q2e - Q2k,-2)*(-4*m6 + m2*pow(Q2k,2)) -
		     8*pow(l1k1,-2)*pow(-k1k2 + l1k1 + l1k2,-1)*
		      pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*(-4*m6 + m2*pow(Q2k,2))) +
		  (G2 + G3)*(16*m2*M2*pow(l1k1,-2) + 8*m2*M2*pow(l1k2,-2) -
		     4*m2*(2*l1k1*M2 + M2*(Q2e + 2*Q2h - Q2k) +
		        Q2h*(Q2h - Q2k + 2*(S - Sk + Sq2)))*pow(k1k2 - l1k1 - l1k2,-1)*
		      pow(l1k2,-2) + 4*m4*(-(M2*Q2h) + 2*S*(Q2h + 2*S))*pow(l1k1,-2)*
		      pow(-k1k2 + l1k1 + l1k2,-2) +
		     4*m4*(-(M2*Q2h) + 2*S*(Q2h + 2*S))*pow(l1k2,-2)*
		      pow(-k1k2 + l1k1 + l1k2,-2) + 64*m2*M2*pow(2*l1k1 + Q2e - Q2k,-2) +
		     16*m4*(-(M2*Q2h) + (Q2k - 2*(S - Sk + Sq2))*
		         (-Q2h + Q2k - 2*(S - Sk + Sq2)))*pow(l1k2,-2)*
		      pow(2*l1k1 + Q2e - Q2k,-2) -
		     64*m4*(M2*Q2h + 2*(S + Sq2)*(Q2h - 2*(S + Sq2)))*
		      pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*pow(2*l1k1 + Q2e - Q2k,-2) -
		     8*m2*(2*k1k2*M2 + M2*(Q2e + 2*Q2h - Q2k) - 2*Q2h*S)*pow(l1k2,-2)*
		      pow(2*l1k1 + Q2e - Q2k,-1) +
		     64*m2*M2*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) -
		     16*m4*(M2*Q2h + (Q2h - Q2k - 2*(S + Sk))*(Q2k + 2*(S + Sk)))*
		      pow(l1k1,-2)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) +
		     16*m2*(2*l1k2*M2 + M2*(Q2e + Q2h) + 2*Q2h*(S + Sq2))*pow(l1k1,-1)*
		      pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) -
		     64*m4*(M2*Q2h + 2*(S + Sq2)*(Q2h - 2*(S + Sq2)))*
		      pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-2)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) \
		- 32*m2*(2*l1k2*M2 + M2*(Q2e + Q2h) + Q2h*(Q2h - Q2k - 2*(S + Sk)))*
		      pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2) \
		- 8*m2*pow(l1k1,-2)*pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
		      (Q2h*Q2k + M2*(Q2h + Q2k) + 2*Q2h*S - 4*Q2k*S +
		        2*m2*(2*M2 + Q2h - 2*Sk) + 2*Q2h*Sk - 2*Q2k*Sk - 8*S*Sk -
		        8*pow(S,2) - 4*pow(Sk,2)) -
		     4*m2*pow(l1k1,-2)*pow(-k1k2 + l1k1 + l1k2,-1)*
		      (M2*(Q2h + Q2k) + 2*m2*(2*M2 + Q2h - 2*Sk) -
		        2*(Q2h*S + 2*Q2k*S + Q2k*Sk + 4*S*Sk + 4*pow(S,2) + 2*pow(Sk,2))) -
		     8*pow(l1k1,-2)*pow(-k1k2 + l1k1 + l1k2,-1)*
		      pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-1)*
		      (2*m4*(M2*Q2k - 4*S*(S + Sk) - Q2k*(2*S + Sk)) +
		        m2*Q2k*(M2*Q2k - Q2k*(2*S + Sk) -
		           2*(2*S*Sk + 2*pow(S,2) + pow(Sk,2)))) +
		     16*pow(l1k1,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      pow(2*k1k2 - 2*l1k2 + Q2h - Q2k,-2)*
		      (2*m4*(2*M2*Q2h - 4*(Q2k + 2*(S + Sk))*(S + Sq2) +
		           Q2h*(Q2k + 2*(2*S + Sk + Sq2))) +
		        m2*Q2h*(2*M2*Q2h - 4*Q2k*S - 4*Q2k*Sk - 8*S*Sk - 8*S*Sq2 +
		           Q2h*(Q2k + 2*(2*S + Sk + Sq2)) - pow(Q2k,2) - 8*pow(S,2) -
		           4*pow(Sk,2) - 4*pow(Sq2,2))) -
		     16*m2*pow(l1k2,-1)*pow(2*l1k1 + Q2e - Q2k,-2)*
		      (Q2h*Q2k + M2*(Q2h + Q2k) - 2*Q2h*S + 4*Q2k*S +
		        2*m2*(2*M2 + Q2h - 2*Sk) + 2*Q2h*Sk - 2*Q2k*Sk + 8*S*Sk -
		        2*Q2h*Sq2 + 4*Q2k*Sq2 - 16*S*Sq2 + 8*Sk*Sq2 - 8*pow(S,2) -
		        4*pow(Sk,2) - 8*pow(Sq2,2)) +
		     32*m2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*pow(2*l1k1 + Q2e - Q2k,-2)*
		      (M2*(Q2h + Q2k) + 4*Q2k*S + 2*m2*(2*M2 + Q2h - 2*Sk) - 2*Q2k*Sk +
		        8*S*Sk + 4*Q2k*Sq2 - 16*S*Sq2 + 8*Sk*Sq2 + 2*Q2h*(S + Sq2) -
		        8*pow(S,2) - 4*pow(Sk,2) - 8*pow(Sq2,2)) +
		     pow(k1k2 - l1k1 - l1k2,-1)*pow(l1k2,-2)*pow(2*l1k1 + Q2e - Q2k,-1)*
		      (-8*m4*(2*M2*Q2h + Q2h*(Q2k - 4*S + 2*Sk - 2*Sq2) +
		           4*S*(Q2k - 2*(S - Sk + Sq2))) +
		        4*m2*Q2h*(-2*M2*Q2h - 4*Q2k*S + 4*Q2k*Sk - 8*S*Sk -
		           Q2h*(Q2k - 4*S + 2*Sk - 2*Sq2) - 4*Q2k*Sq2 + 8*S*Sq2 - 8*Sk*Sq2 +
		           pow(Q2k,2) + 8*pow(S,2) + 4*pow(Sk,2) + 4*pow(Sq2,2))) +
		     32*pow(l1k2,-1)*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		      pow(2*l1k1 + Q2e - Q2k,-2)*
		      (2*m4*(M2*Q2k - 4*(S + Sq2)*(S - Sk + Sq2) +
		           Q2k*(2*S - Sk + 2*Sq2)) +
		        m2*Q2k*(M2*Q2k + Q2k*(2*S - Sk + 2*Sq2) -
		           2*(-2*S*Sk + 4*S*Sq2 - 2*Sk*Sq2 + 2*pow(S,2) + pow(Sk,2) +
		              2*pow(Sq2,2)))));
}

long double Melem::melem2_l1k1_2(const long double l1k1,const long double l1k2,const long double k1k2,const long double Q2e,
		const long double Q2h,const long double Q2k,const long double S,const long double Sk,const long double Sq2,
		const long double f1, const long double f2)const {

	long double p1p2 = M2 + Q2h/2.;
	long double G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
	long double G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
	long double G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;
	long double l2k1 = l1k1 - Q2k/2. + Q2e/2.;
	long double l2k2 = l1k2 - k1k2 + Q2k/2 - Q2h/2;

	return m2*pow(l1k1,-2)*(l2k2*pow(k1k2 - l1k1 - l1k2,-1)*
		     (-4*G1*(4*l1k2*l2k2 + l1k2*Q2h - l2k2*Q2h + l1k2*Q2k - l2k2*Q2k -
		          k1k2*(4*l2k2 + Q2h + Q2k) + l1k1*(4*l2k2 + Q2h + Q2k) +
		          pow(Q2k,2)) + 4*(G2 + G3)*
		        (-4*l1k2*l2k2*M2 - l1k2*M2*Q2h + l2k2*M2*Q2h -
		          l1k2*M2*Q2k + l2k2*M2*Q2k - l1k2*Q2h*Q2k - 2*l1k2*Q2h*S -
		          2*l2k2*Q2h*S + 4*l1k2*Q2k*S - 4*l2k2*Q2k*S - 2*l1k2*Q2h*Sk +
		          2*l1k2*Q2k*Sk - 2*l2k2*Q2k*Sk + 8*l1k2*S*Sk - 8*l2k2*S*Sk +
		          4*Q2k*S*Sk - M2*pow(Q2k,2) + 2*S*pow(Q2k,2) + Sk*pow(Q2k,2) +
		          8*l1k2*pow(S,2) - 8*l2k2*pow(S,2) + 4*Q2k*pow(S,2) +
		          k1k2*(4*l2k2*M2 + Q2h*Q2k + M2*(Q2h + Q2k) + 2*Q2h*S -
		             4*Q2k*S + 2*Q2h*Sk - 2*Q2k*Sk - 8*S*Sk - 8*pow(S,2) -
		             4*pow(Sk,2)) - l1k1*
		           (4*l2k2*M2 + Q2h*Q2k + M2*(Q2h + Q2k) + 2*Q2h*S - 4*Q2k*S +
		             2*Q2h*Sk - 2*Q2k*Sk - 8*S*Sk - 8*pow(S,2) - 4*pow(Sk,2)) +
		          4*l1k2*pow(Sk,2) - 4*l2k2*pow(Sk,2) + 2*Q2k*pow(Sk,2))) +
		    4*l1k1*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (G1*(4*l1k2*Q2e + 4*l1k1*(2*l1k2 + Q2e - Q2h) + 8*pow(l1k1,2) +
		          4*pow(l1k2,2) + pow(Q2e,2) + pow(Q2h,2)) +
		       (G2 + G3)*(4*l1k2*M2*Q2e + 4*l1k2*Q2h*S + 2*Q2e*Q2h*S -
		          4*Q2h*Q2k*S - 4*Q2h*Q2k*Sk - 8*Q2h*S*Sk + 4*l1k2*Q2h*Sq2 +
		          2*Q2e*Q2h*Sq2 - 8*Q2h*S*Sq2 +
		          2*l1k1*(4*l1k2*M2 + 2*M2*(Q2e - Q2h) +
		             Q2h*(-Q2h + Q2k + 2*(2*S + Sk + Sq2))) + 8*M2*pow(l1k1,2) +
		          4*M2*pow(l1k2,2) + M2*pow(Q2e,2) + M2*pow(Q2h,2) +
		          Q2k*pow(Q2h,2) + 2*S*pow(Q2h,2) + 2*Sk*pow(Q2h,2) -
		          Q2h*pow(Q2k,2) - 8*Q2h*pow(S,2) - 4*Q2h*pow(Sk,2) -
		          4*Q2h*pow(Sq2,2))))*pow(pow(l1k1,2) + pow(l2k2,2),-1);
}

long double Melem::melem2_l1k2_2(const long double l1k1,const long double l1k2,const long double k1k2,const long double Q2e,
		const long double Q2h,const long double Q2k,const long double S,const long double Sk,const long double Sq2,
		const long double f1, const long double f2)const {

	long double p1p2 = M2 + Q2h/2.;
	long double G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
	long double G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
	long double G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;
	long double l2k1 = l1k1 - Q2k/2. + Q2e/2.;
	long double l2k2 = l1k2 - k1k2 + Q2k/2 - Q2h/2;

	return 2*m2*pow(l1k2,-2)*(-(l2k1*pow(k1k2 - l1k1 - l1k2,-1)*
		       (2*G1*(2*l1k2*l2k1 - l1k2*Q2e + l2k1*Q2e - 2*l1k2*Q2h +
		            2*l2k1*Q2h + l1k2*Q2k - l2k1*Q2k -
		            k1k2*(2*l1k1 + 2*l1k2 + 2*l2k1 - Q2e - 2*Q2h + Q2k) +
		            l1k1*(4*l2k1 - Q2e - 2*Q2h + Q2k) + 2*pow(k1k2,2) + pow(Q2h,2)) \
		+ (G2 + G3)*(4*l1k2*l2k1*M2 - 2*l1k2*M2*Q2e + 2*l2k1*M2*Q2e -
		            4*l1k2*M2*Q2h + 4*l2k1*M2*Q2h + 2*l1k2*M2*Q2k -
		            2*l2k1*M2*Q2k - 2*l2k1*Q2h*Q2k + 4*l1k2*Q2h*S +
		            4*l2k1*Q2h*S + 4*Q2h*Q2k*S -
		            2*k1k2*(2*l1k1*M2 + 2*l1k2*M2 + 2*l2k1*M2 - M2*Q2e -
		               2*M2*Q2h + M2*Q2k + 2*Q2h*S) +
		            2*l1k1*(4*l2k1*M2 + M2*(-Q2e - 2*Q2h + Q2k) + 2*Q2h*S) -
		            4*l2k1*Q2h*Sk - 4*Q2h*Q2k*Sk + 8*Q2h*S*Sk + 4*l2k1*Q2h*Sq2 +
		            4*Q2h*Q2k*Sq2 - 8*Q2h*S*Sq2 + 8*Q2h*Sk*Sq2 +
		            4*M2*pow(k1k2,2) + 2*l2k1*pow(Q2h,2) + 2*M2*pow(Q2h,2) +
		            Q2k*pow(Q2h,2) - 4*S*pow(Q2h,2) + 2*Sk*pow(Q2h,2) -
		            2*Sq2*pow(Q2h,2) - Q2h*pow(Q2k,2) - 8*Q2h*pow(S,2) -
		            4*Q2h*pow(Sk,2) - 4*Q2h*pow(Sq2,2)))) +
		    2*l1k2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (G1*(8*l1k1*l1k2 + 4*l1k2*(Q2e - Q2h) - Q2e*Q2h - Q2e*Q2k + Q2h*Q2k -
		          2*l1k1*(Q2h + Q2k) + 8*pow(l1k2,2) + pow(Q2h,2) + 2*pow(Q2k,2)) +
		       (G2 + G3)*(4*l1k2*M2*(Q2e - Q2h) - M2*Q2e*Q2h - M2*Q2e*Q2k +
		          M2*Q2h*Q2k - Q2e*Q2h*Q2k + 2*Q2e*Q2h*S - 4*Q2e*Q2k*S +
		          4*Q2h*Q2k*S - 2*Q2e*Q2h*Sk + 2*Q2e*Q2k*Sk - 2*Q2h*Q2k*Sk -
		          8*Q2e*S*Sk + 8*Q2h*S*Sk + 8*Q2k*S*Sk -
		          2*l1k2*Q2h*(Q2k - 4*S + 2*Sk - 4*Sq2) + 2*Q2e*Q2h*Sq2 -
		          4*Q2e*Q2k*Sq2 + 4*Q2h*Q2k*Sq2 + 16*Q2e*S*Sq2 - 16*Q2h*S*Sq2 -
		          16*Q2k*S*Sq2 - 8*Q2e*Sk*Sq2 + 8*Q2h*Sk*Sq2 + 8*Q2k*Sk*Sq2 +
		          8*M2*pow(l1k2,2) + M2*pow(Q2h,2) + Q2k*pow(Q2h,2) -
		          2*S*pow(Q2h,2) + 2*Sk*pow(Q2h,2) - 2*Sq2*pow(Q2h,2) +
		          2*M2*pow(Q2k,2) + 4*S*pow(Q2k,2) - 2*Sk*pow(Q2k,2) +
		          4*Sq2*pow(Q2k,2) + 8*Q2e*pow(S,2) - 8*Q2h*pow(S,2) -
		          8*Q2k*pow(S,2) + 4*Q2e*pow(Sk,2) - 4*Q2h*pow(Sk,2) -
		          4*Q2k*pow(Sk,2) + 8*Q2e*pow(Sq2,2) - 8*Q2h*pow(Sq2,2) -
		          8*Q2k*pow(Sq2,2) + 2*l1k1*
		           (4*l1k2*M2 - Q2h*Q2k - M2*(Q2h + Q2k) + 2*Q2h*S - 4*Q2k*S -
		             2*Q2h*Sk + 2*Q2k*Sk - 8*S*Sk + 2*Q2h*Sq2 - 4*Q2k*Sq2 +
		             16*S*Sq2 - 8*Sk*Sq2 + 8*pow(S,2) + 4*pow(Sk,2) + 8*pow(Sq2,2))))\
		)*pow(pow(l1k2,2) + pow(l2k1,2),-1);
}

long double Melem::melem2_l2k1_2(const long double l1k1,const long double l1k2,const long double k1k2,const long double Q2e,
		const long double Q2h,const long double Q2k,const long double S,const long double Sk,const long double Sq2,
		const long double f1, const long double f2)const {

	long double p1p2 = M2 + Q2h/2.;
	long double G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
	long double G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
	long double G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;
	long double l2k1 = l1k1 - Q2k/2. + Q2e/2.;
	long double l2k2 = l1k2 - k1k2 + Q2k/2 - Q2h/2;

	return 2*m2*pow(l2k1,-2)*(-(l2k1*pow(k1k2 - l1k1 - l1k2,-1)*
		       (2*G1*(2*l1k2*l2k1 - l1k2*Q2e + l2k1*Q2e - 2*l1k2*Q2h +
		            2*l2k1*Q2h + l1k2*Q2k - l2k1*Q2k -
		            k1k2*(2*l1k1 + 2*l1k2 + 2*l2k1 - Q2e - 2*Q2h + Q2k) +
		            l1k1*(4*l2k1 - Q2e - 2*Q2h + Q2k) + 2*pow(k1k2,2) + pow(Q2h,2)) \
		+ (G2 + G3)*(4*l1k2*l2k1*M2 - 2*l1k2*M2*Q2e + 2*l2k1*M2*Q2e -
		            4*l1k2*M2*Q2h + 4*l2k1*M2*Q2h + 2*l1k2*M2*Q2k -
		            2*l2k1*M2*Q2k - 2*l2k1*Q2h*Q2k + 4*l1k2*Q2h*S +
		            4*l2k1*Q2h*S + 4*Q2h*Q2k*S -
		            2*k1k2*(2*l1k1*M2 + 2*l1k2*M2 + 2*l2k1*M2 - M2*Q2e -
		               2*M2*Q2h + M2*Q2k + 2*Q2h*S) +
		            2*l1k1*(4*l2k1*M2 + M2*(-Q2e - 2*Q2h + Q2k) + 2*Q2h*S) -
		            4*l2k1*Q2h*Sk - 4*Q2h*Q2k*Sk + 8*Q2h*S*Sk + 4*l2k1*Q2h*Sq2 +
		            4*Q2h*Q2k*Sq2 - 8*Q2h*S*Sq2 + 8*Q2h*Sk*Sq2 +
		            4*M2*pow(k1k2,2) + 2*l2k1*pow(Q2h,2) + 2*M2*pow(Q2h,2) +
		            Q2k*pow(Q2h,2) - 4*S*pow(Q2h,2) + 2*Sk*pow(Q2h,2) -
		            2*Sq2*pow(Q2h,2) - Q2h*pow(Q2k,2) - 8*Q2h*pow(S,2) -
		            4*Q2h*pow(Sk,2) - 4*Q2h*pow(Sq2,2)))) +
		    2*l1k2*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (G1*(8*l1k1*l1k2 + 4*l1k2*(Q2e - Q2h) - Q2e*Q2h - Q2e*Q2k + Q2h*Q2k -
		          2*l1k1*(Q2h + Q2k) + 8*pow(l1k2,2) + pow(Q2h,2) + 2*pow(Q2k,2)) +
		       (G2 + G3)*(4*l1k2*M2*(Q2e - Q2h) - M2*Q2e*Q2h - M2*Q2e*Q2k +
		          M2*Q2h*Q2k - Q2e*Q2h*Q2k + 2*Q2e*Q2h*S - 4*Q2e*Q2k*S +
		          4*Q2h*Q2k*S - 2*Q2e*Q2h*Sk + 2*Q2e*Q2k*Sk - 2*Q2h*Q2k*Sk -
		          8*Q2e*S*Sk + 8*Q2h*S*Sk + 8*Q2k*S*Sk -
		          2*l1k2*Q2h*(Q2k - 4*S + 2*Sk - 4*Sq2) + 2*Q2e*Q2h*Sq2 -
		          4*Q2e*Q2k*Sq2 + 4*Q2h*Q2k*Sq2 + 16*Q2e*S*Sq2 - 16*Q2h*S*Sq2 -
		          16*Q2k*S*Sq2 - 8*Q2e*Sk*Sq2 + 8*Q2h*Sk*Sq2 + 8*Q2k*Sk*Sq2 +
		          8*M2*pow(l1k2,2) + M2*pow(Q2h,2) + Q2k*pow(Q2h,2) -
		          2*S*pow(Q2h,2) + 2*Sk*pow(Q2h,2) - 2*Sq2*pow(Q2h,2) +
		          2*M2*pow(Q2k,2) + 4*S*pow(Q2k,2) - 2*Sk*pow(Q2k,2) +
		          4*Sq2*pow(Q2k,2) + 8*Q2e*pow(S,2) - 8*Q2h*pow(S,2) -
		          8*Q2k*pow(S,2) + 4*Q2e*pow(Sk,2) - 4*Q2h*pow(Sk,2) -
		          4*Q2k*pow(Sk,2) + 8*Q2e*pow(Sq2,2) - 8*Q2h*pow(Sq2,2) -
		          8*Q2k*pow(Sq2,2) + 2*l1k1*
		           (4*l1k2*M2 - Q2h*Q2k - M2*(Q2h + Q2k) + 2*Q2h*S - 4*Q2k*S -
		             2*Q2h*Sk + 2*Q2k*Sk - 8*S*Sk + 2*Q2h*Sq2 - 4*Q2k*Sq2 +
		             16*S*Sq2 - 8*Sk*Sq2 + 8*pow(S,2) + 4*pow(Sk,2) + 8*pow(Sq2,2))))\
		)*pow(pow(l1k2,2) + pow(l2k1,2),-1);
}

long double Melem::melem2_l2k2_2(const long double l1k1,const long double l1k2,const long double k1k2,const long double Q2e,
		const long double Q2h,const long double Q2k,const long double S,const long double Sk,const long double Sq2,
		const long double f1, const long double f2)const {

	long double p1p2 = M2 + Q2h/2.;
	long double G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
	long double G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
	long double G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;
	long double l2k1 = l1k1 - Q2k/2. + Q2e/2.;
	long double l2k2 = l1k2 - k1k2 + Q2k/2 - Q2h/2;

	return m2*pow(l2k2,-2)*(l2k2*pow(k1k2 - l1k1 - l1k2,-1)*
		     (-4*G1*(4*l1k2*l2k2 + l1k2*Q2h - l2k2*Q2h + l1k2*Q2k - l2k2*Q2k -
		          k1k2*(4*l2k2 + Q2h + Q2k) + l1k1*(4*l2k2 + Q2h + Q2k) +
		          pow(Q2k,2)) + 4*(G2 + G3)*
		        (-4*l1k2*l2k2*M2 - l1k2*M2*Q2h + l2k2*M2*Q2h -
		          l1k2*M2*Q2k + l2k2*M2*Q2k - l1k2*Q2h*Q2k - 2*l1k2*Q2h*S -
		          2*l2k2*Q2h*S + 4*l1k2*Q2k*S - 4*l2k2*Q2k*S - 2*l1k2*Q2h*Sk +
		          2*l1k2*Q2k*Sk - 2*l2k2*Q2k*Sk + 8*l1k2*S*Sk - 8*l2k2*S*Sk +
		          4*Q2k*S*Sk - M2*pow(Q2k,2) + 2*S*pow(Q2k,2) + Sk*pow(Q2k,2) +
		          8*l1k2*pow(S,2) - 8*l2k2*pow(S,2) + 4*Q2k*pow(S,2) +
		          k1k2*(4*l2k2*M2 + Q2h*Q2k + M2*(Q2h + Q2k) + 2*Q2h*S -
		             4*Q2k*S + 2*Q2h*Sk - 2*Q2k*Sk - 8*S*Sk - 8*pow(S,2) -
		             4*pow(Sk,2)) - l1k1*
		           (4*l2k2*M2 + Q2h*Q2k + M2*(Q2h + Q2k) + 2*Q2h*S - 4*Q2k*S +
		             2*Q2h*Sk - 2*Q2k*Sk - 8*S*Sk - 8*pow(S,2) - 4*pow(Sk,2)) +
		          4*l1k2*pow(Sk,2) - 4*l2k2*pow(Sk,2) + 2*Q2k*pow(Sk,2))) +
		    4*l1k1*pow(2*l1k1 + 2*l1k2 + Q2e - Q2h,-1)*
		     (G1*(4*l1k2*Q2e + 4*l1k1*(2*l1k2 + Q2e - Q2h) + 8*pow(l1k1,2) +
		          4*pow(l1k2,2) + pow(Q2e,2) + pow(Q2h,2)) +
		       (G2 + G3)*(4*l1k2*M2*Q2e + 4*l1k2*Q2h*S + 2*Q2e*Q2h*S -
		          4*Q2h*Q2k*S - 4*Q2h*Q2k*Sk - 8*Q2h*S*Sk + 4*l1k2*Q2h*Sq2 +
		          2*Q2e*Q2h*Sq2 - 8*Q2h*S*Sq2 +
		          2*l1k1*(4*l1k2*M2 + 2*M2*(Q2e - Q2h) +
		             Q2h*(-Q2h + Q2k + 2*(2*S + Sk + Sq2))) + 8*M2*pow(l1k1,2) +
		          4*M2*pow(l1k2,2) + M2*pow(Q2e,2) + M2*pow(Q2h,2) +
		          Q2k*pow(Q2h,2) + 2*S*pow(Q2h,2) + 2*Sk*pow(Q2h,2) -
		          Q2h*pow(Q2k,2) - 8*Q2h*pow(S,2) - 4*Q2h*pow(Sk,2) -
		          4*Q2h*pow(Sq2,2))))*pow(pow(l1k1,2) + pow(l2k2,2),-1);
}
