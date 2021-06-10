/*
 * gamma_loop.cpp
 *
 *  Created on: Mar 28, 2017
 *      Author: razvan
 */

#include "config.h"

#include <cmath>
#include "const.h"
#include "gamma_loop.h"
#include "parameters.h"

using namespace constants;
using namespace POLARES;

Gamma_Loop::Gamma_Loop()
:melem_interf(0) {

}

int Gamma_Loop::set_param(const Parameters* param) {
	this->param = param;
	SI.set_param(param);
	m = param->m;
	m2 = pow(m,2.);
	m4 = pow(m,4.);
	m6 = pow(m,6.);
	m8 = pow(m,8.);
	m10 = pow(m,10.);
	m12 = pow(m,12.);
	m14 = pow(m,14.);
	M = param->M;
	M2 = pow(M,2.);
	c = alpha/(4.*pi);
	deltam = - (3.*SI.B0_0mm(m2) + 4.) * m;
	deltaZ1 = (SI.B0_0mm(m2) + 4. + 4.*m2*SI.C0_mm0m0m(m2));

	return 0;
}

long double Gamma_Loop::se_off_shell(const long double Q2h, const long double Q2e, const long double l1k,
		const long double S, const long double U, const long double G1, const long double G23)const {

	long double diag_se_off_shell = G1*(SI.C0_m0M0mm(-2*l1k + m2,m2)*
		      (16*m2 + (4*m4*Q2e - 4*m4*Q2h)*pow(l1k,-2) +
		        (-8*m4 - 8*m2*Q2h)*pow(l1k,-1) + 16*m4*pow(2*l1k + Q2e - Q2h,-1)) \
		+ SI.C0_m0M0mm(2*l1k + m2 + Q2e - Q2h,m2)*
		      (16*m2 - 8*m4*pow(l1k,-1) +
		        (16*m4*Q2e - 16*m4*Q2h)*pow(2*l1k + Q2e - Q2h,-2) +
		        (16*m4 + 16*m2*Q2h)*pow(2*l1k + Q2e - Q2h,-1)) +
		     SI.B0_00m(m2)*(32 + (-8*m4 - 2*m2*Q2e + 6*m2*Q2h)*pow(l1k,-2) +
		        (8*m2 + 4*Q2e + 4*Q2h)*pow(l1k,-1) +
		        (-32*l1k + 8*Q2e - 8*Q2h)*pow(2*l1k - m2,-1) +
		        (-32*m4 - 8*m2*Q2e + 24*m2*Q2h)*pow(2*l1k + Q2e - Q2h,-2) +
		        (-16*m2 - 8*Q2e - 8*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        (-32*l1k - 24*Q2e + 24*Q2h)*pow(2*l1k + m2 + Q2e - Q2h,-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*(32*m4 - 8*pow(Q2h,2)) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (4*Q2e*Q2h - 4*pow(Q2e,2) + 4*pow(Q2h,2)) +
		        pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (8*Q2e*Q2h - 8*pow(Q2e,2) + 8*pow(Q2h,2))) +
		     2*pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-2)*
		      (32*(2 + deltaZ1)*(Q2e - Q2h)*pow(l1k,3) +
		        32*(2 + deltaZ1)*pow(l1k,4) +
		        m2*(4*(1 + deltaZ1)*m2 + Q2e - (3 + 2*deltaZ1)*Q2h)*
		         pow(Q2e - Q2h,2) - 4*pow(l1k,2)*
		         (4*(1 + deltaZ1)*m2*Q2e + (9 + 4*deltaZ1)*Q2e*Q2h -
		           (5 + 3*deltaZ1)*pow(Q2e,2) - 3*(2 + deltaZ1)*pow(Q2h,2)) +
		        2*l1k*(Q2e - Q2h)*(-4*(1 + deltaZ1)*m2*Q2e - Q2e*Q2h +
		           (1 + deltaZ1)*pow(Q2e,2) + (2 + deltaZ1)*pow(Q2h,2))) +
		     SI.B0_M0m(-2*l1k + m2,m2)*(-8 +
		        (36*m4 + 2*m2*Q2e - 10*m2*Q2h)*pow(l1k,-2) +
		        (-20*m2 - 4*Q2e - 12*Q2h)*pow(l1k,-1) +
		        pow(l1k,-3)*(-2*m4*Q2e + 10*m4*Q2h - 16*pow(m,6)) +
		        (32*l1k - 8*Q2e + 8*Q2h)*pow(2*l1k - m2,-1) +
		        (24*m2 + 8*Q2e + 8*Q2h)*pow(2*l1k + Q2e - Q2h,-1) -
		        40*m4*pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1) +
		        pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*Q2e*Q2h + 8*pow(Q2e,2) - 8*pow(Q2h,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (32*pow(m,6) - 8*m2*pow(Q2h,2))) +
		     SI.B0_m0m(m2)*(-32 + (-36*m4 + 8*m2*Q2h)*pow(l1k,-2) +
		        (40*m2 + 8*Q2h)*pow(l1k,-1) +
		        pow(l1k,-3)*(2*m4*Q2e - 10*m4*Q2h + 16*pow(m,6)) +
		        (-16*m4*Q2e + 80*m4*Q2h - 128*pow(m,6))*
		         pow(2*l1k + Q2e - Q2h,-3) +
		        (-144*m4 + 32*m2*Q2h)*pow(2*l1k + Q2e - Q2h,-2) +
		        (-80*m2 - 16*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        80*m4*pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-2)*
		         (64*pow(m,6) - 16*m2*pow(Q2h,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*pow(m,6) + 8*m2*pow(Q2h,2))) +
		     SI.B0_M0m(2*l1k + m2 + Q2e - Q2h,m2)*
		      (-8 + (-12*m2 - 4*Q2e - 4*Q2h)*pow(l1k,-1) +
		        (16*m4*Q2e - 80*m4*Q2h + 128*pow(m,6))*
		         pow(2*l1k + Q2e - Q2h,-3) +
		        (144*m4 + 8*m2*Q2e - 40*m2*Q2h)*pow(2*l1k + Q2e - Q2h,-2) +
		        (40*m2 + 8*Q2e + 24*Q2h)*pow(2*l1k + Q2e - Q2h,-1) -
		        40*m4*pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1) +
		        (32*l1k + 24*Q2e - 24*Q2h)*pow(2*l1k + m2 + Q2e - Q2h,-1) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-4*Q2e*Q2h + 4*pow(Q2e,2) - 4*pow(Q2h,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-2)*
		         (-64*pow(m,6) + 16*m2*pow(Q2h,2)))) +
		  G23*(SI.C0_m0M0mm(2*l1k + m2 + Q2e - Q2h,m2)*
		      (16*m2*M2 + 4*m2*S + 4*m2*U +
		        pow(l1k,-1)*(-8*m4*M2 - 2*m2*Q2e*Q2h + 6*m2*Q2h*S - 4*m4*U +
		           2*m2*Q2e*U + 2*m2*Q2h*U - 4*m2*S*U - 2*m2*pow(Q2h,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (16*m4*M2 + 16*m2*M2*Q2h + 4*m2*Q2h*S + 8*m4*U +
		           4*m2*pow(Q2h,2) - 8*m2*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (16*m4*M2*Q2e - 16*m4*M2*Q2h + 8*m4*Q2e*Q2h - 8*m4*Q2e*S -
		           16*m4*Q2h*S - 8*m4*Q2h*U + 16*m4*S*U + 16*m4*pow(S,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-4*m4*Q2h*S - 12*m4*Q2h*U + 8*m4*S*U - 4*m2*Q2h*S*U +
		           4*m4*pow(Q2h,2) + 6*m2*S*pow(Q2h,2) + 2*m2*U*pow(Q2h,2) -
		           2*m2*pow(Q2h,3) - 4*m2*Q2h*pow(S,2) + 8*m4*pow(U,2))) +
		     SI.B0_00m(m2)*(16*Q2h - 12*S - 12*U +
		        pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (8*M2*Q2e - 16*l1k*Q2h - 16*M2*Q2h - 10*Q2e*Q2h + 12*l1k*S +
		           6*Q2e*S + 18*Q2h*S + 12*l1k*U + 12*Q2e*U - 8*Q2h*U - 12*S*U -
		           12*pow(S,2)) + pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*m2*M2 + 8*M2*Q2e - 32*M2*Q2h - 2*Q2e*Q2h - 4*m2*S +
		           6*Q2e*S - 16*Q2h*S - 4*m2*U + 10*Q2h*U - 12*S*U -
		           4*pow(Q2h,2) + 12*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (-8*m2*M2*Q2e + 24*m2*M2*Q2h - 4*m2*Q2e*Q2h + 4*m2*Q2e*S +
		           24*m2*Q2h*S + 4*m2*Q2h*U - 8*m2*S*U - 24*m2*pow(S,2)) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-16*M2*Q2e*Q2h + 5*Q2e*Q2h*S - 4*Q2e*Q2h*U - 6*Q2e*S*U +
		           8*Q2h*S*U + 4*M2*pow(Q2e,2) - Q2h*pow(Q2e,2) +
		           3*U*pow(Q2e,2) + 16*M2*pow(Q2h,2) - 2*Q2e*pow(Q2h,2) -
		           4*S*pow(Q2h,2) + 4*pow(Q2h,3) - 2*Q2h*pow(S,2)) +
		        pow(l1k,-1)*(8*m2*M2 - 4*M2*Q2e + 16*M2*Q2h + Q2e*Q2h +
		           2*m2*S - 5*Q2h*S + 2*m2*U - 3*Q2e*U + 8*Q2h*U + 6*S*U +
		           2*pow(Q2h,2) - 6*pow(U,2)) +
		        pow(2*l1k - m2,-1)*(-8*M2*Q2e - 16*l1k*Q2h + 16*M2*Q2h +
		           2*Q2e*Q2h + 12*l1k*S - 6*Q2e*S + 2*Q2h*S + 12*l1k*U - 24*Q2h*U +
		           12*S*U + 8*pow(Q2h,2) + 12*pow(U,2)) +
		        pow(l1k,-2)*(-2*m2*M2*Q2e + 6*m2*M2*Q2h - m2*Q2e*Q2h +
		           m2*Q2h*S + m2*Q2e*U + 6*m2*Q2h*U - 2*m2*S*U - 6*m2*pow(U,2)) \
		+ pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*M2*Q2e*Q2h - 8*Q2e*Q2h*S + 10*Q2e*Q2h*U - 12*Q2e*S*U +
		           16*Q2h*S*U + 8*M2*pow(Q2e,2) - 2*Q2h*pow(Q2e,2) +
		           6*S*pow(Q2e,2) + 32*M2*pow(Q2h,2) - 4*Q2e*pow(Q2h,2) -
		           8*U*pow(Q2h,2) + 8*pow(Q2h,3) - 4*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*m2*M2*Q2h + 16*m2*Q2h*S + 16*m2*Q2h*U - 24*m2*S*U -
		           4*Q2h*S*U - 12*m2*pow(Q2h,2) - 8*M2*pow(Q2h,2) -
		           2*pow(Q2h,3) - 4*m2*pow(S,2) + 2*Q2h*pow(S,2) - 4*m2*pow(U,2) +
		           2*Q2h*pow(U,2))) + pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-2)*
		      (16*(4*(2 + deltaZ1)*M2 + S + U)*pow(l1k,4) +
		        4*pow(l1k,3)*(16*(2 + deltaZ1)*M2*(Q2e - Q2h) + 3*Q2e*S +
		           3*Q2h*S + 4*deltaZ1*Q2h*S + 5*Q2e*U - 11*Q2h*U -
		           4*deltaZ1*Q2h*U - 2*pow(S,2) + 2*pow(U,2)) +
		        m2*pow(Q2e - Q2h,2)*(2*M2*(Q2e - (3 + 2*deltaZ1)*Q2h) - Q2h*S +
		           Q2e*(Q2h - U) - 6*Q2h*U - 4*deltaZ1*Q2h*U + 2*S*U + 6*pow(U,2) +
		           4*deltaZ1*pow(U,2)) +
		        2*pow(l1k,2)*(11*Q2e*Q2h*S + 8*deltaZ1*Q2e*Q2h*S - 13*Q2e*Q2h*U -
		           4*deltaZ1*Q2e*Q2h*U - 2*Q2e*S*U - 2*Q2h*S*U +
		           4*m2*(Q2h - S - U)*
		            (Q2e + (3 + 2*deltaZ1)*Q2h - 2*(2 + deltaZ1)*(S + U)) -
		           Q2h*pow(Q2e,2) + S*pow(Q2e,2) + 4*U*pow(Q2e,2) -
		           Q2e*pow(Q2h,2) - 4*S*pow(Q2h,2) - 4*deltaZ1*S*pow(Q2h,2) +
		           17*U*pow(Q2h,2) + 8*deltaZ1*U*pow(Q2h,2) +
		           4*M2*(-((9 + 4*deltaZ1)*Q2e*Q2h) +
		              (5 + 3*deltaZ1)*pow(Q2e,2) + 3*(2 + deltaZ1)*pow(Q2h,2)) -
		           2*Q2e*pow(S,2) - 4*Q2h*pow(S,2) - 4*deltaZ1*Q2h*pow(S,2) +
		           4*Q2e*pow(U,2) - 10*Q2h*pow(U,2) - 4*deltaZ1*Q2h*pow(U,2)) +
		        l1k*(Q2e - Q2h)*(7*Q2e*Q2h*S + 4*deltaZ1*Q2e*Q2h*S - Q2e*Q2h*U -
		           2*Q2e*S*U - 2*Q2h*S*U - Q2h*pow(Q2e,2) + U*pow(Q2e,2) -
		           Q2e*pow(Q2h,2) + S*pow(Q2h,2) + 8*U*pow(Q2h,2) +
		           4*deltaZ1*U*pow(Q2h,2) +
		           4*M2*(-(Q2e*Q2h) + (1 + deltaZ1)*pow(Q2e,2) +
		              (2 + deltaZ1)*pow(Q2h,2)) +
		           2*m2*(Q2e*(2*Q2h - S - 3*U) +
		              2*(S + U)*(S + 7*U + 4*deltaZ1*U) -
		              Q2h*(9*S + 4*deltaZ1*S + 19*U + 12*deltaZ1*U) +
		              (6 + 4*deltaZ1)*pow(Q2h,2)) - 6*Q2h*pow(S,2) -
		           4*deltaZ1*Q2h*pow(S,2) + 2*Q2e*pow(U,2) - 8*Q2h*pow(U,2) -
		           4*deltaZ1*Q2h*pow(U,2))) +
		     SI.B0_M0m(2*l1k + m2 + Q2e - Q2h,m2)*
		      (8*M2 - 8*Q2h + 10*S + 10*U +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (40*m2*M2 + 8*M2*Q2e + 4*m2*Q2h + 24*M2*Q2h + 2*Q2e*Q2h -
		           4*m2*S - 2*Q2e*S + 14*Q2h*S + 12*m2*U + 2*Q2h*U + 4*S*U +
		           2*pow(Q2h,2) - 16*pow(S,2)) +
		        pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-8*M2*Q2e + 16*l1k*Q2h + 16*M2*Q2h + 10*Q2e*Q2h - 12*l1k*S -
		           6*Q2e*S - 18*Q2h*S - 12*l1k*U - 12*Q2e*U + 8*Q2h*U + 12*S*U +
		           12*pow(S,2)) + pow(2*l1k + Q2e - Q2h,-2)*
		         (16*m4*M2 + 8*m2*M2*Q2e + 16*m4*Q2h - 40*m2*M2*Q2h +
		           12*m2*Q2e*Q2h - 32*m4*S - 12*m2*Q2e*S - 36*m2*Q2h*S +
		           8*m4*U - 28*m2*Q2h*U + 24*m2*S*U + 4*m2*pow(Q2h,2) +
		           32*m2*pow(S,2)) + pow(2*l1k + Q2e - Q2h,-3)*
		         (16*m4*M2*Q2e - 80*m4*M2*Q2h + 8*m4*Q2e*Q2h - 8*m4*Q2e*S -
		           80*m4*Q2h*S - 8*m4*Q2h*U + 16*m4*S*U + 80*m4*pow(S,2)) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (16*M2*Q2e*Q2h - 5*Q2e*Q2h*S + 4*Q2e*Q2h*U + 6*Q2e*S*U -
		           8*Q2h*S*U - 4*M2*pow(Q2e,2) + Q2h*pow(Q2e,2) -
		           3*U*pow(Q2e,2) - 16*M2*pow(Q2h,2) + 2*Q2e*pow(Q2h,2) +
		           4*S*pow(Q2h,2) - 4*pow(Q2h,3) + 2*Q2h*pow(S,2)) +
		        pow(l1k,-1)*(-12*m2*M2 + 4*M2*Q2e - 2*m2*Q2h - 16*M2*Q2h -
		           2*Q2e*Q2h + 8*Q2h*S - 8*m2*U + 4*Q2e*U - Q2h*U - 8*S*U -
		           5*pow(Q2h,2) + 2*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*m4*M2 + 8*m4*Q2h + 8*m2*M2*Q2h + 4*m2*Q2h*S -
		           20*m4*U - 18*m2*Q2h*U + 16*m2*S*U - 2*Q2h*S*U +
		           2*m2*pow(Q2h,2) + 5*S*pow(Q2h,2) - U*pow(Q2h,2) - pow(Q2h,3) -
		           4*Q2h*pow(S,2) + 12*m2*pow(U,2) + 2*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-2)*
		         (32*m4*M2*Q2h - 20*m4*Q2h*S - 28*m4*Q2h*U + 40*m4*S*U -
		           4*m2*Q2h*S*U + 20*m4*pow(Q2h,2) + 16*m2*M2*pow(Q2h,2) +
		           14*m2*S*pow(Q2h,2) + 10*m2*U*pow(Q2h,2) - 2*m2*pow(Q2h,3) -
		           12*m2*Q2h*pow(S,2) + 8*m4*pow(U,2) - 8*m2*Q2h*pow(U,2))) +
		     SI.C0_m0M0mm(-2*l1k + m2,m2)*
		      (16*m2*M2 + 4*m2*S + 4*m2*U +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (16*m4*M2 + 4*m2*Q2e*Q2h + 8*m4*S - 4*m2*Q2e*S -
		           4*m2*Q2h*S - 12*m2*Q2h*U + 8*m2*S*U + 4*m2*pow(Q2h,2)) +
		        pow(l1k,-1)*(-8*m4*M2 - 8*m2*M2*Q2h - 4*m4*S - 2*m2*Q2h*U -
		           2*m2*pow(Q2h,2) + 4*m2*pow(U,2)) +
		        pow(l1k,-2)*(4*m4*M2*Q2e - 4*m4*M2*Q2h + 2*m4*Q2e*Q2h -
		           2*m4*Q2h*S - 2*m4*Q2e*U - 4*m4*Q2h*U + 4*m4*S*U +
		           4*m4*pow(U,2)) + pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-12*m4*Q2h*S - 4*m4*Q2h*U + 8*m4*S*U - 4*m2*Q2h*S*U +
		           4*m4*pow(Q2h,2) + 2*m2*S*pow(Q2h,2) + 6*m2*U*pow(Q2h,2) -
		           2*m2*pow(Q2h,3) + 8*m4*pow(S,2) - 4*m2*Q2h*pow(U,2))) +
		     SI.B0_M0m(-2*l1k + m2,m2)*(8*M2 - 8*Q2h + 10*S + 10*U +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (24*m2*M2 - 8*M2*Q2e + 4*m2*Q2h + 32*M2*Q2h + 4*Q2e*Q2h +
		           16*m2*S - 8*Q2e*S + 2*Q2h*S - 16*Q2h*U + 16*S*U +
		           10*pow(Q2h,2) - 4*pow(S,2)) +
		        pow(2*l1k - m2,-1)*(8*M2*Q2e + 16*l1k*Q2h - 16*M2*Q2h -
		           2*Q2e*Q2h - 12*l1k*S + 6*Q2e*S - 2*Q2h*S - 12*l1k*U + 24*Q2h*U -
		           12*S*U - 8*pow(Q2h,2) - 12*pow(U,2)) +
		        pow(l1k,-1)*(-20*m2*M2 - 4*M2*Q2e - 2*m2*Q2h - 12*M2*Q2h -
		           Q2e*Q2h - 6*m2*S - Q2h*S + 2*m2*U + Q2e*U - 7*Q2h*U - 2*S*U -
		           pow(Q2h,2) + 8*pow(U,2)) +
		        pow(l1k,-2)*(4*m4*M2 + 2*m2*M2*Q2e + 4*m4*Q2h -
		           10*m2*M2*Q2h + 3*m2*Q2e*Q2h + 2*m4*S - 7*m2*Q2h*S -
		           8*m4*U - 3*m2*Q2e*U - 9*m2*Q2h*U + 6*m2*S*U +
		           m2*pow(Q2h,2) + 8*m2*pow(U,2)) +
		        pow(l1k,-3)*(-2*m4*M2*Q2e + 10*m4*M2*Q2h - m4*Q2e*Q2h +
		           m4*Q2h*S + m4*Q2e*U + 10*m4*Q2h*U - 2*m4*S*U -
		           10*m4*pow(U,2)) + pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*m4*M2 + 8*m4*Q2h + 8*m2*M2*Q2h - 20*m4*S -
		           18*m2*Q2h*S + 4*m2*Q2h*U + 16*m2*S*U - 2*Q2h*S*U +
		           2*m2*pow(Q2h,2) - S*pow(Q2h,2) + 5*U*pow(Q2h,2) - pow(Q2h,3) +
		           12*m2*pow(S,2) + 2*Q2h*pow(S,2) - 4*Q2h*pow(U,2)) +
		        pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (32*M2*Q2e*Q2h + 8*Q2e*Q2h*S - 10*Q2e*Q2h*U + 12*Q2e*S*U -
		           16*Q2h*S*U - 8*M2*pow(Q2e,2) + 2*Q2h*pow(Q2e,2) -
		           6*S*pow(Q2e,2) - 32*M2*pow(Q2h,2) + 4*Q2e*pow(Q2h,2) +
		           8*U*pow(Q2h,2) - 8*pow(Q2h,3) + 4*Q2h*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*m4*M2*Q2h + 14*m4*Q2h*S + 10*m4*Q2h*U - 20*m4*S*U +
		           2*m2*Q2h*S*U - 10*m4*pow(Q2h,2) - 8*m2*M2*pow(Q2h,2) -
		           5*m2*S*pow(Q2h,2) - 7*m2*U*pow(Q2h,2) + m2*pow(Q2h,3) -
		           4*m4*pow(S,2) + 4*m2*Q2h*pow(S,2) + 6*m2*Q2h*pow(U,2))) +
		     SI.B0_m0m(m2)*(-32*M2 - 8*S - 8*U +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-48*m2*M2 - 8*m2*Q2h - 16*M2*Q2h - 4*Q2e*Q2h - 8*m2*S +
		           4*Q2e*S - 8*m2*U + 12*Q2h*U - 8*S*U - 8*pow(Q2h,2) + 8*pow(S,2)) \
		+ pow(2*l1k + Q2e - Q2h,-2)*(-16*m4*M2 - 16*m4*Q2h + 32*m2*M2*Q2h -
		           8*m2*Q2e*Q2h + 32*m4*S + 8*m2*Q2e*S + 28*m2*Q2h*S - 8*m4*U +
		           24*m2*Q2h*U - 16*m2*S*U - 4*m2*pow(Q2h,2) - 24*m2*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-3)*
		         (-16*m4*M2*Q2e + 80*m4*M2*Q2h - 8*m4*Q2e*Q2h + 8*m4*Q2e*S +
		           80*m4*Q2h*S + 8*m4*Q2h*U - 16*m4*S*U - 80*m4*pow(S,2)) +
		        pow(l1k,-1)*(24*m2*M2 + 4*m2*Q2h + 8*M2*Q2h + 2*Q2e*Q2h +
		           4*m2*S - 6*Q2h*S + 4*m2*U - 2*Q2e*U + 4*S*U + 4*pow(Q2h,2) -
		           4*pow(U,2)) + pow(l1k,-2)*
		         (-4*m4*M2 - 4*m4*Q2h + 8*m2*M2*Q2h - 2*m2*Q2e*Q2h -
		           2*m4*S + 6*m2*Q2h*S + 8*m4*U + 2*m2*Q2e*U + 7*m2*Q2h*U -
		           4*m2*S*U - m2*pow(Q2h,2) - 6*m2*pow(U,2)) +
		        pow(l1k,-3)*(2*m4*M2*Q2e - 10*m4*M2*Q2h + m4*Q2e*Q2h -
		           m4*Q2h*S - m4*Q2e*U - 10*m4*Q2h*U + 2*m4*S*U + 10*m4*pow(U,2)\
		) + pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (16*m4*M2 - 16*m4*Q2h - 16*m2*M2*Q2h + 20*m4*S +
		           6*m2*Q2h*S + 20*m4*U + 6*m2*Q2h*U - 24*m2*S*U + 8*Q2h*S*U -
		           8*S*pow(Q2h,2) - 8*U*pow(Q2h,2) + 4*pow(Q2h,3) - 8*m2*pow(S,2) +
		           4*Q2h*pow(S,2) - 8*m2*pow(U,2) + 4*Q2h*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (16*m4*M2*Q2h - 14*m4*Q2h*S - 10*m4*Q2h*U + 20*m4*S*U -
		           2*m2*Q2h*S*U + 10*m4*pow(Q2h,2) + 8*m2*M2*pow(Q2h,2) +
		           5*m2*S*pow(Q2h,2) + 7*m2*U*pow(Q2h,2) - m2*pow(Q2h,3) +
		           4*m4*pow(S,2) - 4*m2*Q2h*pow(S,2) - 6*m2*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-2)*
		         (-32*m4*M2*Q2h + 20*m4*Q2h*S + 28*m4*Q2h*U - 40*m4*S*U +
		           4*m2*Q2h*S*U - 20*m4*pow(Q2h,2) - 16*m2*M2*pow(Q2h,2) -
		           14*m2*S*pow(Q2h,2) - 10*m2*U*pow(Q2h,2) + 2*m2*pow(Q2h,3) +
		           12*m2*Q2h*pow(S,2) - 8*m4*pow(U,2) + 8*m2*Q2h*pow(U,2))));

	melem_interf = c*diag_se_off_shell/M1gamma(Q2h, Q2e, l1k, S, U, G1, G23);

	return 2.*melem_interf;

}

long double Gamma_Loop::se_on_shell(const long double Q2h, const long double Q2e, const long double l1k,
		const long double S, const long double U, const long double G1, const long double G23)const {

	long double diag_se_on_shell = G1*(SI.B0_M0m(2*l1k + m2 + Q2e - Q2h,m2)*
		      (32 + 8*Q2e*pow(l1k,-1) +
		        (64*m4*Q2h - 128*pow(m,6))*pow(2*l1k + Q2e - Q2h,-3) +
		        (-128*m4 + 16*m2*Q2e + 32*m2*Q2h)*pow(2*l1k + Q2e - Q2h,-2) -
		        8*Q2e*pow(2*l1k + Q2e - Q2h,-1) +
		        32*m4*pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1) +
		        (-48*l1k - 40*Q2e + 32*Q2h)*pow(2*l1k + m2 + Q2e - Q2h,-1) +
		        pow(l1k,-1)*(12*Q2e*Q2h - 8*pow(Q2e,2))*
		         pow(2*l1k + m2 + Q2e - Q2h,-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-2)*
		         (64*pow(m,6) - 16*m2*pow(Q2h,2))) +
		     SI.B0_00m(m2)*(-48 + (-8*m4 + 2*m2*Q2e)*pow(l1k,-2) +
		        (8*m2 - 8*Q2e + 4*Q2h)*pow(l1k,-1) +
		        pow(l1k,-3)*(-4*m4*Q2h + 8*pow(m,6)) +
		        (48*l1k - 16*Q2e + 8*Q2h)*pow(2*l1k - m2,-1) +
		        (32*m4*Q2h - 64*pow(m,6))*pow(2*l1k + Q2e - Q2h,-3) +
		        (-32*m4 + 8*m2*Q2e)*pow(2*l1k + Q2e - Q2h,-2) +
		        (-16*m2 + 16*Q2e - 8*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        pow(2*l1k - m2,-1)*(-24*Q2e*Q2h + 16*pow(Q2e,2))*
		         pow(2*l1k + Q2e - Q2h,-1) +
		        (48*l1k + 40*Q2e - 32*Q2h)*pow(2*l1k + m2 + Q2e - Q2h,-1) +
		        pow(l1k,-1)*(-12*Q2e*Q2h + 8*pow(Q2e,2))*
		         pow(2*l1k + m2 + Q2e - Q2h,-1) +
		        8*pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*pow(Q2h,2) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-2)*
		         (32*pow(m,6) - 8*m2*pow(Q2h,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*pow(m,6) + 4*m2*pow(Q2h,2))) +
		     SI.B0_M0m(-2*l1k + m2,m2)*(32 +
		        (-32*m4 + 4*m2*Q2e + 8*m2*Q2h)*pow(l1k,-2) + 4*Q2e*pow(l1k,-1) +
		        pow(l1k,-3)*(-8*m4*Q2h + 16*pow(m,6)) +
		        (-48*l1k + 16*Q2e - 8*Q2h)*pow(2*l1k - m2,-1) -
		        16*Q2e*pow(2*l1k + Q2e - Q2h,-1) +
		        32*m4*pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1) +
		        pow(2*l1k - m2,-1)*(24*Q2e*Q2h - 16*pow(Q2e,2))*
		         pow(2*l1k + Q2e - Q2h,-1) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*pow(m,6) + 8*m2*pow(Q2h,2))) -
		     2*pow(l1k,-3)*pow(2*l1k + Q2e - Q2h,-3)*
		      (96*(1 + deltaZ1)*(Q2e - Q2h)*pow(l1k,5) +
		        64*(1 + deltaZ1)*pow(l1k,6) +
		        16*(Q2e - Q2h)*pow(l1k,3)*
		         (-4*m4 + deltam*m*(4*m2 - 2*Q2e - Q2h) +
		           m2*(-2*deltaZ1*Q2e + Q2h) +
		           (1 + deltaZ1)*(-(Q2e*Q2h) + pow(Q2e,2) + pow(Q2h,2))) +
		        8*pow(l1k,4)*(-8*m4 + m2*(-4*deltaZ1*Q2e + 2*Q2h) +
		           deltam*(-2*m*(2*Q2e + Q2h) + 8*pow(m,3)) +
		           (1 + deltaZ1)*(-12*Q2e*Q2h + 7*pow(Q2e,2) + 7*pow(Q2h,2))) -
		        l1k*m*pow(Q2e - Q2h,2)*
		         (deltam*(4*m2 - Q2e)*(4*m2 - Q2e - Q2h) +
		           m*(-16*m4 + Q2e*(Q2h + 2*deltaZ1*Q2h) +
		              4*m2*(Q2e - deltaZ1*Q2e + (2 + deltaZ1)*Q2h) - pow(Q2e,2) -
		              2*(1 + deltaZ1)*pow(Q2h,2))) +
		        2*(Q2e - Q2h)*pow(l1k,2)*
		         (4*m4*((-3 + deltaZ1)*Q2e - deltaZ1*Q2h) + 16*pow(m,6) +
		           deltam*m*(-16*m4 + 4*m2*(4*Q2e - Q2h) + Q2e*Q2h -
		              5*pow(Q2e,2) + 2*pow(Q2h,2)) +
		           m2*(Q2e*(Q2h + 2*deltaZ1*Q2h) + (1 - 4*deltaZ1)*pow(Q2e,2) +
		              2*deltaZ1*pow(Q2h,2)) +
		           (1 + deltaZ1)*(-(Q2h*pow(Q2e,2)) + pow(Q2e,3) + Q2e*pow(Q2h,2) -
		              pow(Q2h,3))) + 2*(deltam - m)*(2*m2 - Q2h)*pow(m,3)*
		         pow(-Q2e + Q2h,3))) + G23*
		   (SI.B0_M0m(-2*l1k + m2,m2)*(16*M2 + 8*Q2h - 4*S - 4*U +
		        (4*M2*Q2e + 2*Q2e*Q2h - 2*Q2h*S - 4*m2*U - 2*Q2e*U + 4*S*U)*
		         pow(l1k,-1) + (-8*m4*Q2h - 8*m2*M2*Q2h + 16*m4*S +
		           4*m2*Q2h*S - 4*m2*Q2h*U - 8*m2*S*U)*pow(l1k,-1)*
		         pow(2*l1k + Q2e - Q2h,-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-24*M2*Q2h - 8*m2*S - 8*Q2h*S - 4*pow(Q2h,2) + 8*pow(S,2)) +
		        pow(2*l1k - m2,-1)*(-16*l1k*M2 - 16*l1k*Q2h + 16*M2*Q2h +
		           8*l1k*S - 2*Q2e*S + 10*Q2h*S + 8*l1k*U + 2*Q2e*U - 10*Q2h*U +
		           4*pow(Q2h,2) - 4*pow(S,2) + 4*pow(U,2)) +
		        pow(l1k,-2)*(4*m2*M2*Q2e - 4*m4*Q2h + 8*m2*M2*Q2h +
		           4*m2*Q2h*S + 8*m4*U + 4*m2*Q2h*U - 4*m2*pow(U,2)) +
		        pow(l1k,-3)*(-8*m4*M2*Q2h - 8*m4*Q2h*U + 8*m4*pow(U,2)) +
		        pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*M2*Q2e*Q2h - 6*Q2e*Q2h*S + 6*Q2e*Q2h*U - 4*Q2e*S*U +
		           4*Q2h*S*U - 2*Q2h*pow(Q2e,2) + 2*S*pow(Q2e,2) +
		           24*M2*pow(Q2h,2) - 2*Q2e*pow(Q2h,2) + 8*S*pow(Q2h,2) -
		           2*U*pow(Q2h,2) + 4*pow(Q2h,3) + 4*Q2e*pow(S,2) -
		           8*Q2h*pow(S,2) - 4*Q2h*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (16*m4*M2*Q2h - 8*m4*Q2h*S - 8*m4*Q2h*U + 16*m4*S*U +
		           8*m4*pow(Q2h,2) + 8*m2*M2*pow(Q2h,2) + 4*m2*S*pow(Q2h,2) +
		           4*m2*U*pow(Q2h,2) - 4*m2*Q2h*pow(S,2) - 4*m2*Q2h*pow(U,2))) +
		     SI.B0_00m(m2)*(-16*M2 - 16*Q2h + 8*S + 8*U +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (16*M2*Q2h + 4*Q2e*Q2h - 4*Q2e*S + 8*Q2h*S - 12*Q2h*U + 8*S*U +
		           4*pow(Q2h,2) - 8*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (8*m2*M2*Q2e - 8*m4*Q2h + 16*m4*S - 8*m2*Q2h*S +
		           8*m2*Q2h*U + 8*m2*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-3)*
		         (32*m4*M2*Q2h + 32*m4*Q2h*S - 32*m4*pow(S,2)) +
		        pow(2*l1k - m2,-1)*(16*l1k*M2 + 16*l1k*Q2h - 16*M2*Q2h -
		           8*l1k*S + 2*Q2e*S - 10*Q2h*S - 8*l1k*U - 2*Q2e*U + 10*Q2h*U -
		           4*pow(Q2h,2) + 4*pow(S,2) - 4*pow(U,2)) +
		        pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (16*l1k*M2 + 8*M2*Q2e + 16*l1k*Q2h + 8*M2*Q2h + 8*Q2e*Q2h -
		           8*l1k*S - 2*Q2e*S - 6*Q2h*S - 8*l1k*U - 6*Q2e*U + 14*Q2h*U -
		           4*pow(Q2h,2) + 4*pow(S,2) - 4*pow(U,2)) +
		        pow(l1k,-1)*(-8*M2*Q2h - 2*Q2e*Q2h + 6*Q2h*S + 2*Q2e*U -
		           4*Q2h*U - 4*S*U - 2*pow(Q2h,2) + 4*pow(U,2)) +
		        pow(l1k,-2)*(2*m2*M2*Q2e - 2*m4*Q2h + 2*m2*Q2h*S + 4*m4*U -
		           2*m2*Q2h*U + 2*m2*pow(U,2)) +
		        pow(l1k,-3)*(-4*m4*M2*Q2h - 4*m4*Q2h*U + 4*m4*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*m4*Q2h + 8*m2*M2*Q2h + 8*m4*S - 8*m2*Q2h*S + 8*m4*U -
		           8*m2*Q2h*U + 8*m2*S*U + 8*m2*pow(Q2h,2) + 8*M2*pow(Q2h,2) +
		           4*S*pow(Q2h,2) + 4*U*pow(Q2h,2) - 4*Q2h*pow(S,2) - 4*Q2h*pow(U,2)\
		) + pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (16*M2*Q2e*Q2h + 6*Q2e*Q2h*S - 6*Q2e*Q2h*U + 4*Q2e*S*U -
		           4*Q2h*S*U + 2*Q2h*pow(Q2e,2) - 2*S*pow(Q2e,2) -
		           24*M2*pow(Q2h,2) + 2*Q2e*pow(Q2h,2) - 8*S*pow(Q2h,2) +
		           2*U*pow(Q2h,2) - 4*pow(Q2h,3) - 4*Q2e*pow(S,2) +
		           8*Q2h*pow(S,2) + 4*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (8*M2*Q2e*Q2h - 3*Q2e*Q2h*S + 3*Q2e*Q2h*U + 2*Q2e*S*U -
		           2*Q2h*S*U + Q2h*pow(Q2e,2) - U*pow(Q2e,2) - 12*M2*pow(Q2h,2) +
		           Q2e*pow(Q2h,2) + S*pow(Q2h,2) - 4*U*pow(Q2h,2) - 2*pow(Q2h,3) +
		           2*Q2h*pow(S,2) - 2*Q2e*pow(U,2) + 4*Q2h*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (8*m4*M2*Q2h - 4*m4*Q2h*S - 4*m4*Q2h*U + 8*m4*S*U +
		           4*m4*pow(Q2h,2) + 4*m2*M2*pow(Q2h,2) + 2*m2*S*pow(Q2h,2) +
		           2*m2*U*pow(Q2h,2) - 2*m2*Q2h*pow(S,2) - 2*m2*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-2)*
		         (-16*m4*M2*Q2h + 8*m4*Q2h*S + 8*m4*Q2h*U - 16*m4*S*U -
		           8*m4*pow(Q2h,2) - 8*m2*M2*pow(Q2h,2) - 4*m2*S*pow(Q2h,2) -
		           4*m2*U*pow(Q2h,2) + 4*m2*Q2h*pow(S,2) + 4*m2*Q2h*pow(U,2))) +
		     SI.B0_M0m(2*l1k + m2 + Q2e - Q2h,m2)*
		      (16*M2 + 8*Q2h - 4*S - 4*U +
		        (-8*M2*Q2e - 4*Q2e*Q2h + 8*m2*S + 4*Q2e*S + 4*Q2h*U - 8*S*U)*
		         pow(2*l1k + Q2e - Q2h,-1) +
		        (-8*m4*Q2h - 8*m2*M2*Q2h - 4*m2*Q2h*S + 16*m4*U +
		           4*m2*Q2h*U - 8*m2*S*U)*pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (16*m2*M2*Q2e - 16*m4*Q2h + 32*m2*M2*Q2h + 32*m4*S +
		           16*m2*Q2h*S + 16*m2*Q2h*U - 16*m2*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-3)*
		         (64*m4*M2*Q2h + 64*m4*Q2h*S - 64*m4*pow(S,2)) +
		        pow(l1k,-1)*(12*M2*Q2h + 4*m2*U + 4*Q2h*U + 2*pow(Q2h,2) -
		           4*pow(U,2)) + pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-16*l1k*M2 - 8*M2*Q2e - 16*l1k*Q2h - 8*M2*Q2h - 8*Q2e*Q2h +
		           8*l1k*S + 2*Q2e*S + 6*Q2h*S + 8*l1k*U + 6*Q2e*U - 14*Q2h*U +
		           4*pow(Q2h,2) - 4*pow(S,2) + 4*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-8*M2*Q2e*Q2h + 3*Q2e*Q2h*S - 3*Q2e*Q2h*U - 2*Q2e*S*U +
		           2*Q2h*S*U - Q2h*pow(Q2e,2) + U*pow(Q2e,2) + 12*M2*pow(Q2h,2) -
		           Q2e*pow(Q2h,2) - S*pow(Q2h,2) + 4*U*pow(Q2h,2) + 2*pow(Q2h,3) -
		           2*Q2h*pow(S,2) + 2*Q2e*pow(U,2) - 4*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-2)*
		         (-32*m4*M2*Q2h + 16*m4*Q2h*S + 16*m4*Q2h*U - 32*m4*S*U -
		           16*m4*pow(Q2h,2) - 16*m2*M2*pow(Q2h,2) - 8*m2*S*pow(Q2h,2) -
		           8*m2*U*pow(Q2h,2) + 8*m2*Q2h*pow(S,2) + 8*m2*Q2h*pow(U,2))) -
		     2*pow(l1k,-3)*pow(2*l1k + Q2e - Q2h,-3)*
		      (16*(1 + deltaZ1)*(6*M2*(Q2e - Q2h) + Q2h*(S - U))*pow(l1k,5) +
		        64*(1 + deltaZ1)*M2*pow(l1k,6) +
		        2*(deltam - m)*(M2*Q2h + (Q2h - U)*U)*pow(m,3)*pow(Q2e - Q2h,3) -
		        l1k*m*pow(Q2e - Q2h,2)*
		         (deltam*(Q2e*Q2h*S + Q2e*Q2h*U +
		              M2*(-8*m2*Q2h + Q2e*Q2h + pow(Q2e,2)) +
		              m2*(-(Q2e*(Q2h - 2*U)) + 4*U*(S + 3*U) - 2*Q2h*(S + 8*U) +
		                 3*pow(Q2h,2)) - Q2h*pow(S,2) - Q2e*pow(U,2)) +
		           m*(-(Q2e*Q2h*S) + Q2e*Q2h*U + 2*deltaZ1*Q2e*Q2h*U +
		              M2*(-2*Q2h*(-4*m2 + Q2h + deltaZ1*Q2h) +
		                 Q2e*(Q2h + 2*deltaZ1*Q2h) - pow(Q2e,2)) +
		              m2*(Q2e*(Q2h - 2*U) - 4*U*(S + 3*U) + 2*Q2h*(S + 8*U) -
		                 3*pow(Q2h,2)) - 2*U*pow(Q2h,2) - 2*deltaZ1*U*pow(Q2h,2) +
		              Q2h*pow(S,2) - Q2e*pow(U,2) - 2*deltaZ1*Q2e*pow(U,2) +
		              2*Q2h*pow(U,2) + 2*deltaZ1*Q2h*pow(U,2))) +
		        2*(Q2e - Q2h)*pow(l1k,2)*
		         ((1 + deltaZ1)*(Q2e - Q2h)*
		            (M2*(pow(Q2e,2) + pow(Q2h,2)) +
		              Q2h*(Q2e*S + (Q2h - U)*U - pow(S,2))) +
		           m4*(-8*M2*Q2h - 4*Q2h*S - 22*Q2h*U + 4*S*U +
		              Q2e*(-5*Q2h + 2*S + 8*U) + 7*pow(Q2h,2) + 12*pow(U,2)) +
		           deltam*m*(-3*Q2e*Q2h*S - 3*Q2e*Q2h*U + 2*Q2e*S*U - 2*Q2h*S*U +
		              M2*(-(Q2e*Q2h) + 2*Q2h*(4*m2 + Q2h) - 3*pow(Q2e,2)) +
		              m2*(-4*U*(S + 3*U) + Q2h*(4*S + 22*U) +
		                 Q2e*(5*Q2h - 2*(S + 4*U)) - 7*pow(Q2h,2)) +
		              2*S*pow(Q2h,2) + 2*U*pow(Q2h,2) + Q2h*pow(S,2) +
		              3*Q2e*pow(U,2) - 2*Q2h*pow(U,2)) +
		           m2*(M2*(-(Q2e*(Q2h + 2*deltaZ1*Q2h)) + 3*pow(Q2e,2) +
		                 2*deltaZ1*pow(Q2h,2)) +
		              Q2e*((1 + 2*deltaZ1)*U*(2*S + 3*U) +
		                 Q2h*(S - 2*deltaZ1*S - (5 + 8*deltaZ1)*U) +
		                 2*(1 + deltaZ1)*pow(Q2h,2)) -
		              Q2h*(2*S*(U + 2*deltaZ1*U) - 2*Q2h*(3*U + deltaZ1*(S + 4*U)) +
		                 2*(1 + deltaZ1)*pow(Q2h,2) + pow(S,2) +
		                 2*(2 + 3*deltaZ1)*pow(U,2)))) +
		        4*pow(l1k,3)*(deltam*m*
		            (-((Q2e - Q2h)*(4*Q2h - S - 3*U)*(S + U)) -
		              4*M2*(pow(Q2e,2) - pow(Q2h,2)) +
		              2*m2*(4*Q2e*Q2h - 3*Q2e*S + Q2h*S - 5*Q2e*U + 7*Q2h*U -
		                 4*pow(Q2h,2) + 2*pow(S,2) - 2*pow(U,2))) +
		           2*m4*(Q2e*(-4*Q2h + 3*S + 5*U) - Q2h*(S + 7*U) + 4*pow(Q2h,2) -
		              2*pow(S,2) + 2*pow(U,2)) +
		           m2*(Q2e - Q2h)*(4*M2*(Q2e + Q2h) -
		              2*Q2h*(S + 3*deltaZ1*S + 3*U + 5*deltaZ1*U) +
		              4*(1 + deltaZ1)*pow(Q2h,2) +
		              (1 + 2*deltaZ1)*(4*S*U + pow(S,2) + 3*pow(U,2))) +
		           (1 + deltaZ1)*(Q2e - Q2h)*
		            (4*M2*(-(Q2e*Q2h) + pow(Q2e,2) + pow(Q2h,2)) -
		              Q2h*(Q2h*(S - 3*U) + Q2e*(-3*S + U) + 2*(pow(S,2) + pow(U,2))))\
		) + 8*pow(l1k,4)*(4*m4*(-Q2h + S + U) +
		           deltam*m*(-2*M2*(Q2e + Q2h) + 4*m2*(Q2h - S - U) -
		              (2*Q2h - S - U)*(S + U)) +
		           (1 + deltaZ1)*(M2*(-12*Q2e*Q2h + 7*pow(Q2e,2) + 7*pow(Q2h,2)) -
		              Q2h*(-3*Q2e*S + 2*Q2h*S + 2*Q2e*U - 3*Q2h*U + pow(S,2) +
		                 pow(U,2))) + m2*
		            (2*M2*(Q2e + Q2h) - 2*(1 + 2*deltaZ1)*Q2h*(S + U) +
		              2*(1 + deltaZ1)*pow(Q2h,2) + (1 + 2*deltaZ1)*pow(S + U,2)))));

	melem_interf = c*diag_se_on_shell/M1gamma(Q2h, Q2e, l1k, S, U, G1, G23);

	return 2.*melem_interf;

}

long double Gamma_Loop::vb_on_shell(const long double Q2h, const long double Q2e, const long double l1k,
		const long double S, const long double U, const long double G1, const long double G23)const {

	long double diag_vb_on_shell = G1*(32 + (8*m4 - 8*m2*Q2h)*pow(l1k,-2) +
		     (-16*m2 + 8*Q2e + 9*Q2h)*pow(l1k,-1) +
		     (32*m4 - 32*m2*Q2h)*pow(2*l1k + Q2e - Q2h,-2) +
		     (32*m2 - 16*Q2e - 18*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		     pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		      (-32*m4 + 16*m2*Q2h + 16*pow(Q2h,2)) +
		     pow(l1k,-1)*(-4*m2*pow(Q2e,2) - pow(Q2e,3))*
		      pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		     (8*l1k*Q2e + 4*l1k*Q2h + 2*pow(Q2e,2) + 2*pow(Q2h,2))*
		      pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		     pow(2*l1k + Q2e - Q2h,-1)*(-2*Q2e*pow(Q2h,2) + 2*pow(Q2h,3))*
		      pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		     SI.C0_mMQm0m(2*l1k + m2 + Q2e - Q2h,Q2h,m2)*
		      (-26*l1k + 8*m2 + 4*Q2e - 13*Q2h +
		        (-96*m4 + 16*m2*Q2h + 16*Q2e*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        pow(l1k,-1)*(32*m4 + 20*m2*Q2e - 12*m2*Q2h - (9*Q2e*Q2h)/2. +
		           (13*pow(Q2e,2))/2. - 8*pow(Q2h,2)) +
		        pow(2*l1k + Q2e - Q2h,-2)*(-128*pow(m,6) + 32*m2*pow(Q2h,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (32*m4*Q2h + 64*pow(m,6) - 16*m2*pow(Q2h,2) - 8*pow(Q2h,3)) +
		        (192*m2*Q2e*pow(l1k,3) + 24*Q2e*Q2h*pow(l1k,3) +
		           48*Q2h*pow(l1k,4) - 96*pow(l1k,5) - 96*l1k*m4*pow(Q2e,2) +
		           144*m2*pow(l1k,2)*pow(Q2e,2) - 24*Q2h*pow(l1k,2)*pow(Q2e,2) +
		           168*pow(l1k,3)*pow(Q2e,2) - 72*l1k*m2*pow(Q2e,3) -
		           96*m4*pow(Q2e,3) - 12*l1k*Q2h*pow(Q2e,3) +
		           96*pow(l1k,2)*pow(Q2e,3) - 18*l1k*pow(Q2e,4) -
		           84*m2*pow(Q2e,4) + 3*Q2h*pow(Q2e,4) - 24*pow(Q2e,5))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        pow(l1k,-1)*(-24*m4*pow(Q2e,4) - 18*m2*pow(Q2e,5) +
		           (3*Q2h*pow(Q2e,5))/2. - (9*pow(Q2e,6))/2.)*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        (64*l1k*m4 - 96*l1k*m2*Q2e - 32*m4*Q2e - 20*l1k*Q2e*Q2h +
		           32*m2*pow(l1k,2) + 64*Q2e*pow(l1k,2) - 24*Q2h*pow(l1k,2) +
		           128*pow(l1k,3) - 40*l1k*pow(Q2e,2) - 60*m2*pow(Q2e,2) +
		           2*Q2h*pow(Q2e,2) - 24*pow(Q2e,3))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*(-64*Q2e*pow(m,6) - 16*m4*pow(Q2e,2) -
		           2*m2*pow(Q2e,3) + 3*Q2h*pow(Q2e,3) - 2*pow(Q2e,4))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1)) +
		     SI.B0_M0m(2*l1k + m2 + Q2e - Q2h,m2)*
		      (-23 + (20*m2 + (13*Q2e)/2. + 4*Q2h)*pow(l1k,-1) +
		        (-32*m4 + 32*m2*Q2h)*pow(2*l1k + Q2e - Q2h,-2) +
		        (-32*m2 + 4*Q2e + 8*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        (32*l1k + 24*Q2e - 40*Q2h)*pow(2*l1k + m2 + Q2e - Q2h,-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (16*m4 - 8*m2*Q2h - 8*pow(Q2h,2)) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-16*Q2e*Q2h + 4*pow(Q2e,2) + 8*pow(Q2h,2)) +
		        (96*m2*Q2e*pow(l1k,2) + 72*Q2e*pow(l1k,3) + 48*Q2h*pow(l1k,3) -
		           48*pow(l1k,4) - 48*m4*pow(Q2e,2) - 24*l1k*Q2h*pow(Q2e,2) +
		           120*pow(l1k,2)*pow(Q2e,2) + 12*l1k*pow(Q2e,3) -
		           72*m2*pow(Q2e,3) - 27*pow(Q2e,4))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        pow(l1k,-1)*(-24*m4*pow(Q2e,3) - 24*m2*pow(Q2e,4) +
		           3*Q2h*pow(Q2e,4) - (15*pow(Q2e,5))/2.)*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        (16*l1k*m2 + 64*m4 - 4*l1k*Q2e - 72*m2*Q2e - 28*l1k*Q2h +
		           72*pow(l1k,2) - 26*pow(Q2e,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*m4*Q2e + 128*pow(m,6) + 12*Q2e*pow(Q2h,2) - 8*pow(Q2h,3))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) -
		        8*pow(2*l1k + m2 + Q2e - Q2h,-1)*pow(Q2h,3)*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*(-32*m4*Q2e - 64*pow(m,6) - 8*m2*pow(Q2e,2) +
		           9*Q2h*pow(Q2e,2) - 3*pow(Q2e,3) + 4*Q2e*pow(Q2h,2) +
		           8*pow(Q2h,3))*pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2e,2),-1) + pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-4*Q2e*pow(Q2h,3) + 8*pow(Q2h,4))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1)) +
		     (8*m2*pow(Q2e,2) + 2*pow(Q2e,3))*pow(2*l1k + Q2e - Q2h,-1)*
		      pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		     (-8*l1k*Q2e - 4*l1k*Q2h + 2*Q2e*Q2h - 2*pow(Q2e,2) + 4*pow(Q2h,2))*
		      pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		     pow(l1k,-1)*(Q2e*pow(Q2h,2) - pow(Q2h,3))*
		      pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		     SI.B0_qmm(Q2h,m2)*(-4 + (16*m4 - 12*m2*Q2h)*pow(l1k,-2) +
		        (-24*m2 - 4*Q2e + 19*Q2h)*pow(l1k,-1) +
		        (64*m4 - 48*m2*Q2h)*pow(2*l1k + Q2e - Q2h,-2) +
		        (48*m2 + 8*Q2e - 38*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-64*m4 + 16*m2*Q2h + 24*pow(Q2h,2)) +
		        (-96*Q2e*pow(l1k,3) - 48*Q2h*pow(l1k,3) + 48*l1k*m2*pow(Q2e,2) +
		           24*l1k*Q2h*pow(Q2e,2) - 96*pow(l1k,2)*pow(Q2e,2) +
		           48*m2*pow(Q2e,3) + 24*pow(Q2e,4))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        pow(l1k,-1)*(12*m2*pow(Q2e,4) - 3*Q2h*pow(Q2e,4) + 6*pow(Q2e,5))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        pow(l1k,-1)*(32*m4*Q2e - 4*m2*pow(Q2e,2) - 7*Q2h*pow(Q2e,2) +
		           pow(Q2e,3))*pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2e,2),-1) + (-32*l1k*m2 + 64*l1k*Q2e + 56*l1k*Q2h +
		           12*Q2e*Q2h + 34*pow(Q2e,2) + 18*pow(Q2h,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*(-18*Q2e*pow(Q2h,2) + 18*pow(Q2h,3))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (-24*m2*pow(Q2e,4) + 6*Q2h*pow(Q2e,4) - 12*pow(Q2e,5))*
		         pow(2*l1k + Q2e - Q2h,-1)*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (-72*Q2e*Q2h*pow(l1k,2) + 96*Q2e*pow(l1k,3) + 48*Q2h*pow(l1k,3) -
		           48*l1k*m2*pow(Q2e,2) - 12*l1k*Q2h*pow(Q2e,2) +
		           24*pow(l1k,2)*pow(Q2e,2) - 24*l1k*pow(Q2e,3) +
		           24*m2*pow(Q2e,3) + 6*Q2h*pow(Q2e,3) + 12*pow(Q2e,4) -
		           72*pow(l1k,2)*pow(Q2h,2) + 36*l1k*pow(Q2h,3) +
		           6*Q2e*pow(Q2h,3) - 6*pow(Q2h,4))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (-64*m4*Q2e + 8*m2*pow(Q2e,2) + 14*Q2h*pow(Q2e,2) -
		           2*pow(Q2e,3))*pow(2*l1k + Q2e - Q2h,-1)*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        (32*l1k*m2 - 64*l1k*Q2e + 16*m2*Q2e - 72*l1k*Q2h + 16*Q2e*Q2h +
		           16*pow(l1k,2) + 8*pow(Q2e,2) + 50*pow(Q2h,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*(9*Q2e*pow(Q2h,2) - 9*pow(Q2h,3))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)) +
		     SI.B0_m0m(m2)*((-8*m4 + 8*m2*Q2h)*pow(l1k,-2) +
		        (4*m2 - (17*Q2e)/2. - 14*Q2h)*pow(l1k,-1) +
		        (-32*m4 + 32*m2*Q2h)*pow(2*l1k + Q2e - Q2h,-2) +
		        (-8*m2 + 17*Q2e + 28*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (32*m4 - 16*m2*Q2h - 16*pow(Q2h,2)) +
		        (-96*m2*Q2e*pow(l1k,2) + 24*Q2e*pow(l1k,3) + 48*pow(l1k,4) -
		           48*l1k*m2*pow(Q2e,2) + 48*m4*pow(Q2e,2) -
		           24*pow(l1k,2)*pow(Q2e,2) - 12*l1k*pow(Q2e,3) +
		           24*m2*pow(Q2e,3) + 3*pow(Q2e,4))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        pow(l1k,-1)*(24*m4*pow(Q2e,3) + 12*m2*pow(Q2e,4) +
		           (3*pow(Q2e,5))/2.)*pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2e,2),-2) + pow(l1k,-1)*
		         (64*pow(m,6) + 8*m2*pow(Q2e,2) + 3*pow(Q2e,3))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (16*l1k*m2 - 64*m4 - 52*l1k*Q2e + 72*m2*Q2e - 24*l1k*Q2h -
		           4*Q2e*Q2h - 72*pow(l1k,2) - 2*pow(Q2e,2) - 12*pow(Q2h,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (32*m4*Q2e - 128*pow(m,6) + 8*Q2e*pow(Q2h,2) - 12*pow(Q2h,3))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (-48*m4*pow(Q2e,3) - 24*m2*pow(Q2e,4) - 3*pow(Q2e,5))*
		         pow(2*l1k + Q2e - Q2h,-1)*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (-96*m2*Q2e*pow(l1k,2) + 12*Q2e*Q2h*pow(l1k,2) -
		           24*Q2e*pow(l1k,3) - 96*Q2h*pow(l1k,3) + 48*pow(l1k,4) -
		           48*l1k*m2*pow(Q2e,2) + 48*m4*pow(Q2e,2) +
		           12*l1k*Q2h*pow(Q2e,2) - 12*pow(l1k,2)*pow(Q2e,2) -
		           6*l1k*pow(Q2e,3) + 24*m2*pow(Q2e,3) + 3*Q2h*pow(Q2e,3) +
		           3*pow(Q2e,4) + 6*l1k*Q2e*pow(Q2h,2) +
		           72*pow(l1k,2)*pow(Q2h,2) - 3*pow(Q2e,2)*pow(Q2h,2) -
		           24*l1k*pow(Q2h,3) - 3*Q2e*pow(Q2h,3) + 3*pow(Q2h,4))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (-128*pow(m,6) - 16*m2*pow(Q2e,2) - 6*pow(Q2e,3))*
		         pow(2*l1k + Q2e - Q2h,-1)*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        (-16*l1k*m2 - 64*m4 + 4*l1k*Q2e + 64*m2*Q2e + 104*l1k*Q2h +
		           12*Q2e*Q2h - 80*pow(l1k,2) + 12*pow(Q2e,2) - 44*pow(Q2h,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*(-16*m4*Q2e + 64*pow(m,6) - 4*Q2e*pow(Q2h,2) +
		           6*pow(Q2h,3))*pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2h,2),-1)) + SI.C0_mMQm0m(-2*l1k + m2,Q2h,m2)*
		      (26*l1k + 16*m2 + Q2e - 26*Q2h +
		        (48*m4 - 8*m2*Q2h - 8*Q2e*Q2h)*pow(l1k,-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-64*m4 - 40*m2*Q2e + 24*m2*Q2h + 9*Q2e*Q2h - 13*pow(Q2e,2) +
		           16*pow(Q2h,2)) + pow(l1k,-2)*
		         (-32*pow(m,6) + 8*m2*pow(Q2h,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (32*m4*Q2h + 64*pow(m,6) - 16*m2*pow(Q2h,2) - 8*pow(Q2h,3)) +
		        (48*m4*pow(Q2e,4) + 36*m2*pow(Q2e,5) - 3*Q2h*pow(Q2e,5) +
		           9*pow(Q2e,6))*pow(2*l1k + Q2e - Q2h,-1)*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (-192*m2*Q2e*pow(l1k,3) + 24*Q2e*Q2h*pow(l1k,3) -
		           48*Q2e*pow(l1k,4) - 192*Q2h*pow(l1k,4) + 96*pow(l1k,5) +
		           96*l1k*m4*pow(Q2e,2) - 96*m2*pow(l1k,2)*pow(Q2e,2) +
		           24*Q2h*pow(l1k,2)*pow(Q2e,2) - 24*pow(l1k,3)*pow(Q2e,2) +
		           72*l1k*m2*pow(Q2e,3) - 48*m4*pow(Q2e,3) -
		           6*l1k*Q2h*pow(Q2e,3) + 18*l1k*pow(Q2e,4) - 36*m2*pow(Q2e,4) -
		           6*Q2h*pow(Q2e,4) - 9*pow(Q2e,5) +
		           12*Q2e*pow(l1k,2)*pow(Q2h,2) + 144*pow(l1k,3)*pow(Q2h,2) -
		           6*l1k*pow(Q2e,2)*pow(Q2h,2) + 3*pow(Q2e,3)*pow(Q2h,2) -
		           6*l1k*Q2e*pow(Q2h,3) - 48*pow(l1k,2)*pow(Q2h,3) +
		           6*l1k*pow(Q2h,4))*pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2h,2),-2) + (128*Q2e*pow(m,6) + 32*m4*pow(Q2e,2) +
		           4*m2*pow(Q2e,3) - 6*Q2h*pow(Q2e,3) + 4*pow(Q2e,4))*
		         pow(2*l1k + Q2e - Q2h,-1)*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        (-64*l1k*m4 + 128*l1k*m2*Q2e - 64*m4*Q2e + 24*l1k*Q2e*Q2h +
		           8*Q2e*pow(l1k,2) + 168*Q2h*pow(l1k,2) - 128*pow(l1k,3) +
		           32*l1k*pow(Q2e,2) - 16*m2*pow(Q2e,2) - 4*Q2h*pow(Q2e,2) -
		           13*pow(Q2e,3) - 72*l1k*pow(Q2h,2) - 8*Q2e*pow(Q2h,2) +
		           10*pow(Q2h,3))*pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2h,2),-1)) + SI.B0_00m(m2)*
		      (32 - 8*m4*pow(l1k,-2) + (4*Q2e - 9*Q2h)*pow(l1k,-1) +
		        (-32*l1k + 8*Q2e - 24*Q2h)*pow(2*l1k - m2,-1) -
		        32*m4*pow(2*l1k + Q2e - Q2h,-2) +
		        (-8*Q2e + 18*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        (32*m4 + 16*m2*Q2h)*pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1) +
		        (-32*l1k - 24*Q2e + 40*Q2h)*pow(2*l1k + m2 + Q2e - Q2h,-1) +
		        pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (32*Q2e*Q2h - 8*pow(Q2e,2) - 16*pow(Q2h,2)) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (16*Q2e*Q2h - 4*pow(Q2e,2) - 8*pow(Q2h,2)) +
		        (-8*l1k*Q2e - 4*l1k*Q2h - 8*Q2e*Q2h - 6*pow(Q2e,2) - 6*pow(Q2h,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*(4*m2*pow(Q2e,2) - 2*Q2h*pow(Q2e,2) - pow(Q2e,3) -
		           4*Q2e*pow(Q2h,2) - 8*pow(Q2h,3))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        8*pow(2*l1k + m2 + Q2e - Q2h,-1)*pow(Q2h,3)*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*(-2*Q2e*pow(Q2h,2) + 2*pow(Q2h,3))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (4*Q2e*pow(Q2h,3) - 8*pow(Q2h,4))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (8*l1k*Q2e + 4*l1k*Q2h - 10*Q2e*Q2h - 2*pow(Q2e,2) - 8*pow(Q2h,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*(Q2e*pow(Q2h,2) - pow(Q2h,3))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) -
		        8*pow(2*l1k - m2,-1)*pow(Q2h,3)*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*m2*pow(Q2e,2) + 4*Q2h*pow(Q2e,2) + 2*pow(Q2e,3) +
		           8*Q2e*pow(Q2h,2) + 16*pow(Q2h,3))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (8*Q2e*pow(Q2h,3) - 16*pow(Q2h,4))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)) +
		     SI.B0_M0m(-2*l1k + m2,m2)*(-21 + (-8*m4 + 8*m2*Q2h)*pow(l1k,-2) +
		        (16*m2 - 2*Q2e - 4*Q2h)*pow(l1k,-1) +
		        (32*l1k - 8*Q2e + 24*Q2h)*pow(2*l1k - m2,-1) +
		        (-40*m2 - 13*Q2e - 8*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (16*m4 - 8*m2*Q2h - 8*pow(Q2h,2)) +
		        pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*Q2e*Q2h + 8*pow(Q2e,2) + 16*pow(Q2h,2)) +
		        (48*m4*pow(Q2e,3) + 48*m2*pow(Q2e,4) - 6*Q2h*pow(Q2e,4) +
		           15*pow(Q2e,5))*pow(2*l1k + Q2e - Q2h,-1)*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (96*m2*Q2e*pow(l1k,2) + 60*Q2e*Q2h*pow(l1k,2) -
		           72*Q2e*pow(l1k,3) + 48*Q2h*pow(l1k,3) - 48*pow(l1k,4) +
		           96*l1k*m2*pow(Q2e,2) - 48*m4*pow(Q2e,2) -
		           12*pow(l1k,2)*pow(Q2e,2) + 30*l1k*pow(Q2e,3) -
		           48*m2*pow(Q2e,3) - 9*Q2h*pow(Q2e,3) - 15*pow(Q2e,4) -
		           6*l1k*Q2e*pow(Q2h,2) + 3*pow(Q2e,2)*pow(Q2h,2) -
		           12*l1k*pow(Q2h,3) - 3*Q2e*pow(Q2h,3) + 3*pow(Q2h,4))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (-16*l1k*m2 + 64*m4 + 52*l1k*Q2e - 80*m2*Q2e - 36*l1k*Q2h -
		           18*Q2e*Q2h + 64*pow(l1k,2) - 18*pow(Q2e,2) + 2*pow(Q2h,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (64*m4*Q2e + 128*pow(m,6) + 16*m2*pow(Q2e,2) -
		           18*Q2h*pow(Q2e,2) + 6*pow(Q2e,3) - 8*Q2e*pow(Q2h,2) -
		           16*pow(Q2h,3))*pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2h,2),-1) + 8*pow(2*l1k - m2,-1)*pow(Q2h,3)*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*(16*m4*Q2e - 64*pow(m,6) - 6*Q2e*pow(Q2h,2) +
		           4*pow(Q2h,3))*pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2h,2),-1) + pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*Q2e*pow(Q2h,3) + 16*pow(Q2h,4))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1))) +
		  G23*(32*M2 + 4*Q2h + pow(l1k,-1)*
		      (8*M2*Q2e + 8*M2*Q2h + (Q2e*Q2h)/2. + 6*Q2h*S + 2*Q2h*U -
		        pow(Q2h,2)/2.) + pow(2*l1k + Q2e - Q2h,-1)*
		      (-16*M2*Q2e - 16*M2*Q2h - Q2e*Q2h - 4*Q2h*S - 12*Q2h*U +
		        pow(Q2h,2)) + pow(2*l1k + Q2e - Q2h,-2)*
		      (-32*m2*M2*Q2h - 16*m2*Q2h*S - 4*m2*pow(Q2h,2) - 4*S*pow(Q2h,2) +
		        pow(Q2h,3) + 16*m2*pow(S,2) + 4*Q2h*pow(S,2)) +
		     pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		      (32*m2*M2*Q2h - 8*m2*Q2h*S - 8*m2*Q2h*U + 16*m2*S*U +
		        12*m2*pow(Q2h,2) + 16*M2*pow(Q2h,2) + 6*S*pow(Q2h,2) +
		        6*U*pow(Q2h,2) + pow(Q2h,3) - 6*Q2h*pow(S,2) - 6*Q2h*pow(U,2)) +
		     pow(l1k,-2)*(-8*m2*M2*Q2h - 4*m2*Q2h*U - m2*pow(Q2h,2) -
		        U*pow(Q2h,2) + pow(Q2h,3)/4. + 4*m2*pow(U,2) + Q2h*pow(U,2)) +
		     (-2*l1k*Q2e*Q2h - 4*l1k*Q2e*S + 12*l1k*Q2h*S - 2*Q2e*Q2h*S -
		        4*l1k*Q2e*U - 4*l1k*Q2h*U - 2*Q2e*Q2h*U + 4*Q2e*S*U + 4*Q2h*S*U -
		        8*Q2h*pow(l1k,2) + 2*Q2h*pow(Q2e,2) - 2*S*pow(Q2e,2) -
		        4*U*pow(Q2e,2) - 6*l1k*pow(Q2h,2) - 2*Q2e*pow(Q2h,2) +
		        12*S*pow(Q2h,2) - 2*U*pow(Q2h,2) - 4*pow(Q2h,3) + 4*Q2e*pow(S,2) -
		        4*Q2h*pow(S,2))*pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) +
		        pow(Q2e,2),-1) + pow(2*l1k + Q2e - Q2h,-1)*
		      (-6*Q2e*S*pow(Q2h,2) + 4*S*U*pow(Q2h,2) + 2*Q2e*pow(Q2h,3) +
		        12*S*pow(Q2h,3) - 2*U*pow(Q2h,3) - 4*pow(Q2h,4) +
		        4*Q2e*Q2h*pow(S,2) - 8*pow(Q2h,2)*pow(S,2))*
		      pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		     pow(l1k,-1)*(-(Q2h*S*pow(Q2e,2)) + 2*S*U*pow(Q2e,2) +
		        (Q2h*pow(Q2e,3))/2. - U*pow(Q2e,3) + 2*Q2e*S*pow(Q2h,2) -
		        (pow(Q2e,2)*pow(Q2h,2))/2. - (Q2e*pow(Q2h,3))/2. + 2*S*pow(Q2h,3) -
		        pow(Q2h,4)/2. - 2*pow(Q2h,2)*pow(S,2))*
		      pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		     pow(2*l1k + Q2e - Q2h,-2)*
		      (4*S*pow(Q2h,4) - pow(Q2h,5) - 4*pow(Q2h,3)*pow(S,2))*
		      pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		     pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		      (2*S*pow(Q2h,4) - pow(Q2h,5)/2. - 2*pow(Q2h,3)*pow(S,2))*
		      pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		     SI.B0_M0m(2*l1k + m2 + Q2e - Q2h,m2)*
		      (8*l1k + 8*m2 - 6*M2 - Q2e + Q2h - 5*S + 3*U +
		        pow(l1k,-1)*(12*m2*M2 + 2*m2*Q2e + 10*M2*Q2e + m2*Q2h -
		           4*M2*Q2h + (5*Q2e*Q2h)/2. - 4*m2*S + (3*Q2e*S)/2. -
		           (3*Q2h*S)/2. + 4*m2*U + (3*Q2e*U)/2. - S*U - pow(Q2e,2) -
		           pow(Q2h,2)/2. - 2*pow(S,2)) +
		        pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (16*l1k*Q2h + 8*Q2e*Q2h - 8*l1k*S - 2*Q2e*S - 22*Q2h*S - 8*l1k*U -
		           6*Q2e*U + 6*Q2h*U + 8*S*U + 4*pow(Q2h,2) + 8*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-24*m2*M2 + 12*M2*Q2e - 6*m2*Q2h + 2*M2*Q2h - 8*m2*S +
		           4*Q2e*S - 2*Q2h*S + 8*Q2h*U - 2*pow(Q2h,2) + 10*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (32*m2*M2*Q2h + 16*m2*Q2h*S + 4*m2*pow(Q2h,2) +
		           4*S*pow(Q2h,2) - pow(Q2h,3) - 16*m2*pow(S,2) - 4*Q2h*pow(S,2)) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-3*Q2e*Q2h*S + 3*Q2e*Q2h*U + 2*Q2e*S*U - 6*Q2h*S*U +
		           Q2h*pow(Q2e,2) - U*pow(Q2e,2) - Q2e*pow(Q2h,2) + S*pow(Q2h,2) +
		           2*U*pow(Q2h,2) + 2*Q2h*pow(S,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*m2*M2*Q2h + 4*m2*Q2h*S + 4*m2*Q2h*U - 8*m2*S*U -
		           6*m2*pow(Q2h,2) - 8*M2*pow(Q2h,2) - 2*S*pow(Q2h,2) -
		           4*U*pow(Q2h,2) - pow(Q2h,3)/2. + 2*Q2h*pow(S,2) + 4*Q2h*pow(U,2)) \
		+ pow(l1k,-1)*(12*m2*S*U*pow(Q2e,3) - 6*Q2h*S*U*pow(Q2e,3) -
		           (15*Q2h*S*pow(Q2e,4))/2. - 6*m2*U*pow(Q2e,4) +
		           3*Q2h*U*pow(Q2e,4) + 9*S*U*pow(Q2e,4) + 3*Q2h*pow(Q2e,5) +
		           (3*S*pow(Q2e,5))/2. - (9*U*pow(Q2e,5))/2. - (3*pow(Q2e,6))/4. +
		           9*S*pow(Q2e,3)*pow(Q2h,2) - 3*pow(Q2e,4)*pow(Q2h,2) +
		           3*Q2h*pow(Q2e,3)*pow(S,2) - 6*pow(Q2e,2)*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        (48*l1k*m2*Q2e*S*U + 24*l1k*Q2e*Q2h*S*U -
		           48*m2*Q2e*S*pow(l1k,2) - 48*m2*Q2e*U*pow(l1k,2) -
		           48*Q2e*Q2h*U*pow(l1k,2) + 48*Q2h*S*U*pow(l1k,2) -
		           48*Q2e*Q2h*pow(l1k,3) - 264*Q2e*S*pow(l1k,3) +
		           72*Q2h*S*pow(l1k,3) + 24*Q2e*U*pow(l1k,3) -
		           48*Q2h*U*pow(l1k,3) - 48*S*U*pow(l1k,3) + 168*Q2e*pow(l1k,4) -
		           72*Q2h*pow(l1k,4) - 144*S*pow(l1k,4) + 48*U*pow(l1k,4) +
		           96*pow(l1k,5) - 48*l1k*m2*S*pow(Q2e,2) -
		           84*l1k*Q2h*S*pow(Q2e,2) - 72*l1k*m2*U*pow(Q2e,2) +
		           72*l1k*S*U*pow(Q2e,2) + 48*m2*S*U*pow(Q2e,2) -
		           12*Q2h*S*U*pow(Q2e,2) + 60*Q2h*pow(l1k,2)*pow(Q2e,2) -
		           168*S*pow(l1k,2)*pow(Q2e,2) - 72*U*pow(l1k,2)*pow(Q2e,2) +
		           84*pow(l1k,3)*pow(Q2e,2) + 72*l1k*Q2h*pow(Q2e,3) -
		           36*l1k*S*pow(Q2e,3) - 12*m2*S*pow(Q2e,3) -
		           48*Q2h*S*pow(Q2e,3) - 84*l1k*U*pow(Q2e,3) -
		           36*m2*U*pow(Q2e,3) + 12*Q2h*U*pow(Q2e,3) + 48*S*U*pow(Q2e,3) -
		           12*pow(l1k,2)*pow(Q2e,3) - 24*l1k*pow(Q2e,4) +
		           (51*Q2h*pow(Q2e,4))/2. + 3*S*pow(Q2e,4) - 33*U*pow(Q2e,4) -
		           (15*pow(Q2e,5))/2. + 84*l1k*Q2e*S*pow(Q2h,2) -
		           42*Q2e*pow(l1k,2)*pow(Q2h,2) + 48*S*pow(l1k,2)*pow(Q2h,2) -
		           12*pow(l1k,3)*pow(Q2h,2) - 45*l1k*pow(Q2e,2)*pow(Q2h,2) +
		           48*S*pow(Q2e,2)*pow(Q2h,2) - (39*pow(Q2e,3)*pow(Q2h,2))/2. -
		           6*l1k*Q2e*pow(Q2h,3) + 12*l1k*S*pow(Q2h,3) +
		           6*Q2e*S*pow(Q2h,3) - 6*pow(l1k,2)*pow(Q2h,3) -
		           (3*pow(Q2e,2)*pow(Q2h,3))/2. - 3*l1k*pow(Q2h,4) -
		           (3*Q2e*pow(Q2h,4))/2. + 6*S*pow(Q2h,4) - (3*pow(Q2h,5))/2. +
		           48*l1k*m2*Q2e*pow(S,2) + 12*l1k*Q2e*Q2h*pow(S,2) +
		           120*Q2e*pow(l1k,2)*pow(S,2) + 48*pow(l1k,3)*pow(S,2) +
		           84*l1k*pow(Q2e,2)*pow(S,2) + 24*m2*pow(Q2e,2)*pow(S,2) +
		           12*Q2h*pow(Q2e,2)*pow(S,2) + 18*pow(Q2e,3)*pow(S,2) -
		           36*l1k*pow(Q2h,2)*pow(S,2) - 30*Q2e*pow(Q2h,2)*pow(S,2) -
		           6*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (6*S*pow(Q2h,5) - (3*pow(Q2h,6))/2. - 6*pow(Q2h,4)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (16*l1k*Q2e*Q2h*S + 8*Q2h*S*pow(Q2e,2) - 8*l1k*Q2e*pow(Q2h,2) -
		           40*l1k*S*pow(Q2h,2) - 36*Q2e*S*pow(Q2h,2) - 8*S*U*pow(Q2h,2) -
		           4*pow(Q2e,2)*pow(Q2h,2) + 24*l1k*pow(Q2h,3) +
		           20*Q2e*pow(Q2h,3) + 36*S*pow(Q2h,3) + 8*U*pow(Q2h,3) -
		           24*pow(Q2h,4) - 8*l1k*Q2e*pow(S,2) + 16*l1k*Q2h*pow(S,2) +
		           16*Q2e*Q2h*pow(S,2) - 4*pow(Q2e,2)*pow(S,2) -
		           12*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (4*M2*Q2e*pow(Q2h,2) + 16*Q2e*S*pow(Q2h,2) - 8*S*U*pow(Q2h,2) -
		           2*M2*pow(Q2h,3) - 4*Q2e*pow(Q2h,3) - 20*S*pow(Q2h,3) +
		           4*U*pow(Q2h,3) + (11*pow(Q2h,4))/2. + 64*m4*pow(S,2) -
		           16*m2*Q2e*pow(S,2) - 16*Q2e*Q2h*pow(S,2) +
		           12*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-4*Q2e*S*U*pow(Q2h,2) - 2*Q2e*S*pow(Q2h,3) +
		           4*Q2e*U*pow(Q2h,3) + 8*S*U*pow(Q2h,3) + 4*S*pow(Q2h,4) -
		           8*U*pow(Q2h,4) + 2*Q2e*pow(Q2h,2)*pow(S,2) -
		           4*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-2*S*pow(Q2h,4) + pow(Q2h,5)/2. + 2*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (-4*S*pow(Q2h,4) + pow(Q2h,5) + 4*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (32*l1k*m2*M2 - 40*l1k*m2*Q2e - 24*l1k*M2*Q2e -
		           40*m2*M2*Q2e - 20*l1k*M2*Q2h - 44*l1k*Q2e*Q2h +
		           10*M2*Q2e*Q2h + 16*l1k*m2*S - 32*m4*S + 44*l1k*Q2e*S +
		           28*m2*Q2e*S + 56*l1k*Q2h*S + 38*Q2e*Q2h*S - 32*l1k*m2*U -
		           32*m4*U + 4*l1k*Q2e*U - 20*m2*Q2e*U - 12*l1k*Q2h*U -
		           8*Q2e*Q2h*U - 16*l1k*S*U + 16*m2*S*U - 24*Q2e*S*U +
		           20*Q2h*S*U - 32*m2*pow(l1k,2) + 56*M2*pow(l1k,2) -
		           46*Q2e*pow(l1k,2) - 18*Q2h*pow(l1k,2) + 72*S*pow(l1k,2) -
		           8*U*pow(l1k,2) - 56*pow(l1k,3) + 3*l1k*pow(Q2e,2) -
		           16*m2*pow(Q2e,2) - 46*M2*pow(Q2e,2) -
		           (57*Q2h*pow(Q2e,2))/2. - 2*S*pow(Q2e,2) + 10*U*pow(Q2e,2) +
		           (19*pow(Q2e,3))/2. - 15*l1k*pow(Q2h,2) - 2*M2*pow(Q2h,2) -
		           (Q2e*pow(Q2h,2))/2. + 22*S*pow(Q2h,2) - 4*U*pow(Q2h,2) -
		           (37*pow(Q2h,3))/2. - 48*l1k*pow(S,2) + 24*m2*pow(S,2) -
		           20*Q2e*pow(S,2) - 18*Q2h*pow(S,2) - 8*l1k*pow(U,2) -
		           8*m2*pow(U,2) - 8*Q2e*pow(U,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*(-16*m4*Q2e*U + 32*m4*S*U + 4*m2*Q2e*S*U +
		           6*Q2e*Q2h*S*U - 12*m2*M2*pow(Q2e,2) + 8*M2*Q2h*pow(Q2e,2) +
		           4*m2*S*pow(Q2e,2) + 9*Q2h*S*pow(Q2e,2) - 2*m2*U*pow(Q2e,2) -
		           Q2h*U*pow(Q2e,2) - 8*S*U*pow(Q2e,2) - 2*m2*pow(Q2e,3) -
		           10*M2*pow(Q2e,3) - (11*Q2h*pow(Q2e,3))/2. - 3*S*pow(Q2e,3) +
		           3*U*pow(Q2e,3) + (7*pow(Q2e,4))/4. + 3*Q2e*S*pow(Q2h,2) -
		           4*Q2e*U*pow(Q2h,2) + 8*S*U*pow(Q2h,2) +
		           (3*pow(Q2e,2)*pow(Q2h,2))/2. + (Q2e*pow(Q2h,3))/2. +
		           2*S*pow(Q2h,3) - 8*U*pow(Q2h,3) + pow(Q2h,4)/2. -
		           5*Q2e*Q2h*pow(S,2) + 2*pow(Q2e,2)*pow(S,2) -
		           6*pow(Q2h,2)*pow(S,2) - 4*m2*Q2e*pow(U,2) -
		           2*pow(Q2e,2)*pow(U,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1)) +
		     SI.C0_mMQm0m(2*l1k + m2 + Q2e - Q2h,Q2h,m2)*
		      (16*l1k*m2 - 20*l1k*M2 - 4*l1k*Q2e + 12*m2*Q2e + 16*M2*Q2e +
		        6*l1k*Q2h + 4*m2*Q2h - 12*M2*Q2h + 6*Q2e*Q2h - 6*l1k*S -
		        16*m2*S + 7*Q2e*S - 16*Q2h*S - 2*l1k*U + 8*m2*U - Q2e*U +
		        5*Q2h*U + 10*S*U + 8*pow(l1k,2) - 7*pow(Q2e,2) + pow(Q2h,2) +
		        8*pow(S,2) + pow(2*l1k + Q2e - Q2h,-1)*
		         (16*m2*M2*Q2e - 16*m4*Q2h + 32*m2*M2*Q2h +
		           16*M2*Q2e*Q2h + 32*m4*S + 16*m2*Q2h*S + 8*Q2e*Q2h*S +
		           16*m2*Q2h*U - 4*m2*pow(Q2h,2) - 12*S*pow(Q2h,2) +
		           8*U*pow(Q2h,2) + pow(Q2h,3) - 8*Q2e*pow(S,2) + 12*Q2h*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (64*m4*M2*Q2h + 64*m4*Q2h*S + 32*m2*M2*pow(Q2h,2) +
		           32*m2*S*pow(Q2h,2) - 64*m4*pow(S,2) - 32*m2*Q2h*pow(S,2)) +
		        2*pow(U,2) + pow(l1k,-1)*
		         (12*m2*M2*Q2e - 8*m4*Q2h - 20*m2*M2*Q2h + 3*m2*Q2e*Q2h -
		           4*M2*Q2e*Q2h - 4*m2*Q2e*S - 4*m2*Q2h*S - Q2e*Q2h*S +
		           16*m4*U + 4*m2*Q2e*U + 8*m2*Q2h*U - (3*Q2e*Q2h*U)/2. -
		           8*m2*S*U + 2*Q2e*S*U + Q2h*S*U + 2*m2*pow(Q2e,2) +
		           7*M2*pow(Q2e,2) + 2*Q2h*pow(Q2e,2) + (5*S*pow(Q2e,2))/2. -
		           (3*pow(Q2e,3))/2. - 7*m2*pow(Q2h,2) - 8*M2*pow(Q2h,2) +
		           (Q2e*pow(Q2h,2))/2. - (9*S*pow(Q2h,2))/2. - 4*U*pow(Q2h,2) -
		           2*Q2e*pow(S,2) + Q2e*pow(U,2) + 4*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*m4*M2*Q2h + 16*m4*Q2h*S + 16*m4*Q2h*U - 32*m4*S*U -
		           16*m2*Q2h*S*U - 16*m4*pow(Q2h,2) - 32*m2*M2*pow(Q2h,2) -
		           8*m2*pow(Q2h,3) - 8*M2*pow(Q2h,3) - 4*S*pow(Q2h,3) -
		           4*U*pow(Q2h,3) + 8*m2*Q2h*pow(S,2) + 4*pow(Q2h,2)*pow(S,2) +
		           8*m2*Q2h*pow(U,2) + 4*pow(Q2h,2)*pow(U,2)) +
		        (96*m2*Q2e*S*U*pow(l1k,2) + 48*Q2e*Q2h*S*U*pow(l1k,2) -
		           96*m2*Q2e*S*pow(l1k,3) - 144*Q2e*Q2h*S*pow(l1k,3) -
		           96*m2*Q2e*U*pow(l1k,3) - 72*Q2e*Q2h*U*pow(l1k,3) -
		           96*Q2e*S*U*pow(l1k,3) + 48*Q2h*S*U*pow(l1k,3) -
		           624*Q2e*S*pow(l1k,4) + 144*Q2e*U*pow(l1k,4) -
		           48*Q2h*U*pow(l1k,4) - 96*S*U*pow(l1k,4) + 432*Q2e*pow(l1k,5) -
		           48*Q2h*pow(l1k,5) - 288*S*pow(l1k,5) + 96*U*pow(l1k,5) +
		           192*pow(l1k,6) + 144*l1k*m2*S*U*pow(Q2e,2) -
		           144*m2*S*pow(l1k,2)*pow(Q2e,2) -
		           264*Q2h*S*pow(l1k,2)*pow(Q2e,2) -
		           192*m2*U*pow(l1k,2)*pow(Q2e,2) -
		           24*Q2h*U*pow(l1k,2)*pow(Q2e,2) + 48*S*U*pow(l1k,2)*pow(Q2e,2) +
		           156*Q2h*pow(l1k,3)*pow(Q2e,2) - 504*S*pow(l1k,3)*pow(Q2e,2) +
		           336*pow(l1k,4)*pow(Q2e,2) - 72*l1k*m2*S*pow(Q2e,3) -
		           180*l1k*Q2h*S*pow(Q2e,3) - 144*l1k*m2*U*pow(Q2e,3) +
		           12*l1k*Q2h*U*pow(Q2e,3) + 96*l1k*S*U*pow(Q2e,3) +
		           72*m2*S*U*pow(Q2e,3) - 12*Q2h*S*U*pow(Q2e,3) +
		           192*Q2h*pow(l1k,2)*pow(Q2e,3) - 168*S*pow(l1k,2)*pow(Q2e,3) -
		           120*U*pow(l1k,2)*pow(Q2e,3) + 60*pow(l1k,3)*pow(Q2e,3) +
		           99*l1k*Q2h*pow(Q2e,4) - 6*l1k*S*pow(Q2e,4) -
		           12*m2*S*pow(Q2e,4) - 54*Q2h*S*pow(Q2e,4) -
		           90*l1k*U*pow(Q2e,4) - 48*m2*U*pow(Q2e,4) +
		           9*Q2h*U*pow(Q2e,4) + 42*S*U*pow(Q2e,4) -
		           60*pow(l1k,2)*pow(Q2e,4) - 39*l1k*pow(Q2e,5) +
		           24*Q2h*pow(Q2e,5) + 9*S*pow(Q2e,5) - 27*U*pow(Q2e,5) -
		           9*pow(Q2e,6) + 72*Q2e*S*pow(l1k,2)*pow(Q2h,2) -
		           24*Q2e*pow(l1k,3)*pow(Q2h,2) + 24*S*pow(l1k,3)*pow(Q2h,2) +
		           72*l1k*S*pow(Q2e,2)*pow(Q2h,2) -
		           48*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2) -
		           36*l1k*pow(Q2e,3)*pow(Q2h,2) + 30*S*pow(Q2e,3)*pow(Q2h,2) -
		           12*pow(Q2e,4)*pow(Q2h,2) + 96*m2*Q2e*pow(l1k,2)*pow(S,2) +
		           96*Q2e*Q2h*pow(l1k,2)*pow(S,2) + 240*Q2e*pow(l1k,3)*pow(S,2) +
		           48*Q2h*pow(l1k,3)*pow(S,2) + 96*pow(l1k,4)*pow(S,2) +
		           96*l1k*m2*pow(Q2e,2)*pow(S,2) +
		           72*l1k*Q2h*pow(Q2e,2)*pow(S,2) +
		           216*pow(l1k,2)*pow(Q2e,2)*pow(S,2) +
		           84*l1k*pow(Q2e,3)*pow(S,2) + 24*m2*pow(Q2e,3)*pow(S,2) +
		           24*Q2h*pow(Q2e,3)*pow(S,2) + 12*pow(Q2e,4)*pow(S,2) -
		           36*l1k*Q2e*pow(Q2h,2)*pow(S,2) -
		           24*pow(l1k,2)*pow(Q2h,2)*pow(S,2) -
		           18*pow(Q2e,2)*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        pow(l1k,-1)*(12*m2*S*U*pow(Q2e,4) - 3*Q2h*S*U*pow(Q2e,4) -
		           6*Q2h*S*pow(Q2e,5) - 6*m2*U*pow(Q2e,5) +
		           (3*Q2h*U*pow(Q2e,5))/2. + 6*S*U*pow(Q2e,5) +
		           (9*Q2h*pow(Q2e,6))/4. + (3*S*pow(Q2e,6))/2. - 3*U*pow(Q2e,6) -
		           (3*pow(Q2e,7))/4. + (9*S*pow(Q2e,4)*pow(Q2h,2))/2. -
		           (3*pow(Q2e,5)*pow(Q2h,2))/2. + 3*Q2h*pow(Q2e,4)*pow(S,2) -
		           3*pow(Q2e,3)*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (4*S*pow(Q2h,4) - pow(Q2h,5) - 4*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (-48*l1k*m2*M2*Q2e - 64*l1k*m4*S + 104*l1k*m2*Q2e*S -
		           32*m4*Q2e*S + 88*l1k*Q2e*Q2h*S - 64*l1k*m4*U -
		           24*l1k*m2*Q2e*U - 64*m4*Q2e*U - 8*l1k*Q2e*Q2h*U +
		           64*m4*S*U - 48*l1k*Q2e*S*U - 8*m2*Q2e*S*U + 16*l1k*Q2h*S*U +
		           12*Q2e*Q2h*S*U + 64*m2*M2*pow(l1k,2) -
		           112*m2*Q2e*pow(l1k,2) + 16*M2*Q2e*pow(l1k,2) -
		           16*M2*Q2h*pow(l1k,2) - 60*Q2e*Q2h*pow(l1k,2) +
		           64*m2*S*pow(l1k,2) + 80*Q2e*S*pow(l1k,2) +
		           64*Q2h*S*pow(l1k,2) - 32*m2*U*pow(l1k,2) -
		           8*Q2h*U*pow(l1k,2) - 16*S*U*pow(l1k,2) - 64*m2*pow(l1k,3) +
		           80*M2*pow(l1k,3) - 76*Q2e*pow(l1k,3) - 12*Q2h*pow(l1k,3) +
		           96*S*pow(l1k,3) - 16*U*pow(l1k,3) - 80*pow(l1k,4) -
		           72*l1k*m2*pow(Q2e,2) - 72*l1k*M2*pow(Q2e,2) -
		           64*m2*M2*pow(Q2e,2) - 68*l1k*Q2h*pow(Q2e,2) +
		           12*M2*Q2h*pow(Q2e,2) - 8*l1k*S*pow(Q2e,2) +
		           44*m2*S*pow(Q2e,2) + 42*Q2h*S*pow(Q2e,2) +
		           24*l1k*U*pow(Q2e,2) - 2*Q2h*U*pow(Q2e,2) - 36*S*U*pow(Q2e,2) +
		           24*pow(l1k,2)*pow(Q2e,2) + 50*l1k*pow(Q2e,3) -
		           20*m2*pow(Q2e,3) - 44*M2*pow(Q2e,3) - 29*Q2h*pow(Q2e,3) -
		           20*S*pow(Q2e,3) + 16*U*pow(Q2e,3) + 19*pow(Q2e,4) +
		           4*l1k*S*pow(Q2h,2) + 2*Q2e*S*pow(Q2h,2) -
		           4*pow(l1k,2)*pow(Q2h,2) + 3*pow(Q2e,2)*pow(Q2h,2) -
		           2*l1k*pow(Q2h,3) - Q2e*pow(Q2h,3) + 4*S*pow(Q2h,3) -
		           pow(Q2h,4) - 16*l1k*m2*pow(S,2) + 64*m4*pow(S,2) -
		           44*l1k*Q2e*pow(S,2) - 24*m2*Q2e*pow(S,2) -
		           20*l1k*Q2h*pow(S,2) - 16*Q2e*Q2h*pow(S,2) -
		           56*pow(l1k,2)*pow(S,2) - 4*pow(Q2e,2)*pow(S,2) -
		           6*pow(Q2h,2)*pow(S,2) - 16*l1k*m2*pow(U,2) -
		           12*l1k*Q2e*pow(U,2) - 16*m2*Q2e*pow(U,2) -
		           8*pow(l1k,2)*pow(U,2) - 6*pow(Q2e,2)*pow(U,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*(32*m4*Q2e*S*U - 16*m4*U*pow(Q2e,2) -
		           4*m2*S*U*pow(Q2e,2) + 2*Q2h*S*U*pow(Q2e,2) -
		           12*m2*M2*pow(Q2e,3) + 4*M2*Q2h*pow(Q2e,3) +
		           4*m2*S*pow(Q2e,3) + 7*Q2h*S*pow(Q2e,3) + 2*m2*U*pow(Q2e,3) -
		           8*S*U*pow(Q2e,3) - 2*m2*pow(Q2e,4) - 7*M2*pow(Q2e,4) -
		           (17*Q2h*pow(Q2e,4))/4. - 4*S*pow(Q2e,4) + 3*U*pow(Q2e,4) +
		           (9*pow(Q2e,5))/4. + pow(Q2e,3)*pow(Q2h,2) -
		           3*Q2h*pow(Q2e,2)*pow(S,2) + 2*pow(Q2e,3)*pow(S,2) -
		           Q2e*pow(Q2h,2)*pow(S,2) - 4*m2*pow(Q2e,2)*pow(U,2) -
		           pow(Q2e,3)*pow(U,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1)) +
		     (-6*l1k*Q2e*Q2h + 4*l1k*Q2e*S + 4*l1k*Q2h*S - 2*Q2e*Q2h*S +
		        4*l1k*Q2e*U - 12*l1k*Q2h*U - 10*Q2e*Q2h*U + 4*Q2e*S*U + 4*Q2h*S*U -
		        8*Q2h*pow(l1k,2) + Q2h*pow(Q2e,2) - 2*S*pow(Q2e,2) +
		        14*l1k*pow(Q2h,2) + 4*Q2e*pow(Q2h,2) - 4*S*pow(Q2h,2) +
		        18*U*pow(Q2h,2) - 9*pow(Q2h,3) + 4*Q2e*pow(U,2) - 4*Q2h*pow(U,2))*
		      pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		     pow(2*l1k + Q2e - Q2h,-1)*
		      (2*Q2h*U*pow(Q2e,2) - 4*S*U*pow(Q2e,2) - Q2h*pow(Q2e,3) +
		        2*S*pow(Q2e,3) - 4*Q2e*U*pow(Q2h,2) + pow(Q2e,2)*pow(Q2h,2) +
		        Q2e*pow(Q2h,3) - 4*U*pow(Q2h,3) + pow(Q2h,4) + 4*pow(Q2h,2)*pow(U,2)\
		)*pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		     pow(l1k,-1)*(3*Q2e*U*pow(Q2h,2) - 2*S*U*pow(Q2h,2) - Q2e*pow(Q2h,3) +
		        S*pow(Q2h,3) - 6*U*pow(Q2h,3) + 2*pow(Q2h,4) - 2*Q2e*Q2h*pow(U,2) +
		        4*pow(Q2h,2)*pow(U,2))*
		      pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		     pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		      (2*U*pow(Q2h,4) - pow(Q2h,5)/2. - 2*pow(Q2h,3)*pow(U,2))*
		      pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		     pow(l1k,-2)*(U*pow(Q2h,4) - pow(Q2h,5)/4. - pow(Q2h,3)*pow(U,2))*
		      pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		     SI.B0_qmm(Q2h,m2)*(16*M2 - 8*Q2e - 12*Q2h + 16*S + 16*U +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-4*M2*Q2e - 20*M2*Q2h + 4*Q2e*Q2h + 4*Q2h*S - 4*Q2e*U -
		           20*Q2h*U + 2*pow(Q2e,2) - 8*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (-48*m2*M2*Q2h - 32*m2*Q2h*S - 4*m2*pow(Q2h,2) -
		           4*S*pow(Q2h,2) + pow(Q2h,3) + 32*m2*pow(S,2) + 4*Q2h*pow(S,2)) +
		        pow(l1k,-1)*(2*M2*Q2e + 10*M2*Q2h - 2*Q2e*Q2h + 2*Q2e*S +
		           10*Q2h*S - 2*Q2h*U - pow(Q2e,2) + 4*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (48*m2*M2*Q2h - 16*m2*Q2h*S - 16*m2*Q2h*U + 32*m2*S*U +
		           20*m2*pow(Q2h,2) + 24*M2*pow(Q2h,2) + 10*S*pow(Q2h,2) +
		           10*U*pow(Q2h,2) + pow(Q2h,3) - 10*Q2h*pow(S,2) - 10*Q2h*pow(U,2)) \
		+ pow(l1k,-2)*(-12*m2*M2*Q2h - 8*m2*Q2h*U - m2*pow(Q2h,2) -
		           U*pow(Q2h,2) + pow(Q2h,3)/4. + 8*m2*pow(U,2) + Q2h*pow(U,2)) +
		        pow(l1k,-1)*(6*Q2h*S*U*pow(Q2e,3) + 3*Q2h*S*pow(Q2e,4) -
		           3*Q2h*U*pow(Q2e,4) - 6*S*U*pow(Q2e,4) - (3*Q2h*pow(Q2e,5))/2. +
		           3*U*pow(Q2e,5) - 9*S*pow(Q2e,3)*pow(Q2h,2) +
		           3*pow(Q2e,4)*pow(Q2h,2) + 6*pow(Q2e,2)*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        (-24*l1k*Q2e*Q2h*S*U - 144*Q2e*Q2h*S*pow(l1k,2) +
		           48*Q2e*Q2h*U*pow(l1k,2) - 48*Q2e*S*U*pow(l1k,2) -
		           48*Q2h*S*U*pow(l1k,2) + 120*Q2e*Q2h*pow(l1k,3) +
		           48*Q2e*S*pow(l1k,3) - 144*Q2h*S*pow(l1k,3) +
		           48*Q2e*U*pow(l1k,3) + 48*Q2h*U*pow(l1k,3) + 96*Q2h*pow(l1k,4) -
		           24*l1k*Q2h*S*pow(Q2e,2) - 72*l1k*S*U*pow(Q2e,2) +
		           12*Q2h*S*U*pow(Q2e,2) + 24*Q2h*pow(l1k,2)*pow(Q2e,2) +
		           72*S*pow(l1k,2)*pow(Q2e,2) + 96*U*pow(l1k,2)*pow(Q2e,2) -
		           24*l1k*Q2h*pow(Q2e,3) + 36*l1k*S*pow(Q2e,3) +
		           12*Q2h*S*pow(Q2e,3) + 72*l1k*U*pow(Q2e,3) -
		           12*Q2h*U*pow(Q2e,3) - 36*S*U*pow(Q2e,3) - 12*Q2h*pow(Q2e,4) +
		           6*S*pow(Q2e,4) + 24*U*pow(Q2e,4) - 108*l1k*Q2e*S*pow(Q2h,2) +
		           60*Q2e*pow(l1k,2)*pow(Q2h,2) - 72*S*pow(l1k,2)*pow(Q2h,2) +
		           24*pow(l1k,3)*pow(Q2h,2) + 54*l1k*pow(Q2e,2)*pow(Q2h,2) -
		           54*S*pow(Q2e,2)*pow(Q2h,2) + 21*pow(Q2e,3)*pow(Q2h,2) +
		           12*l1k*Q2e*pow(Q2h,3) - 24*l1k*S*pow(Q2h,3) -
		           12*Q2e*S*pow(Q2h,3) + 12*pow(l1k,2)*pow(Q2h,3) +
		           3*pow(Q2e,2)*pow(Q2h,3) + 6*l1k*pow(Q2h,4) + 3*Q2e*pow(Q2h,4) -
		           12*S*pow(Q2h,4) + 3*pow(Q2h,5) + 48*l1k*Q2e*Q2h*pow(S,2) -
		           48*Q2e*pow(l1k,2)*pow(S,2) + 48*Q2h*pow(l1k,2)*pow(S,2) -
		           48*l1k*pow(Q2e,2)*pow(S,2) + 12*Q2h*pow(Q2e,2)*pow(S,2) -
		           12*pow(Q2e,3)*pow(S,2) + 48*l1k*pow(Q2h,2)*pow(S,2) +
		           36*Q2e*pow(Q2h,2)*pow(S,2) + 12*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-12*S*pow(Q2h,5) + 3*pow(Q2h,6) + 12*pow(Q2h,4)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*M2*Q2e*pow(Q2h,2) - 26*Q2e*S*pow(Q2h,2) +
		           12*S*U*pow(Q2h,2) + 4*M2*pow(Q2h,3) + 6*Q2e*pow(Q2h,3) +
		           24*S*pow(Q2h,3) - 6*U*pow(Q2h,3) - 6*pow(Q2h,4) +
		           20*Q2e*Q2h*pow(S,2) - 20*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (4*S*pow(Q2h,4) - pow(Q2h,5) - 4*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (2*S*pow(Q2h,4) - pow(Q2h,5)/2. - 2*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (24*l1k*M2*Q2e + 24*l1k*M2*Q2h + 26*l1k*Q2e*Q2h -
		           12*M2*Q2e*Q2h + 32*l1k*m2*S - 72*l1k*Q2e*S + 16*m2*Q2e*S -
		           34*Q2e*Q2h*S + 32*l1k*m2*U - 16*l1k*Q2e*U + 32*m2*Q2e*U -
		           2*Q2e*Q2h*U + 16*l1k*S*U - 32*m2*S*U + 24*Q2e*S*U -
		           8*Q2h*S*U - 32*M2*pow(l1k,2) + 56*Q2e*pow(l1k,2) -
		           8*Q2h*pow(l1k,2) - 48*S*pow(l1k,2) + 32*pow(l1k,3) +
		           36*l1k*pow(Q2e,2) + 32*M2*pow(Q2e,2) + 24*Q2h*pow(Q2e,2) -
		           28*S*pow(Q2e,2) - 16*U*pow(Q2e,2) + 10*pow(Q2e,3) -
		           10*l1k*pow(Q2h,2) + 4*M2*pow(Q2h,2) - 2*Q2e*pow(Q2h,2) +
		           26*S*pow(Q2h,2) - 6*U*pow(Q2h,2) - 6*pow(Q2h,3) +
		           24*l1k*pow(S,2) - 32*m2*pow(S,2) + 24*Q2e*pow(S,2) -
		           8*Q2h*pow(S,2) + 8*l1k*pow(U,2) + 8*Q2e*pow(U,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*(-16*m2*Q2e*S*U - 2*Q2e*Q2h*S*U -
		           8*M2*Q2h*pow(Q2e,2) - 8*Q2h*S*pow(Q2e,2) +
		           8*m2*U*pow(Q2e,2) - Q2h*U*pow(Q2e,2) + 8*S*U*pow(Q2e,2) +
		           6*M2*pow(Q2e,3) + (9*Q2h*pow(Q2e,3))/2. - 2*S*pow(Q2e,3) -
		           4*U*pow(Q2e,3) + pow(Q2e,4) - Q2e*S*pow(Q2h,2) -
		           (3*pow(Q2e,2)*pow(Q2h,2))/2. - (Q2e*pow(Q2h,3))/2. +
		           2*S*pow(Q2h,3) - pow(Q2h,4)/2. + 4*Q2e*Q2h*pow(S,2) +
		           2*pow(Q2h,2)*pow(S,2) + 2*pow(Q2e,2)*pow(U,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-12*Q2h*S*U*pow(Q2e,3) + 6*Q2h*S*pow(Q2e,4) -
		           6*Q2h*U*pow(Q2e,4) + 12*S*U*pow(Q2e,4) + 3*Q2h*pow(Q2e,5) -
		           6*S*pow(Q2e,5) + 18*U*pow(Q2e,3)*pow(Q2h,2) -
		           6*pow(Q2e,4)*pow(Q2h,2) - 12*pow(Q2e,2)*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (24*l1k*Q2e*Q2h*S*U + 48*Q2e*Q2h*S*pow(l1k,2) +
		           144*Q2e*Q2h*U*pow(l1k,2) - 48*Q2e*S*U*pow(l1k,2) -
		           48*Q2h*S*U*pow(l1k,2) + 72*Q2e*Q2h*pow(l1k,3) -
		           48*Q2e*S*pow(l1k,3) - 48*Q2h*S*pow(l1k,3) -
		           48*Q2e*U*pow(l1k,3) + 144*Q2h*U*pow(l1k,3) +
		           96*Q2h*pow(l1k,4) - 12*l1k*Q2h*S*pow(Q2e,2) -
		           12*l1k*Q2h*U*pow(Q2e,2) + 24*l1k*S*U*pow(Q2e,2) -
		           12*Q2h*pow(l1k,2)*pow(Q2e,2) + 24*S*pow(l1k,2)*pow(Q2e,2) +
		           6*l1k*Q2h*pow(Q2e,3) - 12*l1k*S*pow(Q2e,3) +
		           6*Q2h*U*pow(Q2e,3) - 12*S*U*pow(Q2e,3) - 3*Q2h*pow(Q2e,4) +
		           6*S*pow(Q2e,4) - 12*l1k*Q2e*S*pow(Q2h,2) -
		           72*l1k*Q2e*U*pow(Q2h,2) + 48*l1k*S*U*pow(Q2h,2) -
		           84*Q2e*pow(l1k,2)*pow(Q2h,2) + 72*S*pow(l1k,2)*pow(Q2h,2) -
		           288*U*pow(l1k,2)*pow(Q2h,2) - 216*pow(l1k,3)*pow(Q2h,2) -
		           12*U*pow(Q2e,2)*pow(Q2h,2) + 3*pow(Q2e,3)*pow(Q2h,2) +
		           30*l1k*Q2e*pow(Q2h,3) - 36*l1k*S*pow(Q2h,3) +
		           204*l1k*U*pow(Q2h,3) + 6*Q2e*U*pow(Q2h,3) - 12*S*U*pow(Q2h,3) +
		           192*pow(l1k,2)*pow(Q2h,3) + 3*pow(Q2e,2)*pow(Q2h,3) -
		           84*l1k*pow(Q2h,4) - 3*Q2e*pow(Q2h,4) + 6*S*pow(Q2h,4) -
		           60*U*pow(Q2h,4) + 18*pow(Q2h,5) + 48*l1k*Q2e*Q2h*pow(U,2) -
		           48*Q2e*pow(l1k,2)*pow(U,2) + 48*Q2h*pow(l1k,2)*pow(U,2) -
		           96*l1k*pow(Q2h,2)*pow(U,2) + 48*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        pow(l1k,-1)*(6*U*pow(Q2h,5) - (3*pow(Q2h,6))/2. -
		           6*pow(Q2h,4)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (-56*l1k*M2*Q2e + 8*l1k*M2*Q2h - 42*l1k*Q2e*Q2h +
		           4*M2*Q2e*Q2h - 32*l1k*m2*S + 16*l1k*Q2e*S + 16*m2*Q2e*S +
		           16*l1k*Q2h*S - 10*Q2e*Q2h*S - 32*l1k*m2*U + 24*l1k*Q2e*U +
		           64*l1k*Q2h*U - 46*Q2e*Q2h*U - 16*l1k*S*U - 32*m2*S*U +
		           16*Q2e*S*U - 32*M2*pow(l1k,2) + 8*Q2e*pow(l1k,2) +
		           40*Q2h*pow(l1k,2) - 16*S*pow(l1k,2) - 64*U*pow(l1k,2) -
		           32*pow(l1k,3) - 4*l1k*pow(Q2e,2) + 12*M2*pow(Q2e,2) +
		           11*Q2h*pow(Q2e,2) - 8*S*pow(Q2e,2) - 4*U*pow(Q2e,2) +
		           2*pow(Q2e,3) - 6*l1k*pow(Q2h,2) + 8*M2*pow(Q2h,2) +
		           22*Q2e*pow(Q2h,2) - 10*S*pow(Q2h,2) + 10*U*pow(Q2h,2) -
		           9*pow(Q2h,3) - 8*l1k*pow(S,2) + 4*Q2e*pow(S,2) +
		           4*Q2h*pow(S,2) - 24*l1k*pow(U,2) - 32*m2*pow(U,2) +
		           12*Q2e*pow(U,2) + 4*Q2h*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (32*m2*Q2e*S*U + 4*Q2e*Q2h*S*U + 16*M2*Q2h*pow(Q2e,2) -
		           16*m2*S*pow(Q2e,2) + 2*Q2h*S*pow(Q2e,2) +
		           16*Q2h*U*pow(Q2e,2) - 16*S*U*pow(Q2e,2) - 12*M2*pow(Q2e,3) -
		           9*Q2h*pow(Q2e,3) + 8*S*pow(Q2e,3) + 4*U*pow(Q2e,3) -
		           2*pow(Q2e,4) + 2*Q2e*U*pow(Q2h,2) + 3*pow(Q2e,2)*pow(Q2h,2) +
		           Q2e*pow(Q2h,3) - 4*U*pow(Q2h,3) + pow(Q2h,4) -
		           4*pow(Q2e,2)*pow(S,2) - 8*Q2e*Q2h*pow(U,2) -
		           4*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*(4*M2*Q2e*pow(Q2h,2) + 13*Q2e*U*pow(Q2h,2) -
		           6*S*U*pow(Q2h,2) - 2*M2*pow(Q2h,3) - 3*Q2e*pow(Q2h,3) +
		           3*S*pow(Q2h,3) - 12*U*pow(Q2h,3) + 3*pow(Q2h,4) -
		           10*Q2e*Q2h*pow(U,2) + 10*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (2*U*pow(Q2h,4) - pow(Q2h,5)/2. - 2*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-2)*(U*pow(Q2h,4) - pow(Q2h,5)/4. - pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)) +
		     SI.C0_mMQm0m(-2*l1k + m2,Q2h,m2)*
		      (-16*l1k*m2 + 4*l1k*M2 + 16*l1k*Q2e + 4*m2*Q2e + 12*M2*Q2e -
		        2*l1k*Q2h + 12*m2*Q2h - 18*M2*Q2h - 4*Q2e*Q2h + 2*l1k*S + Q2e*S +
		        4*Q2h*S - 18*l1k*U - 24*m2*U + 15*Q2e*U - 13*Q2h*U + 10*S*U -
		        16*pow(l1k,2) - 4*pow(Q2e,2) + 4*pow(Q2h,2) + 6*pow(U,2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-24*m2*M2*Q2e + 16*m4*Q2h + 40*m2*M2*Q2h - 6*m2*Q2e*Q2h +
		           8*M2*Q2e*Q2h - 32*m4*S - 8*m2*Q2e*S - 16*m2*Q2h*S +
		           3*Q2e*Q2h*S + 8*m2*Q2e*U + 8*m2*Q2h*U + 2*Q2e*Q2h*U +
		           16*m2*S*U - 4*Q2e*S*U - 2*Q2h*S*U - 4*m2*pow(Q2e,2) -
		           14*M2*pow(Q2e,2) - 4*Q2h*pow(Q2e,2) - 5*U*pow(Q2e,2) +
		           3*pow(Q2e,3) + 14*m2*pow(Q2h,2) + 16*M2*pow(Q2h,2) -
		           Q2e*pow(Q2h,2) + 8*S*pow(Q2h,2) + 9*U*pow(Q2h,2) -
		           2*Q2e*pow(S,2) - 8*Q2h*pow(S,2) + 4*Q2e*pow(U,2)) +
		        pow(l1k,-1)*(-8*m2*M2*Q2e + 8*m4*Q2h - 16*m2*M2*Q2h -
		           8*M2*Q2e*Q2h - 8*m2*Q2h*S - 16*m4*U - 8*m2*Q2h*U -
		           4*Q2e*Q2h*U + 2*m2*pow(Q2h,2) - 4*S*pow(Q2h,2) +
		           6*U*pow(Q2h,2) - pow(Q2h,3)/2. + 4*Q2e*pow(U,2) - 6*Q2h*pow(U,2)) \
		+ pow(l1k,-2)*(16*m4*M2*Q2h + 16*m4*Q2h*U + 8*m2*M2*pow(Q2h,2) +
		           8*m2*U*pow(Q2h,2) - 16*m4*pow(U,2) - 8*m2*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*m4*M2*Q2h + 16*m4*Q2h*S + 16*m4*Q2h*U - 32*m4*S*U -
		           16*m2*Q2h*S*U - 16*m4*pow(Q2h,2) - 32*m2*M2*pow(Q2h,2) -
		           8*m2*pow(Q2h,3) - 8*M2*pow(Q2h,3) - 4*S*pow(Q2h,3) -
		           4*U*pow(Q2h,3) + 8*m2*Q2h*pow(S,2) + 4*pow(Q2h,2)*pow(S,2) +
		           8*m2*Q2h*pow(U,2) + 4*pow(Q2h,2)*pow(U,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-24*m2*S*U*pow(Q2e,4) + 6*Q2h*S*U*pow(Q2e,4) +
		           12*m2*S*pow(Q2e,5) - 3*Q2h*S*pow(Q2e,5) +
		           12*Q2h*U*pow(Q2e,5) - 12*S*U*pow(Q2e,5) -
		           (9*Q2h*pow(Q2e,6))/2. + 6*S*pow(Q2e,6) - 3*U*pow(Q2e,6) +
		           (3*pow(Q2e,7))/2. - 9*U*pow(Q2e,4)*pow(Q2h,2) +
		           3*pow(Q2e,5)*pow(Q2h,2) - 6*Q2h*pow(Q2e,4)*pow(U,2) +
		           6*pow(Q2e,3)*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (96*m2*Q2e*S*U*pow(l1k,2) + 96*m2*Q2e*S*pow(l1k,3) -
		           48*Q2e*Q2h*S*pow(l1k,3) + 96*m2*Q2e*U*pow(l1k,3) -
		           264*Q2e*Q2h*U*pow(l1k,3) + 144*Q2h*S*U*pow(l1k,3) -
		           240*Q2e*Q2h*pow(l1k,4) + 48*Q2e*S*pow(l1k,4) +
		           192*Q2h*S*pow(l1k,4) + 240*Q2e*U*pow(l1k,4) -
		           720*Q2h*U*pow(l1k,4) - 96*S*U*pow(l1k,4) + 144*Q2e*pow(l1k,5) -
		           528*Q2h*pow(l1k,5) - 96*S*pow(l1k,5) + 288*U*pow(l1k,5) +
		           192*pow(l1k,6) - 48*l1k*m2*S*U*pow(Q2e,2) -
		           12*l1k*Q2h*S*U*pow(Q2e,2) - 48*m2*S*pow(l1k,2)*pow(Q2e,2) -
		           12*Q2h*U*pow(l1k,2)*pow(Q2e,2) + 24*S*U*pow(l1k,2)*pow(Q2e,2) +
		           12*Q2h*pow(l1k,3)*pow(Q2e,2) - 24*U*pow(l1k,3)*pow(Q2e,2) -
		           24*pow(l1k,4)*pow(Q2e,2) + 24*l1k*m2*S*pow(Q2e,3) +
		           6*l1k*Q2h*S*pow(Q2e,3) + 12*l1k*Q2h*U*pow(Q2e,3) -
		           24*l1k*S*U*pow(Q2e,3) + 24*m2*S*U*pow(Q2e,3) +
		           6*Q2h*S*U*pow(Q2e,3) - 12*S*pow(l1k,2)*pow(Q2e,3) +
		           12*U*pow(l1k,2)*pow(Q2e,3) + 12*pow(l1k,3)*pow(Q2e,3) -
		           3*l1k*Q2h*pow(Q2e,4) + 12*l1k*S*pow(Q2e,4) -
		           12*m2*S*pow(Q2e,4) - 3*Q2h*S*pow(Q2e,4) - 6*l1k*U*pow(Q2e,4) -
		           9*Q2h*U*pow(Q2e,4) + 12*S*U*pow(Q2e,4) -
		           6*pow(l1k,2)*pow(Q2e,4) + 3*l1k*pow(Q2e,5) + 3*Q2h*pow(Q2e,5) -
		           6*S*pow(Q2e,5) + 3*U*pow(Q2e,5) - (3*pow(Q2e,6))/2. +
		           12*Q2e*S*pow(l1k,2)*pow(Q2h,2) +
		           84*Q2e*U*pow(l1k,2)*pow(Q2h,2) - 72*S*U*pow(l1k,2)*pow(Q2h,2) +
		           144*Q2e*pow(l1k,3)*pow(Q2h,2) - 144*S*pow(l1k,3)*pow(Q2h,2) +
		           696*U*pow(l1k,3)*pow(Q2h,2) + 600*pow(l1k,4)*pow(Q2h,2) +
		           12*l1k*U*pow(Q2e,2)*pow(Q2h,2) +
		           6*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2) -
		           3*l1k*pow(Q2e,3)*pow(Q2h,2) - 6*l1k*Q2e*U*pow(Q2h,3) +
		           12*l1k*S*U*pow(Q2h,3) - 36*Q2e*pow(l1k,2)*pow(Q2h,3) +
		           48*S*pow(l1k,2)*pow(Q2h,3) - 324*U*pow(l1k,2)*pow(Q2h,3) -
		           360*pow(l1k,3)*pow(Q2h,3) - 3*l1k*pow(Q2e,2)*pow(Q2h,3) +
		           3*l1k*Q2e*pow(Q2h,4) - 6*l1k*S*pow(Q2h,4) +
		           72*l1k*U*pow(Q2h,4) + 120*pow(l1k,2)*pow(Q2h,4) -
		           21*l1k*pow(Q2h,5) - 6*U*pow(Q2h,5) + (3*pow(Q2h,6))/2. +
		           96*m2*Q2e*pow(l1k,2)*pow(U,2) -
		           24*Q2e*Q2h*pow(l1k,2)*pow(U,2) + 48*Q2e*pow(l1k,3)*pow(U,2) -
		           240*Q2h*pow(l1k,3)*pow(U,2) + 96*pow(l1k,4)*pow(U,2) -
		           12*l1k*Q2h*pow(Q2e,2)*pow(U,2) + 6*Q2h*pow(Q2e,3)*pow(U,2) +
		           192*pow(l1k,2)*pow(Q2h,2)*pow(U,2) -
		           60*l1k*pow(Q2h,3)*pow(U,2) + 6*pow(Q2h,4)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (112*l1k*m2*M2*Q2e + 32*l1k*M2*Q2e*Q2h + 64*l1k*m4*S -
		           8*l1k*m2*Q2e*S - 32*m4*Q2e*S - 2*l1k*Q2e*Q2h*S +
		           64*l1k*m4*U - 40*l1k*m2*Q2e*U + 78*l1k*Q2e*Q2h*U +
		           64*m4*S*U + 8*l1k*Q2e*S*U - 8*m2*Q2e*S*U - 6*Q2e*Q2h*S*U +
		           64*m2*M2*pow(l1k,2) - 16*m2*Q2e*pow(l1k,2) -
		           48*M2*Q2e*pow(l1k,2) + 24*M2*Q2h*pow(l1k,2) +
		           104*Q2e*Q2h*pow(l1k,2) - 16*Q2e*S*pow(l1k,2) -
		           32*Q2h*S*pow(l1k,2) + 96*m2*U*pow(l1k,2) -
		           120*Q2e*U*pow(l1k,2) + 88*Q2h*U*pow(l1k,2) -
		           16*S*U*pow(l1k,2) + 64*m2*pow(l1k,3) - 16*M2*pow(l1k,3) -
		           100*Q2e*pow(l1k,3) + 28*Q2h*pow(l1k,3) + 16*S*pow(l1k,3) +
		           16*pow(l1k,4) + 8*l1k*m2*pow(Q2e,2) + 28*l1k*M2*pow(Q2e,2) -
		           24*m2*M2*pow(Q2e,2) - 5*l1k*Q2h*pow(Q2e,2) -
		           6*M2*Q2h*pow(Q2e,2) + 4*m2*S*pow(Q2e,2) +
		           3*Q2h*S*pow(Q2e,2) + 16*l1k*U*pow(Q2e,2) + 8*m2*U*pow(Q2e,2) +
		           6*Q2h*U*pow(Q2e,2) - 10*S*U*pow(Q2e,2) +
		           22*pow(l1k,2)*pow(Q2e,2) - 9*l1k*pow(Q2e,3) -
		           4*m2*pow(Q2e,3) - 14*M2*pow(Q2e,3) - 4*Q2h*pow(Q2e,3) +
		           3*S*pow(Q2e,3) - 8*U*pow(Q2e,3) + (9*pow(Q2e,4))/2. -
		           12*l1k*M2*pow(Q2h,2) - 33*l1k*Q2e*pow(Q2h,2) -
		           4*M2*Q2e*pow(Q2h,2) + 20*l1k*S*pow(Q2h,2) +
		           2*Q2e*S*pow(Q2h,2) - 92*l1k*U*pow(Q2h,2) -
		           12*Q2e*U*pow(Q2h,2) + 4*S*U*pow(Q2h,2) -
		           62*pow(l1k,2)*pow(Q2h,2) - pow(Q2e,2)*pow(Q2h,2) +
		           35*l1k*pow(Q2h,3) + 2*M2*pow(Q2h,3) + 3*Q2e*pow(Q2h,3) -
		           4*S*pow(Q2h,3) + 28*U*pow(Q2h,3) - (15*pow(Q2h,4))/2. +
		           16*l1k*m2*pow(S,2) + 4*l1k*Q2e*pow(S,2) - 8*m2*Q2e*pow(S,2) -
		           2*Q2e*Q2h*pow(S,2) - 2*pow(Q2e,2)*pow(S,2) +
		           16*l1k*m2*pow(U,2) + 64*m4*pow(U,2) - 36*l1k*Q2e*pow(U,2) -
		           16*m2*Q2e*pow(U,2) + 68*l1k*Q2h*pow(U,2) +
		           6*Q2e*Q2h*pow(U,2) - 48*pow(l1k,2)*pow(U,2) +
		           4*pow(Q2e,2)*pow(U,2) - 28*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-64*m4*Q2e*S*U + 32*m4*S*pow(Q2e,2) + 8*m2*S*U*pow(Q2e,2) -
		           4*Q2h*S*U*pow(Q2e,2) + 24*m2*M2*pow(Q2e,3) -
		           8*M2*Q2h*pow(Q2e,3) - 4*m2*S*pow(Q2e,3) -
		           8*m2*U*pow(Q2e,3) - 14*Q2h*U*pow(Q2e,3) + 16*S*U*pow(Q2e,3) +
		           4*m2*pow(Q2e,4) + 14*M2*pow(Q2e,4) +
		           (17*Q2h*pow(Q2e,4))/2. - 6*S*pow(Q2e,4) + 8*U*pow(Q2e,4) -
		           (9*pow(Q2e,5))/2. - 2*pow(Q2e,3)*pow(Q2h,2) +
		           8*m2*pow(Q2e,2)*pow(S,2) + 2*pow(Q2e,3)*pow(S,2) +
		           6*Q2h*pow(Q2e,2)*pow(U,2) - 4*pow(Q2e,3)*pow(U,2) +
		           2*Q2e*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*(-2*U*pow(Q2h,4) + pow(Q2h,5)/2. +
		           2*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)) +
		     SI.B0_m0m(m2)*(-8*l1k - 16*m2 - 24*M2 + 15*Q2e - 4*Q2h - 2*S - 8*U +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (32*m2*M2*Q2h + 16*m2*Q2h*S + 4*m2*pow(Q2h,2) +
		           4*S*pow(Q2h,2) - pow(Q2h,3) - 16*m2*pow(S,2) - 4*Q2h*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (48*m2*M2 + 4*m2*Q2e + 20*M2*Q2e + 8*m2*Q2h + 18*M2*Q2h +
		           4*Q2e*Q2h + 16*m2*S - 5*Q2e*S + 10*Q2h*S - 8*m2*U + 7*Q2e*U +
		           9*Q2h*U + 6*S*U - 4*pow(Q2e,2) - 2*pow(Q2h,2) - 10*pow(S,2) -
		           4*pow(U,2)) + pow(l1k,-1)*
		         (-24*m2*M2 - 2*m2*Q2e - 10*M2*Q2e - 4*m2*Q2h - 9*M2*Q2h -
		           2*Q2e*Q2h + 4*m2*S - (7*Q2e*S)/2. - (9*Q2h*S)/2. - 8*m2*U +
		           (5*Q2e*U)/2. - 5*Q2h*U - 3*S*U + 2*pow(Q2e,2) + pow(Q2h,2) +
		           2*pow(S,2) + 5*pow(U,2)) +
		        pow(l1k,-2)*(8*m2*M2*Q2h + 4*m2*Q2h*U + m2*pow(Q2h,2) +
		           U*pow(Q2h,2) - pow(Q2h,3)/4. - 4*m2*pow(U,2) - Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*m2*M2*Q2h + 8*m2*Q2h*S + 8*m2*Q2h*U - 16*m2*S*U -
		           12*m2*pow(Q2h,2) - 16*M2*pow(Q2h,2) - 6*S*pow(Q2h,2) -
		           6*U*pow(Q2h,2) - pow(Q2h,3) + 6*Q2h*pow(S,2) + 6*Q2h*pow(U,2)) +
		        pow(l1k,-1)*(-12*m2*S*U*pow(Q2e,3) + (9*Q2h*S*pow(Q2e,4))/2. +
		           6*m2*U*pow(Q2e,4) - 3*S*U*pow(Q2e,4) - (3*Q2h*pow(Q2e,5))/2. -
		           (3*S*pow(Q2e,5))/2. + (3*U*pow(Q2e,5))/2. + (3*pow(Q2e,6))/4. -
		           3*Q2h*pow(Q2e,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        (-48*l1k*m2*Q2e*S*U + 48*m2*Q2e*S*pow(l1k,2) +
		           144*Q2e*Q2h*S*pow(l1k,2) + 48*m2*Q2e*U*pow(l1k,2) +
		           48*Q2e*S*U*pow(l1k,2) - 72*Q2e*Q2h*pow(l1k,3) +
		           216*Q2e*S*pow(l1k,3) + 72*Q2h*S*pow(l1k,3) -
		           72*Q2e*U*pow(l1k,3) + 48*S*U*pow(l1k,3) - 168*Q2e*pow(l1k,4) -
		           24*Q2h*pow(l1k,4) + 144*S*pow(l1k,4) - 48*U*pow(l1k,4) -
		           96*pow(l1k,5) + 48*l1k*m2*S*pow(Q2e,2) +
		           108*l1k*Q2h*S*pow(Q2e,2) + 72*l1k*m2*U*pow(Q2e,2) -
		           48*m2*S*U*pow(Q2e,2) - 84*Q2h*pow(l1k,2)*pow(Q2e,2) +
		           96*S*pow(l1k,2)*pow(Q2e,2) - 24*U*pow(l1k,2)*pow(Q2e,2) -
		           84*pow(l1k,3)*pow(Q2e,2) - 48*l1k*Q2h*pow(Q2e,3) +
		           12*m2*S*pow(Q2e,3) + 36*Q2h*S*pow(Q2e,3) +
		           12*l1k*U*pow(Q2e,3) + 36*m2*U*pow(Q2e,3) - 12*S*U*pow(Q2e,3) +
		           12*pow(l1k,2)*pow(Q2e,3) + 24*l1k*pow(Q2e,4) -
		           (27*Q2h*pow(Q2e,4))/2. - 9*S*pow(Q2e,4) + 9*U*pow(Q2e,4) +
		           (15*pow(Q2e,5))/2. + 24*l1k*Q2e*S*pow(Q2h,2) -
		           18*Q2e*pow(l1k,2)*pow(Q2h,2) + 24*S*pow(l1k,2)*pow(Q2h,2) -
		           12*pow(l1k,3)*pow(Q2h,2) - 9*l1k*pow(Q2e,2)*pow(Q2h,2) +
		           6*S*pow(Q2e,2)*pow(Q2h,2) - (3*pow(Q2e,3)*pow(Q2h,2))/2. -
		           6*l1k*Q2e*pow(Q2h,3) + 12*l1k*S*pow(Q2h,3) +
		           6*Q2e*S*pow(Q2h,3) - 6*pow(l1k,2)*pow(Q2h,3) -
		           (3*pow(Q2e,2)*pow(Q2h,3))/2. - 3*l1k*pow(Q2h,4) -
		           (3*Q2e*pow(Q2h,4))/2. + 6*S*pow(Q2h,4) - (3*pow(Q2h,5))/2. -
		           48*l1k*m2*Q2e*pow(S,2) - 60*l1k*Q2e*Q2h*pow(S,2) -
		           72*Q2e*pow(l1k,2)*pow(S,2) - 48*Q2h*pow(l1k,2)*pow(S,2) -
		           48*pow(l1k,3)*pow(S,2) - 36*l1k*pow(Q2e,2)*pow(S,2) -
		           24*m2*pow(Q2e,2)*pow(S,2) - 24*Q2h*pow(Q2e,2)*pow(S,2) -
		           6*pow(Q2e,3)*pow(S,2) - 12*l1k*pow(Q2h,2)*pow(S,2) -
		           6*Q2e*pow(Q2h,2)*pow(S,2) - 6*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (6*S*pow(Q2h,5) - (3*pow(Q2h,6))/2. - 6*pow(Q2h,4)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (4*M2*Q2e*pow(Q2h,2) + 16*Q2e*S*pow(Q2h,2) - 8*S*U*pow(Q2h,2) -
		           2*M2*pow(Q2h,3) - 4*Q2e*pow(Q2h,3) - 24*S*pow(Q2h,3) +
		           4*U*pow(Q2h,3) + (13*pow(Q2h,4))/2. - 64*m4*pow(S,2) +
		           16*m2*Q2e*pow(S,2) - 8*Q2e*Q2h*pow(S,2) +
		           24*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-2*S*pow(Q2h,4) + pow(Q2h,5)/2. + 2*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (-4*S*pow(Q2h,4) + pow(Q2h,5) + 4*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (-32*l1k*m2*M2 + 40*l1k*m2*Q2e + 40*m2*M2*Q2e -
		           4*l1k*M2*Q2h + 16*l1k*Q2e*Q2h + 2*M2*Q2e*Q2h -
		           48*l1k*m2*S + 32*m4*S + 24*l1k*Q2e*S - 44*m2*Q2e*S -
		           44*l1k*Q2h*S - 10*Q2e*Q2h*S + 32*m4*U + 8*l1k*Q2e*U -
		           12*m2*Q2e*U + 8*l1k*Q2h*U + 4*Q2e*Q2h*U + 16*m2*S*U +
		           4*Q2e*S*U - 8*Q2h*S*U + 32*m2*pow(l1k,2) -
		           24*M2*pow(l1k,2) - 10*Q2e*pow(l1k,2) + 18*Q2h*pow(l1k,2) -
		           24*S*pow(l1k,2) + 8*U*pow(l1k,2) + 24*pow(l1k,3) -
		           39*l1k*pow(Q2e,2) + 16*m2*pow(Q2e,2) + 14*M2*pow(Q2e,2) +
		           (13*Q2h*pow(Q2e,2))/2. + 28*S*pow(Q2e,2) + 2*U*pow(Q2e,2) -
		           (39*pow(Q2e,3))/2. + 11*l1k*pow(Q2h,2) - 2*M2*pow(Q2h,2) +
		           (5*Q2e*pow(Q2h,2))/2. - 24*S*pow(Q2h,2) + 4*U*pow(Q2h,2) +
		           (13*pow(Q2h,3))/2. + 24*l1k*pow(S,2) + 8*m2*pow(S,2) +
		           22*Q2h*pow(S,2) + 8*m2*pow(U,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*(16*m4*Q2e*U - 32*m4*S*U + 12*m2*Q2e*S*U +
		           12*m2*M2*pow(Q2e,2) - 4*m2*S*pow(Q2e,2) -
		           2*Q2h*S*pow(Q2e,2) - 6*m2*U*pow(Q2e,2) + 2*S*U*pow(Q2e,2) +
		           2*m2*pow(Q2e,3) + 4*M2*pow(Q2e,3) + (3*Q2h*pow(Q2e,3))/2. +
		           5*S*pow(Q2e,3) - (11*pow(Q2e,4))/4. - 2*Q2e*S*pow(Q2h,2) +
		           (pow(Q2e,2)*pow(Q2h,2))/2. + (Q2e*pow(Q2h,3))/2. -
		           2*S*pow(Q2h,3) + pow(Q2h,4)/2. + Q2e*Q2h*pow(S,2) -
		           2*pow(Q2e,2)*pow(S,2) + 2*pow(Q2h,2)*pow(S,2) +
		           4*m2*Q2e*pow(U,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (24*m2*S*U*pow(Q2e,3) - 12*m2*S*pow(Q2e,4) -
		           9*Q2h*U*pow(Q2e,4) + 6*S*U*pow(Q2e,4) + 3*Q2h*pow(Q2e,5) -
		           3*S*pow(Q2e,5) + 3*U*pow(Q2e,5) - (3*pow(Q2e,6))/2. +
		           6*Q2h*pow(Q2e,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (48*l1k*m2*Q2e*S*U + 48*m2*Q2e*S*pow(l1k,2) -
		           24*Q2e*Q2h*S*pow(l1k,2) + 48*m2*Q2e*U*pow(l1k,2) -
		           132*Q2e*Q2h*U*pow(l1k,2) + 72*Q2h*S*U*pow(l1k,2) -
		           120*Q2e*Q2h*pow(l1k,3) + 24*Q2e*S*pow(l1k,3) +
		           96*Q2h*S*pow(l1k,3) + 120*Q2e*U*pow(l1k,3) -
		           360*Q2h*U*pow(l1k,3) - 48*S*U*pow(l1k,3) + 72*Q2e*pow(l1k,4) -
		           264*Q2h*pow(l1k,4) - 48*S*pow(l1k,4) + 144*U*pow(l1k,4) +
		           96*pow(l1k,5) - 24*l1k*m2*S*pow(Q2e,2) -
		           6*l1k*Q2h*U*pow(Q2e,2) + 12*l1k*S*U*pow(Q2e,2) -
		           24*m2*S*U*pow(Q2e,2) - 6*Q2h*S*U*pow(Q2e,2) +
		           6*Q2h*pow(l1k,2)*pow(Q2e,2) - 12*U*pow(l1k,2)*pow(Q2e,2) -
		           12*pow(l1k,3)*pow(Q2e,2) - 6*l1k*S*pow(Q2e,3) +
		           12*m2*S*pow(Q2e,3) + 3*Q2h*S*pow(Q2e,3) + 6*l1k*U*pow(Q2e,3) +
		           6*Q2h*U*pow(Q2e,3) - 6*S*U*pow(Q2e,3) +
		           6*pow(l1k,2)*pow(Q2e,3) - 3*l1k*pow(Q2e,4) -
		           (3*Q2h*pow(Q2e,4))/2. + 3*S*pow(Q2e,4) - 3*U*pow(Q2e,4) +
		           (3*pow(Q2e,5))/2. + 6*l1k*Q2e*S*pow(Q2h,2) +
		           42*l1k*Q2e*U*pow(Q2h,2) - 36*l1k*S*U*pow(Q2h,2) +
		           72*Q2e*pow(l1k,2)*pow(Q2h,2) - 72*S*pow(l1k,2)*pow(Q2h,2) +
		           348*U*pow(l1k,2)*pow(Q2h,2) + 300*pow(l1k,3)*pow(Q2h,2) +
		           3*l1k*pow(Q2e,2)*pow(Q2h,2) + 6*U*pow(Q2e,2)*pow(Q2h,2) -
		           (3*pow(Q2e,3)*pow(Q2h,2))/2. - 18*l1k*Q2e*pow(Q2h,3) +
		           24*l1k*S*pow(Q2h,3) - 162*l1k*U*pow(Q2h,3) -
		           3*Q2e*U*pow(Q2h,3) + 6*S*U*pow(Q2h,3) -
		           180*pow(l1k,2)*pow(Q2h,3) - (3*pow(Q2e,2)*pow(Q2h,3))/2. +
		           60*l1k*pow(Q2h,4) + (3*Q2e*pow(Q2h,4))/2. - 3*S*pow(Q2h,4) +
		           36*U*pow(Q2h,4) - (21*pow(Q2h,5))/2. +
		           48*l1k*m2*Q2e*pow(U,2) - 12*l1k*Q2e*Q2h*pow(U,2) +
		           24*Q2e*pow(l1k,2)*pow(U,2) - 120*Q2h*pow(l1k,2)*pow(U,2) +
		           48*pow(l1k,3)*pow(U,2) - 6*Q2h*pow(Q2e,2)*pow(U,2) +
		           96*l1k*pow(Q2h,2)*pow(U,2) - 30*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        pow(l1k,-1)*(-3*U*pow(Q2h,5) + (3*pow(Q2h,6))/4. +
		           3*pow(Q2h,4)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (32*l1k*m2*M2 - 8*l1k*m2*Q2e - 24*l1k*M2*Q2e +
		           56*m2*M2*Q2e + 12*l1k*M2*Q2h + 52*l1k*Q2e*Q2h +
		           16*M2*Q2e*Q2h + 32*m4*S - 12*l1k*Q2e*S - 12*m2*Q2e*S -
		           16*l1k*Q2h*S + 3*Q2e*Q2h*S + 48*l1k*m2*U + 32*m4*U -
		           60*l1k*Q2e*U - 20*m2*Q2e*U + 44*l1k*Q2h*U + 39*Q2e*Q2h*U +
		           16*m2*S*U - 2*Q2e*S*U - 8*Q2h*S*U + 32*m2*pow(l1k,2) -
		           8*M2*pow(l1k,2) - 50*Q2e*pow(l1k,2) + 14*Q2h*pow(l1k,2) +
		           8*S*pow(l1k,2) + 8*pow(l1k,3) + 11*l1k*pow(Q2e,2) +
		           4*m2*pow(Q2e,2) + 8*M2*pow(Q2e,2) - (5*Q2h*pow(Q2e,2))/2. +
		           3*S*pow(Q2e,2) + 10*U*pow(Q2e,2) - (11*pow(Q2e,3))/2. -
		           31*l1k*pow(Q2h,2) - 6*M2*pow(Q2h,2) -
		           (33*Q2e*pow(Q2h,2))/2. + 10*S*pow(Q2h,2) - 46*U*pow(Q2h,2) +
		           (35*pow(Q2h,3))/2. + 8*m2*pow(S,2) - 24*l1k*pow(U,2) +
		           8*m2*pow(U,2) - 18*Q2e*pow(U,2) + 34*Q2h*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*(-2*M2*Q2e*pow(Q2h,2) - 8*Q2e*U*pow(Q2h,2) +
		           4*S*U*pow(Q2h,2) + M2*pow(Q2h,3) + 2*Q2e*pow(Q2h,3) -
		           2*S*pow(Q2h,3) + 12*U*pow(Q2h,3) - (13*pow(Q2h,4))/4. +
		           32*m4*pow(U,2) - 8*m2*Q2e*pow(U,2) + 4*Q2e*Q2h*pow(U,2) -
		           12*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*m4*Q2e*S + 64*m4*S*U - 24*m2*Q2e*S*U -
		           24*m2*M2*pow(Q2e,2) + 12*m2*S*pow(Q2e,2) +
		           8*m2*U*pow(Q2e,2) + 4*Q2h*U*pow(Q2e,2) - 4*S*U*pow(Q2e,2) -
		           4*m2*pow(Q2e,3) - 8*M2*pow(Q2e,3) - 3*Q2h*pow(Q2e,3) -
		           10*U*pow(Q2e,3) + (11*pow(Q2e,4))/2. + 4*Q2e*U*pow(Q2h,2) -
		           pow(Q2e,2)*pow(Q2h,2) - Q2e*pow(Q2h,3) + 4*U*pow(Q2h,3) -
		           pow(Q2h,4) - 8*m2*Q2e*pow(S,2) - 2*Q2e*Q2h*pow(U,2) +
		           4*pow(Q2e,2)*pow(U,2) - 4*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-2)*(-(U*pow(Q2h,4)) + pow(Q2h,5)/4. + pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-2*U*pow(Q2h,4) + pow(Q2h,5)/2. + 2*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)) +
		     SI.B0_M0m(-2*l1k + m2,m2)*(8*m2 - 2*M2 - 6*Q2e + 3*Q2h - S - 3*U +
		        pow(2*l1k - m2,-1)*(16*l1k*Q2h - 8*l1k*S + 2*Q2e*S - 2*Q2h*S -
		           8*l1k*U - 2*Q2e*U + 26*Q2h*U - 8*S*U - 12*pow(Q2h,2) - 8*pow(U,2)\
		) + pow(l1k,-1)*(12*m2*M2 - 6*M2*Q2e + 3*m2*Q2h - M2*Q2h -
		           4*Q2h*S + 4*m2*U - 2*Q2e*U + Q2h*U + pow(Q2h,2) - 5*pow(U,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-24*m2*M2 - 4*m2*Q2e - 20*M2*Q2e - 2*m2*Q2h + 8*M2*Q2h -
		           5*Q2e*Q2h - 8*m2*S - 3*Q2e*S + 8*m2*U - 3*Q2e*U + 3*Q2h*U +
		           2*S*U + 2*pow(Q2e,2) + pow(Q2h,2) + 4*pow(U,2)) +
		        pow(l1k,-2)*(8*m2*M2*Q2h + 4*m2*Q2h*U + m2*pow(Q2h,2) +
		           U*pow(Q2h,2) - pow(Q2h,3)/4. - 4*m2*pow(U,2) - Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*m2*M2*Q2h + 4*m2*Q2h*S + 4*m2*Q2h*U - 8*m2*S*U -
		           6*m2*pow(Q2h,2) - 8*M2*pow(Q2h,2) - 4*S*pow(Q2h,2) -
		           2*U*pow(Q2h,2) - pow(Q2h,3)/2. + 4*Q2h*pow(S,2) + 2*Q2h*pow(U,2)) \
		+ pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (6*Q2e*Q2h*S - 6*Q2e*Q2h*U + 4*Q2e*S*U - 12*Q2h*S*U +
		           2*Q2h*pow(Q2e,2) - 2*S*pow(Q2e,2) - 2*Q2e*pow(Q2h,2) +
		           4*S*pow(Q2h,2) + 2*U*pow(Q2h,2) + 4*Q2h*pow(U,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-24*m2*S*U*pow(Q2e,3) + 12*Q2h*S*U*pow(Q2e,3) +
		           12*m2*S*pow(Q2e,4) - 6*Q2h*S*pow(Q2e,4) +
		           15*Q2h*U*pow(Q2e,4) - 18*S*U*pow(Q2e,4) - 6*Q2h*pow(Q2e,5) +
		           9*S*pow(Q2e,5) - 3*U*pow(Q2e,5) + (3*pow(Q2e,6))/2. -
		           18*U*pow(Q2e,3)*pow(Q2h,2) + 6*pow(Q2e,4)*pow(Q2h,2) -
		           6*Q2h*pow(Q2e,3)*pow(U,2) + 12*pow(Q2e,2)*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (-48*l1k*m2*Q2e*S*U - 24*l1k*Q2e*Q2h*S*U -
		           48*m2*Q2e*S*pow(l1k,2) - 24*Q2e*Q2h*S*pow(l1k,2) -
		           48*m2*Q2e*U*pow(l1k,2) - 12*Q2e*Q2h*U*pow(l1k,2) +
		           48*Q2e*S*U*pow(l1k,2) - 24*Q2h*S*U*pow(l1k,2) +
		           48*Q2e*Q2h*pow(l1k,3) + 24*Q2e*S*pow(l1k,3) -
		           48*Q2h*S*pow(l1k,3) - 72*Q2e*U*pow(l1k,3) +
		           216*Q2h*U*pow(l1k,3) + 48*S*U*pow(l1k,3) - 72*Q2e*pow(l1k,4) +
		           168*Q2h*pow(l1k,4) + 48*S*pow(l1k,4) - 144*U*pow(l1k,4) -
		           96*pow(l1k,5) + 24*l1k*m2*S*pow(Q2e,2) +
		           12*l1k*Q2h*S*pow(Q2e,2) + 18*l1k*Q2h*U*pow(Q2e,2) -
		           36*l1k*S*U*pow(Q2e,2) + 24*m2*S*U*pow(Q2e,2) +
		           6*Q2h*S*U*pow(Q2e,2) + 6*Q2h*pow(l1k,2)*pow(Q2e,2) -
		           24*S*pow(l1k,2)*pow(Q2e,2) + 12*U*pow(l1k,2)*pow(Q2e,2) +
		           12*pow(l1k,3)*pow(Q2e,2) - 6*l1k*Q2h*pow(Q2e,3) +
		           18*l1k*S*pow(Q2e,3) - 12*m2*S*pow(Q2e,3) -
		           3*Q2h*S*pow(Q2e,3) - 6*l1k*U*pow(Q2e,3) - 12*Q2h*U*pow(Q2e,3) +
		           18*S*U*pow(Q2e,3) - 6*pow(l1k,2)*pow(Q2e,3) +
		           3*l1k*pow(Q2e,4) + (9*Q2h*pow(Q2e,4))/2. - 9*S*pow(Q2e,4) +
		           3*U*pow(Q2e,4) - (3*pow(Q2e,5))/2. + 6*l1k*Q2e*S*pow(Q2h,2) +
		           30*l1k*Q2e*U*pow(Q2h,2) - 12*l1k*S*U*pow(Q2h,2) +
		           12*Q2e*pow(l1k,2)*pow(Q2h,2) - 60*U*pow(l1k,2)*pow(Q2h,2) -
		           84*pow(l1k,3)*pow(Q2h,2) - 3*l1k*pow(Q2e,2)*pow(Q2h,2) +
		           6*U*pow(Q2e,2)*pow(Q2h,2) - (3*pow(Q2e,3)*pow(Q2h,2))/2. -
		           12*l1k*Q2e*pow(Q2h,3) + 12*l1k*S*pow(Q2h,3) -
		           42*l1k*U*pow(Q2h,3) - 3*Q2e*U*pow(Q2h,3) + 6*S*U*pow(Q2h,3) -
		           12*pow(l1k,2)*pow(Q2h,3) - (3*pow(Q2e,2)*pow(Q2h,3))/2. +
		           24*l1k*pow(Q2h,4) + (3*Q2e*pow(Q2h,4))/2. - 3*S*pow(Q2h,4) +
		           24*U*pow(Q2h,4) - (15*pow(Q2h,5))/2. -
		           48*l1k*m2*Q2e*pow(U,2) - 36*l1k*Q2e*Q2h*pow(U,2) +
		           24*Q2e*pow(l1k,2)*pow(U,2) + 72*Q2h*pow(l1k,2)*pow(U,2) -
		           48*pow(l1k,3)*pow(U,2) + 6*Q2h*pow(Q2e,2)*pow(U,2) -
		           18*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        pow(l1k,-1)*(-3*U*pow(Q2h,5) + (3*pow(Q2h,6))/4. +
		           3*pow(Q2h,4)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (-32*l1k*m2*M2 + 8*l1k*m2*Q2e + 80*l1k*M2*Q2e -
		           56*m2*M2*Q2e - 20*l1k*M2*Q2h - 16*l1k*Q2e*Q2h -
		           20*M2*Q2e*Q2h + 32*l1k*m2*S - 32*m4*S - 4*m2*Q2e*S +
		           4*l1k*Q2h*S + Q2e*Q2h*S - 16*l1k*m2*U - 32*m4*U +
		           40*l1k*Q2e*U + 20*m2*Q2e*U - 120*l1k*Q2h*U - 7*Q2e*Q2h*U +
		           16*l1k*S*U + 16*m2*S*U - 10*Q2e*S*U + 12*Q2h*S*U -
		           32*m2*pow(l1k,2) + 40*M2*pow(l1k,2) + 42*Q2e*pow(l1k,2) -
		           62*Q2h*pow(l1k,2) + 8*S*pow(l1k,2) + 64*U*pow(l1k,2) +
		           24*pow(l1k,3) - 7*l1k*pow(Q2e,2) - 4*m2*pow(Q2e,2) -
		           20*M2*pow(Q2e,2) - (15*Q2h*pow(Q2e,2))/2. + 3*S*pow(Q2e,2) -
		           6*U*pow(Q2e,2) + (7*pow(Q2e,3))/2. + 59*l1k*pow(Q2h,2) -
		           2*M2*pow(Q2h,2) + (9*Q2e*pow(Q2h,2))/2. - 8*S*pow(Q2h,2) +
		           66*U*pow(Q2h,2) - (71*pow(Q2h,3))/2. + 8*l1k*pow(S,2) -
		           8*m2*pow(S,2) - 4*Q2e*pow(S,2) - 4*Q2h*pow(S,2) +
		           48*l1k*pow(U,2) + 24*m2*pow(U,2) + 10*Q2e*pow(U,2) -
		           42*Q2h*pow(U,2))*pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2h,2),-1) + pow(l1k,-1)*
		         (-2*M2*Q2e*pow(Q2h,2) - 8*Q2e*U*pow(Q2h,2) + 4*S*U*pow(Q2h,2) +
		           M2*pow(Q2h,3) + 2*Q2e*pow(Q2h,3) - 2*S*pow(Q2h,3) +
		           10*U*pow(Q2h,3) - (11*pow(Q2h,4))/4. - 32*m4*pow(U,2) +
		           8*m2*Q2e*pow(U,2) + 8*Q2e*Q2h*pow(U,2) - 6*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k - m2,-1)*(16*l1k*Q2e*Q2h*U - 8*l1k*Q2e*pow(Q2h,2) -
		           40*l1k*U*pow(Q2h,2) + 8*Q2e*U*pow(Q2h,2) + 8*S*U*pow(Q2h,2) +
		           24*l1k*pow(Q2h,3) - 4*Q2e*pow(Q2h,3) - 8*S*pow(Q2h,3) -
		           16*U*pow(Q2h,3) + 12*pow(Q2h,4) - 8*l1k*Q2e*pow(U,2) +
		           16*l1k*Q2h*pow(U,2) - 4*Q2e*Q2h*pow(U,2) + 4*pow(Q2h,2)*pow(U,2)\
		)*pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (32*m4*Q2e*S - 64*m4*S*U - 8*m2*Q2e*S*U - 12*Q2e*Q2h*S*U +
		           24*m2*M2*pow(Q2e,2) - 16*M2*Q2h*pow(Q2e,2) +
		           4*m2*S*pow(Q2e,2) + 2*Q2h*S*pow(Q2e,2) - 8*m2*U*pow(Q2e,2) -
		           18*Q2h*U*pow(Q2e,2) + 16*S*U*pow(Q2e,2) + 4*m2*pow(Q2e,3) +
		           20*M2*pow(Q2e,3) + 11*Q2h*pow(Q2e,3) - 6*S*pow(Q2e,3) +
		           6*U*pow(Q2e,3) - (7*pow(Q2e,4))/2. + 8*Q2e*S*pow(Q2h,2) -
		           6*Q2e*U*pow(Q2h,2) - 16*S*U*pow(Q2h,2) -
		           3*pow(Q2e,2)*pow(Q2h,2) - Q2e*pow(Q2h,3) + 16*S*pow(Q2h,3) -
		           4*U*pow(Q2h,3) - pow(Q2h,4) + 8*m2*Q2e*pow(S,2) +
		           4*pow(Q2e,2)*pow(S,2) + 10*Q2e*Q2h*pow(U,2) -
		           4*pow(Q2e,2)*pow(U,2) + 12*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*Q2e*S*U*pow(Q2h,2) + 8*Q2e*S*pow(Q2h,3) -
		           4*Q2e*U*pow(Q2h,3) + 16*S*U*pow(Q2h,3) - 16*S*pow(Q2h,4) +
		           8*U*pow(Q2h,4) + 4*Q2e*pow(Q2h,2)*pow(U,2) -
		           8*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-2)*(-(U*pow(Q2h,4)) + pow(Q2h,5)/4. + pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-2*U*pow(Q2h,4) + pow(Q2h,5)/2. + 2*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)) +
		     SI.B0_00m(m2)*(12*Q2h - 8*S - 8*U +
		        pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-16*l1k*Q2h - 8*Q2e*Q2h + 8*l1k*S + 2*Q2e*S + 22*Q2h*S + 8*l1k*U +
		           6*Q2e*U - 6*Q2h*U - 8*S*U - 4*pow(Q2h,2) - 8*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-3*Q2e*Q2h + 4*Q2e*S - 12*Q2h*S + 8*Q2h*U - 8*S*U + 3*pow(Q2h,2) +
		           8*pow(S,2)) + pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (3*Q2e*Q2h*S - 3*Q2e*Q2h*U - 2*Q2e*S*U + 6*Q2h*S*U -
		           Q2h*pow(Q2e,2) + U*pow(Q2e,2) + Q2e*pow(Q2h,2) - S*pow(Q2h,2) -
		           2*U*pow(Q2h,2) - 2*Q2h*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (16*m2*Q2h*S - 4*m2*pow(Q2h,2) - 4*S*pow(Q2h,2) + pow(Q2h,3) -
		           16*m2*pow(S,2) + 4*Q2h*pow(S,2)) +
		        pow(l1k,-1)*((3*Q2e*Q2h)/2. - 4*Q2h*S - 2*Q2e*U + 6*Q2h*U + 4*S*U -
		           (3*pow(Q2h,2))/2. - 4*pow(U,2)) +
		        pow(2*l1k - m2,-1)*(-16*l1k*Q2h + 8*l1k*S - 2*Q2e*S + 2*Q2h*S +
		           8*l1k*U + 2*Q2e*U - 26*Q2h*U + 8*S*U + 12*pow(Q2h,2) + 8*pow(U,2)) \
		+ pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-6*Q2e*Q2h*S + 6*Q2e*Q2h*U - 4*Q2e*S*U + 12*Q2h*S*U -
		           2*Q2h*pow(Q2e,2) + 2*S*pow(Q2e,2) + 2*Q2e*pow(Q2h,2) -
		           4*S*pow(Q2h,2) - 2*U*pow(Q2h,2) - 4*Q2h*pow(U,2)) +
		        pow(l1k,-2)*(4*m2*Q2h*U - m2*pow(Q2h,2) - U*pow(Q2h,2) +
		           pow(Q2h,3)/4. - 4*m2*pow(U,2) + Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (8*m2*Q2h*S + 8*m2*Q2h*U - 16*m2*S*U - 4*m2*pow(Q2h,2) -
		           2*S*pow(Q2h,2) - 2*U*pow(Q2h,2) + pow(Q2h,3) + 2*Q2h*pow(S,2) +
		           2*Q2h*pow(U,2)) + (2*l1k*Q2e*Q2h + 4*l1k*Q2e*S - 12*l1k*Q2h*S +
		           6*Q2e*Q2h*S + 4*l1k*Q2e*U + 4*l1k*Q2h*U + 6*Q2e*Q2h*U -
		           4*Q2e*S*U - 4*Q2h*S*U + 8*Q2h*pow(l1k,2) - 2*Q2h*pow(Q2e,2) +
		           2*S*pow(Q2e,2) + 4*U*pow(Q2e,2) + 14*l1k*pow(Q2h,2) -
		           24*S*pow(Q2h,2) + 6*U*pow(Q2h,2) + 18*pow(Q2h,3) -
		           4*Q2e*pow(S,2) + 4*Q2h*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-6*Q2e*S*pow(Q2h,2) + 4*S*U*pow(Q2h,2) + 2*Q2e*pow(Q2h,3) +
		           20*S*pow(Q2h,3) - 2*U*pow(Q2h,3) - 6*pow(Q2h,4) +
		           4*Q2e*Q2h*pow(S,2) - 16*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*(-4*Q2e*Q2h*S*U + Q2h*S*pow(Q2e,2) +
		           2*Q2h*U*pow(Q2e,2) - 2*S*U*pow(Q2e,2) - (Q2h*pow(Q2e,3))/2. +
		           U*pow(Q2e,3) + 4*Q2e*U*pow(Q2h,2) - 8*S*U*pow(Q2h,2) -
		           (pow(Q2e,2)*pow(Q2h,2))/2. - (Q2e*pow(Q2h,3))/2. -
		           2*S*pow(Q2h,3) + 8*U*pow(Q2h,3) - pow(Q2h,4)/2. +
		           2*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-16*l1k*Q2e*Q2h*S - 8*Q2h*S*pow(Q2e,2) + 8*l1k*Q2e*pow(Q2h,2) +
		           40*l1k*S*pow(Q2h,2) + 36*Q2e*S*pow(Q2h,2) + 8*S*U*pow(Q2h,2) +
		           4*pow(Q2e,2)*pow(Q2h,2) - 24*l1k*pow(Q2h,3) -
		           20*Q2e*pow(Q2h,3) - 36*S*pow(Q2h,3) - 8*U*pow(Q2h,3) +
		           24*pow(Q2h,4) + 8*l1k*Q2e*pow(S,2) - 16*l1k*Q2h*pow(S,2) -
		           16*Q2e*Q2h*pow(S,2) + 4*pow(Q2e,2)*pow(S,2) +
		           12*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (4*S*pow(Q2h,4) - pow(Q2h,5) - 4*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (2*S*pow(Q2h,4) - pow(Q2h,5)/2. - 2*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (4*Q2e*S*U*pow(Q2h,2) + 2*Q2e*S*pow(Q2h,3) - 4*Q2e*U*pow(Q2h,3) -
		           8*S*U*pow(Q2h,3) - 4*S*pow(Q2h,4) + 8*U*pow(Q2h,4) -
		           2*Q2e*pow(Q2h,2)*pow(S,2) + 4*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (6*l1k*Q2e*Q2h - 4*l1k*Q2e*S - 4*l1k*Q2h*S + 6*Q2e*Q2h*S -
		           4*l1k*Q2e*U + 12*l1k*Q2h*U + 14*Q2e*Q2h*U - 4*Q2e*S*U -
		           4*Q2h*S*U + 8*Q2h*pow(l1k,2) - Q2h*pow(Q2e,2) + 2*S*pow(Q2e,2) -
		           22*l1k*pow(Q2h,2) - 10*Q2e*pow(Q2h,2) + 8*S*pow(Q2h,2) -
		           30*U*pow(Q2h,2) + 27*pow(Q2h,3) - 4*Q2e*pow(U,2) + 4*Q2h*pow(U,2)\
		)*pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (8*Q2e*Q2h*S*U - 4*Q2h*S*pow(Q2e,2) - 2*Q2h*U*pow(Q2e,2) +
		           4*S*U*pow(Q2e,2) + Q2h*pow(Q2e,3) - 2*S*pow(Q2e,3) -
		           8*Q2e*S*pow(Q2h,2) + 16*S*U*pow(Q2h,2) + pow(Q2e,2)*pow(Q2h,2) +
		           Q2e*pow(Q2h,3) - 16*S*pow(Q2h,3) + 4*U*pow(Q2h,3) + pow(Q2h,4) -
		           4*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k - m2,-1)*(-16*l1k*Q2e*Q2h*U + 8*l1k*Q2e*pow(Q2h,2) +
		           40*l1k*U*pow(Q2h,2) - 8*Q2e*U*pow(Q2h,2) - 8*S*U*pow(Q2h,2) -
		           24*l1k*pow(Q2h,3) + 4*Q2e*pow(Q2h,3) + 8*S*pow(Q2h,3) +
		           16*U*pow(Q2h,3) - 12*pow(Q2h,4) + 8*l1k*Q2e*pow(U,2) -
		           16*l1k*Q2h*pow(U,2) + 4*Q2e*Q2h*pow(U,2) - 4*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*(3*Q2e*U*pow(Q2h,2) - 2*S*U*pow(Q2h,2) -
		           Q2e*pow(Q2h,3) + S*pow(Q2h,3) - 10*U*pow(Q2h,3) + 3*pow(Q2h,4) -
		           2*Q2e*Q2h*pow(U,2) + 8*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (2*U*pow(Q2h,4) - pow(Q2h,5)/2. - 2*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-2)*(U*pow(Q2h,4) - pow(Q2h,5)/4. - pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (8*Q2e*S*U*pow(Q2h,2) - 8*Q2e*S*pow(Q2h,3) + 4*Q2e*U*pow(Q2h,3) -
		           16*S*U*pow(Q2h,3) + 16*S*pow(Q2h,4) - 8*U*pow(Q2h,4) -
		           4*Q2e*pow(Q2h,2)*pow(U,2) + 8*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)));

	melem_interf = c*(diag_vb_on_shell/M1gamma(Q2h, Q2e, l1k, S, U, G1, G23) - deltaZ1);

	return 2.*melem_interf;

}

long double Gamma_Loop::vb_off_shell(const long double Q2h, const long double Q2e, const long double l1k,
		const long double S, const long double U, const long double G1, const long double G23)const {

	long double diag_vb_off_shell = G1*(-32 + 16*m2*pow(l1k,-1) +
		     SI.C0_m0M0mm(2*l1k + m2 + Q2e - Q2h,m2)*
		      (-16*l1k + 16*m2 - 16*Q2e +
		        pow(l1k,-1)*(32*m4 + 8*m2*Q2e - 8*m2*Q2h - 8*pow(Q2e,2))) -
		     32*Q2h*pow(Q2e - Q2h,-1) - 32*m2*pow(2*l1k + Q2e - Q2h,-1) +
		     SI.C0_mmqm0m(Q2e,m2)*
		      (-64*m2 - 16*Q2e + pow(l1k,-1)*
		         (-32*m4 + 8*m2*Q2e + 16*m2*Q2h + 4*Q2e*Q2h + 4*pow(Q2e,2)) +
		        (64*m4 - 16*m2*Q2e - 32*m2*Q2h - 8*Q2e*Q2h - 8*pow(Q2e,2))*
		         pow(2*l1k + Q2e - Q2h,-1)) +
		     SI.C0_m0M0mm(-2*l1k + m2,m2)*
		      (16*l1k + 16*m2 - 8*Q2e - 8*Q2h +
		        (-64*m4 - 16*m2*Q2e + 16*m2*Q2h + 16*pow(Q2e,2))*
		         pow(2*l1k + Q2e - Q2h,-1)) +
		     SI.B0_qmm(Q2e,m2)*((-4*Q2e + 4*Q2h)*pow(l1k,-1) +
		        32*Q2e*pow(4*m2 + Q2e,-1) +
		        pow(l1k,-1)*(-4*Q2e*Q2h + 4*pow(Q2e,2))*pow(4*m2 + Q2e,-1) +
		        (8*Q2e - 8*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        16*m2*Q2h*pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1) +
		        (8*Q2e*Q2h - 8*pow(Q2e,2))*pow(4*m2 + Q2e,-1)*
		         pow(2*l1k + Q2e - Q2h,-1) - 32*pow(Q2e - Q2h,-2)*pow(Q2h,2)) +
		     SI.C0_0qQmmm(Q2e,Q2h,m2)*
		      (16*Q2e - 16*Q2h - 64*m2*Q2h*pow(Q2e - Q2h,-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*(64*m4 - 24*pow(Q2e,2) - 8*pow(Q2h,2)) +
		        pow(l1k,-1)*(-32*m4 + 12*pow(Q2e,2) + 4*pow(Q2h,2))) +
		     SI.B0_00m(m2)*(-32 + (16*m2 - 8*Q2e + 12*Q2h)*pow(l1k,-1) +
		        (32*l1k - 16*Q2e + 24*Q2h)*pow(2*l1k - m2,-1) +
		        (-32*m2 + 16*Q2e - 24*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        (32*l1k + 32*Q2e - 40*Q2h)*pow(2*l1k + m2 + Q2e - Q2h,-1) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-20*Q2e*Q2h + 8*pow(Q2e,2) + 12*pow(Q2h,2)) +
		        pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-40*Q2e*Q2h + 16*pow(Q2e,2) + 24*pow(Q2h,2))) +
		     SI.D0_mm0QqMm0mm(Q2h,Q2e,-2*l1k + m2,m2)*
		      (64*m4 + 16*l1k*Q2e - 16*m2*Q2e - 32*m2*Q2h - 8*Q2e*Q2h +
		        pow(l1k,-1)*(32*m4*Q2e - 32*m4*Q2h - 16*m2*Q2e*Q2h +
		           64*pow(m,6)) - 8*pow(Q2e,2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-192*m4*Q2e + 128*m4*Q2h + 48*m2*Q2e*Q2h - 128*pow(m,6) +
		           16*pow(Q2e,3) - 16*m2*pow(Q2h,2))) +
		     SI.D0_mm0QqMm0mm(Q2h,Q2e,2*l1k + m2 + Q2e - Q2h,m2)*
		      (64*m4 - 16*l1k*Q2e - 16*m2*Q2e - 32*m2*Q2h - 16*pow(Q2e,2) +
		        (-64*m4*Q2e + 64*m4*Q2h + 32*m2*Q2e*Q2h - 128*pow(m,6))*
		         pow(2*l1k + Q2e - Q2h,-1) +
		        pow(l1k,-1)*(96*m4*Q2e - 64*m4*Q2h - 24*m2*Q2e*Q2h +
		           64*pow(m,6) - 8*pow(Q2e,3) + 8*m2*pow(Q2h,2))) +
		     SI.C0_mMQm0m(2*l1k + m2 + Q2e - Q2h,Q2h,m2)*
		      (-12*l1k + 8*m2 - 12*Q2e + 8*Q2h +
		        pow(l1k,-1)*(48*m4 - 12*m2*Q2h - 4*Q2e*Q2h - 15*pow(Q2e,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*m4 + 16*m2*Q2e - 8*Q2e*Q2h + 8*pow(Q2h,2)) +
		        (32*l1k*m2*Q2e - 32*m4*Q2e + 16*l1k*Q2e*Q2h +
		           32*m2*pow(l1k,2) - 32*Q2e*pow(l1k,2) - 16*pow(l1k,3) +
		           8*l1k*pow(Q2e,2) + 24*m2*pow(Q2e,2) + 24*pow(Q2e,3))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*(-16*m4*pow(Q2e,2) + 8*m2*pow(Q2e,3) -
		           4*Q2h*pow(Q2e,3) + 7*pow(Q2e,4))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1)) +
		     SI.B0_M0m(2*l1k + m2 + Q2e - Q2h,m2)*
		      ((-4*m2 - 2*Q2e - 8*Q2h)*pow(l1k,-1) +
		        (32*m4 - 16*m2*Q2e - 16*m2*Q2h)*pow(2*l1k + Q2e - Q2h,-2) +
		        (32*m2 - 12*Q2e - 20*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        (-32*l1k - 32*Q2e + 40*Q2h)*pow(2*l1k + m2 + Q2e - Q2h,-1) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (20*Q2e*Q2h - 8*pow(Q2e,2) - 12*pow(Q2h,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*m4 + 8*m2*Q2h + 8*pow(Q2h,2)) +
		        pow(l1k,-1)*(-16*m4*Q2e + 4*m2*pow(Q2e,2) - 8*Q2h*pow(Q2e,2) +
		           10*pow(Q2e,3))*pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2e,2),-1) + (16*l1k*m2 + 8*l1k*Q2h + 16*Q2e*Q2h +
		           20*pow(Q2e,2) + 4*pow(Q2h,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*(-4*Q2e*pow(Q2h,2) + 4*pow(Q2h,3))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1)) +
		     SI.C0_mMQm0m(-2*l1k + m2,Q2h,m2)*
		      (4*l1k + 8*m2 - 6*Q2e + 4*Q2h +
		        (-96*m4 + 24*m2*Q2h + 8*Q2e*Q2h + 30*pow(Q2e,2))*
		         pow(2*l1k + Q2e - Q2h,-1) +
		        pow(l1k,-1)*(16*m4 - 8*m2*Q2e + 4*Q2e*Q2h - 4*pow(Q2h,2)) +
		        (32*m4*pow(Q2e,2) - 16*m2*pow(Q2e,3) + 8*Q2h*pow(Q2e,3) -
		           14*pow(Q2e,4))*pow(2*l1k + Q2e - Q2h,-1)*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        (-32*m4*Q2e - 8*l1k*Q2e*Q2h + 32*m2*pow(l1k,2) -
		           8*Q2e*pow(l1k,2) - 64*Q2h*pow(l1k,2) + 48*pow(l1k,3) -
		           28*l1k*pow(Q2e,2) + 16*m2*pow(Q2e,2) + 6*Q2h*pow(Q2e,2) +
		           14*pow(Q2e,3) + 28*l1k*pow(Q2h,2) + 6*Q2e*pow(Q2h,2) -
		           4*pow(Q2h,3))*pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2h,2),-1)) + SI.B0_M0m(-2*l1k + m2,m2)*
		      (2 + (8*m4 - 4*m2*Q2e - 4*m2*Q2h)*pow(l1k,-2) +
		        (-16*m2 + 6*Q2e + 10*Q2h)*pow(l1k,-1) +
		        (-32*l1k + 16*Q2e - 24*Q2h)*pow(2*l1k - m2,-1) +
		        (8*m2 + 4*Q2e + 16*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (40*Q2e*Q2h - 16*pow(Q2e,2) - 24*pow(Q2h,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*m4 + 8*m2*Q2h + 8*pow(Q2h,2)) +
		        (32*m4*Q2e - 8*m2*pow(Q2e,2) + 16*Q2h*pow(Q2e,2) -
		           20*pow(Q2e,3))*pow(2*l1k + Q2e - Q2h,-1)*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        (-16*l1k*m2 - 8*m2*Q2e + 12*Q2e*Q2h - 8*pow(l1k,2) +
		           20*pow(Q2e,2) + 6*pow(Q2h,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*(2*Q2e*pow(Q2h,2) - 2*pow(Q2h,3))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)) +
		     SI.B0_m0m(m2)*(22 + (-8*m4 + 4*m2*Q2e + 4*m2*Q2h)*pow(l1k,-2) +
		        (20*m2 - 2*Q2e - 18*Q2h)*pow(l1k,-1) -
		        32*Q2e*pow(4*m2 + Q2e,-1) +
		        pow(l1k,-1)*(4*Q2e*Q2h - 4*pow(Q2e,2))*pow(4*m2 + Q2e,-1) +
		        (-32*m4 + 16*m2*Q2e + 16*m2*Q2h)*pow(2*l1k + Q2e - Q2h,-2) +
		        (-40*m2 + 4*Q2e + 36*Q2h)*pow(2*l1k + Q2e - Q2h,-1) +
		        (-8*Q2e*Q2h + 8*pow(Q2e,2))*pow(4*m2 + Q2e,-1)*
		         pow(2*l1k + Q2e - Q2h,-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (32*m4 - 16*m2*Q2h - 16*pow(Q2h,2)) +
		        pow(l1k,-1)*(16*m4*Q2e - 12*m2*pow(Q2e,2) - 4*pow(Q2e,3))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (-16*l1k*m2 + 24*l1k*Q2e - 16*m2*Q2e + 8*l1k*Q2h +
		           16*pow(l1k,2) + 4*pow(Q2h,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*(-4*Q2e*pow(Q2h,2) + 4*pow(Q2h,3))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (-32*m4*Q2e + 24*m2*pow(Q2e,2) + 8*pow(Q2e,3))*
		         pow(2*l1k + Q2e - Q2h,-1)*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        (16*l1k*m2 - 8*l1k*Q2e - 8*m2*Q2e - 32*l1k*Q2h +
		           24*pow(l1k,2) - 8*pow(Q2e,2) + 14*pow(Q2h,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*(2*Q2e*pow(Q2h,2) - 2*pow(Q2h,3))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)) +
		     SI.B0_qmm(Q2h,m2)*(8 + (-16*m2 + 10*Q2e)*pow(l1k,-1) +
		        (32*m2 - 20*Q2e)*pow(2*l1k + Q2e - Q2h,-1) -
		        16*m2*Q2h*pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1) +
		        32*pow(Q2e - Q2h,-2)*pow(Q2h,2) +
		        pow(l1k,-1)*(8*m2*pow(Q2e,2) + 8*Q2h*pow(Q2e,2) - 6*pow(Q2e,3))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (-24*l1k*Q2e + 16*m2*Q2e - 16*l1k*Q2h - 16*Q2e*Q2h -
		           16*pow(l1k,2) - 20*pow(Q2e,2) - 8*pow(Q2h,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*(8*Q2e*pow(Q2h,2) - 8*pow(Q2h,3))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (-16*m2*pow(Q2e,2) - 16*Q2h*pow(Q2e,2) + 12*pow(Q2e,3))*
		         pow(2*l1k + Q2e - Q2h,-1)*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        (8*l1k*Q2e + 16*m2*Q2e + 32*l1k*Q2h - 12*Q2e*Q2h -
		           16*pow(l1k,2) - 12*pow(Q2e,2) - 20*pow(Q2h,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*(-4*Q2e*pow(Q2h,2) + 4*pow(Q2h,3))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1))) +
		  G23*(-4*Q2h - 4*S - 4*U +
		     pow(2*l1k + Q2e - Q2h,-1)*
		      (-2*Q2e*Q2h - 4*m2*S + 2*Q2e*S + 24*Q2h*S - 4*m2*U + 6*Q2h*U -
		        4*S*U - 14*pow(Q2h,2) - 12*pow(S,2)) +
		     pow(2*l1k + Q2e - Q2h,-2)*
		      (-4*m2*Q2e*Q2h + 4*m2*Q2e*S + 16*m2*Q2h*S + 12*m2*Q2h*U -
		        8*m2*S*U - 8*m2*pow(Q2h,2) - 8*m2*pow(S,2)) +
		     pow(l1k,-1)*(Q2e*Q2h + 2*m2*S - 3*Q2h*S + 2*m2*U - Q2e*U - 12*Q2h*U +
		        2*S*U + 7*pow(Q2h,2) + 6*pow(U,2)) +
		     pow(Q2e - Q2h,-1)*(-16*Q2h*S - 16*Q2h*U + 32*S*U + 16*pow(S,2) +
		        16*pow(U,2)) + pow(l1k,-2)*
		      (-(m2*Q2e*Q2h) + 3*m2*Q2h*S + m2*Q2e*U + 4*m2*Q2h*U - 2*m2*S*U -
		        2*m2*pow(Q2h,2) - 2*m2*pow(U,2)) +
		     pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		      (-8*m2*S*U + 12*Q2h*S*U + 4*m2*pow(Q2h,2) - 24*S*pow(Q2h,2) -
		        16*U*pow(Q2h,2) + 14*pow(Q2h,3) - 4*m2*pow(S,2) + 10*Q2h*pow(S,2) -
		        4*m2*pow(U,2) + 2*Q2h*pow(U,2)) +
		     pow(l1k,-1)*pow(Q2e - Q2h,-1)*
		      (8*S*pow(Q2h,2) - 8*U*pow(Q2h,2) - 8*Q2h*pow(S,2) + 8*Q2h*pow(U,2)) +
		     pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-2)*
		      (8*m2*Q2h*S*U - 8*m2*S*pow(Q2h,2) - 8*m2*U*pow(Q2h,2) +
		        4*m2*pow(Q2h,3) + 4*m2*Q2h*pow(S,2) + 4*m2*Q2h*pow(U,2)) +
		     pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		      (-4*m2*Q2h*S*U + 4*m2*S*pow(Q2h,2) + 4*m2*U*pow(Q2h,2) -
		        8*S*U*pow(Q2h,2) - 2*m2*pow(Q2h,3) + 8*S*pow(Q2h,3) +
		        8*U*pow(Q2h,3) - 4*pow(Q2h,4) - 2*m2*Q2h*pow(S,2) -
		        4*pow(Q2h,2)*pow(S,2) - 2*m2*Q2h*pow(U,2) - 4*pow(Q2h,2)*pow(U,2)) +
		     pow(l1k,-2)*pow(Q2e - Q2h,-1)*
		      (8*S*U*pow(Q2h,2) - 8*S*pow(Q2h,3) - 8*U*pow(Q2h,3) + 4*pow(Q2h,4) +
		        4*pow(Q2h,2)*pow(S,2) + 4*pow(Q2h,2)*pow(U,2)) +
		     (8*l1k*Q2e*Q2h - 16*l1k*Q2h*S - 8*Q2e*Q2h*S + 8*Q2h*pow(l1k,2) +
		        2*Q2h*pow(Q2e,2) + 4*l1k*pow(Q2h,2) + 2*Q2e*pow(Q2h,2) -
		        8*S*pow(Q2h,2) + 2*pow(Q2h,3) + 8*Q2h*pow(S,2))*
		      pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		     pow(2*l1k + Q2e - Q2h,-1)*
		      (-8*S*pow(Q2h,3) + 2*pow(Q2h,4) + 8*pow(Q2h,2)*pow(S,2))*
		      pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		     (16*l1k*Q2h*U + 8*Q2h*pow(l1k,2) - 12*l1k*pow(Q2h,2) -
		        16*U*pow(Q2h,2) + 6*pow(Q2h,3) + 8*Q2h*pow(U,2))*
		      pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		     pow(l1k,-1)*(4*U*pow(Q2h,3) - pow(Q2h,4) - 4*pow(Q2h,2)*pow(U,2))*
		      pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		     SI.B0_00m(m2)*(-12*Q2h + 12*S + 12*U +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (2*Q2e*Q2h + 4*m2*S - 6*Q2e*S + 16*Q2h*S + 4*m2*U - 2*Q2h*U +
		           12*S*U - 6*pow(Q2h,2) - 12*pow(S,2)) +
		        pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (16*l1k*Q2h + 10*Q2e*Q2h - 12*l1k*S - 6*Q2e*S - 18*Q2h*S -
		           12*l1k*U - 12*Q2e*U + 8*Q2h*U + 12*S*U + 12*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (4*m2*Q2e*Q2h - 4*m2*Q2e*S - 16*m2*Q2h*S - 12*m2*Q2h*U +
		           8*m2*S*U + 8*m2*pow(Q2h,2) + 8*m2*pow(S,2)) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-(Q2e*Q2h*S) + 4*Q2e*Q2h*U + 6*Q2e*S*U - 8*Q2h*S*U +
		           Q2h*pow(Q2e,2) - 3*U*pow(Q2e,2) - 2*Q2e*pow(Q2h,2) +
		           4*S*pow(Q2h,2) - 2*Q2h*pow(S,2)) +
		        pow(2*l1k - m2,-1)*(16*l1k*Q2h - 2*Q2e*Q2h - 12*l1k*S + 6*Q2e*S -
		           2*Q2h*S - 12*l1k*U + 24*Q2h*U - 12*S*U - 8*pow(Q2h,2) -
		           12*pow(U,2)) + pow(l1k,-1)*
		         (-(Q2e*Q2h) - 2*m2*S + Q2h*S - 2*m2*U + 3*Q2e*U - 8*Q2h*U -
		           6*S*U + 3*pow(Q2h,2) + 6*pow(U,2)) +
		        pow(l1k,-2)*(m2*Q2e*Q2h - 3*m2*Q2h*S - m2*Q2e*U - 4*m2*Q2h*U +
		           2*m2*S*U + 2*m2*pow(Q2h,2) + 2*m2*pow(U,2)) +
		        pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (8*Q2e*Q2h*S - 2*Q2e*Q2h*U + 12*Q2e*S*U - 16*Q2h*S*U +
		           2*Q2h*pow(Q2e,2) - 6*S*pow(Q2e,2) - 4*Q2e*pow(Q2h,2) +
		           8*U*pow(Q2h,2) - 4*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (8*m2*S*U + 4*Q2h*S*U - 4*m2*pow(Q2h,2) - 4*S*pow(Q2h,2) -
		           4*U*pow(Q2h,2) + 2*pow(Q2h,3) + 4*m2*pow(S,2) +
		           2*Q2h*pow(S,2) + 4*m2*pow(U,2) + 2*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-2)*
		         (-8*m2*Q2h*S*U + 8*m2*S*pow(Q2h,2) + 8*m2*U*pow(Q2h,2) -
		           4*m2*pow(Q2h,3) - 4*m2*Q2h*pow(S,2) - 4*m2*Q2h*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (4*m2*Q2h*S*U - 4*m2*S*pow(Q2h,2) - 4*m2*U*pow(Q2h,2) +
		           2*m2*pow(Q2h,3) + 2*m2*Q2h*pow(S,2) + 2*m2*Q2h*pow(U,2)) +
		        (-8*l1k*Q2e*Q2h + 16*l1k*Q2h*S + 8*Q2e*Q2h*S - 8*Q2h*pow(l1k,2) -
		           2*Q2h*pow(Q2e,2) - 12*l1k*pow(Q2h,2) - 6*Q2e*pow(Q2h,2) +
		           24*S*pow(Q2h,2) - 14*pow(Q2h,3) - 8*Q2h*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-32*l1k*S*pow(Q2h,2) - 16*Q2e*S*pow(Q2h,2) + 16*l1k*pow(Q2h,3) +
		           8*Q2e*pow(Q2h,3) + 32*S*pow(Q2h,3) - 16*pow(Q2h,4) +
		           16*l1k*Q2h*pow(S,2) + 8*Q2e*Q2h*pow(S,2) -
		           16*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*S*pow(Q2h,3) + 2*pow(Q2h,4) + 8*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (-16*l1k*Q2h*U - 8*Q2h*pow(l1k,2) + 20*l1k*pow(Q2h,2) +
		           32*U*pow(Q2h,2) - 22*pow(Q2h,3) - 8*Q2h*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*(4*U*pow(Q2h,3) - pow(Q2h,4) - 4*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k - m2,-1)*(-32*l1k*U*pow(Q2h,2) + 16*l1k*pow(Q2h,3) -
		           16*U*pow(Q2h,3) + 8*pow(Q2h,4) + 16*l1k*Q2h*pow(U,2) +
		           8*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)) +
		     SI.C0_m0M0mm(2*l1k + m2 + Q2e - Q2h,m2)*
		      (-304*l1k*M2 - 80*m2*M2 - 72*M2*Q2e - 64*M2*Q2h - 4*m2*S -
		        132*Q2h*S - 4*m2*U - 384*Q2h*U + 64*S*U + 244*pow(Q2h,2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-12*m2*Q2h*S - 8*m4*U + 4*m2*pow(Q2h,2) + 8*m2*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (-8*m4*Q2e*Q2h + 8*m4*Q2e*S + 32*m4*Q2h*S + 24*m4*Q2h*U -
		           16*m4*S*U - 16*m4*pow(Q2h,2) - 16*m4*pow(S,2)) +
		        144*pow(U,2) + pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-4*m4*Q2h*S + 4*m4*Q2h*U - 8*m4*S*U - 4*m2*Q2h*S*U +
		           4*m4*pow(Q2h,2) + 6*m2*S*pow(Q2h,2) + 2*m2*U*pow(Q2h,2) -
		           2*m2*pow(Q2h,3) - 4*m2*Q2h*pow(S,2) - 8*m4*pow(U,2)) +
		        pow(l1k,-1)*(-24*m2*M2*Q2e + 2*m2*Q2e*Q2h - 12*M2*Q2e*Q2h -
		           30*m2*Q2h*S - 12*Q2e*Q2h*S + 4*m4*U - 2*m2*Q2e*U -
		           62*m2*Q2h*U - 42*Q2e*Q2h*U + 4*m2*S*U + 52*Q2h*S*U -
		           12*M2*pow(Q2e,2) + 38*m2*pow(Q2h,2) + 26*Q2e*pow(Q2h,2) -
		           74*S*pow(Q2h,2) - 108*U*pow(Q2h,2) + 64*pow(Q2h,3) +
		           12*Q2h*pow(S,2) + 16*m2*pow(U,2) + 16*Q2e*pow(U,2) +
		           44*Q2h*pow(U,2)) + pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-2)*
		         (16*m4*Q2h*S*U - 16*m4*S*pow(Q2h,2) - 16*m4*U*pow(Q2h,2) +
		           8*m4*pow(Q2h,3) + 8*m4*Q2h*pow(S,2) + 8*m4*Q2h*pow(U,2)) +
		        pow(l1k,-1)*(4*S*U*pow(Q2e,5)*pow(Q2h,2) -
		           4*U*pow(Q2e,6)*pow(Q2h,2) - 8*S*U*pow(Q2e,4)*pow(Q2h,3) -
		           4*S*pow(Q2e,5)*pow(Q2h,3) + 8*U*pow(Q2e,5)*pow(Q2h,3) +
		           2*pow(Q2e,6)*pow(Q2h,3) + 4*S*U*pow(Q2e,3)*pow(Q2h,4) +
		           8*S*pow(Q2e,4)*pow(Q2h,4) - 4*U*pow(Q2e,4)*pow(Q2h,4) -
		           4*pow(Q2e,5)*pow(Q2h,4) - 4*S*pow(Q2e,3)*pow(Q2h,5) +
		           2*pow(Q2e,4)*pow(Q2h,5) + 2*pow(Q2e,4)*pow(Q2h,3)*pow(S,2) -
		           4*pow(Q2e,3)*pow(Q2h,4)*pow(S,2) +
		           2*pow(Q2e,2)*pow(Q2h,5)*pow(S,2) + 2*Q2h*pow(Q2e,6)*pow(U,2) -
		           4*pow(Q2e,5)*pow(Q2h,2)*pow(U,2) +
		           2*pow(Q2e,4)*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        (1536*m4*Q2e*Q2h*S*U*pow(l1k,2) - 1024*m4*Q2e*Q2h*S*pow(l1k,3) -
		           1024*m4*Q2e*Q2h*U*pow(l1k,3) + 512*m4*Q2e*S*U*pow(l1k,3) -
		           512*m4*Q2h*S*U*pow(l1k,3) + 7808*m2*Q2e*Q2h*S*U*pow(l1k,3) -
		           4096*m2*Q2e*Q2h*S*pow(l1k,4) - 6144*m2*Q2e*Q2h*U*pow(l1k,4) +
		           2560*m2*Q2e*S*U*pow(l1k,4) + 12032*Q2e*Q2h*S*U*pow(l1k,4) -
		           14336*Q2e*Q2h*S*pow(l1k,5) - 26624*Q2e*Q2h*U*pow(l1k,5) +
		           9728*Q2e*S*U*pow(l1k,5) + 512*l1k*Q2e*Q2h*S*U*pow(m,6) +
		           10368*Q2h*S*U*pow(l1k,3)*pow(Q2e,2) -
		           11264*Q2h*S*pow(l1k,4)*pow(Q2e,2) -
		           24576*Q2h*U*pow(l1k,4)*pow(Q2e,2) +
		           8064*S*U*pow(l1k,4)*pow(Q2e,2) +
		           1984*Q2h*S*U*pow(l1k,2)*pow(Q2e,3) -
		           2624*Q2h*S*pow(l1k,3)*pow(Q2e,3) -
		           8256*Q2h*U*pow(l1k,3)*pow(Q2e,3) +
		           2112*S*U*pow(l1k,3)*pow(Q2e,3) + 256*l1k*Q2h*S*U*pow(Q2e,4) -
		           320*Q2h*S*pow(l1k,2)*pow(Q2e,4) -
		           1632*Q2h*U*pow(l1k,2)*pow(Q2e,4) +
		           288*S*U*pow(l1k,2)*pow(Q2e,4) - 16*l1k*Q2h*S*pow(Q2e,5) -
		           176*l1k*Q2h*U*pow(Q2e,5) + 16*l1k*S*U*pow(Q2e,5) +
		           16*Q2h*S*U*pow(Q2e,5) - 8*Q2h*U*pow(Q2e,6) +
		           1536*l1k*m4*Q2e*S*U*pow(Q2h,2) -
		           2560*m4*Q2e*S*pow(l1k,2)*pow(Q2h,2) -
		           3072*m4*Q2e*U*pow(l1k,2)*pow(Q2h,2) +
		           6528*m2*Q2e*S*U*pow(l1k,2)*pow(Q2h,2) +
		           1024*m4*Q2e*pow(l1k,3)*pow(Q2h,2) -
		           10752*m2*Q2e*S*pow(l1k,3)*pow(Q2h,2) -
		           15872*m2*Q2e*U*pow(l1k,3)*pow(Q2h,2) -
		           2432*m2*S*U*pow(l1k,3)*pow(Q2h,2) -
		           1664*Q2e*S*U*pow(l1k,3)*pow(Q2h,2) +
		           5120*m2*Q2e*pow(l1k,4)*pow(Q2h,2) -
		           14848*Q2e*S*pow(l1k,4)*pow(Q2h,2) -
		           16896*Q2e*U*pow(l1k,4)*pow(Q2h,2) +
		           20480*Q2e*pow(l1k,5)*pow(Q2h,2) -
		           512*l1k*Q2e*S*pow(m,6)*pow(Q2h,2) -
		           512*l1k*Q2e*U*pow(m,6)*pow(Q2h,2) -
		           512*l1k*S*U*pow(m,6)*pow(Q2h,2) +
		           2496*S*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2) -
		           13312*S*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) -
		           15104*U*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) +
		           17920*pow(l1k,4)*pow(Q2e,2)*pow(Q2h,2) +
		           288*l1k*S*U*pow(Q2e,3)*pow(Q2h,2) -
		           2432*S*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2) -
		           2336*U*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2) +
		           5440*pow(l1k,3)*pow(Q2e,3)*pow(Q2h,2) -
		           288*l1k*S*pow(Q2e,4)*pow(Q2h,2) -
		           368*l1k*U*pow(Q2e,4)*pow(Q2h,2) +
		           40*S*U*pow(Q2e,4)*pow(Q2h,2) +
		           976*pow(l1k,2)*pow(Q2e,4)*pow(Q2h,2) +
		           96*l1k*pow(Q2e,5)*pow(Q2h,2) - 16*S*pow(Q2e,5)*pow(Q2h,2) -
		           56*U*pow(Q2e,5)*pow(Q2h,2) + 4*pow(Q2e,6)*pow(Q2h,2) -
		           1664*l1k*m4*Q2e*S*pow(Q2h,3) - 1664*l1k*m4*Q2e*U*pow(Q2h,3) -
		           1280*l1k*m4*S*U*pow(Q2h,3) + 2080*l1k*m2*Q2e*S*U*pow(Q2h,3) +
		           256*m4*Q2e*S*U*pow(Q2h,3) +
		           2048*m4*Q2e*pow(l1k,2)*pow(Q2h,3) +
		           1024*m4*S*pow(l1k,2)*pow(Q2h,3) -
		           8192*m2*Q2e*S*pow(l1k,2)*pow(Q2h,3) +
		           1536*m4*U*pow(l1k,2)*pow(Q2h,3) -
		           10496*m2*Q2e*U*pow(l1k,2)*pow(Q2h,3) -
		           4224*m2*S*U*pow(l1k,2)*pow(Q2h,3) -
		           2368*Q2e*S*U*pow(l1k,2)*pow(Q2h,3) +
		           10752*m2*Q2e*pow(l1k,3)*pow(Q2h,3) +
		           3584*m2*S*pow(l1k,3)*pow(Q2h,3) +
		           3072*Q2e*S*pow(l1k,3)*pow(Q2h,3) +
		           6656*m2*U*pow(l1k,3)*pow(Q2h,3) +
		           5888*Q2e*U*pow(l1k,3)*pow(Q2h,3) +
		           8192*Q2e*pow(l1k,4)*pow(Q2h,3) +
		           256*l1k*Q2e*pow(m,6)*pow(Q2h,3) +
		           512*l1k*S*pow(m,6)*pow(Q2h,3) +
		           512*l1k*U*pow(m,6)*pow(Q2h,3) +
		           96*l1k*S*U*pow(Q2e,2)*pow(Q2h,3) -
		           2816*S*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3) -
		           2816*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3) +
		           8768*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,3) -
		           352*l1k*S*pow(Q2e,3)*pow(Q2h,3) -
		           112*l1k*U*pow(Q2e,3)*pow(Q2h,3) -
		           48*S*U*pow(Q2e,3)*pow(Q2h,3) +
		           1408*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,3) +
		           232*l1k*pow(Q2e,4)*pow(Q2h,3) - 48*S*pow(Q2e,4)*pow(Q2h,3) +
		           56*U*pow(Q2e,4)*pow(Q2h,3) + 32*pow(Q2e,5)*pow(Q2h,3) +
		           896*l1k*m4*Q2e*pow(Q2h,4) + 1408*l1k*m4*S*pow(Q2h,4) -
		           2368*l1k*m2*Q2e*S*pow(Q2h,4) - 256*m4*Q2e*S*pow(Q2h,4) +
		           1408*l1k*m4*U*pow(Q2h,4) - 2688*l1k*m2*Q2e*U*pow(Q2h,4) -
		           256*m4*Q2e*U*pow(Q2h,4) - 1696*l1k*m2*S*U*pow(Q2h,4) -
		           256*m4*S*U*pow(Q2h,4) - 512*l1k*Q2e*S*U*pow(Q2h,4) +
		           256*m2*Q2e*S*U*pow(Q2h,4) - 1280*m4*pow(l1k,2)*pow(Q2h,4) +
		           6272*m2*Q2e*pow(l1k,2)*pow(Q2h,4) +
		           5504*m2*S*pow(l1k,2)*pow(Q2h,4) +
		           2944*Q2e*S*pow(l1k,2)*pow(Q2h,4) +
		           7552*m2*U*pow(l1k,2)*pow(Q2h,4) +
		           3584*Q2e*U*pow(l1k,2)*pow(Q2h,4) -
		           5120*m2*pow(l1k,3)*pow(Q2h,4) -
		           4032*Q2e*pow(l1k,3)*pow(Q2h,4) - 256*l1k*pow(m,6)*pow(Q2h,4) -
		           8*S*U*pow(Q2e,2)*pow(Q2h,4) +
		           1408*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,4) +
		           64*S*pow(Q2e,3)*pow(Q2h,4) + 8*U*pow(Q2e,3)*pow(Q2h,4) -
		           36*pow(Q2e,4)*pow(Q2h,4) - 768*l1k*m4*pow(Q2h,5) +
		           1488*l1k*m2*Q2e*pow(Q2h,5) + 128*m4*Q2e*pow(Q2h,5) +
		           1984*l1k*m2*S*pow(Q2h,5) + 256*m4*S*pow(Q2h,5) +
		           512*l1k*Q2e*S*pow(Q2h,5) - 256*m2*Q2e*S*pow(Q2h,5) +
		           2304*l1k*m2*U*pow(Q2h,5) + 256*m4*U*pow(Q2h,5) +
		           512*l1k*Q2e*U*pow(Q2h,5) - 256*m2*Q2e*U*pow(Q2h,5) -
		           256*m2*S*U*pow(Q2h,5) - 4608*m2*pow(l1k,2)*pow(Q2h,5) -
		           2080*Q2e*pow(l1k,2)*pow(Q2h,5) - 1296*l1k*m2*pow(Q2h,6) -
		           128*m4*pow(Q2h,6) - 256*l1k*Q2e*pow(Q2h,6) +
		           128*m2*Q2e*pow(Q2h,6) + 256*m2*S*pow(Q2h,6) +
		           256*m2*U*pow(Q2h,6) - 128*m2*pow(Q2h,7) +
		           512*m4*Q2e*Q2h*pow(l1k,2)*pow(S,2) +
		           256*m4*Q2e*pow(l1k,3)*pow(S,2) -
		           256*m4*Q2h*pow(l1k,3)*pow(S,2) +
		           2240*m2*Q2e*Q2h*pow(l1k,3)*pow(S,2) +
		           768*m2*Q2e*pow(l1k,4)*pow(S,2) +
		           3456*Q2e*Q2h*pow(l1k,4)*pow(S,2) +
		           2304*Q2e*pow(l1k,5)*pow(S,2) +
		           256*l1k*Q2e*Q2h*pow(m,6)*pow(S,2) +
		           2784*Q2h*pow(l1k,3)*pow(Q2e,2)*pow(S,2) +
		           1600*pow(l1k,4)*pow(Q2e,2)*pow(S,2) +
		           384*Q2h*pow(l1k,2)*pow(Q2e,3)*pow(S,2) +
		           256*pow(l1k,3)*pow(Q2e,3)*pow(S,2) +
		           24*l1k*Q2h*pow(Q2e,4)*pow(S,2) +
		           16*pow(l1k,2)*pow(Q2e,4)*pow(S,2) +
		           768*l1k*m4*Q2e*pow(Q2h,2)*pow(S,2) +
		           256*m4*pow(l1k,2)*pow(Q2h,2)*pow(S,2) +
		           2176*m2*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(S,2) -
		           576*m2*pow(l1k,3)*pow(Q2h,2)*pow(S,2) +
		           448*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(S,2) -
		           256*l1k*pow(m,6)*pow(Q2h,2)*pow(S,2) +
		           1376*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) +
		           176*l1k*pow(Q2e,3)*pow(Q2h,2)*pow(S,2) +
		           12*pow(Q2e,4)*pow(Q2h,2)*pow(S,2) -
		           640*l1k*m4*pow(Q2h,3)*pow(S,2) +
		           880*l1k*m2*Q2e*pow(Q2h,3)*pow(S,2) +
		           128*m4*Q2e*pow(Q2h,3)*pow(S,2) -
		           1152*m2*pow(l1k,2)*pow(Q2h,3)*pow(S,2) -
		           864*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(S,2) +
		           128*l1k*pow(Q2e,2)*pow(Q2h,3)*pow(S,2) +
		           16*pow(Q2e,3)*pow(Q2h,3)*pow(S,2) -
		           688*l1k*m2*pow(Q2h,4)*pow(S,2) - 128*m4*pow(Q2h,4)*pow(S,2) -
		           256*l1k*Q2e*pow(Q2h,4)*pow(S,2) +
		           128*m2*Q2e*pow(Q2h,4)*pow(S,2) -
		           28*pow(Q2e,2)*pow(Q2h,4)*pow(S,2) -
		           128*m2*pow(Q2h,5)*pow(S,2) +
		           1024*m4*Q2e*Q2h*pow(l1k,2)*pow(U,2) +
		           256*m4*Q2e*pow(l1k,3)*pow(U,2) -
		           256*m4*Q2h*pow(l1k,3)*pow(U,2) +
		           5824*m2*Q2e*Q2h*pow(l1k,3)*pow(U,2) +
		           1792*m2*Q2e*pow(l1k,4)*pow(U,2) +
		           8064*Q2e*Q2h*pow(l1k,4)*pow(U,2) +
		           8448*Q2e*pow(l1k,5)*pow(U,2) +
		           256*l1k*Q2e*Q2h*pow(m,6)*pow(U,2) +
		           6496*Q2h*pow(l1k,3)*pow(Q2e,2)*pow(U,2) +
		           8256*pow(l1k,4)*pow(Q2e,2)*pow(U,2) +
		           992*Q2h*pow(l1k,2)*pow(Q2e,3)*pow(U,2) +
		           3072*pow(l1k,3)*pow(Q2e,3)*pow(U,2) +
		           144*l1k*Q2h*pow(Q2e,4)*pow(U,2) +
		           672*pow(l1k,2)*pow(Q2e,4)*pow(U,2) +
		           80*l1k*pow(Q2e,5)*pow(U,2) + 24*Q2h*pow(Q2e,5)*pow(U,2) +
		           4*pow(Q2e,6)*pow(U,2) + 768*l1k*m4*Q2e*pow(Q2h,2)*pow(U,2) -
		           256*m4*pow(l1k,2)*pow(Q2h,2)*pow(U,2) +
		           4352*m2*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(U,2) -
		           2112*m2*pow(l1k,3)*pow(Q2h,2)*pow(U,2) -
		           2112*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(U,2) -
		           256*l1k*pow(m,6)*pow(Q2h,2)*pow(U,2) +
		           1328*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) +
		           96*l1k*pow(Q2e,3)*pow(Q2h,2)*pow(U,2) -
		           20*pow(Q2e,4)*pow(Q2h,2)*pow(U,2) -
		           640*l1k*m4*pow(Q2h,3)*pow(U,2) +
		           1200*l1k*m2*Q2e*pow(Q2h,3)*pow(U,2) +
		           128*m4*Q2e*pow(Q2h,3)*pow(U,2) -
		           3072*m2*pow(l1k,2)*pow(Q2h,3)*pow(U,2) -
		           1504*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(U,2) +
		           8*l1k*pow(Q2e,2)*pow(Q2h,3)*pow(U,2) -
		           8*pow(Q2e,3)*pow(Q2h,3)*pow(U,2) -
		           1008*l1k*m2*pow(Q2h,4)*pow(U,2) -
		           128*m4*pow(Q2h,4)*pow(U,2) - 256*l1k*Q2e*pow(Q2h,4)*pow(U,2) +
		           128*m2*Q2e*pow(Q2h,4)*pow(U,2) - 128*m2*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        (128*l1k*m4*M2*Q2e - 128*l1k*m4*M2*Q2h +
		           352*l1k*m2*M2*Q2e*Q2h + 320*l1k*m4*Q2h*S +
		           552*l1k*m2*Q2e*Q2h*S + 128*m4*Q2e*Q2h*S + 320*l1k*m4*Q2h*U +
		           1144*l1k*m2*Q2e*Q2h*U + 224*m4*Q2e*Q2h*U - 128*l1k*m4*S*U -
		           320*l1k*m2*Q2e*S*U - 64*m4*Q2e*S*U - 1024*l1k*m2*Q2h*S*U -
		           256*m4*Q2h*S*U - 944*l1k*Q2e*Q2h*S*U - 360*m2*Q2e*Q2h*S*U +
		           384*m2*M2*Q2e*pow(l1k,2) - 320*M2*Q2e*Q2h*pow(l1k,2) +
		           1056*m2*Q2h*S*pow(l1k,2) + 1520*Q2e*Q2h*S*pow(l1k,2) +
		           1568*m2*Q2h*U*pow(l1k,2) + 4352*Q2e*Q2h*U*pow(l1k,2) -
		           640*m2*S*U*pow(l1k,2) - 1056*Q2e*S*U*pow(l1k,2) -
		           4224*Q2h*S*U*pow(l1k,2) + 1152*M2*Q2e*pow(l1k,3) +
		           3584*Q2h*S*pow(l1k,3) + 6656*Q2h*U*pow(l1k,3) -
		           2432*S*U*pow(l1k,3) + 80*l1k*M2*Q2h*pow(Q2e,2) +
		           416*l1k*Q2h*S*pow(Q2e,2) + 1584*l1k*Q2h*U*pow(Q2e,2) -
		           256*l1k*S*U*pow(Q2e,2) - 216*Q2h*S*U*pow(Q2e,2) +
		           800*M2*pow(l1k,2)*pow(Q2e,2) + 128*l1k*M2*pow(Q2e,3) +
		           32*M2*Q2h*pow(Q2e,3) + 20*Q2h*S*pow(Q2e,3) +
		           164*Q2h*U*pow(Q2e,3) - 8*S*U*pow(Q2e,3) + 8*M2*pow(Q2e,4) -
		           256*l1k*m4*pow(Q2h,2) - 288*l1k*m2*M2*pow(Q2h,2) -
		           800*l1k*m2*Q2e*pow(Q2h,2) - 144*m4*Q2e*pow(Q2h,2) -
		           128*l1k*M2*Q2e*pow(Q2h,2) + 64*m2*M2*Q2e*pow(Q2h,2) +
		           1368*l1k*m2*S*pow(Q2h,2) + 160*m4*S*pow(Q2h,2) +
		           1136*l1k*Q2e*S*pow(Q2h,2) + 532*m2*Q2e*S*pow(Q2h,2) +
		           1384*l1k*m2*U*pow(Q2h,2) + 128*m4*U*pow(Q2h,2) +
		           1072*l1k*Q2e*U*pow(Q2h,2) + 884*m2*Q2e*U*pow(Q2h,2) -
		           1696*l1k*S*U*pow(Q2h,2) - 232*m2*S*U*pow(Q2h,2) -
		           128*Q2e*S*U*pow(Q2h,2) - 1280*m2*pow(l1k,2)*pow(Q2h,2) -
		           2896*Q2e*pow(l1k,2)*pow(Q2h,2) + 5504*S*pow(l1k,2)*pow(Q2h,2) +
		           7552*U*pow(l1k,2)*pow(Q2h,2) - 5120*pow(l1k,3)*pow(Q2h,2) -
		           992*l1k*pow(Q2e,2)*pow(Q2h,2) - 40*M2*pow(Q2e,2)*pow(Q2h,2) +
		           284*S*pow(Q2e,2)*pow(Q2h,2) + 344*U*pow(Q2e,2)*pow(Q2h,2) -
		           100*pow(Q2e,3)*pow(Q2h,2) - 768*l1k*m2*pow(Q2h,3) -
		           16*m4*pow(Q2h,3) - 64*m2*M2*pow(Q2h,3) -
		           616*l1k*Q2e*pow(Q2h,3) - 540*m2*Q2e*pow(Q2h,3) +
		           1984*l1k*S*pow(Q2h,3) + 148*m2*S*pow(Q2h,3) +
		           108*Q2e*S*pow(Q2h,3) + 2304*l1k*U*pow(Q2h,3) -
		           124*m2*U*pow(Q2h,3) + 40*Q2e*U*pow(Q2h,3) -
		           256*S*U*pow(Q2h,3) - 4608*pow(l1k,2)*pow(Q2h,3) -
		           204*pow(Q2e,2)*pow(Q2h,3) - 1296*l1k*pow(Q2h,4) +
		           116*m2*pow(Q2h,4) + 256*S*pow(Q2h,4) + 256*U*pow(Q2h,4) -
		           128*pow(Q2h,5) - 64*l1k*m4*pow(S,2) -
		           64*l1k*m2*Q2e*pow(S,2) - 368*l1k*m2*Q2h*pow(S,2) -
		           128*m4*Q2h*pow(S,2) - 200*l1k*Q2e*Q2h*pow(S,2) -
		           64*m2*Q2e*Q2h*pow(S,2) - 192*m2*pow(l1k,2)*pow(S,2) -
		           112*Q2e*pow(l1k,2)*pow(S,2) - 1152*Q2h*pow(l1k,2)*pow(S,2) -
		           576*pow(l1k,3)*pow(S,2) - 8*l1k*pow(Q2e,2)*pow(S,2) -
		           24*Q2h*pow(Q2e,2)*pow(S,2) - 688*l1k*pow(Q2h,2)*pow(S,2) -
		           192*m2*pow(Q2h,2)*pow(S,2) - 84*Q2e*pow(Q2h,2)*pow(S,2) -
		           128*pow(Q2h,3)*pow(S,2) - 64*l1k*m4*pow(U,2) -
		           384*l1k*m2*Q2e*pow(U,2) - 64*m4*Q2e*pow(U,2) -
		           592*l1k*m2*Q2h*pow(U,2) - 128*m4*Q2h*pow(U,2) -
		           472*l1k*Q2e*Q2h*pow(U,2) - 344*m2*Q2e*Q2h*pow(U,2) -
		           448*m2*pow(l1k,2)*pow(U,2) - 1584*Q2e*pow(l1k,2)*pow(U,2) -
		           3072*Q2h*pow(l1k,2)*pow(U,2) - 2112*pow(l1k,3)*pow(U,2) -
		           616*l1k*pow(Q2e,2)*pow(U,2) - 136*Q2h*pow(Q2e,2)*pow(U,2) -
		           68*pow(Q2e,3)*pow(U,2) - 1008*l1k*pow(Q2h,2)*pow(U,2) +
		           8*m2*pow(Q2h,2)*pow(U,2) - 40*Q2e*pow(Q2h,2)*pow(U,2) -
		           128*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*(-96*m4*Q2e*Q2h*S*U - 8*Q2h*S*U*pow(Q2e,3) +
		           4*M2*Q2h*pow(Q2e,4) + 6*Q2h*U*pow(Q2e,4) +
		           96*m4*Q2e*S*pow(Q2h,2) + 128*m4*Q2e*U*pow(Q2h,2) +
		           64*m4*S*U*pow(Q2h,2) - 128*m2*Q2e*S*U*pow(Q2h,2) -
		           24*S*U*pow(Q2e,2)*pow(Q2h,2) - 8*M2*pow(Q2e,3)*pow(Q2h,2) +
		           14*S*pow(Q2e,3)*pow(Q2h,2) + 22*U*pow(Q2e,3)*pow(Q2h,2) -
		           4*pow(Q2e,4)*pow(Q2h,2) - 64*m4*Q2e*pow(Q2h,3) -
		           64*m4*S*pow(Q2h,3) + 148*m2*Q2e*S*pow(Q2h,3) -
		           96*m4*U*pow(Q2h,3) + 180*m2*Q2e*U*pow(Q2h,3) +
		           56*m2*S*U*pow(Q2h,3) + 4*M2*pow(Q2e,2)*pow(Q2h,3) +
		           18*S*pow(Q2e,2)*pow(Q2h,3) + 4*U*pow(Q2e,2)*pow(Q2h,3) -
		           12*pow(Q2e,3)*pow(Q2h,3) + 48*m4*pow(Q2h,4) -
		           100*m2*Q2e*pow(Q2h,4) - 76*m2*S*pow(Q2h,4) -
		           108*m2*U*pow(Q2h,4) + 64*m2*pow(Q2h,5) -
		           32*m4*Q2e*Q2h*pow(S,2) + 16*m4*pow(Q2h,2)*pow(S,2) -
		           48*m2*Q2e*pow(Q2h,2)*pow(S,2) -
		           10*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) + 12*m2*pow(Q2h,3)*pow(S,2) -
		           6*Q2e*pow(Q2h,3)*pow(S,2) - 64*m4*Q2e*Q2h*pow(U,2) -
		           10*Q2h*pow(Q2e,3)*pow(U,2) - 2*pow(Q2e,4)*pow(U,2) +
		           48*m4*pow(Q2h,2)*pow(U,2) - 80*m2*Q2e*pow(Q2h,2)*pow(U,2) -
		           4*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) + 44*m2*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1)) +
		     SI.C0_0qQmmm(Q2e,Q2h,m2)*
		      (256*l1k*M2 + 192*m2*M2 + 80*M2*Q2e + 48*M2*Q2h + 248*Q2h*S +
		        472*Q2h*U - 96*S*U - 336*pow(Q2h,2) - 48*pow(S,2) - 176*pow(U,2) +
		        pow(Q2e - Q2h,-1)*(-32*m2*Q2h*S - 32*m2*Q2h*U + 64*m2*S*U +
		           32*m2*pow(S,2) + 32*m2*pow(U,2)) +
		        pow(l1k,-1)*(48*m2*M2*Q2e + 12*M2*Q2e*Q2h + 48*m2*Q2h*S +
		           16*Q2e*Q2h*S + 72*m2*Q2h*U + 58*Q2e*Q2h*U - 68*Q2h*S*U +
		           16*M2*pow(Q2e,2) - 48*m2*pow(Q2h,2) + 4*M2*pow(Q2h,2) -
		           34*Q2e*pow(Q2h,2) + 90*S*pow(Q2h,2) + 136*U*pow(Q2h,2) -
		           76*pow(Q2h,3) - 16*Q2h*pow(S,2) - 16*m2*pow(U,2) -
		           24*Q2e*pow(U,2) - 60*Q2h*pow(U,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-96*m2*M2*Q2e - 24*M2*Q2e*Q2h - 144*m2*Q2h*S -
		           116*Q2e*Q2h*S - 96*m2*Q2h*U - 32*Q2e*Q2h*U + 136*Q2h*S*U -
		           32*M2*pow(Q2e,2) + 96*m2*pow(Q2h,2) - 8*M2*pow(Q2h,2) +
		           68*Q2e*pow(Q2h,2) - 272*S*pow(Q2h,2) - 180*U*pow(Q2h,2) +
		           152*pow(Q2h,3) + 32*m2*pow(S,2) + 48*Q2e*pow(S,2) +
		           120*Q2h*pow(S,2) + 32*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(Q2e - Q2h,-1)*
		         (16*m2*S*pow(Q2h,2) - 16*m2*U*pow(Q2h,2) - 16*m2*Q2h*pow(S,2) +
		           16*m2*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-112*m2*Q2h*S*U + 136*m2*S*pow(Q2h,2) + 120*m2*U*pow(Q2h,2) -
		           152*S*U*pow(Q2h,2) - 72*m2*pow(Q2h,3) + 212*S*pow(Q2h,3) +
		           172*U*pow(Q2h,3) - 116*pow(Q2h,4) - 64*m2*Q2h*pow(S,2) -
		           96*pow(Q2h,2)*pow(S,2) - 48*m2*Q2h*pow(U,2) -
		           56*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-2)*pow(Q2e - Q2h,-1)*
		         (16*m2*S*U*pow(Q2h,2) - 16*m2*S*pow(Q2h,3) -
		           16*m2*U*pow(Q2h,3) + 8*m2*pow(Q2h,4) +
		           8*m2*pow(Q2h,2)*pow(S,2) + 8*m2*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (16*m2*S*U*pow(Q2h,2) - 16*m2*S*pow(Q2h,3) -
		           16*m2*U*pow(Q2h,3) + 40*S*U*pow(Q2h,3) + 8*m2*pow(Q2h,4) -
		           40*S*pow(Q2h,4) - 40*U*pow(Q2h,4) + 20*pow(Q2h,5) +
		           8*m2*pow(Q2h,2)*pow(S,2) + 20*pow(Q2h,3)*pow(S,2) +
		           8*m2*pow(Q2h,2)*pow(U,2) + 20*pow(Q2h,3)*pow(U,2)) +
		        pow(l1k,-1)*(-4*S*U*pow(Q2e,5)*pow(Q2h,2) +
		           4*U*pow(Q2e,6)*pow(Q2h,2) + 8*S*U*pow(Q2e,4)*pow(Q2h,3) +
		           4*S*pow(Q2e,5)*pow(Q2h,3) - 8*U*pow(Q2e,5)*pow(Q2h,3) -
		           2*pow(Q2e,6)*pow(Q2h,3) - 4*S*U*pow(Q2e,3)*pow(Q2h,4) -
		           8*S*pow(Q2e,4)*pow(Q2h,4) + 4*U*pow(Q2e,4)*pow(Q2h,4) +
		           4*pow(Q2e,5)*pow(Q2h,4) + 4*S*pow(Q2e,3)*pow(Q2h,5) -
		           2*pow(Q2e,4)*pow(Q2h,5) - 2*pow(Q2e,4)*pow(Q2h,3)*pow(S,2) +
		           4*pow(Q2e,3)*pow(Q2h,4)*pow(S,2) -
		           2*pow(Q2e,2)*pow(Q2h,5)*pow(S,2) - 2*Q2h*pow(Q2e,6)*pow(U,2) +
		           4*pow(Q2e,5)*pow(Q2h,2)*pow(U,2) -
		           2*pow(Q2e,4)*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        (-2048*m4*Q2e*Q2h*S*U*pow(l1k,2) -
		           6144*m2*Q2e*Q2h*S*U*pow(l1k,3) +
		           4096*m2*Q2e*Q2h*S*pow(l1k,4) + 4096*m2*Q2e*Q2h*U*pow(l1k,4) -
		           2048*m2*Q2e*S*U*pow(l1k,4) - 11264*Q2e*Q2h*S*U*pow(l1k,4) +
		           12288*Q2e*Q2h*S*pow(l1k,5) + 20480*Q2e*Q2h*U*pow(l1k,5) -
		           8192*Q2e*S*U*pow(l1k,5) - 9216*Q2h*S*U*pow(l1k,3)*pow(Q2e,2) +
		           10240*Q2h*S*pow(l1k,4)*pow(Q2e,2) +
		           20480*Q2h*U*pow(l1k,4)*pow(Q2e,2) -
		           7168*S*U*pow(l1k,4)*pow(Q2e,2) -
		           1792*Q2h*S*U*pow(l1k,2)*pow(Q2e,3) +
		           2560*Q2h*S*pow(l1k,3)*pow(Q2e,3) +
		           7680*Q2h*U*pow(l1k,3)*pow(Q2e,3) -
		           2048*S*U*pow(l1k,3)*pow(Q2e,3) - 240*l1k*Q2h*S*U*pow(Q2e,4) +
		           320*Q2h*S*pow(l1k,2)*pow(Q2e,4) +
		           1600*Q2h*U*pow(l1k,2)*pow(Q2e,4) -
		           288*S*U*pow(l1k,2)*pow(Q2e,4) + 16*l1k*Q2h*S*pow(Q2e,5) +
		           176*l1k*Q2h*U*pow(Q2e,5) - 16*l1k*S*U*pow(Q2e,5) -
		           16*Q2h*S*U*pow(Q2e,5) + 8*Q2h*U*pow(Q2e,6) -
		           1024*l1k*m4*Q2e*S*U*pow(Q2h,2) +
		           3072*m4*Q2e*S*pow(l1k,2)*pow(Q2h,2) +
		           3072*m4*Q2e*U*pow(l1k,2)*pow(Q2h,2) -
		           6400*m2*Q2e*S*U*pow(l1k,2)*pow(Q2h,2) +
		           8192*m2*Q2e*S*pow(l1k,3)*pow(Q2h,2) +
		           12288*m2*Q2e*U*pow(l1k,3)*pow(Q2h,2) +
		           2048*m2*S*U*pow(l1k,3)*pow(Q2h,2) +
		           2560*Q2e*S*U*pow(l1k,3)*pow(Q2h,2) -
		           4096*m2*Q2e*pow(l1k,4)*pow(Q2h,2) +
		           14336*Q2e*S*pow(l1k,4)*pow(Q2h,2) +
		           16384*Q2e*U*pow(l1k,4)*pow(Q2h,2) -
		           16384*Q2e*pow(l1k,5)*pow(Q2h,2) -
		           1888*S*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2) +
		           11776*S*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) +
		           12800*U*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) -
		           15360*pow(l1k,4)*pow(Q2e,2)*pow(Q2h,2) -
		           176*l1k*S*U*pow(Q2e,3)*pow(Q2h,2) +
		           2240*S*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2) +
		           1728*U*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2) -
		           5120*pow(l1k,3)*pow(Q2e,3)*pow(Q2h,2) +
		           256*l1k*S*pow(Q2e,4)*pow(Q2h,2) +
		           256*l1k*U*pow(Q2e,4)*pow(Q2h,2) -
		           40*S*U*pow(Q2e,4)*pow(Q2h,2) -
		           960*pow(l1k,2)*pow(Q2e,4)*pow(Q2h,2) -
		           96*l1k*pow(Q2e,5)*pow(Q2h,2) + 24*S*pow(Q2e,5)*pow(Q2h,2) +
		           48*U*pow(Q2e,5)*pow(Q2h,2) - 4*pow(Q2e,6)*pow(Q2h,2) +
		           1024*l1k*m4*Q2e*S*pow(Q2h,3) + 1024*l1k*m4*Q2e*U*pow(Q2h,3) +
		           1024*l1k*m4*S*U*pow(Q2h,3) - 1280*l1k*m2*Q2e*S*U*pow(Q2h,3) -
		           256*m4*Q2e*S*U*pow(Q2h,3) -
		           2048*m4*Q2e*pow(l1k,2)*pow(Q2h,3) -
		           1024*m4*S*pow(l1k,2)*pow(Q2h,3) +
		           8192*m2*Q2e*S*pow(l1k,2)*pow(Q2h,3) -
		           1024*m4*U*pow(l1k,2)*pow(Q2h,3) +
		           9728*m2*Q2e*U*pow(l1k,2)*pow(Q2h,3) +
		           3840*m2*S*U*pow(l1k,2)*pow(Q2h,3) +
		           1536*Q2e*S*U*pow(l1k,2)*pow(Q2h,3) -
		           8192*m2*Q2e*pow(l1k,3)*pow(Q2h,3) -
		           3072*m2*S*pow(l1k,3)*pow(Q2h,3) -
		           4608*Q2e*S*pow(l1k,3)*pow(Q2h,3) -
		           5120*m2*U*pow(l1k,3)*pow(Q2h,3) -
		           6656*Q2e*U*pow(l1k,3)*pow(Q2h,3) -
		           8192*Q2e*pow(l1k,4)*pow(Q2h,3) -
		           80*l1k*S*U*pow(Q2e,2)*pow(Q2h,3) +
		           2048*S*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3) +
		           2048*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3) -
		           7168*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,3) +
		           240*l1k*S*pow(Q2e,3)*pow(Q2h,3) +
		           80*l1k*U*pow(Q2e,3)*pow(Q2h,3) + 48*S*U*pow(Q2e,3)*pow(Q2h,3) -
		           1024*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,3) -
		           160*l1k*pow(Q2e,4)*pow(Q2h,3) + 32*S*pow(Q2e,4)*pow(Q2h,3) -
		           40*U*pow(Q2e,4)*pow(Q2h,3) - 32*pow(Q2e,5)*pow(Q2h,3) -
		           512*l1k*m4*Q2e*pow(Q2h,4) - 1024*l1k*m4*S*pow(Q2h,4) +
		           1408*l1k*m2*Q2e*S*pow(Q2h,4) + 256*m4*Q2e*S*pow(Q2h,4) -
		           1024*l1k*m4*U*pow(Q2h,4) + 1664*l1k*m2*Q2e*U*pow(Q2h,4) +
		           256*m4*Q2e*U*pow(Q2h,4) + 1280*l1k*m2*S*U*pow(Q2h,4) +
		           256*m4*S*U*pow(Q2h,4) + 512*l1k*Q2e*S*U*pow(Q2h,4) -
		           256*m2*Q2e*S*U*pow(Q2h,4) + 1024*m4*pow(l1k,2)*pow(Q2h,4) -
		           5888*m2*Q2e*pow(l1k,2)*pow(Q2h,4) -
		           5120*m2*S*pow(l1k,2)*pow(Q2h,4) -
		           1792*Q2e*S*pow(l1k,2)*pow(Q2h,4) -
		           6656*m2*U*pow(l1k,2)*pow(Q2h,4) -
		           2304*Q2e*U*pow(l1k,2)*pow(Q2h,4) +
		           4096*m2*pow(l1k,3)*pow(Q2h,4) +
		           4608*Q2e*pow(l1k,3)*pow(Q2h,4) + 8*S*U*pow(Q2e,2)*pow(Q2h,4) -
		           1024*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,4) -
		           56*S*pow(Q2e,3)*pow(Q2h,4) - 16*U*pow(Q2e,3)*pow(Q2h,4) +
		           36*pow(Q2e,4)*pow(Q2h,4) + 512*l1k*m4*pow(Q2h,5) -
		           896*l1k*m2*Q2e*pow(Q2h,5) - 128*m4*Q2e*pow(Q2h,5) -
		           1408*l1k*m2*S*pow(Q2h,5) - 256*m4*S*pow(Q2h,5) -
		           512*l1k*Q2e*S*pow(Q2h,5) + 256*m2*Q2e*S*pow(Q2h,5) -
		           1664*l1k*m2*U*pow(Q2h,5) - 256*m4*U*pow(Q2h,5) -
		           512*l1k*Q2e*U*pow(Q2h,5) + 256*m2*Q2e*U*pow(Q2h,5) +
		           256*m2*S*U*pow(Q2h,5) + 4096*m2*pow(l1k,2)*pow(Q2h,5) +
		           1280*Q2e*pow(l1k,2)*pow(Q2h,5) + 896*l1k*m2*pow(Q2h,6) +
		           128*m4*pow(Q2h,6) + 256*l1k*Q2e*pow(Q2h,6) -
		           128*m2*Q2e*pow(Q2h,6) - 256*m2*S*pow(Q2h,6) -
		           256*m2*U*pow(Q2h,6) + 128*m2*pow(Q2h,7) -
		           1024*m4*Q2e*Q2h*pow(l1k,2)*pow(S,2) -
		           1536*m2*Q2e*Q2h*pow(l1k,3)*pow(S,2) -
		           1024*m2*Q2e*pow(l1k,4)*pow(S,2) -
		           3584*Q2e*Q2h*pow(l1k,4)*pow(S,2) -
		           2048*Q2e*pow(l1k,5)*pow(S,2) -
		           2560*Q2h*pow(l1k,3)*pow(Q2e,2)*pow(S,2) -
		           1536*pow(l1k,4)*pow(Q2e,2)*pow(S,2) -
		           384*Q2h*pow(l1k,2)*pow(Q2e,3)*pow(S,2) -
		           256*pow(l1k,3)*pow(Q2e,3)*pow(S,2) -
		           16*l1k*Q2h*pow(Q2e,4)*pow(S,2) -
		           16*pow(l1k,2)*pow(Q2e,4)*pow(S,2) - 4*Q2h*pow(Q2e,5)*pow(S,2) -
		           512*l1k*m4*Q2e*pow(Q2h,2)*pow(S,2) -
		           2432*m2*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(S,2) +
		           512*m2*pow(l1k,3)*pow(Q2h,2)*pow(S,2) +
		           256*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(S,2) -
		           1200*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) -
		           160*l1k*pow(Q2e,3)*pow(Q2h,2)*pow(S,2) -
		           8*pow(Q2e,4)*pow(Q2h,2)*pow(S,2) +
		           512*l1k*m4*pow(Q2h,3)*pow(S,2) -
		           512*l1k*m2*Q2e*pow(Q2h,3)*pow(S,2) -
		           128*m4*Q2e*pow(Q2h,3)*pow(S,2) +
		           1152*m2*pow(l1k,2)*pow(Q2h,3)*pow(S,2) +
		           512*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(S,2) -
		           80*l1k*pow(Q2e,2)*pow(Q2h,3)*pow(S,2) -
		           12*pow(Q2e,3)*pow(Q2h,3)*pow(S,2) +
		           512*l1k*m2*pow(Q2h,4)*pow(S,2) + 128*m4*pow(Q2h,4)*pow(S,2) +
		           256*l1k*Q2e*pow(Q2h,4)*pow(S,2) -
		           128*m2*Q2e*pow(Q2h,4)*pow(S,2) +
		           24*pow(Q2e,2)*pow(Q2h,4)*pow(S,2) +
		           128*m2*pow(Q2h,5)*pow(S,2) -
		           1024*m4*Q2e*Q2h*pow(l1k,2)*pow(U,2) -
		           4608*m2*Q2e*Q2h*pow(l1k,3)*pow(U,2) -
		           1024*m2*Q2e*pow(l1k,4)*pow(U,2) -
		           7680*Q2e*Q2h*pow(l1k,4)*pow(U,2) -
		           6144*Q2e*pow(l1k,5)*pow(U,2) -
		           5632*Q2h*pow(l1k,3)*pow(Q2e,2)*pow(U,2) -
		           6656*pow(l1k,4)*pow(Q2e,2)*pow(U,2) -
		           768*Q2h*pow(l1k,2)*pow(Q2e,3)*pow(U,2) -
		           2816*pow(l1k,3)*pow(Q2e,3)*pow(U,2) -
		           96*l1k*Q2h*pow(Q2e,4)*pow(U,2) -
		           656*pow(l1k,2)*pow(Q2e,4)*pow(U,2) -
		           80*l1k*pow(Q2e,5)*pow(U,2) - 20*Q2h*pow(Q2e,5)*pow(U,2) -
		           4*pow(Q2e,6)*pow(U,2) - 512*l1k*m4*Q2e*pow(Q2h,2)*pow(U,2) -
		           3968*m2*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(U,2) +
		           1536*m2*pow(l1k,3)*pow(Q2h,2)*pow(U,2) +
		           2304*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(U,2) -
		           944*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) -
		           80*l1k*pow(Q2e,3)*pow(Q2h,2)*pow(U,2) +
		           16*pow(Q2e,4)*pow(Q2h,2)*pow(U,2) +
		           512*l1k*m4*pow(Q2h,3)*pow(U,2) -
		           768*l1k*m2*Q2e*pow(Q2h,3)*pow(U,2) -
		           128*m4*Q2e*pow(Q2h,3)*pow(U,2) +
		           2688*m2*pow(l1k,2)*pow(Q2h,3)*pow(U,2) +
		           1024*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(U,2) +
		           4*pow(Q2e,3)*pow(Q2h,3)*pow(U,2) +
		           768*l1k*m2*pow(Q2h,4)*pow(U,2) + 128*m4*pow(Q2h,4)*pow(U,2) +
		           256*l1k*Q2e*pow(Q2h,4)*pow(U,2) -
		           128*m2*Q2e*pow(Q2h,4)*pow(U,2) +
		           4*pow(Q2e,2)*pow(Q2h,4)*pow(U,2) + 128*m2*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (8*S*U*pow(Q2e,5)*pow(Q2h,2) - 8*S*pow(Q2e,6)*pow(Q2h,2) -
		           16*S*U*pow(Q2e,4)*pow(Q2h,3) + 16*S*pow(Q2e,5)*pow(Q2h,3) -
		           8*U*pow(Q2e,5)*pow(Q2h,3) + 4*pow(Q2e,6)*pow(Q2h,3) +
		           8*S*U*pow(Q2e,3)*pow(Q2h,4) - 8*S*pow(Q2e,4)*pow(Q2h,4) +
		           16*U*pow(Q2e,4)*pow(Q2h,4) - 8*pow(Q2e,5)*pow(Q2h,4) -
		           8*U*pow(Q2e,3)*pow(Q2h,5) + 4*pow(Q2e,4)*pow(Q2h,5) +
		           4*Q2h*pow(Q2e,6)*pow(S,2) - 8*pow(Q2e,5)*pow(Q2h,2)*pow(S,2) +
		           4*pow(Q2e,4)*pow(Q2h,3)*pow(S,2) +
		           4*pow(Q2e,4)*pow(Q2h,3)*pow(U,2) -
		           8*pow(Q2e,3)*pow(Q2h,4)*pow(U,2) +
		           4*pow(Q2e,2)*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        (-256*l1k*m2*M2*Q2e*Q2h - 64*l1k*m2*Q2e*Q2h*S -
		           320*m4*Q2e*Q2h*S - 960*l1k*m2*Q2e*Q2h*U - 320*m4*Q2e*Q2h*U +
		           256*l1k*m2*Q2e*S*U + 128*m4*Q2e*S*U + 768*l1k*m2*Q2h*S*U +
		           384*m4*Q2h*S*U + 576*l1k*Q2e*Q2h*S*U + 768*m2*Q2e*Q2h*S*U -
		           512*m2*M2*Q2e*pow(l1k,2) + 256*M2*Q2e*Q2h*pow(l1k,2) -
		           1024*m2*Q2h*S*pow(l1k,2) - 2016*Q2e*Q2h*S*pow(l1k,2) -
		           1024*m2*Q2h*U*pow(l1k,2) - 4448*Q2e*Q2h*U*pow(l1k,2) +
		           512*m2*S*U*pow(l1k,2) + 1152*Q2e*S*U*pow(l1k,2) +
		           3840*Q2h*S*U*pow(l1k,2) - 1024*M2*Q2e*pow(l1k,3) -
		           3072*Q2h*S*pow(l1k,3) - 5120*Q2h*U*pow(l1k,3) +
		           2048*S*U*pow(l1k,3) - 440*l1k*Q2h*S*pow(Q2e,2) -
		           1768*l1k*Q2h*U*pow(Q2e,2) + 320*l1k*S*U*pow(Q2e,2) +
		           248*Q2h*S*U*pow(Q2e,2) - 768*M2*pow(l1k,2)*pow(Q2e,2) -
		           128*l1k*M2*pow(Q2e,3) - 32*M2*Q2h*pow(Q2e,3) -
		           32*Q2h*S*pow(Q2e,3) - 188*Q2h*U*pow(Q2e,3) + 8*S*U*pow(Q2e,3) -
		           8*M2*pow(Q2e,4) + 256*l1k*m2*M2*pow(Q2h,2) +
		           512*l1k*m2*Q2e*pow(Q2h,2) + 256*m4*Q2e*pow(Q2h,2) +
		           128*l1k*M2*Q2e*pow(Q2h,2) - 64*m2*M2*Q2e*pow(Q2h,2) -
		           1024*l1k*m2*S*pow(Q2h,2) - 192*m4*S*pow(Q2h,2) -
		           232*l1k*Q2e*S*pow(Q2h,2) - 1168*m2*Q2e*S*pow(Q2h,2) -
		           1024*l1k*m2*U*pow(Q2h,2) - 192*m4*U*pow(Q2h,2) -
		           376*l1k*Q2e*U*pow(Q2h,2) - 1328*m2*Q2e*U*pow(Q2h,2) +
		           1280*l1k*S*U*pow(Q2h,2) + 192*m2*S*U*pow(Q2h,2) +
		           400*Q2e*S*U*pow(Q2h,2) + 1024*m2*pow(l1k,2)*pow(Q2h,2) +
		           3136*Q2e*pow(l1k,2)*pow(Q2h,2) - 5120*S*pow(l1k,2)*pow(Q2h,2) -
		           6656*U*pow(l1k,2)*pow(Q2h,2) + 4096*pow(l1k,3)*pow(Q2h,2) +
		           1056*l1k*pow(Q2e,2)*pow(Q2h,2) +
		           40*M2*pow(Q2e,2)*pow(Q2h,2) - 340*S*pow(Q2e,2)*pow(Q2h,2) -
		           332*U*pow(Q2e,2)*pow(Q2h,2) + 116*pow(Q2e,3)*pow(Q2h,2) +
		           512*l1k*m2*pow(Q2h,3) + 64*m2*M2*pow(Q2h,3) -
		           32*l1k*Q2e*pow(Q2h,3) + 880*m2*Q2e*pow(Q2h,3) -
		           1408*l1k*S*pow(Q2h,3) - 16*m2*S*pow(Q2h,3) -
		           500*Q2e*S*pow(Q2h,3) - 1664*l1k*U*pow(Q2h,3) +
		           208*m2*U*pow(Q2h,3) - 336*Q2e*U*pow(Q2h,3) +
		           256*S*U*pow(Q2h,3) + 4096*pow(l1k,2)*pow(Q2h,3) +
		           212*pow(Q2e,2)*pow(Q2h,3) + 896*l1k*pow(Q2h,4) -
		           208*m2*pow(Q2h,4) + 208*Q2e*pow(Q2h,4) - 256*S*pow(Q2h,4) -
		           256*U*pow(Q2h,4) + 128*pow(Q2h,5) - 128*l1k*m2*Q2e*pow(S,2) +
		           64*m4*Q2e*pow(S,2) + 256*l1k*m2*Q2h*pow(S,2) +
		           192*m4*Q2h*pow(S,2) - 144*l1k*Q2e*Q2h*pow(S,2) +
		           352*m2*Q2e*Q2h*pow(S,2) + 256*m2*pow(l1k,2)*pow(S,2) +
		           320*Q2e*pow(l1k,2)*pow(S,2) + 1152*Q2h*pow(l1k,2)*pow(S,2) +
		           512*pow(l1k,3)*pow(S,2) + 48*Q2h*pow(Q2e,2)*pow(S,2) +
		           4*pow(Q2e,3)*pow(S,2) + 512*l1k*pow(Q2h,2)*pow(S,2) +
		           160*m2*pow(Q2h,2)*pow(S,2) + 284*Q2e*pow(Q2h,2)*pow(S,2) +
		           128*pow(Q2h,3)*pow(S,2) + 384*l1k*m2*Q2e*pow(U,2) +
		           64*m4*Q2e*pow(U,2) + 512*l1k*m2*Q2h*pow(U,2) +
		           192*m4*Q2h*pow(U,2) + 336*l1k*Q2e*Q2h*pow(U,2) +
		           480*m2*Q2e*Q2h*pow(U,2) + 256*m2*pow(l1k,2)*pow(U,2) +
		           1600*Q2e*pow(l1k,2)*pow(U,2) + 2688*Q2h*pow(l1k,2)*pow(U,2) +
		           1536*pow(l1k,3)*pow(U,2) + 704*l1k*pow(Q2e,2)*pow(U,2) +
		           120*Q2h*pow(Q2e,2)*pow(U,2) + 84*pow(Q2e,3)*pow(U,2) +
		           768*l1k*pow(Q2h,2)*pow(U,2) - 32*m2*pow(Q2h,2)*pow(U,2) +
		           116*Q2e*pow(Q2h,2)*pow(U,2) + 128*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*(8*Q2h*S*U*pow(Q2e,3) - 4*M2*Q2h*pow(Q2e,4) -
		           6*Q2h*U*pow(Q2e,4) + 32*m4*Q2e*S*pow(Q2h,2) -
		           32*m4*Q2e*U*pow(Q2h,2) - 64*m2*Q2e*S*U*pow(Q2h,2) +
		           24*S*U*pow(Q2e,2)*pow(Q2h,2) + 8*M2*pow(Q2e,3)*pow(Q2h,2) -
		           14*S*pow(Q2e,3)*pow(Q2h,2) - 22*U*pow(Q2e,3)*pow(Q2h,2) +
		           4*pow(Q2e,4)*pow(Q2h,2) - 32*m4*S*pow(Q2h,3) +
		           104*m2*Q2e*S*pow(Q2h,3) + 32*m4*U*pow(Q2h,3) +
		           24*m2*Q2e*U*pow(Q2h,3) - 80*Q2e*S*U*pow(Q2h,3) -
		           4*M2*pow(Q2e,2)*pow(Q2h,3) - 18*S*pow(Q2e,2)*pow(Q2h,3) -
		           4*U*pow(Q2e,2)*pow(Q2h,3) + 12*pow(Q2e,3)*pow(Q2h,3) -
		           32*m2*Q2e*pow(Q2h,4) - 40*m2*S*pow(Q2h,4) +
		           80*Q2e*S*pow(Q2h,4) + 40*m2*U*pow(Q2h,4) +
		           80*Q2e*U*pow(Q2h,4) - 40*Q2e*pow(Q2h,5) -
		           32*m4*Q2e*Q2h*pow(S,2) + 32*m4*pow(Q2h,2)*pow(S,2) -
		           72*m2*Q2e*pow(Q2h,2)*pow(S,2) +
		           10*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) +
		           40*m2*pow(Q2h,3)*pow(S,2) - 34*Q2e*pow(Q2h,3)*pow(S,2) +
		           32*m4*Q2e*Q2h*pow(U,2) + 10*Q2h*pow(Q2e,3)*pow(U,2) +
		           2*pow(Q2e,4)*pow(U,2) - 32*m4*pow(Q2h,2)*pow(U,2) +
		           8*m2*Q2e*pow(Q2h,2)*pow(U,2) +
		           4*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) - 40*m2*pow(Q2h,3)*pow(U,2) -
		           40*Q2e*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*Q2h*S*U*pow(Q2e,3) + 8*M2*Q2h*pow(Q2e,4) +
		           12*Q2h*S*pow(Q2e,4) - 48*S*U*pow(Q2e,2)*pow(Q2h,2) -
		           16*M2*pow(Q2e,3)*pow(Q2h,2) + 44*S*pow(Q2e,3)*pow(Q2h,2) +
		           28*U*pow(Q2e,3)*pow(Q2h,2) - 8*pow(Q2e,4)*pow(Q2h,2) +
		           8*M2*pow(Q2e,2)*pow(Q2h,3) + 8*S*pow(Q2e,2)*pow(Q2h,3) +
		           36*U*pow(Q2e,2)*pow(Q2h,3) - 24*pow(Q2e,3)*pow(Q2h,3) -
		           20*Q2h*pow(Q2e,3)*pow(S,2) - 4*pow(Q2e,4)*pow(S,2) -
		           8*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) -
		           20*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) - 12*Q2e*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-2)*(32*m4*Q2e*S*U*pow(Q2h,2) - 32*m4*Q2e*S*pow(Q2h,3) -
		           32*m4*Q2e*U*pow(Q2h,3) - 32*m4*S*U*pow(Q2h,3) +
		           40*m2*Q2e*S*U*pow(Q2h,3) + 16*m4*Q2e*pow(Q2h,4) +
		           32*m4*S*pow(Q2h,4) - 40*m2*Q2e*S*pow(Q2h,4) +
		           32*m4*U*pow(Q2h,4) - 40*m2*Q2e*U*pow(Q2h,4) -
		           40*m2*S*U*pow(Q2h,4) - 16*m4*pow(Q2h,5) +
		           20*m2*Q2e*pow(Q2h,5) + 40*m2*S*pow(Q2h,5) +
		           40*m2*U*pow(Q2h,5) - 20*m2*pow(Q2h,6) +
		           16*m4*Q2e*pow(Q2h,2)*pow(S,2) - 16*m4*pow(Q2h,3)*pow(S,2) +
		           20*m2*Q2e*pow(Q2h,3)*pow(S,2) - 20*m2*pow(Q2h,4)*pow(S,2) +
		           16*m4*Q2e*pow(Q2h,2)*pow(U,2) - 16*m4*pow(Q2h,3)*pow(U,2) +
		           20*m2*Q2e*pow(Q2h,3)*pow(U,2) - 20*m2*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1)) +
		     SI.B0_qmm(Q2e,m2)*(32*M2*Q2e*pow(4*m2 + Q2e,-1) +
		        32*M2*Q2h*pow(Q2e - Q2h,-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (8*M2*Q2e + 8*M2*Q2h + 8*Q2h*S + 8*Q2h*U - 16*pow(Q2h,2)) +
		        pow(l1k,-1)*(-4*M2*Q2e - 4*M2*Q2h - 4*Q2h*S - 4*Q2h*U +
		           8*pow(Q2h,2)) + pow(l1k,-1)*pow(4*m2 + Q2e,-1)*
		         (-4*M2*Q2e*Q2h + 4*Q2e*Q2h*S - 4*Q2e*Q2h*U + 8*Q2h*S*U +
		           4*M2*pow(Q2e,2) - 4*Q2h*pow(S,2) - 4*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (8*Q2h*S*U - 8*M2*pow(Q2h,2) - 44*S*pow(Q2h,2) -
		           32*U*pow(Q2h,2) + 36*pow(Q2h,3) + 12*Q2h*pow(S,2) +
		           4*Q2h*pow(U,2)) + pow(l1k,-1)*pow(Q2e - Q2h,-1)*
		         (12*S*pow(Q2h,2) - 12*U*pow(Q2h,2) - 8*Q2h*pow(S,2) +
		           8*Q2h*pow(U,2)) + pow(4*m2 + Q2e,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (8*M2*Q2e*Q2h + 8*Q2e*Q2h*S - 8*Q2e*Q2h*U - 16*Q2h*S*U -
		           8*M2*pow(Q2e,2) + 8*Q2h*pow(S,2) + 8*Q2h*pow(U,2)) +
		        pow(Q2e - Q2h,-2)*(32*Q2h*S*U - 16*S*pow(Q2h,2) - 16*U*pow(Q2h,2) +
		           16*Q2h*pow(S,2) + 16*Q2h*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-28*S*U*pow(Q2h,2) + 42*S*pow(Q2h,3) + 38*U*pow(Q2h,3) -
		           26*pow(Q2h,4) - 16*pow(Q2h,2)*pow(S,2) - 12*pow(Q2h,2)*pow(U,2)) \
		+ pow(l1k,-1)*pow(4*m2 + Q2e,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (8*S*U*pow(Q2h,2) - 4*pow(Q2h,2)*pow(S,2) - 4*pow(Q2h,2)*pow(U,2)) \
		+ pow(l1k,-2)*pow(Q2e - Q2h,-1)*
		         (12*S*U*pow(Q2h,2) - 18*S*pow(Q2h,3) - 14*U*pow(Q2h,3) +
		           10*pow(Q2h,4) + 8*pow(Q2h,2)*pow(S,2) + 4*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-1)*pow(Q2e - Q2h,-2)*
		         (8*S*pow(Q2h,3) - 8*U*pow(Q2h,3) - 8*pow(Q2h,2)*pow(S,2) +
		           8*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-3)*pow(Q2e - Q2h,-1)*
		         (-4*S*U*pow(Q2h,3) + 4*S*pow(Q2h,4) + 4*U*pow(Q2h,4) -
		           2*pow(Q2h,5) - 2*pow(Q2h,3)*pow(S,2) - 2*pow(Q2h,3)*pow(U,2)) +
		        pow(l1k,-2)*pow(Q2e - Q2h,-2)*
		         (8*S*U*pow(Q2h,3) - 8*S*pow(Q2h,4) - 8*U*pow(Q2h,4) +
		           4*pow(Q2h,5) + 4*pow(Q2h,3)*pow(S,2) + 4*pow(Q2h,3)*pow(U,2)) +
		        pow(l1k,-3)*pow(2*l1k + Q2e - Q2h,-1)*
		         (12*S*U*pow(Q2h,3) - 12*S*pow(Q2h,4) - 12*U*pow(Q2h,4) +
		           6*pow(Q2h,5) + 6*pow(Q2h,3)*pow(S,2) + 6*pow(Q2h,3)*pow(U,2)) +
		        (-16*Q2e*Q2h*S*U + 32*Q2e*S*pow(Q2h,2) + 32*Q2e*U*pow(Q2h,2) -
		           32*Q2e*pow(Q2h,3) - 8*Q2e*Q2h*pow(S,2) - 8*Q2e*Q2h*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*(24*Q2e*S*U*pow(Q2h,2) + 8*U*pow(Q2e,2)*pow(Q2h,2) -
		           40*Q2e*S*pow(Q2h,3) - 40*Q2e*U*pow(Q2h,3) - 8*S*U*pow(Q2h,3) -
		           4*pow(Q2e,2)*pow(Q2h,3) + 28*Q2e*pow(Q2h,4) + 8*S*pow(Q2h,4) +
		           8*U*pow(Q2h,4) - 4*pow(Q2h,5) + 16*Q2e*pow(Q2h,2)*pow(S,2) -
		           4*pow(Q2h,3)*pow(S,2) - 4*Q2h*pow(Q2e,2)*pow(U,2) +
		           12*Q2e*pow(Q2h,2)*pow(U,2) - 4*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (16*Q2e*S*U*pow(Q2h,2) - 16*S*pow(Q2e,2)*pow(Q2h,2) -
		           16*Q2e*S*pow(Q2h,3) - 16*Q2e*U*pow(Q2h,3) + 16*S*U*pow(Q2h,3) +
		           8*pow(Q2e,2)*pow(Q2h,3) + 8*Q2e*pow(Q2h,4) - 16*S*pow(Q2h,4) -
		           16*U*pow(Q2h,4) + 8*pow(Q2h,5) + 8*Q2h*pow(Q2e,2)*pow(S,2) +
		           8*Q2e*pow(Q2h,2)*pow(S,2) + 8*pow(Q2h,3)*pow(S,2) +
		           8*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-2)*(-16*m2*Q2e*S*U*pow(Q2h,2) +
		           24*m2*Q2e*S*pow(Q2h,3) + 24*m2*Q2e*U*pow(Q2h,3) -
		           16*Q2e*S*U*pow(Q2h,3) - 16*m2*Q2e*pow(Q2h,4) -
		           8*m2*S*pow(Q2h,4) + 16*Q2e*S*pow(Q2h,4) - 8*m2*U*pow(Q2h,4) +
		           16*Q2e*U*pow(Q2h,4) + 8*m2*pow(Q2h,5) - 8*Q2e*pow(Q2h,5) -
		           8*m2*Q2e*pow(Q2h,2)*pow(S,2) - 8*Q2e*pow(Q2h,3)*pow(S,2) -
		           8*m2*Q2e*pow(Q2h,2)*pow(U,2) - 8*Q2e*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*S*U*pow(Q2h,4) + 8*S*pow(Q2h,5) + 8*U*pow(Q2h,5) -
		           4*pow(Q2h,6) - 4*pow(Q2h,4)*pow(S,2) - 4*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-3)*(8*m2*Q2e*S*U*pow(Q2h,3) - 8*m2*Q2e*S*pow(Q2h,4) -
		           8*m2*Q2e*U*pow(Q2h,4) - 8*m2*S*U*pow(Q2h,4) +
		           4*m2*Q2e*pow(Q2h,5) + 8*m2*S*pow(Q2h,5) + 8*m2*U*pow(Q2h,5) -
		           4*m2*pow(Q2h,6) + 4*m2*Q2e*pow(Q2h,3)*pow(S,2) -
		           4*m2*pow(Q2h,4)*pow(S,2) + 4*m2*Q2e*pow(Q2h,3)*pow(U,2) -
		           4*m2*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1)) +
		     SI.C0_m0M0mm(-2*l1k + m2,m2)*
		      (48*l1k*M2 - 16*m2*M2 - 16*M2*Q2e - 40*M2*Q2h - 4*m2*S -
		        76*Q2h*S - 4*m2*U + 8*Q2h*U + 36*pow(Q2h,2) + 32*pow(S,2) -
		        16*pow(U,2) + pow(l1k,-1)*
		         (4*m4*S + 6*m2*Q2h*U - 2*m2*pow(Q2h,2) - 4*m2*pow(U,2)) +
		        pow(l1k,-2)*(-2*m4*Q2e*Q2h + 6*m4*Q2h*S + 2*m4*Q2e*U +
		           8*m4*Q2h*U - 4*m4*S*U - 4*m4*pow(Q2h,2) - 4*m4*pow(U,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (48*m2*M2*Q2e - 4*m2*Q2e*Q2h + 24*M2*Q2e*Q2h - 8*m4*S +
		           4*m2*Q2e*S + 124*m2*Q2h*S + 84*Q2e*Q2h*S + 60*m2*Q2h*U +
		           24*Q2e*Q2h*U - 8*m2*S*U - 104*Q2h*S*U + 24*M2*pow(Q2e,2) -
		           76*m2*pow(Q2h,2) - 52*Q2e*pow(Q2h,2) + 216*S*pow(Q2h,2) +
		           148*U*pow(Q2h,2) - 128*pow(Q2h,3) - 32*m2*pow(S,2) -
		           32*Q2e*pow(S,2) - 88*Q2h*pow(S,2) - 24*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (4*m4*Q2h*S - 4*m4*Q2h*U - 8*m4*S*U + 92*m2*Q2h*S*U +
		           4*m4*pow(Q2h,2) - 126*m2*S*pow(Q2h,2) - 90*m2*U*pow(Q2h,2) +
		           128*S*U*pow(Q2h,2) + 62*m2*pow(Q2h,3) - 180*S*pow(Q2h,3) -
		           148*U*pow(Q2h,3) + 100*pow(Q2h,4) - 8*m4*pow(S,2) +
		           64*m2*Q2h*pow(S,2) + 80*pow(Q2h,2)*pow(S,2) +
		           28*m2*Q2h*pow(U,2) + 48*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*m4*Q2h*S*U + 8*m4*S*pow(Q2h,2) + 8*m4*U*pow(Q2h,2) -
		           16*m2*S*U*pow(Q2h,2) - 4*m4*pow(Q2h,3) + 16*m2*S*pow(Q2h,3) +
		           16*m2*U*pow(Q2h,3) - 36*S*U*pow(Q2h,3) - 8*m2*pow(Q2h,4) +
		           36*S*pow(Q2h,4) + 36*U*pow(Q2h,4) - 18*pow(Q2h,5) -
		           4*m4*Q2h*pow(S,2) - 8*m2*pow(Q2h,2)*pow(S,2) -
		           18*pow(Q2h,3)*pow(S,2) - 4*m4*Q2h*pow(U,2) -
		           8*m2*pow(Q2h,2)*pow(U,2) - 18*pow(Q2h,3)*pow(U,2)) +
		        (512*m4*Q2e*Q2h*S*U*pow(l1k,2) + 1024*m4*Q2e*Q2h*S*pow(l1k,3) +
		           1024*m4*Q2e*Q2h*U*pow(l1k,3) - 512*m4*Q2e*S*U*pow(l1k,3) +
		           512*m4*Q2h*S*U*pow(l1k,3) - 1664*m2*Q2e*Q2h*S*U*pow(l1k,3) +
		           2048*m2*Q2e*Q2h*U*pow(l1k,4) - 512*m2*Q2e*S*U*pow(l1k,4) -
		           256*Q2e*Q2h*S*U*pow(l1k,4) + 2048*Q2e*Q2h*S*pow(l1k,5) +
		           6144*Q2e*Q2h*U*pow(l1k,5) - 1536*Q2e*S*U*pow(l1k,5) -
		           512*l1k*Q2e*Q2h*S*U*pow(m,6) -
		           640*Q2h*S*U*pow(l1k,3)*pow(Q2e,2) +
		           1024*Q2h*S*pow(l1k,4)*pow(Q2e,2) +
		           4096*Q2h*U*pow(l1k,4)*pow(Q2e,2) -
		           896*S*U*pow(l1k,4)*pow(Q2e,2) -
		           32*Q2h*S*U*pow(l1k,2)*pow(Q2e,3) +
		           64*Q2h*S*pow(l1k,3)*pow(Q2e,3) +
		           576*Q2h*U*pow(l1k,3)*pow(Q2e,3) -
		           64*S*U*pow(l1k,3)*pow(Q2e,3) + 32*Q2h*U*pow(l1k,2)*pow(Q2e,4) -
		           512*l1k*m4*Q2e*S*U*pow(Q2h,2) -
		           512*m4*Q2e*S*pow(l1k,2)*pow(Q2h,2) +
		           512*m2*Q2e*S*U*pow(l1k,2)*pow(Q2h,2) -
		           1024*m4*Q2e*pow(l1k,3)*pow(Q2h,2) +
		           2560*m2*Q2e*S*pow(l1k,3)*pow(Q2h,2) +
		           3584*m2*Q2e*U*pow(l1k,3)*pow(Q2h,2) +
		           384*m2*S*U*pow(l1k,3)*pow(Q2h,2) -
		           128*Q2e*S*U*pow(l1k,3)*pow(Q2h,2) -
		           1024*m2*Q2e*pow(l1k,4)*pow(Q2h,2) -
		           512*Q2e*S*pow(l1k,4)*pow(Q2h,2) -
		           512*Q2e*U*pow(l1k,4)*pow(Q2h,2) -
		           4096*Q2e*pow(l1k,5)*pow(Q2h,2) +
		           512*l1k*Q2e*S*pow(m,6)*pow(Q2h,2) +
		           512*l1k*Q2e*U*pow(m,6)*pow(Q2h,2) +
		           512*l1k*S*U*pow(m,6)*pow(Q2h,2) +
		           768*S*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) +
		           1024*U*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) -
		           2560*pow(l1k,4)*pow(Q2e,2)*pow(Q2h,2) -
		           16*l1k*S*U*pow(Q2e,3)*pow(Q2h,2) +
		           32*U*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2) -
		           320*pow(l1k,3)*pow(Q2e,3)*pow(Q2h,2) +
		           16*l1k*S*pow(Q2e,4)*pow(Q2h,2) + 8*S*U*pow(Q2e,4)*pow(Q2h,2) -
		           16*pow(l1k,2)*pow(Q2e,4)*pow(Q2h,2) -
		           8*S*pow(Q2e,5)*pow(Q2h,2) + 640*l1k*m4*Q2e*S*pow(Q2h,3) +
		           640*l1k*m4*Q2e*U*pow(Q2h,3) + 256*l1k*m4*S*U*pow(Q2h,3) -
		           544*l1k*m2*Q2e*S*U*pow(Q2h,3) + 128*m4*Q2e*S*U*pow(Q2h,3) -
		           1024*m2*Q2e*S*pow(l1k,2)*pow(Q2h,3) -
		           512*m4*U*pow(l1k,2)*pow(Q2h,3) -
		           256*m2*Q2e*U*pow(l1k,2)*pow(Q2h,3) +
		           256*m2*S*U*pow(l1k,2)*pow(Q2h,3) +
		           576*Q2e*S*U*pow(l1k,2)*pow(Q2h,3) -
		           2560*m2*Q2e*pow(l1k,3)*pow(Q2h,3) -
		           512*m2*S*pow(l1k,3)*pow(Q2h,3) +
		           768*Q2e*S*pow(l1k,3)*pow(Q2h,3) -
		           1536*m2*U*pow(l1k,3)*pow(Q2h,3) -
		           512*Q2e*U*pow(l1k,3)*pow(Q2h,3) +
		           1024*Q2e*pow(l1k,4)*pow(Q2h,3) -
		           256*l1k*Q2e*pow(m,6)*pow(Q2h,3) -
		           512*l1k*S*pow(m,6)*pow(Q2h,3) -
		           512*l1k*U*pow(m,6)*pow(Q2h,3) -
		           576*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,3) +
		           16*l1k*U*pow(Q2e,3)*pow(Q2h,3) - 8*S*U*pow(Q2e,3)*pow(Q2h,3) -
		           8*l1k*pow(Q2e,4)*pow(Q2h,3) + 8*S*pow(Q2e,4)*pow(Q2h,3) -
		           8*U*pow(Q2e,4)*pow(Q2h,3) + 4*pow(Q2e,5)*pow(Q2h,3) -
		           384*l1k*m4*Q2e*pow(Q2h,4) - 384*l1k*m4*S*pow(Q2h,4) +
		           640*l1k*m2*Q2e*S*pow(Q2h,4) - 128*m4*Q2e*S*pow(Q2h,4) -
		           384*l1k*m4*U*pow(Q2h,4) + 576*l1k*m2*Q2e*U*pow(Q2h,4) -
		           128*m4*Q2e*U*pow(Q2h,4) + 160*l1k*m2*S*U*pow(Q2h,4) -
		           128*m4*S*U*pow(Q2h,4) - 128*l1k*Q2e*S*U*pow(Q2h,4) +
		           64*m2*Q2e*S*U*pow(Q2h,4) + 256*m4*pow(l1k,2)*pow(Q2h,4) +
		           384*m2*Q2e*pow(l1k,2)*pow(Q2h,4) -
		           128*m2*S*pow(l1k,2)*pow(Q2h,4) -
		           768*Q2e*S*pow(l1k,2)*pow(Q2h,4) -
		           640*m2*U*pow(l1k,2)*pow(Q2h,4) -
		           640*Q2e*U*pow(l1k,2)*pow(Q2h,4) +
		           1024*m2*pow(l1k,3)*pow(Q2h,4) - 64*Q2e*pow(l1k,3)*pow(Q2h,4) +
		           256*l1k*pow(m,6)*pow(Q2h,4) + 8*U*pow(Q2e,3)*pow(Q2h,4) -
		           4*pow(Q2e,4)*pow(Q2h,4) + 256*l1k*m4*pow(Q2h,5) -
		           336*l1k*m2*Q2e*pow(Q2h,5) + 64*m4*Q2e*pow(Q2h,5) -
		           256*l1k*m2*S*pow(Q2h,5) + 128*m4*S*pow(Q2h,5) +
		           128*l1k*Q2e*S*pow(Q2h,5) - 64*m2*Q2e*S*pow(Q2h,5) -
		           192*l1k*m2*U*pow(Q2h,5) + 128*m4*U*pow(Q2h,5) +
		           128*l1k*Q2e*U*pow(Q2h,5) - 64*m2*Q2e*U*pow(Q2h,5) -
		           64*m2*S*U*pow(Q2h,5) + 256*m2*pow(l1k,2)*pow(Q2h,5) +
		           416*Q2e*pow(l1k,2)*pow(Q2h,5) + 144*l1k*m2*pow(Q2h,6) -
		           64*m4*pow(Q2h,6) - 64*l1k*Q2e*pow(Q2h,6) +
		           32*m2*Q2e*pow(Q2h,6) + 64*m2*S*pow(Q2h,6) +
		           64*m2*U*pow(Q2h,6) - 32*m2*pow(Q2h,7) +
		           512*m4*Q2e*Q2h*pow(l1k,2)*pow(S,2) -
		           256*m4*Q2e*pow(l1k,3)*pow(S,2) +
		           256*m4*Q2h*pow(l1k,3)*pow(S,2) -
		           704*m2*Q2e*Q2h*pow(l1k,3)*pow(S,2) +
		           256*m2*Q2e*pow(l1k,4)*pow(S,2) +
		           384*Q2e*Q2h*pow(l1k,4)*pow(S,2) - 256*Q2e*pow(l1k,5)*pow(S,2) -
		           256*l1k*Q2e*Q2h*pow(m,6)*pow(S,2) -
		           96*Q2h*pow(l1k,3)*pow(Q2e,2)*pow(S,2) -
		           64*pow(l1k,4)*pow(Q2e,2)*pow(S,2) +
		           16*Q2h*pow(l1k,2)*pow(Q2e,3)*pow(S,2) -
		           8*l1k*Q2h*pow(Q2e,4)*pow(S,2) + 4*Q2h*pow(Q2e,5)*pow(S,2) -
		           256*l1k*m4*Q2e*pow(Q2h,2)*pow(S,2) -
		           256*m4*pow(l1k,2)*pow(Q2h,2)*pow(S,2) +
		           576*m2*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(S,2) +
		           64*m2*pow(l1k,3)*pow(Q2h,2)*pow(S,2) -
		           576*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(S,2) +
		           256*l1k*pow(m,6)*pow(Q2h,2)*pow(S,2) -
		           4*pow(Q2e,4)*pow(Q2h,2)*pow(S,2) +
		           128*l1k*m4*pow(Q2h,3)*pow(S,2) -
		           304*l1k*m2*Q2e*pow(Q2h,3)*pow(S,2) +
		           64*m4*Q2e*pow(Q2h,3)*pow(S,2) -
		           64*m2*pow(l1k,2)*pow(Q2h,3)*pow(S,2) +
		           352*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(S,2) +
		           112*l1k*m2*pow(Q2h,4)*pow(S,2) - 64*m4*pow(Q2h,4)*pow(S,2) -
		           64*l1k*Q2e*pow(Q2h,4)*pow(S,2) +
		           32*m2*Q2e*pow(Q2h,4)*pow(S,2) - 32*m2*pow(Q2h,5)*pow(S,2) -
		           256*m4*Q2e*pow(l1k,3)*pow(U,2) +
		           256*m4*Q2h*pow(l1k,3)*pow(U,2) -
		           1216*m2*Q2e*Q2h*pow(l1k,3)*pow(U,2) -
		           768*m2*Q2e*pow(l1k,4)*pow(U,2) -
		           128*Q2e*Q2h*pow(l1k,4)*pow(U,2) -
		           2304*Q2e*pow(l1k,5)*pow(U,2) -
		           256*l1k*Q2e*Q2h*pow(m,6)*pow(U,2) -
		           480*Q2h*pow(l1k,3)*pow(Q2e,2)*pow(U,2) -
		           1600*pow(l1k,4)*pow(Q2e,2)*pow(U,2) -
		           16*Q2h*pow(l1k,2)*pow(Q2e,3)*pow(U,2) -
		           256*pow(l1k,3)*pow(Q2e,3)*pow(U,2) -
		           16*pow(l1k,2)*pow(Q2e,4)*pow(U,2) -
		           256*l1k*m4*Q2e*pow(Q2h,2)*pow(U,2) +
		           256*m4*pow(l1k,2)*pow(Q2h,2)*pow(U,2) -
		           64*m2*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(U,2) +
		           576*m2*pow(l1k,3)*pow(Q2h,2)*pow(U,2) +
		           448*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(U,2) +
		           256*l1k*pow(m,6)*pow(Q2h,2)*pow(U,2) -
		           16*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) +
		           128*l1k*m4*pow(Q2h,3)*pow(U,2) -
		           240*l1k*m2*Q2e*pow(Q2h,3)*pow(U,2) +
		           64*m4*Q2e*pow(Q2h,3)*pow(U,2) +
		           320*m2*pow(l1k,2)*pow(Q2h,3)*pow(U,2) +
		           224*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(U,2) -
		           8*l1k*pow(Q2e,2)*pow(Q2h,3)*pow(U,2) +
		           4*pow(Q2e,3)*pow(Q2h,3)*pow(U,2) +
		           48*l1k*m2*pow(Q2h,4)*pow(U,2) - 64*m4*pow(Q2h,4)*pow(U,2) -
		           64*l1k*Q2e*pow(Q2h,4)*pow(U,2) +
		           32*m2*Q2e*pow(Q2h,4)*pow(U,2) -
		           4*pow(Q2e,2)*pow(Q2h,4)*pow(U,2) - 32*m2*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*S*U*pow(Q2e,5)*pow(Q2h,2) + 8*S*pow(Q2e,6)*pow(Q2h,2) +
		           16*S*U*pow(Q2e,4)*pow(Q2h,3) - 16*S*pow(Q2e,5)*pow(Q2h,3) +
		           8*U*pow(Q2e,5)*pow(Q2h,3) - 4*pow(Q2e,6)*pow(Q2h,3) -
		           8*S*U*pow(Q2e,3)*pow(Q2h,4) + 8*S*pow(Q2e,4)*pow(Q2h,4) -
		           16*U*pow(Q2e,4)*pow(Q2h,4) + 8*pow(Q2e,5)*pow(Q2h,4) +
		           8*U*pow(Q2e,3)*pow(Q2h,5) - 4*pow(Q2e,4)*pow(Q2h,5) -
		           4*Q2h*pow(Q2e,6)*pow(S,2) + 8*pow(Q2e,5)*pow(Q2h,2)*pow(S,2) -
		           4*pow(Q2e,4)*pow(Q2h,3)*pow(S,2) -
		           4*pow(Q2e,4)*pow(Q2h,3)*pow(U,2) +
		           8*pow(Q2e,3)*pow(Q2h,4)*pow(U,2) -
		           4*pow(Q2e,2)*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        (-128*l1k*m4*M2*Q2e + 128*l1k*m4*M2*Q2h -
		           96*l1k*m2*M2*Q2e*Q2h - 320*l1k*m4*Q2h*S -
		           344*l1k*m2*Q2e*Q2h*S + 64*m4*Q2e*Q2h*S - 320*l1k*m4*Q2h*U -
		           264*l1k*m2*Q2e*Q2h*U - 32*m4*Q2e*Q2h*U + 128*l1k*m4*S*U +
		           64*l1k*m2*Q2e*S*U + 256*l1k*m2*Q2h*S*U - 64*m4*Q2h*S*U +
		           240*l1k*Q2e*Q2h*S*U - 264*m2*Q2e*Q2h*S*U +
		           128*m2*M2*Q2e*pow(l1k,2) + 192*M2*Q2e*Q2h*pow(l1k,2) +
		           32*m2*Q2h*S*pow(l1k,2) + 304*Q2e*Q2h*S*pow(l1k,2) -
		           480*m2*Q2h*U*pow(l1k,2) - 320*Q2e*Q2h*U*pow(l1k,2) +
		           128*m2*S*U*pow(l1k,2) + 32*Q2e*S*U*pow(l1k,2) +
		           256*Q2h*S*U*pow(l1k,2) - 128*M2*Q2e*pow(l1k,3) -
		           512*Q2h*S*pow(l1k,3) - 1536*Q2h*U*pow(l1k,3) +
		           384*S*U*pow(l1k,3) - 16*l1k*M2*Q2h*pow(Q2e,2) -
		           24*l1k*Q2h*S*pow(Q2e,2) - 72*l1k*Q2h*U*pow(Q2e,2) -
		           16*Q2h*S*U*pow(Q2e,2) - 32*M2*pow(l1k,2)*pow(Q2e,2) +
		           8*M2*Q2h*pow(Q2e,3) + 12*Q2h*S*pow(Q2e,3) +
		           256*l1k*m4*pow(Q2h,2) + 32*l1k*m2*M2*pow(Q2h,2) +
		           288*l1k*m2*Q2e*pow(Q2h,2) - 16*m4*Q2e*pow(Q2h,2) -
		           64*l1k*M2*Q2e*pow(Q2h,2) + 32*m2*M2*Q2e*pow(Q2h,2) -
		           392*l1k*m2*S*pow(Q2h,2) + 32*m4*S*pow(Q2h,2) -
		           712*l1k*Q2e*S*pow(Q2h,2) + 516*m2*Q2e*S*pow(Q2h,2) -
		           376*l1k*m2*U*pow(Q2h,2) + 64*m4*U*pow(Q2h,2) -
		           328*l1k*Q2e*U*pow(Q2h,2) + 276*m2*Q2e*U*pow(Q2h,2) +
		           160*l1k*S*U*pow(Q2h,2) - 136*m2*S*U*pow(Q2h,2) -
		           320*Q2e*S*U*pow(Q2h,2) + 256*m2*pow(l1k,2)*pow(Q2h,2) -
		           16*Q2e*pow(l1k,2)*pow(Q2h,2) - 128*S*pow(l1k,2)*pow(Q2h,2) -
		           640*U*pow(l1k,2)*pow(Q2h,2) + 1024*pow(l1k,3)*pow(Q2h,2) +
		           48*l1k*pow(Q2e,2)*pow(Q2h,2) - 8*M2*pow(Q2e,2)*pow(Q2h,2) +
		           56*S*pow(Q2e,2)*pow(Q2h,2) + 28*U*pow(Q2e,2)*pow(Q2h,2) -
		           8*pow(Q2e,3)*pow(Q2h,2) + 256*l1k*m2*pow(Q2h,3) -
		           16*m4*pow(Q2h,3) - 32*m2*M2*pow(Q2h,3) +
		           408*l1k*Q2e*pow(Q2h,3) - 268*m2*Q2e*pow(Q2h,3) -
		           256*l1k*S*pow(Q2h,3) + 52*m2*S*pow(Q2h,3) +
		           424*Q2e*S*pow(Q2h,3) - 192*l1k*U*pow(Q2h,3) +
		           148*m2*U*pow(Q2h,3) + 360*Q2e*U*pow(Q2h,3) -
		           64*S*U*pow(Q2h,3) + 256*pow(l1k,2)*pow(Q2h,3) -
		           32*pow(Q2e,2)*pow(Q2h,3) + 144*l1k*pow(Q2h,4) -
		           28*m2*pow(Q2h,4) - 232*Q2e*pow(Q2h,4) + 64*S*pow(Q2h,4) +
		           64*U*pow(Q2h,4) - 32*pow(Q2h,5) + 64*l1k*m4*pow(S,2) +
		           128*l1k*m2*Q2e*pow(S,2) - 32*m4*Q2e*pow(S,2) +
		           80*l1k*m2*Q2h*pow(S,2) - 32*m4*Q2h*pow(S,2) +
		           280*l1k*Q2e*Q2h*pow(S,2) - 240*m2*Q2e*Q2h*pow(S,2) -
		           64*m2*pow(l1k,2)*pow(S,2) - 144*Q2e*pow(l1k,2)*pow(S,2) -
		           64*Q2h*pow(l1k,2)*pow(S,2) + 64*pow(l1k,3)*pow(S,2) +
		           8*l1k*pow(Q2e,2)*pow(S,2) - 24*Q2h*pow(Q2e,2)*pow(S,2) -
		           4*pow(Q2e,3)*pow(S,2) + 112*l1k*pow(Q2h,2)*pow(S,2) -
		           32*m2*pow(Q2h,2)*pow(S,2) - 192*Q2e*pow(Q2h,2)*pow(S,2) -
		           32*pow(Q2h,3)*pow(S,2) + 64*l1k*m4*pow(U,2) +
		           64*l1k*m2*Q2e*pow(U,2) + 32*m4*Q2e*pow(U,2) +
		           112*l1k*m2*Q2h*pow(U,2) - 32*m4*Q2h*pow(U,2) +
		           8*l1k*Q2e*Q2h*pow(U,2) - 40*m2*Q2e*Q2h*pow(U,2) +
		           192*m2*pow(l1k,2)*pow(U,2) + 176*Q2e*pow(l1k,2)*pow(U,2) +
		           320*Q2h*pow(l1k,2)*pow(U,2) + 576*pow(l1k,3)*pow(U,2) +
		           40*l1k*pow(Q2e,2)*pow(U,2) + 48*l1k*pow(Q2h,2)*pow(U,2) -
		           88*m2*pow(Q2h,2)*pow(U,2) - 116*Q2e*pow(Q2h,2)*pow(U,2) -
		           32*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (16*Q2h*S*U*pow(Q2e,3) - 8*M2*Q2h*pow(Q2e,4) -
		           12*Q2h*S*pow(Q2e,4) + 48*S*U*pow(Q2e,2)*pow(Q2h,2) +
		           16*M2*pow(Q2e,3)*pow(Q2h,2) - 44*S*pow(Q2e,3)*pow(Q2h,2) -
		           28*U*pow(Q2e,3)*pow(Q2h,2) + 8*pow(Q2e,4)*pow(Q2h,2) -
		           8*M2*pow(Q2e,2)*pow(Q2h,3) - 8*S*pow(Q2e,2)*pow(Q2h,3) -
		           36*U*pow(Q2e,2)*pow(Q2h,3) + 24*pow(Q2e,3)*pow(Q2h,3) +
		           20*Q2h*pow(Q2e,3)*pow(S,2) + 4*pow(Q2e,4)*pow(S,2) +
		           8*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) +
		           20*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) + 12*Q2e*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*(96*m4*Q2e*Q2h*S*U - 128*m4*Q2e*S*pow(Q2h,2) -
		           96*m4*Q2e*U*pow(Q2h,2) - 64*m4*S*U*pow(Q2h,2) +
		           160*m2*Q2e*S*U*pow(Q2h,2) + 64*m4*Q2e*pow(Q2h,3) +
		           96*m4*S*pow(Q2h,3) - 212*m2*Q2e*S*pow(Q2h,3) +
		           64*m4*U*pow(Q2h,3) - 180*m2*Q2e*U*pow(Q2h,3) -
		           56*m2*S*U*pow(Q2h,3) + 72*Q2e*S*U*pow(Q2h,3) -
		           48*m4*pow(Q2h,4) + 116*m2*Q2e*pow(Q2h,4) +
		           108*m2*S*pow(Q2h,4) - 72*Q2e*S*pow(Q2h,4) +
		           76*m2*U*pow(Q2h,4) - 72*Q2e*U*pow(Q2h,4) - 64*m2*pow(Q2h,5) +
		           36*Q2e*pow(Q2h,5) + 64*m4*Q2e*Q2h*pow(S,2) -
		           48*m4*pow(Q2h,2)*pow(S,2) + 96*m2*Q2e*pow(Q2h,2)*pow(S,2) -
		           44*m2*pow(Q2h,3)*pow(S,2) + 36*Q2e*pow(Q2h,3)*pow(S,2) +
		           32*m4*Q2e*Q2h*pow(U,2) - 16*m4*pow(Q2h,2)*pow(U,2) +
		           64*m2*Q2e*pow(Q2h,2)*pow(U,2) - 12*m2*pow(Q2h,3)*pow(U,2) +
		           36*Q2e*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-2)*(-16*m4*Q2e*S*U*pow(Q2h,2) + 16*m4*Q2e*S*pow(Q2h,3) +
		           16*m4*Q2e*U*pow(Q2h,3) + 16*m4*S*U*pow(Q2h,3) -
		           36*m2*Q2e*S*U*pow(Q2h,3) - 8*m4*Q2e*pow(Q2h,4) -
		           16*m4*S*pow(Q2h,4) + 36*m2*Q2e*S*pow(Q2h,4) -
		           16*m4*U*pow(Q2h,4) + 36*m2*Q2e*U*pow(Q2h,4) +
		           36*m2*S*U*pow(Q2h,4) + 8*m4*pow(Q2h,5) -
		           18*m2*Q2e*pow(Q2h,5) - 36*m2*S*pow(Q2h,5) -
		           36*m2*U*pow(Q2h,5) + 18*m2*pow(Q2h,6) -
		           8*m4*Q2e*pow(Q2h,2)*pow(S,2) + 8*m4*pow(Q2h,3)*pow(S,2) -
		           18*m2*Q2e*pow(Q2h,3)*pow(S,2) + 18*m2*pow(Q2h,4)*pow(S,2) -
		           8*m4*Q2e*pow(Q2h,2)*pow(U,2) + 8*m4*pow(Q2h,3)*pow(U,2) -
		           18*m2*Q2e*pow(Q2h,3)*pow(U,2) + 18*m2*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1)) +
		     SI.D0_mm0QqMm0mm(Q2h,Q2e,2*l1k + m2 + Q2e - Q2h,m2)*
		      (-512*l1k*m2*M2 - 528*l1k*M2*Q2e - 192*m2*M2*Q2e -
		        832*l1k*M2*Q2h - 96*m2*M2*Q2h - 144*M2*Q2e*Q2h -
		        1528*l1k*Q2h*S - 400*m2*Q2h*S - 172*Q2e*Q2h*S - 3928*l1k*Q2h*U -
		        768*m2*Q2h*U - 672*Q2e*Q2h*U + 896*l1k*S*U + 192*m2*S*U +
		        80*Q2e*S*U + 952*Q2h*S*U - 1792*M2*pow(l1k,2) -
		        88*M2*pow(Q2e,2) + 2576*l1k*pow(Q2h,2) + 496*m2*pow(Q2h,2) -
		        144*M2*pow(Q2h,2) + 420*Q2e*pow(Q2h,2) - 1284*S*pow(Q2h,2) -
		        2212*U*pow(Q2h,2) + 1316*pow(Q2h,3) + 80*l1k*pow(S,2) +
		        32*m2*pow(S,2) + 176*Q2h*pow(S,2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (64*m4*M2*Q2h + 32*m2*M2*Q2e*Q2h + 64*m4*Q2h*S +
		           32*m2*Q2e*Q2h*S - 64*m4*pow(S,2) - 32*m2*Q2e*pow(S,2)) +
		        1392*l1k*pow(U,2) + 256*m2*pow(U,2) + 264*Q2e*pow(U,2) +
		        904*Q2h*pow(U,2) + pow(l1k,-1)*
		         (-32*m4*M2*Q2h - 32*m2*M2*Q2e*Q2h - 40*m2*Q2e*Q2h*S -
		           64*m4*Q2h*U - 128*m2*Q2e*Q2h*U + 160*m2*Q2h*S*U +
		           68*Q2e*Q2h*S*U - 40*m2*M2*pow(Q2e,2) -
		           20*M2*Q2h*pow(Q2e,2) - 12*Q2h*S*pow(Q2e,2) -
		           54*Q2h*U*pow(Q2e,2) - 12*M2*pow(Q2e,3) + 16*m4*pow(Q2h,2) +
		           8*m2*M2*pow(Q2h,2) + 72*m2*Q2e*pow(Q2h,2) -
		           192*m2*S*pow(Q2h,2) - 102*Q2e*S*pow(Q2h,2) -
		           272*m2*U*pow(Q2h,2) - 168*Q2e*U*pow(Q2h,2) +
		           200*S*U*pow(Q2h,2) + 34*pow(Q2e,2)*pow(Q2h,2) +
		           152*m2*pow(Q2h,3) + 100*Q2e*pow(Q2h,3) - 232*S*pow(Q2h,3) -
		           296*U*pow(Q2h,3) + 164*pow(Q2h,4) + 40*m2*Q2h*pow(S,2) +
		           12*Q2e*Q2h*pow(S,2) + 68*pow(Q2h,2)*pow(S,2) + 32*m4*pow(U,2) +
		           48*m2*Q2e*pow(U,2) + 120*m2*Q2h*pow(U,2) +
		           68*Q2e*Q2h*pow(U,2) + 20*pow(Q2e,2)*pow(U,2) +
		           132*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-1)*(4*S*U*pow(Q2e,6)*pow(Q2h,2) -
		           4*U*pow(Q2e,7)*pow(Q2h,2) - 8*S*U*pow(Q2e,5)*pow(Q2h,3) -
		           4*S*pow(Q2e,6)*pow(Q2h,3) + 8*U*pow(Q2e,6)*pow(Q2h,3) +
		           2*pow(Q2e,7)*pow(Q2h,3) + 4*S*U*pow(Q2e,4)*pow(Q2h,4) +
		           8*S*pow(Q2e,5)*pow(Q2h,4) - 4*U*pow(Q2e,5)*pow(Q2h,4) -
		           4*pow(Q2e,6)*pow(Q2h,4) - 4*S*pow(Q2e,4)*pow(Q2h,5) +
		           2*pow(Q2e,5)*pow(Q2h,5) + 2*pow(Q2e,5)*pow(Q2h,3)*pow(S,2) -
		           4*pow(Q2e,4)*pow(Q2h,4)*pow(S,2) +
		           2*pow(Q2e,3)*pow(Q2h,5)*pow(S,2) + 2*Q2h*pow(Q2e,7)*pow(U,2) -
		           4*pow(Q2e,6)*pow(Q2h,2)*pow(U,2) +
		           2*pow(Q2e,5)*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        (11264*m4*Q2e*Q2h*S*U*pow(l1k,3) -
		           4096*m4*Q2e*Q2h*S*pow(l1k,4) - 4096*m4*Q2e*Q2h*U*pow(l1k,4) +
		           2048*m4*Q2e*S*U*pow(l1k,4) +
		           52224*m2*Q2e*Q2h*S*U*pow(l1k,4) -
		           20480*m2*Q2e*Q2h*S*pow(l1k,5) -
		           28672*m2*Q2e*Q2h*U*pow(l1k,5) + 12288*m2*Q2e*S*U*pow(l1k,5) +
		           105472*Q2e*Q2h*S*U*pow(l1k,5) - 81920*Q2e*Q2h*S*pow(l1k,6) -
		           139264*Q2e*Q2h*U*pow(l1k,6) + 53248*Q2e*S*U*pow(l1k,6) +
		           2048*Q2e*Q2h*S*U*pow(l1k,2)*pow(m,6) +
		           87808*Q2h*S*U*pow(l1k,4)*pow(Q2e,2) -
		           71680*Q2h*S*pow(l1k,5)*pow(Q2e,2) -
		           141312*Q2h*U*pow(l1k,5)*pow(Q2e,2) +
		           49152*S*U*pow(l1k,5)*pow(Q2e,2) +
		           19456*Q2h*S*U*pow(l1k,3)*pow(Q2e,3) -
		           21760*Q2h*S*pow(l1k,4)*pow(Q2e,3) -
		           57600*Q2h*U*pow(l1k,4)*pow(Q2e,3) +
		           16512*S*U*pow(l1k,4)*pow(Q2e,3) +
		           3072*Q2h*S*U*pow(l1k,2)*pow(Q2e,4) -
		           3904*Q2h*S*pow(l1k,3)*pow(Q2e,4) -
		           14784*Q2h*U*pow(l1k,3)*pow(Q2e,4) +
		           3264*S*U*pow(l1k,3)*pow(Q2e,4) + 320*l1k*Q2h*S*U*pow(Q2e,5) -
		           384*Q2h*S*pow(l1k,2)*pow(Q2e,5) -
		           2336*Q2h*U*pow(l1k,2)*pow(Q2e,5) +
		           352*S*U*pow(l1k,2)*pow(Q2e,5) - 16*l1k*Q2h*S*pow(Q2e,6) -
		           208*l1k*Q2h*U*pow(Q2e,6) + 16*l1k*S*U*pow(Q2e,6) +
		           16*Q2h*S*U*pow(Q2e,6) - 8*Q2h*U*pow(Q2e,7) +
		           12544*m4*Q2e*S*U*pow(l1k,2)*pow(Q2h,2) -
		           16384*m4*Q2e*S*pow(l1k,3)*pow(Q2h,2) -
		           20480*m4*Q2e*U*pow(l1k,3)*pow(Q2h,2) -
		           3072*m4*S*U*pow(l1k,3)*pow(Q2h,2) +
		           59904*m2*Q2e*S*U*pow(l1k,3)*pow(Q2h,2) +
		           4096*m4*Q2e*pow(l1k,4)*pow(Q2h,2) -
		           73728*m2*Q2e*S*pow(l1k,4)*pow(Q2h,2) -
		           108544*m2*Q2e*U*pow(l1k,4)*pow(Q2h,2) -
		           13312*m2*S*U*pow(l1k,4)*pow(Q2h,2) +
		           15616*Q2e*S*U*pow(l1k,4)*pow(Q2h,2) +
		           24576*m2*Q2e*pow(l1k,5)*pow(Q2h,2) -
		           135168*Q2e*S*pow(l1k,5)*pow(Q2h,2) -
		           192512*Q2e*U*pow(l1k,5)*pow(Q2h,2) +
		           110592*Q2e*pow(l1k,6)*pow(Q2h,2) +
		           2048*l1k*Q2e*S*U*pow(m,6)*pow(Q2h,2) -
		           3072*Q2e*S*pow(l1k,2)*pow(m,6)*pow(Q2h,2) -
		           3072*Q2e*U*pow(l1k,2)*pow(m,6)*pow(Q2h,2) +
		           36608*S*U*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) -
		           114432*S*pow(l1k,4)*pow(Q2e,2)*pow(Q2h,2) -
		           158464*U*pow(l1k,4)*pow(Q2e,2)*pow(Q2h,2) +
		           106496*pow(l1k,5)*pow(Q2e,2)*pow(Q2h,2) +
		           4736*S*U*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2) -
		           24320*S*pow(l1k,3)*pow(Q2e,3)*pow(Q2h,2) -
		           30976*U*pow(l1k,3)*pow(Q2e,3)*pow(Q2h,2) +
		           39680*pow(l1k,4)*pow(Q2e,3)*pow(Q2h,2) +
		           512*l1k*S*U*pow(Q2e,4)*pow(Q2h,2) -
		           3648*S*pow(l1k,2)*pow(Q2e,4)*pow(Q2h,2) -
		           4512*U*pow(l1k,2)*pow(Q2e,4)*pow(Q2h,2) +
		           9344*pow(l1k,3)*pow(Q2e,4)*pow(Q2h,2) -
		           352*l1k*S*pow(Q2e,5)*pow(Q2h,2) -
		           624*l1k*U*pow(Q2e,5)*pow(Q2h,2) +
		           56*S*U*pow(Q2e,5)*pow(Q2h,2) +
		           1360*pow(l1k,2)*pow(Q2e,5)*pow(Q2h,2) +
		           112*l1k*pow(Q2e,6)*pow(Q2h,2) - 16*S*pow(Q2e,6)*pow(Q2h,2) -
		           72*U*pow(Q2e,6)*pow(Q2h,2) + 4*pow(Q2e,7)*pow(Q2h,2) +
		           4608*l1k*m4*Q2e*S*U*pow(Q2h,3) -
		           15872*m4*Q2e*S*pow(l1k,2)*pow(Q2h,3) -
		           18432*m4*Q2e*U*pow(l1k,2)*pow(Q2h,3) -
		           7936*m4*S*U*pow(l1k,2)*pow(Q2h,3) +
		           27456*m2*Q2e*S*U*pow(l1k,2)*pow(Q2h,3) +
		           14336*m4*Q2e*pow(l1k,3)*pow(Q2h,3) +
		           5120*m4*S*pow(l1k,3)*pow(Q2h,3) -
		           76288*m2*Q2e*S*pow(l1k,3)*pow(Q2h,3) +
		           7168*m4*U*pow(l1k,3)*pow(Q2h,3) -
		           104960*m2*Q2e*U*pow(l1k,3)*pow(Q2h,3) -
		           33024*m2*S*U*pow(l1k,3)*pow(Q2h,3) -
		           20864*Q2e*S*U*pow(l1k,3)*pow(Q2h,3) +
		           75776*m2*Q2e*pow(l1k,4)*pow(Q2h,3) +
		           20480*m2*S*pow(l1k,4)*pow(Q2h,3) -
		           9472*Q2e*S*pow(l1k,4)*pow(Q2h,3) +
		           34816*m2*U*pow(l1k,4)*pow(Q2h,3) +
		           5888*Q2e*U*pow(l1k,4)*pow(Q2h,3) +
		           112640*Q2e*pow(l1k,5)*pow(Q2h,3) -
		           2048*l1k*Q2e*S*pow(m,6)*pow(Q2h,3) -
		           2048*l1k*Q2e*U*pow(m,6)*pow(Q2h,3) -
		           2048*l1k*S*U*pow(m,6)*pow(Q2h,3) +
		           256*Q2e*S*U*pow(m,6)*pow(Q2h,3) +
		           2048*Q2e*pow(l1k,2)*pow(m,6)*pow(Q2h,3) +
		           1024*S*pow(l1k,2)*pow(m,6)*pow(Q2h,3) +
		           1024*U*pow(l1k,2)*pow(m,6)*pow(Q2h,3) +
		           5440*S*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3) -
		           42240*S*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,3) -
		           48896*U*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,3) +
		           96768*pow(l1k,4)*pow(Q2e,2)*pow(Q2h,3) +
		           128*l1k*S*U*pow(Q2e,3)*pow(Q2h,3) -
		           5440*S*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,3) -
		           5440*U*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,3) +
		           18304*pow(l1k,3)*pow(Q2e,3)*pow(Q2h,3) -
		           608*l1k*S*pow(Q2e,4)*pow(Q2h,3) -
		           144*l1k*U*pow(Q2e,4)*pow(Q2h,3) -
		           64*S*U*pow(Q2e,4)*pow(Q2h,3) +
		           2720*pow(l1k,2)*pow(Q2e,4)*pow(Q2h,3) +
		           376*l1k*pow(Q2e,5)*pow(Q2h,3) - 64*S*pow(Q2e,5)*pow(Q2h,3) +
		           72*U*pow(Q2e,5)*pow(Q2h,3) + 40*pow(Q2e,6)*pow(Q2h,3) -
		           5120*l1k*m4*Q2e*S*pow(Q2h,4) - 5376*l1k*m4*Q2e*U*pow(Q2h,4) -
		           3840*l1k*m4*S*U*pow(Q2h,4) + 5632*l1k*m2*Q2e*S*U*pow(Q2h,4) +
		           640*m4*Q2e*S*U*pow(Q2h,4) +
		           11008*m4*Q2e*pow(l1k,2)*pow(Q2h,4) +
		           10752*m4*S*pow(l1k,2)*pow(Q2h,4) -
		           32192*m2*Q2e*S*pow(l1k,2)*pow(Q2h,4) +
		           13312*m4*U*pow(l1k,2)*pow(Q2h,4) -
		           40128*m2*Q2e*U*pow(l1k,2)*pow(Q2h,4) -
		           20416*m2*S*U*pow(l1k,2)*pow(Q2h,4) -
		           8384*Q2e*S*U*pow(l1k,2)*pow(Q2h,4) -
		           6144*m4*pow(l1k,3)*pow(Q2h,4) +
		           65024*m2*Q2e*pow(l1k,3)*pow(Q2h,4) +
		           44032*m2*S*pow(l1k,3)*pow(Q2h,4) +
		           26752*Q2e*S*pow(l1k,3)*pow(Q2h,4) +
		           65536*m2*U*pow(l1k,3)*pow(Q2h,4) +
		           36992*Q2e*U*pow(l1k,3)*pow(Q2h,4) -
		           27648*m2*pow(l1k,4)*pow(Q2h,4) -
		           11008*Q2e*pow(l1k,4)*pow(Q2h,4) +
		           1024*l1k*Q2e*pow(m,6)*pow(Q2h,4) +
		           2048*l1k*S*pow(m,6)*pow(Q2h,4) -
		           256*Q2e*S*pow(m,6)*pow(Q2h,4) +
		           2048*l1k*U*pow(m,6)*pow(Q2h,4) -
		           256*Q2e*U*pow(m,6)*pow(Q2h,4) - 256*S*U*pow(m,6)*pow(Q2h,4) -
		           1024*pow(l1k,2)*pow(m,6)*pow(Q2h,4) -
		           5440*S*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,4) -
		           5440*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,4) +
		           27264*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,4) -
		           8*S*U*pow(Q2e,3)*pow(Q2h,4) +
		           2720*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,4) +
		           80*S*pow(Q2e,4)*pow(Q2h,4) + 8*U*pow(Q2e,4)*pow(Q2h,4) -
		           44*pow(Q2e,5)*pow(Q2h,4) + 2944*l1k*m4*Q2e*pow(Q2h,5) +
		           4352*l1k*m4*S*pow(Q2h,5) - 6144*l1k*m2*Q2e*S*pow(Q2h,5) -
		           640*m4*Q2e*S*pow(Q2h,5) + 4608*l1k*m4*U*pow(Q2h,5) -
		           7040*l1k*m2*Q2e*U*pow(Q2h,5) - 640*m4*Q2e*U*pow(Q2h,5) -
		           4992*l1k*m2*S*U*pow(Q2h,5) - 640*m4*S*U*pow(Q2h,5) -
		           800*l1k*Q2e*S*U*pow(Q2h,5) + 400*m2*Q2e*S*U*pow(Q2h,5) -
		           8192*m4*pow(l1k,2)*pow(Q2h,5) +
		           23040*m2*Q2e*pow(l1k,2)*pow(Q2h,5) +
		           24384*m2*S*pow(l1k,2)*pow(Q2h,5) +
		           9408*Q2e*S*pow(l1k,2)*pow(Q2h,5) +
		           31296*m2*U*pow(l1k,2)*pow(Q2h,5) +
		           11200*Q2e*U*pow(l1k,2)*pow(Q2h,5) -
		           41984*m2*pow(l1k,3)*pow(Q2h,5) -
		           22656*Q2e*pow(l1k,3)*pow(Q2h,5) -
		           1024*l1k*pow(m,6)*pow(Q2h,5) + 128*Q2e*pow(m,6)*pow(Q2h,5) +
		           256*S*pow(m,6)*pow(Q2h,5) + 256*U*pow(m,6)*pow(Q2h,5) +
		           2720*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,5) -
		           2560*l1k*m4*pow(Q2h,6) + 3776*l1k*m2*Q2e*pow(Q2h,6) +
		           320*m4*Q2e*pow(Q2h,6) + 5504*l1k*m2*S*pow(Q2h,6) +
		           640*m4*S*pow(Q2h,6) + 800*l1k*Q2e*S*pow(Q2h,6) -
		           400*m2*Q2e*S*pow(Q2h,6) + 6400*l1k*m2*U*pow(Q2h,6) +
		           640*m4*U*pow(Q2h,6) + 800*l1k*Q2e*U*pow(Q2h,6) -
		           400*m2*Q2e*U*pow(Q2h,6) - 400*m2*S*U*pow(Q2h,6) -
		           18240*m2*pow(l1k,2)*pow(Q2h,6) -
		           6112*Q2e*pow(l1k,2)*pow(Q2h,6) - 128*pow(m,6)*pow(Q2h,6) -
		           3456*l1k*m2*pow(Q2h,7) - 320*m4*pow(Q2h,7) -
		           400*l1k*Q2e*pow(Q2h,7) + 200*m2*Q2e*pow(Q2h,7) +
		           400*m2*S*pow(Q2h,7) + 400*m2*U*pow(Q2h,7) -
		           200*m2*pow(Q2h,8) + 4096*m4*Q2e*Q2h*pow(l1k,3)*pow(S,2) +
		           1024*m4*Q2e*pow(l1k,4)*pow(S,2) +
		           15872*m2*Q2e*Q2h*pow(l1k,4)*pow(S,2) +
		           4096*m2*Q2e*pow(l1k,5)*pow(S,2) +
		           30208*Q2e*Q2h*pow(l1k,5)*pow(S,2) +
		           14336*Q2e*pow(l1k,6)*pow(S,2) +
		           1024*Q2e*Q2h*pow(l1k,2)*pow(m,6)*pow(S,2) +
		           24192*Q2h*pow(l1k,4)*pow(Q2e,2)*pow(S,2) +
		           11264*pow(l1k,5)*pow(Q2e,2)*pow(S,2) +
		           4384*Q2h*pow(l1k,3)*pow(Q2e,3)*pow(S,2) +
		           2624*pow(l1k,4)*pow(Q2e,3)*pow(S,2) +
		           480*Q2h*pow(l1k,2)*pow(Q2e,4)*pow(S,2) +
		           320*pow(l1k,3)*pow(Q2e,4)*pow(S,2) +
		           24*l1k*Q2h*pow(Q2e,5)*pow(S,2) +
		           16*pow(l1k,2)*pow(Q2e,5)*pow(S,2) +
		           4992*m4*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(S,2) -
		           1024*m4*pow(l1k,3)*pow(Q2h,2)*pow(S,2) +
		           19456*m2*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(S,2) -
		           3584*m2*pow(l1k,4)*pow(Q2h,2)*pow(S,2) +
		           10112*Q2e*pow(l1k,4)*pow(Q2h,2)*pow(S,2) +
		           1024*l1k*Q2e*pow(m,6)*pow(Q2h,2)*pow(S,2) +
		           14976*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) +
		           2176*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2)*pow(S,2) +
		           224*l1k*pow(Q2e,4)*pow(Q2h,2)*pow(S,2) +
		           12*pow(Q2e,5)*pow(Q2h,2)*pow(S,2) +
		           2176*l1k*m4*Q2e*pow(Q2h,3)*pow(S,2) -
		           2688*m4*pow(l1k,2)*pow(Q2h,3)*pow(S,2) +
		           10208*m2*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(S,2) -
		           9344*m2*pow(l1k,3)*pow(Q2h,3)*pow(S,2) -
		           6208*Q2e*pow(l1k,3)*pow(Q2h,3)*pow(S,2) -
		           1024*l1k*pow(m,6)*pow(Q2h,3)*pow(S,2) +
		           128*Q2e*pow(m,6)*pow(Q2h,3)*pow(S,2) +
		           2720*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3)*pow(S,2) +
		           240*l1k*pow(Q2e,3)*pow(Q2h,3)*pow(S,2) +
		           24*pow(Q2e,4)*pow(Q2h,3)*pow(S,2) -
		           1792*l1k*m4*pow(Q2h,4)*pow(S,2) +
		           2368*l1k*m2*Q2e*pow(Q2h,4)*pow(S,2) +
		           320*m4*Q2e*pow(Q2h,4)*pow(S,2) -
		           7200*m2*pow(l1k,2)*pow(Q2h,4)*pow(S,2) -
		           3296*Q2e*pow(l1k,2)*pow(Q2h,4)*pow(S,2) -
		           128*pow(m,6)*pow(Q2h,4)*pow(S,2) -
		           36*pow(Q2e,3)*pow(Q2h,4)*pow(S,2) -
		           2048*l1k*m2*pow(Q2h,5)*pow(S,2) -
		           320*m4*pow(Q2h,5)*pow(S,2) - 400*l1k*Q2e*pow(Q2h,5)*pow(S,2) +
		           200*m2*Q2e*pow(Q2h,5)*pow(S,2) - 200*m2*pow(Q2h,6)*pow(S,2) +
		           7168*m4*Q2e*Q2h*pow(l1k,3)*pow(U,2) +
		           1024*m4*Q2e*pow(l1k,4)*pow(U,2) +
		           38400*m2*Q2e*Q2h*pow(l1k,4)*pow(U,2) +
		           8192*m2*Q2e*pow(l1k,5)*pow(U,2) +
		           79360*Q2e*Q2h*pow(l1k,5)*pow(U,2) +
		           43008*Q2e*pow(l1k,6)*pow(U,2) +
		           1024*Q2e*Q2h*pow(l1k,2)*pow(m,6)*pow(U,2) +
		           64128*Q2h*pow(l1k,4)*pow(Q2e,2)*pow(U,2) +
		           46080*pow(l1k,5)*pow(Q2e,2)*pow(U,2) +
		           13152*Q2h*pow(l1k,3)*pow(Q2e,3)*pow(U,2) +
		           20544*pow(l1k,4)*pow(Q2e,3)*pow(U,2) +
		           1888*Q2h*pow(l1k,2)*pow(Q2e,4)*pow(U,2) +
		           5760*pow(l1k,3)*pow(Q2e,4)*pow(U,2) +
		           256*l1k*Q2h*pow(Q2e,5)*pow(U,2) +
		           992*pow(l1k,2)*pow(Q2e,5)*pow(U,2) +
		           96*l1k*pow(Q2e,6)*pow(U,2) + 32*Q2h*pow(Q2e,6)*pow(U,2) +
		           4*pow(Q2e,7)*pow(U,2) +
		           7552*m4*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(U,2) -
		           2048*m4*pow(l1k,3)*pow(Q2h,2)*pow(U,2) +
		           41984*m2*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(U,2) -
		           10752*m2*pow(l1k,4)*pow(Q2h,2)*pow(U,2) +
		           2944*Q2e*pow(l1k,4)*pow(Q2h,2)*pow(U,2) +
		           1024*l1k*Q2e*pow(m,6)*pow(Q2h,2)*pow(U,2) +
		           21632*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) +
		           2608*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2)*pow(U,2) +
		           128*l1k*pow(Q2e,4)*pow(Q2h,2)*pow(U,2) -
		           28*pow(Q2e,5)*pow(Q2h,2)*pow(U,2) +
		           2432*l1k*m4*Q2e*pow(Q2h,3)*pow(U,2) -
		           5248*m4*pow(l1k,2)*pow(Q2h,3)*pow(U,2) +
		           17376*m2*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(U,2) -
		           25216*m2*pow(l1k,3)*pow(Q2h,3)*pow(U,2) -
		           14912*Q2e*pow(l1k,3)*pow(Q2h,3)*pow(U,2) -
		           1024*l1k*pow(m,6)*pow(Q2h,3)*pow(U,2) +
		           128*Q2e*pow(m,6)*pow(Q2h,3)*pow(U,2) +
		           2720*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3)*pow(U,2) +
		           8*l1k*pow(Q2e,3)*pow(Q2h,3)*pow(U,2) -
		           8*pow(Q2e,4)*pow(Q2h,3)*pow(U,2) -
		           2048*l1k*m4*pow(Q2h,4)*pow(U,2) +
		           3264*l1k*m2*Q2e*pow(Q2h,4)*pow(U,2) +
		           320*m4*Q2e*pow(Q2h,4)*pow(U,2) -
		           13344*m2*pow(l1k,2)*pow(Q2h,4)*pow(U,2) -
		           5088*Q2e*pow(l1k,2)*pow(Q2h,4)*pow(U,2) -
		           128*pow(m,6)*pow(Q2h,4)*pow(U,2) -
		           2944*l1k*m2*pow(Q2h,5)*pow(U,2) -
		           320*m4*pow(Q2h,5)*pow(U,2) - 400*l1k*Q2e*pow(Q2h,5)*pow(U,2) +
		           200*m2*Q2e*pow(Q2h,5)*pow(U,2) - 200*m2*pow(Q2h,6)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        (512*l1k*m4*M2*Q2e*Q2h + 864*l1k*m4*Q2e*Q2h*S +
		           1504*l1k*m4*Q2e*Q2h*U - 384*l1k*m4*Q2e*S*U -
		           1792*l1k*m4*Q2h*S*U - 4416*l1k*m2*Q2e*Q2h*S*U -
		           1056*m4*Q2e*Q2h*S*U + 512*m4*M2*Q2e*pow(l1k,2) +
		           2816*m2*M2*Q2e*Q2h*pow(l1k,2) + 1024*m4*Q2h*S*pow(l1k,2) +
		           4096*m2*Q2e*Q2h*S*pow(l1k,2) + 1024*m4*Q2h*U*pow(l1k,2) +
		           8128*m2*Q2e*Q2h*U*pow(l1k,2) - 256*m4*S*U*pow(l1k,2) -
		           2496*m2*Q2e*S*U*pow(l1k,2) - 7936*m2*Q2h*S*U*pow(l1k,2) -
		           10272*Q2e*Q2h*S*U*pow(l1k,2) + 2048*m2*M2*Q2e*pow(l1k,3) -
		           256*M2*Q2e*Q2h*pow(l1k,3) + 5120*m2*Q2h*S*pow(l1k,3) +
		           13792*Q2e*Q2h*S*pow(l1k,3) + 7168*m2*Q2h*U*pow(l1k,3) +
		           33632*Q2e*Q2h*U*pow(l1k,3) - 3072*m2*S*U*pow(l1k,3) -
		           9216*Q2e*S*U*pow(l1k,3) - 33024*Q2h*S*U*pow(l1k,3) +
		           7168*M2*Q2e*pow(l1k,4) + 20480*Q2h*S*pow(l1k,4) +
		           34816*Q2h*U*pow(l1k,4) - 13312*S*U*pow(l1k,4) +
		           192*Q2e*Q2h*S*pow(m,6) + 192*Q2e*Q2h*U*pow(m,6) -
		           128*Q2e*S*U*pow(m,6) - 128*Q2h*S*U*pow(m,6) -
		           3392*l1k*Q2h*S*U*pow(Q2e,2) +
		           1216*M2*Q2h*pow(l1k,2)*pow(Q2e,2) +
		           5312*Q2h*S*pow(l1k,2)*pow(Q2e,2) +
		           15984*Q2h*U*pow(l1k,2)*pow(Q2e,2) -
		           3424*S*U*pow(l1k,2)*pow(Q2e,2) +
		           5632*M2*pow(l1k,3)*pow(Q2e,2) + 240*l1k*M2*Q2h*pow(Q2e,3) +
		           528*l1k*Q2h*S*pow(Q2e,3) + 2528*l1k*Q2h*U*pow(Q2e,3) -
		           320*l1k*S*U*pow(Q2e,3) - 280*Q2h*S*U*pow(Q2e,3) +
		           1312*M2*pow(l1k,2)*pow(Q2e,3) + 160*l1k*M2*pow(Q2e,4) +
		           48*M2*Q2h*pow(Q2e,4) + 20*Q2h*S*pow(Q2e,4) +
		           204*Q2h*U*pow(Q2e,4) - 8*S*U*pow(Q2e,4) + 8*M2*pow(Q2e,5) -
		           512*l1k*m4*M2*pow(Q2h,2) - 1088*l1k*m4*Q2e*pow(Q2h,2) +
		           1024*l1k*m2*M2*Q2e*pow(Q2h,2) + 64*m4*M2*Q2e*pow(Q2h,2) +
		           2016*l1k*m4*S*pow(Q2h,2) + 6120*l1k*m2*Q2e*S*pow(Q2h,2) +
		           1408*m4*Q2e*S*pow(Q2h,2) + 2016*l1k*m4*U*pow(Q2h,2) +
		           10792*l1k*m2*Q2e*U*pow(Q2h,2) + 1856*m4*Q2e*U*pow(Q2h,2) -
		           2944*l1k*m2*S*U*pow(Q2h,2) - 160*m4*S*U*pow(Q2h,2) -
		           1968*l1k*Q2e*S*U*pow(Q2h,2) - 2136*m2*Q2e*S*U*pow(Q2h,2) -
		           1024*m4*pow(l1k,2)*pow(Q2h,2) -
		           1792*m2*M2*pow(l1k,2)*pow(Q2h,2) -
		           5824*m2*Q2e*pow(l1k,2)*pow(Q2h,2) -
		           1088*M2*Q2e*pow(l1k,2)*pow(Q2h,2) +
		           10752*m2*S*pow(l1k,2)*pow(Q2h,2) +
		           12512*Q2e*S*pow(l1k,2)*pow(Q2h,2) +
		           13312*m2*U*pow(l1k,2)*pow(Q2h,2) +
		           16800*Q2e*U*pow(l1k,2)*pow(Q2h,2) -
		           20416*S*U*pow(l1k,2)*pow(Q2h,2) -
		           6144*m2*pow(l1k,3)*pow(Q2h,2) -
		           23104*Q2e*pow(l1k,3)*pow(Q2h,2) +
		           44032*S*pow(l1k,3)*pow(Q2h,2) + 65536*U*pow(l1k,3)*pow(Q2h,2) -
		           27648*pow(l1k,4)*pow(Q2h,2) - 128*Q2e*pow(m,6)*pow(Q2h,2) +
		           64*S*pow(m,6)*pow(Q2h,2) + 64*U*pow(m,6)*pow(Q2h,2) +
		           4312*l1k*S*pow(Q2e,2)*pow(Q2h,2) +
		           6312*l1k*U*pow(Q2e,2)*pow(Q2h,2) -
		           624*S*U*pow(Q2e,2)*pow(Q2h,2) -
		           10352*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2) -
		           1552*l1k*pow(Q2e,3)*pow(Q2h,2) -
		           56*M2*pow(Q2e,3)*pow(Q2h,2) + 372*S*pow(Q2e,3)*pow(Q2h,2) +
		           568*U*pow(Q2e,3)*pow(Q2h,2) - 124*pow(Q2e,4)*pow(Q2h,2) -
		           1024*l1k*m4*pow(Q2h,3) - 832*l1k*m2*M2*pow(Q2h,3) -
		           64*m4*M2*pow(Q2h,3) - 6688*l1k*m2*Q2e*pow(Q2h,3) -
		           1104*m4*Q2e*pow(Q2h,3) - 288*l1k*M2*Q2e*pow(Q2h,3) +
		           144*m2*M2*Q2e*pow(Q2h,3) + 2824*l1k*m2*S*pow(Q2h,3) -
		           96*m4*S*pow(Q2h,3) + 1944*l1k*Q2e*S*pow(Q2h,3) +
		           2668*m2*Q2e*S*pow(Q2h,3) + 680*l1k*m2*U*pow(Q2h,3) -
		           480*m4*U*pow(Q2h,3) + 1240*l1k*Q2e*U*pow(Q2h,3) +
		           3836*m2*Q2e*U*pow(Q2h,3) - 4992*l1k*S*U*pow(Q2h,3) +
		           312*m2*S*U*pow(Q2h,3) - 8192*m2*pow(l1k,2)*pow(Q2h,3) -
		           9712*Q2e*pow(l1k,2)*pow(Q2h,3) +
		           24384*S*pow(l1k,2)*pow(Q2h,3) + 31296*U*pow(l1k,2)*pow(Q2h,3) -
		           41984*pow(l1k,3)*pow(Q2h,3) - 3728*l1k*pow(Q2e,2)*pow(Q2h,3) +
		           660*S*pow(Q2e,2)*pow(Q2h,3) + 656*U*pow(Q2e,2)*pow(Q2h,3) -
		           332*pow(Q2e,3)*pow(Q2h,3) + 16*l1k*m2*pow(Q2h,4) +
		           368*m4*pow(Q2h,4) - 144*m2*M2*pow(Q2h,4) -
		           520*l1k*Q2e*pow(Q2h,4) - 2228*m2*Q2e*pow(Q2h,4) +
		           5504*l1k*S*pow(Q2h,4) - 644*m2*S*pow(Q2h,4) -
		           64*Q2e*S*pow(Q2h,4) + 6400*l1k*U*pow(Q2h,4) -
		           1572*m2*U*pow(Q2h,4) - 192*Q2e*U*pow(Q2h,4) -
		           400*S*U*pow(Q2h,4) - 18240*pow(l1k,2)*pow(Q2h,4) -
		           328*pow(Q2e,2)*pow(Q2h,4) - 3456*l1k*pow(Q2h,5) +
		           996*m2*pow(Q2h,5) + 128*Q2e*pow(Q2h,5) + 400*S*pow(Q2h,5) +
		           400*U*pow(Q2h,5) - 200*pow(Q2h,6) - 832*l1k*m4*Q2h*pow(S,2) -
		           912*l1k*m2*Q2e*Q2h*pow(S,2) - 384*m4*Q2e*Q2h*pow(S,2) -
		           128*m4*pow(l1k,2)*pow(S,2) - 480*m2*Q2e*pow(l1k,2)*pow(S,2) -
		           2688*m2*Q2h*pow(l1k,2)*pow(S,2) -
		           2432*Q2e*Q2h*pow(l1k,2)*pow(S,2) -
		           1024*m2*pow(l1k,3)*pow(S,2) - 1344*Q2e*pow(l1k,3)*pow(S,2) -
		           9344*Q2h*pow(l1k,3)*pow(S,2) - 3584*pow(l1k,4)*pow(S,2) -
		           64*Q2e*pow(m,6)*pow(S,2) - 64*Q2h*pow(m,6)*pow(S,2) -
		           600*l1k*Q2h*pow(Q2e,2)*pow(S,2) -
		           304*pow(l1k,2)*pow(Q2e,2)*pow(S,2) -
		           8*l1k*pow(Q2e,3)*pow(S,2) - 24*Q2h*pow(Q2e,3)*pow(S,2) -
		           1712*l1k*m2*pow(Q2h,2)*pow(S,2) -
		           192*m4*pow(Q2h,2)*pow(S,2) -
		           1008*l1k*Q2e*pow(Q2h,2)*pow(S,2) -
		           648*m2*Q2e*pow(Q2h,2)*pow(S,2) -
		           7200*pow(l1k,2)*pow(Q2h,2)*pow(S,2) -
		           244*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) -
		           2048*l1k*pow(Q2h,3)*pow(S,2) - 144*m2*pow(Q2h,3)*pow(S,2) -
		           64*Q2e*pow(Q2h,3)*pow(S,2) - 200*pow(Q2h,4)*pow(S,2) -
		           384*l1k*m4*Q2e*pow(U,2) - 960*l1k*m4*Q2h*pow(U,2) -
		           4208*l1k*m2*Q2e*Q2h*pow(U,2) - 736*m4*Q2e*Q2h*pow(U,2) -
		           128*m4*pow(l1k,2)*pow(U,2) -
		           2656*m2*Q2e*pow(l1k,2)*pow(U,2) -
		           5248*m2*Q2h*pow(l1k,2)*pow(U,2) -
		           7328*Q2e*Q2h*pow(l1k,2)*pow(U,2) -
		           2048*m2*pow(l1k,3)*pow(U,2) - 11712*Q2e*pow(l1k,3)*pow(U,2) -
		           25216*Q2h*pow(l1k,3)*pow(U,2) - 10752*pow(l1k,4)*pow(U,2) -
		           64*Q2e*pow(m,6)*pow(U,2) - 64*Q2h*pow(m,6)*pow(U,2) -
		           2616*l1k*Q2h*pow(Q2e,2)*pow(U,2) -
		           5904*pow(l1k,2)*pow(Q2e,2)*pow(U,2) -
		           1016*l1k*pow(Q2e,3)*pow(U,2) - 232*Q2h*pow(Q2e,3)*pow(U,2) -
		           84*pow(Q2e,4)*pow(U,2) - 656*l1k*m2*pow(Q2h,2)*pow(U,2) +
		           96*m4*pow(Q2h,2)*pow(U,2) - 704*l1k*Q2e*pow(Q2h,2)*pow(U,2) -
		           1616*m2*Q2e*pow(Q2h,2)*pow(U,2) -
		           13344*pow(l1k,2)*pow(Q2h,2)*pow(U,2) -
		           328*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) -
		           2944*l1k*pow(Q2h,3)*pow(U,2) + 584*m2*pow(Q2h,3)*pow(U,2) +
		           64*Q2e*pow(Q2h,3)*pow(U,2) - 200*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*(-64*Q2e*Q2h*S*U*pow(m,6) - 8*Q2h*S*U*pow(Q2e,4) +
		           4*M2*Q2h*pow(Q2e,5) + 6*Q2h*U*pow(Q2e,5) -
		           304*m4*Q2e*S*U*pow(Q2h,2) + 32*Q2e*S*pow(m,6)*pow(Q2h,2) +
		           96*Q2e*U*pow(m,6)*pow(Q2h,2) + 64*S*U*pow(m,6)*pow(Q2h,2) -
		           32*S*U*pow(Q2e,3)*pow(Q2h,2) - 8*M2*pow(Q2e,4)*pow(Q2h,2) +
		           14*S*pow(Q2e,4)*pow(Q2h,2) + 30*U*pow(Q2e,4)*pow(Q2h,2) -
		           4*pow(Q2e,5)*pow(Q2h,2) + 328*m4*Q2e*S*pow(Q2h,3) +
		           408*m4*Q2e*U*pow(Q2h,3) + 176*m4*S*U*pow(Q2h,3) -
		           328*m2*Q2e*S*U*pow(Q2h,3) - 32*Q2e*pow(m,6)*pow(Q2h,3) -
		           32*S*pow(m,6)*pow(Q2h,3) - 96*U*pow(m,6)*pow(Q2h,3) +
		           4*M2*pow(Q2e,3)*pow(Q2h,3) + 26*S*pow(Q2e,3)*pow(Q2h,3) +
		           4*U*pow(Q2e,3)*pow(Q2h,3) - 16*pow(Q2e,4)*pow(Q2h,3) -
		           216*m4*Q2e*pow(Q2h,4) - 200*m4*S*pow(Q2h,4) +
		           360*m2*Q2e*S*pow(Q2h,4) - 280*m4*U*pow(Q2h,4) +
		           424*m2*Q2e*U*pow(Q2h,4) + 200*m2*S*U*pow(Q2h,4) +
		           32*pow(m,6)*pow(Q2h,4) + 152*m4*pow(Q2h,5) -
		           228*m2*Q2e*pow(Q2h,5) - 232*m2*S*pow(Q2h,5) -
		           296*m2*U*pow(Q2h,5) + 164*m2*pow(Q2h,6) -
		           112*m4*Q2e*pow(Q2h,2)*pow(S,2) -
		           10*pow(Q2e,3)*pow(Q2h,2)*pow(S,2) + 48*m4*pow(Q2h,3)*pow(S,2) -
		           132*m2*Q2e*pow(Q2h,3)*pow(S,2) -
		           10*pow(Q2e,2)*pow(Q2h,3)*pow(S,2) + 68*m2*pow(Q2h,4)*pow(S,2) -
		           64*Q2e*Q2h*pow(m,6)*pow(U,2) - 14*Q2h*pow(Q2e,4)*pow(U,2) -
		           2*pow(Q2e,5)*pow(U,2) - 192*m4*Q2e*pow(Q2h,2)*pow(U,2) +
		           64*pow(m,6)*pow(Q2h,2)*pow(U,2) -
		           4*pow(Q2e,3)*pow(Q2h,2)*pow(U,2) + 128*m4*pow(Q2h,3)*pow(U,2) -
		           196*m2*Q2e*pow(Q2h,3)*pow(U,2) + 132*m2*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1)) +
		     SI.C0_mmqm0m(Q2e,m2)*
		      (256*l1k*M2 + 64*m2*M2 + 48*M2*Q2e + 192*M2*Q2h + 184*Q2h*S +
		        360*Q2h*U - 48*S*U - 288*pow(Q2h,2) - 24*pow(S,2) - 120*pow(U,2) +
		        pow(l1k,-1)*(32*m2*M2*Q2e + 48*m2*M2*Q2h + 28*M2*Q2e*Q2h +
		           32*m2*Q2h*S + 8*Q2e*Q2h*S + 72*m2*Q2h*U + 42*Q2e*Q2h*U -
		           60*Q2h*S*U + 8*M2*pow(Q2e,2) - 56*m2*pow(Q2h,2) +
		           24*M2*pow(Q2h,2) - 30*Q2e*pow(Q2h,2) + 106*S*pow(Q2h,2) +
		           194*U*pow(Q2h,2) - 126*pow(Q2h,3) - 8*Q2h*pow(S,2) -
		           16*m2*pow(U,2) - 12*Q2e*pow(U,2) - 68*Q2h*pow(U,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-64*m2*M2*Q2e - 96*m2*M2*Q2h - 56*M2*Q2e*Q2h -
		           144*m2*Q2h*S - 84*Q2e*Q2h*S - 64*m2*Q2h*U - 16*Q2e*Q2h*U +
		           120*Q2h*S*U - 16*M2*pow(Q2e,2) + 112*m2*pow(Q2h,2) -
		           48*M2*pow(Q2h,2) + 60*Q2e*pow(Q2h,2) - 388*S*pow(Q2h,2) -
		           212*U*pow(Q2h,2) + 252*pow(Q2h,3) + 32*m2*pow(S,2) +
		           24*Q2e*pow(S,2) + 136*Q2h*pow(S,2) + 16*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-112*m2*Q2h*S*U + 32*m2*M2*pow(Q2h,2) +
		           264*m2*S*pow(Q2h,2) + 200*m2*U*pow(Q2h,2) -
		           252*S*U*pow(Q2h,2) - 184*m2*pow(Q2h,3) + 24*M2*pow(Q2h,3) +
		           478*S*pow(Q2h,3) + 334*U*pow(Q2h,3) - 286*pow(Q2h,4) -
		           80*m2*Q2h*pow(S,2) - 192*pow(Q2h,2)*pow(S,2) -
		           32*m2*Q2h*pow(U,2) - 76*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (128*m2*S*U*pow(Q2h,2) - 168*m2*S*pow(Q2h,3) -
		           152*m2*U*pow(Q2h,3) + 160*S*U*pow(Q2h,3) + 96*m2*pow(Q2h,4) -
		           206*S*pow(Q2h,4) - 178*U*pow(Q2h,4) + 112*pow(Q2h,5) +
		           72*m2*pow(Q2h,2)*pow(S,2) + 94*pow(Q2h,3)*pow(S,2) +
		           56*m2*pow(Q2h,2)*pow(U,2) + 66*pow(Q2h,3)*pow(U,2)) +
		        pow(l1k,-3)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*m2*S*U*pow(Q2h,3) + 32*m2*S*pow(Q2h,4) +
		           32*m2*U*pow(Q2h,4) - 32*S*U*pow(Q2h,4) - 16*m2*pow(Q2h,5) +
		           32*S*pow(Q2h,5) + 32*U*pow(Q2h,5) - 16*pow(Q2h,6) -
		           16*m2*pow(Q2h,3)*pow(S,2) - 16*pow(Q2h,4)*pow(S,2) -
		           16*m2*pow(Q2h,3)*pow(U,2) - 16*pow(Q2h,4)*pow(U,2)) +
		        pow(l1k,-1)*(-4*S*U*pow(Q2e,5)*pow(Q2h,2) +
		           4*U*pow(Q2e,6)*pow(Q2h,2) + 4*S*U*pow(Q2e,4)*pow(Q2h,3) +
		           4*S*pow(Q2e,5)*pow(Q2h,3) - 4*U*pow(Q2e,5)*pow(Q2h,3) -
		           2*pow(Q2e,6)*pow(Q2h,3) - 4*S*pow(Q2e,4)*pow(Q2h,4) +
		           2*pow(Q2e,5)*pow(Q2h,4) - 2*pow(Q2e,4)*pow(Q2h,3)*pow(S,2) +
		           2*pow(Q2e,3)*pow(Q2h,4)*pow(S,2) - 2*Q2h*pow(Q2e,6)*pow(U,2) +
		           2*pow(Q2e,5)*pow(Q2h,2)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (8*S*U*pow(Q2e,5)*pow(Q2h,2) - 8*S*pow(Q2e,6)*pow(Q2h,2) -
		           8*S*U*pow(Q2e,4)*pow(Q2h,3) + 8*S*pow(Q2e,5)*pow(Q2h,3) -
		           8*U*pow(Q2e,5)*pow(Q2h,3) + 4*pow(Q2e,6)*pow(Q2h,3) +
		           8*U*pow(Q2e,4)*pow(Q2h,4) - 4*pow(Q2e,5)*pow(Q2h,4) +
		           4*Q2h*pow(Q2e,6)*pow(S,2) - 4*pow(Q2e,5)*pow(Q2h,2)*pow(S,2) +
		           4*pow(Q2e,4)*pow(Q2h,3)*pow(U,2) -
		           4*pow(Q2e,3)*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        (-2560*m4*Q2e*Q2h*S*U*pow(l1k,2) -
		           7168*m2*Q2e*Q2h*S*U*pow(l1k,3) +
		           4096*m2*Q2e*Q2h*S*pow(l1k,4) + 4096*m2*Q2e*Q2h*U*pow(l1k,4) -
		           2048*m2*Q2e*S*U*pow(l1k,4) - 15360*Q2e*Q2h*S*U*pow(l1k,4) +
		           12288*Q2e*Q2h*S*pow(l1k,5) + 20480*Q2e*Q2h*U*pow(l1k,5) -
		           8192*Q2e*S*U*pow(l1k,5) - 12288*Q2h*S*U*pow(l1k,3)*pow(Q2e,2) +
		           10240*Q2h*S*pow(l1k,4)*pow(Q2e,2) +
		           20480*Q2h*U*pow(l1k,4)*pow(Q2e,2) -
		           7168*S*U*pow(l1k,4)*pow(Q2e,2) -
		           2400*Q2h*S*U*pow(l1k,2)*pow(Q2e,3) +
		           2560*Q2h*S*pow(l1k,3)*pow(Q2e,3) +
		           7680*Q2h*U*pow(l1k,3)*pow(Q2e,3) -
		           2048*S*U*pow(l1k,3)*pow(Q2e,3) - 288*l1k*Q2h*S*U*pow(Q2e,4) +
		           320*Q2h*S*pow(l1k,2)*pow(Q2e,4) +
		           1600*Q2h*U*pow(l1k,2)*pow(Q2e,4) -
		           288*S*U*pow(l1k,2)*pow(Q2e,4) + 16*l1k*Q2h*S*pow(Q2e,5) +
		           176*l1k*Q2h*U*pow(Q2e,5) - 16*l1k*S*U*pow(Q2e,5) -
		           16*Q2h*S*U*pow(Q2e,5) + 8*Q2h*U*pow(Q2e,6) -
		           1024*l1k*m4*Q2e*S*U*pow(Q2h,2) +
		           4096*m4*Q2e*S*pow(l1k,2)*pow(Q2h,2) +
		           4096*m4*Q2e*U*pow(l1k,2)*pow(Q2h,2) +
		           512*m4*S*U*pow(l1k,2)*pow(Q2h,2) -
		           10496*m2*Q2e*S*U*pow(l1k,2)*pow(Q2h,2) +
		           9216*m2*Q2e*S*pow(l1k,3)*pow(Q2h,2) +
		           15360*m2*Q2e*U*pow(l1k,3)*pow(Q2h,2) +
		           2048*m2*S*U*pow(l1k,3)*pow(Q2h,2) -
		           1536*Q2e*S*U*pow(l1k,3)*pow(Q2h,2) -
		           4096*m2*Q2e*pow(l1k,4)*pow(Q2h,2) +
		           20480*Q2e*S*pow(l1k,4)*pow(Q2h,2) +
		           28672*Q2e*U*pow(l1k,4)*pow(Q2h,2) -
		           16384*Q2e*pow(l1k,5)*pow(Q2h,2) -
		           512*Q2e*S*U*pow(m,6)*pow(Q2h,2) -
		           5504*S*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2) +
		           15872*S*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) +
		           23040*U*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) -
		           15360*pow(l1k,4)*pow(Q2e,2)*pow(Q2h,2) -
		           656*l1k*S*U*pow(Q2e,3)*pow(Q2h,2) +
		           2944*S*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2) +
		           4480*U*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2) -
		           5120*pow(l1k,3)*pow(Q2e,3)*pow(Q2h,2) +
		           304*l1k*S*pow(Q2e,4)*pow(Q2h,2) +
		           656*l1k*U*pow(Q2e,4)*pow(Q2h,2) -
		           72*S*U*pow(Q2e,4)*pow(Q2h,2) -
		           960*pow(l1k,2)*pow(Q2e,4)*pow(Q2h,2) -
		           96*l1k*pow(Q2e,5)*pow(Q2h,2) + 24*S*pow(Q2e,5)*pow(Q2h,2) +
		           72*U*pow(Q2e,5)*pow(Q2h,2) - 4*pow(Q2e,6)*pow(Q2h,2) +
		           1280*l1k*m4*Q2e*S*pow(Q2h,3) + 1792*l1k*m4*Q2e*U*pow(Q2h,3) +
		           1024*l1k*m4*S*U*pow(Q2h,3) - 2816*l1k*m2*Q2e*S*U*pow(Q2h,3) -
		           1408*m4*Q2e*S*U*pow(Q2h,3) -
		           3072*m4*Q2e*pow(l1k,2)*pow(Q2h,3) -
		           1024*m4*S*pow(l1k,2)*pow(Q2h,3) +
		           14336*m2*Q2e*S*pow(l1k,2)*pow(Q2h,3) -
		           1024*m4*U*pow(l1k,2)*pow(Q2h,3) +
		           17920*m2*Q2e*U*pow(l1k,2)*pow(Q2h,3) +
		           4864*m2*S*U*pow(l1k,2)*pow(Q2h,3) +
		           2176*Q2e*S*U*pow(l1k,2)*pow(Q2h,3) -
		           10240*m2*Q2e*pow(l1k,3)*pow(Q2h,3) -
		           3072*m2*S*pow(l1k,3)*pow(Q2h,3) -
		           512*Q2e*S*pow(l1k,3)*pow(Q2h,3) -
		           5120*m2*U*pow(l1k,3)*pow(Q2h,3) -
		           512*Q2e*U*pow(l1k,3)*pow(Q2h,3) -
		           17408*Q2e*pow(l1k,4)*pow(Q2h,3) +
		           512*Q2e*S*pow(m,6)*pow(Q2h,3) +
		           512*Q2e*U*pow(m,6)*pow(Q2h,3) + 512*S*U*pow(m,6)*pow(Q2h,3) -
		           768*l1k*S*U*pow(Q2e,2)*pow(Q2h,3) +
		           6528*S*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3) +
		           8064*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3) -
		           14336*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,3) +
		           768*l1k*S*pow(Q2e,3)*pow(Q2h,3) +
		           768*l1k*U*pow(Q2e,3)*pow(Q2h,3) - 8*S*U*pow(Q2e,3)*pow(Q2h,3) -
		           2752*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,3) -
		           384*l1k*pow(Q2e,4)*pow(Q2h,3) + 72*S*pow(Q2e,4)*pow(Q2h,3) +
		           16*U*pow(Q2e,4)*pow(Q2h,3) - 44*pow(Q2e,5)*pow(Q2h,3) -
		           1024*l1k*m4*Q2e*pow(Q2h,4) - 1280*l1k*m4*S*pow(Q2h,4) +
		           2944*l1k*m2*Q2e*S*pow(Q2h,4) + 1664*m4*Q2e*S*pow(Q2h,4) -
		           1792*l1k*m4*U*pow(Q2h,4) + 4992*l1k*m2*Q2e*U*pow(Q2h,4) +
		           1664*m4*Q2e*U*pow(Q2h,4) + 2816*l1k*m2*S*U*pow(Q2h,4) +
		           896*m4*S*U*pow(Q2h,4) + 1728*l1k*Q2e*S*U*pow(Q2h,4) -
		           1248*m2*Q2e*S*U*pow(Q2h,4) + 1024*m4*pow(l1k,2)*pow(Q2h,4) -
		           11776*m2*Q2e*pow(l1k,2)*pow(Q2h,4) -
		           6656*m2*S*pow(l1k,2)*pow(Q2h,4) -
		           2432*Q2e*S*pow(l1k,2)*pow(Q2h,4) -
		           9728*m2*U*pow(l1k,2)*pow(Q2h,4) -
		           4992*Q2e*U*pow(l1k,2)*pow(Q2h,4) +
		           4096*m2*pow(l1k,3)*pow(Q2h,4) +
		           2048*Q2e*pow(l1k,3)*pow(Q2h,4) - 256*Q2e*pow(m,6)*pow(Q2h,4) -
		           512*S*pow(m,6)*pow(Q2h,4) - 512*U*pow(m,6)*pow(Q2h,4) +
		           768*l1k*S*pow(Q2e,2)*pow(Q2h,4) +
		           768*l1k*U*pow(Q2e,2)*pow(Q2h,4) -
		           4544*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,4) -
		           384*l1k*pow(Q2e,3)*pow(Q2h,4) + 1024*l1k*m4*pow(Q2h,5) -
		           2688*l1k*m2*Q2e*pow(Q2h,5) - 960*m4*Q2e*pow(Q2h,5) -
		           3200*l1k*m2*S*pow(Q2h,5) - 1152*m4*S*pow(Q2h,5) -
		           1984*l1k*Q2e*S*pow(Q2h,5) + 1376*m2*Q2e*S*pow(Q2h,5) -
		           4736*l1k*m2*U*pow(Q2h,5) - 1152*m4*U*pow(Q2h,5) -
		           2240*l1k*Q2e*U*pow(Q2h,5) + 1504*m2*Q2e*U*pow(Q2h,5) +
		           864*m2*S*U*pow(Q2h,5) + 6400*m2*pow(l1k,2)*pow(Q2h,5) +
		           2880*Q2e*pow(l1k,2)*pow(Q2h,5) + 256*pow(m,6)*pow(Q2h,5) -
		           384*l1k*pow(Q2e,2)*pow(Q2h,5) + 2688*l1k*m2*pow(Q2h,6) +
		           704*m4*pow(Q2h,6) + 1248*l1k*Q2e*pow(Q2h,6) -
		           816*m2*Q2e*pow(Q2h,6) - 992*m2*S*pow(Q2h,6) -
		           1120*m2*U*pow(Q2h,6) + 624*m2*pow(Q2h,7) -
		           1280*m4*Q2e*Q2h*pow(l1k,2)*pow(S,2) -
		           1536*m2*Q2e*Q2h*pow(l1k,3)*pow(S,2) -
		           1024*m2*Q2e*pow(l1k,4)*pow(S,2) -
		           4608*Q2e*Q2h*pow(l1k,4)*pow(S,2) -
		           2048*Q2e*pow(l1k,5)*pow(S,2) -
		           3072*Q2h*pow(l1k,3)*pow(Q2e,2)*pow(S,2) -
		           1536*pow(l1k,4)*pow(Q2e,2)*pow(S,2) -
		           432*Q2h*pow(l1k,2)*pow(Q2e,3)*pow(S,2) -
		           256*pow(l1k,3)*pow(Q2e,3)*pow(S,2) -
		           16*l1k*Q2h*pow(Q2e,4)*pow(S,2) -
		           16*pow(l1k,2)*pow(Q2e,4)*pow(S,2) - 4*Q2h*pow(Q2e,5)*pow(S,2) -
		           256*l1k*m4*Q2e*pow(Q2h,2)*pow(S,2) +
		           256*m4*pow(l1k,2)*pow(Q2h,2)*pow(S,2) -
		           3968*m2*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(S,2) +
		           512*m2*pow(l1k,3)*pow(Q2h,2)*pow(S,2) -
		           256*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(S,2) -
		           256*Q2e*pow(m,6)*pow(Q2h,2)*pow(S,2) -
		           1984*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) -
		           208*l1k*pow(Q2e,3)*pow(Q2h,2)*pow(S,2) -
		           12*pow(Q2e,4)*pow(Q2h,2)*pow(S,2) +
		           256*l1k*m4*pow(Q2h,3)*pow(S,2) -
		           512*l1k*m2*Q2e*pow(Q2h,3)*pow(S,2) -
		           704*m4*Q2e*pow(Q2h,3)*pow(S,2) +
		           1408*m2*pow(l1k,2)*pow(Q2h,3)*pow(S,2) +
		           64*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(S,2) +
		           256*pow(m,6)*pow(Q2h,3)*pow(S,2) -
		           384*l1k*pow(Q2e,2)*pow(Q2h,3)*pow(S,2) -
		           32*pow(Q2e,3)*pow(Q2h,3)*pow(S,2) +
		           768*l1k*m2*pow(Q2h,4)*pow(S,2) + 448*m4*pow(Q2h,4)*pow(S,2) +
		           736*l1k*Q2e*pow(Q2h,4)*pow(S,2) -
		           560*m2*Q2e*pow(Q2h,4)*pow(S,2) + 368*m2*pow(Q2h,5)*pow(S,2) -
		           1280*m4*Q2e*Q2h*pow(l1k,2)*pow(U,2) -
		           5632*m2*Q2e*Q2h*pow(l1k,3)*pow(U,2) -
		           1024*m2*Q2e*pow(l1k,4)*pow(U,2) -
		           11776*Q2e*Q2h*pow(l1k,4)*pow(U,2) -
		           6144*Q2e*pow(l1k,5)*pow(U,2) -
		           9216*Q2h*pow(l1k,3)*pow(Q2e,2)*pow(U,2) -
		           6656*pow(l1k,4)*pow(Q2e,2)*pow(U,2) -
		           1840*Q2h*pow(l1k,2)*pow(Q2e,3)*pow(U,2) -
		           2816*pow(l1k,3)*pow(Q2e,3)*pow(U,2) -
		           272*l1k*Q2h*pow(Q2e,4)*pow(U,2) -
		           656*pow(l1k,2)*pow(Q2e,4)*pow(U,2) -
		           80*l1k*pow(Q2e,5)*pow(U,2) - 32*Q2h*pow(Q2e,5)*pow(U,2) -
		           4*pow(Q2e,6)*pow(U,2) - 768*l1k*m4*Q2e*pow(Q2h,2)*pow(U,2) +
		           256*m4*pow(l1k,2)*pow(Q2h,2)*pow(U,2) -
		           6784*m2*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(U,2) +
		           1536*m2*pow(l1k,3)*pow(Q2h,2)*pow(U,2) -
		           768*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(U,2) -
		           256*Q2e*pow(m,6)*pow(Q2h,2)*pow(U,2) -
		           3520*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) -
		           384*l1k*pow(Q2e,3)*pow(Q2h,2)*pow(U,2) -
		           8*pow(Q2e,4)*pow(Q2h,2)*pow(U,2) +
		           768*l1k*m4*pow(Q2h,3)*pow(U,2) -
		           2304*l1k*m2*Q2e*pow(Q2h,3)*pow(U,2) -
		           704*m4*Q2e*pow(Q2h,3)*pow(U,2) +
		           3712*m2*pow(l1k,2)*pow(Q2h,3)*pow(U,2) +
		           2112*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(U,2) +
		           256*pow(m,6)*pow(Q2h,3)*pow(U,2) -
		           384*l1k*pow(Q2e,2)*pow(Q2h,3)*pow(U,2) -
		           4*pow(Q2e,3)*pow(Q2h,3)*pow(U,2) +
		           2048*l1k*m2*pow(Q2h,4)*pow(U,2) +
		           448*m4*pow(Q2h,4)*pow(U,2) + 992*l1k*Q2e*pow(Q2h,4)*pow(U,2) -
		           688*m2*Q2e*pow(Q2h,4)*pow(U,2) + 496*m2*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        (-256*l1k*m2*M2*Q2e*Q2h - 128*m4*M2*Q2e*Q2h -
		           160*l1k*m2*Q2e*Q2h*S - 224*m4*Q2e*Q2h*S -
		           864*l1k*m2*Q2e*Q2h*U - 224*m4*Q2e*Q2h*U +
		           256*l1k*m2*Q2e*S*U + 64*m4*Q2e*S*U + 1024*l1k*m2*Q2h*S*U +
		           448*m4*Q2h*S*U + 928*l1k*Q2e*Q2h*S*U + 720*m2*Q2e*Q2h*S*U -
		           512*m2*M2*Q2e*pow(l1k,2) - 256*M2*Q2e*Q2h*pow(l1k,2) -
		           1024*m2*Q2h*S*pow(l1k,2) - 1760*Q2e*Q2h*S*pow(l1k,2) -
		           1024*m2*Q2h*U*pow(l1k,2) - 4000*Q2e*Q2h*U*pow(l1k,2) +
		           512*m2*S*U*pow(l1k,2) + 960*Q2e*S*U*pow(l1k,2) +
		           4864*Q2h*S*U*pow(l1k,2) - 1024*M2*Q2e*pow(l1k,3) -
		           3072*Q2h*S*pow(l1k,3) - 5120*Q2h*U*pow(l1k,3) +
		           2048*S*U*pow(l1k,3) - 256*l1k*M2*Q2h*pow(Q2e,2) -
		           376*l1k*Q2h*S*pow(Q2e,2) - 1480*l1k*Q2h*U*pow(Q2e,2) +
		           224*l1k*S*U*pow(Q2e,2) + 256*Q2h*S*U*pow(Q2e,2) -
		           768*M2*pow(l1k,2)*pow(Q2e,2) - 128*l1k*M2*pow(Q2e,3) -
		           56*M2*Q2h*pow(Q2e,3) - 32*Q2h*S*pow(Q2e,3) -
		           156*Q2h*U*pow(Q2e,3) + 8*S*U*pow(Q2e,3) - 8*M2*pow(Q2e,4) +
		           256*l1k*m2*M2*pow(Q2h,2) + 128*m4*M2*pow(Q2h,2) +
		           512*l1k*m2*Q2e*pow(Q2h,2) + 192*m4*Q2e*pow(Q2h,2) +
		           384*l1k*M2*Q2e*pow(Q2h,2) - 320*m2*M2*Q2e*pow(Q2h,2) -
		           1280*l1k*m2*S*pow(Q2h,2) - 480*m4*S*pow(Q2h,2) -
		           432*l1k*Q2e*S*pow(Q2h,2) - 1464*m2*Q2e*S*pow(Q2h,2) -
		           1792*l1k*m2*U*pow(Q2h,2) - 480*m4*U*pow(Q2h,2) -
		           1808*l1k*Q2e*U*pow(Q2h,2) - 1640*m2*Q2e*U*pow(Q2h,2) +
		           2816*l1k*S*U*pow(Q2h,2) + 848*m2*S*U*pow(Q2h,2) +
		           864*Q2e*S*U*pow(Q2h,2) + 1024*m2*pow(l1k,2)*pow(Q2h,2) +
		           2944*Q2e*pow(l1k,2)*pow(Q2h,2) - 6656*S*pow(l1k,2)*pow(Q2h,2) -
		           9728*U*pow(l1k,2)*pow(Q2h,2) + 4096*pow(l1k,3)*pow(Q2h,2) +
		           960*l1k*pow(Q2e,2)*pow(Q2h,2) - 420*S*pow(Q2e,2)*pow(Q2h,2) -
		           664*U*pow(Q2e,2)*pow(Q2h,2) + 108*pow(Q2e,3)*pow(Q2h,2) +
		           1024*l1k*m2*pow(Q2h,3) + 256*m4*pow(Q2h,3) +
		           192*m2*M2*pow(Q2h,3) + 704*l1k*Q2e*pow(Q2h,3) +
		           1280*m2*Q2e*pow(Q2h,3) - 3200*l1k*S*pow(Q2h,3) -
		           968*m2*S*pow(Q2h,3) - 1352*Q2e*S*pow(Q2h,3) -
		           4736*l1k*U*pow(Q2h,3) - 792*m2*U*pow(Q2h,3) -
		           1016*Q2e*U*pow(Q2h,3) + 864*S*U*pow(Q2h,3) +
		           6400*pow(l1k,2)*pow(Q2h,3) + 432*pow(Q2e,2)*pow(Q2h,3) +
		           2688*l1k*pow(Q2h,4) + 416*m2*pow(Q2h,4) + 752*Q2e*pow(Q2h,4) -
		           992*S*pow(Q2h,4) - 1120*U*pow(Q2h,4) + 624*pow(Q2h,5) -
		           64*l1k*m2*Q2e*pow(S,2) + 32*m4*Q2e*pow(S,2) +
		           256*l1k*m2*Q2h*pow(S,2) + 224*m4*Q2h*pow(S,2) -
		           160*l1k*Q2e*Q2h*pow(S,2) + 376*m2*Q2e*Q2h*pow(S,2) +
		           256*m2*pow(l1k,2)*pow(S,2) + 224*Q2e*pow(l1k,2)*pow(S,2) +
		           1408*Q2h*pow(l1k,2)*pow(S,2) + 512*pow(l1k,3)*pow(S,2) +
		           60*Q2h*pow(Q2e,2)*pow(S,2) + 4*pow(Q2e,3)*pow(S,2) +
		           768*l1k*pow(Q2h,2)*pow(S,2) + 424*m2*pow(Q2h,2)*pow(S,2) +
		           544*Q2e*pow(Q2h,2)*pow(S,2) + 368*pow(Q2h,3)*pow(S,2) +
		           320*l1k*m2*Q2e*pow(U,2) + 32*m4*Q2e*pow(U,2) +
		           768*l1k*m2*Q2h*pow(U,2) + 224*m4*Q2h*pow(U,2) +
		           896*l1k*Q2e*Q2h*pow(U,2) + 504*m2*Q2e*Q2h*pow(U,2) +
		           256*m2*pow(l1k,2)*pow(U,2) + 1376*Q2e*pow(l1k,2)*pow(U,2) +
		           3712*Q2h*pow(l1k,2)*pow(U,2) + 1536*pow(l1k,3)*pow(U,2) +
		           544*l1k*pow(Q2e,2)*pow(U,2) + 252*Q2h*pow(Q2e,2)*pow(U,2) +
		           60*pow(Q2e,3)*pow(U,2) + 2048*l1k*pow(Q2h,2)*pow(U,2) +
		           328*m2*pow(Q2h,2)*pow(U,2) + 320*Q2e*pow(Q2h,2)*pow(U,2) +
		           496*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*Q2h*S*U*pow(Q2e,3) + 8*M2*Q2h*pow(Q2e,4) +
		           12*Q2h*S*pow(Q2e,4) - 80*S*U*pow(Q2e,2)*pow(Q2h,2) -
		           8*M2*pow(Q2e,3)*pow(Q2h,2) + 72*S*pow(Q2e,3)*pow(Q2h,2) +
		           28*U*pow(Q2e,3)*pow(Q2h,2) - 8*pow(Q2e,4)*pow(Q2h,2) -
		           80*Q2e*S*U*pow(Q2h,3) + 80*S*pow(Q2e,2)*pow(Q2h,3) +
		           80*U*pow(Q2e,2)*pow(Q2h,3) - 40*pow(Q2e,3)*pow(Q2h,3) +
		           80*Q2e*S*pow(Q2h,4) + 80*Q2e*U*pow(Q2h,4) - 80*S*U*pow(Q2h,4) -
		           40*pow(Q2e,2)*pow(Q2h,4) - 40*Q2e*pow(Q2h,5) +
		           80*S*pow(Q2h,5) + 80*U*pow(Q2h,5) - 40*pow(Q2h,6) -
		           32*Q2h*pow(Q2e,3)*pow(S,2) - 4*pow(Q2e,4)*pow(S,2) -
		           40*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) -
		           40*Q2e*pow(Q2h,3)*pow(S,2) - 40*pow(Q2h,4)*pow(S,2) -
		           20*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) -
		           40*Q2e*pow(Q2h,3)*pow(U,2) - 40*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*(8*Q2h*S*U*pow(Q2e,3) - 4*M2*Q2h*pow(Q2e,4) -
		           6*Q2h*U*pow(Q2e,4) + 64*m4*Q2e*S*pow(Q2h,2) -
		           64*m4*Q2e*U*pow(Q2h,2) - 256*m2*Q2e*S*U*pow(Q2h,2) +
		           40*S*U*pow(Q2e,2)*pow(Q2h,2) + 4*M2*pow(Q2e,3)*pow(Q2h,2) -
		           14*S*pow(Q2e,3)*pow(Q2h,2) - 36*U*pow(Q2e,3)*pow(Q2h,2) +
		           4*pow(Q2e,4)*pow(Q2h,2) - 32*m4*S*pow(Q2h,3) +
		           480*m2*Q2e*S*pow(Q2h,3) + 32*m4*U*pow(Q2h,3) +
		           160*m2*Q2e*U*pow(Q2h,3) - 280*Q2e*S*U*pow(Q2h,3) -
		           40*S*pow(Q2e,2)*pow(Q2h,3) - 40*U*pow(Q2e,2)*pow(Q2h,3) +
		           20*pow(Q2e,3)*pow(Q2h,3) - 192*m2*Q2e*pow(Q2h,4) -
		           88*m2*S*pow(Q2h,4) + 372*Q2e*S*pow(Q2h,4) +
		           88*m2*U*pow(Q2h,4) + 316*Q2e*U*pow(Q2h,4) +
		           40*S*U*pow(Q2h,4) + 20*pow(Q2e,2)*pow(Q2h,4) -
		           204*Q2e*pow(Q2h,5) - 40*S*pow(Q2h,5) - 40*U*pow(Q2h,5) +
		           20*pow(Q2h,6) - 48*m4*Q2e*Q2h*pow(S,2) +
		           16*m4*pow(Q2h,2)*pow(S,2) - 260*m2*Q2e*pow(Q2h,2)*pow(S,2) +
		           10*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) +
		           60*m2*pow(Q2h,3)*pow(S,2) - 168*Q2e*pow(Q2h,3)*pow(S,2) +
		           20*pow(Q2h,4)*pow(S,2) + 48*m4*Q2e*Q2h*pow(U,2) +
		           16*Q2h*pow(Q2e,3)*pow(U,2) + 2*pow(Q2e,4)*pow(U,2) -
		           16*m4*pow(Q2h,2)*pow(U,2) + 4*m2*Q2e*pow(Q2h,2)*pow(U,2) +
		           20*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) -
		           60*m2*pow(Q2h,3)*pow(U,2) - 112*Q2e*pow(Q2h,3)*pow(U,2) +
		           20*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-2)*(128*m4*Q2e*S*U*pow(Q2h,2) -
		           168*m4*Q2e*S*pow(Q2h,3) - 152*m4*Q2e*U*pow(Q2h,3) -
		           64*m4*S*U*pow(Q2h,3) + 224*m2*Q2e*S*U*pow(Q2h,3) +
		           96*m4*Q2e*pow(Q2h,4) + 104*m4*S*pow(Q2h,4) -
		           270*m2*Q2e*S*pow(Q2h,4) + 88*m4*U*pow(Q2h,4) -
		           242*m2*Q2e*U*pow(Q2h,4) - 96*m2*S*U*pow(Q2h,4) +
		           64*Q2e*S*U*pow(Q2h,4) - 64*m4*pow(Q2h,5) +
		           144*m2*Q2e*pow(Q2h,5) + 142*m2*S*pow(Q2h,5) -
		           64*Q2e*S*pow(Q2h,5) + 114*m2*U*pow(Q2h,5) -
		           64*Q2e*U*pow(Q2h,5) - 80*m2*pow(Q2h,6) + 32*Q2e*pow(Q2h,6) +
		           72*m4*Q2e*pow(Q2h,2)*pow(S,2) - 40*m4*pow(Q2h,3)*pow(S,2) +
		           126*m2*Q2e*pow(Q2h,3)*pow(S,2) - 62*m2*pow(Q2h,4)*pow(S,2) +
		           32*Q2e*pow(Q2h,4)*pow(S,2) + 56*m4*Q2e*pow(Q2h,2)*pow(U,2) -
		           24*m4*pow(Q2h,3)*pow(U,2) + 98*m2*Q2e*pow(Q2h,3)*pow(U,2) -
		           34*m2*pow(Q2h,4)*pow(U,2) + 32*Q2e*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (40*S*U*pow(Q2h,5) - 40*S*pow(Q2h,6) - 40*U*pow(Q2h,6) +
		           20*pow(Q2h,7) + 20*pow(Q2h,5)*pow(S,2) + 20*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-3)*(-32*m4*Q2e*S*U*pow(Q2h,3) + 32*m4*Q2e*S*pow(Q2h,4) +
		           32*m4*Q2e*U*pow(Q2h,4) + 32*m4*S*U*pow(Q2h,4) -
		           32*m2*Q2e*S*U*pow(Q2h,4) - 16*m4*Q2e*pow(Q2h,5) -
		           32*m4*S*pow(Q2h,5) + 32*m2*Q2e*S*pow(Q2h,5) -
		           32*m4*U*pow(Q2h,5) + 32*m2*Q2e*U*pow(Q2h,5) +
		           32*m2*S*U*pow(Q2h,5) + 16*m4*pow(Q2h,6) -
		           16*m2*Q2e*pow(Q2h,6) - 32*m2*S*pow(Q2h,6) -
		           32*m2*U*pow(Q2h,6) + 16*m2*pow(Q2h,7) -
		           16*m4*Q2e*pow(Q2h,3)*pow(S,2) + 16*m4*pow(Q2h,4)*pow(S,2) -
		           16*m2*Q2e*pow(Q2h,4)*pow(S,2) + 16*m2*pow(Q2h,5)*pow(S,2) -
		           16*m4*Q2e*pow(Q2h,3)*pow(U,2) + 16*m4*pow(Q2h,4)*pow(U,2) -
		           16*m2*Q2e*pow(Q2h,4)*pow(U,2) + 16*m2*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1)) +
		     SI.D0_mm0QqMm0mm(Q2h,Q2e,-2*l1k + m2,m2)*
		      (16*l1k*M2*Q2e - 64*m2*M2*Q2e + 64*l1k*M2*Q2h -
		        96*m2*M2*Q2h - 56*M2*Q2e*Q2h + 120*l1k*Q2h*S - 272*m2*Q2h*S -
		        100*Q2e*Q2h*S - 424*l1k*Q2h*U - 224*m2*Q2h*U - 48*Q2e*Q2h*U +
		        64*l1k*S*U + 64*m2*S*U + 152*Q2h*S*U - 256*M2*pow(l1k,2) -
		        16*M2*pow(Q2e,2) + 112*l1k*pow(Q2h,2) + 208*m2*pow(Q2h,2) -
		        48*M2*pow(Q2h,2) + 76*Q2e*pow(Q2h,2) - 444*S*pow(Q2h,2) -
		        268*U*pow(Q2h,2) + 284*pow(Q2h,3) - 80*l1k*pow(S,2) +
		        96*m2*pow(S,2) + 40*Q2e*pow(S,2) + 168*Q2h*pow(S,2) +
		        208*l1k*pow(U,2) + 64*m2*pow(U,2) + 16*Q2e*pow(U,2) +
		        48*Q2h*pow(U,2) + pow(l1k,-1)*
		         (-32*m4*M2*Q2h - 16*m2*M2*Q2e*Q2h - 32*m4*Q2h*U -
		           16*m2*Q2e*Q2h*U + 32*m4*pow(U,2) + 16*m2*Q2e*pow(U,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (64*m4*M2*Q2h + 64*m2*M2*Q2e*Q2h + 128*m4*Q2h*S +
		           256*m2*Q2e*Q2h*S + 80*m2*Q2e*Q2h*U - 320*m2*Q2h*S*U -
		           136*Q2e*Q2h*S*U + 80*m2*M2*pow(Q2e,2) +
		           40*M2*Q2h*pow(Q2e,2) + 108*Q2h*S*pow(Q2e,2) +
		           24*Q2h*U*pow(Q2e,2) + 24*M2*pow(Q2e,3) - 32*m4*pow(Q2h,2) -
		           16*m2*M2*pow(Q2h,2) - 144*m2*Q2e*pow(Q2h,2) +
		           544*m2*S*pow(Q2h,2) + 336*Q2e*S*pow(Q2h,2) +
		           384*m2*U*pow(Q2h,2) + 204*Q2e*U*pow(Q2h,2) -
		           400*S*U*pow(Q2h,2) - 68*pow(Q2e,2)*pow(Q2h,2) -
		           304*m2*pow(Q2h,3) - 200*Q2e*pow(Q2h,3) + 592*S*pow(Q2h,3) +
		           464*U*pow(Q2h,3) - 328*pow(Q2h,4) - 64*m4*pow(S,2) -
		           96*m2*Q2e*pow(S,2) - 240*m2*Q2h*pow(S,2) -
		           136*Q2e*Q2h*pow(S,2) - 40*pow(Q2e,2)*pow(S,2) -
		           264*pow(Q2h,2)*pow(S,2) - 80*m2*Q2h*pow(U,2) -
		           24*Q2e*Q2h*pow(U,2) - 136*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (64*m4*Q2h*S*U - 96*m4*S*pow(Q2h,2) - 32*m4*U*pow(Q2h,2) +
		           304*m2*S*U*pow(Q2h,2) + 32*m4*pow(Q2h,3) -
		           408*m2*S*pow(Q2h,3) - 328*m2*U*pow(Q2h,3) +
		           328*S*U*pow(Q2h,3) + 216*m2*pow(Q2h,4) - 424*S*pow(Q2h,4) -
		           360*U*pow(Q2h,4) + 228*pow(Q2h,5) + 64*m4*Q2h*pow(S,2) +
		           192*m2*pow(Q2h,2)*pow(S,2) + 196*pow(Q2h,3)*pow(S,2) +
		           112*m2*pow(Q2h,2)*pow(U,2) + 132*pow(Q2h,3)*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-64*m2*S*U*pow(Q2h,3) + 64*m2*S*pow(Q2h,4) +
		           64*m2*U*pow(Q2h,4) - 64*S*U*pow(Q2h,4) - 32*m2*pow(Q2h,5) +
		           64*S*pow(Q2h,5) + 64*U*pow(Q2h,5) - 32*pow(Q2h,6) -
		           32*m2*pow(Q2h,3)*pow(S,2) - 32*pow(Q2h,4)*pow(S,2) -
		           32*m2*pow(Q2h,3)*pow(U,2) - 32*pow(Q2h,4)*pow(U,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*S*U*pow(Q2e,6)*pow(Q2h,2) + 8*S*pow(Q2e,7)*pow(Q2h,2) +
		           16*S*U*pow(Q2e,5)*pow(Q2h,3) - 16*S*pow(Q2e,6)*pow(Q2h,3) +
		           8*U*pow(Q2e,6)*pow(Q2h,3) - 4*pow(Q2e,7)*pow(Q2h,3) -
		           8*S*U*pow(Q2e,4)*pow(Q2h,4) + 8*S*pow(Q2e,5)*pow(Q2h,4) -
		           16*U*pow(Q2e,5)*pow(Q2h,4) + 8*pow(Q2e,6)*pow(Q2h,4) +
		           8*U*pow(Q2e,4)*pow(Q2h,5) - 4*pow(Q2e,5)*pow(Q2h,5) -
		           4*Q2h*pow(Q2e,7)*pow(S,2) + 8*pow(Q2e,6)*pow(Q2h,2)*pow(S,2) -
		           4*pow(Q2e,5)*pow(Q2h,3)*pow(S,2) -
		           4*pow(Q2e,5)*pow(Q2h,3)*pow(U,2) +
		           8*pow(Q2e,4)*pow(Q2h,4)*pow(U,2) -
		           4*pow(Q2e,3)*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        (1024*m4*Q2e*Q2h*S*U*pow(l1k,3) - 4096*m4*Q2e*Q2h*S*pow(l1k,4) -
		           4096*m4*Q2e*Q2h*U*pow(l1k,4) + 2048*m4*Q2e*S*U*pow(l1k,4) +
		           11264*m2*Q2e*Q2h*S*U*pow(l1k,4) -
		           4096*m2*Q2e*Q2h*S*pow(l1k,5) -
		           12288*m2*Q2e*Q2h*U*pow(l1k,5) + 4096*m2*Q2e*S*U*pow(l1k,5) +
		           9216*Q2e*Q2h*S*U*pow(l1k,5) - 16384*Q2e*Q2h*S*pow(l1k,6) -
		           40960*Q2e*Q2h*U*pow(l1k,6) + 12288*Q2e*S*U*pow(l1k,6) +
		           2048*Q2e*Q2h*S*U*pow(l1k,2)*pow(m,6) +
		           8448*Q2h*S*U*pow(l1k,4)*pow(Q2e,2) -
		           10240*Q2h*S*pow(l1k,5)*pow(Q2e,2) -
		           30720*Q2h*U*pow(l1k,5)*pow(Q2e,2) +
		           8192*S*U*pow(l1k,5)*pow(Q2e,2) +
		           768*Q2h*S*U*pow(l1k,3)*pow(Q2e,3) -
		           1280*Q2h*S*pow(l1k,4)*pow(Q2e,3) -
		           6400*Q2h*U*pow(l1k,4)*pow(Q2e,3) +
		           1152*S*U*pow(l1k,4)*pow(Q2e,3) +
		           32*Q2h*S*U*pow(l1k,2)*pow(Q2e,4) -
		           64*Q2h*S*pow(l1k,3)*pow(Q2e,4) -
		           704*Q2h*U*pow(l1k,3)*pow(Q2e,4) +
		           64*S*U*pow(l1k,3)*pow(Q2e,4) - 32*Q2h*U*pow(l1k,2)*pow(Q2e,5) +
		           4352*m4*Q2e*S*U*pow(l1k,2)*pow(Q2h,2) -
		           4096*m4*Q2e*U*pow(l1k,3)*pow(Q2h,2) -
		           1024*m4*S*U*pow(l1k,3)*pow(Q2h,2) +
		           2560*m2*Q2e*S*U*pow(l1k,3)*pow(Q2h,2) +
		           4096*m4*Q2e*pow(l1k,4)*pow(Q2h,2) -
		           16384*m2*Q2e*S*pow(l1k,4)*pow(Q2h,2) -
		           26624*m2*Q2e*U*pow(l1k,4)*pow(Q2h,2) -
		           3072*m2*S*U*pow(l1k,4)*pow(Q2h,2) -
		           2304*Q2e*S*U*pow(l1k,4)*pow(Q2h,2) +
		           8192*m2*Q2e*pow(l1k,5)*pow(Q2h,2) -
		           8192*Q2e*S*pow(l1k,5)*pow(Q2h,2) -
		           16384*Q2e*U*pow(l1k,5)*pow(Q2h,2) +
		           28672*Q2e*pow(l1k,6)*pow(Q2h,2) -
		           3072*Q2e*S*pow(l1k,2)*pow(m,6)*pow(Q2h,2) -
		           3072*Q2e*U*pow(l1k,2)*pow(m,6)*pow(Q2h,2) +
		           1280*S*U*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) -
		           9984*S*pow(l1k,4)*pow(Q2e,2)*pow(Q2h,2) -
		           17152*U*pow(l1k,4)*pow(Q2e,2)*pow(Q2h,2) +
		           20480*pow(l1k,5)*pow(Q2e,2)*pow(Q2h,2) +
		           64*S*U*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2) -
		           768*S*pow(l1k,3)*pow(Q2e,3)*pow(Q2h,2) -
		           1280*U*pow(l1k,3)*pow(Q2e,3)*pow(Q2h,2) +
		           3840*pow(l1k,4)*pow(Q2e,3)*pow(Q2h,2) -
		           16*l1k*S*U*pow(Q2e,4)*pow(Q2h,2) -
		           64*S*pow(l1k,2)*pow(Q2e,4)*pow(Q2h,2) -
		           32*U*pow(l1k,2)*pow(Q2e,4)*pow(Q2h,2) +
		           384*pow(l1k,3)*pow(Q2e,4)*pow(Q2h,2) +
		           16*l1k*S*pow(Q2e,5)*pow(Q2h,2) + 8*S*U*pow(Q2e,5)*pow(Q2h,2) +
		           16*pow(l1k,2)*pow(Q2e,5)*pow(Q2h,2) -
		           8*S*pow(Q2e,6)*pow(Q2h,2) - 1024*l1k*m4*Q2e*S*U*pow(Q2h,3) -
		           5632*m4*Q2e*S*pow(l1k,2)*pow(Q2h,3) -
		           6144*m4*Q2e*U*pow(l1k,2)*pow(Q2h,3) -
		           1792*m4*S*U*pow(l1k,2)*pow(Q2h,3) +
		           3776*m2*Q2e*S*U*pow(l1k,2)*pow(Q2h,3) +
		           2048*m4*Q2e*pow(l1k,3)*pow(Q2h,3) +
		           1024*m4*S*pow(l1k,3)*pow(Q2h,3) -
		           512*m2*Q2e*S*pow(l1k,3)*pow(Q2h,3) +
		           3072*m4*U*pow(l1k,3)*pow(Q2h,3) -
		           8704*m2*Q2e*U*pow(l1k,3)*pow(Q2h,3) -
		           3840*m2*S*U*pow(l1k,3)*pow(Q2h,3) -
		           3712*Q2e*S*U*pow(l1k,3)*pow(Q2h,3) +
		           18432*m2*Q2e*pow(l1k,4)*pow(Q2h,3) +
		           4096*m2*S*pow(l1k,4)*pow(Q2h,3) +
		           768*Q2e*S*pow(l1k,4)*pow(Q2h,3) +
		           10240*m2*U*pow(l1k,4)*pow(Q2h,3) +
		           9984*Q2e*U*pow(l1k,4)*pow(Q2h,3) +
		           6144*Q2e*pow(l1k,5)*pow(Q2h,3) +
		           256*Q2e*S*U*pow(m,6)*pow(Q2h,3) +
		           2048*Q2e*pow(l1k,2)*pow(m,6)*pow(Q2h,3) +
		           1024*S*pow(l1k,2)*pow(m,6)*pow(Q2h,3) +
		           1024*U*pow(l1k,2)*pow(m,6)*pow(Q2h,3) +
		           64*S*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3) -
		           1280*S*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,3) -
		           1792*U*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,3) +
		           9728*pow(l1k,4)*pow(Q2e,2)*pow(Q2h,3) -
		           64*S*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,3) -
		           64*U*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,3) +
		           640*pow(l1k,3)*pow(Q2e,3)*pow(Q2h,3) +
		           16*l1k*U*pow(Q2e,4)*pow(Q2h,3) - 8*S*U*pow(Q2e,4)*pow(Q2h,3) +
		           32*pow(l1k,2)*pow(Q2e,4)*pow(Q2h,3) -
		           8*l1k*pow(Q2e,5)*pow(Q2h,3) + 8*S*pow(Q2e,5)*pow(Q2h,3) -
		           8*U*pow(Q2e,5)*pow(Q2h,3) + 4*pow(Q2e,6)*pow(Q2h,3) +
		           1280*l1k*m4*Q2e*S*pow(Q2h,4) + 1024*l1k*m4*Q2e*U*pow(Q2h,4) +
		           256*l1k*m4*S*U*pow(Q2h,4) - 896*l1k*m2*Q2e*S*U*pow(Q2h,4) +
		           256*m4*Q2e*S*U*pow(Q2h,4) +
		           3840*m4*Q2e*pow(l1k,2)*pow(Q2h,4) +
		           2560*m4*S*pow(l1k,2)*pow(Q2h,4) -
		           5056*m2*Q2e*S*pow(l1k,2)*pow(Q2h,4) +
		           3072*m4*U*pow(l1k,2)*pow(Q2h,4) -
		           4800*m2*Q2e*U*pow(l1k,2)*pow(Q2h,4) -
		           1344*m2*S*U*pow(l1k,2)*pow(Q2h,4) +
		           832*Q2e*S*U*pow(l1k,2)*pow(Q2h,4) -
		           2048*m4*pow(l1k,3)*pow(Q2h,4) +
		           3584*m2*Q2e*pow(l1k,3)*pow(Q2h,4) +
		           4096*m2*S*pow(l1k,3)*pow(Q2h,4) +
		           5248*Q2e*S*pow(l1k,3)*pow(Q2h,4) +
		           9216*m2*U*pow(l1k,3)*pow(Q2h,4) +
		           5248*Q2e*U*pow(l1k,3)*pow(Q2h,4) -
		           7168*m2*pow(l1k,4)*pow(Q2h,4) -
		           4864*Q2e*pow(l1k,4)*pow(Q2h,4) -
		           256*Q2e*S*pow(m,6)*pow(Q2h,4) -
		           256*Q2e*U*pow(m,6)*pow(Q2h,4) - 256*S*U*pow(m,6)*pow(Q2h,4) -
		           1024*pow(l1k,2)*pow(m,6)*pow(Q2h,4) -
		           64*S*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,4) -
		           64*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,4) +
		           896*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,4) +
		           32*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,4) +
		           8*U*pow(Q2e,4)*pow(Q2h,4) - 4*pow(Q2e,5)*pow(Q2h,4) -
		           640*l1k*m4*Q2e*pow(Q2h,5) - 512*l1k*m4*S*pow(Q2h,5) +
		           1024*l1k*m2*Q2e*S*pow(Q2h,5) - 256*m4*Q2e*S*pow(Q2h,5) -
		           256*l1k*m4*U*pow(Q2h,5) + 896*l1k*m2*Q2e*U*pow(Q2h,5) -
		           256*m4*Q2e*U*pow(Q2h,5) + 256*l1k*m2*S*U*pow(Q2h,5) -
		           256*m4*S*U*pow(Q2h,5) - 160*l1k*Q2e*S*U*pow(Q2h,5) +
		           80*m2*Q2e*S*U*pow(Q2h,5) - 2048*m4*pow(l1k,2)*pow(Q2h,5) +
		           3072*m2*Q2e*pow(l1k,2)*pow(Q2h,5) +
		           1856*m2*S*pow(l1k,2)*pow(Q2h,5) -
		           1088*Q2e*S*pow(l1k,2)*pow(Q2h,5) +
		           2112*m2*U*pow(l1k,2)*pow(Q2h,5) -
		           832*Q2e*U*pow(l1k,2)*pow(Q2h,5) -
		           5120*m2*pow(l1k,3)*pow(Q2h,5) -
		           3456*Q2e*pow(l1k,3)*pow(Q2h,5) + 128*Q2e*pow(m,6)*pow(Q2h,5) +
		           256*S*pow(m,6)*pow(Q2h,5) + 256*U*pow(m,6)*pow(Q2h,5) +
		           32*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,5) + 256*l1k*m4*pow(Q2h,6) -
		           512*l1k*m2*Q2e*pow(Q2h,6) + 128*m4*Q2e*pow(Q2h,6) -
		           384*l1k*m2*S*pow(Q2h,6) + 256*m4*S*pow(Q2h,6) +
		           160*l1k*Q2e*S*pow(Q2h,6) - 80*m2*Q2e*S*pow(Q2h,6) -
		           256*l1k*m2*U*pow(Q2h,6) + 256*m4*U*pow(Q2h,6) +
		           160*l1k*Q2e*U*pow(Q2h,6) - 80*m2*Q2e*U*pow(Q2h,6) -
		           80*m2*S*U*pow(Q2h,6) - 1344*m2*pow(l1k,2)*pow(Q2h,6) +
		           544*Q2e*pow(l1k,2)*pow(Q2h,6) - 128*pow(m,6)*pow(Q2h,6) +
		           192*l1k*m2*pow(Q2h,7) - 128*m4*pow(Q2h,7) -
		           80*l1k*Q2e*pow(Q2h,7) + 40*m2*Q2e*pow(Q2h,7) +
		           80*m2*S*pow(Q2h,7) + 80*m2*U*pow(Q2h,7) - 40*m2*pow(Q2h,8) -
		           1024*m4*Q2e*Q2h*pow(l1k,3)*pow(S,2) +
		           1024*m4*Q2e*pow(l1k,4)*pow(S,2) +
		           3584*m2*Q2e*Q2h*pow(l1k,4)*pow(S,2) +
		           512*Q2e*Q2h*pow(l1k,5)*pow(S,2) +
		           2048*Q2e*pow(l1k,6)*pow(S,2) +
		           1024*Q2e*Q2h*pow(l1k,2)*pow(m,6)*pow(S,2) +
		           1408*Q2h*pow(l1k,4)*pow(Q2e,2)*pow(S,2) +
		           1024*pow(l1k,5)*pow(Q2e,2)*pow(S,2) +
		           32*Q2h*pow(l1k,3)*pow(Q2e,3)*pow(S,2) +
		           64*pow(l1k,4)*pow(Q2e,3)*pow(S,2) +
		           16*Q2h*pow(l1k,2)*pow(Q2e,4)*pow(S,2) -
		           8*l1k*Q2h*pow(Q2e,5)*pow(S,2) + 4*Q2h*pow(Q2e,6)*pow(S,2) +
		           1920*m4*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(S,2) -
		           1536*m2*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(S,2) -
		           512*m2*pow(l1k,4)*pow(Q2h,2)*pow(S,2) +
		           1664*Q2e*pow(l1k,4)*pow(Q2h,2)*pow(S,2) +
		           384*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) +
		           32*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2)*pow(S,2) -
		           4*pow(Q2e,5)*pow(Q2h,2)*pow(S,2) -
		           640*l1k*m4*Q2e*pow(Q2h,3)*pow(S,2) -
		           640*m4*pow(l1k,2)*pow(Q2h,3)*pow(S,2) +
		           2080*m2*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(S,2) -
		           384*m2*pow(l1k,3)*pow(Q2h,3)*pow(S,2) -
		           1984*Q2e*pow(l1k,3)*pow(Q2h,3)*pow(S,2) +
		           128*Q2e*pow(m,6)*pow(Q2h,3)*pow(S,2) +
		           32*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3)*pow(S,2) +
		           256*l1k*m4*pow(Q2h,4)*pow(S,2) -
		           512*l1k*m2*Q2e*pow(Q2h,4)*pow(S,2) +
		           128*m4*Q2e*pow(Q2h,4)*pow(S,2) -
		           608*m2*pow(l1k,2)*pow(Q2h,4)*pow(S,2) +
		           544*Q2e*pow(l1k,2)*pow(Q2h,4)*pow(S,2) -
		           128*pow(m,6)*pow(Q2h,4)*pow(S,2) +
		           192*l1k*m2*pow(Q2h,5)*pow(S,2) - 128*m4*pow(Q2h,5)*pow(S,2) -
		           80*l1k*Q2e*pow(Q2h,5)*pow(S,2) +
		           40*m2*Q2e*pow(Q2h,5)*pow(S,2) - 40*m2*pow(Q2h,6)*pow(S,2) +
		           2048*m4*Q2e*Q2h*pow(l1k,3)*pow(U,2) +
		           1024*m4*Q2e*pow(l1k,4)*pow(U,2) +
		           9728*m2*Q2e*Q2h*pow(l1k,4)*pow(U,2) +
		           4096*m2*Q2e*pow(l1k,5)*pow(U,2) +
		           8704*Q2e*Q2h*pow(l1k,5)*pow(U,2) +
		           14336*Q2e*pow(l1k,6)*pow(U,2) +
		           1024*Q2e*Q2h*pow(l1k,2)*pow(m,6)*pow(U,2) +
		           7552*Q2h*pow(l1k,4)*pow(Q2e,2)*pow(U,2) +
		           11264*pow(l1k,5)*pow(Q2e,2)*pow(U,2) +
		           608*Q2h*pow(l1k,3)*pow(Q2e,3)*pow(U,2) +
		           2624*pow(l1k,4)*pow(Q2e,3)*pow(U,2) +
		           16*Q2h*pow(l1k,2)*pow(Q2e,4)*pow(U,2) +
		           320*pow(l1k,3)*pow(Q2e,4)*pow(U,2) +
		           16*pow(l1k,2)*pow(Q2e,5)*pow(U,2) +
		           2432*m4*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(U,2) -
		           1024*m4*pow(l1k,3)*pow(Q2h,2)*pow(U,2) +
		           4608*m2*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(U,2) -
		           3584*m2*pow(l1k,4)*pow(Q2h,2)*pow(U,2) -
		           4480*Q2e*pow(l1k,4)*pow(Q2h,2)*pow(U,2) +
		           896*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) +
		           16*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2)*pow(U,2) -
		           384*l1k*m4*Q2e*pow(Q2h,3)*pow(U,2) -
		           1152*m4*pow(l1k,2)*pow(Q2h,3)*pow(U,2) +
		           1824*m2*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(U,2) -
		           3968*m2*pow(l1k,3)*pow(Q2h,3)*pow(U,2) -
		           1984*Q2e*pow(l1k,3)*pow(Q2h,3)*pow(U,2) +
		           128*Q2e*pow(m,6)*pow(Q2h,3)*pow(U,2) +
		           32*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3)*pow(U,2) -
		           8*l1k*pow(Q2e,3)*pow(Q2h,3)*pow(U,2) +
		           4*pow(Q2e,4)*pow(Q2h,3)*pow(U,2) -
		           384*l1k*m2*Q2e*pow(Q2h,4)*pow(U,2) +
		           128*m4*Q2e*pow(Q2h,4)*pow(U,2) -
		           864*m2*pow(l1k,2)*pow(Q2h,4)*pow(U,2) +
		           288*Q2e*pow(l1k,2)*pow(Q2h,4)*pow(U,2) -
		           128*pow(m,6)*pow(Q2h,4)*pow(U,2) -
		           4*pow(Q2e,3)*pow(Q2h,4)*pow(U,2) +
		           64*l1k*m2*pow(Q2h,5)*pow(U,2) - 128*m4*pow(Q2h,5)*pow(U,2) -
		           80*l1k*Q2e*pow(Q2h,5)*pow(U,2) +
		           40*m2*Q2e*pow(Q2h,5)*pow(U,2) - 40*m2*pow(Q2h,6)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (16*Q2h*S*U*pow(Q2e,4) - 8*M2*Q2h*pow(Q2e,5) -
		           12*Q2h*S*pow(Q2e,5) + 64*S*U*pow(Q2e,3)*pow(Q2h,2) +
		           16*M2*pow(Q2e,4)*pow(Q2h,2) - 60*S*pow(Q2e,4)*pow(Q2h,2) -
		           28*U*pow(Q2e,4)*pow(Q2h,2) + 8*pow(Q2e,5)*pow(Q2h,2) -
		           8*M2*pow(Q2e,3)*pow(Q2h,3) - 8*S*pow(Q2e,3)*pow(Q2h,3) -
		           52*U*pow(Q2e,3)*pow(Q2h,3) + 32*pow(Q2e,4)*pow(Q2h,3) +
		           28*Q2h*pow(Q2e,4)*pow(S,2) + 4*pow(Q2e,5)*pow(S,2) +
		           8*pow(Q2e,3)*pow(Q2h,2)*pow(S,2) +
		           20*pow(Q2e,3)*pow(Q2h,2)*pow(U,2) +
		           20*pow(Q2e,2)*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        (-480*l1k*m4*Q2e*Q2h*S + 160*l1k*m4*Q2e*Q2h*U +
		           128*l1k*m4*Q2e*S*U + 640*l1k*m2*Q2e*Q2h*S*U -
		           608*m4*Q2e*Q2h*S*U + 512*m4*M2*Q2e*pow(l1k,2) +
		           768*m2*M2*Q2e*Q2h*pow(l1k,2) + 1024*m4*Q2h*S*pow(l1k,2) +
		           1536*m2*Q2e*Q2h*S*pow(l1k,2) + 1024*m4*Q2h*U*pow(l1k,2) +
		           1856*m2*Q2e*Q2h*U*pow(l1k,2) - 256*m4*S*U*pow(l1k,2) -
		           448*m2*Q2e*S*U*pow(l1k,2) - 1792*m2*Q2h*S*U*pow(l1k,2) -
		           928*Q2e*Q2h*S*U*pow(l1k,2) - 768*M2*Q2e*Q2h*pow(l1k,3) +
		           1024*m2*Q2h*S*pow(l1k,3) + 32*Q2e*Q2h*S*pow(l1k,3) +
		           3072*m2*Q2h*U*pow(l1k,3) + 4256*Q2e*Q2h*U*pow(l1k,3) -
		           1024*m2*S*U*pow(l1k,3) - 768*Q2e*S*U*pow(l1k,3) -
		           3840*Q2h*S*U*pow(l1k,3) + 1024*M2*Q2e*pow(l1k,4) +
		           4096*Q2h*S*pow(l1k,4) + 10240*Q2h*U*pow(l1k,4) -
		           3072*S*U*pow(l1k,4) + 192*Q2e*Q2h*S*pow(m,6) +
		           192*Q2e*Q2h*U*pow(m,6) - 128*Q2e*S*U*pow(m,6) -
		           128*Q2h*S*U*pow(m,6) - 32*l1k*Q2h*S*U*pow(Q2e,2) +
		           64*M2*Q2h*pow(l1k,2)*pow(Q2e,2) +
		           224*Q2h*S*pow(l1k,2)*pow(Q2e,2) +
		           1328*Q2h*U*pow(l1k,2)*pow(Q2e,2) -
		           160*S*U*pow(l1k,2)*pow(Q2e,2) +
		           512*M2*pow(l1k,3)*pow(Q2e,2) - 16*l1k*M2*Q2h*pow(Q2e,3) -
		           24*l1k*Q2h*S*pow(Q2e,3) + 72*l1k*Q2h*U*pow(Q2e,3) -
		           16*Q2h*S*U*pow(Q2e,3) + 32*M2*pow(l1k,2)*pow(Q2e,3) +
		           8*M2*Q2h*pow(Q2e,4) + 12*Q2h*S*pow(Q2e,4) +
		           64*l1k*m4*Q2e*pow(Q2h,2) - 256*l1k*m2*M2*Q2e*pow(Q2h,2) +
		           64*m4*M2*Q2e*pow(Q2h,2) + 32*l1k*m4*S*pow(Q2h,2) -
		           2120*l1k*m2*Q2e*S*pow(Q2h,2) + 1056*m4*Q2e*S*pow(Q2h,2) +
		           32*l1k*m4*U*pow(Q2h,2) - 328*l1k*m2*Q2e*U*pow(Q2h,2) +
		           800*m4*Q2e*U*pow(Q2h,2) + 320*l1k*m2*S*U*pow(Q2h,2) -
		           224*m4*S*U*pow(Q2h,2) + 1232*l1k*Q2e*S*U*pow(Q2h,2) -
		           1176*m2*Q2e*S*U*pow(Q2h,2) - 1024*m4*pow(l1k,2)*pow(Q2h,2) -
		           256*m2*M2*pow(l1k,2)*pow(Q2h,2) -
		           1600*m2*Q2e*pow(l1k,2)*pow(Q2h,2) +
		           320*M2*Q2e*pow(l1k,2)*pow(Q2h,2) +
		           2560*m2*S*pow(l1k,2)*pow(Q2h,2) +
		           2720*Q2e*S*pow(l1k,2)*pow(Q2h,2) +
		           3072*m2*U*pow(l1k,2)*pow(Q2h,2) +
		           1184*Q2e*U*pow(l1k,2)*pow(Q2h,2) -
		           1344*S*U*pow(l1k,2)*pow(Q2h,2) -
		           2048*m2*pow(l1k,3)*pow(Q2h,2) -
		           1984*Q2e*pow(l1k,3)*pow(Q2h,2) + 4096*S*pow(l1k,3)*pow(Q2h,2) +
		           9216*U*pow(l1k,3)*pow(Q2h,2) - 7168*pow(l1k,4)*pow(Q2h,2) -
		           128*Q2e*pow(m,6)*pow(Q2h,2) + 64*S*pow(m,6)*pow(Q2h,2) +
		           64*U*pow(m,6)*pow(Q2h,2) - 112*l1k*S*pow(Q2e,2)*pow(Q2h,2) +
		           48*l1k*U*pow(Q2e,2)*pow(Q2h,2) - 80*S*U*pow(Q2e,2)*pow(Q2h,2) -
		           720*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2) -
		           16*l1k*pow(Q2e,3)*pow(Q2h,2) - 8*M2*pow(Q2e,3)*pow(Q2h,2) +
		           72*S*pow(Q2e,3)*pow(Q2h,2) + 28*U*pow(Q2e,3)*pow(Q2h,2) -
		           8*pow(Q2e,4)*pow(Q2h,2) + 64*l1k*m2*M2*pow(Q2h,3) -
		           64*m4*M2*pow(Q2h,3) + 928*l1k*m2*Q2e*pow(Q2h,3) -
		           624*m4*Q2e*pow(Q2h,3) - 96*l1k*M2*Q2e*pow(Q2h,3) +
		           48*m2*M2*Q2e*pow(Q2h,3) - 392*l1k*m2*S*pow(Q2h,3) -
		           2328*l1k*Q2e*S*pow(Q2h,3) + 1860*m2*Q2e*S*pow(Q2h,3) -
		           680*l1k*m2*U*pow(Q2h,3) + 64*m4*U*pow(Q2h,3) -
		           1592*l1k*Q2e*U*pow(Q2h,3) + 1412*m2*Q2e*U*pow(Q2h,3) +
		           256*l1k*S*U*pow(Q2h,3) - 104*m2*S*U*pow(Q2h,3) -
		           736*Q2e*S*U*pow(Q2h,3) - 2048*m2*pow(l1k,2)*pow(Q2h,3) -
		           1552*Q2e*pow(l1k,2)*pow(Q2h,3) + 1856*S*pow(l1k,2)*pow(Q2h,3) +
		           2112*U*pow(l1k,2)*pow(Q2h,3) - 5120*pow(l1k,3)*pow(Q2h,3) +
		           48*l1k*pow(Q2e,2)*pow(Q2h,3) + 80*S*pow(Q2e,2)*pow(Q2h,3) +
		           80*U*pow(Q2e,2)*pow(Q2h,3) - 40*pow(Q2e,3)*pow(Q2h,3) +
		           368*l1k*m2*pow(Q2h,4) + 80*m4*pow(Q2h,4) -
		           48*m2*M2*pow(Q2h,4) + 1352*l1k*Q2e*pow(Q2h,4) -
		           1052*m2*Q2e*pow(Q2h,4) - 384*l1k*S*pow(Q2h,4) -
		           188*m2*S*pow(Q2h,4) + 928*Q2e*S*pow(Q2h,4) -
		           256*l1k*U*pow(Q2h,4) - 12*m2*U*pow(Q2h,4) +
		           800*Q2e*U*pow(Q2h,4) - 80*S*U*pow(Q2h,4) -
		           1344*pow(l1k,2)*pow(Q2h,4) - 40*pow(Q2e,2)*pow(Q2h,4) +
		           192*l1k*pow(Q2h,5) + 156*m2*pow(Q2h,5) - 496*Q2e*pow(Q2h,5) +
		           80*S*pow(Q2h,5) + 80*U*pow(Q2h,5) - 40*pow(Q2h,6) +
		           256*l1k*m4*Q2e*pow(S,2) + 64*l1k*m4*Q2h*pow(S,2) +
		           976*l1k*m2*Q2e*Q2h*pow(S,2) - 448*m4*Q2e*Q2h*pow(S,2) -
		           128*m4*pow(l1k,2)*pow(S,2) - 480*m2*Q2e*pow(l1k,2)*pow(S,2) -
		           640*m2*Q2h*pow(l1k,2)*pow(S,2) -
		           992*Q2e*Q2h*pow(l1k,2)*pow(S,2) + 320*Q2e*pow(l1k,3)*pow(S,2) -
		           384*Q2h*pow(l1k,3)*pow(S,2) - 512*pow(l1k,4)*pow(S,2) -
		           64*Q2e*pow(m,6)*pow(S,2) - 64*Q2h*pow(m,6)*pow(S,2) +
		           72*l1k*Q2h*pow(Q2e,2)*pow(S,2) -
		           16*pow(l1k,2)*pow(Q2e,2)*pow(S,2) + 8*l1k*pow(Q2e,3)*pow(S,2) -
		           32*Q2h*pow(Q2e,3)*pow(S,2) - 4*pow(Q2e,4)*pow(S,2) +
		           176*l1k*m2*pow(Q2h,2)*pow(S,2) - 64*m4*pow(Q2h,2)*pow(S,2) +
		           992*l1k*Q2e*pow(Q2h,2)*pow(S,2) -
		           816*m2*Q2e*pow(Q2h,2)*pow(S,2) -
		           608*pow(l1k,2)*pow(Q2h,2)*pow(S,2) -
		           40*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) +
		           192*l1k*pow(Q2h,3)*pow(S,2) + 40*m2*pow(Q2h,3)*pow(S,2) -
		           432*Q2e*pow(Q2h,3)*pow(S,2) - 40*pow(Q2h,4)*pow(S,2) -
		           128*l1k*m4*Q2e*pow(U,2) - 64*l1k*m4*Q2h*pow(U,2) -
		           272*l1k*m2*Q2e*Q2h*pow(U,2) - 224*m4*Q2e*Q2h*pow(U,2) -
		           128*m4*pow(l1k,2)*pow(U,2) - 608*m2*Q2e*pow(l1k,2)*pow(U,2) -
		           1152*m2*Q2h*pow(l1k,2)*pow(U,2) -
		           192*Q2e*Q2h*pow(l1k,2)*pow(U,2) -
		           1024*m2*pow(l1k,3)*pow(U,2) - 1856*Q2e*pow(l1k,3)*pow(U,2) -
		           3968*Q2h*pow(l1k,3)*pow(U,2) - 3584*pow(l1k,4)*pow(U,2) -
		           64*Q2e*pow(m,6)*pow(U,2) - 64*Q2h*pow(m,6)*pow(U,2) -
		           56*l1k*Q2h*pow(Q2e,2)*pow(U,2) -
		           624*pow(l1k,2)*pow(Q2e,2)*pow(U,2) -
		           40*l1k*pow(Q2e,3)*pow(U,2) + 208*l1k*m2*pow(Q2h,2)*pow(U,2) -
		           96*m4*pow(Q2h,2)*pow(U,2) + 368*l1k*Q2e*pow(Q2h,2)*pow(U,2) -
		           424*m2*Q2e*pow(Q2h,2)*pow(U,2) -
		           864*pow(l1k,2)*pow(Q2h,2)*pow(U,2) -
		           20*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) +
		           64*l1k*pow(Q2h,3)*pow(U,2) - 80*m2*pow(Q2h,3)*pow(U,2) -
		           304*Q2e*pow(Q2h,3)*pow(U,2) - 40*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*(64*Q2e*Q2h*S*U*pow(m,6) +
		           304*m4*Q2e*S*U*pow(Q2h,2) - 96*Q2e*S*pow(m,6)*pow(Q2h,2) -
		           32*Q2e*U*pow(m,6)*pow(Q2h,2) - 64*S*U*pow(m,6)*pow(Q2h,2) -
		           408*m4*Q2e*S*pow(Q2h,3) - 328*m4*Q2e*U*pow(Q2h,3) -
		           176*m4*S*U*pow(Q2h,3) + 456*m2*Q2e*S*U*pow(Q2h,3) +
		           32*Q2e*pow(m,6)*pow(Q2h,3) + 96*S*pow(m,6)*pow(Q2h,3) +
		           32*U*pow(m,6)*pow(Q2h,3) + 216*m4*Q2e*pow(Q2h,4) +
		           280*m4*S*pow(Q2h,4) - 552*m2*Q2e*S*pow(Q2h,4) +
		           200*m4*U*pow(Q2h,4) - 488*m2*Q2e*U*pow(Q2h,4) -
		           200*m2*S*U*pow(Q2h,4) + 128*Q2e*S*U*pow(Q2h,4) -
		           32*pow(m,6)*pow(Q2h,4) - 152*m4*pow(Q2h,5) +
		           292*m2*Q2e*pow(Q2h,5) + 296*m2*S*pow(Q2h,5) -
		           128*Q2e*S*pow(Q2h,5) + 232*m2*U*pow(Q2h,5) -
		           128*Q2e*U*pow(Q2h,5) - 164*m2*pow(Q2h,6) + 64*Q2e*pow(Q2h,6) +
		           64*Q2e*Q2h*pow(m,6)*pow(S,2) +
		           192*m4*Q2e*pow(Q2h,2)*pow(S,2) -
		           64*pow(m,6)*pow(Q2h,2)*pow(S,2) -
		           128*m4*pow(Q2h,3)*pow(S,2) + 260*m2*Q2e*pow(Q2h,3)*pow(S,2) -
		           132*m2*pow(Q2h,4)*pow(S,2) + 64*Q2e*pow(Q2h,4)*pow(S,2) +
		           112*m4*Q2e*pow(Q2h,2)*pow(U,2) - 48*m4*pow(Q2h,3)*pow(U,2) +
		           196*m2*Q2e*pow(Q2h,3)*pow(U,2) - 68*m2*pow(Q2h,4)*pow(U,2) +
		           64*Q2e*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-2)*(-64*m4*Q2e*S*U*pow(Q2h,3) + 64*m4*Q2e*S*pow(Q2h,4) +
		           64*m4*Q2e*U*pow(Q2h,4) + 64*m4*S*U*pow(Q2h,4) -
		           64*m2*Q2e*S*U*pow(Q2h,4) - 32*m4*Q2e*pow(Q2h,5) -
		           64*m4*S*pow(Q2h,5) + 64*m2*Q2e*S*pow(Q2h,5) -
		           64*m4*U*pow(Q2h,5) + 64*m2*Q2e*U*pow(Q2h,5) +
		           64*m2*S*U*pow(Q2h,5) + 32*m4*pow(Q2h,6) -
		           32*m2*Q2e*pow(Q2h,6) - 64*m2*S*pow(Q2h,6) -
		           64*m2*U*pow(Q2h,6) + 32*m2*pow(Q2h,7) -
		           32*m4*Q2e*pow(Q2h,3)*pow(S,2) + 32*m4*pow(Q2h,4)*pow(S,2) -
		           32*m2*Q2e*pow(Q2h,4)*pow(S,2) + 32*m2*pow(Q2h,5)*pow(S,2) -
		           32*m4*Q2e*pow(Q2h,3)*pow(U,2) + 32*m4*pow(Q2h,4)*pow(U,2) -
		           32*m2*Q2e*pow(Q2h,4)*pow(U,2) + 32*m2*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1)) +
		     SI.B0_M0m(2*l1k + m2 + Q2e - Q2h,m2)*
		      (14*l1k - 8*M2 + 7*Q2e - 2*Q2h - 24*S - 10*U +
		        pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (-16*l1k*Q2h - 10*Q2e*Q2h + 12*l1k*S + 6*Q2e*S + 18*Q2h*S +
		           12*l1k*U + 12*Q2e*U - 8*Q2h*U - 12*S*U - 12*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*m2*M2 - 16*M2*Q2e + 4*m2*Q2h - 16*M2*Q2h +
		           2*Q2e*Q2h + 12*m2*S - 2*Q2e*S - 34*Q2h*S - 12*m2*U -
		           22*Q2h*U + 4*S*U + 18*pow(Q2h,2) + 16*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (-16*m2*M2*Q2e - 16*m4*Q2h - 16*m2*M2*Q2h -
		           12*m2*Q2e*Q2h + 32*m4*S + 12*m2*Q2e*S + 36*m2*Q2h*S -
		           8*m4*U + 20*m2*Q2h*U - 24*m2*S*U - 20*m2*pow(Q2h,2) -
		           16*m2*pow(S,2)) + pow(2*l1k + Q2e - Q2h,-3)*
		         (-8*m4*Q2e*Q2h + 8*m4*Q2e*S + 32*m4*Q2h*S + 24*m4*Q2h*U -
		           16*m4*S*U - 16*m4*pow(Q2h,2) - 16*m4*pow(S,2)) +
		        pow(l1k,-1)*pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (Q2e*Q2h*S - 4*Q2e*Q2h*U - 6*Q2e*S*U + 8*Q2h*S*U -
		           Q2h*pow(Q2e,2) + 3*U*pow(Q2e,2) + 2*Q2e*pow(Q2h,2) -
		           4*S*pow(Q2h,2) + 2*Q2h*pow(S,2)) +
		        pow(l1k,-1)*(-6*m2*Q2h + 4*M2*Q2h + Q2e*Q2h - 7*Q2e*S +
		           2*Q2h*S + 12*m2*U - 6*Q2e*U + 15*Q2h*U + 12*S*U +
		           2*pow(Q2e,2) - 12*pow(Q2h,2) + 6*pow(S,2) - 2*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*m4*Q2h + 24*m2*M2*Q2h + 12*m2*Q2h*S + 20*m4*U +
		           30*m2*Q2h*U + 8*m2*S*U - 18*Q2h*S*U - 18*m2*pow(Q2h,2) +
		           8*M2*pow(Q2h,2) + 37*S*pow(Q2h,2) + 27*U*pow(Q2h,2) -
		           20*pow(Q2h,3) - 16*Q2h*pow(S,2) - 12*m2*pow(U,2) -
		           10*Q2h*pow(U,2)) + pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-2)*
		         (-52*m4*Q2h*S - 12*m4*Q2h*U + 24*m4*S*U + 20*m2*Q2h*S*U +
		           20*m4*pow(Q2h,2) - 26*m2*S*pow(Q2h,2) - 14*m2*U*pow(Q2h,2) +
		           10*m2*pow(Q2h,3) + 32*m4*pow(S,2) + 16*m2*Q2h*pow(S,2) -
		           8*m4*pow(U,2) + 4*m2*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-3)*
		         (16*m4*Q2h*S*U - 16*m4*S*pow(Q2h,2) - 16*m4*U*pow(Q2h,2) +
		           8*m4*pow(Q2h,3) + 8*m4*Q2h*pow(S,2) + 8*m4*Q2h*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*m4*Q2h*S - 32*m4*Q2h*U + 32*m4*S*U + 32*m2*Q2h*S*U +
		           16*m4*pow(Q2h,2) - 40*m2*S*pow(Q2h,2) - 40*m2*U*pow(Q2h,2) +
		           12*S*U*pow(Q2h,2) + 24*m2*pow(Q2h,3) - 16*S*pow(Q2h,3) -
		           12*U*pow(Q2h,3) + 8*pow(Q2h,4) + 16*m4*pow(S,2) +
		           16*m2*Q2h*pow(S,2) + 8*pow(Q2h,2)*pow(S,2) + 16*m4*pow(U,2) +
		           16*m2*Q2h*pow(U,2) + 4*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-3)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*m2*S*U*pow(Q2h,2) + 8*m2*S*pow(Q2h,3) + 8*m2*U*pow(Q2h,3) -
		           2*S*U*pow(Q2h,3) - 4*m2*pow(Q2h,4) + 2*S*pow(Q2h,4) +
		           2*U*pow(Q2h,4) - pow(Q2h,5) - 4*m2*pow(Q2h,2)*pow(S,2) -
		           pow(Q2h,3)*pow(S,2) - 4*m2*pow(Q2h,2)*pow(U,2) -
		           pow(Q2h,3)*pow(U,2)) +
		        (-288*Q2e*Q2h*S*pow(l1k,2) + 192*Q2e*Q2h*pow(l1k,3) +
		           384*Q2e*S*pow(l1k,3) - 192*Q2h*S*pow(l1k,3) -
		           240*Q2e*pow(l1k,4) + 96*Q2h*pow(l1k,4) + 192*S*pow(l1k,4) -
		           96*pow(l1k,5) - 144*l1k*Q2h*S*pow(Q2e,2) +
		           144*Q2h*pow(l1k,2)*pow(Q2e,2) + 288*S*pow(l1k,2)*pow(Q2e,2) -
		           240*pow(l1k,3)*pow(Q2e,2) + 48*l1k*Q2h*pow(Q2e,3) +
		           96*l1k*S*pow(Q2e,3) - 24*Q2h*S*pow(Q2e,3) -
		           120*pow(l1k,2)*pow(Q2e,3) - 30*l1k*pow(Q2e,4) +
		           6*Q2h*pow(Q2e,4) + 12*S*pow(Q2e,4) - 3*pow(Q2e,5) +
		           96*l1k*Q2e*Q2h*pow(S,2) - 144*Q2e*pow(l1k,2)*pow(S,2) +
		           96*Q2h*pow(l1k,2)*pow(S,2) - 96*pow(l1k,3)*pow(S,2) -
		           72*l1k*pow(Q2e,2)*pow(S,2) + 24*Q2h*pow(Q2e,2)*pow(S,2) -
		           12*pow(Q2e,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        (-32*l1k*m2*M2 - 32*l1k*M2*Q2e - 16*m2*M2*Q2e -
		           116*l1k*Q2e*Q2h - 16*l1k*m2*S + 76*l1k*Q2e*S - 8*m2*Q2e*S -
		           408*l1k*Q2h*S + 204*Q2e*Q2h*S - 16*l1k*m2*U + 8*l1k*Q2e*U -
		           16*m2*Q2e*U - 8*l1k*Q2h*U - 20*Q2e*Q2h*U + 16*m2*S*U -
		           8*Q2e*S*U + 8*Q2h*S*U - 32*M2*pow(l1k,2) -
		           64*Q2e*pow(l1k,2) + 272*Q2h*pow(l1k,2) + 136*S*pow(l1k,2) -
		           96*pow(l1k,3) - 32*l1k*pow(Q2e,2) - 8*M2*pow(Q2e,2) +
		           104*Q2h*pow(Q2e,2) + 34*S*pow(Q2e,2) + 8*U*pow(Q2e,2) -
		           12*pow(Q2e,3) + 252*l1k*pow(Q2h,2) - 284*Q2e*pow(Q2h,2) -
		           472*S*pow(Q2h,2) - 4*U*pow(Q2h,2) + 338*pow(Q2h,3) -
		           32*l1k*pow(S,2) + 16*m2*pow(S,2) - 20*Q2e*pow(S,2) +
		           144*Q2h*pow(S,2))*pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2e,2),-1) + pow(2*l1k + Q2e - Q2h,-1)*
		         (-12*Q2e*S*pow(Q2h,2) + 8*S*U*pow(Q2h,2) + 4*Q2e*pow(Q2h,3) -
		           4*U*pow(Q2h,3) + 8*Q2e*Q2h*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*(8*m2*Q2e*S*U + 16*Q2e*Q2h*S*U - 6*Q2h*S*pow(Q2e,2) -
		           4*m2*U*pow(Q2e,2) - 8*Q2h*U*pow(Q2e,2) - 4*S*U*pow(Q2e,2) +
		           Q2h*pow(Q2e,3) + 7*S*pow(Q2e,3) + 2*U*pow(Q2e,3) -
		           2*pow(Q2e,4) - 20*Q2e*S*pow(Q2h,2) + 7*pow(Q2e,2)*pow(Q2h,2) -
		           Q2e*pow(Q2h,3) + 4*S*pow(Q2h,3) - pow(Q2h,4) +
		           8*Q2e*Q2h*pow(S,2) - 6*pow(Q2e,2)*pow(S,2) +
		           12*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + m2 + Q2e - Q2h,-1)*
		         (32*l1k*S*pow(Q2h,2) + 16*Q2e*S*pow(Q2h,2) - 16*l1k*pow(Q2h,3) -
		           8*Q2e*pow(Q2h,3) - 32*S*pow(Q2h,3) + 16*pow(Q2h,4) -
		           16*l1k*Q2h*pow(S,2) - 8*Q2e*Q2h*pow(S,2) +
		           16*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (4*S*pow(Q2h,4) - pow(Q2h,5) - 4*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (-16*l1k*m2*Q2e*Q2h + 16*m4*Q2e*Q2h + 32*l1k*m2*Q2e*S +
		           32*l1k*m2*Q2h*S + 216*l1k*Q2e*Q2h*S - 48*m2*Q2e*Q2h*S -
		           32*l1k*Q2e*Q2h*U - 32*m2*Q2e*Q2h*U - 128*l1k*m2*S*U -
		           64*m4*S*U - 48*m2*Q2e*S*U + 64*m2*Q2h*S*U + 16*Q2e*Q2h*S*U -
		           16*m2*Q2e*pow(l1k,2) - 16*m2*Q2h*pow(l1k,2) -
		           120*Q2e*Q2h*pow(l1k,2) - 16*Q2e*S*pow(l1k,2) +
		           80*Q2h*S*pow(l1k,2) - 48*Q2h*pow(l1k,3) - 32*S*pow(l1k,3) +
		           16*pow(l1k,4) + 124*l1k*Q2h*pow(Q2e,2) - 72*l1k*S*pow(Q2e,2) -
		           220*Q2h*S*pow(Q2e,2) - 16*Q2h*U*pow(Q2e,2) +
		           44*pow(l1k,2)*pow(Q2e,2) - 36*l1k*pow(Q2e,3) -
		           130*Q2h*pow(Q2e,3) + 60*S*pow(Q2e,3) + 25*pow(Q2e,4) +
		           32*l1k*m2*pow(Q2h,2) + 16*m4*pow(Q2h,2) -
		           144*l1k*Q2e*pow(Q2h,2) + 68*m2*Q2e*pow(Q2h,2) -
		           96*l1k*S*pow(Q2h,2) - 48*m2*S*pow(Q2h,2) +
		           224*Q2e*S*pow(Q2h,2) - 16*Q2e*U*pow(Q2h,2) +
		           68*pow(l1k,2)*pow(Q2h,2) + 275*pow(Q2e,2)*pow(Q2h,2) +
		           60*l1k*pow(Q2h,3) - 4*m2*pow(Q2h,3) - 219*Q2e*pow(Q2h,3) -
		           112*S*pow(Q2h,3) + 81*pow(Q2h,4) - 64*l1k*m2*pow(S,2) -
		           32*m4*pow(S,2) - 8*m2*Q2e*pow(S,2) - 32*l1k*Q2h*pow(S,2) +
		           48*m2*Q2h*pow(S,2) - 44*Q2e*Q2h*pow(S,2) +
		           16*pow(l1k,2)*pow(S,2) + 28*pow(Q2e,2)*pow(S,2) +
		           32*pow(Q2h,2)*pow(S,2) - 64*l1k*m2*pow(U,2) -
		           32*m4*pow(U,2) + 16*l1k*Q2e*pow(U,2) - 24*m2*Q2e*pow(U,2) +
		           8*pow(Q2e,2)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*(64*m4*Q2e*Q2h*S + 32*m4*Q2e*Q2h*U -
		           32*m4*Q2e*S*U - 64*m4*Q2h*S*U - 56*m2*Q2e*Q2h*S*U -
		           32*m4*Q2e*pow(Q2h,2) + 32*m4*S*pow(Q2h,2) +
		           80*m2*Q2e*S*pow(Q2h,2) + 64*m4*U*pow(Q2h,2) +
		           48*m2*Q2e*U*pow(Q2h,2) + 8*m2*S*U*pow(Q2h,2) -
		           16*Q2e*S*U*pow(Q2h,2) - 8*U*pow(Q2e,2)*pow(Q2h,2) -
		           16*m4*pow(Q2h,3) - 36*m2*Q2e*pow(Q2h,3) -
		           8*m2*S*pow(Q2h,3) + 24*Q2e*S*pow(Q2h,3) + 8*m2*U*pow(Q2h,3) +
		           24*Q2e*U*pow(Q2h,3) + 4*pow(Q2e,2)*pow(Q2h,3) -
		           4*m2*pow(Q2h,4) - 16*Q2e*pow(Q2h,4) - 32*m4*Q2e*pow(S,2) -
		           16*m4*Q2h*pow(S,2) - 40*m2*Q2e*Q2h*pow(S,2) +
		           8*m2*pow(Q2h,2)*pow(S,2) - 16*Q2e*pow(Q2h,2)*pow(S,2) +
		           4*pow(Q2h,3)*pow(S,2) - 48*m4*Q2h*pow(U,2) -
		           16*m2*Q2e*Q2h*pow(U,2) + 4*Q2h*pow(Q2e,2)*pow(U,2) -
		           8*Q2e*pow(Q2h,2)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-2)*(32*m4*Q2e*Q2h*S*U - 32*Q2e*Q2h*S*pow(m,6) -
		           32*Q2e*Q2h*U*pow(m,6) + 32*Q2e*S*U*pow(m,6) -
		           32*Q2h*S*U*pow(m,6) - 40*m4*Q2e*S*pow(Q2h,2) -
		           40*m4*Q2e*U*pow(Q2h,2) - 16*m4*S*U*pow(Q2h,2) +
		           28*m2*Q2e*S*U*pow(Q2h,2) + 16*Q2e*pow(m,6)*pow(Q2h,2) +
		           32*S*pow(m,6)*pow(Q2h,2) + 32*U*pow(m,6)*pow(Q2h,2) +
		           24*m4*Q2e*pow(Q2h,3) + 24*m4*S*pow(Q2h,3) -
		           32*m2*Q2e*S*pow(Q2h,3) + 24*m4*U*pow(Q2h,3) -
		           28*m2*Q2e*U*pow(Q2h,3) - 8*m2*S*U*pow(Q2h,3) +
		           4*Q2e*S*U*pow(Q2h,3) - 16*pow(m,6)*pow(Q2h,3) -
		           16*m4*pow(Q2h,4) + 16*m2*Q2e*pow(Q2h,4) +
		           12*m2*S*pow(Q2h,4) - 4*Q2e*S*pow(Q2h,4) + 8*m2*U*pow(Q2h,4) -
		           4*Q2e*U*pow(Q2h,4) - 6*m2*pow(Q2h,5) + 2*Q2e*pow(Q2h,5) +
		           16*m4*Q2e*Q2h*pow(S,2) + 16*Q2e*pow(m,6)*pow(S,2) -
		           16*Q2h*pow(m,6)*pow(S,2) - 8*m4*pow(Q2h,2)*pow(S,2) +
		           16*m2*Q2e*pow(Q2h,2)*pow(S,2) - 6*m2*pow(Q2h,3)*pow(S,2) +
		           2*Q2e*pow(Q2h,3)*pow(S,2) + 16*m4*Q2e*Q2h*pow(U,2) +
		           16*Q2e*pow(m,6)*pow(U,2) - 16*Q2h*pow(m,6)*pow(U,2) -
		           8*m4*pow(Q2h,2)*pow(U,2) + 12*m2*Q2e*pow(Q2h,2)*pow(U,2) -
		           2*m2*pow(Q2h,3)*pow(U,2) + 2*Q2e*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-3)*(-8*m4*Q2e*S*U*pow(Q2h,2) + 8*m4*Q2e*S*pow(Q2h,3) +
		           8*m4*Q2e*U*pow(Q2h,3) + 8*m4*S*U*pow(Q2h,3) -
		           2*m2*Q2e*S*U*pow(Q2h,3) - 4*m4*Q2e*pow(Q2h,4) -
		           8*m4*S*pow(Q2h,4) + 2*m2*Q2e*S*pow(Q2h,4) -
		           8*m4*U*pow(Q2h,4) + 2*m2*Q2e*U*pow(Q2h,4) +
		           2*m2*S*U*pow(Q2h,4) + 4*m4*pow(Q2h,5) - m2*Q2e*pow(Q2h,5) -
		           2*m2*S*pow(Q2h,5) - 2*m2*U*pow(Q2h,5) + m2*pow(Q2h,6) -
		           4*m4*Q2e*pow(Q2h,2)*pow(S,2) + 4*m4*pow(Q2h,3)*pow(S,2) -
		           m2*Q2e*pow(Q2h,3)*pow(S,2) + m2*pow(Q2h,4)*pow(S,2) -
		           4*m4*Q2e*pow(Q2h,2)*pow(U,2) + 4*m4*pow(Q2h,3)*pow(U,2) -
		           m2*Q2e*pow(Q2h,3)*pow(U,2) + m2*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        (-128*m2*Q2e*S*pow(l1k,3) + 736*Q2e*Q2h*S*pow(l1k,3) +
		           64*m2*Q2e*pow(l1k,4) - 480*Q2e*Q2h*pow(l1k,4) -
		           320*Q2e*S*pow(l1k,4) - 320*Q2h*S*pow(l1k,4) +
		           192*Q2e*pow(l1k,5) + 192*Q2h*pow(l1k,5) + 128*S*pow(l1k,5) -
		           64*pow(l1k,6) + 248*l1k*Q2h*S*pow(Q2e,3) +
		           196*l1k*Q2h*pow(Q2e,4) - 168*l1k*S*pow(Q2e,4) +
		           220*Q2h*S*pow(Q2e,4) - 64*l1k*pow(Q2e,5) + 130*Q2h*pow(Q2e,5) -
		           60*S*pow(Q2e,5) - 25*pow(Q2e,6) +
		           256*Q2e*S*pow(l1k,2)*pow(Q2h,2) - 80*Q2e*pow(l1k,3)*pow(Q2h,2) +
		           384*S*pow(l1k,3)*pow(Q2h,2) - 272*pow(l1k,4)*pow(Q2h,2) +
		           352*l1k*S*pow(Q2e,2)*pow(Q2h,2) - 92*l1k*pow(Q2e,3)*pow(Q2h,2) -
		           256*S*pow(Q2e,3)*pow(Q2h,2) - 267*pow(Q2e,4)*pow(Q2h,2) -
		           416*l1k*Q2e*S*pow(Q2h,3) - 116*Q2e*pow(l1k,2)*pow(Q2h,3) +
		           448*S*pow(l1k,2)*pow(Q2h,3) - 240*pow(l1k,3)*pow(Q2h,3) -
		           368*l1k*pow(Q2e,2)*pow(Q2h,3) + 96*S*pow(Q2e,2)*pow(Q2h,3) +
		           243*pow(Q2e,3)*pow(Q2h,3) + 324*l1k*Q2e*pow(Q2h,4) -
		           324*pow(l1k,2)*pow(Q2h,4) - 81*pow(Q2e,2)*pow(Q2h,4) +
		           64*m2*Q2e*pow(l1k,2)*pow(S,2) -
		           272*Q2e*Q2h*pow(l1k,2)*pow(S,2) + 128*Q2e*pow(l1k,3)*pow(S,2) +
		           128*Q2h*pow(l1k,3)*pow(S,2) - 64*pow(l1k,4)*pow(S,2) -
		           16*l1k*Q2h*pow(Q2e,2)*pow(S,2) - 96*l1k*pow(Q2e,3)*pow(S,2) +
		           60*Q2h*pow(Q2e,3)*pow(S,2) - 28*pow(Q2e,4)*pow(S,2) +
		           112*l1k*Q2e*pow(Q2h,2)*pow(S,2) -
		           128*pow(l1k,2)*pow(Q2h,2)*pow(S,2) -
		           40*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) - 16*l1k*pow(Q2h,3)*pow(S,2) +
		           8*Q2e*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1)*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1)) +
		     SI.C0_mMQm0m(2*l1k + m2 + Q2e - Q2h,Q2h,m2)*
		      (-296*l1k*M2 - 64*m2*M2 + 28*l1k*Q2e - 68*M2*Q2e - 10*l1k*Q2h -
		        24*m2*Q2h - 120*M2*Q2h + 5*Q2e*Q2h - 64*l1k*S + 8*m2*S -
		        34*Q2e*S - 122*Q2h*S - 4*l1k*U + 8*m2*U - 8*Q2e*U - 384*Q2h*U +
		        68*S*U + 44*pow(l1k,2) + 11*pow(Q2e,2) + 236*pow(Q2h,2) +
		        16*pow(S,2) + pow(2*l1k + Q2e - Q2h,-1)*
		         (16*m2*M2*Q2h - 8*M2*Q2e*Q2h - 16*m2*Q2h*S - 8*Q2e*Q2h*S +
		           16*m2*pow(Q2h,2) + 8*M2*pow(Q2h,2) - 16*S*pow(Q2h,2) +
		           8*pow(Q2h,3) - 16*m2*pow(S,2) + 8*Q2e*pow(S,2) + 8*Q2h*pow(S,2)) \
		+ pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (32*m4*Q2h*S - 32*m2*Q2h*S*U - 16*m4*pow(Q2h,2) +
		           88*m2*S*pow(Q2h,2) + 16*m2*U*pow(Q2h,2) - 8*S*U*pow(Q2h,2) -
		           32*m2*pow(Q2h,3) + 24*S*pow(Q2h,3) + 4*U*pow(Q2h,3) -
		           8*pow(Q2h,4) - 48*m2*Q2h*pow(S,2) - 16*pow(Q2h,2)*pow(S,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (32*m4*Q2h*S*U - 48*m4*S*pow(Q2h,2) - 16*m4*U*pow(Q2h,2) +
		           24*m2*S*U*pow(Q2h,2) + 16*m4*pow(Q2h,3) -
		           36*m2*S*pow(Q2h,3) - 12*m2*U*pow(Q2h,3) + 4*S*U*pow(Q2h,3) +
		           12*m2*pow(Q2h,4) - 6*S*pow(Q2h,4) - 2*U*pow(Q2h,4) +
		           2*pow(Q2h,5) + 32*m4*Q2h*pow(S,2) +
		           24*m2*pow(Q2h,2)*pow(S,2) + 4*pow(Q2h,3)*pow(S,2)) +
		        144*pow(U,2) + pow(l1k,-1)*
		         (-24*m2*M2*Q2e - 48*m2*M2*Q2h - 8*m2*Q2e*Q2h -
		           24*M2*Q2e*Q2h + 4*m2*Q2h*S - 13*Q2e*Q2h*S + 4*m2*Q2e*U -
		           56*m2*Q2h*U - 38*Q2e*Q2h*U - 8*m2*S*U + 6*Q2e*S*U +
		           44*Q2h*S*U - 12*M2*pow(Q2e,2) + Q2h*pow(Q2e,2) -
		           7*S*pow(Q2e,2) - 3*U*pow(Q2e,2) + 2*pow(Q2e,3) +
		           32*m2*pow(Q2h,2) - 24*M2*pow(Q2h,2) + 22*Q2e*pow(Q2h,2) -
		           74*S*pow(Q2h,2) - 158*U*pow(Q2h,2) + 98*pow(Q2h,3) +
		           6*Q2e*pow(S,2) + 10*Q2h*pow(S,2) + 16*m2*pow(U,2) +
		           16*Q2e*pow(U,2) + 60*Q2h*pow(U,2)) +
		        (-384*Q2e*Q2h*S*pow(l1k,3) + 240*Q2e*Q2h*pow(l1k,4) +
		           960*Q2e*S*pow(l1k,4) - 192*Q2h*S*pow(l1k,4) -
		           576*Q2e*pow(l1k,5) + 96*Q2h*pow(l1k,5) + 384*S*pow(l1k,5) -
		           192*pow(l1k,6) - 288*Q2h*S*pow(l1k,2)*pow(Q2e,2) +
		           240*Q2h*pow(l1k,3)*pow(Q2e,2) + 960*S*pow(l1k,3)*pow(Q2e,2) -
		           720*pow(l1k,4)*pow(Q2e,2) - 96*l1k*Q2h*S*pow(Q2e,3) +
		           120*Q2h*pow(l1k,2)*pow(Q2e,3) + 480*S*pow(l1k,2)*pow(Q2e,3) -
		           480*pow(l1k,3)*pow(Q2e,3) + 30*l1k*Q2h*pow(Q2e,4) +
		           120*l1k*S*pow(Q2e,4) - 12*Q2h*S*pow(Q2e,4) -
		           180*pow(l1k,2)*pow(Q2e,4) - 36*l1k*pow(Q2e,5) +
		           3*Q2h*pow(Q2e,5) + 12*S*pow(Q2e,5) - 3*pow(Q2e,6) +
		           144*Q2e*Q2h*pow(l1k,2)*pow(S,2) - 384*Q2e*pow(l1k,3)*pow(S,2) +
		           96*Q2h*pow(l1k,3)*pow(S,2) - 192*pow(l1k,4)*pow(S,2) +
		           72*l1k*Q2h*pow(Q2e,2)*pow(S,2) -
		           288*pow(l1k,2)*pow(Q2e,2)*pow(S,2) -
		           96*l1k*pow(Q2e,3)*pow(S,2) + 12*Q2h*pow(Q2e,3)*pow(S,2) -
		           12*pow(Q2e,4)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        (-64*l1k*m2*M2*Q2e - 32*l1k*m2*Q2e*S + 476*l1k*Q2e*Q2h*S -
		           48*l1k*m2*Q2e*U - 16*l1k*Q2e*Q2h*U + 32*l1k*m2*S*U -
		           40*l1k*Q2e*S*U + 32*m2*Q2e*S*U + 16*Q2e*Q2h*S*U -
		           64*m2*M2*pow(l1k,2) - 48*M2*Q2e*pow(l1k,2) -
		           272*Q2e*Q2h*pow(l1k,2) - 32*m2*S*pow(l1k,2) +
		           312*Q2e*S*pow(l1k,2) - 760*Q2h*S*pow(l1k,2) -
		           32*m2*U*pow(l1k,2) + 48*Q2e*U*pow(l1k,2) - 16*S*U*pow(l1k,2) -
		           32*M2*pow(l1k,3) - 224*Q2e*pow(l1k,3) + 528*Q2h*pow(l1k,3) +
		           416*S*pow(l1k,3) + 16*U*pow(l1k,3) - 256*pow(l1k,4) -
		           24*l1k*M2*pow(Q2e,2) - 16*m2*M2*pow(Q2e,2) +
		           168*l1k*Q2h*pow(Q2e,2) + 156*l1k*S*pow(Q2e,2) -
		           8*m2*S*pow(Q2e,2) - 322*Q2h*S*pow(Q2e,2) +
		           48*l1k*U*pow(Q2e,2) - 24*m2*U*pow(Q2e,2) -
		           16*Q2h*U*pow(Q2e,2) - 28*S*U*pow(Q2e,2) -
		           128*pow(l1k,2)*pow(Q2e,2) - 56*l1k*pow(Q2e,3) -
		           4*M2*pow(Q2e,3) - 140*Q2h*pow(Q2e,3) + 50*S*pow(Q2e,3) +
		           20*U*pow(Q2e,3) - 16*pow(Q2e,4) - 584*l1k*Q2e*pow(Q2h,2) -
		           896*l1k*S*pow(Q2h,2) + 1112*Q2e*S*pow(Q2h,2) +
		           480*pow(l1k,2)*pow(Q2h,2) + 672*pow(Q2e,2)*pow(Q2h,2) +
		           648*l1k*pow(Q2h,3) - 1352*Q2e*pow(Q2h,3) - 1296*S*pow(Q2h,3) +
		           1320*pow(Q2h,4) + 32*l1k*m2*pow(S,2) - 88*l1k*Q2e*pow(S,2) +
		           16*m2*Q2e*pow(S,2) + 256*l1k*Q2h*pow(S,2) -
		           188*Q2e*Q2h*pow(S,2) - 144*pow(l1k,2)*pow(S,2) -
		           36*pow(Q2e,2)*pow(S,2) + 384*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*(8*m2*S*U*pow(Q2e,2) + 8*Q2h*S*U*pow(Q2e,2) +
		           Q2h*S*pow(Q2e,3) - 4*m2*U*pow(Q2e,3) - 4*Q2h*U*pow(Q2e,3) -
		           6*S*U*pow(Q2e,3) - Q2h*pow(Q2e,4) + 7*S*pow(Q2e,4) +
		           3*U*pow(Q2e,4) - 2*pow(Q2e,5) - 12*S*pow(Q2e,2)*pow(Q2h,2) +
		           4*pow(Q2e,3)*pow(Q2h,2) + 2*Q2h*pow(Q2e,2)*pow(S,2) -
		           6*pow(Q2e,3)*pow(S,2) + 8*Q2e*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (2560*m4*Q2e*Q2h*S*U*pow(l1k,2) - 1024*m4*Q2e*Q2h*S*pow(l1k,3) -
		           1024*m4*Q2e*Q2h*U*pow(l1k,3) + 512*m4*Q2e*S*U*pow(l1k,3) +
		           512*m4*Q2h*S*U*pow(l1k,3) + 9344*m2*Q2e*Q2h*S*U*pow(l1k,3) -
		           4096*m2*Q2e*Q2h*S*pow(l1k,4) - 6144*m2*Q2e*Q2h*U*pow(l1k,4) +
		           2560*m2*Q2e*S*U*pow(l1k,4) + 16128*Q2e*Q2h*S*U*pow(l1k,4) -
		           14336*Q2e*Q2h*S*pow(l1k,5) - 26624*Q2e*Q2h*U*pow(l1k,5) +
		           9728*Q2e*S*U*pow(l1k,5) + 512*l1k*Q2e*Q2h*S*U*pow(m,6) +
		           12992*Q2h*S*U*pow(l1k,3)*pow(Q2e,2) -
		           11264*Q2h*S*pow(l1k,4)*pow(Q2e,2) -
		           24576*Q2h*U*pow(l1k,4)*pow(Q2e,2) +
		           8064*S*U*pow(l1k,4)*pow(Q2e,2) +
		           2304*Q2h*S*U*pow(l1k,2)*pow(Q2e,3) -
		           2624*Q2h*S*pow(l1k,3)*pow(Q2e,3) -
		           8256*Q2h*U*pow(l1k,3)*pow(Q2e,3) +
		           2112*S*U*pow(l1k,3)*pow(Q2e,3) + 272*l1k*Q2h*S*U*pow(Q2e,4) -
		           320*Q2h*S*pow(l1k,2)*pow(Q2e,4) -
		           1632*Q2h*U*pow(l1k,2)*pow(Q2e,4) +
		           288*S*U*pow(l1k,2)*pow(Q2e,4) - 16*l1k*Q2h*S*pow(Q2e,5) -
		           176*l1k*Q2h*U*pow(Q2e,5) + 16*l1k*S*U*pow(Q2e,5) +
		           16*Q2h*S*U*pow(Q2e,5) - 8*Q2h*U*pow(Q2e,6) +
		           2432*l1k*m4*Q2e*S*U*pow(Q2h,2) -
		           3072*m4*Q2e*S*pow(l1k,2)*pow(Q2h,2) -
		           4608*m4*Q2e*U*pow(l1k,2)*pow(Q2h,2) -
		           1024*m4*S*U*pow(l1k,2)*pow(Q2h,2) +
		           9472*m2*Q2e*S*U*pow(l1k,2)*pow(Q2h,2) +
		           1024*m4*Q2e*pow(l1k,3)*pow(Q2h,2) -
		           12800*m2*Q2e*S*pow(l1k,3)*pow(Q2h,2) -
		           20992*m2*Q2e*U*pow(l1k,3)*pow(Q2h,2) -
		           2432*m2*S*U*pow(l1k,3)*pow(Q2h,2) +
		           1280*Q2e*S*U*pow(l1k,3)*pow(Q2h,2) +
		           5120*m2*Q2e*pow(l1k,4)*pow(Q2h,2) -
		           19968*Q2e*S*pow(l1k,4)*pow(Q2h,2) -
		           32256*Q2e*U*pow(l1k,4)*pow(Q2h,2) +
		           20480*Q2e*pow(l1k,5)*pow(Q2h,2) -
		           1024*l1k*Q2e*S*pow(m,6)*pow(Q2h,2) -
		           1024*l1k*Q2e*U*pow(m,6)*pow(Q2h,2) +
		           512*l1k*S*U*pow(m,6)*pow(Q2h,2) +
		           512*Q2e*S*U*pow(m,6)*pow(Q2h,2) +
		           5280*S*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2) -
		           16512*S*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) -
		           25984*U*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) +
		           17920*pow(l1k,4)*pow(Q2e,2)*pow(Q2h,2) +
		           592*l1k*S*U*pow(Q2e,3)*pow(Q2h,2) -
		           2784*S*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2) -
		           4288*U*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2) +
		           5440*pow(l1k,3)*pow(Q2e,3)*pow(Q2h,2) -
		           304*l1k*S*pow(Q2e,4)*pow(Q2h,2) -
		           560*l1k*U*pow(Q2e,4)*pow(Q2h,2) +
		           56*S*U*pow(Q2e,4)*pow(Q2h,2) +
		           976*pow(l1k,2)*pow(Q2e,4)*pow(Q2h,2) +
		           96*l1k*pow(Q2e,5)*pow(Q2h,2) - 16*S*pow(Q2e,5)*pow(Q2h,2) -
		           64*U*pow(Q2e,5)*pow(Q2h,2) + 4*pow(Q2e,6)*pow(Q2h,2) -
		           3200*l1k*m4*Q2e*S*pow(Q2h,3) - 3968*l1k*m4*Q2e*U*pow(Q2h,3) -
		           1408*l1k*m4*S*U*pow(Q2h,3) + 4032*l1k*m2*Q2e*S*U*pow(Q2h,3) +
		           768*m4*Q2e*S*U*pow(Q2h,3) +
		           3072*m4*Q2e*pow(l1k,2)*pow(Q2h,3) +
		           1024*m4*S*pow(l1k,2)*pow(Q2h,3) -
		           11776*m2*Q2e*S*pow(l1k,2)*pow(Q2h,3) +
		           1536*m4*U*pow(l1k,2)*pow(Q2h,3) -
		           18176*m2*Q2e*U*pow(l1k,2)*pow(Q2h,3) -
		           5248*m2*S*U*pow(l1k,2)*pow(Q2h,3) -
		           3200*Q2e*S*U*pow(l1k,2)*pow(Q2h,3) +
		           14336*m2*Q2e*pow(l1k,3)*pow(Q2h,3) +
		           3584*m2*S*pow(l1k,3)*pow(Q2h,3) -
		           384*Q2e*S*pow(l1k,3)*pow(Q2h,3) +
		           6656*m2*U*pow(l1k,3)*pow(Q2h,3) +
		           2944*Q2e*U*pow(l1k,3)*pow(Q2h,3) +
		           18432*Q2e*pow(l1k,4)*pow(Q2h,3) +
		           768*l1k*Q2e*pow(m,6)*pow(Q2h,3) -
		           512*Q2e*S*pow(m,6)*pow(Q2h,3) -
		           512*Q2e*U*pow(m,6)*pow(Q2h,3) - 512*S*U*pow(m,6)*pow(Q2h,3) +
		           784*l1k*S*U*pow(Q2e,2)*pow(Q2h,3) -
		           6144*S*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3) -
		           7936*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3) +
		           15808*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,3) -
		           688*l1k*S*pow(Q2e,3)*pow(Q2h,3) -
		           784*l1k*U*pow(Q2e,3)*pow(Q2h,3) +
		           16*S*U*pow(Q2e,3)*pow(Q2h,3) +
		           2560*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,3) +
		           336*l1k*pow(Q2e,4)*pow(Q2h,3) - 64*S*pow(Q2e,4)*pow(Q2h,3) -
		           16*U*pow(Q2e,4)*pow(Q2h,3) + 36*pow(Q2e,5)*pow(Q2h,3) +
		           2432*l1k*m4*Q2e*pow(Q2h,4) + 1920*l1k*m4*S*pow(Q2h,4) -
		           4704*l1k*m2*Q2e*S*pow(Q2h,4) - 896*m4*Q2e*S*pow(Q2h,4) +
		           2688*l1k*m4*U*pow(Q2h,4) - 6432*l1k*m2*Q2e*U*pow(Q2h,4) -
		           1024*m4*Q2e*U*pow(Q2h,4) - 2944*l1k*m2*S*U*pow(Q2h,4) -
		           640*m4*S*U*pow(Q2h,4) - 1184*l1k*Q2e*S*U*pow(Q2h,4) +
		           736*m2*Q2e*S*U*pow(Q2h,4) - 1280*m4*pow(l1k,2)*pow(Q2h,4) +
		           11136*m2*Q2e*pow(l1k,2)*pow(Q2h,4) +
		           6784*m2*S*pow(l1k,2)*pow(Q2h,4) +
		           4032*Q2e*S*pow(l1k,2)*pow(Q2h,4) +
		           11392*m2*U*pow(l1k,2)*pow(Q2h,4) +
		           6208*Q2e*U*pow(l1k,2)*pow(Q2h,4) -
		           5120*m2*pow(l1k,3)*pow(Q2h,4) -
		           2752*Q2e*pow(l1k,3)*pow(Q2h,4) - 256*l1k*pow(m,6)*pow(Q2h,4) +
		           256*Q2e*pow(m,6)*pow(Q2h,4) + 512*S*pow(m,6)*pow(Q2h,4) +
		           512*U*pow(m,6)*pow(Q2h,4) - 800*l1k*S*pow(Q2e,2)*pow(Q2h,4) -
		           800*l1k*U*pow(Q2e,2)*pow(Q2h,4) +
		           4480*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,4) +
		           400*l1k*pow(Q2e,3)*pow(Q2h,4) - 8*S*pow(Q2e,3)*pow(Q2h,4) +
		           4*pow(Q2e,4)*pow(Q2h,4) - 1664*l1k*m4*pow(Q2h,5) +
		           3696*l1k*m2*Q2e*pow(Q2h,5) + 576*m4*Q2e*pow(Q2h,5) +
		           3488*l1k*m2*S*pow(Q2h,5) + 768*m4*S*pow(Q2h,5) +
		           1312*l1k*Q2e*S*pow(Q2h,5) - 800*m2*Q2e*S*pow(Q2h,5) +
		           4960*l1k*m2*U*pow(Q2h,5) + 896*m4*U*pow(Q2h,5) +
		           1696*l1k*Q2e*U*pow(Q2h,5) - 992*m2*Q2e*U*pow(Q2h,5) -
		           672*m2*S*U*pow(Q2h,5) - 80*Q2e*S*U*pow(Q2h,5) -
		           7168*m2*pow(l1k,2)*pow(Q2h,5) -
		           3808*Q2e*pow(l1k,2)*pow(Q2h,5) - 256*pow(m,6)*pow(Q2h,5) +
		           400*l1k*pow(Q2e,2)*pow(Q2h,5) - 2896*l1k*m2*pow(Q2h,6) -
		           512*m4*pow(Q2h,6) - 912*l1k*Q2e*pow(Q2h,6) +
		           528*m2*Q2e*pow(Q2h,6) + 736*m2*S*pow(Q2h,6) +
		           80*Q2e*S*pow(Q2h,6) + 928*m2*U*pow(Q2h,6) +
		           80*Q2e*U*pow(Q2h,6) - 496*m2*pow(Q2h,7) - 40*Q2e*pow(Q2h,7) +
		           768*m4*Q2e*Q2h*pow(l1k,2)*pow(S,2) +
		           256*m4*Q2e*pow(l1k,3)*pow(S,2) +
		           256*m4*Q2h*pow(l1k,3)*pow(S,2) +
		           2496*m2*Q2e*Q2h*pow(l1k,3)*pow(S,2) +
		           768*m2*Q2e*pow(l1k,4)*pow(S,2) +
		           3968*Q2e*Q2h*pow(l1k,4)*pow(S,2) +
		           2304*Q2e*pow(l1k,5)*pow(S,2) +
		           256*l1k*Q2e*Q2h*pow(m,6)*pow(S,2) +
		           3072*Q2h*pow(l1k,3)*pow(Q2e,2)*pow(S,2) +
		           1600*pow(l1k,4)*pow(Q2e,2)*pow(S,2) +
		           400*Q2h*pow(l1k,2)*pow(Q2e,3)*pow(S,2) +
		           256*pow(l1k,3)*pow(Q2e,3)*pow(S,2) +
		           24*l1k*Q2h*pow(Q2e,4)*pow(S,2) +
		           16*pow(l1k,2)*pow(Q2e,4)*pow(S,2) +
		           832*l1k*m4*Q2e*pow(Q2h,2)*pow(S,2) -
		           512*m4*pow(l1k,2)*pow(Q2h,2)*pow(S,2) +
		           2624*m2*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(S,2) -
		           576*m2*pow(l1k,3)*pow(Q2h,2)*pow(S,2) +
		           1024*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(S,2) +
		           256*l1k*pow(m,6)*pow(Q2h,2)*pow(S,2) +
		           256*Q2e*pow(m,6)*pow(Q2h,2)*pow(S,2) +
		           1824*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) +
		           200*l1k*pow(Q2e,3)*pow(Q2h,2)*pow(S,2) +
		           12*pow(Q2e,4)*pow(Q2h,2)*pow(S,2) -
		           320*l1k*m4*pow(Q2h,3)*pow(S,2) +
		           1312*l1k*m2*Q2e*pow(Q2h,3)*pow(S,2) +
		           320*m4*Q2e*pow(Q2h,3)*pow(S,2) -
		           1280*m2*pow(l1k,2)*pow(Q2h,3)*pow(S,2) -
		           832*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(S,2) -
		           256*pow(m,6)*pow(Q2h,3)*pow(S,2) +
		           352*l1k*pow(Q2e,2)*pow(Q2h,3)*pow(S,2) +
		           28*pow(Q2e,3)*pow(Q2h,3)*pow(S,2) -
		           896*l1k*m2*pow(Q2h,4)*pow(S,2) - 256*m4*pow(Q2h,4)*pow(S,2) -
		           400*l1k*Q2e*pow(Q2h,4)*pow(S,2) +
		           272*m2*Q2e*pow(Q2h,4)*pow(S,2) +
		           4*pow(Q2e,2)*pow(Q2h,4)*pow(S,2) -
		           240*m2*pow(Q2h,5)*pow(S,2) - 40*Q2e*pow(Q2h,5)*pow(S,2) +
		           1792*m4*Q2e*Q2h*pow(l1k,2)*pow(U,2) +
		           256*m4*Q2e*pow(l1k,3)*pow(U,2) +
		           256*m4*Q2h*pow(l1k,3)*pow(U,2) +
		           7616*m2*Q2e*Q2h*pow(l1k,3)*pow(U,2) +
		           1792*m2*Q2e*pow(l1k,4)*pow(U,2) +
		           13696*Q2e*Q2h*pow(l1k,4)*pow(U,2) +
		           8448*Q2e*pow(l1k,5)*pow(U,2) +
		           256*l1k*Q2e*Q2h*pow(m,6)*pow(U,2) +
		           10624*Q2h*pow(l1k,3)*pow(Q2e,2)*pow(U,2) +
		           8256*pow(l1k,4)*pow(Q2e,2)*pow(U,2) +
		           1808*Q2h*pow(l1k,2)*pow(Q2e,3)*pow(U,2) +
		           3072*pow(l1k,3)*pow(Q2e,3)*pow(U,2) +
		           232*l1k*Q2h*pow(Q2e,4)*pow(U,2) +
		           672*pow(l1k,2)*pow(Q2e,4)*pow(U,2) +
		           80*l1k*pow(Q2e,5)*pow(U,2) + 28*Q2h*pow(Q2e,5)*pow(U,2) +
		           4*pow(Q2e,6)*pow(U,2) + 1600*l1k*m4*Q2e*pow(Q2h,2)*pow(U,2) -
		           512*m4*pow(l1k,2)*pow(Q2h,2)*pow(U,2) +
		           7360*m2*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(U,2) -
		           2112*m2*pow(l1k,3)*pow(Q2h,2)*pow(U,2) -
		           512*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(U,2) +
		           256*l1k*pow(m,6)*pow(Q2h,2)*pow(U,2) +
		           256*Q2e*pow(m,6)*pow(Q2h,2)*pow(U,2) +
		           3472*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) +
		           376*l1k*pow(Q2e,3)*pow(Q2h,2)*pow(U,2) +
		           12*pow(Q2e,4)*pow(Q2h,2)*pow(U,2) -
		           1088*l1k*m4*pow(Q2h,3)*pow(U,2) +
		           2784*l1k*m2*Q2e*pow(Q2h,3)*pow(U,2) +
		           448*m4*Q2e*pow(Q2h,3)*pow(U,2) -
		           4480*m2*pow(l1k,2)*pow(Q2h,3)*pow(U,2) -
		           2496*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(U,2) -
		           256*pow(m,6)*pow(Q2h,3)*pow(U,2) +
		           400*l1k*pow(Q2e,2)*pow(Q2h,3)*pow(U,2) -
		           2112*l1k*m2*pow(Q2h,4)*pow(U,2) -
		           384*m4*pow(Q2h,4)*pow(U,2) - 784*l1k*Q2e*pow(Q2h,4)*pow(U,2) +
		           464*m2*Q2e*pow(Q2h,4)*pow(U,2) - 432*m2*pow(Q2h,5)*pow(U,2) -
		           40*Q2e*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        pow(l1k,-1)*(4*S*U*pow(Q2e,5)*pow(Q2h,2) -
		           4*U*pow(Q2e,6)*pow(Q2h,2) + 128*Q2e*S*U*pow(m,6)*pow(Q2h,3) -
		           4*S*U*pow(Q2e,4)*pow(Q2h,3) - 4*S*pow(Q2e,5)*pow(Q2h,3) +
		           4*U*pow(Q2e,5)*pow(Q2h,3) + 2*pow(Q2e,6)*pow(Q2h,3) +
		           128*m4*Q2e*S*U*pow(Q2h,4) - 128*Q2e*S*pow(m,6)*pow(Q2h,4) -
		           128*Q2e*U*pow(m,6)*pow(Q2h,4) - 128*S*U*pow(m,6)*pow(Q2h,4) +
		           4*S*pow(Q2e,4)*pow(Q2h,4) - 2*pow(Q2e,5)*pow(Q2h,4) -
		           128*m4*Q2e*S*pow(Q2h,5) - 128*m4*Q2e*U*pow(Q2h,5) -
		           128*m4*S*U*pow(Q2h,5) + 40*m2*Q2e*S*U*pow(Q2h,5) +
		           64*Q2e*pow(m,6)*pow(Q2h,5) + 128*S*pow(m,6)*pow(Q2h,5) +
		           128*U*pow(m,6)*pow(Q2h,5) + 64*m4*Q2e*pow(Q2h,6) +
		           128*m4*S*pow(Q2h,6) - 40*m2*Q2e*S*pow(Q2h,6) +
		           128*m4*U*pow(Q2h,6) - 40*m2*Q2e*U*pow(Q2h,6) -
		           40*m2*S*U*pow(Q2h,6) - 64*pow(m,6)*pow(Q2h,6) -
		           64*m4*pow(Q2h,7) + 20*m2*Q2e*pow(Q2h,7) +
		           40*m2*S*pow(Q2h,7) + 40*m2*U*pow(Q2h,7) - 20*m2*pow(Q2h,8) +
		           64*Q2e*pow(m,6)*pow(Q2h,3)*pow(S,2) +
		           2*pow(Q2e,4)*pow(Q2h,3)*pow(S,2) +
		           64*m4*Q2e*pow(Q2h,4)*pow(S,2) -
		           64*pow(m,6)*pow(Q2h,4)*pow(S,2) -
		           2*pow(Q2e,3)*pow(Q2h,4)*pow(S,2) - 64*m4*pow(Q2h,5)*pow(S,2) +
		           20*m2*Q2e*pow(Q2h,5)*pow(S,2) - 20*m2*pow(Q2h,6)*pow(S,2) +
		           2*Q2h*pow(Q2e,6)*pow(U,2) - 2*pow(Q2e,5)*pow(Q2h,2)*pow(U,2) +
		           64*Q2e*pow(m,6)*pow(Q2h,3)*pow(U,2) +
		           64*m4*Q2e*pow(Q2h,4)*pow(U,2) -
		           64*pow(m,6)*pow(Q2h,4)*pow(U,2) - 64*m4*pow(Q2h,5)*pow(U,2) +
		           20*m2*Q2e*pow(Q2h,5)*pow(U,2) - 20*m2*pow(Q2h,6)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        pow(l1k,-2)*(32*Q2e*Q2h*S*U*pow(m,6) +
		           24*m4*Q2e*S*U*pow(Q2h,2) - 48*Q2e*S*pow(m,6)*pow(Q2h,2) -
		           16*Q2e*U*pow(m,6)*pow(Q2h,2) - 32*S*U*pow(m,6)*pow(Q2h,2) -
		           36*m4*Q2e*S*pow(Q2h,3) - 12*m4*Q2e*U*pow(Q2h,3) -
		           24*m4*S*U*pow(Q2h,3) + 4*m2*Q2e*S*U*pow(Q2h,3) +
		           16*Q2e*pow(m,6)*pow(Q2h,3) + 48*S*pow(m,6)*pow(Q2h,3) +
		           16*U*pow(m,6)*pow(Q2h,3) + 12*m4*Q2e*pow(Q2h,4) +
		           36*m4*S*pow(Q2h,4) - 6*m2*Q2e*S*pow(Q2h,4) +
		           12*m4*U*pow(Q2h,4) - 2*m2*Q2e*U*pow(Q2h,4) -
		           4*m2*S*U*pow(Q2h,4) - 16*pow(m,6)*pow(Q2h,4) -
		           12*m4*pow(Q2h,5) + 2*m2*Q2e*pow(Q2h,5) + 6*m2*S*pow(Q2h,5) +
		           2*m2*U*pow(Q2h,5) - 2*m2*pow(Q2h,6) +
		           32*Q2e*Q2h*pow(m,6)*pow(S,2) +
		           24*m4*Q2e*pow(Q2h,2)*pow(S,2) -
		           32*pow(m,6)*pow(Q2h,2)*pow(S,2) - 24*m4*pow(Q2h,3)*pow(S,2) +
		           4*m2*Q2e*pow(Q2h,3)*pow(S,2) - 4*m2*pow(Q2h,4)*pow(S,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        (128*l1k*m4*M2*Q2e + 128*l1k*m4*M2*Q2h +
		           32*l1k*m4*Q2e*Q2h + 480*l1k*m2*M2*Q2e*Q2h +
		           128*m4*M2*Q2e*Q2h + 192*l1k*m4*Q2h*S +
		           616*l1k*m2*Q2e*Q2h*S + 16*m4*Q2e*Q2h*S + 192*l1k*m4*Q2h*U +
		           1112*l1k*m2*Q2e*Q2h*U + 176*m4*Q2e*Q2h*U - 128*l1k*m4*S*U -
		           320*l1k*m2*Q2e*S*U - 64*m4*Q2e*S*U - 1472*l1k*m2*Q2h*S*U -
		           256*m4*Q2h*S*U - 1088*l1k*Q2e*Q2h*S*U - 408*m2*Q2e*Q2h*S*U +
		           384*m2*M2*Q2e*pow(l1k,2) - 32*m2*Q2e*Q2h*pow(l1k,2) -
		           64*M2*Q2e*Q2h*pow(l1k,2) + 64*m2*Q2e*S*pow(l1k,2) +
		           1056*m2*Q2h*S*pow(l1k,2) + 1888*Q2e*Q2h*S*pow(l1k,2) +
		           1504*m2*Q2h*U*pow(l1k,2) + 4352*Q2e*Q2h*U*pow(l1k,2) -
		           640*m2*S*U*pow(l1k,2) - 1056*Q2e*S*U*pow(l1k,2) -
		           5248*Q2h*S*U*pow(l1k,2) - 32*m2*Q2e*pow(l1k,3) +
		           1152*M2*Q2e*pow(l1k,3) - 32*m2*Q2h*pow(l1k,3) -
		           240*Q2e*Q2h*pow(l1k,3) - 32*Q2e*S*pow(l1k,3) +
		           3744*Q2h*S*pow(l1k,3) + 6656*Q2h*U*pow(l1k,3) -
		           2432*S*U*pow(l1k,3) - 96*Q2h*pow(l1k,4) - 64*S*pow(l1k,4) +
		           32*pow(l1k,5) + 224*l1k*M2*Q2h*pow(Q2e,2) -
		           24*l1k*Q2h*S*pow(Q2e,2) + 1584*l1k*Q2h*U*pow(Q2e,2) -
		           256*l1k*S*U*pow(Q2e,2) - 224*Q2h*S*U*pow(Q2e,2) +
		           800*M2*pow(l1k,2)*pow(Q2e,2) +
		           248*Q2h*pow(l1k,2)*pow(Q2e,2) - 144*S*pow(l1k,2)*pow(Q2e,2) +
		           88*pow(l1k,3)*pow(Q2e,2) + 128*l1k*M2*pow(Q2e,3) -
		           260*l1k*Q2h*pow(Q2e,3) + 40*M2*Q2h*pow(Q2e,3) +
		           120*l1k*S*pow(Q2e,3) + 480*Q2h*S*pow(Q2e,3) +
		           164*Q2h*U*pow(Q2e,3) - 8*S*U*pow(Q2e,3) -
		           72*pow(l1k,2)*pow(Q2e,3) + 50*l1k*pow(Q2e,4) +
		           8*M2*pow(Q2e,4) + 226*Q2h*pow(Q2e,4) - 84*S*pow(Q2e,4) -
		           32*pow(Q2e,5) - 224*l1k*m4*pow(Q2h,2) -
		           288*l1k*m2*M2*pow(Q2h,2) - 128*m4*M2*pow(Q2h,2) -
		           792*l1k*m2*Q2e*pow(Q2h,2) - 80*m4*Q2e*pow(Q2h,2) -
		           160*l1k*M2*Q2e*pow(Q2h,2) + 160*m2*M2*Q2e*pow(Q2h,2) +
		           1848*l1k*m2*S*pow(Q2h,2) + 368*m4*S*pow(Q2h,2) +
		           1920*l1k*Q2e*S*pow(Q2h,2) + 416*m2*Q2e*S*pow(Q2h,2) +
		           2696*l1k*m2*U*pow(Q2h,2) + 464*m4*U*pow(Q2h,2) +
		           2072*l1k*Q2e*U*pow(Q2h,2) + 1184*m2*Q2e*U*pow(Q2h,2) -
		           2944*l1k*S*U*pow(Q2h,2) - 568*m2*S*U*pow(Q2h,2) -
		           248*Q2e*S*U*pow(Q2h,2) - 1216*m2*pow(l1k,2)*pow(Q2h,2) -
		           3184*Q2e*pow(l1k,2)*pow(Q2h,2) + 6592*S*pow(l1k,2)*pow(Q2h,2) +
		           11392*U*pow(l1k,2)*pow(Q2h,2) - 4984*pow(l1k,3)*pow(Q2h,2) -
		           458*l1k*pow(Q2e,2)*pow(Q2h,2) + 8*M2*pow(Q2e,2)*pow(Q2h,2) -
		           688*S*pow(Q2e,2)*pow(Q2h,2) + 536*U*pow(Q2e,2)*pow(Q2h,2) -
		           794*pow(Q2e,3)*pow(Q2h,2) - 1672*l1k*m2*pow(Q2h,3) -
		           272*m4*pow(Q2h,3) - 128*m2*M2*pow(Q2h,3) -
		           1734*l1k*Q2e*pow(Q2h,3) - 652*m2*Q2e*pow(Q2h,3) -
		           48*M2*Q2e*pow(Q2h,3) + 3264*l1k*S*pow(Q2h,3) +
		           640*m2*S*pow(Q2h,3) + 1144*Q2e*S*pow(Q2h,3) +
		           4960*l1k*U*pow(Q2h,3) + 512*m2*U*pow(Q2h,3) +
		           284*Q2e*U*pow(Q2h,3) - 672*S*U*pow(Q2h,3) -
		           7048*pow(l1k,2)*pow(Q2h,3) + 848*pow(Q2e,2)*pow(Q2h,3) -
		           2734*l1k*pow(Q2h,4) - 276*m2*pow(Q2h,4) -
		           1122*Q2e*pow(Q2h,4) + 412*S*pow(Q2h,4) + 928*U*pow(Q2h,4) -
		           166*pow(Q2h,5) - 64*l1k*m4*pow(S,2) -
		           96*l1k*m2*Q2e*pow(S,2) - 464*l1k*m2*Q2h*pow(S,2) -
		           64*m4*Q2h*pow(S,2) - 376*l1k*Q2e*Q2h*pow(S,2) -
		           8*m2*Q2e*Q2h*pow(S,2) - 192*m2*pow(l1k,2)*pow(S,2) -
		           80*Q2e*pow(l1k,2)*pow(S,2) - 1344*Q2h*pow(l1k,2)*pow(S,2) -
		           544*pow(l1k,3)*pow(S,2) + 48*l1k*pow(Q2e,2)*pow(S,2) +
		           160*Q2h*pow(Q2e,2)*pow(S,2) - 48*pow(Q2e,3)*pow(S,2) -
		           832*l1k*pow(Q2h,2)*pow(S,2) - 216*m2*pow(Q2h,2)*pow(S,2) -
		           292*Q2e*pow(Q2h,2)*pow(S,2) - 144*pow(Q2h,3)*pow(S,2) -
		           64*l1k*m4*pow(U,2) - 384*l1k*m2*Q2e*pow(U,2) -
		           64*m4*Q2e*pow(U,2) - 1104*l1k*m2*Q2h*pow(U,2) -
		           192*m4*Q2h*pow(U,2) - 880*l1k*Q2e*Q2h*pow(U,2) -
		           512*m2*Q2e*Q2h*pow(U,2) - 448*m2*pow(l1k,2)*pow(U,2) -
		           1584*Q2e*pow(l1k,2)*pow(U,2) - 4480*Q2h*pow(l1k,2)*pow(U,2) -
		           2112*pow(l1k,3)*pow(U,2) - 616*l1k*pow(Q2e,2)*pow(U,2) -
		           216*Q2h*pow(Q2e,2)*pow(U,2) - 68*pow(Q2e,3)*pow(U,2) -
		           2112*l1k*pow(Q2h,2)*pow(U,2) - 240*m2*pow(Q2h,2)*pow(U,2) -
		           152*Q2e*pow(Q2h,2)*pow(U,2) - 432*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*(-176*m4*Q2e*Q2h*S*U + 32*Q2e*Q2h*S*pow(m,6) +
		           32*Q2e*Q2h*U*pow(m,6) - 128*Q2h*S*U*pow(m,6) -
		           8*Q2h*S*U*pow(Q2e,3) + 4*M2*Q2h*pow(Q2e,4) +
		           6*Q2h*U*pow(Q2e,4) + 32*m4*M2*Q2e*pow(Q2h,2) +
		           304*m4*Q2e*S*pow(Q2h,2) + 304*m4*Q2e*U*pow(Q2h,2) -
		           176*m4*S*U*pow(Q2h,2) - 252*m2*Q2e*S*U*pow(Q2h,2) -
		           32*Q2e*pow(m,6)*pow(Q2h,2) + 96*S*pow(m,6)*pow(Q2h,2) +
		           96*U*pow(m,6)*pow(Q2h,2) - 32*S*U*pow(Q2e,2)*pow(Q2h,2) -
		           4*M2*pow(Q2e,3)*pow(Q2h,2) + 14*S*pow(Q2e,3)*pow(Q2h,2) +
		           28*U*pow(Q2e,3)*pow(Q2h,2) - 4*pow(Q2e,4)*pow(Q2h,2) -
		           32*m4*M2*pow(Q2h,3) - 216*m4*Q2e*pow(Q2h,3) +
		           24*m2*M2*Q2e*pow(Q2h,3) + 128*m4*S*pow(Q2h,3) +
		           354*m2*Q2e*S*pow(Q2h,3) + 96*m4*U*pow(Q2h,3) +
		           434*m2*Q2e*U*pow(Q2h,3) - 84*m2*S*U*pow(Q2h,3) -
		           48*Q2e*S*U*pow(Q2h,3) - 32*pow(m,6)*pow(Q2h,3) +
		           32*S*pow(Q2e,2)*pow(Q2h,3) + 40*U*pow(Q2e,2)*pow(Q2h,3) -
		           16*pow(Q2e,3)*pow(Q2h,3) - 24*m4*pow(Q2h,4) -
		           24*m2*M2*pow(Q2h,4) - 270*m2*Q2e*pow(Q2h,4) +
		           42*m2*S*pow(Q2h,4) + 52*Q2e*S*pow(Q2h,4) -
		           30*m2*U*pow(Q2h,4) + 44*Q2e*U*pow(Q2h,4) - 40*S*U*pow(Q2h,4) -
		           20*pow(Q2e,2)*pow(Q2h,4) + 38*m2*pow(Q2h,5) -
		           24*Q2e*pow(Q2h,5) + 40*S*pow(Q2h,5) + 40*U*pow(Q2h,5) -
		           20*pow(Q2h,6) - 112*m4*Q2e*Q2h*pow(S,2) -
		           64*Q2h*pow(m,6)*pow(S,2) - 80*m4*pow(Q2h,2)*pow(S,2) -
		           116*m2*Q2e*pow(Q2h,2)*pow(S,2) -
		           10*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) -
		           48*m2*pow(Q2h,3)*pow(S,2) - 24*Q2e*pow(Q2h,3)*pow(S,2) -
		           20*pow(Q2h,4)*pow(S,2) - 96*m4*Q2e*Q2h*pow(U,2) -
		           64*Q2h*pow(m,6)*pow(U,2) - 12*Q2h*pow(Q2e,3)*pow(U,2) -
		           2*pow(Q2e,4)*pow(U,2) - 64*m4*pow(Q2h,2)*pow(U,2) -
		           168*m2*Q2e*pow(Q2h,2)*pow(U,2) -
		           20*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) - 4*m2*pow(Q2h,3)*pow(U,2) -
		           20*Q2e*pow(Q2h,3)*pow(U,2) - 20*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        (-256*m2*Q2e*S*pow(l1k,4) + 1472*Q2e*Q2h*S*pow(l1k,4) +
		           128*m2*Q2e*pow(l1k,5) - 960*Q2e*Q2h*pow(l1k,5) -
		           640*Q2e*S*pow(l1k,5) - 640*Q2h*S*pow(l1k,5) +
		           384*Q2e*pow(l1k,6) + 384*Q2h*pow(l1k,6) + 256*S*pow(l1k,6) -
		           128*pow(l1k,7) - 728*l1k*Q2h*S*pow(Q2e,4) -
		           388*l1k*Q2h*pow(Q2e,5) + 216*l1k*S*pow(Q2e,5) -
		           460*Q2h*S*pow(Q2e,5) + 78*l1k*pow(Q2e,6) - 226*Q2h*pow(Q2e,6) +
		           84*S*pow(Q2e,6) + 32*pow(Q2e,7) +
		           512*Q2e*S*pow(l1k,3)*pow(Q2h,2) -
		           160*Q2e*pow(l1k,4)*pow(Q2h,2) + 768*S*pow(l1k,4)*pow(Q2h,2) -
		           544*pow(l1k,5)*pow(Q2h,2) + 448*l1k*S*pow(Q2e,3)*pow(Q2h,2) +
		           690*l1k*pow(Q2e,4)*pow(Q2h,2) + 992*S*pow(Q2e,4)*pow(Q2h,2) +
		           694*pow(Q2e,5)*pow(Q2h,2) + 592*Q2e*S*pow(l1k,2)*pow(Q2h,3) -
		           232*Q2e*pow(l1k,3)*pow(Q2h,3) + 896*S*pow(l1k,3)*pow(Q2h,3) -
		           480*pow(l1k,4)*pow(Q2h,3) + 1360*l1k*S*pow(Q2e,2)*pow(Q2h,3) -
		           186*l1k*pow(Q2e,3)*pow(Q2h,3) - 940*S*pow(Q2e,3)*pow(Q2h,3) -
		           1168*pow(Q2e,4)*pow(Q2h,3) - 1296*l1k*Q2e*S*pow(Q2h,4) -
		           640*Q2e*pow(l1k,2)*pow(Q2h,4) + 1296*S*pow(l1k,2)*pow(Q2h,4) -
		           648*pow(l1k,3)*pow(Q2h,4) - 1514*l1k*pow(Q2e,2)*pow(Q2h,4) +
		           324*S*pow(Q2e,2)*pow(Q2h,4) + 998*pow(Q2e,3)*pow(Q2h,4) +
		           1320*l1k*Q2e*pow(Q2h,5) - 1320*pow(l1k,2)*pow(Q2h,5) -
		           330*pow(Q2e,2)*pow(Q2h,5) + 128*m2*Q2e*pow(l1k,3)*pow(S,2) -
		           544*Q2e*Q2h*pow(l1k,3)*pow(S,2) + 256*Q2e*pow(l1k,4)*pow(S,2) +
		           256*Q2h*pow(l1k,4)*pow(S,2) - 128*pow(l1k,5)*pow(S,2) -
		           232*l1k*Q2h*pow(Q2e,3)*pow(S,2) + 136*l1k*pow(Q2e,4)*pow(S,2) -
		           184*Q2h*pow(Q2e,4)*pow(S,2) + 48*pow(Q2e,5)*pow(S,2) -
		           288*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(S,2) -
		           256*pow(l1k,3)*pow(Q2h,2)*pow(S,2) -
		           272*l1k*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) +
		           224*pow(Q2e,3)*pow(Q2h,2)*pow(S,2) +
		           368*l1k*Q2e*pow(Q2h,3)*pow(S,2) -
		           384*pow(l1k,2)*pow(Q2h,3)*pow(S,2) -
		           88*pow(Q2e,2)*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1)*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1)) +
		     SI.B0_qmm(Q2h,m2)*(-16*M2 + 12*Q2h - 8*S - 8*U -
		        32*M2*Q2h*pow(Q2e - Q2h,-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-24*M2*Q2h - 8*Q2e*Q2h + 4*Q2e*S + 16*Q2h*S + 28*Q2h*U -
		           8*S*U - 18*pow(Q2h,2) + 8*pow(S,2)) +
		        pow(l1k,-1)*(12*M2*Q2h + 4*Q2e*Q2h - 14*Q2h*S - 2*Q2e*U -
		           8*Q2h*U + 4*S*U + 9*pow(Q2h,2) - 4*pow(U,2)) +
		        pow(Q2e - Q2h,-2)*(-32*Q2h*S*U + 16*S*pow(Q2h,2) +
		           16*U*pow(Q2h,2) - 16*Q2h*pow(S,2) - 16*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(Q2e - Q2h,-1)*
		         (-12*S*pow(Q2h,2) + 12*U*pow(Q2h,2) + 8*Q2h*pow(S,2) -
		           8*Q2h*pow(U,2)) + pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (24*Q2h*S*U + 8*M2*pow(Q2h,2) - 4*S*pow(Q2h,2) -
		           16*U*pow(Q2h,2) - 4*pow(Q2h,3) + 4*Q2h*pow(S,2) + 12*Q2h*pow(U,2)\
		) + pow(l1k,-1)*pow(Q2e - Q2h,-2)*
		         (-8*S*pow(Q2h,3) + 8*U*pow(Q2h,3) + 8*pow(Q2h,2)*pow(S,2) -
		           8*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-2)*pow(Q2e - Q2h,-1)*
		         (-12*S*U*pow(Q2h,2) + 18*S*pow(Q2h,3) + 14*U*pow(Q2h,3) -
		           10*pow(Q2h,4) - 8*pow(Q2h,2)*pow(S,2) - 4*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (12*S*U*pow(Q2h,2) - 26*S*pow(Q2h,3) - 22*U*pow(Q2h,3) +
		           18*pow(Q2h,4) + 8*pow(Q2h,2)*pow(S,2) + 4*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-3)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-12*S*U*pow(Q2h,3) + 12*S*pow(Q2h,4) + 12*U*pow(Q2h,4) -
		           6*pow(Q2h,5) - 6*pow(Q2h,3)*pow(S,2) - 6*pow(Q2h,3)*pow(U,2)) +
		        pow(l1k,-2)*pow(Q2e - Q2h,-2)*
		         (-8*S*U*pow(Q2h,3) + 8*S*pow(Q2h,4) + 8*U*pow(Q2h,4) -
		           4*pow(Q2h,5) - 4*pow(Q2h,3)*pow(S,2) - 4*pow(Q2h,3)*pow(U,2)) +
		        pow(l1k,-3)*pow(Q2e - Q2h,-1)*
		         (4*S*U*pow(Q2h,3) - 4*S*pow(Q2h,4) - 4*U*pow(Q2h,4) +
		           2*pow(Q2h,5) + 2*pow(Q2h,3)*pow(S,2) + 2*pow(Q2h,3)*pow(U,2)) +
		        (288*Q2e*Q2h*S*pow(l1k,2) - 192*Q2e*Q2h*pow(l1k,3) +
		           192*Q2h*S*pow(l1k,3) - 96*Q2h*pow(l1k,4) +
		           144*l1k*Q2h*S*pow(Q2e,2) - 144*Q2h*pow(l1k,2)*pow(Q2e,2) -
		           48*l1k*Q2h*pow(Q2e,3) + 24*Q2h*S*pow(Q2e,3) -
		           6*Q2h*pow(Q2e,4) - 96*l1k*Q2e*Q2h*pow(S,2) -
		           96*Q2h*pow(l1k,2)*pow(S,2) - 24*Q2h*pow(Q2e,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        (32*l1k*M2*Q2e + 16*l1k*Q2e*S + 248*l1k*Q2h*S - 32*Q2e*Q2h*S +
		           24*l1k*Q2e*U + 16*l1k*Q2h*U + 24*Q2e*Q2h*U - 16*l1k*S*U -
		           16*Q2e*S*U - 16*Q2h*S*U + 32*M2*pow(l1k,2) -
		           144*Q2h*pow(l1k,2) + 16*S*pow(l1k,2) + 16*U*pow(l1k,2) +
		           8*M2*pow(Q2e,2) - 52*Q2h*pow(Q2e,2) + 4*S*pow(Q2e,2) +
		           12*U*pow(Q2e,2) - 28*l1k*pow(Q2h,2) + 30*Q2e*pow(Q2h,2) +
		           120*S*pow(Q2h,2) + 8*U*pow(Q2h,2) + 2*pow(Q2h,3) -
		           16*l1k*pow(S,2) - 8*Q2e*pow(S,2) - 96*Q2h*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (24*Q2e*S*pow(Q2h,2) - 16*S*U*pow(Q2h,2) - 8*Q2e*pow(Q2h,3) +
		           8*S*pow(Q2h,3) + 8*U*pow(Q2h,3) - 2*pow(Q2h,4) -
		           16*Q2e*Q2h*pow(S,2) - 8*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*(-16*Q2e*Q2h*S*U + 14*Q2h*S*pow(Q2e,2) +
		           8*Q2h*U*pow(Q2e,2) - 4*S*U*pow(Q2e,2) - 4*Q2h*pow(Q2e,3) +
		           2*U*pow(Q2e,3) + 16*Q2e*S*pow(Q2h,2) -
		           6*pow(Q2e,2)*pow(Q2h,2) + 2*Q2e*pow(Q2h,3) - 8*S*pow(Q2h,3) +
		           2*pow(Q2h,4) - 12*Q2e*Q2h*pow(S,2) - 8*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*S*pow(Q2h,4) + 2*pow(Q2h,5) + 8*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (-192*Q2h*U*pow(l1k,3) - 96*Q2h*pow(l1k,4) +
		           288*U*pow(l1k,2)*pow(Q2h,2) + 192*pow(l1k,3)*pow(Q2h,2) -
		           144*l1k*U*pow(Q2h,3) - 144*pow(l1k,2)*pow(Q2h,3) +
		           48*l1k*pow(Q2h,4) + 24*U*pow(Q2h,4) - 6*pow(Q2h,5) -
		           96*Q2h*pow(l1k,2)*pow(U,2) + 96*l1k*pow(Q2h,2)*pow(U,2) -
		           24*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (-32*l1k*M2*Q2h + 16*l1k*Q2e*Q2h - 8*l1k*Q2e*S - 32*l1k*Q2h*S +
		           20*Q2e*Q2h*S - 8*l1k*Q2h*U + 36*Q2e*Q2h*U + 16*l1k*S*U -
		           8*Q2e*S*U - 24*Q2h*S*U + 32*M2*pow(l1k,2) +
		           16*Q2h*pow(l1k,2) + 16*S*pow(l1k,2) + 16*U*pow(l1k,2) -
		           8*Q2h*pow(Q2e,2) + 4*S*pow(Q2e,2) - 20*l1k*pow(Q2h,2) +
		           8*M2*pow(Q2h,2) - 28*Q2e*pow(Q2h,2) + 20*S*pow(Q2h,2) -
		           8*U*pow(Q2h,2) + 24*pow(Q2h,3) + 16*l1k*pow(U,2) -
		           8*Q2h*pow(U,2))*pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2h,2),-1) + pow(l1k,-1)*
		         (-12*Q2e*U*pow(Q2h,2) + 8*S*U*pow(Q2h,2) + 4*Q2e*pow(Q2h,3) -
		           4*S*pow(Q2h,3) - 4*U*pow(Q2h,3) + pow(Q2h,4) +
		           8*Q2e*Q2h*pow(U,2) + 4*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (32*Q2e*Q2h*S*U - 16*Q2h*S*pow(Q2e,2) - 28*Q2h*U*pow(Q2e,2) +
		           8*S*U*pow(Q2e,2) + 8*Q2h*pow(Q2e,3) - 4*S*pow(Q2e,3) -
		           32*Q2e*U*pow(Q2h,2) + 12*pow(Q2e,2)*pow(Q2h,2) -
		           4*Q2e*pow(Q2h,3) + 16*U*pow(Q2h,3) - 4*pow(Q2h,4) +
		           24*Q2e*Q2h*pow(U,2) + 16*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*U*pow(Q2h,4) + 2*pow(Q2h,5) + 8*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (16*Q2e*S*pow(Q2h,3) - 16*S*U*pow(Q2h,3) - 8*Q2e*pow(Q2h,4) +
		           16*S*pow(Q2h,4) + 16*U*pow(Q2h,4) - 8*pow(Q2h,5) -
		           8*Q2e*pow(Q2h,2)*pow(S,2) - 8*pow(Q2h,3)*pow(S,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        (-128*l1k*Q2e*Q2h*S + 32*m2*Q2e*Q2h*S + 32*m2*Q2e*Q2h*U -
		           64*m2*Q2h*S*U - 16*Q2e*Q2h*S*U + 64*Q2e*Q2h*pow(l1k,2) -
		           32*Q2e*S*pow(l1k,2) - 32*Q2h*S*pow(l1k,2) -
		           32*Q2e*U*pow(l1k,2) - 32*Q2h*U*pow(l1k,2) -
		           28*l1k*Q2h*pow(Q2e,2) + 32*l1k*S*pow(Q2e,2) +
		           72*Q2h*S*pow(Q2e,2) - 16*pow(l1k,2)*pow(Q2e,2) +
		           12*l1k*pow(Q2e,3) + 30*Q2h*pow(Q2e,3) - 24*S*pow(Q2e,3) -
		           8*pow(Q2e,4) + 36*l1k*Q2e*pow(Q2h,2) - 32*m2*Q2e*pow(Q2h,2) +
		           64*l1k*S*pow(Q2h,2) + 32*m2*S*pow(Q2h,2) -
		           48*Q2e*S*pow(Q2h,2) + 32*l1k*U*pow(Q2h,2) +
		           32*m2*U*pow(Q2h,2) + 40*Q2e*U*pow(Q2h,2) -
		           16*pow(l1k,2)*pow(Q2h,2) - 38*pow(Q2e,2)*pow(Q2h,2) -
		           20*l1k*pow(Q2h,3) - 6*Q2e*pow(Q2h,3) + 32*S*pow(Q2h,3) -
		           8*U*pow(Q2h,3) + 6*pow(Q2h,4) + 16*l1k*Q2e*pow(S,2) +
		           16*l1k*Q2h*pow(S,2) - 32*m2*Q2h*pow(S,2) +
		           32*Q2e*Q2h*pow(S,2) - 16*pow(Q2e,2)*pow(S,2) -
		           24*pow(Q2h,2)*pow(S,2) - 16*l1k*Q2e*pow(U,2) -
		           16*l1k*Q2h*pow(U,2) - 32*m2*Q2h*pow(U,2) -
		           16*Q2e*Q2h*pow(U,2) + 8*pow(Q2h,2)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*(16*Q2e*S*pow(Q2h,3) + 8*Q2e*U*pow(Q2h,3) +
		           8*S*U*pow(Q2h,3) - 12*Q2e*pow(Q2h,4) - 8*S*pow(Q2h,4) -
		           8*U*pow(Q2h,4) + 4*pow(Q2h,5) + 4*Q2e*pow(Q2h,2)*pow(U,2) +
		           4*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-2)*(-8*m2*Q2e*S*pow(Q2h,3) - 8*m2*Q2e*U*pow(Q2h,3) +
		           16*m2*S*U*pow(Q2h,3) + 16*Q2e*S*U*pow(Q2h,3) +
		           8*m2*Q2e*pow(Q2h,4) - 8*m2*S*pow(Q2h,4) -
		           16*Q2e*S*pow(Q2h,4) - 8*m2*U*pow(Q2h,4) -
		           16*Q2e*U*pow(Q2h,4) + 8*Q2e*pow(Q2h,5) +
		           8*m2*pow(Q2h,3)*pow(S,2) + 8*Q2e*pow(Q2h,3)*pow(S,2) +
		           8*m2*pow(Q2h,3)*pow(U,2) + 8*Q2e*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (8*S*U*pow(Q2h,4) - 8*S*pow(Q2h,5) - 8*U*pow(Q2h,5) +
		           4*pow(Q2h,6) + 4*pow(Q2h,4)*pow(S,2) + 4*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-3)*(-8*m2*Q2e*S*U*pow(Q2h,3) + 8*m2*Q2e*S*pow(Q2h,4) +
		           8*m2*Q2e*U*pow(Q2h,4) + 8*m2*S*U*pow(Q2h,4) -
		           4*m2*Q2e*pow(Q2h,5) - 8*m2*S*pow(Q2h,5) -
		           8*m2*U*pow(Q2h,5) + 4*m2*pow(Q2h,6) -
		           4*m2*Q2e*pow(Q2h,3)*pow(S,2) + 4*m2*pow(Q2h,4)*pow(S,2) -
		           4*m2*Q2e*pow(Q2h,3)*pow(U,2) + 4*m2*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        (-384*Q2e*Q2h*S*pow(l1k,3) + 224*Q2e*Q2h*pow(l1k,4) +
		           128*Q2e*S*pow(l1k,4) + 128*Q2h*S*pow(l1k,4) -
		           64*Q2e*pow(l1k,5) - 64*Q2h*pow(l1k,5) -
		           32*l1k*Q2h*S*pow(Q2e,3) - 28*l1k*Q2h*pow(Q2e,4) +
		           64*l1k*S*pow(Q2e,4) - 72*Q2h*S*pow(Q2e,4) + 20*l1k*pow(Q2e,5) -
		           30*Q2h*pow(Q2e,5) + 24*S*pow(Q2e,5) + 8*pow(Q2e,6) +
		           64*Q2e*S*pow(l1k,2)*pow(Q2h,2) - 96*Q2e*pow(l1k,3)*pow(Q2h,2) -
		           256*S*pow(l1k,3)*pow(Q2h,2) + 160*pow(l1k,4)*pow(Q2h,2) -
		           160*l1k*S*pow(Q2e,2)*pow(Q2h,2) -
		           40*l1k*pow(Q2e,3)*pow(Q2h,2) + 64*S*pow(Q2e,3)*pow(Q2h,2) +
		           38*pow(Q2e,4)*pow(Q2h,2) + 96*l1k*Q2e*S*pow(Q2h,3) -
		           64*Q2e*pow(l1k,2)*pow(Q2h,3) - 128*S*pow(l1k,2)*pow(Q2h,3) +
		           32*pow(l1k,3)*pow(Q2h,3) + 56*l1k*pow(Q2e,2)*pow(Q2h,3) -
		           16*S*pow(Q2e,2)*pow(Q2h,3) - 16*pow(Q2e,3)*pow(Q2h,3) +
		           160*Q2e*Q2h*pow(l1k,2)*pow(S,2) - 64*Q2e*pow(l1k,3)*pow(S,2) -
		           64*Q2h*pow(l1k,3)*pow(S,2) + 32*l1k*Q2h*pow(Q2e,2)*pow(S,2) +
		           48*l1k*pow(Q2e,3)*pow(S,2) - 32*Q2h*pow(Q2e,3)*pow(S,2) +
		           16*pow(Q2e,4)*pow(S,2) - 64*l1k*Q2e*pow(Q2h,2)*pow(S,2) +
		           96*pow(l1k,2)*pow(Q2h,2)*pow(S,2) +
		           24*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) +
		           16*l1k*pow(Q2h,3)*pow(S,2) - 8*Q2e*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1)*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        (-128*Q2e*Q2h*U*pow(l1k,3) - 96*Q2e*Q2h*pow(l1k,4) +
		           128*Q2e*U*pow(l1k,4) + 128*Q2h*U*pow(l1k,4) +
		           64*Q2e*pow(l1k,5) + 64*Q2h*pow(l1k,5) +
		           32*Q2e*pow(l1k,3)*pow(Q2h,2) - 256*U*pow(l1k,3)*pow(Q2h,2) -
		           160*pow(l1k,4)*pow(Q2h,2) + 32*l1k*U*pow(Q2e,2)*pow(Q2h,2) +
		           64*l1k*Q2e*U*pow(Q2h,3) + 192*U*pow(l1k,2)*pow(Q2h,3) +
		           160*pow(l1k,3)*pow(Q2h,3) - 16*l1k*pow(Q2e,2)*pow(Q2h,3) +
		           16*U*pow(Q2e,2)*pow(Q2h,3) - 28*l1k*Q2e*pow(Q2h,4) -
		           64*l1k*U*pow(Q2h,4) - 24*Q2e*U*pow(Q2h,4) -
		           96*pow(l1k,2)*pow(Q2h,4) - 8*pow(Q2e,2)*pow(Q2h,4) +
		           36*l1k*pow(Q2h,5) + 14*Q2e*pow(Q2h,5) + 8*U*pow(Q2h,5) -
		           6*pow(Q2h,6) - 32*Q2e*Q2h*pow(l1k,2)*pow(U,2) +
		           64*Q2e*pow(l1k,3)*pow(U,2) + 64*Q2h*pow(l1k,3)*pow(U,2) -
		           16*l1k*Q2h*pow(Q2e,2)*pow(U,2) -
		           48*l1k*Q2e*pow(Q2h,2)*pow(U,2) -
		           96*pow(l1k,2)*pow(Q2h,2)*pow(U,2) -
		           8*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) + 32*l1k*pow(Q2h,3)*pow(U,2) +
		           8*Q2e*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1)) +
		     SI.B0_M0m(-2*l1k + m2,m2)*(2*l1k - 12*M2 + 4*Q2e + Q2h - 12*S - 10*U +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (12*m2*Q2h - 8*M2*Q2h - 2*Q2e*Q2h - 24*m2*S + 12*Q2e*S -
		           30*Q2h*S + 14*Q2e*U - 4*Q2h*U - 24*S*U - 4*pow(Q2e,2) +
		           24*pow(Q2h,2) + 4*pow(S,2) - 12*pow(U,2)) +
		        pow(l1k,-1)*(8*m2*M2 + 8*M2*Q2e - 2*m2*Q2h + 8*M2*Q2h -
		           Q2e*Q2h + 6*m2*S + 11*Q2h*S - 6*m2*U + Q2e*U + 17*Q2h*U -
		           2*S*U - 9*pow(Q2h,2) - 8*pow(U,2)) +
		        pow(2*l1k - m2,-1)*(-16*l1k*Q2h + 2*Q2e*Q2h + 12*l1k*S - 6*Q2e*S +
		           2*Q2h*S + 12*l1k*U - 24*Q2h*U + 12*S*U + 8*pow(Q2h,2) +
		           12*pow(U,2)) + pow(l1k,-2)*
		         (-4*m2*M2*Q2e - 4*m4*Q2h - 4*m2*M2*Q2h - 3*m2*Q2e*Q2h -
		           2*m4*S + 5*m2*Q2h*S + 8*m4*U + 3*m2*Q2e*U + 9*m2*Q2h*U -
		           6*m2*S*U - 5*m2*pow(Q2h,2) - 4*m2*pow(U,2)) +
		        pow(l1k,-3)*(m4*Q2e*Q2h - 3*m4*Q2h*S - m4*Q2e*U - 4*m4*Q2h*U +
		           2*m4*S*U + 2*m4*pow(Q2h,2) + 2*m4*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*m4*Q2h + 24*m2*M2*Q2h + 20*m4*S + 62*m2*Q2h*S +
		           12*m2*Q2h*U - 24*m2*S*U - 26*Q2h*S*U - 18*m2*pow(Q2h,2) +
		           8*M2*pow(Q2h,2) + 59*S*pow(Q2h,2) + 37*U*pow(Q2h,2) -
		           32*pow(Q2h,3) - 44*m2*pow(S,2) - 26*Q2h*pow(S,2) -
		           8*Q2h*pow(U,2)) + pow(2*l1k - m2,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*Q2e*Q2h*S + 2*Q2e*Q2h*U - 12*Q2e*S*U + 16*Q2h*S*U -
		           2*Q2h*pow(Q2e,2) + 6*S*pow(Q2e,2) + 4*Q2e*pow(Q2h,2) -
		           8*U*pow(Q2h,2) + 4*Q2h*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-26*m4*Q2h*S - 6*m4*Q2h*U + 20*m4*S*U + 6*m2*Q2h*S*U +
		           6*m4*pow(Q2h,2) - 17*m2*S*pow(Q2h,2) - 11*m2*U*pow(Q2h,2) +
		           20*S*U*pow(Q2h,2) + 11*m2*pow(Q2h,3) - 24*S*pow(Q2h,3) -
		           20*U*pow(Q2h,3) + 12*pow(Q2h,4) + 20*m4*pow(S,2) +
		           6*m2*Q2h*pow(S,2) + 12*pow(Q2h,2)*pow(S,2) +
		           8*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-3)*pow(2*l1k + Q2e - Q2h,-1)*
		         (4*m4*Q2h*S*U - 4*m4*S*pow(Q2h,2) - 4*m4*U*pow(Q2h,2) -
		           8*m2*S*U*pow(Q2h,2) + 2*m4*pow(Q2h,3) + 8*m2*S*pow(Q2h,3) +
		           8*m2*U*pow(Q2h,3) - 2*S*U*pow(Q2h,3) - 4*m2*pow(Q2h,4) +
		           2*S*pow(Q2h,4) + 2*U*pow(Q2h,4) - pow(Q2h,5) +
		           2*m4*Q2h*pow(S,2) - 4*m2*pow(Q2h,2)*pow(S,2) -
		           pow(Q2h,3)*pow(S,2) + 2*m4*Q2h*pow(U,2) -
		           4*m2*pow(Q2h,2)*pow(U,2) - pow(Q2h,3)*pow(U,2)) +
		        (-192*Q2h*U*pow(l1k,3) - 144*Q2h*pow(l1k,4) + 192*U*pow(l1k,4) +
		           96*pow(l1k,5) + 48*pow(l1k,3)*pow(Q2h,2) +
		           48*l1k*U*pow(Q2h,3) + 24*pow(l1k,2)*pow(Q2h,3) -
		           18*l1k*pow(Q2h,4) - 12*U*pow(Q2h,4) + 3*pow(Q2h,5) -
		           48*Q2h*pow(l1k,2)*pow(U,2) + 96*pow(l1k,3)*pow(U,2) -
		           24*l1k*pow(Q2h,2)*pow(U,2) + 12*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (32*l1k*m2*M2 + 16*l1k*M2*Q2h + 20*l1k*Q2e*Q2h +
		           16*l1k*m2*S - 8*l1k*Q2e*S - 8*m2*Q2e*S - 12*Q2e*Q2h*S +
		           16*l1k*m2*U - 36*l1k*Q2e*U + 72*l1k*Q2h*U + 10*Q2e*Q2h*U +
		           16*m2*S*U - 8*Q2e*S*U + 8*Q2h*S*U - 16*M2*pow(l1k,2) -
		           16*Q2e*pow(l1k,2) + 48*Q2h*pow(l1k,2) + 8*S*pow(l1k,2) -
		           48*U*pow(l1k,2) - 32*pow(l1k,3) + 8*l1k*pow(Q2e,2) -
		           2*Q2h*pow(Q2e,2) + 4*S*pow(Q2e,2) + 14*U*pow(Q2e,2) -
		           4*pow(Q2e,3) - 52*l1k*pow(Q2h,2) - 4*M2*pow(Q2h,2) +
		           12*Q2e*pow(Q2h,2) - 6*S*pow(Q2h,2) - 80*U*pow(Q2h,2) +
		           32*pow(Q2h,3) - 32*l1k*pow(U,2) + 16*m2*pow(U,2) -
		           20*Q2e*pow(U,2) + 32*Q2h*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*(6*Q2e*U*pow(Q2h,2) - 4*S*U*pow(Q2h,2) -
		           2*Q2e*pow(Q2h,3) + 2*S*pow(Q2h,3) - 4*Q2e*Q2h*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*m2*Q2e*S*U - 32*Q2e*Q2h*S*U + 8*m2*S*pow(Q2e,2) +
		           16*Q2h*S*pow(Q2e,2) + 12*Q2h*U*pow(Q2e,2) + 8*S*U*pow(Q2e,2) -
		           2*Q2h*pow(Q2e,3) - 4*S*pow(Q2e,3) - 14*U*pow(Q2e,3) +
		           4*pow(Q2e,4) + 40*Q2e*U*pow(Q2h,2) - 14*pow(Q2e,2)*pow(Q2h,2) +
		           2*Q2e*pow(Q2h,3) - 8*U*pow(Q2h,3) + 2*pow(Q2h,4) -
		           16*Q2e*Q2h*pow(U,2) + 12*pow(Q2e,2)*pow(U,2) -
		           24*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k - m2,-1)*(32*l1k*U*pow(Q2h,2) - 16*l1k*pow(Q2h,3) +
		           16*U*pow(Q2h,3) - 8*pow(Q2h,4) - 16*l1k*Q2h*pow(U,2) -
		           8*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (4*U*pow(Q2h,4) - pow(Q2h,5) - 4*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        (-16*l1k*m2*Q2e*Q2h + 16*m4*Q2e*Q2h + 32*l1k*Q2e*Q2h*S -
		           96*m2*Q2e*Q2h*S - 32*l1k*m2*Q2e*U - 32*l1k*m2*Q2h*U -
		           40*l1k*Q2e*Q2h*U - 80*m2*Q2e*Q2h*U + 128*l1k*m2*S*U -
		           64*m4*S*U + 80*m2*Q2e*S*U + 32*Q2e*Q2h*S*U -
		           16*m2*Q2e*pow(l1k,2) - 16*m2*Q2h*pow(l1k,2) -
		           8*Q2e*Q2h*pow(l1k,2) + 16*Q2e*U*pow(l1k,2) -
		           16*Q2h*U*pow(l1k,2) - 16*Q2h*pow(l1k,3) + 32*U*pow(l1k,3) +
		           16*pow(l1k,4) - 16*Q2h*U*pow(Q2e,2) + 16*m4*pow(Q2h,2) -
		           4*l1k*Q2e*pow(Q2h,2) + 60*m2*Q2e*pow(Q2h,2) -
		           96*Q2e*S*pow(Q2h,2) - 8*l1k*U*pow(Q2h,2) -
		           16*m2*U*pow(Q2h,2) - 28*Q2e*U*pow(Q2h,2) +
		           8*pow(Q2e,2)*pow(Q2h,2) + 4*m2*pow(Q2h,3) +
		           47*Q2e*pow(Q2h,3) - 4*U*pow(Q2h,3) + pow(Q2h,4) +
		           64*l1k*m2*pow(S,2) - 32*m4*pow(S,2) - 16*l1k*Q2e*pow(S,2) +
		           72*m2*Q2e*pow(S,2) - 32*m2*Q2h*pow(S,2) +
		           40*Q2e*Q2h*pow(S,2) + 64*l1k*m2*pow(U,2) - 32*m4*pow(U,2) +
		           32*l1k*Q2e*pow(U,2) + 24*m2*Q2e*pow(U,2) +
		           16*m2*Q2h*pow(U,2) - 4*Q2e*Q2h*pow(U,2) +
		           16*pow(l1k,2)*pow(U,2) + 8*pow(Q2e,2)*pow(U,2) -
		           4*pow(Q2h,2)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*(96*m4*Q2e*Q2h*S + 64*m4*Q2e*Q2h*U -
		           96*m4*Q2e*S*U + 64*m4*Q2h*S*U - 40*m2*Q2e*Q2h*S*U -
		           32*m4*Q2e*pow(Q2h,2) - 64*m4*S*pow(Q2h,2) +
		           80*m2*Q2e*S*pow(Q2h,2) - 32*m4*U*pow(Q2h,2) +
		           48*m2*Q2e*U*pow(Q2h,2) - 8*m2*S*U*pow(Q2h,2) -
		           40*Q2e*S*U*pow(Q2h,2) + 16*m4*pow(Q2h,3) -
		           44*m2*Q2e*pow(Q2h,3) - 8*m2*S*pow(Q2h,3) +
		           48*Q2e*S*pow(Q2h,3) + 8*m2*U*pow(Q2h,3) +
		           40*Q2e*U*pow(Q2h,3) + 4*m2*pow(Q2h,4) - 24*Q2e*pow(Q2h,4) -
		           64*m4*Q2e*pow(S,2) + 48*m4*Q2h*pow(S,2) -
		           32*m2*Q2e*Q2h*pow(S,2) - 24*Q2e*pow(Q2h,2)*pow(S,2) -
		           32*m4*Q2e*pow(U,2) + 16*m4*Q2h*pow(U,2) -
		           8*m2*Q2e*Q2h*pow(U,2) - 8*m2*pow(Q2h,2)*pow(U,2) -
		           16*Q2e*pow(Q2h,2)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*Q2e*S*U*pow(Q2h,2) + 16*S*pow(Q2e,2)*pow(Q2h,2) +
		           16*Q2e*U*pow(Q2h,3) - 8*pow(Q2e,2)*pow(Q2h,3) -
		           8*Q2h*pow(Q2e,2)*pow(S,2) - 8*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-2)*(16*m4*Q2e*Q2h*S*U - 32*Q2e*Q2h*S*pow(m,6) -
		           32*Q2e*Q2h*U*pow(m,6) + 32*Q2e*S*U*pow(m,6) -
		           32*Q2h*S*U*pow(m,6) - 24*m4*Q2e*S*pow(Q2h,2) -
		           24*m4*Q2e*U*pow(Q2h,2) + 36*m2*Q2e*S*U*pow(Q2h,2) +
		           16*Q2e*pow(m,6)*pow(Q2h,2) + 32*S*pow(m,6)*pow(Q2h,2) +
		           32*U*pow(m,6)*pow(Q2h,2) + 16*m4*Q2e*pow(Q2h,3) +
		           8*m4*S*pow(Q2h,3) - 40*m2*Q2e*S*pow(Q2h,3) +
		           8*m4*U*pow(Q2h,3) - 36*m2*Q2e*U*pow(Q2h,3) -
		           16*m2*S*U*pow(Q2h,3) + 4*Q2e*S*U*pow(Q2h,3) -
		           16*pow(m,6)*pow(Q2h,3) - 8*m4*pow(Q2h,4) +
		           20*m2*Q2e*pow(Q2h,4) + 20*m2*S*pow(Q2h,4) -
		           4*Q2e*S*pow(Q2h,4) + 16*m2*U*pow(Q2h,4) - 4*Q2e*U*pow(Q2h,4) -
		           10*m2*pow(Q2h,5) + 2*Q2e*pow(Q2h,5) + 8*m4*Q2e*Q2h*pow(S,2) +
		           16*Q2e*pow(m,6)*pow(S,2) - 16*Q2h*pow(m,6)*pow(S,2) +
		           20*m2*Q2e*pow(Q2h,2)*pow(S,2) - 10*m2*pow(Q2h,3)*pow(S,2) +
		           2*Q2e*pow(Q2h,3)*pow(S,2) + 8*m4*Q2e*Q2h*pow(U,2) +
		           16*Q2e*pow(m,6)*pow(U,2) - 16*Q2h*pow(m,6)*pow(U,2) +
		           16*m2*Q2e*pow(Q2h,2)*pow(U,2) - 6*m2*pow(Q2h,3)*pow(U,2) +
		           2*Q2e*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-3)*(-8*m4*Q2e*S*U*pow(Q2h,2) + 8*m4*Q2e*S*pow(Q2h,3) +
		           8*m4*Q2e*U*pow(Q2h,3) + 8*m4*S*U*pow(Q2h,3) -
		           2*m2*Q2e*S*U*pow(Q2h,3) - 4*m4*Q2e*pow(Q2h,4) -
		           8*m4*S*pow(Q2h,4) + 2*m2*Q2e*S*pow(Q2h,4) -
		           8*m4*U*pow(Q2h,4) + 2*m2*Q2e*U*pow(Q2h,4) +
		           2*m2*S*U*pow(Q2h,4) + 4*m4*pow(Q2h,5) - m2*Q2e*pow(Q2h,5) -
		           2*m2*S*pow(Q2h,5) - 2*m2*U*pow(Q2h,5) + m2*pow(Q2h,6) -
		           4*m4*Q2e*pow(Q2h,2)*pow(S,2) + 4*m4*pow(Q2h,3)*pow(S,2) -
		           m2*Q2e*pow(Q2h,3)*pow(S,2) + m2*pow(Q2h,4)*pow(S,2) -
		           4*m4*Q2e*pow(Q2h,2)*pow(U,2) + 4*m4*pow(Q2h,3)*pow(U,2) -
		           m2*Q2e*pow(Q2h,3)*pow(U,2) + m2*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        (128*m2*Q2e*U*pow(l1k,3) + 96*Q2e*Q2h*U*pow(l1k,3) +
		           64*m2*Q2e*pow(l1k,4) + 32*Q2e*Q2h*pow(l1k,4) -
		           64*Q2e*U*pow(l1k,4) + 192*Q2h*U*pow(l1k,4) +
		           128*Q2h*pow(l1k,5) - 128*U*pow(l1k,5) - 64*pow(l1k,6) -
		           16*Q2e*pow(l1k,3)*pow(Q2h,2) - 64*U*pow(l1k,3)*pow(Q2h,2) -
		           80*pow(l1k,4)*pow(Q2h,2) - 8*l1k*Q2e*U*pow(Q2h,3) -
		           4*Q2e*pow(l1k,2)*pow(Q2h,3) + 16*pow(l1k,3)*pow(Q2h,3) -
		           8*l1k*U*pow(Q2h,4) - 4*Q2e*U*pow(Q2h,4) -
		           4*pow(l1k,2)*pow(Q2h,4) + 4*l1k*pow(Q2h,5) + Q2e*pow(Q2h,5) +
		           4*U*pow(Q2h,5) - pow(Q2h,6) + 64*m2*Q2e*pow(l1k,2)*pow(U,2) +
		           48*Q2e*Q2h*pow(l1k,2)*pow(U,2) - 64*Q2e*pow(l1k,3)*pow(U,2) +
		           64*Q2h*pow(l1k,3)*pow(U,2) - 64*pow(l1k,4)*pow(U,2) +
		           16*l1k*Q2e*pow(Q2h,2)*pow(U,2) + 4*Q2e*pow(Q2h,3)*pow(U,2) -
		           4*pow(Q2h,4)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1)) +
		     SI.B0_m0m(m2)*(-16*l1k + 36*M2 - 11*Q2e + Q2h + 32*S + 16*U -
		        32*M2*Q2e*pow(4*m2 + Q2e,-1) +
		        pow(2*l1k + Q2e - Q2h,-2)*
		         (16*m2*M2*Q2e + 16*m4*Q2h + 16*m2*M2*Q2h + 8*m2*Q2e*Q2h -
		           32*m4*S - 8*m2*Q2e*S - 20*m2*Q2h*S + 8*m4*U - 8*m2*Q2h*U +
		           16*m2*S*U + 12*m2*pow(Q2h,2) + 8*m2*pow(S,2)) +
		        pow(2*l1k + Q2e - Q2h,-3)*
		         (8*m4*Q2e*Q2h - 8*m4*Q2e*S - 32*m4*Q2h*S - 24*m4*Q2h*U +
		           16*m4*S*U + 16*m4*pow(Q2h,2) + 16*m4*pow(S,2)) +
		        pow(l1k,-1)*(-8*m2*M2 - 4*M2*Q2e + 8*m2*Q2h - 20*M2*Q2h -
		           3*Q2e*Q2h - 4*m2*S + 7*Q2e*S + 4*Q2h*S - 4*m2*U + 4*Q2e*U -
		           12*Q2h*U - 8*S*U - 2*pow(Q2e,2) + pow(Q2h,2) - 6*pow(S,2) +
		           8*pow(U,2)) + pow(2*l1k + Q2e - Q2h,-1)*
		         (16*m2*M2 + 8*M2*Q2e - 16*m2*Q2h + 40*M2*Q2h + 6*Q2e*Q2h +
		           8*m2*S - 8*Q2e*S + 24*Q2h*S + 8*m2*U - 14*Q2e*U - 8*Q2h*U +
		           16*S*U + 4*pow(Q2e,2) - 2*pow(Q2h,2) - 16*pow(S,2) + 12*pow(U,2)) \
		+ pow(l1k,-2)*(4*m2*M2*Q2e + 4*m4*Q2h + 4*m2*M2*Q2h +
		           2*m2*Q2e*Q2h + 2*m4*S - 2*m2*Q2h*S - 8*m4*U - 2*m2*Q2e*U -
		           5*m2*Q2h*U + 4*m2*S*U + 3*m2*pow(Q2h,2) + 2*m2*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (16*m4*Q2h - 48*m2*M2*Q2h - 20*m4*S - 74*m2*Q2h*S -
		           20*m4*U - 42*m2*Q2h*U + 8*m2*S*U + 8*Q2h*S*U +
		           40*m2*pow(Q2h,2) - 16*M2*pow(Q2h,2) - 44*S*pow(Q2h,2) -
		           12*U*pow(Q2h,2) + 18*pow(Q2h,3) + 40*m2*pow(S,2) +
		           24*Q2h*pow(S,2) + 8*m2*pow(U,2)) +
		        pow(l1k,-3)*(-(m4*Q2e*Q2h) + 3*m4*Q2h*S + m4*Q2e*U +
		           4*m4*Q2h*U - 2*m4*S*U - 2*m4*pow(Q2h,2) - 2*m4*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-2)*
		         (52*m4*Q2h*S + 12*m4*Q2h*U - 24*m4*S*U - 12*m2*Q2h*S*U -
		           20*m4*pow(Q2h,2) + 18*m2*S*pow(Q2h,2) + 6*m2*U*pow(Q2h,2) -
		           6*m2*pow(Q2h,3) - 32*m4*pow(S,2) - 12*m2*Q2h*pow(S,2) +
		           8*m4*pow(U,2)) + pow(4*m2 + Q2e,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*M2*Q2e*Q2h - 8*Q2e*Q2h*S + 8*Q2e*Q2h*U + 16*Q2h*S*U +
		           8*M2*pow(Q2e,2) - 8*Q2h*pow(S,2) - 8*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(4*m2 + Q2e,-1)*
		         (4*M2*Q2e*Q2h - 4*Q2e*Q2h*S + 4*Q2e*Q2h*U - 8*Q2h*S*U -
		           4*M2*pow(Q2e,2) + 4*Q2h*pow(S,2) + 4*Q2h*pow(U,2)) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-3)*
		         (-16*m4*Q2h*S*U + 16*m4*S*pow(Q2h,2) + 16*m4*U*pow(Q2h,2) -
		           8*m4*pow(Q2h,3) - 8*m4*Q2h*pow(S,2) - 8*m4*Q2h*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (58*m4*Q2h*S + 38*m4*Q2h*U - 52*m4*S*U - 42*m2*Q2h*S*U -
		           22*m4*pow(Q2h,2) + 61*m2*S*pow(Q2h,2) + 55*m2*U*pow(Q2h,2) -
		           16*S*U*pow(Q2h,2) - 37*m2*pow(Q2h,3) + 24*S*pow(Q2h,3) +
		           16*U*pow(Q2h,3) - 12*pow(Q2h,4) - 36*m4*pow(S,2) -
		           24*m2*Q2h*pow(S,2) - 12*pow(Q2h,2)*pow(S,2) - 16*m4*pow(U,2) -
		           18*m2*Q2h*pow(U,2) - 4*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-1)*pow(4*m2 + Q2e,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*S*U*pow(Q2h,2) + 4*pow(Q2h,2)*pow(S,2) + 4*pow(Q2h,2)*pow(U,2)) \
		+ pow(l1k,-3)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-4*m4*Q2h*S*U + 4*m4*S*pow(Q2h,2) + 4*m4*U*pow(Q2h,2) +
		           16*m2*S*U*pow(Q2h,2) - 2*m4*pow(Q2h,3) - 16*m2*S*pow(Q2h,3) -
		           16*m2*U*pow(Q2h,3) + 4*S*U*pow(Q2h,3) + 8*m2*pow(Q2h,4) -
		           4*S*pow(Q2h,4) - 4*U*pow(Q2h,4) + 2*pow(Q2h,5) -
		           2*m4*Q2h*pow(S,2) + 8*m2*pow(Q2h,2)*pow(S,2) +
		           2*pow(Q2h,3)*pow(S,2) - 2*m4*Q2h*pow(U,2) +
		           8*m2*pow(Q2h,2)*pow(U,2) + 2*pow(Q2h,3)*pow(U,2)) +
		        (-384*Q2e*S*pow(l1k,3) + 240*Q2e*pow(l1k,4) - 192*S*pow(l1k,4) +
		           96*pow(l1k,5) - 288*S*pow(l1k,2)*pow(Q2e,2) +
		           240*pow(l1k,3)*pow(Q2e,2) - 96*l1k*S*pow(Q2e,3) +
		           120*pow(l1k,2)*pow(Q2e,3) + 30*l1k*pow(Q2e,4) -
		           12*S*pow(Q2e,4) + 3*pow(Q2e,5) + 144*Q2e*pow(l1k,2)*pow(S,2) +
		           96*pow(l1k,3)*pow(S,2) + 72*l1k*pow(Q2e,2)*pow(S,2) +
		           12*pow(Q2e,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-2) +
		        (32*l1k*m2*M2 + 16*m2*M2*Q2e + 124*l1k*Q2e*Q2h +
		           16*l1k*m2*S - 92*l1k*Q2e*S + 8*m2*Q2e*S + 144*l1k*Q2h*S -
		           180*Q2e*Q2h*S + 16*l1k*m2*U - 32*l1k*Q2e*U + 16*m2*Q2e*U -
		           8*l1k*Q2h*U - 4*Q2e*Q2h*U + 16*l1k*S*U - 16*m2*S*U +
		           24*Q2e*S*U + 8*Q2h*S*U + 64*Q2e*pow(l1k,2) -
		           120*Q2h*pow(l1k,2) - 152*S*pow(l1k,2) - 16*U*pow(l1k,2) +
		           96*pow(l1k,3) + 32*l1k*pow(Q2e,2) - 50*Q2h*pow(Q2e,2) -
		           38*S*pow(Q2e,2) - 20*U*pow(Q2e,2) + 12*pow(Q2e,3) -
		           212*l1k*pow(Q2h,2) + 260*Q2e*pow(Q2h,2) + 328*S*pow(Q2h,2) -
		           4*U*pow(Q2h,2) - 326*pow(Q2h,3) + 48*l1k*pow(S,2) -
		           16*m2*pow(S,2) + 28*Q2e*pow(S,2) - 40*Q2h*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-12*Q2e*S*pow(Q2h,2) + 8*S*U*pow(Q2h,2) + 4*Q2e*pow(Q2h,3) -
		           4*U*pow(Q2h,3) + 8*Q2e*Q2h*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*(-8*m2*Q2e*S*U - 8*Q2h*S*pow(Q2e,2) +
		           4*m2*U*pow(Q2e,2) + 8*S*U*pow(Q2e,2) + 3*Q2h*pow(Q2e,3) -
		           7*S*pow(Q2e,3) - 4*U*pow(Q2e,3) + 2*pow(Q2e,4) +
		           4*Q2e*S*pow(Q2h,2) - pow(Q2e,2)*pow(Q2h,2) - Q2e*pow(Q2h,3) +
		           4*S*pow(Q2h,3) - pow(Q2h,4) + 4*Q2e*Q2h*pow(S,2) +
		           6*pow(Q2e,2)*pow(S,2) - 4*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (4*S*pow(Q2h,4) - pow(Q2h,5) - 4*pow(Q2h,3)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1) +
		        (384*Q2h*U*pow(l1k,3) + 240*Q2h*pow(l1k,4) - 192*U*pow(l1k,4) -
		           96*pow(l1k,5) - 288*U*pow(l1k,2)*pow(Q2h,2) -
		           240*pow(l1k,3)*pow(Q2h,2) + 96*l1k*U*pow(Q2h,3) +
		           120*pow(l1k,2)*pow(Q2h,3) - 30*l1k*pow(Q2h,4) -
		           12*U*pow(Q2h,4) + 3*pow(Q2h,5) + 144*Q2h*pow(l1k,2)*pow(U,2) -
		           96*pow(l1k,3)*pow(U,2) - 72*l1k*pow(Q2h,2)*pow(U,2) +
		           12*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (-32*l1k*m2*M2 + 16*l1k*M2*Q2h - 36*l1k*Q2e*Q2h -
		           16*l1k*m2*S + 16*l1k*Q2e*S + 8*m2*Q2e*S + 32*l1k*Q2h*S -
		           8*Q2e*Q2h*S - 16*l1k*m2*U + 36*l1k*Q2e*U - 48*l1k*Q2h*U -
		           46*Q2e*Q2h*U - 16*l1k*S*U - 16*m2*S*U + 16*Q2e*S*U +
		           16*Q2h*S*U - 16*M2*pow(l1k,2) + 16*Q2e*pow(l1k,2) -
		           56*Q2h*pow(l1k,2) - 24*S*pow(l1k,2) + 32*U*pow(l1k,2) +
		           32*pow(l1k,3) - 8*l1k*pow(Q2e,2) + 10*Q2h*pow(Q2e,2) -
		           8*S*pow(Q2e,2) - 14*U*pow(Q2e,2) + 4*pow(Q2e,3) +
		           52*l1k*pow(Q2h,2) - 4*M2*pow(Q2h,2) + 16*Q2e*pow(Q2h,2) -
		           14*S*pow(Q2h,2) + 56*U*pow(Q2h,2) - 34*pow(Q2h,3) +
		           16*l1k*pow(U,2) - 16*m2*pow(U,2) + 20*Q2e*pow(U,2) -
		           16*Q2h*pow(U,2))*pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) +
		           pow(Q2h,2),-1) + pow(l1k,-1)*
		         (6*Q2e*U*pow(Q2h,2) - 4*S*U*pow(Q2h,2) - 2*Q2e*pow(Q2h,3) +
		           2*S*pow(Q2h,3) - 4*Q2e*Q2h*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (16*m2*Q2e*S*U - 8*m2*S*pow(Q2e,2) + 16*Q2h*U*pow(Q2e,2) -
		           16*S*U*pow(Q2e,2) - 6*Q2h*pow(Q2e,3) + 8*S*pow(Q2e,3) +
		           14*U*pow(Q2e,3) - 4*pow(Q2e,4) - 8*Q2e*U*pow(Q2h,2) +
		           2*pow(Q2e,2)*pow(Q2h,2) + 2*Q2e*pow(Q2h,3) - 8*U*pow(Q2h,3) +
		           2*pow(Q2h,4) - 8*Q2e*Q2h*pow(U,2) - 12*pow(Q2e,2)*pow(U,2) +
		           8*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (4*U*pow(Q2h,4) - pow(Q2h,5) - 4*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        (32*l1k*m2*Q2e*Q2h - 32*m4*Q2e*Q2h - 32*l1k*m2*Q2e*S -
		           32*l1k*m2*Q2h*S - 120*l1k*Q2e*Q2h*S + 112*m2*Q2e*Q2h*S +
		           32*l1k*m2*Q2e*U + 32*l1k*m2*Q2h*U + 72*l1k*Q2e*Q2h*U +
		           80*m2*Q2e*Q2h*U + 128*m4*S*U - 32*m2*Q2e*S*U -
		           16*Q2e*Q2h*S*U + 32*m2*Q2e*pow(l1k,2) +
		           32*m2*Q2h*pow(l1k,2) + 64*Q2e*Q2h*pow(l1k,2) +
		           48*Q2e*S*pow(l1k,2) - 48*Q2h*S*pow(l1k,2) +
		           16*Q2e*U*pow(l1k,2) + 48*Q2h*U*pow(l1k,2) + 64*Q2h*pow(l1k,3) +
		           32*S*pow(l1k,3) - 32*U*pow(l1k,3) - 32*pow(l1k,4) -
		           96*l1k*Q2h*pow(Q2e,2) + 40*l1k*S*pow(Q2e,2) +
		           148*Q2h*S*pow(Q2e,2) + 32*Q2h*U*pow(Q2e,2) -
		           28*pow(l1k,2)*pow(Q2e,2) + 24*l1k*pow(Q2e,3) +
		           100*Q2h*pow(Q2e,3) - 36*S*pow(Q2e,3) - 17*pow(Q2e,4) -
		           32*l1k*m2*pow(Q2h,2) - 32*m4*pow(Q2h,2) +
		           112*l1k*Q2e*pow(Q2h,2) - 96*m2*Q2e*pow(Q2h,2) +
		           32*l1k*S*pow(Q2h,2) + 16*m2*S*pow(Q2h,2) -
		           112*Q2e*S*pow(Q2h,2) - 24*l1k*U*pow(Q2h,2) -
		           16*m2*U*pow(Q2h,2) - 28*Q2e*U*pow(Q2h,2) -
		           52*pow(l1k,2)*pow(Q2h,2) - 245*pow(Q2e,2)*pow(Q2h,2) -
		           40*l1k*pow(Q2h,3) + 210*Q2e*pow(Q2h,3) + 80*S*pow(Q2h,3) +
		           12*U*pow(Q2h,3) - 88*pow(Q2h,4) + 64*m4*pow(S,2) -
		           64*m2*Q2e*pow(S,2) + 16*l1k*Q2h*pow(S,2) +
		           16*m2*Q2h*pow(S,2) - 20*Q2e*Q2h*pow(S,2) -
		           16*pow(l1k,2)*pow(S,2) - 12*pow(Q2e,2)*pow(S,2) -
		           8*pow(Q2h,2)*pow(S,2) + 64*m4*pow(U,2) - 32*l1k*Q2e*pow(U,2) +
		           16*l1k*Q2h*pow(U,2) + 16*m2*Q2h*pow(U,2) +
		           28*Q2e*Q2h*pow(U,2) - 16*pow(l1k,2)*pow(U,2) -
		           16*pow(Q2e,2)*pow(U,2) - 4*pow(Q2h,2)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*(-160*m4*Q2e*Q2h*S - 96*m4*Q2e*Q2h*U +
		           128*m4*Q2e*S*U + 96*m2*Q2e*Q2h*S*U + 64*m4*Q2e*pow(Q2h,2) +
		           32*m4*S*pow(Q2h,2) - 160*m2*Q2e*S*pow(Q2h,2) -
		           32*m4*U*pow(Q2h,2) - 96*m2*Q2e*U*pow(Q2h,2) +
		           32*Q2e*S*U*pow(Q2h,2) + 80*m2*Q2e*pow(Q2h,3) +
		           16*m2*S*pow(Q2h,3) - 48*Q2e*S*pow(Q2h,3) -
		           16*m2*U*pow(Q2h,3) - 32*Q2e*U*pow(Q2h,3) + 24*Q2e*pow(Q2h,4) +
		           96*m4*Q2e*pow(S,2) - 32*m4*Q2h*pow(S,2) +
		           72*m2*Q2e*Q2h*pow(S,2) - 8*m2*pow(Q2h,2)*pow(S,2) +
		           24*Q2e*pow(Q2h,2)*pow(S,2) + 32*m4*Q2e*pow(U,2) +
		           32*m4*Q2h*pow(U,2) + 24*m2*Q2e*Q2h*pow(U,2) +
		           8*m2*pow(Q2h,2)*pow(U,2) + 8*Q2e*pow(Q2h,2)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-2)*(-48*m4*Q2e*Q2h*S*U + 64*Q2e*Q2h*S*pow(m,6) +
		           64*Q2e*Q2h*U*pow(m,6) - 64*Q2e*S*U*pow(m,6) +
		           64*Q2h*S*U*pow(m,6) + 64*m4*Q2e*S*pow(Q2h,2) +
		           64*m4*Q2e*U*pow(Q2h,2) + 16*m4*S*U*pow(Q2h,2) -
		           48*m2*Q2e*S*U*pow(Q2h,2) - 32*Q2e*pow(m,6)*pow(Q2h,2) -
		           64*S*pow(m,6)*pow(Q2h,2) - 64*U*pow(m,6)*pow(Q2h,2) -
		           40*m4*Q2e*pow(Q2h,3) - 32*m4*S*pow(Q2h,3) +
		           56*m2*Q2e*S*pow(Q2h,3) - 32*m4*U*pow(Q2h,3) +
		           48*m2*Q2e*U*pow(Q2h,3) + 8*m2*S*U*pow(Q2h,3) -
		           8*Q2e*S*U*pow(Q2h,3) + 32*pow(m,6)*pow(Q2h,3) +
		           24*m4*pow(Q2h,4) - 28*m2*Q2e*pow(Q2h,4) -
		           16*m2*S*pow(Q2h,4) + 8*Q2e*S*pow(Q2h,4) - 8*m2*U*pow(Q2h,4) +
		           8*Q2e*U*pow(Q2h,4) + 8*m2*pow(Q2h,5) - 4*Q2e*pow(Q2h,5) -
		           24*m4*Q2e*Q2h*pow(S,2) - 32*Q2e*pow(m,6)*pow(S,2) +
		           32*Q2h*pow(m,6)*pow(S,2) + 8*m4*pow(Q2h,2)*pow(S,2) -
		           28*m2*Q2e*pow(Q2h,2)*pow(S,2) + 8*m2*pow(Q2h,3)*pow(S,2) -
		           4*Q2e*pow(Q2h,3)*pow(S,2) - 24*m4*Q2e*Q2h*pow(U,2) -
		           32*Q2e*pow(m,6)*pow(U,2) + 32*Q2h*pow(m,6)*pow(U,2) +
		           8*m4*pow(Q2h,2)*pow(U,2) - 20*m2*Q2e*pow(Q2h,2)*pow(U,2) -
		           4*Q2e*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-3)*(16*m4*Q2e*S*U*pow(Q2h,2) - 16*m4*Q2e*S*pow(Q2h,3) -
		           16*m4*Q2e*U*pow(Q2h,3) - 16*m4*S*U*pow(Q2h,3) +
		           4*m2*Q2e*S*U*pow(Q2h,3) + 8*m4*Q2e*pow(Q2h,4) +
		           16*m4*S*pow(Q2h,4) - 4*m2*Q2e*S*pow(Q2h,4) +
		           16*m4*U*pow(Q2h,4) - 4*m2*Q2e*U*pow(Q2h,4) -
		           4*m2*S*U*pow(Q2h,4) - 8*m4*pow(Q2h,5) +
		           2*m2*Q2e*pow(Q2h,5) + 4*m2*S*pow(Q2h,5) +
		           4*m2*U*pow(Q2h,5) - 2*m2*pow(Q2h,6) +
		           8*m4*Q2e*pow(Q2h,2)*pow(S,2) - 8*m4*pow(Q2h,3)*pow(S,2) +
		           2*m2*Q2e*pow(Q2h,3)*pow(S,2) - 2*m2*pow(Q2h,4)*pow(S,2) +
		           8*m4*Q2e*pow(Q2h,2)*pow(U,2) - 8*m4*pow(Q2h,3)*pow(U,2) +
		           2*m2*Q2e*pow(Q2h,3)*pow(U,2) - 2*m2*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        (128*m2*Q2e*S*pow(l1k,3) - 352*Q2e*Q2h*S*pow(l1k,3) -
		           64*m2*Q2e*pow(l1k,4) + 256*Q2e*Q2h*pow(l1k,4) +
		           192*Q2e*S*pow(l1k,4) + 192*Q2h*S*pow(l1k,4) -
		           128*Q2e*pow(l1k,5) - 128*Q2h*pow(l1k,5) - 128*S*pow(l1k,5) +
		           64*pow(l1k,6) - 216*l1k*Q2h*S*pow(Q2e,3) -
		           168*l1k*Q2h*pow(Q2e,4) + 104*l1k*S*pow(Q2e,4) -
		           148*Q2h*S*pow(Q2e,4) + 44*l1k*pow(Q2e,5) - 100*Q2h*pow(Q2e,5) +
		           36*S*pow(Q2e,5) + 17*pow(Q2e,6) -
		           320*Q2e*S*pow(l1k,2)*pow(Q2h,2) +
		           176*Q2e*pow(l1k,3)*pow(Q2h,2) - 128*S*pow(l1k,3)*pow(Q2h,2) +
		           112*pow(l1k,4)*pow(Q2h,2) - 192*l1k*S*pow(Q2e,2)*pow(Q2h,2) +
		           132*l1k*pow(Q2e,3)*pow(Q2h,2) + 192*S*pow(Q2e,3)*pow(Q2h,2) +
		           229*pow(Q2e,4)*pow(Q2h,2) + 320*l1k*Q2e*S*pow(Q2h,3) +
		           180*Q2e*pow(l1k,2)*pow(Q2h,3) - 320*S*pow(l1k,2)*pow(Q2h,3) +
		           208*pow(l1k,3)*pow(Q2h,3) + 312*l1k*pow(Q2e,2)*pow(Q2h,3) -
		           80*S*pow(Q2e,2)*pow(Q2h,3) - 227*pow(Q2e,3)*pow(Q2h,3) -
		           324*l1k*Q2e*pow(Q2h,4) + 324*pow(l1k,2)*pow(Q2h,4) +
		           81*pow(Q2e,2)*pow(Q2h,4) - 64*m2*Q2e*pow(l1k,2)*pow(S,2) +
		           112*Q2e*Q2h*pow(l1k,2)*pow(S,2) - 64*Q2e*pow(l1k,3)*pow(S,2) -
		           64*Q2h*pow(l1k,3)*pow(S,2) + 64*pow(l1k,4)*pow(S,2) -
		           16*l1k*Q2h*pow(Q2e,2)*pow(S,2) + 48*l1k*pow(Q2e,3)*pow(S,2) -
		           28*Q2h*pow(Q2e,3)*pow(S,2) + 12*pow(Q2e,4)*pow(S,2) -
		           48*l1k*Q2e*pow(Q2h,2)*pow(S,2) +
		           32*pow(l1k,2)*pow(Q2h,2)*pow(S,2) +
		           16*pow(Q2e,2)*pow(Q2h,2)*pow(S,2))*
		         pow(4*l1k*Q2e + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2e,2),-1)*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        (-128*m2*Q2e*U*pow(l1k,3) + 32*Q2e*Q2h*U*pow(l1k,3) -
		           64*m2*Q2e*pow(l1k,4) + 64*Q2e*Q2h*pow(l1k,4) -
		           64*Q2e*U*pow(l1k,4) - 320*Q2h*U*pow(l1k,4) - 64*Q2e*pow(l1k,5) -
		           192*Q2h*pow(l1k,5) + 128*U*pow(l1k,5) + 64*pow(l1k,6) -
		           16*Q2e*pow(l1k,3)*pow(Q2h,2) + 320*U*pow(l1k,3)*pow(Q2h,2) +
		           240*pow(l1k,4)*pow(Q2h,2) - 32*l1k*U*pow(Q2e,2)*pow(Q2h,2) -
		           56*l1k*Q2e*U*pow(Q2h,3) + 4*Q2e*pow(l1k,2)*pow(Q2h,3) -
		           192*U*pow(l1k,2)*pow(Q2h,3) - 176*pow(l1k,3)*pow(Q2h,3) +
		           16*l1k*pow(Q2e,2)*pow(Q2h,3) - 16*U*pow(Q2e,2)*pow(Q2h,3) +
		           28*l1k*Q2e*pow(Q2h,4) + 72*l1k*U*pow(Q2h,4) +
		           28*Q2e*U*pow(Q2h,4) + 100*pow(l1k,2)*pow(Q2h,4) +
		           8*pow(Q2e,2)*pow(Q2h,4) - 40*l1k*pow(Q2h,5) -
		           15*Q2e*pow(Q2h,5) - 12*U*pow(Q2h,5) + 7*pow(Q2h,6) -
		           64*m2*Q2e*pow(l1k,2)*pow(U,2) -
		           16*Q2e*Q2h*pow(l1k,2)*pow(U,2) - 128*Q2h*pow(l1k,3)*pow(U,2) +
		           64*pow(l1k,4)*pow(U,2) + 16*l1k*Q2h*pow(Q2e,2)*pow(U,2) +
		           32*l1k*Q2e*pow(Q2h,2)*pow(U,2) +
		           96*pow(l1k,2)*pow(Q2h,2)*pow(U,2) +
		           8*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) - 32*l1k*pow(Q2h,3)*pow(U,2) -
		           12*Q2e*pow(Q2h,3)*pow(U,2) + 4*pow(Q2h,4)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1)) +
		     SI.C0_mMQm0m(-2*l1k + m2,Q2h,m2)*
		      (56*l1k*M2 - 8*l1k*Q2e - 16*M2*Q2e + 6*l1k*Q2h - 8*m2*Q2h -
		        32*M2*Q2h + 10*Q2e*Q2h + 12*l1k*S + 8*m2*S - 8*Q2e*S - 80*Q2h*S +
		        8*l1k*U + 8*m2*U - 18*Q2e*U + 16*Q2h*U + 8*S*U - 4*pow(l1k,2) +
		        4*pow(Q2e,2) + 26*pow(Q2h,2) + 32*pow(S,2) - 12*pow(U,2) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (48*m2*M2*Q2e + 96*m2*M2*Q2h + 16*m2*Q2e*Q2h +
		           48*M2*Q2e*Q2h - 8*m2*Q2e*S + 112*m2*Q2h*S + 76*Q2e*Q2h*S -
		           8*m2*Q2h*U + 26*Q2e*Q2h*U + 16*m2*S*U - 12*Q2e*S*U -
		           88*Q2h*S*U + 24*M2*pow(Q2e,2) - 2*Q2h*pow(Q2e,2) +
		           6*S*pow(Q2e,2) + 14*U*pow(Q2e,2) - 4*pow(Q2e,3) -
		           64*m2*pow(Q2h,2) + 48*M2*pow(Q2h,2) - 44*Q2e*pow(Q2h,2) +
		           316*S*pow(Q2h,2) + 148*U*pow(Q2h,2) - 196*pow(Q2h,3) -
		           32*m2*pow(S,2) - 32*Q2e*pow(S,2) - 120*Q2h*pow(S,2) -
		           12*Q2e*pow(U,2) - 20*Q2h*pow(U,2)) +
		        pow(l1k,-1)*(-8*m2*M2*Q2h + 4*M2*Q2e*Q2h + 8*m2*Q2h*U +
		           4*Q2e*Q2h*U - 8*m2*pow(Q2h,2) - 4*M2*pow(Q2h,2) +
		           8*U*pow(Q2h,2) - 4*pow(Q2h,3) + 8*m2*pow(U,2) - 4*Q2e*pow(U,2) -
		           4*Q2h*pow(U,2)) + pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*m4*Q2h*S + 80*m2*Q2h*S*U + 16*m4*pow(Q2h,2) -
		           32*m2*M2*pow(Q2h,2) - 256*m2*S*pow(Q2h,2) -
		           120*m2*U*pow(Q2h,2) + 196*S*U*pow(Q2h,2) + 152*m2*pow(Q2h,3) -
		           24*M2*pow(Q2h,3) - 406*S*pow(Q2h,3) - 258*U*pow(Q2h,3) +
		           238*pow(Q2h,4) + 96*m2*Q2h*pow(S,2) + 168*pow(Q2h,2)*pow(S,2) +
		           52*pow(Q2h,2)*pow(U,2)) +
		        pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-32*m4*Q2h*S*U + 48*m4*S*pow(Q2h,2) + 16*m4*U*pow(Q2h,2) -
		           120*m2*S*U*pow(Q2h,2) - 16*m4*pow(Q2h,3) +
		           172*m2*S*pow(Q2h,3) + 132*m2*U*pow(Q2h,3) -
		           140*S*U*pow(Q2h,3) - 92*m2*pow(Q2h,4) + 188*S*pow(Q2h,4) +
		           156*U*pow(Q2h,4) - 102*pow(Q2h,5) - 32*m4*Q2h*pow(S,2) -
		           80*m2*pow(Q2h,2)*pow(S,2) - 86*pow(Q2h,3)*pow(S,2) -
		           40*m2*pow(Q2h,2)*pow(U,2) - 54*pow(Q2h,3)*pow(U,2)) +
		        pow(l1k,-3)*pow(2*l1k + Q2e - Q2h,-1)*
		         (32*m2*S*U*pow(Q2h,3) - 32*m2*S*pow(Q2h,4) -
		           32*m2*U*pow(Q2h,4) + 32*S*U*pow(Q2h,4) + 16*m2*pow(Q2h,5) -
		           32*S*pow(Q2h,5) - 32*U*pow(Q2h,5) + 16*pow(Q2h,6) +
		           16*m2*pow(Q2h,3)*pow(S,2) + 16*pow(Q2h,4)*pow(S,2) +
		           16*m2*pow(Q2h,3)*pow(U,2) + 16*pow(Q2h,4)*pow(U,2)) +
		        (768*Q2h*U*pow(l1k,4) + 480*Q2h*pow(l1k,5) - 384*U*pow(l1k,5) -
		           192*pow(l1k,6) - 576*U*pow(l1k,3)*pow(Q2h,2) -
		           480*pow(l1k,4)*pow(Q2h,2) + 192*U*pow(l1k,2)*pow(Q2h,3) +
		           240*pow(l1k,3)*pow(Q2h,3) - 24*l1k*U*pow(Q2h,4) -
		           60*pow(l1k,2)*pow(Q2h,4) + 6*l1k*pow(Q2h,5) +
		           288*Q2h*pow(l1k,3)*pow(U,2) - 192*pow(l1k,4)*pow(U,2) -
		           144*pow(l1k,2)*pow(Q2h,2)*pow(U,2) + 24*l1k*pow(Q2h,3)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-2) +
		        (16*l1k*m2*Q2e*S - 16*l1k*Q2e*Q2h*S - 92*l1k*Q2e*Q2h*U -
		           32*l1k*m2*S*U + 24*l1k*Q2e*S*U + 16*m2*Q2e*S*U +
		           32*l1k*Q2h*S*U + 4*Q2e*Q2h*S*U - 64*m2*M2*pow(l1k,2) +
		           32*M2*Q2h*pow(l1k,2) - 72*Q2e*Q2h*pow(l1k,2) -
		           32*m2*S*pow(l1k,2) + 32*Q2e*S*pow(l1k,2) +
		           64*Q2h*S*pow(l1k,2) - 32*m2*U*pow(l1k,2) +
		           72*Q2e*U*pow(l1k,2) - 96*Q2h*U*pow(l1k,2) - 32*S*U*pow(l1k,2) -
		           32*M2*pow(l1k,3) + 32*Q2e*pow(l1k,3) - 112*Q2h*pow(l1k,3) -
		           48*S*pow(l1k,3) + 64*U*pow(l1k,3) + 64*pow(l1k,4) +
		           20*l1k*Q2h*pow(Q2e,2) - 12*l1k*S*pow(Q2e,2) -
		           8*m2*S*pow(Q2e,2) - 2*Q2h*S*pow(Q2e,2) - 28*l1k*U*pow(Q2e,2) +
		           16*Q2h*U*pow(Q2e,2) - 12*S*U*pow(Q2e,2) -
		           16*pow(l1k,2)*pow(Q2e,2) + 8*l1k*pow(Q2e,3) - 6*Q2h*pow(Q2e,3) +
		           6*S*pow(Q2e,3) + 14*U*pow(Q2e,3) - 4*pow(Q2e,4) -
		           8*l1k*M2*pow(Q2h,2) + 32*l1k*Q2e*pow(Q2h,2) -
		           28*l1k*S*pow(Q2h,2) + 112*l1k*U*pow(Q2h,2) +
		           4*Q2e*U*pow(Q2h,2) - 8*S*U*pow(Q2h,2) +
		           104*pow(l1k,2)*pow(Q2h,2) + 2*pow(Q2e,2)*pow(Q2h,2) -
		           68*l1k*pow(Q2h,3) - 2*Q2e*pow(Q2h,3) + 4*S*pow(Q2h,3) -
		           72*U*pow(Q2h,3) + 34*pow(Q2h,4) - 32*l1k*m2*pow(U,2) +
		           40*l1k*Q2e*pow(U,2) - 32*l1k*Q2h*pow(U,2) -
		           16*Q2e*Q2h*pow(U,2) + 32*pow(l1k,2)*pow(U,2) -
		           12*pow(Q2e,2)*pow(U,2) + 40*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-16*m2*S*U*pow(Q2e,2) - 16*Q2h*S*U*pow(Q2e,2) +
		           8*m2*S*pow(Q2e,3) + 8*Q2h*S*pow(Q2e,3) - 2*Q2h*U*pow(Q2e,3) +
		           12*S*U*pow(Q2e,3) + 2*Q2h*pow(Q2e,4) - 6*S*pow(Q2e,4) -
		           14*U*pow(Q2e,4) + 4*pow(Q2e,5) + 24*U*pow(Q2e,2)*pow(Q2h,2) -
		           8*pow(Q2e,3)*pow(Q2h,2) - 4*Q2h*pow(Q2e,2)*pow(U,2) +
		           12*pow(Q2e,3)*pow(U,2) - 16*Q2e*pow(Q2h,2)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (-8*S*U*pow(Q2e,5)*pow(Q2h,2) + 8*S*pow(Q2e,6)*pow(Q2h,2) +
		           8*S*U*pow(Q2e,4)*pow(Q2h,3) - 8*S*pow(Q2e,5)*pow(Q2h,3) +
		           8*U*pow(Q2e,5)*pow(Q2h,3) - 4*pow(Q2e,6)*pow(Q2h,3) -
		           8*U*pow(Q2e,4)*pow(Q2h,4) + 4*pow(Q2e,5)*pow(Q2h,4) -
		           4*Q2h*pow(Q2e,6)*pow(S,2) + 4*pow(Q2e,5)*pow(Q2h,2)*pow(S,2) -
		           4*pow(Q2e,4)*pow(Q2h,3)*pow(U,2) +
		           4*pow(Q2e,3)*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        (1024*m4*Q2e*Q2h*S*pow(l1k,3) + 1024*m4*Q2e*Q2h*U*pow(l1k,3) -
		           512*m4*Q2e*S*U*pow(l1k,3) - 512*m4*Q2h*S*U*pow(l1k,3) -
		           2176*m2*Q2e*Q2h*S*U*pow(l1k,3) +
		           2048*m2*Q2e*Q2h*U*pow(l1k,4) - 512*m2*Q2e*S*U*pow(l1k,4) -
		           1280*Q2e*Q2h*S*U*pow(l1k,4) + 2048*Q2e*Q2h*S*pow(l1k,5) +
		           6144*Q2e*Q2h*U*pow(l1k,5) - 1536*Q2e*S*U*pow(l1k,5) -
		           512*l1k*Q2e*Q2h*S*U*pow(m,6) -
		           1216*Q2h*S*U*pow(l1k,3)*pow(Q2e,2) +
		           1024*Q2h*S*pow(l1k,4)*pow(Q2e,2) +
		           4096*Q2h*U*pow(l1k,4)*pow(Q2e,2) -
		           896*S*U*pow(l1k,4)*pow(Q2e,2) -
		           64*Q2h*S*U*pow(l1k,2)*pow(Q2e,3) +
		           64*Q2h*S*pow(l1k,3)*pow(Q2e,3) +
		           576*Q2h*U*pow(l1k,3)*pow(Q2e,3) - 64*S*U*pow(l1k,3)*pow(Q2e,3) +
		           32*Q2h*U*pow(l1k,2)*pow(Q2e,4) -
		           1408*l1k*m4*Q2e*S*U*pow(Q2h,2) -
		           1024*m4*Q2e*S*pow(l1k,2)*pow(Q2h,2) +
		           512*m4*Q2e*U*pow(l1k,2)*pow(Q2h,2) +
		           512*m4*S*U*pow(l1k,2)*pow(Q2h,2) +
		           384*m2*Q2e*S*U*pow(l1k,2)*pow(Q2h,2) -
		           1024*m4*Q2e*pow(l1k,3)*pow(Q2h,2) +
		           3584*m2*Q2e*S*pow(l1k,3)*pow(Q2h,2) +
		           5632*m2*Q2e*U*pow(l1k,3)*pow(Q2h,2) +
		           384*m2*S*U*pow(l1k,3)*pow(Q2h,2) -
		           512*Q2e*S*U*pow(l1k,3)*pow(Q2h,2) -
		           1024*m2*Q2e*pow(l1k,4)*pow(Q2h,2) +
		           512*Q2e*S*pow(l1k,4)*pow(Q2h,2) +
		           4608*Q2e*U*pow(l1k,4)*pow(Q2h,2) -
		           4096*Q2e*pow(l1k,5)*pow(Q2h,2) +
		           1024*l1k*Q2e*S*pow(m,6)*pow(Q2h,2) +
		           1024*l1k*Q2e*U*pow(m,6)*pow(Q2h,2) -
		           512*l1k*S*U*pow(m,6)*pow(Q2h,2) -
		           384*S*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2) +
		           1408*S*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) +
		           4224*U*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,2) -
		           2560*pow(l1k,4)*pow(Q2e,2)*pow(Q2h,2) -
		           32*l1k*S*U*pow(Q2e,3)*pow(Q2h,2) +
		           32*S*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2) +
		           384*U*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,2) -
		           320*pow(l1k,3)*pow(Q2e,3)*pow(Q2h,2) +
		           16*l1k*S*pow(Q2e,4)*pow(Q2h,2) +
		           16*l1k*U*pow(Q2e,4)*pow(Q2h,2) + 8*S*U*pow(Q2e,4)*pow(Q2h,2) -
		           16*pow(l1k,2)*pow(Q2e,4)*pow(Q2h,2) -
		           8*S*pow(Q2e,5)*pow(Q2h,2) + 1920*l1k*m4*Q2e*S*pow(Q2h,3) +
		           2176*l1k*m4*Q2e*U*pow(Q2h,3) + 384*l1k*m4*S*U*pow(Q2h,3) -
		           1472*l1k*m2*Q2e*S*U*pow(Q2h,3) + 512*m4*Q2e*S*U*pow(Q2h,3) -
		           1536*m2*Q2e*S*pow(l1k,2)*pow(Q2h,3) -
		           512*m4*U*pow(l1k,2)*pow(Q2h,3) +
		           1280*m2*Q2e*U*pow(l1k,2)*pow(Q2h,3) +
		           512*m2*S*U*pow(l1k,2)*pow(Q2h,3) +
		           1280*Q2e*S*U*pow(l1k,2)*pow(Q2h,3) -
		           4096*m2*Q2e*pow(l1k,3)*pow(Q2h,3) -
		           512*m2*S*pow(l1k,3)*pow(Q2h,3) +
		           1664*Q2e*S*pow(l1k,3)*pow(Q2h,3) -
		           1536*m2*U*pow(l1k,3)*pow(Q2h,3) -
		           1152*Q2e*U*pow(l1k,3)*pow(Q2h,3) -
		           2048*Q2e*pow(l1k,4)*pow(Q2h,3) -
		           768*l1k*Q2e*pow(m,6)*pow(Q2h,3) -
		           32*l1k*S*U*pow(Q2e,2)*pow(Q2h,3) +
		           384*S*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3) +
		           640*U*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,3) -
		           2496*pow(l1k,3)*pow(Q2e,2)*pow(Q2h,3) +
		           32*l1k*S*pow(Q2e,3)*pow(Q2h,3) +
		           32*l1k*U*pow(Q2e,3)*pow(Q2h,3) -
		           192*pow(l1k,2)*pow(Q2e,3)*pow(Q2h,3) -
		           16*l1k*pow(Q2e,4)*pow(Q2h,3) - 8*U*pow(Q2e,4)*pow(Q2h,3) +
		           4*pow(Q2e,5)*pow(Q2h,3) - 1408*l1k*m4*Q2e*pow(Q2h,4) -
		           640*l1k*m4*S*pow(Q2h,4) + 2080*l1k*m2*Q2e*S*pow(Q2h,4) -
		           640*m4*Q2e*S*pow(Q2h,4) - 896*l1k*m4*U*pow(Q2h,4) +
		           1888*l1k*m2*Q2e*U*pow(Q2h,4) - 512*m4*Q2e*U*pow(Q2h,4) +
		           384*l1k*m2*S*U*pow(Q2h,4) - 128*m4*S*U*pow(Q2h,4) -
		           416*l1k*Q2e*S*U*pow(Q2h,4) + 448*m2*Q2e*S*U*pow(Q2h,4) +
		           256*m4*pow(l1k,2)*pow(Q2h,4) -
		           128*m2*Q2e*pow(l1k,2)*pow(Q2h,4) -
		           384*m2*S*pow(l1k,2)*pow(Q2h,4) -
		           1984*Q2e*S*pow(l1k,2)*pow(Q2h,4) -
		           1920*m2*U*pow(l1k,2)*pow(Q2h,4) -
		           1856*Q2e*U*pow(l1k,2)*pow(Q2h,4) +
		           1024*m2*pow(l1k,3)*pow(Q2h,4) + 192*Q2e*pow(l1k,3)*pow(Q2h,4) +
		           256*l1k*pow(m,6)*pow(Q2h,4) + 32*l1k*S*pow(Q2e,2)*pow(Q2h,4) +
		           32*l1k*U*pow(Q2e,2)*pow(Q2h,4) -
		           320*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,4) -
		           16*l1k*pow(Q2e,3)*pow(Q2h,4) + 640*l1k*m4*pow(Q2h,5) -
		           1264*l1k*m2*Q2e*pow(Q2h,5) + 320*m4*Q2e*pow(Q2h,5) -
		           608*l1k*m2*S*pow(Q2h,5) + 256*m4*S*pow(Q2h,5) +
		           544*l1k*Q2e*S*pow(Q2h,5) - 512*m2*Q2e*S*pow(Q2h,5) -
		           672*l1k*m2*U*pow(Q2h,5) + 128*m4*U*pow(Q2h,5) +
		           416*l1k*Q2e*U*pow(Q2h,5) - 448*m2*Q2e*U*pow(Q2h,5) -
		           128*m2*S*U*pow(Q2h,5) + 80*Q2e*S*U*pow(Q2h,5) +
		           1024*m2*pow(l1k,2)*pow(Q2h,5) +
		           1312*Q2e*pow(l1k,2)*pow(Q2h,5) - 16*l1k*pow(Q2e,2)*pow(Q2h,5) +
		           464*l1k*m2*pow(Q2h,6) - 128*m4*pow(Q2h,6) -
		           272*l1k*Q2e*pow(Q2h,6) + 256*m2*Q2e*pow(Q2h,6) +
		           192*m2*S*pow(Q2h,6) - 80*Q2e*S*pow(Q2h,6) +
		           128*m2*U*pow(Q2h,6) - 80*Q2e*U*pow(Q2h,6) - 96*m2*pow(Q2h,7) +
		           40*Q2e*pow(Q2h,7) + 512*m4*Q2e*Q2h*pow(l1k,2)*pow(S,2) -
		           256*m4*Q2e*pow(l1k,3)*pow(S,2) -
		           256*m4*Q2h*pow(l1k,3)*pow(S,2) -
		           960*m2*Q2e*Q2h*pow(l1k,3)*pow(S,2) +
		           256*m2*Q2e*pow(l1k,4)*pow(S,2) +
		           384*Q2e*Q2h*pow(l1k,4)*pow(S,2) - 256*Q2e*pow(l1k,5)*pow(S,2) -
		           256*l1k*Q2e*Q2h*pow(m,6)*pow(S,2) -
		           128*Q2h*pow(l1k,3)*pow(Q2e,2)*pow(S,2) -
		           64*pow(l1k,4)*pow(Q2e,2)*pow(S,2) +
		           16*Q2h*pow(l1k,2)*pow(Q2e,3)*pow(S,2) -
		           8*l1k*Q2h*pow(Q2e,4)*pow(S,2) + 4*Q2h*pow(Q2e,5)*pow(S,2) -
		           576*l1k*m4*Q2e*pow(Q2h,2)*pow(S,2) +
		           256*m4*pow(l1k,2)*pow(Q2h,2)*pow(S,2) +
		           1024*m2*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(S,2) +
		           64*m2*pow(l1k,3)*pow(Q2h,2)*pow(S,2) -
		           896*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(S,2) -
		           256*l1k*pow(m,6)*pow(Q2h,2)*pow(S,2) -
		           16*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) -
		           8*l1k*pow(Q2e,3)*pow(Q2h,2)*pow(S,2) +
		           64*l1k*m4*pow(Q2h,3)*pow(S,2) -
		           864*l1k*m2*Q2e*pow(Q2h,3)*pow(S,2) +
		           320*m4*Q2e*pow(Q2h,3)*pow(S,2) -
		           64*m2*pow(l1k,2)*pow(Q2h,3)*pow(S,2) +
		           768*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(S,2) -
		           16*l1k*pow(Q2e,2)*pow(Q2h,3)*pow(S,2) +
		           192*l1k*m2*pow(Q2h,4)*pow(S,2) - 128*m4*pow(Q2h,4)*pow(S,2) -
		           272*l1k*Q2e*pow(Q2h,4)*pow(S,2) +
		           256*m2*Q2e*pow(Q2h,4)*pow(S,2) - 96*m2*pow(Q2h,5)*pow(S,2) +
		           40*Q2e*pow(Q2h,5)*pow(S,2) -
		           512*m4*Q2e*Q2h*pow(l1k,2)*pow(U,2) -
		           256*m4*Q2e*pow(l1k,3)*pow(U,2) -
		           256*m4*Q2h*pow(l1k,3)*pow(U,2) -
		           1984*m2*Q2e*Q2h*pow(l1k,3)*pow(U,2) -
		           768*m2*Q2e*pow(l1k,4)*pow(U,2) -
		           2176*Q2e*Q2h*pow(l1k,4)*pow(U,2) -
		           2304*Q2e*pow(l1k,5)*pow(U,2) -
		           256*l1k*Q2e*Q2h*pow(m,6)*pow(U,2) -
		           1792*Q2h*pow(l1k,3)*pow(Q2e,2)*pow(U,2) -
		           1600*pow(l1k,4)*pow(Q2e,2)*pow(U,2) -
		           176*Q2h*pow(l1k,2)*pow(Q2e,3)*pow(U,2) -
		           256*pow(l1k,3)*pow(Q2e,3)*pow(U,2) -
		           8*l1k*Q2h*pow(Q2e,4)*pow(U,2) -
		           16*pow(l1k,2)*pow(Q2e,4)*pow(U,2) -
		           832*l1k*m4*Q2e*pow(Q2h,2)*pow(U,2) +
		           256*m4*pow(l1k,2)*pow(Q2h,2)*pow(U,2) -
		           896*m2*Q2e*pow(l1k,2)*pow(Q2h,2)*pow(U,2) +
		           576*m2*pow(l1k,3)*pow(Q2h,2)*pow(U,2) +
		           640*Q2e*pow(l1k,3)*pow(Q2h,2)*pow(U,2) -
		           256*l1k*pow(m,6)*pow(Q2h,2)*pow(U,2) -
		           320*pow(l1k,2)*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) -
		           8*l1k*pow(Q2e,3)*pow(Q2h,2)*pow(U,2) +
		           320*l1k*m4*pow(Q2h,3)*pow(U,2) -
		           672*l1k*m2*Q2e*pow(Q2h,3)*pow(U,2) +
		           192*m4*Q2e*pow(Q2h,3)*pow(U,2) +
		           832*m2*pow(l1k,2)*pow(Q2h,3)*pow(U,2) +
		           640*Q2e*pow(l1k,2)*pow(Q2h,3)*pow(U,2) -
		           16*l1k*pow(Q2e,2)*pow(Q2h,3)*pow(U,2) +
		           4*pow(Q2e,3)*pow(Q2h,3)*pow(U,2) +
		           256*l1k*m2*pow(Q2h,4)*pow(U,2) -
		           144*l1k*Q2e*pow(Q2h,4)*pow(U,2) +
		           192*m2*Q2e*pow(Q2h,4)*pow(U,2) - 32*m2*pow(Q2h,5)*pow(U,2) +
		           40*Q2e*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        pow(l1k,-1)*(-128*Q2e*S*U*pow(m,6)*pow(Q2h,3) -
		           128*m4*Q2e*S*U*pow(Q2h,4) + 128*Q2e*S*pow(m,6)*pow(Q2h,4) +
		           128*Q2e*U*pow(m,6)*pow(Q2h,4) + 128*S*U*pow(m,6)*pow(Q2h,4) +
		           128*m4*Q2e*S*pow(Q2h,5) + 128*m4*Q2e*U*pow(Q2h,5) +
		           128*m4*S*U*pow(Q2h,5) - 40*m2*Q2e*S*U*pow(Q2h,5) -
		           64*Q2e*pow(m,6)*pow(Q2h,5) - 128*S*pow(m,6)*pow(Q2h,5) -
		           128*U*pow(m,6)*pow(Q2h,5) - 64*m4*Q2e*pow(Q2h,6) -
		           128*m4*S*pow(Q2h,6) + 40*m2*Q2e*S*pow(Q2h,6) -
		           128*m4*U*pow(Q2h,6) + 40*m2*Q2e*U*pow(Q2h,6) +
		           40*m2*S*U*pow(Q2h,6) + 64*pow(m,6)*pow(Q2h,6) +
		           64*m4*pow(Q2h,7) - 20*m2*Q2e*pow(Q2h,7) -
		           40*m2*S*pow(Q2h,7) - 40*m2*U*pow(Q2h,7) + 20*m2*pow(Q2h,8) -
		           64*Q2e*pow(m,6)*pow(Q2h,3)*pow(S,2) -
		           64*m4*Q2e*pow(Q2h,4)*pow(S,2) +
		           64*pow(m,6)*pow(Q2h,4)*pow(S,2) + 64*m4*pow(Q2h,5)*pow(S,2) -
		           20*m2*Q2e*pow(Q2h,5)*pow(S,2) + 20*m2*pow(Q2h,6)*pow(S,2) -
		           64*Q2e*pow(m,6)*pow(Q2h,3)*pow(U,2) -
		           64*m4*Q2e*pow(Q2h,4)*pow(U,2) +
		           64*pow(m,6)*pow(Q2h,4)*pow(U,2) + 64*m4*pow(Q2h,5)*pow(U,2) -
		           20*m2*Q2e*pow(Q2h,5)*pow(U,2) + 20*m2*pow(Q2h,6)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-2) +
		        (-128*l1k*m4*M2*Q2e - 128*l1k*m4*M2*Q2h -
		           32*l1k*m4*Q2e*Q2h - 224*l1k*m2*M2*Q2e*Q2h -
		           192*l1k*m4*Q2h*S - 376*l1k*m2*Q2e*Q2h*S + 144*m4*Q2e*Q2h*S -
		           192*l1k*m4*Q2h*U - 264*l1k*m2*Q2e*Q2h*U - 16*m4*Q2e*Q2h*U +
		           128*l1k*m4*S*U + 64*l1k*m2*Q2e*S*U + 448*l1k*m2*Q2h*S*U -
		           64*m4*Q2h*S*U + 256*l1k*Q2e*Q2h*S*U - 248*m2*Q2e*Q2h*S*U +
		           128*m2*M2*Q2e*pow(l1k,2) + 192*M2*Q2e*Q2h*pow(l1k,2) -
		           32*m2*Q2h*S*pow(l1k,2) + 304*Q2e*Q2h*S*pow(l1k,2) +
		           64*m2*Q2e*U*pow(l1k,2) - 480*m2*Q2h*U*pow(l1k,2) -
		           304*Q2e*Q2h*U*pow(l1k,2) + 128*m2*S*U*pow(l1k,2) +
		           32*Q2e*S*U*pow(l1k,2) + 512*Q2h*S*U*pow(l1k,2) +
		           32*m2*Q2e*pow(l1k,3) - 128*M2*Q2e*pow(l1k,3) +
		           32*m2*Q2h*pow(l1k,3) - 512*Q2h*S*pow(l1k,3) +
		           32*Q2e*U*pow(l1k,3) - 1440*Q2h*U*pow(l1k,3) +
		           384*S*U*pow(l1k,3) + 32*Q2e*pow(l1k,4) + 64*Q2h*pow(l1k,4) -
		           64*U*pow(l1k,4) - 32*pow(l1k,5) - 32*l1k*M2*Q2h*pow(Q2e,2) -
		           24*l1k*Q2h*S*pow(Q2e,2) - 72*l1k*Q2h*U*pow(Q2e,2) -
		           16*Q2h*S*U*pow(Q2e,2) - 32*M2*pow(l1k,2)*pow(Q2e,2) +
		           8*M2*Q2h*pow(Q2e,3) + 12*Q2h*S*pow(Q2e,3) +
		           224*l1k*m4*pow(Q2h,2) + 32*l1k*m2*M2*pow(Q2h,2) +
		           280*l1k*m2*Q2e*pow(Q2h,2) - 48*m4*Q2e*pow(Q2h,2) -
		           160*l1k*M2*Q2e*pow(Q2h,2) + 128*m2*M2*Q2e*pow(Q2h,2) -
		           616*l1k*m2*S*pow(Q2h,2) + 48*m4*S*pow(Q2h,2) -
		           944*l1k*Q2e*S*pow(Q2h,2) + 848*m2*Q2e*S*pow(Q2h,2) -
		           920*l1k*m2*U*pow(Q2h,2) - 48*m4*U*pow(Q2h,2) -
		           592*l1k*Q2e*U*pow(Q2h,2) + 240*m2*Q2e*U*pow(Q2h,2) +
		           384*l1k*S*U*pow(Q2h,2) - 120*m2*S*U*pow(Q2h,2) -
		           456*Q2e*S*U*pow(Q2h,2) + 224*m2*pow(l1k,2)*pow(Q2h,2) -
		           16*Q2e*pow(l1k,2)*pow(Q2h,2) - 384*S*pow(l1k,2)*pow(Q2h,2) -
		           1968*U*pow(l1k,2)*pow(Q2h,2) + 976*pow(l1k,3)*pow(Q2h,2) +
		           48*l1k*pow(Q2e,2)*pow(Q2h,2) + 68*S*pow(Q2e,2)*pow(Q2h,2) -
		           24*U*pow(Q2e,2)*pow(Q2h,2) - 8*pow(Q2e,3)*pow(Q2h,2) +
		           648*l1k*m2*pow(Q2h,3) + 16*m4*pow(Q2h,3) -
		           32*m2*M2*pow(Q2h,3) + 670*l1k*Q2e*pow(Q2h,3) -
		           428*m2*Q2e*pow(Q2h,3) + 48*M2*Q2e*pow(Q2h,3) -
		           608*l1k*S*pow(Q2h,3) + 176*m2*S*pow(Q2h,3) +
		           940*Q2e*S*pow(Q2h,3) - 648*l1k*U*pow(Q2h,3) +
		           144*m2*U*pow(Q2h,3) + 636*Q2e*U*pow(Q2h,3) -
		           128*S*U*pow(Q2h,3) + 1048*pow(l1k,2)*pow(Q2h,3) -
		           16*pow(Q2e,2)*pow(Q2h,3) + 450*l1k*pow(Q2h,4) -
		           100*m2*pow(Q2h,4) - 564*Q2e*pow(Q2h,4) + 192*S*pow(Q2h,4) +
		           112*U*pow(Q2h,4) - 88*pow(Q2h,5) + 64*l1k*m4*pow(S,2) +
		           128*l1k*m2*Q2e*pow(S,2) - 32*m4*Q2e*pow(S,2) +
		           208*l1k*m2*Q2h*pow(S,2) - 96*m4*Q2h*pow(S,2) +
		           368*l1k*Q2e*Q2h*pow(S,2) - 344*m2*Q2e*Q2h*pow(S,2) -
		           64*m2*pow(l1k,2)*pow(S,2) - 144*Q2e*pow(l1k,2)*pow(S,2) -
		           64*Q2h*pow(l1k,2)*pow(S,2) + 64*pow(l1k,3)*pow(S,2) +
		           8*l1k*pow(Q2e,2)*pow(S,2) - 28*Q2h*pow(Q2e,2)*pow(S,2) -
		           4*pow(Q2e,3)*pow(S,2) + 192*l1k*pow(Q2h,2)*pow(S,2) -
		           120*m2*pow(Q2h,2)*pow(S,2) - 404*Q2e*pow(Q2h,2)*pow(S,2) -
		           96*pow(Q2h,3)*pow(S,2) + 64*l1k*m4*pow(U,2) +
		           96*l1k*m2*Q2e*pow(U,2) + 32*m4*Q2e*pow(U,2) +
		           336*l1k*m2*Q2h*pow(U,2) + 32*m4*Q2h*pow(U,2) +
		           136*l1k*Q2e*Q2h*pow(U,2) + 48*m2*Q2e*Q2h*pow(U,2) +
		           192*m2*pow(l1k,2)*pow(U,2) + 176*Q2e*pow(l1k,2)*pow(U,2) +
		           864*Q2h*pow(l1k,2)*pow(U,2) + 544*pow(l1k,3)*pow(U,2) +
		           40*l1k*pow(Q2e,2)*pow(U,2) + 28*Q2h*pow(Q2e,2)*pow(U,2) +
		           248*l1k*pow(Q2h,2)*pow(U,2) - 16*m2*pow(Q2h,2)*pow(U,2) -
		           144*Q2e*pow(Q2h,2)*pow(U,2) - 24*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*(176*m4*Q2e*Q2h*S*U - 32*Q2e*Q2h*S*pow(m,6) -
		           32*Q2e*Q2h*U*pow(m,6) + 128*Q2h*S*U*pow(m,6) -
		           32*m4*M2*Q2e*pow(Q2h,2) - 368*m4*Q2e*S*pow(Q2h,2) -
		           240*m4*Q2e*U*pow(Q2h,2) + 176*m4*S*U*pow(Q2h,2) +
		           444*m2*Q2e*S*U*pow(Q2h,2) + 32*Q2e*pow(m,6)*pow(Q2h,2) -
		           96*S*pow(m,6)*pow(Q2h,2) - 96*U*pow(m,6)*pow(Q2h,2) +
		           32*m4*M2*pow(Q2h,3) + 216*m4*Q2e*pow(Q2h,3) -
		           24*m2*M2*Q2e*pow(Q2h,3) - 96*m4*S*pow(Q2h,3) -
		           754*m2*Q2e*S*pow(Q2h,3) - 128*m4*U*pow(Q2h,3) -
		           546*m2*Q2e*U*pow(Q2h,3) + 84*m2*S*U*pow(Q2h,3) +
		           280*Q2e*S*U*pow(Q2h,3) + 32*pow(m,6)*pow(Q2h,3) +
		           24*m4*pow(Q2h,4) + 24*m2*M2*pow(Q2h,4) +
		           430*m2*Q2e*pow(Q2h,4) + 30*m2*S*pow(Q2h,4) -
		           376*Q2e*S*pow(Q2h,4) - 42*m2*U*pow(Q2h,4) -
		           312*Q2e*U*pow(Q2h,4) - 38*m2*pow(Q2h,5) + 204*Q2e*pow(Q2h,5) +
		           160*m4*Q2e*Q2h*pow(S,2) + 64*Q2h*pow(m,6)*pow(S,2) +
		           64*m4*pow(Q2h,2)*pow(S,2) + 328*m2*Q2e*pow(Q2h,2)*pow(S,2) +
		           4*m2*pow(Q2h,3)*pow(S,2) + 172*Q2e*pow(Q2h,3)*pow(S,2) +
		           48*m4*Q2e*Q2h*pow(U,2) + 64*Q2h*pow(m,6)*pow(U,2) +
		           80*m4*pow(Q2h,2)*pow(U,2) + 148*m2*Q2e*pow(Q2h,2)*pow(U,2) +
		           48*m2*pow(Q2h,3)*pow(U,2) + 108*Q2e*pow(Q2h,3)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(2*l1k + Q2e - Q2h,-1)*
		         (16*Q2h*S*U*pow(Q2e,3) - 8*M2*Q2h*pow(Q2e,4) -
		           12*Q2h*S*pow(Q2e,4) + 64*S*U*pow(Q2e,2)*pow(Q2h,2) +
		           8*M2*pow(Q2e,3)*pow(Q2h,2) - 56*S*pow(Q2e,3)*pow(Q2h,2) -
		           28*U*pow(Q2e,3)*pow(Q2h,2) + 8*pow(Q2e,4)*pow(Q2h,2) +
		           80*Q2e*S*U*pow(Q2h,3) - 80*S*pow(Q2e,2)*pow(Q2h,3) -
		           64*U*pow(Q2e,2)*pow(Q2h,3) + 32*pow(Q2e,3)*pow(Q2h,3) -
		           80*Q2e*S*pow(Q2h,4) - 80*Q2e*U*pow(Q2h,4) + 80*S*U*pow(Q2h,4) +
		           40*pow(Q2e,2)*pow(Q2h,4) + 40*Q2e*pow(Q2h,5) - 80*S*pow(Q2h,5) -
		           80*U*pow(Q2h,5) + 40*pow(Q2h,6) + 24*Q2h*pow(Q2e,3)*pow(S,2) +
		           4*pow(Q2e,4)*pow(S,2) + 40*pow(Q2e,2)*pow(Q2h,2)*pow(S,2) +
		           40*Q2e*pow(Q2h,3)*pow(S,2) + 40*pow(Q2h,4)*pow(S,2) +
		           20*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) + 32*Q2e*pow(Q2h,3)*pow(U,2) +
		           40*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-2)*(-32*Q2e*Q2h*S*U*pow(m,6) -
		           120*m4*Q2e*S*U*pow(Q2h,2) + 48*Q2e*S*pow(m,6)*pow(Q2h,2) +
		           16*Q2e*U*pow(m,6)*pow(Q2h,2) + 32*S*U*pow(m,6)*pow(Q2h,2) +
		           172*m4*Q2e*S*pow(Q2h,3) + 132*m4*Q2e*U*pow(Q2h,3) +
		           56*m4*S*U*pow(Q2h,3) - 204*m2*Q2e*S*U*pow(Q2h,3) -
		           16*Q2e*pow(m,6)*pow(Q2h,3) - 48*S*pow(m,6)*pow(Q2h,3) -
		           16*U*pow(m,6)*pow(Q2h,3) - 92*m4*Q2e*pow(Q2h,4) -
		           108*m4*S*pow(Q2h,4) + 252*m2*Q2e*S*pow(Q2h,4) -
		           68*m4*U*pow(Q2h,4) + 220*m2*Q2e*U*pow(Q2h,4) +
		           76*m2*S*U*pow(Q2h,4) - 64*Q2e*S*U*pow(Q2h,4) +
		           16*pow(m,6)*pow(Q2h,4) + 60*m4*pow(Q2h,5) -
		           134*m2*Q2e*pow(Q2h,5) - 124*m2*S*pow(Q2h,5) +
		           64*Q2e*S*pow(Q2h,5) - 92*m2*U*pow(Q2h,5) +
		           64*Q2e*U*pow(Q2h,5) + 70*m2*pow(Q2h,6) - 32*Q2e*pow(Q2h,6) -
		           32*Q2e*Q2h*pow(m,6)*pow(S,2) - 80*m4*Q2e*pow(Q2h,2)*pow(S,2) +
		           32*pow(m,6)*pow(Q2h,2)*pow(S,2) + 48*m4*pow(Q2h,3)*pow(S,2) -
		           118*m2*Q2e*pow(Q2h,3)*pow(S,2) + 54*m2*pow(Q2h,4)*pow(S,2) -
		           32*Q2e*pow(Q2h,4)*pow(S,2) - 40*m4*Q2e*pow(Q2h,2)*pow(U,2) +
		           8*m4*pow(Q2h,3)*pow(U,2) - 86*m2*Q2e*pow(Q2h,3)*pow(U,2) +
		           22*m2*pow(Q2h,4)*pow(U,2) - 32*Q2e*pow(Q2h,4)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-1)*pow(2*l1k + Q2e - Q2h,-1)*
		         (-40*S*U*pow(Q2h,5) + 40*S*pow(Q2h,6) + 40*U*pow(Q2h,6) -
		           20*pow(Q2h,7) - 20*pow(Q2h,5)*pow(S,2) - 20*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        pow(l1k,-3)*(32*m4*Q2e*S*U*pow(Q2h,3) - 32*m4*Q2e*S*pow(Q2h,4) -
		           32*m4*Q2e*U*pow(Q2h,4) - 32*m4*S*U*pow(Q2h,4) +
		           32*m2*Q2e*S*U*pow(Q2h,4) + 16*m4*Q2e*pow(Q2h,5) +
		           32*m4*S*pow(Q2h,5) - 32*m2*Q2e*S*pow(Q2h,5) +
		           32*m4*U*pow(Q2h,5) - 32*m2*Q2e*U*pow(Q2h,5) -
		           32*m2*S*U*pow(Q2h,5) - 16*m4*pow(Q2h,6) +
		           16*m2*Q2e*pow(Q2h,6) + 32*m2*S*pow(Q2h,6) +
		           32*m2*U*pow(Q2h,6) - 16*m2*pow(Q2h,7) +
		           16*m4*Q2e*pow(Q2h,3)*pow(S,2) - 16*m4*pow(Q2h,4)*pow(S,2) +
		           16*m2*Q2e*pow(Q2h,4)*pow(S,2) - 16*m2*pow(Q2h,5)*pow(S,2) +
		           16*m4*Q2e*pow(Q2h,3)*pow(U,2) - 16*m4*pow(Q2h,4)*pow(U,2) +
		           16*m2*Q2e*pow(Q2h,4)*pow(U,2) - 16*m2*pow(Q2h,5)*pow(U,2))*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1) +
		        (-256*m2*Q2e*U*pow(l1k,4) + 64*Q2e*Q2h*U*pow(l1k,4) -
		           128*m2*Q2e*pow(l1k,5) + 128*Q2e*Q2h*pow(l1k,5) -
		           128*Q2e*U*pow(l1k,5) - 640*Q2h*U*pow(l1k,5) -
		           128*Q2e*pow(l1k,6) - 384*Q2h*pow(l1k,6) + 256*U*pow(l1k,6) +
		           128*pow(l1k,7) - 32*Q2e*pow(l1k,4)*pow(Q2h,2) +
		           640*U*pow(l1k,4)*pow(Q2h,2) + 480*pow(l1k,5)*pow(Q2h,2) +
		           16*Q2e*U*pow(l1k,2)*pow(Q2h,3) + 8*Q2e*pow(l1k,3)*pow(Q2h,3) -
		           384*U*pow(l1k,3)*pow(Q2h,3) - 352*pow(l1k,4)*pow(Q2h,3) +
		           32*l1k*U*pow(Q2e,2)*pow(Q2h,3) + 56*l1k*Q2e*U*pow(Q2h,4) -
		           8*Q2e*pow(l1k,2)*pow(Q2h,4) + 208*U*pow(l1k,2)*pow(Q2h,4) +
		           200*pow(l1k,3)*pow(Q2h,4) - 16*l1k*pow(Q2e,2)*pow(Q2h,4) +
		           16*U*pow(Q2e,2)*pow(Q2h,4) - 30*l1k*Q2e*pow(Q2h,5) -
		           88*l1k*U*pow(Q2h,5) - 32*Q2e*U*pow(Q2h,5) -
		           112*pow(l1k,2)*pow(Q2h,5) - 8*pow(Q2e,2)*pow(Q2h,5) +
		           46*l1k*pow(Q2h,6) + 16*Q2e*pow(Q2h,6) + 16*U*pow(Q2h,6) -
		           8*pow(Q2h,7) - 128*m2*Q2e*pow(l1k,3)*pow(U,2) -
		           32*Q2e*Q2h*pow(l1k,3)*pow(U,2) - 256*Q2h*pow(l1k,4)*pow(U,2) +
		           128*pow(l1k,5)*pow(U,2) + 192*pow(l1k,3)*pow(Q2h,2)*pow(U,2) -
		           16*l1k*pow(Q2e,2)*pow(Q2h,2)*pow(U,2) -
		           24*l1k*Q2e*pow(Q2h,3)*pow(U,2) -
		           96*pow(l1k,2)*pow(Q2h,3)*pow(U,2) -
		           8*pow(Q2e,2)*pow(Q2h,3)*pow(U,2) + 40*l1k*pow(Q2h,4)*pow(U,2) +
		           16*Q2e*pow(Q2h,4)*pow(U,2) - 8*pow(Q2h,5)*pow(U,2))*
		         pow(-4*l1k*Q2h + 4*m2*Q2h + 4*pow(l1k,2) + pow(Q2h,2),-1)*
		         pow(-2*l1k*Q2e*Q2h + 2*m2*Q2e*Q2h + 4*Q2e*pow(l1k,2) +
		           2*l1k*pow(Q2e,2) - m2*pow(Q2e,2) - m2*pow(Q2h,2),-1)));

	melem_interf = c*diag_vb_off_shell/M1gamma(Q2h, Q2e, l1k, S, U, G1, G23);

	return 2.*melem_interf;
}

long double Gamma_Loop::M1gamma(const long double Q2h, const long double Q2e, const long double l1k,
		const long double S, const long double U, const long double G1, const long double G23)const {

	return -4*pow(l1k,-2)*pow(2*l1k + Q2e - Q2h,-2)*
			   (G1*(16*(Q2e - Q2h)*pow(l1k,3) +
			        16*pow(l1k,4) +
			        m2*(2*m2 - Q2h)*pow(Q2e - Q2h,2) +
			        l1k*(Q2e - Q2h)*
			         (-4*m2*Q2e + pow(Q2e,2) + pow(Q2h,2)) +
			        pow(l1k,2)*(-8*m2*Q2e - 8*Q2e*Q2h +
			           6*pow(Q2e,2) + 6*pow(Q2h,2))) +
			     G23*(4*(4*M2*(Q2e - Q2h) +
			           Q2h*(S - U))*pow(l1k,3) +
			        16*M2*pow(l1k,4) -
			        m2*(M2*Q2h + (Q2h - U)*U)*
			         pow(Q2e - Q2h,2) +
			        l1k*(Q2e - Q2h)*
			         (2*m2*(Q2h - 2*U)*(Q2h - S - U) +
			           M2*(pow(Q2e,2) + pow(Q2h,2)) +
			           Q2h*(Q2e*S + (Q2h - U)*U - pow(S,2)))\
			         + pow(l1k,2)*
			         (M2*(-8*Q2e*Q2h + 6*pow(Q2e,2) +
			              6*pow(Q2h,2)) -
			           2*Q2h*(Q2h*S - 2*Q2h*U +
			              Q2e*(-2*S + U) + pow(S,2) + pow(U,2))
						  + 4*m2*pow(-Q2h + S + U,2))));
}
