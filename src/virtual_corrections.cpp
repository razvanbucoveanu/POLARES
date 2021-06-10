/*
 * virtual_corrections.cpp
 *
 *  Created on: Dec 2, 2015
 *      Author: razvan
 */

#include "virtual_corrections.h"
#include "parameters.h"
#include "interpolation.h"
#include <cmath>
#include "const.h"
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_zeta.h>
#include <iostream>

//using namespace constants;
using namespace POLARES;

Virtual_Corrections::Virtual_Corrections(): interpolation(0), param(0), v(0), S(0),
		f1e(0), f2e(0), p1p2(0), G1(0), G2(0), G3(0), m(constants::me), m2(constants::me2),
		m4(constants::me4), M(constants::mpr), M2(constants::mpr2),
		pi(constants::pi), pi2(constants::pi2), m_me(constants::me),
		m_mu(constants::m_mu), m_tau(constants::m_tau), alpha(constants::alpha),
		sw2(constants::sw2_msbar), mw2(constants::mw2), mz2(constants::mz2) {
}

int Virtual_Corrections::set_param(const Parameters* param, const Interpolation* interpolation) {
	this->param = param;
	this->interpolation = interpolation;
	SI.set_param(param);
	GL.set_param(param);
	m = param->m;
	m2 = pow(m,2.);
	m4 = pow(m,4.);
	M = param->M;
	M2 = pow(M,2.);

	return 0;
}

double Virtual_Corrections::running_alpha(const double Q2)const {

	long double run_alpha = alpha;

	if (param->flag[param->vac_pol] == 1)
		run_alpha = alpha/(1.-d_vac_1st(Q2,m_me)/2.);
	if (param->flag[param->vac_pol] == 2)
		run_alpha = alpha/(1.-(d_vac_1st(Q2,m_me)/2.
				+ d_vac_1st(Q2, m_mu)/2. + d_vac_1st(Q2,m_tau)/2.));
	if (param->flag[param->vac_pol] >= 3)
		run_alpha = alpha/(1.-interpolation->d_vac_hadr(Q2)/2.);

	return run_alpha;
}

//vacuum polarization
double Virtual_Corrections::d_vac_1st(const double Q2, const double m_lepton)const {

	v = sqrt(1. + 4.*pow(m_lepton,2.)/Q2);

	return 2.*alpha/(3.*pi)*((pow(v,2.)-8./3.) + v*(3.-pow(v,2.))*log((v+1.)/(v-1.))/2.);
}

double Virtual_Corrections::d_vac_2nd(const double Q2, const double m_lepton)const {

	double v_log = log(Q2/pow(m_lepton,2.));

	return 2.*pow(alpha/pi,2.)*(v_log/4.+gsl_sf_zeta_int(3)-5./24.);
}

double Virtual_Corrections::f1e_1st(const double Q2)const {

	S = param->en * M;
	v = sqrt(1. + 4.*m2/Q2);

	return  alpha/pi/2. * ((pow(v,2.)+1.)*log((v+1.)/(v-1.))*log((pow(v,2.)-1.)/(4.*pow(v,2.)))/(4.*v)
			+ (2.*pow(v,2.)+1.)*log((v+1.)/(v-1.))/(2.*v)
			- 2. + (pow(v,2.)+1.)*(gsl_sf_dilog((v+1.)/(2.*v)) - gsl_sf_dilog((v-1.)/(2.*v)))/(2.*v));

}

double Virtual_Corrections::f2e_1st(const double Q2)const {

	S = param->en * M;
	v = sqrt(1. + 4.*m2/Q2);

	return	alpha/4./pi * (pow(v,2.) - 1.) * log((v+1.)/(v-1)) / v;

}

double Virtual_Corrections::f1e_2nd(const double Q2)const {

	const double L = log(Q2/m2);
	const double L2 = pow(L,2.);
	const double L3 = pow(L,3.);
	const double L4 = pow(L,4.);
	const double zeta3 = gsl_sf_zeta_int(3);

	const double total = 1171./216. - pi2*log(2.)/2. + 13./32.*pi2 - 59./1440.*pi2*pi2
			- 9./4.*zeta3 - 1627./864.*L - 13./144.*pi2*L
			+ 3./2.*zeta3*L + 229./288.*L2 - pi2*L2/48. -31./144.*L3 + L4/32.;

	const double f1e_2nd_ee = -pow(L,3.)/36. + 19.*pow(L,2.)/72. - L*(265./216. + pi2/36.);

	return pow(alpha/pi,2.) *  (total - f1e_2nd_ee);
//	return pow(alpha/pi,2.) * total;
}

//vertex correction
double Virtual_Corrections::d_vert(const double Q2, const double f1, const double f2)const {

	S = param->en * M;
	v = sqrt(1. + 4.*m2/Q2);
	f1e = f1e_1st(Q2);
	f2e = f2e_1st(Q2);

	p1p2 = M2 + Q2/2.;
	G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
	G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
	G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;

	return ((-(f2e*Q2*(6*G1 + (G2 + G3)*(4*M2 + Q2))) +
       4*f1e*(G1*(2*m2 - Q2) -
          (G2 + G3)*(M2*Q2 + 2*(Q2 - 2*S)*S)))
			*pow(G1*(2*m2 - Q2) -
       (G2 + G3)*(M2*Q2 + 2*(Q2 - 2*S)*S),-1))/2.;
}

double Virtual_Corrections::d_vert_carbon(const double Q2)const {

	S = param->en * M;
	v = sqrt(1. + 4.*m2/Q2);
	f1e = f1e_1st(Q2);
	f2e = f2e_1st(Q2);

	return 2.*f1e + f2e*Q2*(4*M2+Q2)/2./((M2*Q2+2.*S*(Q2-2.*S)));
}

double Virtual_Corrections::d_vert_quad(const double Q2, const double f1, const double f2)const {

	S = param->en * M;
	v = sqrt(1. + 4.*m2/Q2);
	f1e = f1e_1st(Q2);
	f2e = f2e_1st(Q2);

	p1p2 = M2 + Q2/2.;
	G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
	G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
	G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;

	return (pow(m,-2)*(-8*f1e*f2e*m2*Q2*
        (6*G1 + (G2 + G3)*(4*M2 + Q2)) +
       16*m2*(G1*(2*m2 - Q2) -
          (G2 + G3)*(M2*Q2 + 2*(Q2 - 2*S)*S))
         *pow(f1e,2) + Q2*pow(f2e,2)*
        (2*G1*(-8*m2 + Q2) -
          (G2 + G3)*
           (4*m2*(4*M2 + Q2) -
             pow(Q2 - 4*S,2))))*
     pow(G1*(2*m2 - Q2) -
       (G2 + G3)*(M2*Q2 + 2*(Q2 - 2*S)*S),-1))/16.;
}

double Virtual_Corrections::d_vert_pol_g(const double Q2, const double f1, const double f2,
		const double f1z, const double f2z, const double gae, const double ga, const double gv)const {

	S = param->en * M;
	v = sqrt(1. + 4.*m2/Q2);
	f1e = f1e_1st(Q2);
	f2e = f2e_1st(Q2);

	return (4*f1*M2*(Q2*(4*ga*
            (f1e*(f1z + f2z)*m2 +
              f1z*(f1e + f2e)*M2) -
           ((2*f1e*(f1z + f2z) +
                 f2e*(2*f1z + 3*f2z))*ga +
              2*(f1e + f2e)*gae*gv)*Q2) +
        8*(f1e*f1z*ga + (f1e + f2e)*gae*gv)*Q2*
         S - 16*f1e*f1z*ga*pow(S,2)) +
     f2*Q2*(f2e*Q2*(-4*
            (3*f1z*ga + 2*f2z*ga + 2*gae*gv)*M2\
            + f2z*ga*Q2) + 32*f2e*gae*gv*M2*S +
        4*f1e*(M2*(4*(f1z + f2z)*ga*m2 -
              (2*f1z*ga + f2z*ga + 2*gae*gv)*Q2)\
            + 2*(4*gae*gv*M2 + f2z*ga*Q2)*S -
           4*f2z*ga*pow(S,2))))*
   pow(4*M2*Q2*(4*ga*
         ((f1 + f2)*(f1z + f2z)*m2 + f1*f1z*M2)
          - ((2*f1*(f1z + f2z) +
              f2*(2*f1z + f2z))*ga +
           2*(f1 + f2)*gae*gv)*Q2) +
     8*Q2*(4*(f1*f1z*ga + (f1 + f2)*gae*gv)*
         M2 + f2*f2z*ga*Q2)*S -
     16*ga*(4*f1*f1z*M2 + f2*f2z*Q2)*pow(S,2),-1);;
}

double Virtual_Corrections::d_vert_pol_Z(const double Q2, const double f1, const double f2,
		const double f1z, const double f2z, const double gae, const double ga, const double gv)const {

//	double f1e, f2e, S, fa, B0q2, C0q2;

	double S = param->en * M;

	double deltaZA = -alpha/4./pi*(SI.B0_0mm(m2) + 2. + 4.*m2*SI.C0_mmqm0m(Q2, m2));
	double deltaZ1 = alpha/4./pi * (SI.B0_0mm(m2) + 4. + 4.*m2*SI.C0_mm0m0m(m2));
	double result = (pow(pi,-1)*pow(4*m2 + Q2,-1)*
		     (-4*alpha*SI.B0_m0m(m2)*
		        (2*m2*(3*M2*Q2 -
		             2*(Q2 - 2*S)*(Q2 - S)) +
		          (8*M2 - 2*Q2)*m4 +
		          Q2*(M2*Q2 - pow(Q2 - 2*S,2))) +
		       alpha*SI.B0_qmm(Q2,m2)*
		        (4*(4*M2 - Q2)*m4 -
		          3*Q2*(-(M2*Q2) +
		             pow(Q2 - 2*S,2)) +
		          m2*(16*M2*Q2 + 32*Q2*S -
		             11*pow(Q2,2) - 16*pow(S,2))) +
		       2*(4*m2 + Q2)*
		        (alpha*(m2*(4*M2 - Q2) + M2*Q2 -
		             pow(Q2 - 2*S,2)) +
		          alpha*(2*m2 + Q2)*
		           SI.C0_mmqm0m(Q2,m2)*
		           (-(M2*Q2) + m2*(-4*M2 + Q2) +
		             pow(Q2 - 2*S,2)) +
		          pi*(deltaZ1*Q2*(Q2 - 4*S) +
		             deltaZA*
		              (-2*M2*Q2 +
		                2*m2*(-4*M2 + Q2) -
		                4*Q2*S + pow(Q2,2) + 8*pow(S,2)
		                ))))*
		     pow(-(M2*Q2) + m2*(-4*M2 + Q2) +
		       pow(Q2 - 2*S,2),-1))/4.;

	return result;
}

double Virtual_Corrections::d_vert_pol(const double Q2, const double f1, const double f2,
		const double f1z, const double f2z, const double gae, const double ga, const double gv)const {

	return d_vert_pol_g(Q2, f1, f2, f1z, f2z, gae, ga, gv) + d_vert_pol_Z(Q2, f1, f2, f1z, f2z, gae, ga, gv);
}


double Virtual_Corrections::d_vert_pol_quad_gZ(const double Q2, const double f1, const double f2,
		const double f1z, const double f2z, const double gae, const double ga, const double gv)const {

	double fa, B0q2, C0q2;

	S = param->en * M;
	v = sqrt(1. + 4.*m2/Q2);
	f1e = f1e_1st(Q2);
	f2e = f2e_1st(Q2);

	B0q2 = SI.B0_qmm(Q2, m2);
	C0q2 = SI.C0_mmqm0m(Q2, m2);

	fa = (alpha*pow(pi,-1)*(-4*B0q2*m2 + 16*C0q2*m4 + 4*Q2 - 3*B0q2*Q2 + 12*C0q2*m2*Q2 -
		       (4*m2 + 3*Q2)*log(m2) + 2*C0q2*pow(Q2,2))*pow(4*m2 + Q2,-1))/4.;

	return (4*f1*M2*(8*Q2*S*(f1e*f1z*fa*ga + gae*gv*pow(f1e + f2e,2)) +
	        Q2*(-(Q2*((2*f1e*(f1z + f2z) + f2e*(2*f1z + 3*f2z))*fa*ga +
	                2*gae*gv*pow(f1e + f2e,2))) +
	           4*fa*ga*(f1e*(f1z + f2z)*m2 + f1z*(f1e + f2e)*M2)) -
	        16*f1e*f1z*fa*ga*pow(S,2)) +
	     f2*Q2*(-8*gae*gv*(Q2 - 4*S)*pow(f1e,2)*M2 +
	        f2e*(f2z*fa*ga*Q2*(Q2 - 8*M2) - 12*f1z*fa*ga*Q2*M2 -
	           8*f2e*gae*gv*(Q2 - 4*S)*M2) +
	        4*f1e*((-(((2*f1z + f2z)*fa*ga + 4*f2e*gae*gv)*Q2) +
	              4*(f1z + f2z)*fa*ga*m2)*M2 +
	           2*S*(f2z*fa*ga*Q2 + 8*f2e*gae*gv*M2) - 4*f2z*fa*ga*pow(S,2))))*
	   pow(8*Q2*S*(f2*f2z*ga*Q2 + 4*(f1*f1z*ga + (f1 + f2)*gae*gv)*M2) +
	     4*Q2*M2*(-(((2*f1*(f1z + f2z) + f2*(2*f1z + f2z))*ga +
	             2*(f1 + f2)*gae*gv)*Q2) +
	        4*ga*((f1 + f2)*(f1z + f2z)*m2 + f1*f1z*M2)) -
	     16*ga*(f2*f2z*Q2 + 4*f1*f1z*M2)*pow(S,2),-1);
}

//soft-bremsstrahlung
double Virtual_Corrections::d_brems_ee(const double Q2)const {

	double B_l1l2, B_l1l1, B_l2l2, l2, l1l2, en1;

	en1 = param->en - Q2/(2.*M);
	l2 = sqrt(pow(en1,2.)-m2);
	l1l2 = m2 + Q2/2.;

//	v = sqrt(1. + 4.*m2/Q2);
//	double d_IR = alpha/pi*log(m2/param->lambda)*(1.-(pow(v,2.)+1.)/2./v*log((v+1.)/(v-1.)));

	B_l1l1 = log(2.*param->min[param->E_gamma]/m) + param->en * log(m/(param->en + param->l1)) / param->l1;
	B_l2l2 = log(2.*param->min[param->E_gamma]/m) + en1 * log(m/(en1 + l2)) / l2;

	double beta, v, l, xi;
	beta = (l1l2 + sqrt(pow(l1l2,2.) - m4)) / m2;
	v = (beta * l1l2 - m2) / (beta*param->en - en1);
	l = beta * param->en - en1;
	xi = sqrt(1. + 4.*m2/Q2);
// (pow(xi,2.)+1.) * log((xi+1.)/(xi-1.))/(2.*xi) ~ log(Q2/m2)
	B_l1l2 = log(4.*pow(param->min[param->E_gamma],2.)/m2) * (pow(xi,2.)+1.) * log((xi+1.)/(xi-1.))/(2.*xi)
			+ beta * l1l2 * ( (pow(log((param->en-param->l1)/(param->en+param->l1)),2.)
			- pow(log((en1-l2)/(en1+l2)),2.))/4.
			+  gsl_sf_dilog(1. - beta*(param->en-param->l1)/v) -  gsl_sf_dilog(1. - (en1-l2)/v)
			+  gsl_sf_dilog(1. - beta*(param->en+param->l1)/v) -  gsl_sf_dilog(1. - (en1+l2)/v) )/ (v*l);

	return - alpha*(B_l1l1 - B_l1l2 + B_l2l2)/pi;
}

double Virtual_Corrections::d_brems_ee(const double Q2, const double en1)const {

	double B_l1l2, B_l1l1, B_l2l2, l2, l1l2;

//	en1 = param->en - Q2/(2.*M);
	l2 = sqrt(pow(en1,2.)-m2);
	l1l2 = m2 + Q2/2.;

//	v = sqrt(1. + 4.*m2/Q2);
//	double d_IR = alpha/pi*log(m2/param->lambda)*(1.-(pow(v,2.)+1.)/2./v*log((v+1.)/(v-1.)));

	B_l1l1 = log(2.*param->min[param->E_gamma]/m) + param->en * log(m/(param->en + param->l1)) / param->l1;
	B_l2l2 = log(2.*param->min[param->E_gamma]/m) + en1 * log(m/(en1 + l2)) / l2;

	double beta, v, l, xi;
	beta = (l1l2 + sqrt(pow(l1l2,2.) - m4)) / m2;
	v = (beta * l1l2 - m2) / (beta*param->en - en1);
	l = beta * param->en - en1;
	xi = sqrt(1. + 4.*m2/Q2);
// (pow(xi,2.)+1.) * log((xi+1.)/(xi-1.))/(2.*xi) ~ log(Q2/m2)
	B_l1l2 = log(4.*pow(param->min[param->E_gamma],2.)/m2) * (pow(xi,2.)+1.) * log((xi+1.)/(xi-1.))/(2.*xi)
			+ beta * l1l2 * ( (pow(log((param->en-param->l1)/(param->en+param->l1)),2.)
			- pow(log((en1-l2)/(en1+l2)),2.))/4.
			+  gsl_sf_dilog(1. - beta*(param->en-param->l1)/v) -  gsl_sf_dilog(1. - (en1-l2)/v)
			+  gsl_sf_dilog(1. - beta*(param->en+param->l1)/v) -  gsl_sf_dilog(1. - (en1+l2)/v) )/ (v*l);

	return - alpha*(B_l1l1 - B_l1l2 + B_l2l2)/pi;
}

double Virtual_Corrections::d_brems_hadr(const double Q2)const {

	double en = param->en;
	double delta_eg = param->min[param->E_gamma];
	double rho = sqrt(Q2 + 4.*M2);
	double x = pow(sqrt(Q2)+rho,2.)/4./M2;
	double en1 = en - Q2/(2.*M);
	double ep = M + en - en1;
	double p2 = sqrt(pow(ep,2.)-M2);
	double eta = en/en1;

//	double delta1 = log(4.*pow(delta_eg,2.)/Q2/x)*log(eta)
//									+ gsl_sf_dilog(1.-eta/x) - gsl_sf_dilog(1.-1./eta/x);

	double delta2 = log(4.*pow(delta_eg,2.)/M2)*(ep/p2*log(x)-1.) + 1.
									+ ep/p2*(-pow(log(x),2.)/2.-log(x)*log(pow(rho,2.)/M2)
									+ log(x) - gsl_sf_dilog(1.-1/pow(x,2.))+2.*gsl_sf_dilog(-1./x)
									+ pi2/6.);

	return alpha/pi*pow(param->Z_target,2.)*delta2;
}

double Virtual_Corrections::d_brems_hadr_interf(const double Q2)const {

	double en = param->en;
	double delta_eg = param->min[param->E_gamma];
	double rho = sqrt(Q2 + 4.*M2);
	double x = pow(sqrt(Q2)+rho,2.)/4./M2;
	double en1 = en - Q2/(2.*M);
	double ep = M + en - en1;
	double p2 = sqrt(pow(ep,2.)-M2);
	double eta = en/en1;

	double delta1 = log(4.*pow(delta_eg,2.)/Q2/x)*log(eta)
									+ gsl_sf_dilog(1.-eta/x) - gsl_sf_dilog(1.-1./eta/x);

//	double delta2 = log(4.*pow(delta_eg,2.)/M2)*(ep/p2*log(x)-1.) + 1.
//									+ ep/p2*(-pow(log(x),2.)/2.-log(x)*log(pow(rho,2.)/M2)
//									+ log(x) - gsl_sf_dilog(1.-1/pow(x,2.))+2.*gsl_sf_dilog(-1./x)
//									+ pi2/6.);

	return alpha/pi*param->Z_lepton*param->Z_target*2.*delta1;
}

double Virtual_Corrections::d_brems_ee_test(const double Q2)const {

//	double B_l1l2, B_l1l1, B_l2l2, l2, l1l2;
	double en1;

	en1 = param->en - Q2/(2.*M);
	double L = log(Q2/m2);

	double thl = acos((-Q2 * (1. + param->en/M) + 2.*pow(param->en,2.) - 2.*m2) /
				(2.*param->l1*sqrt(pow(param->en-Q2/(2.*M),2.)-m2)));

	//static limit
//	return alpha/4./pi * (8.*(L - 1.)*log(param->min[param->E_gamma]/param->en)
//			+ 2.*pow(L,2.) + 4.*gsl_sf_dilog(pow(cos(param->thl/2.),2.)) - 4.*pi2/3.);

	//Van
	return alpha/pi * (log(pow(param->min[param->E_gamma],2.)/(param->en*en1)) * (L - 1.)
			- 1./2.*pow(log(param->en/en1),2.) + 1./2.*pow(L,2.) - pow(pi,2.)/3.
			+ gsl_sf_dilog(pow(cos(thl/2.),2.)));

}

double Virtual_Corrections::d_box_Feshbach(const double Q2)const {

	double en1 = param->en - Q2/(2.*M);
	double thl = acos((-Q2 * (1. + param->en/M) + 2.*pow(param->en,2.) - 2.*m2) /
				(2.*param->l1*sqrt(pow(param->en-Q2/(2.*M),2.)-m2)));

	double delta_F = alpha*pi*(sin(thl/2.)-pow(sin(thl/2.),2.))/pow(cos(thl/2.),2.);

//	double delta_log = alpha/pi*Q2/M/param->en*log(sqrt(Q2)/param->en)
//						*(log(sqrt(Q2)/2./param->en)+1.);
//	double delta_IR = 2.*alpha/pi*log(Q2/lambda2)*log(std::abs(u-M2)/(s-M2));

	return -param->Z_lepton*param->Z_target*delta_F;
}

//double Virtual_Corrections::d_box_MTj(const double Q2)const {
//
//	double en1 = param->en - Q2/(2.*M);
//
//	double delta = log(param->en/en1)*log(pow(Q2/2./M,2.)/param->en/en1)
//						+ 2.*(gsl_sf_dilog(1.-M/2./param->en)
//								- gsl_sf_dilog(1.-M/2./en1));
//	double delta_IR = alpha/pi*Q2/M/param->en*log(param->lambda/Q2);
//
//	return alpha/pi*delta + d_box_Feshbach(Q2);
//}

long double Virtual_Corrections::d_gamma_loop(const long double Q2h, const long double Q2e, const long double l1k,
		const long double S, const long double U, const long double f1, const long double f2)const {

	long double p1p2 = M2 + Q2h/2.;
	long double G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
	long double G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
	long double G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;
	long double G23 = G2 + G3;

	long double gamma_loop;

	gamma_loop = GL.vb_on_shell(Q2h, Q2e, l1k, S, U, G1, G23) + GL.se_on_shell(Q2h, Q2e, l1k, S, U, G1, G23)
				+ GL.se_off_shell(Q2h, Q2e, l1k, S, U, G1, G23) + GL.vb_off_shell(Q2h, Q2e, l1k, S, U, G1, G23);

	return gamma_loop;
//	return SI.B0_M0m(-2*l1k + m2, m2);
}

double Virtual_Corrections::d_2nd_total(const double Q2, const double f1, const double f2)const {
//R. Hill, arXiv:1605.02613

	double vert_1st = d_vert(Q2, f1, f2);
	double vert_quad = d_vert_quad(Q2, f1, f2);
	double sp = d_brems_ee(Q2);
	double vert_2nd = 2.*f1e_2nd(Q2);

	return vert_quad + vert_1st*sp + pow(sp,2.)/2. + vert_2nd;

}

double Virtual_Corrections::d_2nd_pol_total(const double Q2, const double f1, const double f2,
		const double f1z, const double f2z, const double gae, const double ga, const double gv)const {
//R. Hill, arXiv:1605.02613

	double vert_1st = d_vert_pol(Q2, f1, f2, f1z, f2z, gae, ga, gv);
	double vert_quad = d_vert_pol_quad_gZ(Q2, f1, f2, f1z, f2z, gae, ga, gv);
	double sp = d_brems_ee(Q2);
	double vert_2nd = 2.*f1e_2nd(Q2);

	return vert_quad + vert_1st*sp + pow(sp,2.)/2. + vert_2nd;

}

double Virtual_Corrections::kappa_weak(const double Q2)const {
//arXiv : 1107.4683, formulas extracted from Jegerlehner's code alphaQED.
//via A. Weber, H. Spiesberger, 28.12.2013, Mathematica notebook

	double kappa_tot(0.), kappa_bosonic(0.), kappa_fermionic(0.);
	double cw2 = 1.-sw2;

	double Jkapbkl = -alpha/2./pi/sw2*(-(42.*cw2 + 1.)/12.*log(cw2) +
			1./18. + (6.*cw2 + 7.)/18.) - (444.*cw2 + 43.)/720.*Q2/mw2 +
					(530.*cw2 + 71.)/8400.*pow(Q2/mw2,2.);
	double Jkapbgr = -alpha/2./pi/sw2*(-(42.*cw2 + 1.)/12.*log(cw2) +
			1./18. - (1./2.*sqrt(1. + 4.*mw2/Q2)*
					log((sqrt(1. + 4.*mw2/Q2) + 1.)/(sqrt(1. + 4.*mw2/Q2) - 1.)) - 1.)
					*((7. - 4.*mw2/Q2)*cw2 + 1./6.*(1. + 4.*mw2/Q2)) -
					mw2/Q2 * (0.75 - mw2/Q2 + (mw2/Q2 - 1.5)*sqrt(1. + 4.*mw2/Q2)*
							log((sqrt(1. + 4.*mw2/Q2) + 1.)/(sqrt(1. + 4.*mw2/Q2) - 1.)) +
							mw2/Q2*(2. - mw2/Q2)*pow(log((sqrt(1. + 4.*mw2/Q2) + 1.)/
									(sqrt(1. + 4.*mw2/Q2) - 1.)),2.)));

	if (Q2 < mz2 / 1e6) kappa_bosonic = Jkapbkl;
	else kappa_bosonic = Jkapbgr;

	enum T3_enum {T3l, T3uq, T3dq, T3sq, T3cq, T3bq, T3tq};
	enum mq_enum {muq, mdq, msq, mcq, mbq, mtq};
	enum Q_enum {Ql, Quq, Qdq, Qsq, Qcq, Qbq, Qtq};
	double T3[] = {-1./2., 1./2., -1./2., -1./2., 1./2., -1./2., 1./2.};
	double mq[] = {0.1, 0.1, 0.1, 1.275, 4.18, 173.5};
	double Q[] = {-1., 2./3., -1./3., -1./3., 2./3., -1./3., 2./3.};

	if (Q2 < pow(m_me,2.) / 1e6) kappa_fermionic += Jkapflkl(T3[T3l], Q[Ql], Q2, m);
	else kappa_fermionic += Jkapflgr(T3[T3l], Q[Ql], Q2, m_me);

	if (Q2 < pow(m_mu,2.) / 1e6) kappa_fermionic += Jkapflkl(T3[T3l], Q[Ql], Q2, m_mu);
	else kappa_fermionic += Jkapflgr(T3[T3l], Q[Ql], Q2, m_mu);

	if (Q2 < pow(m_tau,2.) / 1e6) kappa_fermionic += Jkapflkl(T3[T3l], Q[Ql], Q2, m_tau);
	else kappa_fermionic += Jkapflgr(T3[T3l], Q[Ql], Q2, m_tau);

	if (Q2 < pow(mq[muq],2.) / 1e6) kappa_fermionic += Jkapfqkl(T3[T3uq], Q[Quq], Q2, mq[muq]);
	else kappa_fermionic += Jkapfqgr(T3[T3uq], Q[Quq], Q2, mq[muq]);

	if (Q2 < pow(mq[mdq],2.) / 1e6) kappa_fermionic += Jkapfqkl(T3[T3dq], Q[Qdq], Q2, mq[mdq]);
	else kappa_fermionic += Jkapfqgr(T3[T3dq], Q[Qdq], Q2, mq[mdq]);

	if (Q2 < pow(mq[msq],2.) / 1e6) kappa_fermionic += Jkapfqkl(T3[T3sq], Q[Qsq], Q2, mq[msq]);
	else kappa_fermionic += Jkapfqgr(T3[T3sq], Q[Qsq], Q2, mq[msq]);

	if (Q2 < pow(mq[mcq],2.) / 1e6) kappa_fermionic += Jkapfqkl(T3[T3cq], Q[Qcq], Q2, mq[mcq]);
	else kappa_fermionic += Jkapfqgr(T3[T3cq], Q[Qcq], Q2, mq[mcq]);

	if (Q2 < pow(mq[mbq],2.) / 1e6) kappa_fermionic += Jkapfqkl(T3[T3bq], Q[Qbq], Q2, mq[mbq]);
	else kappa_fermionic += Jkapfqgr(T3[T3bq], Q[Qbq], Q2, mq[mbq]);

	if (Q2 < pow(mq[mtq],2.) / 1e6) kappa_fermionic += Jkapfqkl(T3[T3tq], Q[Qtq], Q2, mq[mtq]);
	else kappa_fermionic += Jkapfqgr(T3[T3tq], Q[Qtq], Q2, mq[mtq]);

	kappa_fermionic *= - alpha/6./pi/sw2;

	kappa_tot = 1. + kappa_fermionic + kappa_bosonic;

	return kappa_tot;
}

double Virtual_Corrections::Jkapflkl(const double T3l, const double Ql, const double Q2, const double ml)const {

	return (T3l*Ql - 2.*sw2*pow(Ql,2.)) * (log(pow(ml,2.)/mz2) + 1./5.*Q2/pow(ml,2.) -
		    3./140.*pow(Q2,2.)/pow(ml,4.));
}

double Virtual_Corrections::Jkapflgr(const double T3l, const double Ql, const double Q2, const double ml)const {

	return (T3l*Ql - 2.*sw2*pow(Ql,2.))*(log(pow(ml,2.)/mz2) - 5./3. +
			   4.*pow(ml,2.)/Q2 + (1. - 2.*pow(ml,2.)/Q2)*sqrt(1. + 4.*pow(ml,2.)/Q2)*
			    log((sqrt(1. + 4.*pow(ml,2.)/Q2) + 1.)/(sqrt(1. + 4.*pow(ml,2.)/Q2) - 1.)));
}

double Virtual_Corrections::Jkapfqkl(const double T3q, const double Qq, const double Q2, const double mq)const {

	return 3.*(T3q*Qq - 2.*sw2*pow(Qq,2.))*(log(pow(mq,2.)/mz2) + 1./5.*Q2/pow(mq,2.) -
		    3./140.*pow(Q2,2.)/pow(mq,4.));
}

double Virtual_Corrections::Jkapfqgr(const double T3q, const double Qq, const double Q2, const double mq)const {

	return 3.*(T3q*Qq - 2.*sw2*pow(Qq,2.))*(log(pow(mq,2.)/mz2) - 5./3. +
			   4.*pow(mq,2.)/Q2 + (1. - 2.*pow(mq,2.)/Q2)*sqrt(1. + 4.*pow(mq,2.)/Q2)*
			    log((sqrt(1. + 4.*pow(mq,2.)/Q2) + 1.)/(sqrt(1. + 4.*pow(mq,2.)/Q2) - 1.)));
}
