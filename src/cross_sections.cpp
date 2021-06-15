/*
 * cross_sections.cpp
 *
 *  Created on: Dec 4., 2015
 *      Author: razvan
 */

#include "cross_sections.h"
#include <math.h>
#include "const.h"
#include "interpolation.h"
#include <iostream>
#include <iomanip>
#include "config.h"

//using namespace constants;
using namespace POLARES;

//initialization

Cross_Sections::Cross_Sections(): interpolation(0), param(0), s(0), melem_interf(0),
		melem2(0), sigma_born(0), gpe(0), gpm(0), gpze(0), gpzm(0), gae(0), f1(0),
		f2(0), f1z(0), f2z(0), tau(0), eps(0), thl(0), l2(0), a(0), b(0), l1k(0),
		cospsi(0), x(0), phig(0), gamma_loop(0), m(constants::me), m2(constants::me2),
		m4(constants::me4), m_me(constants::me), m_mu(constants::m_mu), m_tau(constants::m_tau),
		M(constants::mpr), M2(constants::mpr2), M4(constants::mpr4),
		alpha(constants::alpha), pi(constants::pi), pi2(constants::pi2), gf(constants::gf) {
	}

int Cross_Sections::set_param(const Parameters* param, const Interpolation* interpolation) {

	this->param = param;
	this->interpolation = interpolation;
	ff.change_flag(param->flag[param->form_factors]);
	VC.set_param(param, interpolation);
	melem_2nd.set_param(param);
	melem_2nd_pol.set_param(param);
	m = param->m;
	m2 = pow(m,2.);
	m4 = pow(m,4.);
	M = param->M;
	M2 = pow(M,2.);

//	fout = fopen("events_test.txt", "w");
	return 0;
}

//Born cross section
long double Cross_Sections::crsect_born (const long double Q2)const {

	ff.ffactp (Q2, gpe, gpm);

//Pauli and Dirac Form Factors
	tau = Q2/(4.*M2);
	f2 = (gpm - gpe) / (1. + tau);
	f1 = (gpe + tau*gpm) / (1. + tau);

	long double S = param->en * M;

//Matrix element squared
	//F1F2 contribution
	melem2 = 4.*f1*f2*Q2*(Q2 - 2.*m2);
	//F1^2. contribution
	melem2 -= 2.*pow(f1,2.)*
			(4.*Q2*S + 2.*Q2*m2 + 2.*Q2*M2 - pow(Q2,2.) - 8*pow(S,2.));
	//F2^2. contribution
	melem2 += Q2*pow(f2,2.)/M2*(-2.*Q2*S + Q2*M2 - 4.*m2*M2 + 4.*pow(S,2.));

	melem2 *= pow(4.*pi*alpha,2.)/pow(Q2,2.);

//Cross Section
	s = melem2 / (64.*M2*pi*(pow(param->en,2.)-m2));

	return s;
}

//Elastic Cross section
long double Cross_Sections::crsect_elastic (const long double Q2)const {

	ff.ffactp (Q2, gpe, gpm);

//Pauli and Dirac Form Factors
	tau = Q2/(4.*M2);
	f2 = (gpm - gpe) / (1. + tau);
	f1 = (gpe + tau*gpm) / (1. + tau);

	long double S = param->en * M;

//Matrix element squared
	//F1F2 contribution
	melem2 = 4.*f1*f2*Q2*(Q2 - 2.*m2)*M2;
	//F1^2. contribution
	melem2 -= 2.*pow(f1,2.)*
			M2*(4.*Q2*S + 2.*Q2*m2 + 2.*Q2*M2 - pow(Q2,2.) - 8*pow(S,2.));
	//F2^2. contribution
	melem2 += Q2*pow(f2,2.)*(-2.*Q2*S + Q2*M2 - 4.*m2*M2 + 4.*pow(S,2.));

	melem2 *= pow(4.*pi*VC.running_alpha(Q2),2.)/pow(Q2,2.)/M2;
//	melem2 *= pow(4.*pi*alpha,2.)/pow(Q2,2.)/M2;

//Cross Section
	sigma_born = melem2 / (64.*M2*pi*(pow(param->en,2.)-m2));

	s = sigma_born;

	if (param->flag[param->order] != 0) {
		s += sigma_born * (VC.d_vert(Q2, f1, f2) + VC.d_brems_ee(Q2));
//		s = sigma_born * VC.d_vert(Q2, f1, f2);

//		if (param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_1st(Q2,m_me);
//		if (param->flag[param->vac_pol] == 2.)
//			s += sigma_born * (VC.d_vac_1st(Q2,m_me) + VC.d_vac_1st(Q2,m_mu) + VC.d_vac_1st(Q2,m_tau));
//		if (param->flag[param->vac_pol] >= 3.)
//			s += sigma_born * interpolation->d_vac_hadr(Q2);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_2nd(Q2,m_me);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] > 1)
//			s += sigma_born * (VC.d_vac_2nd(Q2,m_me) + VC.d_vac_2nd(Q2,m_mu) + VC.d_vac_2nd(Q2,m_tau));

		if (param->flag[param->hadr_corr] == 2. ||
				param->flag[param->hadr_corr] == 3.) {
			s += sigma_born * VC.d_brems_hadr(Q2);
		}

//		if (param->flag[param->tpe] == 1) {
//			//		en1 = param->en - Q2/(2.*M);
//			thl = acos((-Q2 * (1. + param->en/M) + 2.*pow(param->en,2.) - 2.*m2) /
//					(2.*param->l1*sqrt(pow(param->en-Q2/(2.*M),2.)-m2)));
//			//		eps = 1. / (1. - 2.*(1. + tau) * (2.*m2 - Q2) / (4.*param->en*en1 - Q2));
//			eps = 1. / (1. + 2.*(1. + tau)*pow(tan(thl/2.),2.));
//			if (param->flag[param->lepton] == 0)
//				s += sigma_born * interpolation->tpe(Q2, eps);
//			if (param->flag[param->lepton] == 1)
//				s -= sigma_born * interpolation->tpe(Q2, eps);
//		}
	}

	if (param->flag[param->order] == 2.) {
		s += sigma_born * VC.d_2nd_total(Q2, f1, f2);
	}


	if (param->flag[param->hadr_corr] == 1 ||
			param->flag[param->hadr_corr] == 3.) {
		s += sigma_born * VC.d_brems_hadr_interf(Q2);
		if (param->flag[param->tpe] == 1) {
			s += sigma_born * VC.d_box_Feshbach(Q2);
		}

		if (param->flag[param->tpe] == 2.) {
			s += sigma_born * interpolation->tpe_Tomalak(Q2);
		}
	}

	return s;
}

//Differential born cross section in thl
long double Cross_Sections::crsect_born_thl (const long double Q2)const {

	long double mott, en1, eps;

	en1 = param->en - Q2/(2.*M);
	l2 = sqrt(pow(en1,2.) - m2);
	tau = Q2/(4.*M2);
	eps = 1./(1. - 2.*(1.+tau)*(2.*m2 - Q2) / (4.*param->en*en1 - Q2));
//	eps = 1./(1. + 2.*(1.+tau)*pow(tan(param->thl/2.),2.));
	mott = pow(alpha,2.)*(1. - Q2/(4.*param->en*en1)) * M * en1 * param->l1 * l2
	 / (4. * pow(param->en,3.) * pow(Q2/(4.*param->en*en1),2.) * (M*param->en*en1 + m2*(en1-param->en-M)));

//	mott = pow(alpha,2.) * pow(en1/param->en,2.) / (Q2*pow(tan(param->thl/2.),2.));

	ff.ffactp(Q2, gpe, gpm);

	s = (eps*pow(gpe,2.) + tau*pow(gpm,2.))*mott/(eps*(1.+tau));

//	long double thl = acos((- Q2/2. - m2 + param->en*en1) / (param->l1*l2));

	return 2.*pi*s;
}

//Differential elastic cross section in thl
long double Cross_Sections::crsect_elastic_thl (const long double Q2)const {

	long double mott, en1, eps;

	en1 = param->en - Q2/(2.*M);
	l2 = sqrt(pow(en1,2.) - m2);
	tau = Q2/(4.*M2);
	eps = 1./(1. - 2.*(1.+tau)*(2.*m2 - Q2) / (4.*param->en*en1 - Q2));
	//	eps = 1./(1. + 2.*(1.+tau)*pow(tan(param->thl/2.),2.));
	mott = pow(VC.running_alpha(Q2),2.)*(1. - Q2/(4.*param->en*en1)) * M * en1 * param->l1 * l2
			/ (4. * pow(param->en,3.) * pow(Q2/(4.*param->en*en1),2.) * (M*param->en*en1 + m2*(en1-param->en-M)));

	//	mott = pow(alpha,2.) * pow(en1/param->en,2.) / (Q2*pow(tan(param->thl/2.),2.));

	ff.ffactp(Q2, gpe, gpm);

	//Pauli and Dirac Form Factors
	tau = Q2/(4.*M2);
	f2 = (gpm - gpe) / (1. + tau);
	f1 = (gpe + tau*gpm) / (1. + tau);

	sigma_born = (eps*pow(gpe,2.) + tau*pow(gpm,2.))*mott/(eps*(1.+tau));

	s = sigma_born;

	if (param->flag[param->order] != 0) {
		s += sigma_born * (VC.d_vert(Q2, f1, f2) + VC.d_brems_ee(Q2));

//		if (param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_1st(Q2,me);
//		if (param->flag[param->vac_pol] == 2.)
//			s += sigma_born * (VC.d_vac_1st(Q2,me) + VC.d_vac_1st(Q2,m_mu) + VC.d_vac_1st(Q2,m_tau));
//		if (param->flag[param->vac_pol] >= 3.)
//			s += sigma_born * interpolation->d_vac_hadr(Q2);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_2nd(Q2,me);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] > 1)
//			s += sigma_born * (VC.d_vac_2nd(Q2,me) + VC.d_vac_2nd(Q2,m_mu) + VC.d_vac_2nd(Q2,m_tau));

//		if (param->flag[param->tpe] == 1) {
//			//		en1 = param->en - Q2/(2.*M);
//			thl = acos((-Q2 * (1. + param->en/M) + 2.*pow(param->en,2.) - 2.*m2) /
//					(2.*param->l1*sqrt(pow(param->en-Q2/(2.*M),2.)-m2)));
//			//		eps = 1. / (1. - 2.*(1. + tau) * (2.*m2 - Q2) / (4.*param->en*en1 - Q2));
//			eps = 1. / (1. + 2.*(1. + tau)*pow(tan(thl/2.),2.));
//			if (param->flag[param->lepton] == 0)
//				s += sigma_born * interpolation->tpe(Q2, eps);
//			if (param->flag[param->lepton] == 1)
//				s -= sigma_born * interpolation->tpe(Q2, eps);
//		}
	}

	if (param->flag[param->order] == 2.) {
		s +=  sigma_born * VC.d_2nd_total(Q2, f1, f2);
	}


	if (param->flag[param->hadr_corr] == 1 ||
			param->flag[param->hadr_corr] == 3.) {
		s += sigma_born * VC.d_brems_hadr_interf(Q2);

		if (param->flag[param->tpe] == 1) {
			s += sigma_born * VC.d_box_Feshbach(Q2);
		}

		if (param->flag[param->tpe] == 2.) {
			s += sigma_born * interpolation->tpe_Tomalak(Q2);
		}
	}

	if (param->flag[param->hadr_corr] == 2. ||
			param->flag[param->hadr_corr] == 3.) {
		s += sigma_born * VC.d_brems_hadr(Q2);
	}

	return 2.*pi*s;
}

// Born Interference Cross Section
long double Cross_Sections::interf_born (const long double Q2)const {

	long double en1 = param->en - Q2/(2.*M);          // energy of the scattered electron
	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = VC.kappa_weak(Q2);
		gv = -1./2. + 2.*kappa*param->sw2;

		ff.ffz (Q2, gpze, gpzm, kappa, param->sw2);
	}
	else ff.ffz (Q2, gpze, gpzm, param->sw2);
	ff.ffgae (Q2, gae);
	ff.ffactp (Q2, gpe, gpm);

// Pauli and Dirac Form Factors
	tau = Q2/(4.*M2);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);
	f2z=(gpzm-gpze)/(1.+tau);
	f1z=(gpze+tau*gpzm)/(1.+tau);

	long double sp1, sl2, sp2, S;
	l2 = sqrt(pow(en1,2.) - m2);
	thl = acos((param->en*en1 - Q2/2. - m2) / (param->l1*l2));
	S = param->en*M;
	sp1 = M*param->l1/m;
	sl2 = (param->l1*en1 - l2*param->en*cos(thl))/m;
	sp2 = (param->l1 * (param->en - en1 + M) - param->en * (param->l1 - l2*cos(thl)))/m;

// Matrix Element for the interference term
	melem_interf = 4.*f1*M2*(Q2*(-(f2z*ga*(sl2 + sp1 - sp2)) + 2.*gae*gv*(sp1 + sp2)) +
			2.*f1z*ga*(Q2*sp2 - 2.*S*(sp1 + sp2) + 2.*sl2*M2));
	melem_interf += f2*Q2*(-4.*f1z*ga*(sl2 + sp1 - sp2)*M2 + 8*gae*gv*(sp1 + sp2)*M2 +
					f2z*ga*(Q2*(sl2 + sp1 + sp2) - 4.*(S*(sp1 + sp2) + (sp1 - sp2)*M2)));
	melem_interf *= m/2./M2;

	melem_interf *= - 2.*gf*(4.*alpha*pi)/(sqrt(2.)*Q2);

//Cross Section
	s = melem_interf / (64.*M2*pi*(pow(param->en,2.)-m2));

	return s;
}

// Elastic Interference Cross Section
long double Cross_Sections::interf_elastic (const long double Q2)const {
	long double en1 = param->en - Q2/(2.*M);          // energy of the scattered electron
	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = VC.kappa_weak(Q2);
		gv = -1./2. + 2.*kappa*param->sw2;

		ff.ffz (Q2, gpze, gpzm, kappa, param->sw2);
	}
	else ff.ffz (Q2, gpze, gpzm, param->sw2);
	ff.ffgae (Q2, gae);
	ff.ffactp (Q2, gpe, gpm);

// Pauli and Dirac Form Factors
	tau = Q2/(4.*M2);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);
	f2z=(gpzm-gpze)/(1.+tau);
	f1z=(gpze+tau*gpzm)/(1.+tau);

	long double sp1, sl2, sp2, S;
	l2 = sqrt(pow(en1,2.) - m2);
	thl = acos((param->en*en1 - Q2/2. - m2) / (param->l1*l2));
	S = param->en*M;
	sp1 = M*param->l1/m;
	sl2 = (param->l1*en1 - l2*param->en*cos(thl))/m;
	sp2 = (param->l1 * (param->en - en1 + M) - param->en * (param->l1 - l2*cos(thl)))/m;

// Matrix Element for the interference term
	melem_interf = 4.*f1*M2*(Q2*(-(f2z*ga*(sl2 + sp1 - sp2)) + 2.*gae*gv*(sp1 + sp2)) +
			2.*f1z*ga*(Q2*sp2 - 2.*S*(sp1 + sp2) + 2.*sl2*M2));
	melem_interf += f2*Q2*(-4.*f1z*ga*(sl2 + sp1 - sp2)*M2 + 8*gae*gv*(sp1 + sp2)*M2 +
					f2z*ga*(Q2*(sl2 + sp1 + sp2) - 4.*(S*(sp1 + sp2) + (sp1 - sp2)*M2)));
	melem_interf *= m/2./M2;

	melem_interf *= - 2.*gf*(4.*VC.running_alpha(Q2)*pi)/(sqrt(2.)*Q2);

//Cross Section
	sigma_born = melem_interf / (64.*M2*pi*(pow(param->en,2.)-m2));

	s = sigma_born;

	if (param->flag[param->order] != 0) {

		s += sigma_born * (VC.d_vert_pol(Q2, f1, f2, f1z, f2z, gae, ga, gv) + VC.d_brems_ee(Q2));

//		if (param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_1st(Q2,me) / 2.;
//		if (param->flag[param->vac_pol] == 2.)
//			s += sigma_born * (VC.d_vac_1st(Q2,me) + VC.d_vac_1st(Q2,m_mu) + VC.d_vac_1st(Q2,m_tau)) / 2.;
//		if (param->flag[param->vac_pol] >= 3.)
//			s += sigma_born * interpolation->d_vac_hadr(Q2) / 2.;
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_2nd(Q2,me)/2.;
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] > 1)
//			s += sigma_born * (VC.d_vac_2nd(Q2,me) + VC.d_vac_2nd(Q2,m_mu) + VC.d_vac_2nd(Q2,m_tau)) / 2.;


//		if (param->flag[param->tpe] == 1) {
//			//		en1 = param->en - Q2/(2.*M);
//			tau = Q2/(4.*pow(M,2.));
//			//		eps = 1. / (1. - 2.*(1. + tau) * (2.*m2 - Q2) / (4.*param->en*en1 - Q2));
//			eps = 1. / (1. + 2.*(1. + tau)*pow(tan(thl/2.),2.));
//			if (param->flag[param->lepton] == 0)
//				s += sigma_born * interpolation->tpe(Q2, eps) / 2.;
//			if (param->flag[param->lepton] == 1)
//				s -= sigma_born * interpolation->tpe(Q2, eps) / 2.;
//		}
	}

	if (param->flag[param->order] == 2.) {
		s +=  sigma_born * VC.d_2nd_pol_total(Q2, f1, f2, f1z, f2z, gae, ga, gv);
	}

	return s;
}

long double Cross_Sections::asymm_born_test (const long double Q2)const {

	long double en1, eps, eps1, F1, F;
	long double gv = -1./2. + 2.*param->sw2;

	en1 = param->en - Q2/(2.*M);          // energy of the scattered electron

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = VC.kappa_weak(Q2);
		gv = -1./2. + 2.*kappa*param->sw2;

		ff.ffz (Q2, gpze, gpzm, kappa, param->sw2);
	}
	else ff.ffz (Q2, gpze, gpzm, param->sw2);
	ff.ffgae (Q2, gae);
	ff.ffactp (Q2, gpe, gpm);

// Pauli and Dirac Form Factors
	tau = Q2/(4.*M2);
	eps = 1. / (1. + 2.*(1.+tau)*pow(tan(param->thl/2.),2.));
	eps1 = sqrt(tau*(1.+tau)*(1.-pow(eps,2.)));
	F1 = eps*gpe*gpze + tau*gpm*gpzm - (1.-4.*param->sw2)*eps1*gpm*gae;
	F = eps*gpe*gpe + tau*gpm*gpm;

	long double melem_gamma = 64. * pow(alpha*pi*M,2.) * F / (Q2 * pow(tan(param->thl/2.),2.) * eps * (1.+tau));
	long double melem_interf = -16. * gf * alpha * pi * M2 * F1 / (sqrt(2.) * pow(tan(param->thl/2.),2.) * eps * (1.+tau));

	return melem_interf / melem_gamma;
}

//leading order interference cross section differential in thl
long double Cross_Sections::interf_born_thl (const long double Q2) const {

	long double en1 = param->en - Q2/(2.*M);          // energy of the scattered electron
	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = VC.kappa_weak(Q2);
		gv = -1./2. + 2.*kappa*param->sw2;

		ff.ffz (Q2, gpze, gpzm, kappa, param->sw2);
	}
	else ff.ffz (Q2, gpze, gpzm, param->sw2);
	ff.ffgae (Q2, gae);
	ff.ffactp (Q2, gpe, gpm);

// Pauli and Dirac Form Factors
	tau = Q2/(4.*M2);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);
	f2z=(gpzm-gpze)/(1.+tau);
	f1z=(gpze+tau*gpzm)/(1.+tau);

	long double sp1, sl2, sp2, S;
	l2 = sqrt(pow(en1,2.) - m2);
	thl = acos((param->en*en1 - Q2/2. - m2) / (param->l1*l2));
	S = param->en*M;
	sp1 = M*param->l1/m;
	sl2 = (param->l1*en1 - l2*param->en*cos(thl))/m;
	sp2 = (param->l1 * (param->en - en1 + M) - param->en * (param->l1 - l2*cos(thl)))/m;

// Matrix Element for the interference term
	melem_interf = 4.*f1*M2*(Q2*(-(f2z*ga*(sl2 + sp1 - sp2)) + 2.*gae*gv*(sp1 + sp2)) +
			2.*f1z*ga*(Q2*sp2 - 2.*S*(sp1 + sp2) + 2.*sl2*M2));
	melem_interf += f2*Q2*(-4.*f1z*ga*(sl2 + sp1 - sp2)*M2 + 8*gae*gv*(sp1 + sp2)*M2 +
					f2z*ga*(Q2*(sl2 + sp1 + sp2) - 4.*(S*(sp1 + sp2) + (sp1 - sp2)*M2)));
	melem_interf *= m/2./M2;

	melem_interf *= - 8.*gf*(4.*alpha*pi)/(4.*sqrt(2.)*Q2);

	s = melem_interf*pow(en1,2.)/(64.*M2*pow(pi,2.)*pow(param->l1,2.));

	return 2.*pi*s;
}

//leading order interference cross section differential in thl
long double Cross_Sections::interf_elastic_thl (const long double Q2) const {

	long double en1 = param->en - Q2/(2.*M);          // energy of the scattered electron
	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = VC.kappa_weak(Q2);
		gv = -1./2. + 2.*kappa*param->sw2;

		ff.ffz (Q2, gpze, gpzm, kappa, param->sw2);
	}
	else ff.ffz (Q2, gpze, gpzm, param->sw2);
	ff.ffgae (Q2, gae);
	ff.ffactp (Q2, gpe, gpm);

// Pauli and Dirac Form Factors
	tau = Q2/(4.*M2);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);
	f2z=(gpzm-gpze)/(1.+tau);
	f1z=(gpze+tau*gpzm)/(1.+tau);

	long double sp1, sl2, sp2, S;
	l2 = sqrt(pow(en1,2.) - m2);
	thl = acos((param->en*en1 - Q2/2. - m2) / (param->l1*l2));
	S = param->en*M;
	sp1 = M*param->l1/m;
	sl2 = (param->l1*en1 - l2*param->en*cos(thl))/m;
	sp2 = (param->l1 * (param->en - en1 + M) - param->en * (param->l1 - l2*cos(thl)))/m;

// Matrix Element for the interference term
	melem_interf = 4.*f1*M2*(Q2*(-(f2z*ga*(sl2 + sp1 - sp2)) + 2.*gae*gv*(sp1 + sp2)) +
			2.*f1z*ga*(Q2*sp2 - 2.*S*(sp1 + sp2) + 2.*sl2*M2));
	melem_interf += f2*Q2*(-4.*f1z*ga*(sl2 + sp1 - sp2)*M2 + 8*gae*gv*(sp1 + sp2)*M2 +
					f2z*ga*(Q2*(sl2 + sp1 + sp2) - 4.*(S*(sp1 + sp2) + (sp1 - sp2)*M2)));
	melem_interf *= m/2./M2;

	melem_interf *= - 8.*gf*(4.*VC.running_alpha(Q2)*pi)/(4.*sqrt(2.)*Q2);

	sigma_born = melem_interf*pow(en1,2.)/(64.*M2*pow(pi,2.)*pow(param->l1,2.));

	s = sigma_born;

	if (param->flag[param->order] == 1 || param->flag[param->order] == 2.) {
		s += sigma_born * (VC.d_vert_pol(Q2, f1, f2, f1z, f2z, gae, ga, gv) + VC.d_brems_ee(Q2));

//		if (param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_1st(Q2,me) / 2.;
//		if (param->flag[param->vac_pol] == 2.)
//			s += sigma_born * (VC.d_vac_1st(Q2,me) + VC.d_vac_1st(Q2,m_mu) + VC.d_vac_1st(Q2,m_tau)) / 2.;
//		if (param->flag[param->vac_pol] >= 3.)
//			s += sigma_born * interpolation->d_vac_hadr(Q2) / 2.;
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_2nd(Q2,me)/2.;
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] > 1)
//			s += sigma_born * (VC.d_vac_2nd(Q2,me) + VC.d_vac_2nd(Q2,m_mu) + VC.d_vac_2nd(Q2,m_tau)) / 2.;


//		if (param->flag[param->tpe] == 1) {
//			//		en1 = param->en - Q2/(2.*M);
//			tau = Q2/(4.*pow(M,2.));
//			//		eps = 1. / (1. - 2.*(1. + tau) * (2.*m2 - Q2) / (4.*param->en*en1 - Q2));
//			eps = 1. / (1. + 2.*(1. + tau)*pow(tan(thl/2.),2.));
//			if (param->flag[param->lepton] == 0)
//				s += sigma_born * interpolation->tpe(Q2, eps) / 2.;
//			if (param->flag[param->lepton] == 1)
//				s -= sigma_born * interpolation->tpe(Q2, eps) / 2.;
//		}
	}

	if (param->flag[param->order] == 2.) {
		s +=  sigma_born * VC.d_2nd_pol_total(Q2, f1, f2, f1z, f2z, gae, ga, gv);
	}

	return 2.*pi*s;
}

// Bremsstrahlung Cross Section
long double Cross_Sections::crsect_brems_1st(const long double en1, const long double thl, const long double eg,
		const long double thg, const long double phig) const{

	long double q2, Q2, Q2e, S1, S2;

	l2 = sqrt(pow(en1,2.)-m2);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S1 = en1 * M;
	S2 = M * (param->en - en1);

	long double S = 2.*param->en*M;
	long double U = -2.*en1*M;

	l1k = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	q2 = 2.*M*(en1 + eg - param->en);  //calculates the momentum transfer squared
	Q2 = -q2;

	tau = Q2/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	//matrix element squared in l1k, Q'2. (Q2), Q2 (Q2e), Sq (S2) and S (S1)
	//F1F2 contribution
	melem2 = 4.*f1*f2*Q2*M2*(-16.*(Q2 - Q2e)*pow(l1k,3.) + 16.*pow(l1k,4.) -
			l1k*(Q2 - Q2e)*(Q2e*(Q2e - 4.*m2) + pow(Q2,2.)) +
			m2*(-Q2 + 2.*m2)*pow(Q2 - Q2e,2.) +
			pow(l1k,2.)*(-8.*Q2*Q2e - 8.*Q2e*m2 + 6.*pow(Q2,2.) + 6.*pow(Q2e,2.)));
    //F2^2. contribution
	melem2 -= Q2*pow(f2,2.)*(-16.*pow(l1k,4.)*M2 +
			8.*pow(l1k,3.)*(Q2*(2.*S1 + S2) + 2.*(Q2 - Q2e)*M2) +
			m2*(2.*S1*(Q2 + 2.*S1) + Q2*M2 - 4.*m2*M2)*pow(Q2 - Q2e,2.) +
			2.*pow(l1k,2.)*(M2*(4.*Q2*Q2e - 3.*pow(Q2,2.) - 3.*pow(Q2e,2.)) -
			2.*Q2*(3.*Q2*S1 - 3.*Q2e*S1 + Q2*S2 - 2.*Q2e*S2 + 4.*S1*S2 + 4.*pow(S1,2.) + 2.*pow(S2,2.)) +
			2.*m2*(-4.*Q2*S2 + pow(Q2,2.) + 4.*(Q2e*M2 + pow(S2,2.)))) -
			l1k*(Q2 - Q2e)*(2.*m2*(4.*Q2*S1 - 2.*Q2*S2 - 8.*S1*S2 + 4.*Q2e*M2 +
			pow(Q2,2.)) - M2*(pow(Q2,2.) + pow(Q2e,2.)) -
			2.*Q2*(Q2*S1 - Q2e*(S1 + S2) + 2.*(2.*S1*S2 + 2.*pow(S1,2.) + pow(S2,2.)))));
	//F1^2. contribution
	melem2 -= 2.*pow(f1,2.)*M2*(16.*pow(l1k,4.)*(-Q2 + 2.*M2) -
			16.*pow(l1k,3.)*(-(Q2*(Q2 - Q2e + 2.*S1 + S2)) + 2.*(Q2 - Q2e)*M2) +
			m2*pow(Q2 - Q2e,2.)*(4.*Q2*S1 - 2.*Q2*m2 - 2.*Q2*M2 + pow(Q2,2.) +
			8.*pow(S1,2.)) + l1k*(Q2 - Q2e)*
			(-4.*m2*(Q2*(Q2e + 4.*S1 - 2.*S2) - 8.*S1*S2 + pow(Q2,2.)) -
			2.*M2*(pow(Q2,2.) + pow(Q2e,2.)) +
			Q2*(4.*Q2*S1 - 4.*Q2e*(S1 + S2) + pow(Q2,2.) + pow(Q2e,2.) +
			8.*(2.*S1*S2 + 2.*pow(S1,2.) + pow(S2,2.)))) +
			2.*pow(l1k,2.)*(M2*(-8.*Q2*Q2e + 6.*pow(Q2,2.) + 6.*pow(Q2e,2.)) +
			4.*m2*(Q2*(Q2e - 4.*S2) + pow(Q2,2.) + 4.*pow(S2,2.)) -
			Q2*(-4.*Q2*(Q2e - 3.*S1 - S2) - 4.*Q2e*(3.*S1 + 2.*S2) + 3.*pow(Q2,2.) + 3.*pow(Q2e,2.) +
			8.*(2.*S1*S2 + 2.*pow(S1,2.) + pow(S2,2.)))));

	melem2 /= pow(l1k, 2.) * M2 * pow(2.*l1k - Q2 + Q2e, 2.) * pow(q2,2.);

	melem2 *= pow(4.*pi*alpha,3.); //final matrix element squared

//check if melem2 is bigger than 0
//	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
	if (melem2 < 0.) melem2 = 0.;

//	std::cout.precision(18);
//	if (eg < 1e-1.2) std::cout << melem2 <<"\t"<< f1 <<"\t"<< f2 <<"\t"<<
//			Q2 <<"\t"<< Q2e <<"\t"<< S1 <<"\t"<< S2 <<"\t"<< l1k <<"\n";

	x = std::abs(sin(thl)*sin(thg)*sin(phig));

	s =  melem2 / (32. * pow(2.*pi,4.) * M * param->l1 * x); //final cross section formula calculated with 4.29

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

	long double sigma_brems = s;

	if (param->flag[param->order] == 2.) {
		long double v = sqrt(1. + 4.*m2/Q2e);
		long double d_IR = alpha/pi*log(m2/param->lambda)*(1.-(pow(v,2.)+1.)/2./v*log((v+1.)/(v-1.)));
		if (param->flag[param->GL] == 0) {
//		long double S = 2.*M*param->en;
//		long double U = -2.*M*en1;
			long double gamma_loop = VC.d_gamma_loop(Q2, Q2e, l1k, S, U, f1, f2);
			if (std::abs(gamma_loop) > 0.5) gamma_loop = VC.d_vert(Q2, f1, f2);
			s *=  1. + VC.d_brems_ee(Q2, en1) + gamma_loop - d_IR;
//			if (std::abs(gamma_loop) > 1e4) {
//				std::cout.precision(18);
//				long double p1p2 = M2 + Q2/2.;
//				long double G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
//				long double G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
//				long double G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;
//				long double G23 = G2 + G3;
//				std::cout << "\n invariants: "<< Q2 <<"\t"<< Q2e <<"\t"<< l1k
//						 <<"\t"<< S <<"\t"<< U <<"\t"<< f1 <<"\t"<< f2 <<"\n";
//				std::cout <<"gamma_loop: " << gamma_loop
//						<<"\t"<< VC.GL.se_on_shell(Q2, Q2e, l1k, S, U, G1, G23)
//						<<"\t"<< VC.GL.se_off_shell(Q2, Q2e, l1k, S, U, G1, G23)
//						<<"\t"<< VC.GL.vb_on_shell(Q2, Q2e, l1k, S, U, G1, G23)
//						<<"\t"<< VC.GL.vb_off_shell(Q2, Q2e, l1k, S, U, G1, G23)
//						<<"\t"<< VC.GL.M1gamma(Q2, Q2e, l1k, S, U, G1, G23) <<"\n";
//			}
//			if (std::abs(gamma_loop) > 1.) {
//				std::cout << gamma_loop << "\n";
//				std::cout << en1 <<"\t"<< eg <<"\t"<< thl <<"\t"<< thg << "\n";
//			}
		}
		if (param->flag[param->GL] == 1) {
			s *=  1. + VC.d_brems_ee(Q2e, en1) + VC.d_vert(Q2, f1, f2);
		}

//		if (param->flag[param->vac_pol] == 1)
//			s += sigma_brems * VC.d_vac_1st(Q2,me);
//		if (param->flag[param->vac_pol] == 2.)
//			s += sigma_brems * (VC.d_vac_1st(Q2,me) + VC.d_vac_1st(Q2,m_mu) + VC.d_vac_1st(Q2,m_tau));
//		if (param->flag[param->vac_pol] >= 3.)
//			s += sigma_brems * interpolation->d_vac_hadr(Q2);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_2nd(Q2,me);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] > 1)
//			s += sigma_born * (VC.d_vac_2nd(Q2,me) + VC.d_vac_2nd(Q2,m_mu) + VC.d_vac_2nd(Q2,m_tau));

//		if (param->flag[param->tpe] == 1) {
//		//		en1 = param->en - Q2/(2.*M);
//			tau = Q2/(4.*pow(M,2.));
//		//		eps = 1. / (1. - 2.*(1. + tau) * (2.*m2 - Q2) / (4.*param->en*en1 - Q2));
//			eps = 1. / (1. + 2.*(1. + tau)*pow(tan(thl/2.),2.));
//			if (param->flag[param->lepton] == 0)
//			s += sigma_brems * interpolation->tpe(Q2, eps);
//			if (param->flag[param->lepton] == 1)
//			s -= sigma_brems * interpolation->tpe(Q2, eps);
//		}
	}

	//the cross section is multiplied by two to account for the range pi to 2.*pi of the azymuthal angle of the photon
	return 2.*s;
}

// Bremsstrahlung Cross Section
long double Cross_Sections::crsect_brems_1st_test(const long double en1, const long double thl, const long double eg,
		const long double thg, const long double phig) const{

	long double q2, Q2, Q2e, S1, S2;

	l2 = sqrt(pow(en1,2.)-m2);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S1 = en1 * M;
	S2 = M * (param->en - en1);

	long double S = 2.*param->en*M;
	long double U = -2.*en1*M;

	l1k = param->en*eg - param->l1*eg*sin(thg)*sin(thl)*cos(phig) - param->l1*eg*cos(thg)*cos(thl); //the product between the 4.-momenta of the initial electron and final photon
	q2 = 2.*M*(en1 + eg - param->en);  //calculates the momentum transfer squared
	Q2 = -q2;

	tau = Q2/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	//matrix element squared in l1k, Q'2. (Q2), Q2 (Q2e), Sq (S2) and S (S1)
	//F1F2 contribution
	melem2 = 4.*f1*f2*Q2*M2*(-16.*(Q2 - Q2e)*pow(l1k,3.) + 16.*pow(l1k,4.) -
			l1k*(Q2 - Q2e)*(Q2e*(Q2e - 4.*m2) + pow(Q2,2.)) +
			m2*(-Q2 + 2.*m2)*pow(Q2 - Q2e,2.) +
			pow(l1k,2.)*(-8.*Q2*Q2e - 8.*Q2e*m2 + 6.*pow(Q2,2.) + 6.*pow(Q2e,2.)));
    //F2^2. contribution
	melem2 -= Q2*pow(f2,2.)*(-16.*pow(l1k,4.)*M2 +
			8.*pow(l1k,3.)*(Q2*(2.*S1 + S2) + 2.*(Q2 - Q2e)*M2) +
			m2*(2.*S1*(Q2 + 2.*S1) + Q2*M2 - 4.*m2*M2)*pow(Q2 - Q2e,2.) +
			2.*pow(l1k,2.)*(M2*(4.*Q2*Q2e - 3.*pow(Q2,2.) - 3.*pow(Q2e,2.)) -
			2.*Q2*(3.*Q2*S1 - 3.*Q2e*S1 + Q2*S2 - 2.*Q2e*S2 + 4.*S1*S2 + 4.*pow(S1,2.) + 2.*pow(S2,2.)) +
			2.*m2*(-4.*Q2*S2 + pow(Q2,2.) + 4.*(Q2e*M2 + pow(S2,2.)))) -
			l1k*(Q2 - Q2e)*(2.*m2*(4.*Q2*S1 - 2.*Q2*S2 - 8.*S1*S2 + 4.*Q2e*M2 +
			pow(Q2,2.)) - M2*(pow(Q2,2.) + pow(Q2e,2.)) -
			2.*Q2*(Q2*S1 - Q2e*(S1 + S2) + 2.*(2.*S1*S2 + 2.*pow(S1,2.) + pow(S2,2.)))));
	//F1^2. contribution
	melem2 -= 2.*pow(f1,2.)*M2*(16.*pow(l1k,4.)*(-Q2 + 2.*M2) -
			16.*pow(l1k,3.)*(-(Q2*(Q2 - Q2e + 2.*S1 + S2)) + 2.*(Q2 - Q2e)*M2) +
			m2*pow(Q2 - Q2e,2.)*(4.*Q2*S1 - 2.*Q2*m2 - 2.*Q2*M2 + pow(Q2,2.) +
			8.*pow(S1,2.)) + l1k*(Q2 - Q2e)*
			(-4.*m2*(Q2*(Q2e + 4.*S1 - 2.*S2) - 8.*S1*S2 + pow(Q2,2.)) -
			2.*M2*(pow(Q2,2.) + pow(Q2e,2.)) +
			Q2*(4.*Q2*S1 - 4.*Q2e*(S1 + S2) + pow(Q2,2.) + pow(Q2e,2.) +
			8.*(2.*S1*S2 + 2.*pow(S1,2.) + pow(S2,2.)))) +
			2.*pow(l1k,2.)*(M2*(-8.*Q2*Q2e + 6.*pow(Q2,2.) + 6.*pow(Q2e,2.)) +
			4.*m2*(Q2*(Q2e - 4.*S2) + pow(Q2,2.) + 4.*pow(S2,2.)) -
			Q2*(-4.*Q2*(Q2e - 3.*S1 - S2) - 4.*Q2e*(3.*S1 + 2.*S2) + 3.*pow(Q2,2.) + 3.*pow(Q2e,2.) +
			8.*(2.*S1*S2 + 2.*pow(S1,2.) + pow(S2,2.)))));

	melem2 /= pow(l1k, 2.) * M2 * pow(2.*l1k - Q2 + Q2e, 2.) * pow(q2,2.);

	melem2 *= pow(4.*pi*alpha,3.); //final matrix element squared

//check if melem2 is bigger than 0
//	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
	if (melem2 < 0.) melem2 = 0.;

//	std::cout.precision(18);
//	if (eg < 1e-1.2) std::cout << melem2 <<"\t"<< f1 <<"\t"<< f2 <<"\t"<<
//			Q2 <<"\t"<< Q2e <<"\t"<< S1 <<"\t"<< S2 <<"\t"<< l1k <<"\n";

	x = std::abs(param->l1*sin(thl)*sin(thg)*sin(phig));

	s =  melem2 * l2 / (32. * pow(2.*pi,4.) * M * param->l1 * x); //final cross section formula calculated with 4.29

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

	long double sigma_brems = s;

	if (param->flag[param->order] == 2.) {
		long double v = sqrt(1. + 4.*m2/Q2e);
		long double d_IR = alpha/pi*log(m2/param->lambda)*(1.-(pow(v,2.)+1.)/2./v*log((v+1.)/(v-1.)));
		if (param->flag[param->GL] == 0) {
//		long double S = 2.*M*param->en;
//		long double U = -2.*M*en1;
			long double gamma_loop = VC.d_gamma_loop(Q2, Q2e, l1k, S, U, f1, f2);
			if (std::abs(gamma_loop) > 0.5) gamma_loop = VC.d_vert(Q2, f1, f2);
			s *=  1. + VC.d_brems_ee(Q2, en1) + gamma_loop - d_IR;
//			if (std::abs(gamma_loop) > 1e4) {
//				std::cout.precision(18);
//				long double p1p2 = M2 + Q2/2.;
//				long double G1 = 4.*pow(f1+f2,2.)*(M2-p1p2);
//				long double G2 = pow(f2,2.)*p1p2/M2 - f2*(4.*f1+3.*f2);
//				long double G3 = pow(2.*f1+f2,2.) + pow(f2,2.)*p1p2/M2;
//				long double G23 = G2 + G3;
//				std::cout << "\n invariants: "<< Q2 <<"\t"<< Q2e <<"\t"<< l1k
//						 <<"\t"<< S <<"\t"<< U <<"\t"<< f1 <<"\t"<< f2 <<"\n";
//				std::cout <<"gamma_loop: " << gamma_loop
//						<<"\t"<< VC.GL.se_on_shell(Q2, Q2e, l1k, S, U, G1, G23)
//						<<"\t"<< VC.GL.se_off_shell(Q2, Q2e, l1k, S, U, G1, G23)
//						<<"\t"<< VC.GL.vb_on_shell(Q2, Q2e, l1k, S, U, G1, G23)
//						<<"\t"<< VC.GL.vb_off_shell(Q2, Q2e, l1k, S, U, G1, G23)
//						<<"\t"<< VC.GL.M1gamma(Q2, Q2e, l1k, S, U, G1, G23) <<"\n";
//			}
//			if (std::abs(gamma_loop) > 1.) {
//				std::cout << gamma_loop << "\n";
//				std::cout << en1 <<"\t"<< eg <<"\t"<< thl <<"\t"<< thg << "\n";
//			}
		}
		if (param->flag[param->GL] == 1) {
			s *=  1. + VC.d_brems_ee(Q2e, en1) + VC.d_vert(Q2, f1, f2);
		}

//		if (param->flag[param->vac_pol] == 1)
//			s += sigma_brems * VC.d_vac_1st(Q2,me);
//		if (param->flag[param->vac_pol] == 2.)
//			s += sigma_brems * (VC.d_vac_1st(Q2,me) + VC.d_vac_1st(Q2,m_mu) + VC.d_vac_1st(Q2,m_tau));
//		if (param->flag[param->vac_pol] >= 3.)
//			s += sigma_brems * interpolation->d_vac_hadr(Q2);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_2nd(Q2,me);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] > 1)
//			s += sigma_born * (VC.d_vac_2nd(Q2,me) + VC.d_vac_2nd(Q2,m_mu) + VC.d_vac_2nd(Q2,m_tau));

//		if (param->flag[param->tpe] == 1) {
//		//		en1 = param->en - Q2/(2.*M);
//			tau = Q2/(4.*pow(M,2.));
//		//		eps = 1. / (1. - 2.*(1. + tau) * (2.*m2 - Q2) / (4.*param->en*en1 - Q2));
//			eps = 1. / (1. + 2.*(1. + tau)*pow(tan(thl/2.),2.));
//			if (param->flag[param->lepton] == 0)
//			s += sigma_brems * interpolation->tpe(Q2, eps);
//			if (param->flag[param->lepton] == 1)
//			s -= sigma_brems * interpolation->tpe(Q2, eps);
//		}
	}

	//the cross section is multiplied by two to account for the range pi to 2.*pi of the azymuthal angle of the photon
	return 2.*s;
}

long double Cross_Sections::crsect_brems_1st_hadr(const long double en1, const long double thl, const long double eg,
		const long double thg, const long double phig) const{

	long double q2, Q2, Q2e, S1, S2;

	l2 = sqrt(pow(en1,2.)-m2);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S1 = en1 * M;
	S2 = M * (param->en - en1);

	l1k = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	q2 = 2.*M*(en1 + eg - param->en);  //calculates the momentum transfer squared
	Q2 = -q2;

	tau = Q2e/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2e, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	long double S = 2.*param->en*M;
	long double U = -2.*en1*M;

//	long double F12 = f1;
//	long double F22 = f2;
//	long double kfli = l1k;
//	long double kflf = l1k+Q2e/2.-Q2/2.;
//	long double kfpi = S/2.+U/2.-Q2/2.;
//	long double kfpf = S/2.+U/2.-Q2e/2.;
//	long double lfpi = -U/2.;
//	long double lfpf = -U/2.+Q2/2.-l1k;
//	long double lilf = m2+Q2e/2.;
//	long double lipi = S/2.;
//	long double lipf = S/2.-Q2e/2.-l1k;
//	long double pipf = M2 + Q2/2.;

//	long double F10 = 1.;
//	long double F20 = 1.79;
	long double f10 = 1.;
	long double f20 = 1.79;

	melem2 = (pow(-Q2 + S + U,-2.)*pow(-Q2e + S + U,-2.)/M4*
		     (2.*f1*f2*M2*(4.*pow(f10,2.)*M2*
		           (2.*m2*(4.*Q2e*M2*pow(-Q2 + Q2e,2.) +
		                (Q2e - S - U)*(-Q2 + S + U)*pow(Q2 + Q2e,2.)) +
		             Q2e*(-4.*Q2e*M2*pow(-Q2 + Q2e,2.) -
		                (Q2e - S - U)*(-Q2 + S + U)*
		                 (4.*l1k*(-Q2 + Q2e) + 8*pow(l1k,2.) + pow(Q2 + Q2e,2.)))) +
		          f10*f20*(Q2e - S - U)*(-Q2 + S + U)*
		           (6*l1k*Q2*Q2e*S - 14*l1k*Q2*Q2e*U - 8*l1k*Q2*S*U + 8*l1k*Q2e*S*U -
		             2.*Q2*Q2e*S*U + 8*Q2e*S*pow(l1k,2.) + 8*Q2e*U*pow(l1k,2.) +
		             16*S*U*pow(l1k,2.) + 4.*(3.*Q2 + Q2e)*(Q2 - S - U)*(Q2e - S - U)*
		              m2 - Q2e*S*pow(Q2,2.) + 4.*Q2e*U*pow(Q2,2.) + 6*l1k*S*pow(Q2e,2.) +
		             3.*Q2*S*pow(Q2e,2.) + 2.*l1k*U*pow(Q2e,2.) - Q2*U*pow(Q2e,2.) -
		             2.*S*U*pow(Q2e,2.) - pow(Q2,2.)*pow(Q2e,2.) - Q2*pow(Q2e,3.) +
		             2.*S*pow(Q2e,3.) + U*pow(Q2e,3.) - 2.*Q2e*M2*pow(-Q2 + Q2e,2.) -
		             8*l1k*Q2*pow(S,2.) - 4.*l1k*Q2e*pow(S,2.) - 2.*Q2*Q2e*pow(S,2.) +
		             8*pow(l1k,2.)*pow(S,2.) + 2.*pow(Q2,2.)*pow(S,2.) - 2.*pow(Q2e,2.)*pow(S,2.) +
		             12*l1k*Q2e*pow(U,2.) - 4.*Q2*Q2e*pow(U,2.) + 8*pow(l1k,2.)*pow(U,2.) +
		             2.*pow(Q2e,2.)*pow(U,2.)) +
		          (Q2e - S - U)*(-Q2 + S + U)*pow(f20,2.)*
		           (-2.*Q2*Q2e*S*U + 8*pow(l1k,2.)*
		              ((S + U)*(Q2e + S + U) + 4.*Q2e*M2) - Q2e*S*pow(Q2,2.) +
		             4.*Q2e*U*pow(Q2,2.) + 2.*Q2e*M2*pow(Q2,2.) + 3.*Q2*S*pow(Q2e,2.) -
		             Q2*U*pow(Q2e,2.) - 2.*S*U*pow(Q2e,2.) - 4.*Q2*M2*pow(Q2e,2.) -
		             pow(Q2,2.)*pow(Q2e,2.) - Q2*pow(Q2e,3.) + 2.*S*pow(Q2e,3.) + U*pow(Q2e,3.) +
		             2.*M2*pow(Q2e,3.) -
		             4.*m2*(-((3.*Q2 + Q2e)*(Q2 - S - U)*(Q2e - S - U)) +
		                2.*M2*pow(-Q2 + Q2e,2.)) - 2.*Q2*Q2e*pow(S,2.) +
		             2.*pow(Q2,2.)*pow(S,2.) - 2.*pow(Q2e,2.)*pow(S,2.) - 4.*Q2*Q2e*pow(U,2.) +
		             2.*pow(Q2e,2.)*pow(U,2.) +
		             2.*l1k*(-4.*Q2*S*(S + U) + 8*Q2e*(-Q2 + Q2e)*M2 +
		                (3.*S + U)*pow(Q2e,2.) +
		                Q2e*(3.*Q2*S - 7*Q2*U + 4.*S*U - 2.*pow(S,2.) + 6*pow(U,2.))))) +
		       4.*pow(f1,2.)*M2*(-2.*f10*f20*(Q2e - S - U)*(-Q2 + S + U)*M2*
		           (4.*l1k*Q2e*(-Q2 + Q2e) + 8*Q2e*pow(l1k,2.) -
		             (-Q2e + 2.*m2)*pow(-Q2 + Q2e,2.)) +
		          (Q2e - S - U)*(-Q2 + S + U)*pow(f20,2.)*
		           (4.*l1k*(-(Q2*S*(S + U)) + Q2e*(Q2*(S - U) + U*(S + U)) +
		                Q2e*(-Q2 + Q2e)*M2) - Q2e*S*pow(Q2,2.) + Q2e*U*pow(Q2,2.) +
		             Q2*S*pow(Q2e,2.) - Q2*U*pow(Q2e,2.) -
		             2.*m2*(2.*Q2*(Q2e - S - U)*(-Q2 + S + U) +
		                M2*pow(-Q2 + Q2e,2.)) - Q2*Q2e*pow(S,2.) +
		             pow(Q2,2.)*pow(S,2.) - Q2*Q2e*pow(U,2.) + pow(Q2e,2.)*pow(U,2.) +
		             4.*pow(l1k,2.)*(2.*Q2e*M2 + pow(S + U,2.))) -
		          2.*pow(f10,2.)*M2*(8*Q2*Q2e*S*U*M2 +
		             4.*l1k*(Q2e - S - U)*(Q2e*(Q2 - S - U)*(Q2 - Q2e + S - U) +
		                2.*(-Q2 + Q2e)*(Q2e - 2.*S)*M2) +
		             8*(Q2e - S - U)*pow(l1k,2.)*
		              (Q2e*(-Q2 + S + U) + 2.*(Q2e - S - U)*M2) -
		             4.*Q2e*S*U*pow(Q2,2.) - 8*Q2e*S*M2*pow(Q2,2.) -
		             4.*Q2e*U*M2*pow(Q2,2.) - 4.*Q2e*M4*pow(Q2,2.) +
		             Q2e*S*pow(Q2,3.) + Q2e*U*pow(Q2,3.) - 4.*Q2*S*U*pow(Q2e,2.) +
		             4.*Q2*S*M2*pow(Q2e,2.) - 4.*Q2*U*M2*pow(Q2e,2.) +
		             8*Q2*M4*pow(Q2e,2.) + S*pow(Q2,2.)*pow(Q2e,2.) +
		             3.*U*pow(Q2,2.)*pow(Q2e,2.) + 6*M2*pow(Q2,2.)*pow(Q2e,2.) -
		             pow(Q2,3.)*pow(Q2e,2.) + 3.*Q2*S*pow(Q2e,3.) + Q2*U*pow(Q2e,3.) -
		             4.*S*U*pow(Q2e,3.) - 4.*Q2*M2*pow(Q2e,3.) -
		             4.*S*M2*pow(Q2e,3.) - 4.*M4*pow(Q2e,3.) - Q2*pow(Q2e,4.) +
		             S*pow(Q2e,4.) + U*pow(Q2e,4.) + 2.*M2*pow(Q2e,4.) +
		             4.*Q2*Q2e*U*pow(S,2.) - 4.*Q2*Q2e*M2*pow(S,2.) -
		             Q2e*pow(Q2,2.)*pow(S,2.) + 4.*M2*pow(Q2,2.)*pow(S,2.) -
		             4.*Q2*pow(Q2e,2.)*pow(S,2.) + 6*U*pow(Q2e,2.)*pow(S,2.) +
		             4.*M2*pow(Q2e,2.)*pow(S,2.) - 3.*pow(Q2e,3.)*pow(S,2.) +
		             2.*Q2*Q2e*pow(S,3.) - 4.*Q2e*U*pow(S,3.) + 4.*pow(Q2e,2.)*pow(S,3.) -
		             2.*Q2e*pow(S,4.) + 6*Q2*Q2e*S*pow(U,2.) + 4.*Q2*Q2e*M2*pow(U,2.) -
		             3.*Q2e*pow(Q2,2.)*pow(U,2.) - 4.*Q2*pow(Q2e,2.)*pow(U,2.) +
		             4.*S*pow(Q2e,2.)*pow(U,2.) - pow(Q2e,3.)*pow(U,2.) -
		             4.*Q2e*pow(S,2.)*pow(U,2.) + 4.*Q2*Q2e*pow(U,3.) - 4.*Q2e*S*pow(U,3.) +
		             2.*pow(Q2e,2.)*pow(U,3.) - 2.*Q2e*pow(U,4.) -
		             2.*m2*(2.*Q2e*M2*pow(-Q2 + Q2e,2.) +
		                (Q2e - S - U)*(-Q2 + S + U)*
		                 (2.*Q2e*(Q2 - S - U) - 2.*Q2*(S + U) + pow(Q2,2.) + pow(Q2e,2.) +
		                   2.*pow(S + U,2.))))) -
		       pow(f2,2.)*(-2.*f10*f20*(Q2e - S - U)*(-Q2 + S + U)*M2*
		           (-2.*Q2*Q2e*S*U + 8*pow(l1k,2.)*
		              ((S + U)*(Q2e + S + U) + 4.*Q2e*M2) - Q2e*S*pow(Q2,2.) +
		             4.*Q2e*U*pow(Q2,2.) + 2.*Q2e*M2*pow(Q2,2.) + 3.*Q2*S*pow(Q2e,2.) -
		             Q2*U*pow(Q2e,2.) - 2.*S*U*pow(Q2e,2.) - 4.*Q2*M2*pow(Q2e,2.) -
		             pow(Q2,2.)*pow(Q2e,2.) - Q2*pow(Q2e,3.) + 2.*S*pow(Q2e,3.) + U*pow(Q2e,3.) +
		             2.*M2*pow(Q2e,3.) -
		             4.*m2*(-((3.*Q2 + Q2e)*(Q2 - S - U)*(Q2e - S - U)) +
		                2.*M2*pow(-Q2 + Q2e,2.)) - 2.*Q2*Q2e*pow(S,2.) +
		             2.*pow(Q2,2.)*pow(S,2.) - 2.*pow(Q2e,2.)*pow(S,2.) - 4.*Q2*Q2e*pow(U,2.) +
		             2.*pow(Q2e,2.)*pow(U,2.) +
		             2.*l1k*(-4.*Q2*S*(S + U) + 8*Q2e*(-Q2 + Q2e)*M2 +
		                (3.*S + U)*pow(Q2e,2.) +
		                Q2e*(3.*Q2*S - 7*Q2*U + 4.*S*U - 2.*pow(S,2.) + 6*pow(U,2.)))) +
		          (Q2 - S - U)*(Q2e - S - U)*pow(f20,2.)*
		           (-(Q2e*S*U*pow(Q2,2.)) - 2.*Q2e*S*M2*pow(Q2,2.) +
		             2.*Q2e*U*M2*pow(Q2,2.) + 2.*Q2e*M4*pow(Q2,2.) -
		             2.*Q2*S*U*pow(Q2e,2.) + 4.*Q2*S*M2*pow(Q2e,2.) -
		             4.*S*U*M2*pow(Q2e,2.) - 4.*Q2*M4*pow(Q2e,2.) +
		             U*pow(Q2,2.)*pow(Q2e,2.) + Q2*S*pow(Q2e,3.) - S*U*pow(Q2e,3.) -
		             2.*Q2*M2*pow(Q2e,3.) + 2.*S*M2*pow(Q2e,3.) +
		             2.*U*M2*pow(Q2e,3.) + 2.*M4*pow(Q2e,3.) +
		             2.*l1k*(-(Q2e*(Q2e - S - U)*(S - U)*(-Q2 + S + U)) +
		                8*Q2e*(-Q2 + Q2e)*M4 +
		                2.*M2*(-2.*Q2*S*(S + U) + 2.*Q2e*(Q2*(S - U) + U*(S + U)) -
		                   Q2*pow(Q2e,2.) + pow(Q2e,3.))) + 2.*Q2*Q2e*U*pow(S,2.) -
		             2.*Q2*Q2e*M2*pow(S,2.) + 2.*M2*pow(Q2,2.)*pow(S,2.) -
		             2.*Q2*pow(Q2e,2.)*pow(S,2.) + 3.*U*pow(Q2e,2.)*pow(S,2.) -
		             2.*M2*pow(Q2e,2.)*pow(S,2.) - pow(Q2e,3.)*pow(S,2.) +
		             Q2*Q2e*pow(S,3.) - 2.*Q2e*U*pow(S,3.) + 2.*pow(Q2e,2.)*pow(S,3.) -
		             Q2e*pow(S,4.) + 3.*Q2*Q2e*S*pow(U,2.) - 2.*Q2*Q2e*M2*pow(U,2.) -
		             Q2e*pow(Q2,2.)*pow(U,2.) - 2.*Q2*pow(Q2e,2.)*pow(U,2.) +
		             2.*S*pow(Q2e,2.)*pow(U,2.) - 2.*Q2e*pow(S,2.)*pow(U,2.) + 2.*Q2*Q2e*pow(U,3.) -
		             2.*Q2e*S*pow(U,3.) + pow(Q2e,2.)*pow(U,3.) - Q2e*pow(U,4.) +
		             8*pow(l1k,2.)*M2*
		              (4.*Q2e*M2 + pow(Q2e,2.) + pow(S + U,2.)) -
		             2.*m2*(-4.*(Q2 + Q2e)*(Q2 - S - U)*(Q2e - S - U)*M2 +
		                4.*M4*pow(-Q2 + Q2e,2.) -
		                pow(-Q2 + S + U,2.)*pow(-Q2e + S + U,2.))) +
		          2.*pow(f10,2.)*M2*(Q2e*
		              ((Q2e - S - U)*(-Q2 + S + U)*
		                 (Q2*(2.*U*(S + U) - Q2*(S + 2.*U)) +
		                   Q2e*(Q2*(3.*S + U) - 2.*U*(3.*S + U) + pow(Q2,2.)) -
		                   3.*(Q2 - U)*pow(Q2e,2.)) + 4.*Q2e*M4*pow(-Q2 + Q2e,2.) -
		                4.*M2*(S*pow(Q2e,3.) +
		                   pow(Q2e,2.)*(-(Q2*(3.*S + U)) + pow(Q2,2.) - pow(S,2.)) -
		                   pow(Q2,2.)*pow(S,2.) + Q2*Q2e*(2.*S*U + U*(-Q2 + U) + 3.*pow(S,2.))))\
		              + 2.*l1k*Q2e*(Q2e - S - U)*
		              (4.*(-Q2 + Q2e)*(Q2e - 2.*S)*M2 -
		                (Q2 - S - U)*(Q2*S + 3.*Q2*U + Q2e*(-2.*Q2 - 5*S + U) + 2.*pow(Q2e,2.) +
		                   2.*pow(S,2.) - 2.*pow(U,2.))) +
		             8*Q2e*pow(l1k,2.)*(-Q2 + S + U + 2.*M2)*pow(-Q2e + S + U,2.) -
		             4.*m2*(4.*Q2*Q2e*(Q2e - S - U)*(-Q2 + S + U)*M2 +
		                4.*Q2e*M4*pow(-Q2 + Q2e,2.) -
		                Q2*pow(-Q2 + S + U,2.)*pow(-Q2e + S + U,2.))))))/2.;

	melem2 *= pow(4.*pi*alpha,3.)/pow(Q2e,2.); //final matrix element squared

//check if melem2 is bigger than 0
//	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
	if (melem2 < 0.) melem2 = 0.;

	x = sin(thl) * sin(thg) * sin(phig);

	s =  melem2 / (32. * pow(2.*pi,4.) * M * param->l1 * x); //final cross section formula calculated with 4.29

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

	//the cross section is multiplied by two to account for the range pi to 2.*pi of the azymuhtal angle of the photon
	return 2.*s;
}

long double Cross_Sections::crsect_brems_1st_hadr_interf(const long double en1, const long double thl, const long double eg,
		const long double thg, const long double phig) const{

	long double q2, Q2, Q2e, S1, S2;

	l2 = sqrt(pow(en1,2.)-m2);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S1 = en1 * M;
	S2 = M * (param->en - en1);

	l1k = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	q2 = 2.*M*(en1 + eg - param->en);  //calculates the momentum transfer squared
	Q2 = -q2;

	tau = Q2/(4.*M2);
	ff.ffactp (Q2, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	long double tau_e = Q2e/(4.*M2);
	ff.ffactp (Q2e, gpe, gpm);
	long double f2e=(gpm-gpe)/(1.+tau);
	long double f1e=(gpe+tau*gpm)/(1.+tau);

	long double S = 2.*param->en*M;
	long double U = -2.*en1*M;

	long double f10 = 1.;
	long double f20 = 1.79;

	melem2 = (pow(l1k,-1.)/M2*pow(2.*l1k - Q2 + Q2e,-1.)*pow(Q2 - S - U,-1.)*
		     pow(Q2e - S - U,-1.)*(f1*(f20*(-2.*f1e*
		              ((-Q2 + Q2e)*(-(Q2*S) + Q2e*(Q2 - U)) -
		                2.*l1k*(Q2*(S + U) + Q2e*(-2.*Q2 + S + U)))*M2*
		              (4.*l1k*(-Q2 + Q2e) + 8*pow(l1k,2.) + pow(-Q2 + Q2e,2.)) +
		             f2e*(16*pow(l1k,3.)*(Q2*(Q2e - S - U)*(S + U) +
		                   2.*Q2e*(-Q2 + Q2e)*M2) +
		                4.*pow(l1k,2.)*(2.*Q2e*(-Q2 + Q2e)*(-3.*Q2 + 3.*Q2e - 2.*S + 2.*U)*
		                    M2 + 2.*Q2*(S + U)*
		                    (2.*Q2*S + Q2*U + pow(S,2.) - pow(U,2.)) -
		                   Q2*Q2e*(4.*Q2*S + 2.*Q2*U + 6*S*U + 7*pow(S,2.) - pow(U,2.)) +
		                   pow(Q2e,2.)*(6*Q2*S - pow(S,2.) + pow(U,2.))) -
		                (-Q2 + Q2e)*(4.*(Q2e - 2.*S)*m2*
		                    (Q2*(Q2e - S - U)*(-Q2 + S + U) + M2*pow(-Q2 + Q2e,2.)) +
		                   (-(Q2*S) + Q2e*(Q2 - U))*
		                    (-(S*pow(Q2e,2.)) + M2*pow(-Q2 + Q2e,2.) -
		                      Q2*(-2.*S*U + Q2*(S + U) + pow(S,2.) - pow(U,2.)) +
		                      Q2e*(Q2*(S - 2.*U) + pow(Q2,2.) + pow(S,2.) + pow(U,2.)))) -
		                2.*l1k*(6*Q2e*S*U*pow(Q2,2.) - 5*Q2e*S*pow(Q2,3.) - 3.*Q2e*U*pow(Q2,3.) +
		                   4.*S*U*pow(Q2,3.) + 4.*S*pow(Q2,2.)*pow(Q2e,2.) -
		                   4.*U*pow(Q2,2.)*pow(Q2e,2.) + 2.*pow(Q2,3.)*pow(Q2e,2.) -
		                   5*Q2*S*pow(Q2e,3.) + Q2*U*pow(Q2e,3.) + 2.*S*U*pow(Q2e,3.) -
		                   M2*(Q2*(S + U) + Q2e*(-4.*Q2 - 3.*S + 5*U) + 2.*pow(Q2e,2.))*
		                    pow(-Q2 + Q2e,2.) +
		                   8*(Q2e - S - U)*m2*
		                    (Q2*(Q2e - S - U)*(-Q2 + S + U) + M2*pow(-Q2 + Q2e,2.)) -
		                   4.*Q2*Q2e*U*pow(S,2.) - 7*Q2e*pow(Q2,2.)*pow(S,2.) -
		                   U*pow(Q2,2.)*pow(S,2.) + 3.*pow(Q2,3.)*pow(S,2.) +
		                   9*Q2*pow(Q2e,2.)*pow(S,2.) - U*pow(Q2e,2.)*pow(S,2.) +
		                   pow(Q2e,3.)*pow(S,2.) - 4.*Q2*Q2e*pow(S,3.) + 3.*pow(Q2,2.)*pow(S,3.) -
		                   pow(Q2e,2.)*pow(S,3.) + 3.*Q2e*pow(Q2,2.)*pow(U,2.) -
		                   5*S*pow(Q2,2.)*pow(U,2.) + pow(Q2,3.)*pow(U,2.) +
		                   3.*Q2*pow(Q2e,2.)*pow(U,2.) - S*pow(Q2e,2.)*pow(U,2.) -
		                   pow(Q2e,3.)*pow(U,2.) - pow(Q2,2.)*pow(U,3.) - pow(Q2e,2.)*pow(U,3.))))\
		           + 2.*f10*M2*(f2e*
		              (-8*Q2e*(-Q2 + Q2e)*(3.*Q2 - 2.*(S + 2.*U))*pow(l1k,2.) -
		                32*Q2e*(Q2 - S - U)*pow(l1k,3.) -
		                (-Q2 + Q2e)*(2.*Q2e*m2*
		                    (-2.*Q2e*(2.*Q2 + S - U) + Q2*(-Q2 + 6*S + 2.*U) + pow(Q2e,2.)) +
		                   (-(Q2*S) + Q2e*(Q2 - U))*(Q2*Q2e + pow(Q2,2.) + 2.*pow(Q2e,2.))) +
		                2.*l1k*(2.*Q2e*(-2.*Q2 + S + 3.*U)*pow(Q2,2.) + (S + U)*pow(Q2,3.) +
		                   Q2*(2.*Q2 + 3.*S - 5*U)*pow(Q2e,2.) -
		                   4.*Q2e*m2*(-4.*Q2*Q2e + Q2*(-Q2 + 4.*(S + U)) + pow(Q2e,2.)) +
		                   (-6*Q2 + 2.*S + 6*U)*pow(Q2e,3.))) +
		             2.*f1e*(16*pow(l1k,3.)*(Q2*(S + U) + Q2e*(-2.*Q2 + S + U) +
		                   4.*(Q2e - S - U)*M2) +
		                8*pow(l1k,2.)*(2.*(-Q2 + Q2e)*(3.*Q2e - 2.*(2.*S + U))*M2 +
		                   (-3.*Q2 + S + 2.*U)*pow(Q2e,2.) -
		                   Q2*(2.*Q2*S + Q2*U + pow(S,2.) - pow(U,2.)) +
		                   Q2e*(3.*Q2*S - 3.*Q2*U + 3.*pow(Q2,2.) - pow(S,2.) + pow(U,2.))) +
		                (-Q2 + Q2e)*(-((-(Q2*S) + Q2e*(Q2 - U))*
		                      (-2.*Q2e*S - 2.*Q2*U - 2.*(Q2 + Q2e)*M2 + pow(Q2,2.) +
		                        pow(Q2e,2.) + 2.*pow(S,2.) + 2.*pow(U,2.))) +
		                   2.*m2*(Q2e*(3.*Q2 - 4.*U)*(Q2 - S - U) -
		                      (3.*S + 2.*U)*pow(Q2,2.) + (Q2 - U)*pow(Q2e,2.) +
		                      2.*Q2*(4.*S*U + pow(S,2.) + 3.*pow(U,2.)) - 4.*U*pow(S + U,2.))) +
		                2.*l1k*(-5*Q2e*S*pow(Q2,2.) + 9*Q2e*U*pow(Q2,2.) - 4.*S*U*pow(Q2,2.) -
		                   4.*Q2e*pow(Q2,3.) + 3.*S*pow(Q2,3.) + U*pow(Q2,3.) +
		                   9*Q2*S*pow(Q2e,2.) - 5*Q2*U*pow(Q2e,2.) - 4.*S*U*pow(Q2e,2.) +
		                   4.*pow(Q2,2.)*pow(Q2e,2.) - 4.*Q2*pow(Q2e,3.) + S*pow(Q2e,3.) +
		                   3.*U*pow(Q2e,3.) +
		                   2.*M2*(2.*Q2*Q2e*(2.*Q2 + 3.*S - U) - (5*S + U)*pow(Q2,2.) -
		                      (2.*Q2 + 5*S + U)*pow(Q2e,2.) + 2.*pow(Q2e,3.)) -
		                   8*Q2*Q2e*pow(S,2.) + 2.*Q2*U*pow(S,2.) + 2.*Q2e*U*pow(S,2.) +
		                   2.*pow(Q2,2.)*pow(S,2.) - 2.*pow(Q2e,2.)*pow(S,2.) + 2.*Q2*pow(S,3.) +
		                   2.*Q2e*pow(S,3.) - 8*Q2*Q2e*pow(U,2.) + 2.*Q2*S*pow(U,2.) +
		                   2.*Q2e*S*pow(U,2.) - 2.*pow(Q2,2.)*pow(U,2.) + 2.*pow(Q2e,2.)*pow(U,2.) +
		                   2.*Q2*pow(U,3.) + 2.*Q2e*pow(U,3.) +
		                   2.*m2*((2.*Q2 - S - U)*pow(Q2e,2.) +
		                      2.*Q2e*(-5*Q2*(S + U) + 3.*pow(Q2,2.) + 2.*pow(S + U,2.)) -
		                      (S + U)*(-8*Q2*(S + U) + 5*pow(Q2,2.) + 4.*pow(S + U,2.))))))) +
		       f2*(f20*(f2e*((-Q2 + Q2e)*(-3.*Q2e*S*U*pow(Q2,2.) + 2.*Q2e*S*pow(Q2,3.) +
		                   2.*Q2e*U*pow(Q2,3.) - 2.*S*U*pow(Q2,3.) + U*pow(Q2,2.)*pow(Q2e,2.) -
		                   pow(Q2,3.)*pow(Q2e,2.) + Q2*S*pow(Q2e,3.) - S*U*pow(Q2e,3.) -
		                   4.*(Q2 + Q2e - 2.*(S + U))*m2*
		                    (Q2*(Q2e - S - U)*(-Q2 + S + U) + M2*pow(-Q2 + Q2e,2.)) +
		                   2.*U*pow(Q2,2.)*pow(S,2.) - pow(Q2,3.)*pow(S,2.) -
		                   2.*Q2*pow(Q2e,2.)*pow(S,2.) + U*pow(Q2e,2.)*pow(S,2.) +
		                   Q2*Q2e*pow(S,3.) + Q2*Q2e*S*pow(U,2.) - 4.*Q2e*pow(Q2,2.)*pow(U,2.) +
		                   2.*S*pow(Q2,2.)*pow(U,2.) + Q2*pow(Q2e,2.)*pow(U,2.) +
		                   2.*Q2*Q2e*pow(U,3.) - pow(Q2e,2.)*pow(U,3.)) -
		                2.*l1k*(6*Q2e*S*U*pow(Q2,2.) - 2.*Q2e*S*pow(Q2,3.) - 8*Q2e*U*pow(Q2,3.) +
		                   4.*S*U*pow(Q2,3.) + S*pow(Q2,2.)*pow(Q2e,2.) +
		                   U*pow(Q2,2.)*pow(Q2e,2.) + 2.*pow(Q2,3.)*pow(Q2e,2.) -
		                   5*Q2*S*pow(Q2e,3.) + Q2*U*pow(Q2e,3.) + 2.*S*U*pow(Q2e,3.) -
		                   2.*Q2e*(Q2 + Q2e - 2.*(S + U))*M2*pow(-Q2 + Q2e,2.) +
		                   8*(Q2 + Q2e - 2.*(S + U))*m2*
		                    (Q2*(Q2e - S - U)*(-Q2 + S + U) + M2*pow(-Q2 + Q2e,2.)) -
		                   7*Q2*Q2e*U*pow(S,2.) - 5*Q2e*pow(Q2,2.)*pow(S,2.) +
		                   2.*U*pow(Q2,2.)*pow(S,2.) + 10*Q2*pow(Q2e,2.)*pow(S,2.) -
		                   U*pow(Q2e,2.)*pow(S,2.) + pow(Q2e,3.)*pow(S,2.) - 5*Q2*Q2e*pow(S,3.) +
		                   4.*pow(Q2,2.)*pow(S,3.) - pow(Q2e,2.)*pow(S,3.) -
		                   7*Q2*Q2e*S*pow(U,2.) + 13*Q2e*pow(Q2,2.)*pow(U,2.) -
		                   4.*S*pow(Q2,2.)*pow(U,2.) + 2.*pow(Q2,3.)*pow(U,2.) -
		                   8*Q2*pow(Q2e,2.)*pow(U,2.) + 5*S*pow(Q2e,2.)*pow(U,2.) -
		                   pow(Q2e,3.)*pow(U,2.) - 5*Q2*Q2e*pow(U,3.) - 2.*pow(Q2,2.)*pow(U,3.) +
		                   5*pow(Q2e,2.)*pow(U,3.)) +
		                16*pow(l1k,3.)*(2.*Q2e*(Q2 + Q2e - 2.*(S + U))*M2 -
		                   (S + U)*(-(Q2*Q2e) + pow(S + U,2.))) +
		                4.*pow(l1k,2.)*(Q2*(S + U)*(Q2*(-S + U) + 6*S*(S + U)) +
		                   6*Q2e*(-Q2 + Q2e)*(Q2 + Q2e - 2.*(S + U))*M2 +
		                   pow(Q2e,2.)*(6*Q2*S - pow(S,2.) + pow(U,2.)) -
		                   2.*Q2e*(3.*U*pow(Q2,2.) + 5*Q2*(pow(S,2.) - pow(U,2.)) +
		                      3.*U*pow(S + U,2.)))) -
		             f1e*(16*(S + U)*pow(l1k,3.)*
		                 (-((Q2 - S - U)*(S + U)) + 2.*(-Q2 + Q2e)*M2) +
		                (-Q2 + Q2e)*(4.*(Q2 - 2.*U)*m2*
		                    (Q2*(Q2e - S - U)*(-Q2 + S + U) + M2*pow(-Q2 + Q2e,2.)) +
		                   (-(Q2*S) + Q2e*(Q2 - U))*
		                    (M2*pow(-Q2 + Q2e,2.) +
		                      Q2*(-(Q2*U) + pow(S,2.) + pow(U,2.)) -
		                      Q2e*(Q2*S - 2.*Q2*U + 2.*pow(U,2.)))) +
		                2.*l1k*(8*(Q2 - S - U)*m2*
		                    (Q2*(Q2e - S - U)*(-Q2 + S + U) + M2*pow(-Q2 + Q2e,2.)) +
		                   (-Q2 + Q2e)*((-Q2 + Q2e)*(-(Q2*(5*S + U)) + Q2e*(2.*Q2 - S + 3.*U))*
		                       M2 +
		                      Q2e*((-3.*S + 5*U)*pow(Q2,2.) + Q2*(pow(S,2.) - 11*pow(U,2.)) +
		                         6*(S + U)*pow(U,2.)) +
		                      Q2*(3.*Q2*pow(S,2.) - 3.*U*pow(S,2.) - pow(S,3.) - Q2*pow(U,2.) -
		                         S*pow(U,2.) + pow(U,3.)))) +
		                4.*pow(l1k,2.)*(4.*(-Q2 + Q2e)*(-(Q2*(2.*S + U)) + Q2e*(S + 2.*U))*
		                    M2 + Q2*(S + U)*
		                    (Q2*(5*S + U) - 2.*(3.*S*U + 2.*pow(S,2.) + pow(U,2.))) +
		                   Q2e*(-4.*(S - U)*pow(Q2,2.) +
		                      3.*Q2*(-2.*S*U + pow(S,2.) - 3.*pow(U,2.)) + 6*U*pow(S + U,2.))))) +
		          f10*(2.*f1e*M2*(-8*(-Q2 + Q2e)*
		                 (Q2e*(6*Q2 - S - 2.*U) - 3.*Q2*(2.*S + U))*pow(l1k,2.) +
		                16*(3.*Q2*(S + U) + Q2e*(-4.*Q2 + S + U))*pow(l1k,3.) +
		                (-Q2 + Q2e)*(2.*Q2*m2*
		                    (Q2e*(4.*Q2 - 2.*S - 6*U) - Q2*(Q2 + 2.*S - 2.*U) + pow(Q2e,2.)) -
		                   (-(Q2*S) + Q2e*(Q2 - U))*(Q2*Q2e + 2.*pow(Q2,2.) + pow(Q2e,2.))) +
		                2.*l1k*(Q2e*(-8*Q2 - 9*S + 5*U)*pow(Q2,2.) + 2.*(4.*S + U)*pow(Q2,3.) +
		                   2.*Q2*(3.*Q2 + 4.*S - U)*pow(Q2e,2.) +
		                   4.*Q2*m2*(4.*Q2e*(Q2 - S - U) - pow(Q2,2.) + pow(Q2e,2.)) +
		                   (-6*Q2 + S + 3.*U)*pow(Q2e,3.))) +
		             f2e*(16*pow(l1k,3.)*((Q2e - S - U)*(S + U)*(-Q2 + S + U) +
		                   2.*Q2e*(-Q2 + Q2e)*M2) +
		                (-Q2 + Q2e)*(-4.*m2*
		                    (Q2*(Q2 - 2.*U)*(Q2 - S - U)*(-Q2e + S + U) +
		                      M2*(-3.*Q2*Q2e*(Q2 - 2.*(S + U)) +
		                         (Q2 - 2.*U)*pow(Q2,2.) - (3.*Q2 + 2.*S)*pow(Q2e,2.) + pow(Q2e,3.))
		                      ) - (-(Q2*S) + Q2e*(Q2 - U))*
		                    (-2.*Q2e*S*U + (-Q2 + U)*pow(Q2e,2.) +
		                      M2*pow(Q2 + Q2e,2.) + Q2*(-(Q2*U) + pow(S,2.) + pow(U,2.))
		                      )) + 4.*pow(l1k,2.)*
		                 (2.*Q2e*(-Q2 + Q2e)*(-3.*Q2 + 3.*Q2e - 2.*S + 2.*U)*M2 +
		                   pow(Q2e,2.)*(6*S*U + U*(-6*Q2 + 5*U) + pow(S,2.)) +
		                   Q2*(S + U)*(-(Q2*(5*S + U)) +
		                      2.*(3.*S*U + 2.*pow(S,2.) + pow(U,2.))) +
		                   2.*Q2e*(3.*S*pow(Q2,2.) - 2.*Q2*(pow(S,2.) - pow(U,2.)) -
		                      (S + 2.*U)*pow(S + U,2.))) -
		                2.*l1k*(2.*Q2e*S*U*pow(Q2,2.) + 3.*Q2e*S*pow(Q2,3.) - 3.*Q2e*U*pow(Q2,3.) -
		                   4.*Q2*S*U*pow(Q2e,2.) - S*pow(Q2,2.)*pow(Q2e,2.) -
		                   U*pow(Q2,2.)*pow(Q2e,2.) + 6*Q2*U*pow(Q2e,3.) - 2.*S*U*pow(Q2e,3.) -
		                   2.*pow(Q2,2.)*pow(Q2e,3.) -
		                   M2*(-(Q2e*(4.*Q2 + S - 7*U)*pow(Q2,2.)) +
		                      (S + U)*pow(Q2,3.) + Q2*(2.*Q2 + 11*S - 5*U)*pow(Q2e,2.) +
		                      (-8*Q2 - 3.*S + 5*U)*pow(Q2e,3.) + 2.*pow(Q2e,4.)) -
		                   5*Q2*Q2e*U*pow(S,2.) + 2.*Q2e*pow(Q2,2.)*pow(S,2.) +
		                   3.*U*pow(Q2,2.)*pow(S,2.) - 3.*pow(Q2,3.)*pow(S,2.) +
		                   3.*Q2*pow(Q2e,2.)*pow(S,2.) + 4.*U*pow(Q2e,2.)*pow(S,2.) -
		                   3.*Q2*Q2e*pow(S,3.) + pow(Q2,2.)*pow(S,3.) - 5*Q2*Q2e*S*pow(U,2.) +
		                   6*Q2e*pow(Q2,2.)*pow(U,2.) + S*pow(Q2,2.)*pow(U,2.) +
		                   pow(Q2,3.)*pow(U,2.) - Q2*pow(Q2e,2.)*pow(U,2.) +
		                   6*S*pow(Q2e,2.)*pow(U,2.) - 4.*pow(Q2e,3.)*pow(U,2.) -
		                   3.*Q2*Q2e*pow(U,3.) - pow(Q2,2.)*pow(U,3.) + 2.*pow(Q2e,2.)*pow(U,3.) +
		                   8*m2*(M2*
		                       (-3.*Q2*Q2e*(Q2 - 2.*(S + U)) + (Q2 - S - U)*pow(Q2,2.) -
		                         (3.*Q2 + S + U)*pow(Q2e,2.) + pow(Q2e,3.)) -
		                      Q2*(Q2e - S - U)*pow(-Q2 + S + U,2.))))))))/4.;

	melem2 *= pow(4.*pi*alpha,3.)/(Q2e*Q2); //final matrix element squared
	melem2 *= -param->Z_lepton*param->Z_target; //to take into account the charge numbers

	x = sin(thl) * sin(thg) * sin(phig);

	s =  melem2 / (32. * pow(2.*pi,4.) * M * param->l1 * x); //final cross section formula calculated with 4.29

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

	//the cross section is multiplied by two to account for the range pi to 2.*pi of the azymuhtal angle of the photon
	return 2.*s;
}

long double Cross_Sections::crsect_brems_1st_ps2(const long double en1, const long double thl, const long double eg,
		const long double thg, const long double phig) const{

	long double q2, Q2, Q2e, S1, S2;

	l2 = sqrt(pow(en1,2.)-m2);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S1 = en1 * M;
	S2 = M * (param->en - en1);

//	long double S = 2.*param->en*M;
//	long double U = -2.*en1*M;

	l1k = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
//	long double l2k = l2k = en1*eg - eg*l2*cospsi;

	q2 = 2.*M*(en1 + eg - param->en);  //calculates the momentum transfer squared
	Q2 = -q2;

	tau = Q2/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	//matrix element squared in l1k, Q'2. (Q2), Q2 (Q2e), Sq (S2) and S (S1)
	//F1F2 contribution
	melem2 = 4.*f1*f2*Q2*M2*(-16.*(Q2 - Q2e)*pow(l1k,3.) + 16.*pow(l1k,4.) -
			l1k*(Q2 - Q2e)*(Q2e*(Q2e - 4.*m2) + pow(Q2,2.)) +
			m2*(-Q2 + 2.*m2)*pow(Q2 - Q2e,2.) +
			pow(l1k,2.)*(-8.*Q2*Q2e - 8.*Q2e*m2 + 6.*pow(Q2,2.) + 6.*pow(Q2e,2.)));
    //F2^2. contribution
	melem2 -= Q2*pow(f2,2.)*(-16.*pow(l1k,4.)*M2 +
			8.*pow(l1k,3.)*(Q2*(2.*S1 + S2) + 2.*(Q2 - Q2e)*M2) +
			m2*(2.*S1*(Q2 + 2.*S1) + Q2*M2 - 4.*m2*M2)*pow(Q2 - Q2e,2.) +
			2.*pow(l1k,2.)*(M2*(4.*Q2*Q2e - 3.*pow(Q2,2.) - 3.*pow(Q2e,2.)) -
			2.*Q2*(3.*Q2*S1 - 3.*Q2e*S1 + Q2*S2 - 2.*Q2e*S2 + 4.*S1*S2 + 4.*pow(S1,2.) + 2.*pow(S2,2.)) +
			2.*m2*(-4.*Q2*S2 + pow(Q2,2.) + 4.*(Q2e*M2 + pow(S2,2.)))) -
			l1k*(Q2 - Q2e)*(2.*m2*(4.*Q2*S1 - 2.*Q2*S2 - 8.*S1*S2 + 4.*Q2e*M2 +
			pow(Q2,2.)) - M2*(pow(Q2,2.) + pow(Q2e,2.)) -
			2.*Q2*(Q2*S1 - Q2e*(S1 + S2) + 2.*(2.*S1*S2 + 2.*pow(S1,2.) + pow(S2,2.)))));
	//F1^2. contribution
	melem2 -= 2.*pow(f1,2.)*M2*(16.*pow(l1k,4.)*(-Q2 + 2.*M2) -
			16.*pow(l1k,3.)*(-(Q2*(Q2 - Q2e + 2.*S1 + S2)) + 2.*(Q2 - Q2e)*M2) +
			m2*pow(Q2 - Q2e,2.)*(4.*Q2*S1 - 2.*Q2*m2 - 2.*Q2*M2 + pow(Q2,2.) +
			8.*pow(S1,2.)) + l1k*(Q2 - Q2e)*
			(-4.*m2*(Q2*(Q2e + 4.*S1 - 2.*S2) - 8.*S1*S2 + pow(Q2,2.)) -
			2.*M2*(pow(Q2,2.) + pow(Q2e,2.)) +
			Q2*(4.*Q2*S1 - 4.*Q2e*(S1 + S2) + pow(Q2,2.) + pow(Q2e,2.) +
			8.*(2.*S1*S2 + 2.*pow(S1,2.) + pow(S2,2.)))) +
			2.*pow(l1k,2.)*(M2*(-8.*Q2*Q2e + 6.*pow(Q2,2.) + 6.*pow(Q2e,2.)) +
			4.*m2*(Q2*(Q2e - 4.*S2) + pow(Q2,2.) + 4.*pow(S2,2.)) -
			Q2*(-4.*Q2*(Q2e - 3.*S1 - S2) - 4.*Q2e*(3.*S1 + 2.*S2) + 3.*pow(Q2,2.) + 3.*pow(Q2e,2.) +
			8.*(2.*S1*S2 + 2.*pow(S1,2.) + pow(S2,2.)))));

	melem2 /= pow(l1k, 2.) * M2 * pow(2.*l1k - Q2 + Q2e, 2.) * pow(q2,2.);

	melem2 *= pow(4.*pi*alpha,3.); //final matrix element squared

//check if melem2 is bigger than 0
//	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
	if (melem2 < 0.) melem2 = 0.;

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double A = param->l1*cos(thl)-eg*cospsi;
	long double B = param->en+M-eg;

	x = std::abs(A*en1-B*l2)/eg/pow(l2,2.);

	s =  melem2 / (32. * pow(2.*pi,4.) * M * param->l1 * x); //final cross section formula calculated with 4.29

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

	long double sigma_brems = s;

	if (param->flag[param->order] == 2.) {
		long double v = sqrt(1. + 4.*m2/Q2e);
		long double d_IR = alpha/pi*log(m2/param->lambda)*(1.-(pow(v,2.)+1.)/2./v*log((v+1.)/(v-1.)));
		if (param->flag[param->GL] == 0) {
//			long double S = 2.*M*param->en;
//			long double U = -2.*M*en1;
			gamma_loop = VC.d_gamma_loop(Q2, Q2e, l1k, S1, S2, f1, f2);
			s *=  1. + VC.d_brems_ee(Q2, en1) + gamma_loop - d_IR;
			if (std::abs(gamma_loop) > 0.5) {
//				std::cout << gamma_loop << "\n";
				s = 0.;
				}
		}
		if (param->flag[param->GL] == 1) {
			s *=  1. + VC.d_brems_ee(Q2e, en1) + VC.d_vert(Q2, f1, f2);
		}

//		if (param->flag[param->vac_pol] == 1)
//			s += sigma_brems * VC.d_vac_1st(Q2,me);
//		if (param->flag[param->vac_pol] == 2.)
//			s += sigma_brems * (VC.d_vac_1st(Q2,me) + VC.d_vac_1st(Q2,m_mu) + VC.d_vac_1st(Q2,m_tau));
//		if (param->flag[param->vac_pol] >= 3.)
//			s += sigma_brems * interpolation->d_vac_hadr(Q2);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_2nd(Q2,me);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] > 1)
//			s += sigma_born * (VC.d_vac_2nd(Q2,me) + VC.d_vac_2nd(Q2,m_mu) + VC.d_vac_2nd(Q2,m_tau));


//		if (param->flag[param->tpe] == 1) {
//		//		en1 = param->en - Q2/(2.*M);
//			tau = Q2/(4.*pow(M,2.));
//		//		eps = 1. / (1. - 2.*(1. + tau) * (2.*m2 - Q2) / (4.*param->en*en1 - Q2));
//			eps = 1. / (1. + 2.*(1. + tau)*pow(tan(thl/2.),2.));
//			if (param->flag[param->lepton] == 0)
//			s += sigma_brems * interpolation->tpe(Q2, eps);
//			if (param->flag[param->lepton] == 1)
//			s -= sigma_brems * interpolation->tpe(Q2, eps);
//		}
	}

	return s;
}

long double Cross_Sections::crsect_brems_1st_sg_diff(const long double en1, const long double thl, const long double eg,
		const long double thg, const long double phig) const{

	long double q2, Q2, Q2e, S1, S2;

	l2 = sqrt(pow(en1,2.)-m2);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S1 = en1 * M;
	S2 = M * (param->en - en1);

	l1k = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
//	long double l2k = l2k = en1*eg - eg*l2*cospsi;

	q2 = 2.*M*(en1 + eg - param->en);  //calculates the momentum transfer squared
	Q2 = -q2;
	long double Q2h = Q2;

	tau = Q2/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	//matrix element squared in l1k, Q'2. (Q2), Q2 (Q2e), Sq (S2) and S (S1)
	melem2 = -(pow(l1k,-1.)*pow(M,-2.)*pow(2.*l1k + Q2e - Q2h,-1.)*
		      (8*f1*f2*M2*(4.*(Q2e - 2.*Q2h)*pow(l1k,2.) -
		           (Q2h + m2)*pow(Q2e - Q2h,2.) +
		           2.*l1k*(-3.*Q2e*Q2h + pow(Q2e,2.) + 2.*pow(Q2h,2.))) +
		        pow(f2,2.)*((Q2e - Q2h)*Q2h*(2.*Q2h*S1 + Q2e*(Q2h + 2.*S1) -
		              2.*Q2h*S2 - 8*S1*S2 - 2.*(Q2e - Q2h)*M2) +
		           4.*pow(l1k,2.)*(-(Q2h*(Q2e - 4.*S2)) +
		              4.*(Q2e - Q2h)*M2) + 2.*l1k*(Q2h*
		               (Q2h*(4.*S1 - 2.*S2) - 8*S2*(2.*S1 + S2) +
		                 Q2e*(Q2h + 4.*S1 + 6*S2) - pow(Q2e,2.)) +
		              4.*M2*pow(Q2e - Q2h,2.)) +
		           m2*(-4.*M2*pow(Q2e - Q2h,2.) +
		              Q2h*pow(Q2e + Q2h - 4.*S2,2.))) +
		        4.*pow(f1,2.)*M2*(-2.*Q2e*Q2h*S2 - 8*Q2e*S1*S2 +
		           8*Q2h*S1*S2 + 4.*Q2e*Q2h*m2 -
		           8*Q2e*S2*m2 - 8*Q2h*S2*m2 +
		           8*pow(l1k,2.)*(-Q2h + 2.*S2 + 2.*M2) +
		           2.*S1*pow(Q2e,2.) + 2.*M2*pow(Q2e - Q2h,2.) +
		           Q2e*pow(Q2h,2.) - 2.*S1*pow(Q2h,2.) + 2.*S2*pow(Q2h,2.) -
		           pow(Q2h,3.) + 4.*l1k*(-(Q2e*Q2h) + 2.*Q2e*S1 + 2.*Q2h*S1 +
		        3.*Q2e*S2 - Q2h*S2 - 8*S1*S2 + 2.*(Q2e - Q2h)*M2 +
		        pow(Q2h,2.) - 4.*pow(S2,2.)) + 16*m2*pow(S2,2.))))/2.;

	melem2 *= pow(4.*pi*alpha,3.)/pow(q2,2.); //final matrix element squared

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double A = param->l1*cos(thl)-eg*cospsi;
	long double B = param->en+M-eg;

	x = std::abs(A*en1-B*l2)/eg/pow(l2,2.);

	s =  melem2 / (32. * pow(2.*pi,4.) * M * param->l1 * x); //final cross section formula calculated with 4.29

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

	return s;
}

long double Cross_Sections::crsect_brems_2nd (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared
	Sk = Sq2 - Q2k/2. - eg*M;

	tau = Q2h/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2h, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	melem2 = melem_2nd.melem2(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = melem_2nd.melem2_inv2(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);

	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared

//check if melem2 is bigger than 0
	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
//	if (melem2 < 0.) melem2 = 0.;
//	if (std::isinf(melem2) != 0) std::cout<<"Warning! the matrix element squared is infinite\n";

//	std::cout.precision(18);
//	if (eg < 1e-1.0) std::cout <<"melem2 = "<< melem2 <<"\t"<< f1 <<"\t"<< f2 <<"\t"<<
//			Q2h <<"\t"<< Q2e <<"\t"<< Q2k <<"\t"<< S <<"\t"<< Sq2 <<"\t"<<
//			Sk <<"\t"<< l1k1 <<"\t"<< l1k2 <<"\t"<< k1k2 <<"\n\n";

	x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));

	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);

//	if (s < 0.) s = 0.;

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

//		long double test = Q2h - Q2e + 2.*(-l1k1-l1k2+l2k1+l2k2+k1k2);

// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
	return s/2.;
}

long double Cross_Sections::crsect_brems_2nd_phig1 (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared

	Sk = Sq2 - Q2k/2. - eg*M;

	tau = Q2h/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2h, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	melem2 = melem_2nd.melem2(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = melem_2nd.melem2_inv2(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);

	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared

//check if melem2 is bigger than 0
	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
//	if (melem2 < 0.) melem2 = 0.;
//	if (std::isinf(melem2) != 0) std::cout<<"Warning! the matrix element squared is infinite\n";

	x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));

	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);

//	if (s < 0.) s = 0.;

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
	return s/2.;
}

long double Cross_Sections::crsect_brems_2nd_l1k1_1 (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared

	Sk = Sq2 - Q2k/2. - eg*M;

	tau = Q2h/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2h, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	melem2 = melem_2nd.melem2_l1k1_1(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = melem_2nd.melem2_inv2(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);

	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared

//check if melem2 is bigger than 0
//	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
//	if (melem2 < 0.) melem2 = 0.;
//	if (std::isinf(melem2) != 0) std::cout<<"Warning! the matrix element squared is infinite\n";

	x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));

	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);

//	if (s < 0.) s = 0.;

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
	return s/2.;
}

long double Cross_Sections::crsect_brems_2nd_l1k2_1 (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared
	Sk = Sq2 - Q2k/2. - eg*M;

	tau = Q2h/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2h, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	melem2 = melem_2nd.melem2_l1k2_1(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = melem_2nd.melem2_inv2(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);

	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared

//check if melem2 is bigger than 0
//	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
//	if (melem2 < 0.) melem2 = 0.;
//	if (std::isinf(melem2) != 0) std::cout<<"Warning! the matrix element squared is infinite\n";

//	std::cout.precision(18);
//	if (eg < 1e-1.0) std::cout <<"melem2 = "<< melem2 <<"\t"<< f1 <<"\t"<< f2 <<"\t"<<
//			Q2h <<"\t"<< Q2e <<"\t"<< Q2k <<"\t"<< S <<"\t"<< Sq2 <<"\t"<<
//			Sk <<"\t"<< l1k1 <<"\t"<< l1k2 <<"\t"<< k1k2 <<"\n\n";

	x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));

	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);

//	if (s < 0.) s = 0.;

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

//		long double test = Q2h - Q2e + 2.*(-l1k1-l1k2+l2k1+l2k2+k1k2);

// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
	return s/2.;
}

long double Cross_Sections::crsect_brems_2nd_l2k1_1 (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared

	Sk = Sq2 - Q2k/2. - eg*M;

	tau = Q2h/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2h, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	melem2 = melem_2nd.melem2_l2k1_1(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = melem_2nd.melem2_inv2(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);

	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared

//check if melem2 is bigger than 0
//	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
//	if (melem2 < 0.) melem2 = 0.;
//	if (std::isinf(melem2) != 0) std::cout<<"Warning! the matrix element squared is infinite\n";

	x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));

	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);

//	if (s < 0.) s = 0.;

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
	return s/2.;
}

long double Cross_Sections::crsect_brems_2nd_l2k2_1 (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared
	Sk = Sq2 - Q2k/2. - eg*M;

	tau = Q2h/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2h, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	melem2 = melem_2nd.melem2_l2k2_1(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = melem_2nd.melem2_inv2(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);

	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared

//check if melem2 is bigger than 0
//	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
//	if (melem2 < 0.) melem2 = 0.;
//	if (std::isinf(melem2) != 0) std::cout<<"Warning! the matrix element squared is infinite\n";

//	std::cout.precision(18);
//	if (eg < 1e-1.0) std::cout <<"melem2 = "<< melem2 <<"\t"<< f1 <<"\t"<< f2 <<"\t"<<
//			Q2h <<"\t"<< Q2e <<"\t"<< Q2k <<"\t"<< S <<"\t"<< Sq2 <<"\t"<<
//			Sk <<"\t"<< l1k1 <<"\t"<< l1k2 <<"\t"<< k1k2 <<"\n\n";

	x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));

	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);

//	if (s < 0.) s = 0.;

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

//		long double test = Q2h - Q2e + 2.*(-l1k1-l1k2+l2k1+l2k2+k1k2);

// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
	return s/2.;
}

//long double Cross_Sections::crsect_brems_2nd_l2k2_1 (const long double en1, const long double thl, const long double eg, const long double thg,
//		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{
//
//	long double Q2h, Q2k, Q2e, S, Sq2, Sk;
//
//	l2 = sqrt(pow(en1,2.)-m2);
//
//	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
//	long double alpha2 = -param->l1*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
//
//	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
//	S = en1 * M;
//	Sq2 = M * (param->en - en1);
//
//	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
//			+ cos(thg1)*cos(thg);
//	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
//	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);
//
//	long double l2k1 = en1*eg - l2*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
//	long double l2k2 = en1*eg1 - l2*eg1*cos(thg1);
//	long double l1k1 = param->en*eg - eg*param->l1*cospsi;
//	long double l1k2 = param->en*eg1 - eg1*param->l1*cospsi1;
//	long double k1k2 = eg1 * eg * (1. - cosxi);
//
//	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
//	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared
//	Sk = Sq2 - Q2k/2. - eg*M;
//
//	tau = Q2h/(4.*M2);
//
////Pauli and Dirac Form Factors
//	ff.ffactp (Q2h, gpe, gpm);
//	f2=(gpm-gpe)/(1.+tau);
//	f1=(gpe+tau*gpm)/(1.+tau);
//
//	melem2 = melem_2nd.melem2_l2k2_1(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
////	melem2 = melem_2nd.melem2_inv2(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);
//
//	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared
//
////check if melem2 is bigger than 0
////	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
////	if (melem2 < 0.) melem2 = 0.;
////	if (std::isinf(melem2) != 0) std::cout<<"Warning! the matrix element squared is infinite\n";
//
////	std::cout.precision(18);
////	if (eg < 1e-1.0) std::cout <<"melem2 = "<< melem2 <<"\t"<< f1 <<"\t"<< f2 <<"\t"<<
////			Q2h <<"\t"<< Q2e <<"\t"<< Q2k <<"\t"<< S <<"\t"<< Sq2 <<"\t"<<
////			Sk <<"\t"<< l1k1 <<"\t"<< l1k2 <<"\t"<< k1k2 <<"\n\n";
//
////	std::cout << melem2 << "\n";
//
//	x = std::abs(-alpha1*cos(phig) + alpha2*sin(phig));
//
//	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);
//
////	std::cout << s << "\n";
//
////	if (s < 0.) s = 0.;
//
////check if s is a real number
//	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";
//
////		long double test = Q2h - Q2e + 2.*(-l1k1-l1k2+l2k1+l2k2+k1k2);
//
//// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
//	return s/2.;
//}


long double Cross_Sections::crsect_brems_2nd_l1k1_2 (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared

	Sk = Sq2 - Q2k/2. - eg*M;

	tau = Q2h/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2h, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	melem2 = melem_2nd.melem2_l1k1_2(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = melem_2nd.melem2_inv2(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);

	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared

//check if melem2 is bigger than 0
//	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
//	if (melem2 < 0.) melem2 = 0.;
//	if (std::isinf(melem2) != 0) std::cout<<"Warning! the matrix element squared is infinite\n";

	x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));

	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);

//	if (s < 0.) s = 0.;

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
	return s/2.;
}

long double Cross_Sections::crsect_brems_2nd_l1k2_2 (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared
	Sk = Sq2 - Q2k/2. - eg*M;

	tau = Q2h/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2h, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	melem2 = melem_2nd.melem2_l1k2_2(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = melem_2nd.melem2_inv2(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);

	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared

//check if melem2 is bigger than 0
//	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
//	if (melem2 < 0.) melem2 = 0.;
//	if (std::isinf(melem2) != 0) std::cout<<"Warning! the matrix element squared is infinite\n";

//	std::cout.precision(18);
//	if (eg < 1e-1.0) std::cout <<"melem2 = "<< melem2 <<"\t"<< f1 <<"\t"<< f2 <<"\t"<<
//			Q2h <<"\t"<< Q2e <<"\t"<< Q2k <<"\t"<< S <<"\t"<< Sq2 <<"\t"<<
//			Sk <<"\t"<< l1k1 <<"\t"<< l1k2 <<"\t"<< k1k2 <<"\n\n";

	x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));

	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);

//	if (s < 0.) s = 0.;

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

//		long double test = Q2h - Q2e + 2.*(-l1k1-l1k2+l2k1+l2k2+k1k2);

// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
	return s/2.;
}

long double Cross_Sections::crsect_brems_2nd_l2k1_2 (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared

	Sk = Sq2 - Q2k/2. - eg*M;

	tau = Q2h/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2h, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	melem2 = melem_2nd.melem2_l2k1_2(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = melem_2nd.melem2_inv2(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);

	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared

//check if melem2 is bigger than 0
//	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
//	if (melem2 < 0.) melem2 = 0.;
//	if (std::isinf(melem2) != 0) std::cout<<"Warning! the matrix element squared is infinite\n";

	x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));

	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);

//	if (s < 0.) s = 0.;

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
	return s/2.;
}

long double Cross_Sections::crsect_brems_2nd_l2k2_2 (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared
	Sk = Sq2 - Q2k/2. - eg*M;

	tau = Q2h/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2h, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

	melem2 = melem_2nd.melem2_l2k2_2(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = melem_2nd.melem2_inv2(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);

	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared

//check if melem2 is bigger than 0
//	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
//	if (melem2 < 0.) melem2 = 0.;
//	if (std::isinf(melem2) != 0) std::cout<<"Warning! the matrix element squared is infinite\n";

//	std::cout.precision(18);
//	if (eg < 1e-1.0) std::cout <<"melem2 = "<< melem2 <<"\t"<< f1 <<"\t"<< f2 <<"\t"<<
//			Q2h <<"\t"<< Q2e <<"\t"<< Q2k <<"\t"<< S <<"\t"<< Sq2 <<"\t"<<
//			Sk <<"\t"<< l1k1 <<"\t"<< l1k2 <<"\t"<< k1k2 <<"\n\n";

	x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));

	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);

//	if (s < 0.) s = 0.;

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

//		long double test = Q2h - Q2e + 2.*(-l1k1-l1k2+l2k1+l2k2+k1k2);

// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
	return s/2.;
}

//long double Cross_Sections::crsect_brems_2nd_l2k2_2 (const long double en1, const long double thl, const long double eg, const long double thg,
//		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{
//
//	long double Q2h, Q2k, Q2e, S, Sq2, Sk;
//
//	l2 = sqrt(pow(en1,2.)-m2);
//
//	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
//	long double alpha2 = -param->l1*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);
//
//	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
//	S = en1 * M;
//	Sq2 = M * (param->en - en1);
//
//	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
//			+ cos(thg1)*cos(thg);
//	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
//	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);
//
//	long double l2k1 = en1*eg - l2*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
//	long double l2k2 = en1*eg1 - l2*eg1*cos(thg1);
//	long double k1k2 = eg1 * eg * (1. - cosxi);
//	long double l1k1 = param->en*eg - eg*param->l1*cospsi;
//	long double l1k2 = param->en*eg1 - eg1*param->l1*cospsi1;
//
//	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
//	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared
//	Sk = Sq2 - Q2k/2. - eg*M;
//
//	tau = Q2h/(4.*M2);
//
////Pauli and Dirac Form Factors
//	ff.ffactp (Q2h, gpe, gpm);
//	f2=(gpm-gpe)/(1.+tau);
//	f1=(gpe+tau*gpm)/(1.+tau);
//
//	melem2 = melem_2nd.melem2_l2k2_2(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
////	melem2 = melem_2nd.melem2_inv2(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);
//
//	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared
//
////check if melem2 is bigger than 0
////	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
////	if (melem2 < 0.) melem2 = 0.;
////	if (std::isinf(melem2) != 0) std::cout<<"Warning! the matrix element squared is infinite\n";
//
////	std::cout.precision(18);
////	if (eg < 1e-1.0) std::cout <<"melem2 = "<< melem2 <<"\t"<< f1 <<"\t"<< f2 <<"\t"<<
////			Q2h <<"\t"<< Q2e <<"\t"<< Q2k <<"\t"<< S <<"\t"<< Sq2 <<"\t"<<
////			Sk <<"\t"<< l1k1 <<"\t"<< l1k2 <<"\t"<< k1k2 <<"\n\n";
//
//	x = std::abs(-alpha1*cos(phig) + alpha2*sin(phig));
//
//	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);
//
////	if (s < 0.) s = 0.;
//
////check if s is a real number
//	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";
//
////		long double test = Q2h - Q2e + 2.*(-l1k1-l1k2+l2k1+l2k2+k1k2);
//
//// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
//	return s/2.;
//}

long double Cross_Sections::crsect_brems_2nd_add (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	long double p20 = param->en + M - en1 - eg - eg1;
	long double p21 = -l2*sin(thl) - eg*sin(thg)*cos(phig) - eg1*sin(thg1)*cos(phig1);
	long double p22 = - eg*sin(thg)*sin(phig) - eg1*sin(thg1)*sin(phig1);
	long double p23 = param->l1 - l2*cos(thl) - eg*cos(thg) - eg1*cos(thg1);

	long double l1k1_eg1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2_eg1 = param->en - param->l1*cos(thg1);
	long double k1k2_eg1 = eg * (1. - cosxi);
	long double l2k1_eg1 = en1*eg - eg*l2*cospsi;
	long double l2k2_eg1 = en1 - l2*cospsi1;
	long double p1l1 = M*param->en;
	long double p1l2 = M*en1;
	long double p1k1 = M*eg;
	long double p1k2 = M;
	long double p1p2 = M*p20;
	long double l1p2 = param->en*p20 - param->l1*p23;
	long double l2p2 = en1*p20 - l2*sin(thl)*p21 - l2*cos(thl)*p23;
	long double k1p2 = eg*p20 - eg*sin(thg)*cos(phig)*p21 - eg*sin(thg)*sin(phig)*p22
			      - eg*cos(thg)*p23;
	long double p2k2 = p20 - sin(thg1)*cos(phig1)*p21 - sin(thg1)*sin(phig1)*p22
		      - cos(thg1)*p23;
	long double l1l2 = param->en*en1 - param->l1*l2*cos(thl);

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared
	Sk = Sq2 - Q2k/2. - eg*M;
	tau = Q2h/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2h, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

//	melem2 = 2.*melem_2nd.melem2_Rinterfl1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2)
//			 + melem_2nd.melem2_add(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = 2.*melem_2nd.melem2_Rinterfl1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);
//	melem2 = melem_2nd.melem2_add(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = melem_2nd.melem2_R2(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2)
//			 + 2.*melem_2nd.melem2_Rinterf(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
	melem2 = melem_2nd.melem2_R2(p1l1, p1l2, p1k1, p1k2, p1p2, l1p2, l2p2, k1p2, p2k2, l1l2, l1k1_eg1,
			l2k1_eg1, l1k2_eg1, l2k2_eg1, k1k2_eg1, eg1, f1, f2)/eg1
				+ 2.*melem_2nd.melem2_Rinterf(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);

	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared

//	std::cout << melem2 << "\n";

	x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));

	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);

//check if s is a real number
//	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";
	if (s != s )  s = 0.;

// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
	return s/2.;
}

long double Cross_Sections::crsect_brems_2nd_add_phig1 (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig);
	long double alpha2 = l2*eg1*sin(thl)*sin(thg1) + eg*eg1*sin(thg)*sin(thg1)*cos(phig);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

		long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
		long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
		long double k1k2 = eg1 * eg * (1. - cosxi);
		long double l2k1 = en1*eg - eg*l2*cospsi;
		long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

		long double p20 = param->en + M - en1 - eg - eg1;
		long double p21 = -l2*sin(thl) - eg*sin(thg)*cos(phig) - eg1*sin(thg1)*cos(phig1);
		long double p22 = - eg*sin(thg)*sin(phig) - eg1*sin(thg1)*sin(phig1);
		long double p23 = param->l1 - l2*cos(thl) - eg*cos(thg) - eg1*cos(thg1);

		long double l1k1_eg1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
		long double l1k2_eg1 = param->en - param->l1*cos(thg1);
		long double k1k2_eg1 = eg * (1. - cosxi);
		long double l2k1_eg1 = en1*eg - eg*l2*cospsi;
		long double l2k2_eg1 = en1 - l2*cospsi1;
		long double p1l1 = M*param->en;
		long double p1l2 = M*en1;
		long double p1k1 = M*eg;
		long double p1k2 = M;
		long double p1p2 = M*p20;
		long double l1p2 = param->en*p20 - param->l1*p23;
		long double l2p2 = en1*p20 - l2*sin(thl)*p21 - l2*cos(thl)*p23;
		long double k1p2 = eg*p20 - eg*sin(thg)*cos(phig)*p21 - eg*sin(thg)*sin(phig)*p22
				      - eg*cos(thg)*p23;
		long double p2k2 = p20 - sin(thg1)*cos(phig1)*p21 - sin(thg1)*sin(phig1)*p22
			      - cos(thg1)*p23;
		long double l1l2 = param->en*en1 - param->l1*l2*cos(thl);

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared
	Sk = Sq2 - Q2k/2. - eg*M;
	tau = Q2h/(4.*M2);

//Pauli and Dirac Form Factors
	ff.ffactp (Q2h, gpe, gpm);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);

//	melem2 = 2.*melem_2nd.melem2_Rinterfl1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2)
//				 + melem_2nd.melem2_add(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = melem_2nd.melem2_add(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
//	melem2 = 2.*melem_2nd.melem2_Rinterfl1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);
	//	melem2 = melem_2nd.melem2_R2(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2)
	//			 + 2.*melem_2nd.melem2_Rinterf(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2);
	melem2 = melem_2nd.melem2_R2(p1l1, p1l2, p1k1, p1k2, p1p2, l1p2, l2p2, k1p2, p2k2, l1l2, l1k1_eg1,
			l2k1_eg1, l1k2_eg1, l2k2_eg1, k1k2_eg1, eg1, f1, f2)/eg1
				+ 2.*melem_2nd.melem2_Rinterf(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2);

//	std::cout.precision(14);
//	std::cout <<"melem2 = "<< melem_2nd.melem2_Rinterfl1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2) <<"\t"<<
//			l1k1 <<"\t"<< l1k2 <<"\t"<< l2k1 <<"\t"<< l2k2 <<"\t"<< Q2e <<"\t"<< Q2h <<"\t"<< S
//			<<"\t"<< Sk <<"\t"<< Sq2 <<"\t"<< f1 <<"\t"<< f2 <<"\n\n";

	melem2 *= pow(4.*pi*alpha,4.)/pow(Q2h,2.)/4.; //final matrix element squared

	x = std::abs(- alpha1*cos(phig1) + alpha2*sin(phig1));

	s = 2. * melem2 * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);

//check if s is a real number
//	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";
	if (s != s )  s = 0.;

// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
	return s/2.;
}

// Bremsstrahlung Cross Section for the interference term between weak and gamma exchange
long double Cross_Sections::interf_brems_1st (const long double en1, const long double thl, const long double eg, const long double thg) const{

	double q2, Q2, l2k, Q2e, S1, S2;

	l2 = sqrt(pow(en1,2)-m2);

	a = M*(param->en-en1-eg)-param->en*en1+param->l1*l2*cos(thl)-
			param->en*eg+param->l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);
	phig = a/b;       //calculates the cosine of the phi angle of the photon using 4.2.7 and 4.2.8

	phig = acos(phig);

	l1k = param->en*eg - param->l1*eg*cos(thg); //the product between the 4-momenta of the initial electron and final photon

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S1 = en1 * M;
	S2 = M * (param->en - en1);

	q2 = 2.*M*(en1-param->en+eg);  //calculates the momentum transfer squared
	Q2 = -q2;
	tau = Q2/(4.*M2);

	double ga = -1./2.;
	double gv = -1./2. + 2.*param->sw2;

	if (param->flag[param->kappa_weak] == 1) {
		double kappa = VC.kappa_weak(Q2);
		gv = -1./2. + 2.*kappa*param->sw2;

		ff.ffz (Q2, gpze, gpzm, kappa, param->sw2);
	}
	else ff.ffz (Q2, gpze, gpzm, param->sw2);
	ff.ffgae (Q2, gae);
	ff.ffactp (Q2, gpe, gpm);

	//Pauli and Dirac Form Factors
	f2 = (gpm-gpe)/(1.+tau);
	f1 = (gpe+tau*gpm)/(1.+tau);
	f2z = (gpzm-gpze)/(1.+tau);
	f1z = (gpze+tau*gpzm)/(1.+tau);

	double sp1, sp2, sl2, sk;

	sp1 = M*param->l1/m;
	sl2 = (param->l1*en1 - l2*param->en*cos(thl))/m;
	sk = (param->l1 * eg - param->en * eg * cos(thg))/m;
	sp2 = (param->l1 * (param->en - en1 - eg + M) - param->en * (param->l1 - l2*cos(thl) - eg * cos(thg)))/m;

	//matrix element squared expressed in l1k, l2k, q2, eg, param->en and en1 (found in "bremsstrahlung_matrix_element.pdf" document)
	melem_interf = (m*pow(l1k,-2)*pow(M,-2)*pow(2*l1k - Q2 + Q2e,-2)*
		     (f2*(4*gae*gv*Q2*M2*(-2*l1k*(Q2 - Q2e)*
		              (3*Q2*sk + 4*S1*sk - 4*S2*sk - Q2*sl2 - Q2e*sl2 + 4*S2*sl2 -
		                2*Q2e*sp1 - 2*Q2e*sp2) +
		             4*(Q2*(2*sk - sl2) + 4*S2*(-sk + sl2) - Q2e*(sl2 + 2*(sp1 + sp2)))*
		              pow(l1k,2) + (Q2*sk + 2*(2*S1*sk + (sp1 + sp2)*m2))*
		              pow(Q2 - Q2e,2)) + 4*f1z*ga*M2*
		           (2*l1k*Q2*(Q2 - Q2e)*(Q2*(3*sk - sl2) - Q2e*(sk + sp1 - sp2)) +
		             4*(Q2*(3*sk - sl2 + sp1 - sp2) + Q2e*(sk + sl2 - sp1 + sp2))*
		              pow(l1k,3) - Q2*(Q2*sk + (sk + sl2 + sp1 - sp2)*m2)*
		              pow(Q2 - Q2e,2) + 2*pow(l1k,2)*
		              (2*Q2*Q2e*(2*sk - sl2 + 2*sp1 - 2*sp2) +
		                (-7*sk + 3*sl2 - sp1 + sp2)*pow(Q2,2) +
		                (sk + sl2 - sp1 + sp2)*pow(Q2e,2))) +
		          f2z*ga*(4*pow(l1k,3)*(-(Q2*
		                   (-4*S2*sl2 + Q2*(sk + sl2 + sp1 + sp2) +
		                     Q2e*(sk + sl2 + sp1 + sp2))) +
		                4*(Q2*(sk - sl2 + sp1 - sp2) + Q2e*(sk + sl2 - sp1 + sp2))*M2
		                ) + Q2*(-2*sk*(2*S1*(Q2 + 2*S1) + Q2*M2) +
		                m2*(Q2*(sk + sl2 - sp1 - sp2) -
		                   4*(S1*(sp1 + sp2) + (sp1 - sp2)*M2)))*pow(Q2 - Q2e,2) +
		             2*pow(l1k,2)*(-4*M2*
		                 (Q2*Q2e*(3*sl2 - 4*sp1 + 4*sp2) +
		                   (3*sk - 2*sl2 + sp1 - sp2)*pow(Q2,2) -
		                   (sk + sl2 - sp1 + sp2)*pow(Q2e,2)) +
		                Q2*(-16*S1*S2*sl2 +
		                   4*Q2*(S2*(-2*sl2 + sp1 + sp2) + S1*(-2*sk + sl2 + sp1 + sp2)) +
		                   4*Q2e*(S2*sl2 + S1*(2*sk + sl2 + sp1 + sp2)) +
		                   (sk + sl2 + sp1 + sp2)*pow(Q2,2) -
		                   (sk + sl2 + sp1 + sp2)*pow(Q2e,2))) +
		             4*l1k*Q2*(Q2 - Q2e)*(-2*Q2e*S1*sk - Q2e*S1*sl2 + 4*S1*S2*sl2 -
		                Q2e*S1*sp1 - Q2e*S1*sp2 - 2*S2*sp1*m2 - 2*S2*sp2*m2 +
		                Q2*(S1*(4*sk - sl2 - sp1 - sp2) + S2*(sl2 - sp1 - sp2) +
		                   (sp1 + sp2)*m2) +
		                (Q2*(3*sk - sl2) + Q2e*(sl2 - 2*sp1 + 2*sp2))*M2 +
		                4*sk*pow(S1,2)))) -
		       4*f1*M2*(-(gae*gv*Q2*
		             (-2*l1k*(Q2 - Q2e)*(3*Q2*sk + 4*S1*sk - 4*S2*sk - Q2*sl2 - Q2e*sl2 +
		                  4*S2*sl2 - 2*Q2e*sp1 - 2*Q2e*sp2) +
		               4*(Q2*(2*sk - sl2) + 4*S2*(-sk + sl2) - Q2e*(sl2 + 2*(sp1 + sp2)))*
		                pow(l1k,2) + (Q2*sk + 2*(2*S1*sk + (sp1 + sp2)*m2))*
		                pow(Q2 - Q2e,2))) +
		          f2z*ga*(-2*l1k*Q2*(Q2 - Q2e)*(Q2*(3*sk - sl2) - Q2e*(sk + sp1 - sp2)) -
		             4*(Q2*(3*sk - sl2 + sp1 - sp2) + Q2e*(sk + sl2 - sp1 + sp2))*
		              pow(l1k,3) + Q2*(Q2*sk + (sk + sl2 + sp1 - sp2)*m2)*
		              pow(Q2 - Q2e,2) + 2*pow(l1k,2)*
		              (2*Q2*Q2e*(-2*sk + sl2 - 2*sp1 + 2*sp2) +
		                (7*sk - 3*sl2 + sp1 - sp2)*pow(Q2,2) -
		                (sk + sl2 - sp1 + sp2)*pow(Q2e,2))) +
		          f1z*ga*(8*pow(l1k,3)*(-2*S2*sl2 + Q2e*sp1 + Q2*(-sk + sl2 + sp2) +
		                4*sk*M2) -
		             4*pow(l1k,2)*(4*Q2e*S1*sk + 2*Q2e*S1*sl2 + 2*Q2e*S2*sl2 - 8*S1*S2*sl2 +
		                2*Q2e*S1*sp1 + 2*Q2e*S1*sp2 +
		                Q2*(Q2e*(2*sk - sl2 + 2*sp1 - 2*sp2) +
		                   2*(S2*(-2*sl2 + sp1 + sp2) + S1*(-2*sk + sl2 + sp1 + sp2))) +
		                (Q2*(8*sk - 2*sl2) - 2*Q2e*(4*sk + sl2))*M2 +
		                (-3*sk + 2*sl2 + sp2)*pow(Q2,2) - sp1*pow(Q2e,2)) +
		             pow(Q2 - Q2e,2)*(m2*
		                 (2*Q2*sp1 + 4*S1*(sp1 + sp2) - 4*(sk + sl2)*M2) +
		                sk*(4*Q2*S1 - 2*Q2*M2 + pow(Q2,2) + 8*pow(S1,2))) +
		             2*l1k*(Q2 - Q2e)*(Q2*(Q2e*(sk + sp1 - sp2) +
		                   2*(S2*(-sl2 + sp1 + sp2) + S1*(-4*sk + sl2 + sp1 + sp2) -
		                      (sp1 + sp2)*m2)) +
		                (Q2*(6*sk - 2*sl2) - 2*Q2e*(2*sk + sl2))*M2 +
		                (-3*sk + sl2)*pow(Q2,2) +
		                2*(Q2e*S1*(2*sk + sl2 + sp1 + sp2) +
		                   2*(-2*S1*S2*sl2 + S2*(sp1 + sp2)*m2 - 2*sk*pow(S1,2))))))))/2.;

	melem_interf *= - 2.*gf*pow(4.*pi*alpha,2.)/(sqrt(2.)*q2);

	x = sin(thl) * sin(thg) * sin(phig);

	s =  melem_interf / (32. * pow(2.*pi,4.) * M * param->l1 * x); //final cross section formula calculated with 4.2.9

	if (param->flag[param->order] == 2) {
		s *=  1. + VC.d_brems_ee(Q2) + VC.d_vert_pol(Q2, f1, f2, f1z, f2z, gae, ga, gv);
	}

	if (s != s ) std::cout<<"Warning! The Cross Section is not a real number"<<"\n"; //checks if s is a real number

	return 2.*s;
}

// Bremsstrahlung Cross Section for the interference term between weak and gamma exchange
long double Cross_Sections::interf_brems_1st_ps2 (const long double en1, const long double thl, const long double eg, const long double thg, const long double phig) const{

	long double q2, Q2, l2k, Q2e, S1, S2;

	l2 = sqrt(pow(en1,2.)-m2);

	l1k = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S1 = en1 * M;
	S2 = M * (param->en - en1);

	q2 = 2.*M*(en1-param->en+eg);  //calculates the momentum transfer squared
	Q2 = -q2;
	tau = Q2/(4.*M2);

	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = VC.kappa_weak(Q2);
		gv = -1./2. + 2.*kappa*param->sw2;

		ff.ffz (Q2, gpze, gpzm, kappa, param->sw2);
	}
	else ff.ffz (Q2, gpze, gpzm, param->sw2);
	ff.ffgae (Q2, gae);
	ff.ffactp (Q2, gpe, gpm);

	//Pauli and Dirac Form Factors
	f2 = (gpm-gpe)/(1.+tau);
	f1 = (gpe+tau*gpm)/(1.+tau);
	f2z = (gpzm-gpze)/(1.+tau);
	f1z = (gpze+tau*gpzm)/(1.+tau);

	long double sp1, sp2, sl2, sk;

	sp1 = M*param->l1/m;
	sl2 = (param->l1*en1 - l2*param->en*cos(thl))/m;
	sk = (param->l1 * eg - param->en * eg * cos(thg))/m;
	sp2 = (param->l1 * (param->en - en1 - eg + M) - param->en * (param->l1 - l2*cos(thl) - eg * cos(thg)))/m;

	//matrix element squared expressed in l1k, l2k, q2, eg, param->en and en1 (found in "bremsstrahlung_matrix_element.pdf" document)
	melem_interf = (m*pow(l1k,-2.)*pow(M,-2.)*pow(2.*l1k - Q2 + Q2e,-2.)*
		     (f2*(4.*gae*gv*Q2*M2*(-2.*l1k*(Q2 - Q2e)*
		              (3.*Q2*sk + 4.*S1*sk - 4.*S2*sk - Q2*sl2 - Q2e*sl2 + 4.*S2*sl2 -
		                2.*Q2e*sp1 - 2.*Q2e*sp2) +
		             4.*(Q2*(2.*sk - sl2) + 4.*S2*(-sk + sl2) - Q2e*(sl2 + 2.*(sp1 + sp2)))*
		              pow(l1k,2.) + (Q2*sk + 2.*(2.*S1*sk + (sp1 + sp2)*m2))*
		              pow(Q2 - Q2e,2.)) + 4.*f1z*ga*M2*
		           (2.*l1k*Q2*(Q2 - Q2e)*(Q2*(3.*sk - sl2) - Q2e*(sk + sp1 - sp2)) +
		             4.*(Q2*(3.*sk - sl2 + sp1 - sp2) + Q2e*(sk + sl2 - sp1 + sp2))*
		              pow(l1k,3.) - Q2*(Q2*sk + (sk + sl2 + sp1 - sp2)*m2)*
		              pow(Q2 - Q2e,2.) + 2.*pow(l1k,2.)*
		              (2.*Q2*Q2e*(2.*sk - sl2 + 2.*sp1 - 2.*sp2) +
		                (-7*sk + 3.*sl2 - sp1 + sp2)*pow(Q2,2.) +
		                (sk + sl2 - sp1 + sp2)*pow(Q2e,2.))) +
		          f2z*ga*(4.*pow(l1k,3.)*(-(Q2*
		                   (-4.*S2*sl2 + Q2*(sk + sl2 + sp1 + sp2) +
		                     Q2e*(sk + sl2 + sp1 + sp2))) +
		                4.*(Q2*(sk - sl2 + sp1 - sp2) + Q2e*(sk + sl2 - sp1 + sp2))*M2
		                ) + Q2*(-2.*sk*(2.*S1*(Q2 + 2.*S1) + Q2*M2) +
		                m2*(Q2*(sk + sl2 - sp1 - sp2) -
		                   4.*(S1*(sp1 + sp2) + (sp1 - sp2)*M2)))*pow(Q2 - Q2e,2.) +
		             2.*pow(l1k,2.)*(-4.*M2*
		                 (Q2*Q2e*(3.*sl2 - 4.*sp1 + 4.*sp2) +
		                   (3.*sk - 2.*sl2 + sp1 - sp2)*pow(Q2,2.) -
		                   (sk + sl2 - sp1 + sp2)*pow(Q2e,2.)) +
		                Q2*(-1.6*S1*S2*sl2 +
		                   4.*Q2*(S2*(-2.*sl2 + sp1 + sp2) + S1*(-2.*sk + sl2 + sp1 + sp2)) +
		                   4.*Q2e*(S2*sl2 + S1*(2.*sk + sl2 + sp1 + sp2)) +
		                   (sk + sl2 + sp1 + sp2)*pow(Q2,2.) -
		                   (sk + sl2 + sp1 + sp2)*pow(Q2e,2.))) +
		             4.*l1k*Q2*(Q2 - Q2e)*(-2.*Q2e*S1*sk - Q2e*S1*sl2 + 4.*S1*S2*sl2 -
		                Q2e*S1*sp1 - Q2e*S1*sp2 - 2.*S2*sp1*m2 - 2.*S2*sp2*m2 +
		                Q2*(S1*(4.*sk - sl2 - sp1 - sp2) + S2*(sl2 - sp1 - sp2) +
		                   (sp1 + sp2)*m2) +
		                (Q2*(3.*sk - sl2) + Q2e*(sl2 - 2.*sp1 + 2.*sp2))*M2 +
		                4.*sk*pow(S1,2.)))) -
		       4.*f1*M2*(-(gae*gv*Q2*
		             (-2.*l1k*(Q2 - Q2e)*(3.*Q2*sk + 4.*S1*sk - 4.*S2*sk - Q2*sl2 - Q2e*sl2 +
		                  4.*S2*sl2 - 2.*Q2e*sp1 - 2.*Q2e*sp2) +
		               4.*(Q2*(2.*sk - sl2) + 4.*S2*(-sk + sl2) - Q2e*(sl2 + 2.*(sp1 + sp2)))*
		                pow(l1k,2.) + (Q2*sk + 2.*(2.*S1*sk + (sp1 + sp2)*m2))*
		                pow(Q2 - Q2e,2.))) +
		          f2z*ga*(-2.*l1k*Q2*(Q2 - Q2e)*(Q2*(3.*sk - sl2) - Q2e*(sk + sp1 - sp2)) -
		             4.*(Q2*(3.*sk - sl2 + sp1 - sp2) + Q2e*(sk + sl2 - sp1 + sp2))*
		              pow(l1k,3.) + Q2*(Q2*sk + (sk + sl2 + sp1 - sp2)*m2)*
		              pow(Q2 - Q2e,2.) + 2.*pow(l1k,2.)*
		              (2.*Q2*Q2e*(-2.*sk + sl2 - 2.*sp1 + 2.*sp2) +
		                (7*sk - 3.*sl2 + sp1 - sp2)*pow(Q2,2.) -
		                (sk + sl2 - sp1 + sp2)*pow(Q2e,2.))) +
		          f1z*ga*(8*pow(l1k,3.)*(-2.*S2*sl2 + Q2e*sp1 + Q2*(-sk + sl2 + sp2) +
		                4.*sk*M2) -
		             4.*pow(l1k,2.)*(4.*Q2e*S1*sk + 2.*Q2e*S1*sl2 + 2.*Q2e*S2*sl2 - 8*S1*S2*sl2 +
		                2.*Q2e*S1*sp1 + 2.*Q2e*S1*sp2 +
		                Q2*(Q2e*(2.*sk - sl2 + 2.*sp1 - 2.*sp2) +
		                   2.*(S2*(-2.*sl2 + sp1 + sp2) + S1*(-2.*sk + sl2 + sp1 + sp2))) +
		                (Q2*(8*sk - 2.*sl2) - 2.*Q2e*(4.*sk + sl2))*M2 +
		                (-3.*sk + 2.*sl2 + sp2)*pow(Q2,2.) - sp1*pow(Q2e,2.)) +
		             pow(Q2 - Q2e,2.)*(m2*
		                 (2.*Q2*sp1 + 4.*S1*(sp1 + sp2) - 4.*(sk + sl2)*M2) +
		                sk*(4.*Q2*S1 - 2.*Q2*M2 + pow(Q2,2.) + 8*pow(S1,2.))) +
		             2.*l1k*(Q2 - Q2e)*(Q2*(Q2e*(sk + sp1 - sp2) +
		                   2.*(S2*(-sl2 + sp1 + sp2) + S1*(-4.*sk + sl2 + sp1 + sp2) -
		                      (sp1 + sp2)*m2)) +
		                (Q2*(6*sk - 2.*sl2) - 2.*Q2e*(2.*sk + sl2))*M2 +
		                (-3.*sk + sl2)*pow(Q2,2.) +
		                2.*(Q2e*S1*(2.*sk + sl2 + sp1 + sp2) +
		                   2.*(-2.*S1*S2*sl2 + S2*(sp1 + sp2)*m2 - 2.*sk*pow(S1,2.))))))))/2.;

	melem_interf *= - 8.*gf*pow(4.*pi*alpha,2.)/(4.*sqrt(2.)*q2);

	long double cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double A = param->l1*cos(thl)-eg*cospsi;
	long double B = param->en+M-eg;

	x = std::abs(A*en1-B*l2)/eg/pow(l2,2.);

	s =  melem_interf / (32. * pow(2.*pi,4.) * M * param->l1 * x);  //final cross section formula calculated with 4.29

	if (s != s ) std::cout<<"Warning! The Cross Section is not a real number"<<"\n"; //checks if s is a real number

	if (param->flag[param->order] == 2.) {
		s *=  1. + VC.d_brems_ee(Q2) + VC.d_vert(Q2, f1, f2)/2. + VC.d_vert_pol_Z(Q2, f1, f2, f1z, f2z, gae, ga, gv);
	}

	return s;
}

// Bremsstrahlung Cross Section for the interference term between weak and gamma exchange
long double Cross_Sections::interf_brems_1st_test (const long double en1, const long double thl, const long double eg, const long double thg) const{

	long double q2, Q2, l2k, Q2e, S1, S2;

	l2 = sqrt(pow(en1,2.)-m2);

	a = M*(param->en-en1-eg)-param->en*en1+param->l1*l2*cos(thl)-
			param->en*eg+param->l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);
	phig = a/b;       //calculates the cosine of the phi angle of the photon using 4.27 and 4.28

	phig = acos(phig);

	l1k = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S1 = en1 * M;
	S2 = M * (param->en - en1);

	q2 = 2.*M*(en1-param->en+eg);  //calculates the momentum transfer squared
	Q2 = -q2;
	tau = Q2/(4.*M2);

	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = VC.kappa_weak(Q2);
		gv = -1./2. + 2.*kappa*param->sw2;

		ff.ffz (Q2, gpze, gpzm, kappa, param->sw2);
	}
	else ff.ffz (Q2, gpze, gpzm, param->sw2);
	ff.ffgae (Q2, gae);
	ff.ffactp (Q2, gpe, gpm);

	//Pauli and Dirac Form Factors
	f2 = (gpm-gpe)/(1.+tau);
	f1 = (gpe+tau*gpm)/(1.+tau);
	f2z = (gpzm-gpze)/(1.+tau);
	f1z = (gpze+tau*gpzm)/(1.+tau);

	//matrix element squared expressed in l1k, l2k, q2, eg, param->en and en1
	melem_interf = -(pow(l1k,-2.)*pow(M,-2.)*pow(2.*l1k - Q2 + Q2e,-2.)*
		     (f2*Q2*(-4.*gae*gv*M2*
		           (8*(Q2 + Q2e - 2.*S2)*pow(l1k,3.) +
		             l1k*(Q2 - Q2e)*(4.*Q2*S1 + Q2e*(-Q2e + 4.*(S1 + S2)) +
		                2.*(3.*Q2 + 2.*Q2e + 4.*S1 - 8*S2)*m2 + pow(Q2,2.)) -
		             2.*pow(l1k,2.)*(4.*Q2*(S1 - S2) + Q2e*(-3.*Q2e + 4.*S1 + 8*S2) +
		                2.*(3.*Q2 + Q2e - 8*S2)*m2 + 3.*pow(Q2,2.)) -
		             (Q2 + Q2e - 4.*S2)*m2*pow(Q2 - Q2e,2.)) -
		          4.*f1z*ga*M2*(16*pow(l1k,4.) -
		             16*pow(l1k,3.)*(Q2 - Q2e + m2) +
		             m2*pow(Q2 - Q2e,3.) +
		             pow(l1k,2.)*(-8*Q2*Q2e + 4.*(7*Q2 - 6*Q2e)*m2 +
		                6*pow(Q2,2.) + 6*pow(Q2e,2.)) +
		             l1k*(4.*(Q2 - Q2e)*m4 + Q2e*pow(Q2,2.) - pow(Q2,3.) -
		                Q2*pow(Q2e,2.) -
		                2.*m2*(-1.2*Q2*Q2e + 7*pow(Q2,2.) + 5*pow(Q2e,2.)) +
		                pow(Q2e,3.))) +
		          f2z*ga*(-32*pow(l1k,4.)*M2 +
		             16*pow(l1k,3.)*(Q2*(2.*S1 + S2) + 2.*(Q2 - Q2e)*M2 +
		                m2*(-S2 + 2.*M2)) -
		             m2*(2.*S1*(Q2e - 4.*S2) + Q2*(Q2e + 2.*S1 - 2.*S2) +
		                2.*(Q2 - Q2e)*M2)*pow(Q2 - Q2e,2.) -
		             4.*pow(l1k,2.)*(-2.*m2*
		                 (Q2*(-Q2e - 3.*S1 + S2) + S1*(Q2e + 4.*S2) -
		                   7*(Q2 - Q2e)*M2 + pow(Q2,2.)) +
		                M2*(-4.*Q2*Q2e + 3.*pow(Q2,2.) + 3.*pow(Q2e,2.)) +
		                2.*Q2*(-3.*Q2e*S1 - 2.*Q2e*S2 + 4.*S1*S2 + Q2*(3.*S1 + S2) +
		                   4.*pow(S1,2.) + 2.*pow(S2,2.))) -
		             2.*l1k*(Q2 - Q2e)*(8*m4*M2 -
		                M2*(pow(Q2,2.) + pow(Q2e,2.)) +
		                m2*(2.*Q2e*S2 + Q2*(-3.*Q2e - 4.*S1 + 6*S2) -
		                   2.*(7*Q2 - 6*Q2e)*M2 + 2.*pow(Q2,2.) -
		                   8*pow(S1,2.) - 8*pow(S2,2.)) -
		                2.*Q2*(Q2*S1 - Q2e*(S1 + S2) +
		                   2.*(2.*S1*S2 + 2.*pow(S1,2.) + pow(S2,2.)))))) +
		       4.*f1*M2*(Q2*(gae*gv*
		              (-8*(Q2 + Q2e - 2.*S2)*pow(l1k,3.) -
		                l1k*(Q2 - Q2e)*
		                 (4.*Q2*S1 + Q2e*(-Q2e + 4.*(S1 + S2)) +
		                   2.*(3.*Q2 + 2.*Q2e + 4.*S1 - 8*S2)*m2 + pow(Q2,2.)) +
		                2.*pow(l1k,2.)*(4.*Q2*(S1 - S2) +
		                   Q2e*(-3.*Q2e + 4.*S1 + 8*S2) +
		                   2.*(3.*Q2 + Q2e - 8*S2)*m2 + 3.*pow(Q2,2.)) +
		                (Q2 + Q2e - 4.*S2)*m2*pow(Q2 - Q2e,2.)) +
		             f2z*ga*(-1.6*pow(l1k,4.) +
		                16*pow(l1k,3.)*(Q2 - Q2e + m2) -
		                m2*pow(Q2 - Q2e,3.) +
		                pow(l1k,2.)*(8*Q2*Q2e - 4.*(7*Q2 - 6*Q2e)*m2 -
		                   6*pow(Q2,2.) - 6*pow(Q2e,2.)) +
		                l1k*(-4.*(Q2 - Q2e)*m4 - Q2e*pow(Q2,2.) + pow(Q2,3.) +
		                   Q2*pow(Q2e,2.) +
		                   2.*m2*(-1.2*Q2*Q2e + 7*pow(Q2,2.) + 5*pow(Q2e,2.)) -
		                   pow(Q2e,3.)))) +
		          f1z*ga*(16*pow(l1k,4.)*(-Q2 + 2.*M2) -
		             16*pow(l1k,3.)*(-(Q2*(Q2 - Q2e + 2.*S1 + S2)) +
		                2.*(Q2 - Q2e)*M2 +
		                m2*(-Q2 + S2 + 2.*M2)) -
		             m2*(2.*S1*(Q2e - 4.*S2) + 2.*Q2*(S1 - S2) -
		                2.*(Q2 - Q2e)*M2 + pow(Q2,2.))*pow(Q2 - Q2e,2.) +
		             l1k*(Q2 - Q2e)*(-4.*Q2*m4 -
		                2.*M2*(pow(Q2,2.) + pow(Q2e,2.)) +
		                2.*m2*(-2.*Q2e*S2 - 2.*Q2*(Q2e - 2.*S1 + 3.*S2) +
		                   (-1.4*Q2 + 8*Q2e)*M2 + 5*pow(Q2,2.) +
		                   8*pow(S1,2.) + 8*pow(S2,2.)) +
		                Q2*(4.*Q2*S1 - 4.*Q2e*(S1 + S2) + pow(Q2,2.) + pow(Q2e,2.) +
		                   8*(2.*S1*S2 + 2.*pow(S1,2.) + pow(S2,2.)))) +
		             2.*pow(l1k,2.)*(2.*m2*
		                 (2.*Q2*(2.*Q2e - 3.*S1 + S2) + 2.*S1*(Q2e + 4.*S2) +
		                   2.*(7*Q2 - 5*Q2e)*M2 - 5*pow(Q2,2.)) +
		                M2*(-8*Q2*Q2e + 6*pow(Q2,2.) + 6*pow(Q2e,2.)) -
		                Q2*(-4.*Q2*(Q2e - 3.*S1 - S2) - 4.*Q2e*(3.*S1 + 2.*S2) +
		                   3.*pow(Q2,2.) + 3.*pow(Q2e,2.) +
		                   8*(2.*S1*S2 + 2.*pow(S1,2.) + pow(S2,2.))))))))/2.;

	melem_interf *= - 8.*gf*pow(4.*pi*alpha,2.)/(4.*sqrt(2.)*q2);

	x = sin(thl) * sin(thg) * sin(phig);

	s =  melem_interf / (32. * pow(2.*pi,4.) * M * param->l1 * x); //final cross section formula calculated with 4.29

	if (param->flag[param->order] == 2.) {
		s *=  1. + VC.d_brems_ee(Q2e, en1) + VC.d_vert_pol(Q2, f1, f2, f1z, f2z, gae, ga, gv);
	}

	if (s != s ) std::cout<<"Warning! The Cross Section is not a real number"<<"\n"; //checks if s is a real number

	return 2.*s;
}

// Bremsstrahlung Cross Section for the interference term between weak and gamma exchange
long double Cross_Sections::interf_brems_2nd (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared
	long double q2h = - Q2h;
	Sk = Sq2 - Q2k/2. - eg*M;

	tau = Q2h/(4.*M2);

	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = VC.kappa_weak(Q2h);
		gv = -1./2. + 2.*kappa*param->sw2;

		ff.ffz (Q2h, gpze, gpzm, kappa, param->sw2);
	}
	else ff.ffz (Q2h, gpze, gpzm, param->sw2);
	ff.ffgae (Q2e, gae);
	ff.ffactp (Q2e, gpe, gpm);

	//Pauli and Dirac Form Factors
	f2 = (gpm-gpe)/(1.+tau);
	f1 = (gpe+tau*gpm)/(1.+tau);
	f2z = (gpzm-gpze)/(1.+tau);
	f1z = (gpze+tau*gpzm)/(1.+tau);

	//matrix element squared expressed in l1k, l2k, q2, eg, param->en and en1 (found in "bremsstrahlung_matrix_element.pdf" document)
	melem_interf = melem_2nd_pol.melem_interf(l1k1, l1k2, k1k2, Q2e, Q2h, Q2k, S, Sk, Sq2, f1, f2, f1z, f2z, gae, ga, gv);

	melem_interf *= 8.*gf*pow(4.*pi*alpha,3.)/(4.*sqrt(2.)*q2h);

	x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));

	s =  melem_interf * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x); //final cross section formula

	if (s != s ) std::cout<<"Warning! The Cross Section is not a real number"<<"\n"; //checks if s is a real number

	return s;
}

long double Cross_Sections::crsect_brems_2nd_pol_add (const long double en1, const long double thl, const long double eg, const long double thg,
		const long double phig, const long double eg1, const long double thg1, const long double phig1) const{

	long double Q2h, Q2k, Q2e, S, Sq2, Sk;

	l2 = sqrt(pow(en1,2.)-m2);

	long double alpha1 = eg*eg1*sin(thg)*sin(thg1)*sin(phig1);
	long double alpha2 = l2*eg*sin(thl)*sin(thg) + eg*eg1*sin(thg)*sin(thg1)*cos(phig1);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S = en1 * M;
	Sq2 = M * (param->en - en1);

	long double cosxi = sin(thg1)*sin(thg)*cos(phig1)*cos(phig) + sin(thg1)*sin(thg)*sin(phig1)*sin(phig)
			+ cos(thg1)*cos(thg);
	cospsi = sin(thl)*sin(thg)*cos(phig) + cos(thl)*cos(thg);
	long double cospsi1 = sin(thl)*sin(thg1)*cos(phig1) + cos(thl)*cos(thg1);

	long double l1k1 = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	long double l1k2 = param->en*eg1 - param->l1*eg1*cos(thg1);
	long double k1k2 = eg1 * eg * (1. - cosxi);
	long double l2k1 = en1*eg - eg*l2*cospsi;
	long double l2k2 = en1*eg1 - eg1*l2*cospsi1;

	Q2k = Q2e + 2.*l1k1 - 2.*l2k1;
	Q2h = 2.*M*(param->en - en1 - eg - eg1);  //calculates the momentum transfer squared
	Sk = Sq2 - Q2k/2. - eg*M;
	tau = Q2h/(4.*M2);

	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = VC.kappa_weak(Q2h);
		gv = -1./2. + 2.*kappa*param->sw2;

		ff.ffz (Q2h, gpze, gpzm, kappa, param->sw2);
	}
	else
	ff.ffz (Q2h, gpze, gpzm, param->sw2);
	ff.ffgae (Q2e, gae);
	ff.ffactp (Q2e, gpe, gpm);

	//Pauli and Dirac Form Factors
	f2 = (gpm-gpe)/(1.+tau);
	f1 = (gpe+tau*gpm)/(1.+tau);
	f2z = (gpzm-gpze)/(1.+tau);
	f1z = (gpze+tau*gpzm)/(1.+tau);

	melem_interf = melem_2nd_pol.melem2_pol_add_interf_g(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2, f1z, f2z, gae, ga, gv)
					+ melem_2nd_pol.melem2_pol_add_interf_Z(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2, f1z, f2z, gae, ga, gv)
					+ melem_2nd_pol.melem2_pol_add(l1k1, l1k2, l2k1, l2k2, Q2e, Q2h, S, Sk, Sq2, f1, f2, f1z, f2z, gae, ga, gv);

	melem_interf *= 8.*gf*pow(4.*pi*alpha,3.)/(-4.*sqrt(2.)*Q2h); //final matrix element squared

	x = std::abs(- alpha1*cos(phig) + alpha2*sin(phig));

	s = melem_interf * eg * eg1 * l2 / (64. * pow(2.*pi,7.) * M * param->l1 * x);

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

// we need to divide by two to account for the symmetry factor, for having two identical particles in final state
	return s/2.;
}

//For carbon-1.2
//Born cross section
long double Cross_Sections::crsect_born_carbon (const long double Q2)const {

	long double fc = ff.ff_carbon12(Q2);
	long double S = param->en * M;

//Matrix element squared
	melem2 = -4.*(M2*Q2 + 2.*(Q2 - 2.*S)*S)*pow(fc,2.);

	melem2 *= pow(4.*pi*alpha,2.)/pow(Q2,2.);

	melem2 *= pow(param->Z_target,2.);

//Cross Section
	s = melem2 / (64.*M2*pi*(pow(param->en,2.)-m2));

	return s;
}

//Elastic Cross section
long double Cross_Sections::crsect_elastic_carbon (const long double Q2)const {

	long double fc = ff.ff_carbon12(Q2);
	long double S = param->en*M;

	alpha = VC.running_alpha(Q2);

//Matrix element squared
	melem2 = -4.*(M2*Q2 + 2.*(Q2 - 2.*S)*S)*pow(fc,2.);

	melem2 *= pow(4.*pi*alpha,2.)/pow(Q2,2.);

	melem2 *= pow(param->Z_target,2.);

//Cross Section
	sigma_born = melem2 / (64.*M2*pi*(pow(param->en,2.)-m2));

	s = sigma_born;

	if (param->flag[param->order] == 1) {
		s += sigma_born * (VC.d_vert_carbon(Q2) + VC.d_brems_ee(Q2));

//		if (param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_1st(Q2,me);
//		if (param->flag[param->vac_pol] == 2.)
//			s += sigma_born * (VC.d_vac_1st(Q2,me) + VC.d_vac_1st(Q2,m_mu) + VC.d_vac_1st(Q2,m_tau));
//		if (param->flag[param->vac_pol] >= 3.)
//			s += sigma_born * interpolation->d_vac_hadr(Q2);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_2nd(Q2,me);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] > 1)
//			s += sigma_born * (VC.d_vac_2nd(Q2,me) + VC.d_vac_2nd(Q2,m_mu) + VC.d_vac_2nd(Q2,m_tau));

	}

	return s;
}

// Born Interference Cross Section
long double Cross_Sections::interf_born_carbon (const long double Q2)const {

	long double S = param->en*M;
          // energy of the scattered electron

	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	long double fc = ff.ff_carbon12(Q2);
	long double fcPV = ff.ffz_carbon12(Q2);

	long double Qw;

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = param->kappa_weak;
		Qw = -24.*kappa*param->sw2;
		gv = -1./2. + 2.*kappa*param->sw2;
	}
	else Qw = -24.*param->sw2;

// Matrix Element for the interference term
	melem_interf = 8.*fc*fcPV*ga*Qw*(M2*Q2 + 2.*(Q2 - 2.*S)*S);

	melem_interf *= -gf*(4.*alpha*pi)/(sqrt(2.)*Q2);

	melem_interf *= param->Z_target;


//Cross Section
	s = melem_interf / (64.*M2*pi*(pow(param->en,2.)-m2));

	return s;
}

// Elastic Interference Cross Section
long double Cross_Sections::interf_elastic_carbon (const long double Q2)const {

	long double S = param->en*M;
          // energy of the scattered electron

	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	long double fc = ff.ff_carbon12(Q2);
	long double fcPV = ff.ffz_carbon12(Q2);

	alpha = VC.running_alpha(Q2);

	long double Qw;

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = param->kappa_weak;
		Qw = -24.*kappa*param->sw2;
		gv = -1./2. + 2.*kappa*param->sw2;
	}
	else Qw = -24.*param->sw2;

// Matrix Element for the interference term
	melem_interf = 8.*fc*fcPV*ga*Qw*(M2*Q2 + 2.*(Q2 - 2.*S)*S);

	melem_interf *= -gf*(4.*alpha*pi)/(sqrt(2.)*Q2);

	melem_interf *= param->Z_target;


//Cross Section
	sigma_born = melem_interf / (64.*M2*pi*(pow(param->en,2.)-m2));

	s = sigma_born;

	if (param->flag[param->order] == 1) {

		s += sigma_born * (VC.d_vert_carbon(Q2) + VC.d_brems_ee(Q2));

//		if (param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_1st(Q2,me) / 2.;
//		if (param->flag[param->vac_pol] == 2.)
//			s += sigma_born * (VC.d_vac_1st(Q2,me) + VC.d_vac_1st(Q2,m_mu) + VC.d_vac_1st(Q2,m_tau)) / 2.;
//		if (param->flag[param->vac_pol] >= 3.)
//			s += sigma_born * interpolation->d_vac_hadr(Q2) / 2.;
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_2nd(Q2,me);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] > 1)
//			s += sigma_born * (VC.d_vac_2nd(Q2,me) + VC.d_vac_2nd(Q2,m_mu) + VC.d_vac_2nd(Q2,m_tau)) / 2.;

	}

	return s;
}

// Bremsstrahlung Cross Section
long double Cross_Sections::crsect_brems_1st_carbon(const long double en1, const long double thl, const long double eg,
		const long double thg, const long double phig) const{

	long double q2, Q2, Q2e, S1, S2;

	l2 = sqrt(pow(en1,2.)-m2);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));

	long double S = 2.*param->en*M;
	long double U = -2.*en1*M;

	l1k = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	q2 = 2.*M*(en1 + eg - param->en);  //calculates the momentum transfer squared
	Q2 = -q2;

	long double fc = ff.ff_carbon12(Q2);

	//matrix element squared in l1k, Q'2. (Q2), Q2 (Q2e), Sq (S2) and S (S1)
	melem2 = -4.*pow(fc,2.)*pow(l1k,-2.)*pow(2.*l1k - Q2 + Q2e,-2.)*
			   (4.*(4.*M2*(-Q2 + Q2e) + Q2*(S - U))*
			      pow(l1k,3.) + 16*M2*pow(l1k,4.) -
			     m2*(M2*Q2 + (Q2 - U)*U)*pow(-Q2 + Q2e,2.) +
			     l1k*(-Q2 + Q2e)*
			      (2.*m2*(Q2 - 2.*U)*(Q2 - S - U) +
			        M2*(pow(Q2,2.) + pow(Q2e,2.)) +
			        Q2*(Q2e*S + (Q2 - U)*U - pow(S,2.))) +
			     pow(l1k,2.)*(M2*(-8*Q2*Q2e + 6*pow(Q2,2.) +
			           6*pow(Q2e,2.)) -
			        2.*Q2*(Q2*S - 2.*Q2*U + Q2e*(-2.*S + U) +
			           pow(S,2.) + pow(U,2.)) +
			        4.*m2*pow(-Q2 + S + U,2.)));

	melem2 *= pow(4.*pi*alpha,3.)/pow(q2,2.); //final matrix element squared

	melem2 *= pow(param->Z_target,2.);

//check if melem2 is bigger than 0
//	if (melem2 < 0.) std::cout<<"Warning! the matrix element squared is negative\n";
	if (melem2 < 0.) melem2 = 0.;

	x = sin(thl) * sin(thg) * sin(phig);

	s =  melem2 / (32. * pow(2.*pi,4.) * M * param->l1 * x); //final cross section formula calculated with 4.29

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

	//the cross section is multiplied by two to account for the range pi to 2.*pi of the azymuhtal angle of the photon
	return 2.*s;
}

// Bremsstrahlung Cross Section for the interference term between weak and gamma exchange
long double Cross_Sections::interf_brems_1st_carbon (const long double en1, const long double thl, const long double eg, const long double thg) const{

	long double q2, Q2, l2k, Q2e, S1, S2;

	l2 = sqrt(pow(en1,2.)-m2);

	a = M*(param->en-en1-eg)-param->en*en1+param->l1*l2*cos(thl)-
			param->en*eg+param->l1*eg*cos(thg)+en1*eg-eg*l2*cos(thl)*cos(thg)+m2;
	b = eg*l2*sin(thl)*sin(thg);
	phig = a/b;       //calculates the cosine of the phi angle of the photon using 4.27 and 4.28

	phig = acos(phig);

	l1k = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	long double S = 2.*param->en*M;
	long double U = -2.*en1*M;

	q2 = 2.*M*(en1-param->en+eg);  //calculates the momentum transfer squared
	Q2 = -q2;

	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	long double fc = ff.ff_carbon12(Q2);
	long double fcPV = ff.ffz_carbon12(Q2);

	long double Qw;

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = param->kappa_weak;
		Qw = -24.*kappa*param->sw2;
		gv = -1./2. + 2.*kappa*param->sw2;
	}
	else Qw = -24.*param->sw2;

	//matrix element squared expressed in l1k, l2k, q2, eg, param->en and en1 (found in "bremsstrahlung_matrix_element.pdf" document)
//	melem_interf = -4.*fc*fcPV*ga*Qw*pow(l1k,-2.)*
//			   pow(2.*l1k - Q2 + Q2e,-2.)*
//			   (-8*(-4.*M2*(-Q2 + Q2e) + Q2*(-S + U) +
//			        m2*(4.*M2 + S + U))*pow(l1k,3.) +
//			     32*M2*pow(l1k,4.) +
//			     2.*l1k*(-Q2 + Q2e)*
//			      (M2*(pow(Q2,2.) + pow(Q2e,2.)) +
//			        Q2*(Q2e*S + (Q2 - U)*U - pow(S,2.)) +
//			        m2*(M2*(14*Q2 - 8*Q2e) + 3.*Q2*S +
//			           5*Q2*U - 4.*S*U + Q2e*(-3.*Q2 + S + U) +
//			           2.*pow(Q2,2.) - 2.*pow(S,2.) - 4.*pow(U,2.)))
//			       - m2*pow(-Q2 + Q2e,2.)*
//			      (2.*M2*(-Q2 + Q2e) - Q2*S + Q2e*(Q2 - U) -
//			        2.*Q2*U + 2.*S*U + 2.*pow(U,2.)) -
//			     4.*pow(l1k,2.)*(M2*
//			         (4.*Q2*Q2e - 3.*pow(Q2,2.) - 3.*pow(Q2e,2.))\
//			         + Q2*(Q2*S - 2.*Q2e*S - 2.*Q2*U + Q2e*U +
//			           pow(S,2.) + pow(U,2.)) +
//			        m2*(2.*Q2*Q2e + 2.*M2*(-7*Q2 + 5*Q2e) -
//			           Q2*S - 4.*Q2*U + Q2e*U + 2.*S*U -
//			           2.*pow(Q2,2.) + 2.*pow(U,2.))));

	long double sp1 = M*param->l1/m;
	long double sl2 = (param->l1*en1 - l2*param->en*cos(thl))/m;
	long double sk = (param->l1 * eg - param->en * eg * cos(thg))/m;
	long double sp2 = (param->l1 * (param->en - en1 - eg + M) - param->en * (param->l1 - l2*cos(thl) - eg * cos(thg)))/m;

	melem_interf = -4*fc*fcPV*ga*m*Qw*pow(l1k,-2)*
			   pow(2*l1k - Q2 + Q2e,-2)*
			   (4*(8*M2*sk + Q2*sk + Q2*sl2 - 2*S*sl2 +
			        Q2*sp1 + Q2*sp2 +
			        Q2e*(sk + sl2 + sp1 + sp2) - 2*sl2*U)*
			      pow(l1k,3) - (m2*
			         (4*M2*(sk + sl2) +
			           Q2*(sk + sl2 - sp1 - sp2) +
			           2*(sp1 + sp2)*U) +
			        2*sk*(M2*Q2 + (Q2 - U)*U))*
			      pow(-Q2 + Q2e,2) +
			     2*pow(l1k,2)*(4*Q2*S*sl2 +
			        4*M2*(Q2*(-4*sk + sl2) +
			           Q2e*(4*sk + sl2)) - 2*Q2*S*sp1 -
			        2*Q2*S*sp2 - 4*Q2*sk*U + 6*Q2*sl2*U -
			        4*S*sl2*U + Q2e*
			         (-2*S*sl2 + 2*(2*sk + sp1 + sp2)*U) -
			        sk*pow(Q2,2) - sl2*pow(Q2,2) -
			        sp1*pow(Q2,2) - sp2*pow(Q2,2) +
			        (sk + sl2 + sp1 + sp2)*pow(Q2e,2) -
			        4*sl2*pow(U,2)) +
			     2*l1k*(-Q2 + Q2e)*
			      (2*M2*(Q2*(-3*sk + sl2) +
			           Q2e*(2*sk + sl2)) - 2*m2*S*sp1 -
			        2*m2*S*sp2 + 2*Q2e*sk*U + Q2e*sl2*U -
			        2*S*sl2*U - 2*m2*sp1*U + Q2e*sp1*U -
			        2*m2*sp2*U + Q2e*sp2*U +
			        Q2*(S*(sl2 - sp1 - sp2) +
			           2*(m2*(sp1 + sp2) + (-2*sk + sl2)*U))
			         + 2*sk*pow(U,2) - 2*sl2*pow(U,2)));

//	std::cout << melem_interf <<"\t"<< test << "\n";

	melem_interf *= -gf*pow(4.*alpha*pi,2.)/(sqrt(2.)*q2);

	melem_interf *= param->Z_target;

	x = sin(thl) * sin(thg) * sin(phig);

	s =  melem_interf / (32. * pow(2.*pi,4.) * M * param->l1 * x); //final cross section formula calculated with 4.29

	if (s != s ) std::cout<<"Warning! The Cross Section is not a real number"<<"\n"; //checks if s is a real number

	return 2.*s;
}

//For carbon-1.2
//Born cross section
long double Cross_Sections::crsect_born_carbon_thl (const long double Q2)const {

	long double fc = ff.ff_carbon12(Q2);
	long double S = param->en * M;
	long double en1 = param->en - Q2/(2.*M);

//Matrix element squared
	melem2 = -4.*(M2*Q2 + 2.*(Q2 - 2.*S)*S)*pow(fc,2.);
	melem2 *= pow(4.*pi*alpha,2.)/pow(Q2,2.);
	melem2 *= pow(param->Z_target,2.);

//Cross Section
	s = melem2/(64.*M2*pi2)*pow(en1/param->en,2.);

	return 2.*pi*s;
}

//Elastic Cross section
long double Cross_Sections::crsect_elastic_carbon_thl (const long double Q2)const {

	long double fc = ff.ff_carbon12(Q2);
	long double S = param->en * M;
	long double en1 = param->en - Q2/(2.*M);

//Matrix element squared
	melem2 = -4.*(M2*Q2 + 2.*(Q2 - 2.*S)*S)*pow(fc,2.);
	melem2 *= pow(4.*pi*VC.running_alpha(Q2),2.)/pow(Q2,2.);
	melem2 *= pow(param->Z_target,2.);

//Cross Section
	sigma_born = melem2/(64.*M2*pi2)*pow(en1/param->en,2.);
	s = sigma_born;

	if (param->flag[param->order] == 1) {
		s += sigma_born * (VC.d_vert_carbon(Q2) + VC.d_brems_ee(Q2));

//		if (param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_1st(Q2,me);
//		if (param->flag[param->vac_pol] == 2.)
//			s += sigma_born * (VC.d_vac_1st(Q2,me) + VC.d_vac_1st(Q2,m_mu) + VC.d_vac_1st(Q2,m_tau));
//		if (param->flag[param->vac_pol] >= 3.)
//			s += sigma_born * interpolation->d_vac_hadr(Q2);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_2nd(Q2,me);
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] > 1)
//			s += sigma_born * (VC.d_vac_2nd(Q2,me) + VC.d_vac_2nd(Q2,m_mu) + VC.d_vac_2nd(Q2,m_tau));

	}


	return 2.*pi*s;
}

// Born Interference Cross Section
long double Cross_Sections::interf_born_carbon_thl (const long double Q2)const {

	long double en1 = param->en - Q2/(2.*M);
	long double S = param->en*M;

	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	long double fc = ff.ff_carbon12(Q2);
	long double fcPV = ff.ffz_carbon12(Q2);

	long double Qw;

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = param->kappa_weak;
		Qw = -24.*kappa*param->sw2;
		gv = -1./2. + 2.*kappa*param->sw2;
	}
	else Qw = -24.*param->sw2;

// Matrix Element for the interference term
	melem_interf = 8.*fc*fcPV*ga*Qw*(M2*Q2 + 2.*(Q2 - 2.*S)*S);
	melem_interf *= -gf*(4.*alpha*pi)/(sqrt(2.)*Q2);
	melem_interf *= param->Z_target;

//Cross Section
	s = melem_interf/(64.*M2*pi2)*pow(en1/param->en,2.);

	return 2.*pi*s;
}


// Elastic Interference Cross Section
long double Cross_Sections::interf_elastic_carbon_thl (const long double Q2)const {

	long double en1 = param->en - Q2/(2.*M);
	long double S = param->en*M;

	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	long double fc = ff.ff_carbon12(Q2);
	long double fcPV = ff.ffz_carbon12(Q2);

	long double Qw;

	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = param->kappa_weak;
		Qw = -24.*kappa*param->sw2;
		gv = -1./2. + 2.*kappa*param->sw2;
	}
	else Qw = -24.*param->sw2;

// Matrix Element for the interference term
	melem_interf = 8.*fc*fcPV*ga*Qw*(M2*Q2 + 2.*(Q2 - 2.*S)*S);
	melem_interf *= -gf*(4.*VC.running_alpha(Q2)*pi)/(sqrt(2.)*Q2);
	melem_interf *= param->Z_target;

//Cross Section
	sigma_born = melem_interf/(64.*M2*pi2)*pow(en1/param->en,2.);

	s = sigma_born;

	if (param->flag[param->order] == 1) {

		s += sigma_born * (VC.d_vert_carbon(Q2) + VC.d_brems_ee(Q2));

//		if (param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_1st(Q2,me) / 2.;
//		if (param->flag[param->vac_pol] == 2.)
//			s += sigma_born * (VC.d_vac_1st(Q2,me) + VC.d_vac_1st(Q2,m_mu) + VC.d_vac_1st(Q2,m_tau)) / 2.;
//		if (param->flag[param->vac_pol] >= 3.)
//			s += sigma_born * interpolation->d_vac_hadr(Q2) / 2.;
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] == 1)
//			s += sigma_born * VC.d_vac_2nd(Q2,me) / 2.;
//		if (param->flag[param->order] == 2. && param->flag[param->vac_pol] > 1)
//			s += sigma_born * (VC.d_vac_2nd(Q2,me) + VC.d_vac_2nd(Q2,m_mu) + VC.d_vac_2nd(Q2,m_tau)) / 2.;

	}

	return 2.*pi*s;
}

long double Cross_Sections::asymm_born_carbon_test (const long double Q2)const {

	long double en1 = param->en - Q2/(2.*M);
	long double S = param->en*M;
	l2 = sqrt(pow(en1,2.)-m2);

	long double ga = -1./2.;
	long double gv = -1./2. + 2.*param->sw2;

	long double sp1 = M*param->l1/m;
	long double sl2 = (param->l1*en1 - l2*param->en*cos(param->thl))/m;
	long double sp2 = (param->l1 * (param->en - en1 + M) - param->en * (param->l1 - l2*cos(param->thl)))/m;

//	long double sp1 = M*param->en/m;
//	long double sl2 = (param->en*en1 - l2*param->en*cos(thl))/m;
//	long double sp2 = (param->en * (param->en - en1 + M) - param->l1 * (param->l1 - l2*cos(thl)))/m;

	long double fc = ff.ff_carbon12(Q2);
	long double fcPV = ff.ffz_carbon12(Q2);

	long double Qw;
	if (param->flag[param->kappa_weak] == 1) {
		long double kappa = param->kappa_weak;
		Qw = -24.*kappa*param->sw2;
		gv = -1./2. + 2.*kappa*param->sw2;
	}
	else Qw = -24.*param->sw2;

//	double asymm = -gf*Q2*Qw*fcPV/(4.*sqrt(2.)*pi*alpha*param->Z_target*fc);

	double asymm = fcPV*ga*gf*Qw*m*Q2*(4.*M2*sl2+Q2*(sl2+sp1+sp2)-4.*S*(sp1+sp2))
							/(4.*sqrt(2.)*pi*alpha*param->Z_target*fc*(M2*Q2+2*S*(Q2-2*S)));
//	double asymm = fcPV*ga*gf*Qw*Q2*(m2*(4*M2+Q2)+2.*M2*Q2+4*S*(Q2-2*S))
//							/(4.*sqrt(2.)*pi*alpha*param->Z_target*fc*(M2*Q2+2*S*(Q2-2*S)));

//	std::cout << m*(4.*M2*sl2+Q2*(sl2+sp1+sp2)-4.*S*(sp1+sp2)) <<"\t"<< 2.*(M2*Q2+2*S*(Q2-2*S)) << "\n";

//	asymm *= -gf*(4.*alpha*pi)/(sqrt(2.)*Q2);
//	asymm *= param->Z_target;
//	asymm /= (64.*M2*pi2)*pow(en1/param->en,2.);
//	asymm *= 2.*pi;

	return asymm;
}

//Born cross section for electron target
long double Cross_Sections::crsect_born_e (const long double Q2)const {

	long double S = param->en * M;

	melem2 = 2.*(-4.*Q2*S - 2.*Q2*m2 - 2.*Q2*M2 +
     pow(Q2,2.) + 8*pow(S,2.));

	melem2 *= pow(4.*pi*alpha,2.)/pow(Q2,2.);

//Cross Section
	s = melem2 / (64.*M2*pi*(pow(param->en,2.)-m2));

	return s;
}

// Bremsstrahlung Cross Section for electron target
long double Cross_Sections::crsect_brems_1st_e(const long double en1, const long double thl, const long double eg,
		const long double thg, const long double phig) const{

	long double q2, Q2, Q2e, S1, S2;
	long double me = M;
	long double mu = m;

	l2 = sqrt(pow(en1,2.)-m2);

	Q2e = - 2. * (m2 - param->en*en1 + param->l1 * l2 * cos(thl));
	S1 = en1 * M;
	S2 = M * (param->en - en1);

//	long double S = 2.*param->en*M;
//	long double U = -2.*en1*M;

	long double p1k = param->en*eg - param->l1*eg*cos(thg); //the product between the 4.-momenta of the initial electron and final photon
	q2 = 2.*M*(en1 + eg - param->en);  //calculates the momentum transfer squared
	Q2 = -q2;
	long double Q2h = Q2;

	melem2 = 4.*pow(Q2e - 2.*S2,-2.)*pow(Q2h - 2.*S2,-2.)*
   (4.*Q2e*m4*pow(Q2e - Q2h,2.) -
     2.*m2*(-32*Q2e*Q2h*S1*S2 +
        4.*p1k*(Q2e - Q2h)*(Q2e - 2.*S2)*
         (Q2e - 4.*(S1 + S2)) +
        8*Q2h*S1*pow(Q2e,2.) +
        4.*Q2h*S2*pow(Q2e,2.) +
        16*S1*S2*pow(Q2e,2.) - 2.*Q2h*pow(Q2e,3.) -
        4.*S1*pow(Q2e,3.) - 4.*S2*pow(Q2e,3.) +
        pow(Q2e,4.) -
        2.*Q2e*pow(mu,2.)*pow(Q2e - Q2h,2.) -
        4.*Q2e*S1*pow(Q2h,2.) -
        8*Q2e*S2*pow(Q2h,2.) +
        16*S1*S2*pow(Q2h,2.) +
        3.*pow(Q2e,2.)*pow(Q2h,2.) -
        16*Q2e*Q2h*pow(S1,2.) +
        8*pow(Q2e,2.)*pow(S1,2.) +
        8*pow(Q2h,2.)*pow(S1,2.) +
        8*pow(p1k,2.)*pow(Q2e - 2.*S2,2.) -
        8*Q2e*Q2h*pow(S2,2.) +
        8*pow(Q2e,2.)*pow(S2,2.) +
        8*pow(Q2h,2.)*pow(S2,2.)) +
     (Q2e - 2.*S2)*(Q2h - 2.*S2)*
      (4.*p1k*Q2e*(Q2e - Q2h - 2.*(2.*S1 + S2)) +
        8*Q2e*pow(p1k,2.) -
        2.*pow(mu,2.)*
         (2.*Q2e*(Q2h - 2.*S2) - 4.*Q2h*S2 +
           pow(Q2e,2.) + pow(Q2h,2.) + 8*pow(S2,2.))
          + Q2e*(4.*Q2h*S1 - 4.*Q2e*(S1 + S2) +
           pow(Q2e,2.) + pow(Q2h,2.) +
           8*(2.*S1*S2 + 2.*pow(S1,2.) + pow(S2,2.)))));

	melem2 *= pow(4.*pi*alpha,3.)/pow(Q2e,2.); //final matrix element squared

//check if melem2 is bigger than 0
	if (melem2 < 0.) melem2 = 0.;

	x = std::abs(sin(thl)*sin(thg)*sin(phig));

	s =  melem2 / (32. * pow(2.*pi,4.) * M * param->l1 * x); //final cross section formula calculated with 4.29

//check if s is a real number
	if (s != s )  std::cout<<"Warning! The Bremsstralung Cross Section is not a number\n";

	long double sigma_brems = s;

	//the cross section is multiplied by two to account for the range pi to 2.*pi of the azymuthal angle of the photon
	return 2.*s;
}
