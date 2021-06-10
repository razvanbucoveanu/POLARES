/*
 * form_factors.cpp
 *
 *  Created on: Dec 1, 2015
 *      Author: razvan
 */

#include <cmath>
#include "form_factors.h"
#include "const.h"
#include <stdlib.h>
#include <cstdlib>
#include <iostream>


using namespace constants;
using namespace POLARES;

const double Form_factors::pEa10 = 1.041;
const double Form_factors::pEa11 = 0.765;
const double Form_factors::pEa20 = -0.041;
const double Form_factors::pEa21 = 6.2;
const double Form_factors::pEab = -0.23;
const double Form_factors::pEQb = 0.07;
const double Form_factors::pEsigmab = 0.27;

const double Form_factors::pMa10 = 1.002;
const double Form_factors::pMa11 = 0.749;
const double Form_factors::pMa20 = -0.002;
const double Form_factors::pMa21 = 6.0;
const double Form_factors::pMab = -0.13;
const double Form_factors::pMQb = 0.35;
const double Form_factors::pMsigmab = 0.21;

const double Form_factors::a = 2.28409;
const double Form_factors::b = 4.41942;

const double Form_factors::aes = 0.32267;
const double Form_factors::bes = 4.686;

const double Form_factors::ams = 0.044;
const double Form_factors::bms = 0.93;

const double Form_factors::Gap = -1.267;
const double Form_factors::M_a = 1.014;


bool Form_factors::change_flag(const int flag){
	if( (flag<0)||(flag>=FF_WRONG) ){
		std::cerr << "Wrong flag for Form Factors -> exit(1)"<<std::endl;
		exit(1);
		return 1;
	}else{
		this->flag =flag;
		return 0;
	}
}

const double Form_factors::ae[8] = {-0.4980, 5.4592, -34.7281, 114.3173, -262.9808, 329.1395, -227.3306, 66.6980};
const double Form_factors::am[8] = {0.2472, -4.9123, 29.7509, -84.0430, 129.3256, -111.1068, 49.9753, -9.1659};
const double Form_factors::amn[10] = {-1.9147, 6.47767, -17.32918, 31.80021, -37.18707, 27.52359, -12.81713, 3.63457, -0.57277, 0.03843};

void Form_factors::ffactp(const long double Q2, double& ge, double& gm)const {

	if (flag == FF_simple_dipole)
//Standard Dipole parameterization
	{
		ge = pow((1.+Q2/0.71),-2);
		gm = muep*ge;
	}

	if (flag == FF_pol_dipole_bern)
	// Polynomial x Dipole parameterization (Bernauer's dissertation, pp 181)
		{
			ge = 0.; gm = 0.;
			const double* itr_ae(ae);
			const double* itr_am(am);
			double Q2i(Q2);
			for (int i = 0; i <= 7; i++)
			{
	//			ge = ge+ae_L[i]*pow(Q2,i+1);
	//			gm = gm+am_L[i]*pow(Q2,i+1);
				ge += (*itr_ae)*Q2i;
				gm += (*itr_am)*Q2i;
				++itr_ae;
				++itr_am;
				Q2i *= Q2;
			}
			const double g(pow((1.+Q2/0.71),-2));
			ge = g*(1.+ge);
			gm = muep*g*(1.+gm);
		}

		if (flag == FF_FW)
	// Friedrich Walcher parameterization (J. Friedrich and Th. Walcher, Eur. Phys. J. A 17 (2003) 607-623)
		{
			ge = pEa10/pow((1.+Q2/pEa11),2)+pEa20/pow((1.+Q2/pEa21),2) +
			     pEab*Q2*(exp(-0.5*pow(((sqrt(Q2)-pEQb)/pEsigmab),2))+
		             exp(-0.5*pow(((sqrt(Q2)+pEQb)/pEsigmab),2)));

			gm = muep*(pMa10/pow((1.+Q2/pMa11),2)+pMa20/pow((1.+Q2/pMa21),2) +
			     pMab*Q2*(exp(-0.5*pow(((sqrt(Q2)-pMQb)/pMsigmab),2))+
			     exp(-0.5*pow(((sqrt(Q2)+pMQb)/pMsigmab),2))));
		}

		if (flag == FF_static_limit)
	// static limit of the form factors
		{
			ge = 1.;
			gm = muep;
//			gm = 0.;
		}

		if (flag == FF_user_defined)
	// user-defined form factors
		{
			std::cerr << "\nForm factors have not been defined by the user "
					"-> Go to 'form_factors.cpp' line 114\n";
			exit(1);
		}
}

// Neutron Form Factors (S. Galster et al., Nucl. Phys. B 32 (1971) 221-237)
void Form_factors::ffactn (const long double Q2, double& gen, double& gmn)const {

	const double gsd(pow((1.+Q2/0.71),-2));
	const double t(Q2/(4.*mpr2));

	gen = a*t*gsd/(1.+b*t);
	gmn = 0.;
	const double* itr_amn(amn);
	double Q2i(1.);

	for (int i = 0; i <= 9; i++)
	{
		gmn += (*itr_amn)*Q2i;
		Q2i*=Q2;
		++itr_amn;
	}
}

// Form Factors for strangeness (P. Young et al., Phys. Rev. C 79 (2009) 065202)
void Form_factors::sff (const long double Q2, double& gse, double& gsm)const {

	const double t(Q2/(4.*mpr2));
	const double gsd(pow((1.+Q2/0.71),-2.));
	gse = aes*t*gsd / (1.+bes*t);
	gsm = ams + bms*Q2;
}

// G_A^Z Form Factor (K. Hagiwara et al., Phys. Rev. D66, 010001 (2002))
void Form_factors::ffgae (const long double Q2,double& gae)const {

	gae = Gap*pow(1.+Q2/pow(M_a,2.),-2.);
}

// Form Factors for Z exchange (using isospin symmetry)
void Form_factors::ffz (const long double Q2, double& gpze, double& gpzm,
		const double sw2)const {

	double gpe, gpm, gne, gnm, gse, gsm;

	sff (Q2, gse, gsm);
	ffactp (Q2, gpe, gpm);
	ffactn (Q2, gne, gnm);

	gpze = (1.-4.*sw2)*gpe-gne-gse;
	gpzm = (1.-4.*sw2)*gpm-gnm-gsm;
//	gpze *= rho_ew;
//	gpzm *= rho_ew;
}

// Form Factors for Z exchange (using isospin symmetry)
void Form_factors::ffz (const long double Q2, double& gpze, double& gpzm,
		const double kappa, const double sw2)const {

	double gpe, gpm, gne, gnm, gse, gsm;

	sff (Q2, gse, gsm);
	ffactp (Q2, gpe, gpm);
	ffactn (Q2, gne, gnm);

	gpze = (1.-4.*kappa*sw2)*gpe-gne-gse;
	gpzm = (1.-4.*kappa*sw2)*gpm-gnm-gsm;
//	gpze *= rho_ew;
//	gpzm *= rho_ew;
}

double Form_factors::ff_carbon12(const long double Q2) const {

	double fc = 0.;

//	Symmetrized Fermi Form Factor for carbon-12 (J. Piekarewicz et al.,
//	Phys. Rev. C 94 (2016) 034316)
	if (flag == FF_carbon) {
		const double a = 0.4505686556360125 * fm;
		const double c = 2.3425250733227583 * fm;
		const double q = sqrt(Q2);

		fc = 3.*(pi*q*a/sinh(pi*q*a));
		fc *= pi*q*a*sin(q*c)/tanh(pi*q*a) - q*c*cos(q*c);
		fc /= q*c*(pow(q*c,2.)+pow(pi*q*a,2.));
	}

	if (flag == FF_carbon_user_defined)
// user-defined form factors
	{
		std::cerr << "\nCarbon Form factors have not been defined by the user "
				"-> Go to 'form_factors.cpp' line 189\n";
		exit(1);
	}

	return fc;
}

double Form_factors::ffz_carbon12(const long double Q2) const {

	return ff_carbon12(Q2);
}
