/*
 * form_factors.h
 *
 *  Created on: Dec 1, 2015
 *      Author: razvan
 */

#ifndef SRC_FORM_FACTORS_H_
#define SRC_FORM_FACTORS_H_

namespace POLARES {

class Form_factors {
public:

	enum FF_Type {FF_simple_dipole, FF_pol_dipole_bern, FF_FW, FF_static_limit,
		FF_user_defined, FF_carbon, FF_carbon_user_defined, FF_WRONG};

protected:
	int flag;

	//Friedrich-Walcher parametrization of the proton form factors
	//(J. Friedrich and Th. Walcher, Eur. Phys. J. A 17 (2003) 607-623)
	//For G^p_E
	static const double pEa10;// = 1.041;
	static const double pEa11;// = 0.765;
	static const double pEa20;// = -0.041;
	static const double pEa21;// = 6.2;
	static const double pEab;// = -0.23;
	static const double pEQb;// = 0.07;
	static const double pEsigmab;// = 0.27;

	//For G^p_M
	static const double pMa10;// = 1.002;
	static const double pMa11;// = 0.749;
	static const double pMa20;// = -0.002;
	static const double pMa21;// = 6.0;
	static const double pMab;// = -0.13;
	static const double pMQb;// = 0.35;
	static const double pMsigmab;// = 0.21;

	//Polynomial X Dipole Parametrization (Bernauer's dissertation, pp 181)
	//parameters for G^p_E
	static const double ae[8];// = {-0.4980, 5.4592, -34.7281, 114.3173, -262.9808, 329.1395, -227.3306, 66.6980};

	//parameters for G^p_M
	static double const am[8];// = {0.2472, -4.9123, 29.7509, -84.0430, 129.3256, -111.1068, 49.9753, -9.1659};

	//parameters for G^n_M
	static double const amn[10];// = {-1.9147, 6.47767, -17.32918, 31.80021, -37.18707, 27.52359, -12.81713, 3.63457, -0.57277, 0.03843};

	//parameters for G^n_E
	static const double a, b;// = 2.28409, b = 4.41942;

	//parameters for G^s_E
	static const double aes, bes;// = 0.32267, bes = 4.686;

	//parameters for G^s_M
	static const double ams, bms;// = 0.044, bms = 0.93;

	//parameters for G^A_Z
	static const double Gap, M_a;// = -1.267, M_a = 1.014;

public:

	bool change_flag(const int flag);

	void ffactp (const long double Q2, double& ge, double& gm)const;
	void ffactn (const long double Q2, double& gen, double& gmn)const;
	void sff (const long double Q2, double& gse, double& gsm)const;
	void ffgae (const long double Q2,double& gae)const;
	void ffz (const long double Q2, double& gpze, double& gpzm, const double sw2)const;
	void ffz (const long double Q2, double& gpze, double& gpzm, const double kappa, const double sw2)const;

	double ff_carbon12 (const long double Q2)const;
	double ffz_carbon12 (const long double Q2)const;

	Form_factors()
	:flag(FF_simple_dipole){
	}

	Form_factors (const int flag_ff)
	:flag(flag_ff){
	}
};

} //namespace POLARES



#endif /* SRC_FORM_FACTORS_H_ */
