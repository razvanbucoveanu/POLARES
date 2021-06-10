/*
 * cuba_param.h
 *
 *  Created on: Dec 1, 2015
 *      Author: razvan
 */

#ifndef SRC_CUBA_PARAM_H_
#define SRC_CUBA_PARAM_H_

//#include "parameters.h"

namespace POLARES {

class Parameters;

class Cuba_parameters {

protected:
	const Parameters* param;

public:
	static const int NDIM_ELASTIC;
	static const int NDIM_brems_1st;
	static const int NDIM_brems_2nd;
	static const int NDIM_brems_2nd_1diff;
	static const int NDIM_brems_2nd_2diff;
	static const int NDIM_brems_1st_1diff;
	static const int NDIM_brems_1st_2diff;
	static const int NCOMP;
	static const int NVEC;
	static const int KEY;
	static const double EPSABS;
	static const int GRIDNO;

	char *STATEFILE_test;
	char *STATEFILE;

	int GRIDNO_elastic;
	int GRIDNO_brems;
	int GRIDNO_brems_hadr;
	int GRIDNO_brems_hadr_interf;
	int GRIDNO_brems_1st;
	int GRIDNO_brems_test;
	int GRIDNO_brems_interf;
	int GRIDNO_brems_2nd;
	int GRIDNO_brems_2nd_l1k1;
	int GRIDNO_brems_2nd_l1k2;
	int GRIDNO_brems_2nd_l2k1;
	int GRIDNO_brems_2nd_l2k2;
	int GRIDNO_brems_2nd_add_l1k1;
	int GRIDNO_brems_2nd_add_l1k2;
	int GRIDNO_brems_2nd_add_l2k1;
	int GRIDNO_brems_2nd_add_l2k2;
	int GRIDNO_brems_interf_2nd;
	int GRIDNO_brems_interf_2nd_1;
	int GRIDNO_brems_interf_2nd_2;
	int flags;
	int flags_brems;
	double SEED;
	long long int MINEVAL;
	long long int MAXEVAL_LO;
	long long int MAXEVAL_1st;
	long long int MAXEVAL_2nd;
	long long int MAXEVAL_2nd_add;
	double EPSREL;
	int no_cores;
	int *SPIN;
//	Vegas
	long long int NSTART, NINCREASE, NBATCH;
//	Suave
	long long int NNEW, NMIN;
	double FLATNESS;

	void set_param(const Parameters* param);

	//variables for numerical integration
	long long int neval;
	int comp, fail, nregions;

	Cuba_parameters();

};
}

#endif /* SRC_CUBA_PARAM_H_ */
