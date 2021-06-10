#include "cuba_param.h"
#include "parameters.h"

namespace POLARES {

const int Cuba_parameters::NDIM_ELASTIC = 1;
const int Cuba_parameters::NDIM_brems_1st = 4;
const int Cuba_parameters::NDIM_brems_2nd = 7;
const int Cuba_parameters::NDIM_brems_2nd_1diff = 6;
const int Cuba_parameters::NDIM_brems_2nd_2diff = 5;
const int Cuba_parameters::NDIM_brems_1st_1diff = 3;
const int Cuba_parameters::NDIM_brems_1st_2diff = 2;
const int Cuba_parameters::NCOMP = 1;
const int Cuba_parameters::NVEC = 1;
const int Cuba_parameters::KEY = 0;
const double Cuba_parameters::EPSABS = 1e-200;
const int Cuba_parameters::GRIDNO = 0;

void Cuba_parameters::set_param(const Parameters* param) {

		this->param = param;
		 SEED = param->SEED;
		 flags_brems = param->flag[param->int_output];
		 MINEVAL = param->MINEVAL;
		 MAXEVAL_LO = param->MAXEVAL_LO;
		 MAXEVAL_1st = param->MAXEVAL_1st;
		 MAXEVAL_2nd = param->MAXEVAL_2nd;
		 MAXEVAL_2nd_add = param->MAXEVAL_2nd_add;
		 NSTART = param->NSTART;
		 NINCREASE = param->NINCREASE;
		 NBATCH = param->NBATCH;
		 NNEW = param->NNEW;
		 NMIN = param->NMIN;
		 NBATCH = param->NBATCH;
		 FLATNESS = param->FLATNESS;
		 EPSREL = param->EPSREL;
		 no_cores = param->no_cores;
		 GRIDNO_elastic = param->GRIDNO_elastic;
		 GRIDNO_brems = param->GRIDNO_brems;
		 GRIDNO_brems_hadr = param->GRIDNO_brems_hadr;
		 GRIDNO_brems_hadr_interf = param->GRIDNO_brems_hadr_interf;
		 GRIDNO_brems_1st = param->GRIDNO_brems_1st;
		 GRIDNO_brems_test = param->GRIDNO_brems_test;
		 GRIDNO_brems_interf = param->GRIDNO_brems_interf;
		 GRIDNO_brems_2nd = param->GRIDNO_brems_2nd;
		 GRIDNO_brems_2nd_l1k1 = param->GRIDNO_brems_2nd_l1k1;
		 GRIDNO_brems_2nd_l1k2 = param->GRIDNO_brems_2nd_l1k2;
		 GRIDNO_brems_2nd_l2k1 = param->GRIDNO_brems_2nd_l2k1;
		 GRIDNO_brems_2nd_l2k2 = param->GRIDNO_brems_2nd_l2k2;
		 GRIDNO_brems_2nd_add_l1k1 = param->GRIDNO_brems_2nd_add_l1k1;
		 GRIDNO_brems_2nd_add_l1k2 = param->GRIDNO_brems_2nd_add_l1k2;
		 GRIDNO_brems_2nd_add_l2k1 = param->GRIDNO_brems_2nd_add_l2k1;
		 GRIDNO_brems_2nd_add_l2k2 = param->GRIDNO_brems_2nd_add_l2k2;
		 GRIDNO_brems_interf_2nd = param->GRIDNO_brems_interf_2nd;
		 GRIDNO_brems_interf_2nd_1 = param->GRIDNO_brems_interf_2nd_1;
		 GRIDNO_brems_interf_2nd_2 = param->GRIDNO_brems_interf_2nd_2;

	}

Cuba_parameters::Cuba_parameters():
				GRIDNO_elastic(0),
				GRIDNO_brems(0),
				GRIDNO_brems_hadr(0),
				GRIDNO_brems_hadr_interf(0),
				GRIDNO_brems_1st(0),
				GRIDNO_brems_test(0),
				GRIDNO_brems_interf(0),
				GRIDNO_brems_2nd(0),
				GRIDNO_brems_2nd_l1k1(0),
				GRIDNO_brems_2nd_l1k2(0),
				GRIDNO_brems_2nd_l2k1(0),
				GRIDNO_brems_2nd_l2k2(0),
				GRIDNO_brems_2nd_add_l1k1(0),
				GRIDNO_brems_2nd_add_l1k2(0),
				GRIDNO_brems_2nd_add_l2k1(0),
				GRIDNO_brems_2nd_add_l2k2(0),
				GRIDNO_brems_interf_2nd(0),
				GRIDNO_brems_interf_2nd_1(0),
				GRIDNO_brems_interf_2nd_2(0),
				flags(0),
				STATEFILE(0),
				STATEFILE_test("test_file"),
				flags_brems(0),
				SEED(0),
				MINEVAL(0),
				MAXEVAL_LO(0),
				MAXEVAL_1st(0),
				MAXEVAL_2nd(0),
				MAXEVAL_2nd_add(0),
				EPSREL(0),
				no_cores(0),
				SPIN(0),
				NSTART(0),
				NINCREASE(0),
				NBATCH(0),
				NNEW(0),
				NMIN(0),
				FLATNESS(0),
				comp(0),
				neval(0),
				fail(0),
				nregions(0),
				param(0){
		};

}
