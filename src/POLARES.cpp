//------------------------------------------------------------------------------------------------
// POLARES (Radiative Corrections to Polarized Electron-Proton Scattering)
// is a free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by the
// Free Software Foundation (see <http://www.gnu.org/licenses/>).
// If you do so please cite
//	R.-D. Bucoveanu and H. Spiesberger, Eur. Phys. J. A (2019) 55: 57
//	arXiv:1811.04970 [hep-ph].
// Copyright (c) Razvan Bucoveanu, 2019. E-mail: rabucove@uni-mainz.de
//
// Modified: 24.02.-05.03.2022 by HS
//
//================================================================================================
#include "POLARES.h" // Main header file
#include <cuba.h>
#include <iostream>
#include <fstream>
#include "const.h"
#include <stdlib.h>
#include "config.h"
//------------------------------------------------------------------------------------------------

using namespace POLARES;
using namespace constants;

PES::PES()
:seed_is_set(false),
 m(me), M(mpr), m2(me2), M2(mpr2) {
	FS.seed = rand.uniform_int(1e7);
}

PES::~PES(){
}

void PES::set_input(const Input& input){
	this->input = input;
	param.read_input(input);

	m = param.m;
	m2 = pow(m,2.);
	M = param.M;
	M2 = pow(M,2.);

// Prepare to include hadronic vacuum polarization
	if (param.flag[param.vac_pol] == 3) {
		std::ifstream ifvpol("share/vp_Ignatov.dat");
		if (!ifvpol) {
			std::cerr<<"Error: 'vp_Ignatov.dat' cannot be found! "
					"Please copy it in the folder named 'share'.'.\n\n";
			exit (EXIT_FAILURE);
		}
		interpolation.init_d_vac_hadr("share/vp_Ignatov.dat");
	}
	if (param.flag[param.vac_pol] == 4) {
		std::ifstream ifvpol("share/vp_Jeger.dat");
		if (!ifvpol) {
			std::cerr<<"Error: 'vp_Jeger.dat' cannot be found! "
					"Please copy it in the folder named 'share'.'.\n\n";
			exit (EXIT_FAILURE);
		}
		interpolation.init_d_vac_hadr("share/vp_Jeger.dat");
	}
	if (param.flag[param.vac_pol] == 5) {
		std::ifstream ifvpol("share/vp_knt18.dat");
		if (!ifvpol) {
			std::cerr<<"Error: 'vp_knt18.dat' cannot be found! "
					"Please copy it in the folder named 'share'.'.\n\n";
			exit (EXIT_FAILURE);
		}
		interpolation.init_d_vac_hadr("share/vp_knt18.dat");
	}
// Prepare to include two-photon-exchange corrections
	if (param.flag[param.tpe] == 2) {
			std::ifstream iftpe("share/tpe_fwd.dat");
			if (!iftpe) {
				std::cerr<<"Error: 'tpe_fwd.dat' cannot be found! "
						"Please copy it in the folder named 'share'.\n\n";
				exit (EXIT_FAILURE);
			}
			interpolation.init_tpe_Tomalak("share/tpe_fwd.dat");
		}
	if (param.flag[param.tpe] == 3) {
			std::ifstream iftpe("share/tpe_bwd.dat");
			if (!iftpe) {
				std::cerr<<"Error: 'tpe_bwd.dat' cannot be found! "
						"Please copy it in the folder named 'share'.\n\n";
				exit (EXIT_FAILURE);
			}
			interpolation.init_tpe_Tomalak("share/tpe_bwd.dat");
		}
}


int PES::initialization(){

	param.final_param(input);
	CP.set_param(&param);
	integrands.set_param(&param, &interpolation);

	void *USERDATA = &integrands;

	cubacores(CP.no_cores, 10000);

	double integral[CP.NCOMP], error[CP.NCOMP], prob[CP.NCOMP];
	double sigma_hp_2nd_l1k1, sigma_hp_2nd_l1k2, sigma_hp_2nd_l2k1, sigma_hp_2nd_l2k2;
	double sigma_hp_2nd_l1k1_error, sigma_hp_2nd_l1k2_error, sigma_hp_2nd_l2k1_error,
	sigma_hp_2nd_l2k2_error;


// Proton target
	if (param.flag[param.target] == 0) {


// Leading order only
	if (param.flag[param.LO] == 1) {
		std::cout << "Integration for the LO unpolarized cross section ... ";
		if (param.flag[param.int_method] == 0) {
		llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_born,
				USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
				CP.MINEVAL, CP.MAXEVAL_LO, CP.NSTART, CP.NINCREASE, CP.NBATCH,
				CP.GRIDNO, CP.STATEFILE, CP.SPIN,
				&CP.neval, &CP.fail, integral, error, prob);
		}
		if (param.flag[param.int_method] == 1) {
				llSuave(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_born,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_LO, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
		}
		std::cout <<"completed\n";
		output.sigma_unpol_born = integral[CP.comp]*nb;
		errors.sigma_unpol_born = error[CP.comp]*nb;

        if (param.flag[param.asymmetry] == 1) {
            std::cout << "Integration for the LO polarized cross section ... ";
            llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_interf_born,
                    USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
                    CP.MINEVAL, CP.MAXEVAL_LO, CP.NSTART, CP.NINCREASE, CP.NBATCH,
                    CP.GRIDNO, CP.STATEFILE, CP.SPIN,
                    &CP.neval, &CP.fail, integral, error, prob);
            output.sigma_pol_born = param.P*integral[CP.comp]*nb;
            errors.sigma_pol_born = param.P*error[CP.comp]*nb;
            output.asymm_born = output.sigma_pol_born / output.sigma_unpol_born;
            errors.asymm_born = abs(output.sigma_pol_born/output.sigma_unpol_born) *
                    sqrt(pow(errors.sigma_unpol_born/output.sigma_unpol_born,2.)
                            + pow(errors.sigma_pol_born/output.sigma_pol_born,2.));
            output.sigma_born = output.sigma_pol_born + output.sigma_unpol_born;
            errors.sigma_born = sqrt(pow(errors.sigma_pol_born,2.) + 
                                     pow(errors.sigma_unpol_born,2.));
            std::cout <<"completed\n";
        }
    }
// End section for leading order only


// Non-radiative (elastic) cross section: 
// leading order + loop + soft photon  
// for both first order or first+second order 
// and including vacpol and hadronic corrections if requested
	if (param.flag[param.order] >= 0) {
		std::cout << "Integration for soft-photon + loop unpolarized cross section ... ";
		llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_elastic,
				USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
				CP.MINEVAL, CP.MAXEVAL_LO, CP.NSTART, CP.NINCREASE, CP.NBATCH,
				CP.GRIDNO_elastic, CP.STATEFILE, CP.SPIN,
				&CP.neval, &CP.fail, integral, error, prob);
		if (param.flag[param.order] == 0 || param.flag[param.order] == 1) {
			output.sigma_unpol_elastic_1st = integral[CP.comp]*nb;
            errors.sigma_unpol_elastic_1st = error[CP.comp]*nb;
		}
		if (param.flag[param.order] == 2) {
			output.sigma_unpol_elastic_2nd = integral[CP.comp]*nb;
            errors.sigma_unpol_elastic_2nd = error[CP.comp]*nb;
// ... but we want to keep track of the 1st order corrections: 
// reset order and repeat the calculation 
            param.flag[param.order] = 1;
            llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_elastic,
                    USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
                    CP.MINEVAL, CP.MAXEVAL_LO, CP.NSTART, CP.NINCREASE, CP.NBATCH,
                    CP.GRIDNO_elastic, CP.STATEFILE, CP.SPIN,
                    &CP.neval, &CP.fail, integral, error, prob);
            output.sigma_unpol_elastic_1st = integral[CP.comp]*nb;
            errors.sigma_unpol_elastic_1st = error[CP.comp]*nb;
            param.flag[param.order] = 2;
// done recalculation of 1st order elastic cross section 
		}
		if (param.flag[param.brems] == 0 && param.flag[param.brems_hadr] == 0) {
			output.sigma_unpol_1st = output.sigma_unpol_elastic_1st;
			errors.sigma_unpol_1st = errors.sigma_unpol_elastic_1st;
			output.sigma_unpol_2nd = output.sigma_unpol_elastic_2nd;
			errors.sigma_unpol_2nd = errors.sigma_unpol_elastic_2nd;
		}
		std::cout <<"completed\n";

        if (param.flag[param.asymmetry] == 1) {
			std::cout << "Integration for soft-photon + loop polarized cross section ... ";
			llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_interf_elastic,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_LO, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			if (param.flag[param.order] == 0 || param.flag[param.order] == 1) {
				output.sigma_pol_elastic_1st = param.P*integral[CP.comp]*nb;
                errors.sigma_pol_elastic_1st = param.P*error[CP.comp]*nb;
            }
			if (param.flag[param.order] == 2) {
				output.sigma_pol_elastic_2nd = param.P*integral[CP.comp]*nb;
                errors.sigma_pol_elastic_2nd = param.P*error[CP.comp]*nb; 
// keep track of the 1st order corrections: 
// reset order and repeat the calculation 
                param.flag[param.order] = 1;
                llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_interf_elastic,
                        USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
                        CP.MINEVAL, CP.MAXEVAL_LO, CP.NSTART, CP.NINCREASE, CP.NBATCH,
                        CP.GRIDNO, CP.STATEFILE, CP.SPIN,
                        &CP.neval, &CP.fail, integral, error, prob);
                output.sigma_pol_elastic_1st = param.P*integral[CP.comp]*nb;
                errors.sigma_pol_elastic_1st = param.P*error[CP.comp]*nb;
                param.flag[param.order] = 2;
// done recalculation of 1st order elastic cross section 
			}
			if (param.flag[param.brems] == 0) {
				output.sigma_pol_1st = output.sigma_pol_elastic_1st;
				errors.sigma_pol_1st = errors.sigma_pol_elastic_1st;
				output.sigma_pol_2nd = output.sigma_pol_elastic_2nd;
				errors.sigma_pol_2nd = errors.sigma_pol_elastic_2nd;
				output.sigma_1st = output.sigma_pol_1st + output.sigma_unpol_1st;
				errors.sigma_1st = sqrt(pow(errors.sigma_pol_1st,2.) + 
                                        pow(errors.sigma_unpol_1st,2.));
				output.sigma_2nd = output.sigma_pol_2nd + output.sigma_unpol_2nd;
				errors.sigma_2nd = sqrt(pow(errors.sigma_pol_2nd,2.) + 
                                        pow(errors.sigma_unpol_2nd,2.));
                output.asymm_1st = output.sigma_pol_1st / output.sigma_unpol_1st;
                errors.asymm_1st = abs(output.sigma_pol_1st/output.sigma_unpol_1st) *
                        sqrt(pow(errors.sigma_unpol_1st/output.sigma_unpol_1st,2.)
                                + pow(errors.sigma_pol_1st/output.sigma_pol_1st,2.));
                output.rel_asymm_1st = - 100.*(output.asymm_1st - output.asymm_born) / 
                                              output.asymm_born;
                output.asymm_2nd = output.sigma_pol_2nd / output.sigma_unpol_2nd;
                errors.asymm_2nd = abs(output.sigma_pol_2nd/output.sigma_unpol_2nd) *
                        sqrt(pow(errors.sigma_unpol_2nd/output.sigma_unpol_2nd,2.)
                                + pow(errors.sigma_pol_2nd/output.sigma_pol_2nd,2.));
                output.rel_asymm_2nd = - 100.*(output.asymm_2nd - output.asymm_born) / 
                                              output.asymm_born;
			}
            std::cout <<"completed\n";
        }
	}
// End section for leptonic non-radiative cross section


// Radiative contributions: One hard photon radiation 
// At first order: just one-photon bremsstrahlung 
// at second order: including loop and soft-photon correction to one-photon bremsstrahlung
	if (param.flag[param.brems] == 1 || param.flag[param.brems] == 2) {
		std::cout << "Integration for one hard-photon unpolarized cross section ... ";
		if (param.flag[param.int_method] == 0) {
			if (param.flag[param.PS] == 0) {
				llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e7, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}
			if (param.flag[param.PS] == 1) {
				llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_ps2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e7, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_ps2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}
		}
		if (param.flag[param.int_method] == 1) {
			if (param.flag[param.PS] == 0) {
				llSuave(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}

			if (param.flag[param.PS] == 1) {
				llSuave(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_ps2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}
		}
		if (param.flag[param.int_method] == 2) {
			if (param.flag[param.PS] == 0) {
				llCuhre(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
						CP.MINEVAL, CP.MAXEVAL_1st, 9,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}
			if (param.flag[param.PS] == 1) {
				llCuhre(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_ps2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
						CP.MINEVAL, CP.MAXEVAL_1st, 9,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}
		}
		if (param.flag[param.order] == 0) {
			output.sigma_unpol_inelastic_1st = integral[CP.comp]*nb;
			errors.sigma_unpol_inelastic_1st = error[CP.comp]*nb;
		}
		if (param.flag[param.order] == 1 || param.flag[param.order] == 2) {
			output.sigma_unpol_inelastic_1st = integral[CP.comp]*nb;
			errors.sigma_unpol_inelastic_1st = error[CP.comp]*nb;
			output.sigma_unpol_1st = output.sigma_unpol_inelastic_1st + 
                                     output.sigma_unpol_elastic_1st;
			errors.sigma_unpol_1st = sqrt(pow(errors.sigma_unpol_inelastic_1st,2.) + 
                                          pow(errors.sigma_unpol_elastic_1st,2.));
		}
		if (param.flag[param.order] == 2) {
			output.sigma_unpol_inelastic_loop = integral[CP.comp]*nb;
			errors.sigma_unpol_inelastic_loop = error[CP.comp]*nb;
		}
		std::cout <<"completed\n";
	}
// End section for one-photon bremsstrahlung


// Hadronic bremsstrahlung corrections (1st order only)
	if (param.flag[param.brems_hadr] == 1 || param.flag[param.brems_hadr] == 3) {
		std::cout << 
        "Integration for hadronic interference hard-photon unpolarized cross section ... ";
		if (param.flag[param.int_method] == 0) {
			llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_hadr_interf,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e7, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_hadr_interf, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_hadr_interf,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_hadr_interf, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
		}
		if (param.flag[param.int_method] == 1) {
			llSuave(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_hadr_interf,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);

		}
		output.sigma_unpol_inelastic_1st_hadr_interf = integral[CP.comp]*nb;
		errors.sigma_unpol_inelastic_1st_hadr_interf = error[CP.comp]*nb;
		if (param.flag[param.order] == 0 || param.flag[param.order] == 1) {
			output.sigma_unpol_inelastic_1st += output.sigma_unpol_inelastic_1st_hadr_interf;
			errors.sigma_unpol_inelastic_1st = sqrt(pow(errors.sigma_unpol_inelastic_1st,2.) + 
                                            pow(errors.sigma_unpol_inelastic_1st_hadr_interf,2.));
			output.sigma_unpol_1st = output.sigma_unpol_inelastic_1st + 
                                     output.sigma_unpol_elastic_1st;
			errors.sigma_unpol_1st = sqrt(pow(errors.sigma_unpol_inelastic_1st,2.) + 
                                          pow(errors.sigma_unpol_elastic_1st,2.));
		}
		if (param.flag[param.order] == 2) {
			output.sigma_unpol_inelastic_loop += output.sigma_unpol_inelastic_1st_hadr_interf;
			errors.sigma_unpol_inelastic_loop = sqrt(pow(errors.sigma_unpol_inelastic_1st,2.) + 
                                        pow(errors.sigma_unpol_inelastic_1st_hadr_interf,2.));
		}
		std::cout <<"completed\n";
	}

	if (param.flag[param.brems_hadr] == 2 || param.flag[param.brems_hadr] == 3) {
		std::cout << "Integration for purely hadronic one hard-photon unpolarized cross section ... ";
		if (param.flag[param.int_method] == 0) {
			llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_hadr,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e7, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_hadr, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_hadr,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_hadr, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
		}
		if (param.flag[param.int_method] == 1) {
			llSuave(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_hadr,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
		}
		output.sigma_unpol_inelastic_1st_hadr = integral[CP.comp]*nb;
		errors.sigma_unpol_inelastic_1st_hadr = error[CP.comp]*nb;
		if (param.flag[param.order] == 0 || param.flag[param.order] == 1) {
			output.sigma_unpol_inelastic_1st += output.sigma_unpol_inelastic_1st_hadr;
			errors.sigma_unpol_inelastic_1st = sqrt(pow(errors.sigma_unpol_inelastic_1st,2.) + 
                                                pow(errors.sigma_unpol_inelastic_1st_hadr,2.));
			output.sigma_unpol_1st = output.sigma_unpol_inelastic_1st + 
                                     output.sigma_unpol_elastic_1st;
			errors.sigma_unpol_1st = sqrt(pow(errors.sigma_unpol_inelastic_1st,2.) + 
                                          pow(errors.sigma_unpol_elastic_1st,2.));
		}
		if (param.flag[param.order] == 2) {
			output.sigma_unpol_inelastic_loop += output.sigma_unpol_inelastic_1st_hadr;
			errors.sigma_unpol_inelastic_loop = sqrt(pow(errors.sigma_unpol_inelastic_1st,2.) + 
                                                pow(errors.sigma_unpol_inelastic_1st_hadr,2.));
		}
		std::cout <<"completed\n";
	}
// End section for hadronic corrections 


// One hard photon radiation for polarized cross section 
// for second order including loop+soft-photon correction
            if (param.flag[param.asymmetry] == 1 && 
                (param.flag[param.brems] == 1 || 
                 param.flag[param.brems] == 2)) {
                std::cout << "Integration for one hard photon polarized cross section ... ";
                if (param.flag[param.int_method] == 0) {
                    llVegas(CP.NDIM_brems_1st, CP.NCOMP, 
                            integrands.cuba_integrand_interf_brems_1st,
                            USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
                            CP.MINEVAL, 1e7, CP.NSTART, CP.NINCREASE, CP.NBATCH,
                            CP.GRIDNO_brems_interf, CP.STATEFILE, CP.SPIN,
                            &CP.neval, &CP.fail, integral, error, prob);
                    llVegas(CP.NDIM_brems_1st, CP.NCOMP, 
                            integrands.cuba_integrand_interf_brems_1st,
                            USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
                            CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
                            CP.GRIDNO_brems_interf, CP.STATEFILE, CP.SPIN,
                            &CP.neval, &CP.fail, integral, error, prob);
                }
                if (param.flag[param.int_method] == 1) {
                    llSuave(CP.NDIM_brems_1st, CP.NCOMP, 
                            integrands.cuba_integrand_interf_brems_1st,
                            USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
                            CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
                            CP.STATEFILE, CP.SPIN,
                            &CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
                }
                if (param.flag[param.int_method] == 2) {
                    llCuhre(CP.NDIM_brems_1st, CP.NCOMP, 
                            integrands.cuba_integrand_interf_brems_1st,
                            USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
                            CP.MINEVAL, CP.MAXEVAL_1st, 9,
                            CP.STATEFILE, CP.SPIN,
                            &CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
                }
                if (param.flag[param.order] == 1 || param.flag[param.order] == 2) {
                    output.sigma_pol_inelastic_1st = param.P*integral[CP.comp]*nb;
                    errors.sigma_pol_inelastic_1st = param.P*error[CP.comp]*nb;
                    output.sigma_pol_1st = output.sigma_pol_inelastic_1st + 
                                           output.sigma_pol_elastic_1st;
                    errors.sigma_pol_1st = errors.sigma_pol_inelastic_1st + 
                                           errors.sigma_pol_elastic_1st;
                    output.asymm_1st = output.sigma_pol_1st / output.sigma_unpol_1st;
                    errors.asymm_1st = abs(output.asymm_1st) * 
                                       sqrt(pow(errors.sigma_pol_1st / output.sigma_pol_1st, 2.) + 
                                            pow(errors.sigma_unpol_1st / output.sigma_unpol_1st, 2.));
                    output.rel_asymm_1st = - 100.*(output.asymm_1st - output.asymm_born) / 
                                                  output.asymm_born;
                    output.sigma_1st = output.sigma_pol_1st + output.sigma_unpol_1st;
                    errors.sigma_1st = sqrt(pow(output.sigma_pol_1st,2.) + 
                                            pow(output.sigma_unpol_1st,2.));
                }
                if (param.flag[param.order] == 2) {
                    output.sigma_pol_inelastic_loop = param.P*integral[CP.comp]*nb;
                    errors.sigma_pol_inelastic_loop = param.P*error[CP.comp]*nb;
                }
                std::cout << "completed\n";
                if (param.flag[param.order] == 2 && param.flag[param.brems_add] == 1) {
                    std::cout << "Integration for one hard + one soft photon finite contribution ... ";
                    llSuave(CP.NDIM_brems_2nd, CP.NCOMP, 
                            integrands.cuba_integrand_brems_2nd_pol_add,
                            USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
                            CP.MINEVAL, CP.MAXEVAL_2nd_add, CP.NNEW, CP.NMIN, CP.FLATNESS,
                            CP.STATEFILE, CP.SPIN,
                            &CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
                    output.sigma_pol_2nd_add = param.P*integral[CP.comp]*nb;
                    errors.sigma_pol_2nd_add = param.P*error[CP.comp]*nb;
                    output.sigma_pol_inelastic_loop += output.sigma_pol_2nd_add;
                    errors.sigma_pol_inelastic_loop = sqrt(pow(errors.sigma_pol_2nd_add,2.) + 
                                                           pow(errors.sigma_pol_inelastic_loop,2.));
                    std::cout <<"completed\n";
                }
            }
// End section for one-photon bremsstrahlung for polarized cross section


// Second order corrections
        if (param.flag[param.brems] >= 2 || param.flag[param.brems_add] >= 1) {
            std::cout <<"Start integration for 2nd order bremsstrahlung\n";
            std::cout <<"***** This may take a while . . . \n";
        }
        
// Second order correction: finite soft*hard terms
		if (param.flag[param.order] == 2 && param.flag[param.brems_add] == 1) {
			std::cout << "Integration for one hard + one soft photon finite contribution ... ";
			double sigma_add_l1k1, sigma_add_l1k2, sigma_add_l2k1, sigma_add_l2k2;
			double error_sigma_add_l1k1, error_sigma_add_l1k2, error_sigma_add_l2k1, 
                   error_sigma_add_l2k2;
// Vegas integration disabled. Always use Suave, independent of input for 
// integration method. 
//			if (param.flag[param.int_method] == 0) {
//				llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_l1k1,
//						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
//						CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
//						CP.GRIDNO_brems_2nd_add_l1k1, CP.STATEFILE, CP.SPIN,
//						&CP.neval, &CP.fail, integral, error, prob);
//
//				llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_l1k1,
//						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
//						CP.MINEVAL, CP.MAXEVAL_2nd_add, CP.NSTART, CP.NINCREASE, CP.NBATCH,
//						CP.GRIDNO_brems_2nd_add_l1k1, CP.STATEFILE, CP.SPIN,
//						&CP.neval, &CP.fail, integral, error, prob);
//
//				sigma_add_l1k1 = integral[CP.comp]*nb;
//				error_sigma_add_l1k1 = error[CP.comp]*nb;
//
//				llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_l2k1,
//						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
//						CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
//						CP.GRIDNO_brems_2nd_add_l2k1, CP.STATEFILE, CP.SPIN,
//						&CP.neval, &CP.fail, integral, error, prob);
//
//				llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_l2k1,
//						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
//						CP.MINEVAL, CP.MAXEVAL_2nd_add, CP.NSTART, CP.NINCREASE, CP.NBATCH,
//						CP.GRIDNO_brems_2nd_add_l2k1, CP.STATEFILE, CP.SPIN,
//						&CP.neval, &CP.fail, integral, error, prob);
//
//				sigma_add_l2k1 = integral[CP.comp]*nb;
//				error_sigma_add_l2k1 = error[CP.comp]*nb;
//
//				llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_l1k2,
//						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
//						CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
//						CP.GRIDNO_brems_2nd_add_l1k2, CP.STATEFILE, CP.SPIN,
//						&CP.neval, &CP.fail, integral, error, prob);
//
//				llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_l1k2,
//						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
//						CP.MINEVAL, CP.MAXEVAL_2nd_add, CP.NSTART, CP.NINCREASE, CP.NBATCH,
//						CP.GRIDNO_brems_2nd_add_l1k2, CP.STATEFILE, CP.SPIN,
//						&CP.neval, &CP.fail, integral, error, prob);
//
//				sigma_add_l1k2 = integral[CP.comp]*nb;
//				error_sigma_add_l1k2 = error[CP.comp]*nb;
//
//				llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_l2k2,
//						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
//						CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
//						CP.GRIDNO_brems_2nd_add_l2k2, CP.STATEFILE, CP.SPIN,
//						&CP.neval, &CP.fail, integral, error, prob);
//
//				llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_l2k2,
//						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
//						CP.MINEVAL, CP.MAXEVAL_2nd_add, CP.NSTART, CP.NINCREASE, CP.NBATCH,
//						CP.GRIDNO_brems_2nd_add_l2k2, CP.STATEFILE, CP.SPIN,
//						&CP.neval, &CP.fail, integral, error, prob);
//
//				sigma_add_l2k2 = integral[CP.comp]*nb;
//				error_sigma_add_l2k2 = error[CP.comp]*nb;
//			}
//			if (param.flag[param.int_method] == 1) {
				llSuave(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_l1k1,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_2nd_add, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
				sigma_add_l1k1 = integral[CP.comp]*nb;
				error_sigma_add_l1k1 = error[CP.comp]*nb;
				llSuave(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_l2k1,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_2nd_add, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
				sigma_add_l2k1 = integral[CP.comp]*nb;
				error_sigma_add_l2k1 = error[CP.comp]*nb;
				llSuave(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_l1k2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_2nd_add, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
				sigma_add_l1k2 = integral[CP.comp]*nb;
				error_sigma_add_l1k2 = error[CP.comp]*nb;
				llSuave(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_l2k2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_2nd_add, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
				sigma_add_l2k2 = integral[CP.comp]*nb;
				error_sigma_add_l2k2 = error[CP.comp]*nb;
//			}
			double sigma_add = sigma_add_l1k1 + sigma_add_l2k1 + sigma_add_l1k2 + sigma_add_l2k2;
			double error_sigma_add = sqrt(pow(error_sigma_add_l1k1,2.) + 
                                          pow(error_sigma_add_l2k1,2.) + 
                                          pow(error_sigma_add_l1k2,2.) + 
                                          pow(error_sigma_add_l2k2,2.));
			output.sigma_unpol_2nd_add = sigma_add;
			errors.sigma_unpol_2nd_add = error_sigma_add;
			output.sigma_unpol_inelastic_loop += output.sigma_unpol_2nd_add;
			errors.sigma_unpol_inelastic_loop = sqrt(pow(errors.sigma_unpol_2nd_add,2.)
													+ pow(errors.sigma_unpol_inelastic_loop,2.));
			std::cout <<"completed\n";
		}
// End section for 'additional' finite terms


//        
    if (param.flag[param.order] == 2 && param.flag[param.brems] == 1) {
        std::cout << "***** WARNING: ***** \n"
                  << "***** 2nd order loop is combined with 1st order bremsstrahlung \n" 
                  << "***** Result is incomplete! \n\n"; 
        output.sigma_unpol_2nd = output.sigma_unpol_elastic_2nd + 
                                 output.sigma_unpol_inelastic_loop;
        errors.sigma_unpol_2nd = sqrt(pow(errors.sigma_unpol_elastic_2nd,2.) + 
                                      pow(errors.sigma_unpol_inelastic_loop,2.));
		}
//


// 2nd order bremsstrahlung: two hard photons
	if (param.flag[param.brems] == 2 || param.flag[param.brems] == 3) {
		std::cout << "Integration for two-hard-photon unpolarized cross section ... ";
		if (param.flag[param.int_method] == 0) {
			llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l1k1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l1k1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l1k1 = integral[CP.comp]*nb;
			sigma_hp_2nd_l1k1_error = error[CP.comp]*nb;
			llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l1k2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l1k2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l1k2 = integral[CP.comp]*nb;
			sigma_hp_2nd_l1k2_error = error[CP.comp]*nb;
			llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l2k1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l2k1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l2k1 = integral[CP.comp]*nb;
			sigma_hp_2nd_l2k1_error = error[CP.comp]*nb;
			llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l2k2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l2k2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l2k2 = integral[CP.comp]*nb;
			sigma_hp_2nd_l2k2_error = error[CP.comp]*nb;
		}
		if (param.flag[param.int_method] == 1) {
			llSuave(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l1k1 = integral[CP.comp]*nb;
			sigma_hp_2nd_l1k1_error = error[CP.comp]*nb;
			llSuave(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l1k2 = integral[CP.comp]*nb;
			sigma_hp_2nd_l1k2_error = error[CP.comp]*nb;
			llSuave(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l2k1 = integral[CP.comp]*nb;
			sigma_hp_2nd_l2k1_error = error[CP.comp]*nb;
			llSuave(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l2k2 = integral[CP.comp]*nb;
			sigma_hp_2nd_l2k2_error = error[CP.comp]*nb;
		}
		if (param.flag[param.int_method] == 2) {
			llCuhre(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.KEY,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l1k1 = integral[CP.comp]*nb;
			sigma_hp_2nd_l1k1_error = error[CP.comp]*nb;
			llCuhre(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.KEY,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l1k2 = integral[CP.comp]*nb;
			sigma_hp_2nd_l1k2_error = error[CP.comp]*nb;
			llCuhre(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.KEY,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l2k1 = integral[CP.comp]*nb;
			sigma_hp_2nd_l2k1_error = error[CP.comp]*nb;
			llCuhre(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.KEY,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l2k2 = integral[CP.comp]*nb;
			sigma_hp_2nd_l2k2_error = error[CP.comp]*nb;
		}
		output.sigma_unpol_inelastic_2nd = sigma_hp_2nd_l1k1 + 
                                           sigma_hp_2nd_l1k2 + 
                                           sigma_hp_2nd_l2k1 + 
                                           sigma_hp_2nd_l2k2;
		errors.sigma_unpol_inelastic_2nd = sqrt(pow(sigma_hp_2nd_l1k1_error,2.) + 
                                                pow(sigma_hp_2nd_l1k2_error,2.) + 
                                                pow(sigma_hp_2nd_l2k1_error,2.) + 
                                                pow(sigma_hp_2nd_l2k2_error,2.));
		if (param.flag[param.order] == 2) {
			output.sigma_unpol_2nd = output.sigma_unpol_elastic_2nd + 
                                     output.sigma_unpol_inelastic_loop + 
                                     output.sigma_unpol_inelastic_2nd;
			errors.sigma_unpol_2nd = sqrt(pow(errors.sigma_unpol_elastic_2nd,2.) + 
                                          pow(errors.sigma_unpol_inelastic_loop,2.) + 
                                          pow(errors.sigma_unpol_inelastic_2nd,2.));
		}
		std::cout <<"completed\n";
	}
// End section for 2-hard-photon bremsstrahlung 


// Hard radiation for differential cross sections
// brems=4: 1st order 2-fold differential in theta_l and theta_gamma 
	if (param.flag[param.brems] == 4) {
		if (param.flag[param.int_method] == 0) {
			if (param.flag[param.PS] == 0) {
				llVegas(CP.NDIM_brems_1st_2diff, CP.NCOMP, 
                        integrands.cuba_integrand_brems_1st_2diff_v1,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e7, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				llVegas(CP.NDIM_brems_1st_2diff, CP.NCOMP, 
                        integrands.cuba_integrand_brems_1st_2diff_v1,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}
			if (param.flag[param.PS] == 1) {
				llVegas(CP.NDIM_brems_1st_2diff, CP.NCOMP, 
                        integrands.cuba_integrand_brems_1st_ps2_2diff,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e7, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				llVegas(CP.NDIM_brems_1st_2diff, CP.NCOMP, 
                        integrands.cuba_integrand_brems_1st_ps2_2diff,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}
		}
		if (param.flag[param.int_method] == 1) {
			if (param.flag[param.PS] == 0) {
				llSuave(CP.NDIM_brems_1st_2diff, CP.NCOMP, 
                        integrands.cuba_integrand_brems_1st_2diff_v1,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}
			if (param.flag[param.PS] == 1) {
				llSuave(CP.NDIM_brems_1st_2diff, CP.NCOMP, 
                        integrands.cuba_integrand_brems_1st_ps2_2diff,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}
		}
		if (param.flag[param.int_method] == 2) {
			if (param.flag[param.PS] == 0) {
				llCuhre(CP.NDIM_brems_1st_2diff, CP.NCOMP, 
                        integrands.cuba_integrand_brems_1st_2diff_v1,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.KEY,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}
			if (param.flag[param.PS] == 1) {
				llCuhre(CP.NDIM_brems_1st_2diff, CP.NCOMP, 
                        integrands.cuba_integrand_brems_1st_ps2_2diff,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.KEY,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}
		}
		if (param.flag[param.order] == 1) {
			output.sigma_unpol_inelastic_1st = integral[CP.comp]*nb;
			errors.sigma_unpol_inelastic_1st = error[CP.comp]*nb;
		}
		if (param.flag[param.order] == 2) {
			output.sigma_unpol_inelastic_loop = integral[CP.comp]*nb;
			errors.sigma_unpol_inelastic_loop = error[CP.comp]*nb;
		}
	}
// End section bremsstrahlung type = 4

// brems=5: 2nd order 2-fold differential in theta_l and theta_gamma 
	if (param.flag[param.brems] == 5) {
		if (param.flag[param.int_method] == 0) {
			llVegas(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l1k1_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l1k1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			llVegas(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l1k1_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l1k1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l1k1 = integral[CP.comp]*nb;
			sigma_hp_2nd_l1k1_error = error[CP.comp]*nb;
			llVegas(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l1k2_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l1k2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			llVegas(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l1k2_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l1k2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l1k2 = integral[CP.comp]*nb;
			sigma_hp_2nd_l1k2_error = error[CP.comp]*nb;
			llVegas(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l2k1_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l2k1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			llVegas(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l2k1_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l2k1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l2k1 = integral[CP.comp]*nb;
			sigma_hp_2nd_l2k1_error = error[CP.comp]*nb;
			llVegas(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l2k2_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l2k2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			llVegas(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l2k2_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l2k2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l2k2 = integral[CP.comp]*nb;
			sigma_hp_2nd_l2k2_error = error[CP.comp]*nb;
		}
		if (param.flag[param.int_method] == 1) {
			llSuave(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l1k1_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l1k1 = integral[CP.comp]*nb;
			sigma_hp_2nd_l1k1_error = error[CP.comp]*nb;
			llSuave(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l1k2_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l1k2 = integral[CP.comp]*nb;
			sigma_hp_2nd_l1k2_error = error[CP.comp]*nb;
			llSuave(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l2k1_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l2k1 = integral[CP.comp]*nb;
			sigma_hp_2nd_l2k1_error = error[CP.comp]*nb;
			llSuave(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l2k2_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l2k2 = integral[CP.comp]*nb;
			sigma_hp_2nd_l2k2_error = error[CP.comp]*nb;
		}
		if (param.flag[param.int_method] == 2) {
			llCuhre(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l1k1_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.KEY,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l1k1 = integral[CP.comp]*nb;
			sigma_hp_2nd_l1k1_error = error[CP.comp]*nb;
			llCuhre(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l1k2_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.KEY,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l1k2 = integral[CP.comp]*nb;
			sigma_hp_2nd_l1k2_error = error[CP.comp]*nb;
			llCuhre(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l2k1_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.KEY,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l2k1 = integral[CP.comp]*nb;
			sigma_hp_2nd_l2k1_error = error[CP.comp]*nb;
			llCuhre(CP.NDIM_brems_2nd_2diff, CP.NCOMP, 
                    integrands.cuba_integrand_brems_2nd_l2k2_2diff,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.KEY,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_hp_2nd_l2k2 = integral[CP.comp]*nb;
			sigma_hp_2nd_l2k2_error = error[CP.comp]*nb;
		}
		output.sigma_unpol_inelastic_2nd = sigma_hp_2nd_l1k1 + 
                                           sigma_hp_2nd_l1k2 + 
                                           sigma_hp_2nd_l2k1 + 
                                           sigma_hp_2nd_l2k2;
		errors.sigma_unpol_inelastic_2nd = sqrt(pow(sigma_hp_2nd_l1k1_error,2.) + 
                                                pow(sigma_hp_2nd_l1k2_error,2.) + 
                                                pow(sigma_hp_2nd_l2k1_error,2.) + 
                                                pow(sigma_hp_2nd_l2k2_error,2.));

	}
// End section bremsstrahlung type = 5

// brems=6: 2nd order 1-fold differential in theta'_gamma 
	if (param.flag[param.brems] == 6) {
		llVegas(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_1diff,
				USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
				CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
				CP.GRIDNO_brems_2nd, CP.STATEFILE, CP.SPIN,
				&CP.neval, &CP.fail, integral, error, prob);
		llVegas(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_1diff,
				USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
				CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
				CP.GRIDNO_brems_2nd, CP.STATEFILE, CP.SPIN,
				&CP.neval, &CP.fail, integral, error, prob);
		output.sigma_unpol_inelastic_2nd = integral[CP.comp]*nb;
		errors.sigma_unpol_inelastic_2nd = error[CP.comp]*nb;
	}
// End section bremsstrahlung type = 6

// brems=7: 1st order 2-fold differential in theta_l and E' 
	if (param.flag[param.brems] == 7) {
		llVegas(CP.NDIM_brems_1st_2diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_2diff_v2,
				USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
				CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
				CP.GRIDNO, CP.STATEFILE, CP.SPIN,
				&CP.neval, &CP.fail, integral, error, prob);
		output.sigma_unpol_inelastic_1st = integral[CP.comp]*nb;
		errors.sigma_unpol_inelastic_1st = error[CP.comp]*nb;
	}
// End section bremsstrahlung type = 7

// brems=8: 1st order 3-fold differential in theta_l, theta_gamma and phi_gamma 
	if (param.flag[param.brems] == 8) {
		llSuave(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_brems_1st_ps2_3diff,
				USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
				CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
				CP.STATEFILE, CP.SPIN,
				&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
		output.sigma_unpol_inelastic_1st = integral[CP.comp]*nb;
		errors.sigma_unpol_inelastic_1st = error[CP.comp]*nb;
	}
// End section bremsstrahlung type = 8


// Two-photon bremsstrahlung for the polarized cross section 
	if (param.flag[param.asymmetry] == 1 && 
        (param.flag[param.brems] == 2 || 
         param.flag[param.brems] == 3)) {
		std::cout << "Integration for two hard photon polarized cross section ... ";
		double sigma_pol_2hg_1, error_sigma_pol_2hg_1, sigma_pol_2hg_2, error_sigma_pol_2hg_2;
		if (param.flag[param.int_method] == 0) {
			llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_interf_2nd_1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_interf_2nd_1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			sigma_pol_2hg_1 = integral[CP.comp]*nb;
			error_sigma_pol_2hg_1 = error[CP.comp]*nb;
			llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_interf_2nd_2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_interf_2nd_2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			sigma_pol_2hg_2 = integral[CP.comp]*nb;
			error_sigma_pol_2hg_2 = error[CP.comp]*nb;
		}
		if (param.flag[param.int_method] == 1) {
			llSuave(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_pol_2hg_1 = integral[CP.comp]*nb;
			error_sigma_pol_2hg_1 = error[CP.comp]*nb;
			llSuave(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_pol_2hg_2 = integral[CP.comp]*nb;
			error_sigma_pol_2hg_2 = error[CP.comp]*nb;
		}
		if (param.flag[param.int_method] == 2) {
			llCuhre(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.KEY,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_pol_2hg_1 = integral[CP.comp]*nb;
			error_sigma_pol_2hg_1 = error[CP.comp]*nb;
			llCuhre(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.KEY,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			sigma_pol_2hg_2 = integral[CP.comp]*nb;
			error_sigma_pol_2hg_2 = error[CP.comp]*nb;
		}
		output.sigma_pol_inelastic_2nd = param.P*(sigma_pol_2hg_1 + sigma_pol_2hg_2);
		errors.sigma_pol_inelastic_2nd = param.P*sqrt(pow(error_sigma_pol_2hg_1,2.) + 
                                                      pow(error_sigma_pol_2hg_2,2.));
		if (param.flag[param.order] == 2) {
			output.sigma_pol_2nd = output.sigma_pol_elastic_2nd + 
                                   output.sigma_pol_inelastic_loop + 
                                   output.sigma_pol_inelastic_2nd;
			errors.sigma_pol_2nd = sqrt(pow(errors.sigma_pol_elastic_2nd,2.) +
                                        pow(errors.sigma_pol_inelastic_loop,2.) + 
                                        pow(errors.sigma_pol_inelastic_2nd,2.));
			output.asymm_2nd = output.sigma_pol_2nd / output.sigma_unpol_2nd;
			errors.asymm_2nd = abs(output.asymm_2nd) * 
                               sqrt(pow(errors.sigma_pol_2nd / output.sigma_pol_2nd, 2.) + 
                                    pow(errors.sigma_unpol_2nd / output.sigma_unpol_2nd, 2.));
			output.rel_asymm_2nd = - 100.*(output.asymm_2nd - output.asymm_born) / 
                                         output.asymm_born;
			errors.sigma_2nd = sqrt(pow(output.sigma_pol_2nd,2.) + 
                                    pow(output.sigma_unpol_2nd,2.));
			output.sigma_2nd = output.sigma_pol_2nd + output.sigma_unpol_2nd;
			errors.sigma_2nd = sqrt(pow(output.sigma_pol_2nd,2.) + 
                                    pow(output.sigma_unpol_2nd,2.));
		}
		std::cout << "completed\n";
	}
// End section for two-photon bremsstrahlung for polarized cross section 


// To keep track of 1st-order results if 2nd order was requested: 
// Recalculate one-photon bremsstrahlung without loop+soft-photon
	if (param.flag[param.brems] == 2 && param.flag[param.order] == 2) {
        param.flag[param.order] = 1;
		param.MAXEVAL_1st = param.maxeval_1st_aux;
		llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_elastic,
				USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
				CP.MINEVAL, 1e7, CP.NSTART, CP.NINCREASE, CP.NBATCH,
				CP.GRIDNO, CP.STATEFILE, CP.SPIN,
				&CP.neval, &CP.fail, integral, error, prob);
        output.sigma_unpol_elastic_1st = integral[CP.comp]*nb;
        errors.sigma_unpol_elastic_1st = error[CP.comp]*nb;
		if (param.flag[param.int_method] == 0) {
			if (param.flag[param.PS] == 0) {
				llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e7, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}
			if (param.flag[param.PS] == 1) {
				llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_ps2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e7, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_ps2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}
		}
        if (param.flag[param.int_method] == 1) {
			if (param.flag[param.PS] == 0) {
				llSuave(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}
			if (param.flag[param.PS] == 1) {
				llSuave(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_ps2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}
        }
		output.sigma_unpol_inelastic_1st = integral[CP.comp]*nb;
		errors.sigma_unpol_inelastic_1st = error[CP.comp]*nb;
        output.sigma_unpol_1st = output.sigma_unpol_inelastic_1st + output.sigma_unpol_elastic_1st;
		errors.sigma_unpol_1st = errors.sigma_unpol_inelastic_1st + errors.sigma_unpol_elastic_1st;
// add hadronic contribution if non-zero
        output.sigma_unpol_inelastic_1st += output.sigma_unpol_inelastic_1st_hadr_interf;
        errors.sigma_unpol_inelastic_1st = sqrt(pow(errors.sigma_unpol_inelastic_1st,2.) + 
                                            pow(errors.sigma_unpol_inelastic_1st_hadr_interf,2.));
        output.sigma_unpol_1st += output.sigma_unpol_inelastic_1st_hadr_interf;
        errors.sigma_unpol_1st = sqrt(pow(errors.sigma_unpol_1st,2.) + 
                                      pow(errors.sigma_unpol_inelastic_1st_hadr_interf,2.));
        output.sigma_unpol_inelastic_1st += output.sigma_unpol_inelastic_1st_hadr;
        errors.sigma_unpol_inelastic_1st = sqrt(pow(errors.sigma_unpol_inelastic_1st,2.) + 
                                            pow(errors.sigma_unpol_inelastic_1st_hadr,2.));
        output.sigma_unpol_1st += output.sigma_unpol_inelastic_1st_hadr;
        errors.sigma_unpol_1st = sqrt(pow(errors.sigma_unpol_1st,2.) + 
                                      pow(errors.sigma_unpol_inelastic_1st_hadr,2.));

        if (param.flag[param.asymmetry] == 1) {
			llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_interf_elastic,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			output.sigma_pol_elastic_1st = param.P*integral[CP.comp]*nb;
            errors.sigma_pol_elastic_1st = param.P*error[CP.comp]*nb;
			llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_interf_brems_1st,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_interf, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			output.sigma_pol_inelastic_1st = param.P*integral[CP.comp]*nb;
			errors.sigma_pol_inelastic_1st = param.P*error[CP.comp]*nb;
			output.sigma_pol_1st = output.sigma_pol_inelastic_1st + output.sigma_pol_elastic_1st;
			errors.sigma_pol_1st = errors.sigma_pol_inelastic_1st + errors.sigma_pol_elastic_1st;
			output.asymm_1st = output.sigma_pol_1st / output.sigma_unpol_1st;
			errors.asymm_1st = abs(output.asymm_1st) * 
                               sqrt(pow(errors.sigma_pol_1st / output.sigma_pol_1st, 2.) + 
                                    pow(errors.sigma_unpol_1st / output.sigma_unpol_1st, 2.));
			output.sigma_1st = output.sigma_pol_1st + output.sigma_unpol_1st;
			errors.sigma_1st = sqrt(pow(output.sigma_pol_1st,2.) + 
                                    pow(output.sigma_unpol_1st,2.));
		}
		param.flag[param.order] = 2;
		param.MAXEVAL_1st = param.MAXEVAL_gamma_loop;
	}
// End section to repeat one-photon bremsstrahlung at 1st order 

        
// param.en <= 10 ?
	if (param.en <= 10.) {
		if (param.flag[param.brems] == 0) {
			output.sigma_unpol_1st_vect.push_back(output.sigma_unpol_1st);
			output.ev_brems_1st.push_back(0.);
			output.sigma_unpol_1st_elastic_vect.push_back(output.sigma_unpol_elastic_1st);
			if (param.flag[param.asymmetry] == 1) {
				output.sigma_pol_1st_vect.push_back(output.sigma_pol_1st);
				output.sigma_1st_vect.push_back(output.sigma_1st);
			}
		}
		if (param.flag[param.brems] == 1) {
			output.ev_brems_1st.push_back(output.sigma_unpol_inelastic_1st / 
                                          output.sigma_unpol_1st);
			output.sigma_unpol_1st_vect.push_back(output.sigma_unpol_1st);
			output.sigma_unpol_1st_elastic_vect.push_back(output.sigma_unpol_elastic_1st);
			output.sigma_unpol_1st_inelastic_vect.push_back(output.sigma_unpol_inelastic_1st);
			if (param.flag[param.asymmetry] == 1) {
				output.sigma_pol_1st_vect.push_back(output.sigma_pol_1st);
				output.sigma_1st_vect.push_back(output.sigma_1st);
			}
		}
		if (param.flag[param.brems] == 2) {
			output.ev_brems_1st.push_back(output.sigma_unpol_inelastic_1st / 
                                          output.sigma_unpol_2nd);
			output.sigma_unpol_1st_vect.push_back(output.sigma_unpol_1st);
			output.ev_brems_2nd.push_back(output.sigma_unpol_inelastic_2nd / 
                                          output.sigma_unpol_2nd);
			output.ev_brems_2nd_l1k1.push_back(sigma_hp_2nd_l1k1 / 
                                               output.sigma_unpol_inelastic_2nd);
			output.ev_brems_2nd_l1k2.push_back(sigma_hp_2nd_l1k2 / 
                                               output.sigma_unpol_inelastic_2nd);
			output.ev_brems_2nd_l2k1.push_back(sigma_hp_2nd_l2k1 / 
                                               output.sigma_unpol_inelastic_2nd);
			output.ev_brems_2nd_l2k2.push_back(sigma_hp_2nd_l2k2 / 
                                               output.sigma_unpol_inelastic_2nd);
			output.sigma_unpol_2nd_vect.push_back(output.sigma_unpol_2nd);
			if (param.flag[param.asymmetry] == 1) {
				output.sigma_pol_1st_vect.push_back(output.sigma_pol_1st);
				output.sigma_1st_vect.push_back(output.sigma_1st);
				output.sigma_pol_2nd_vect.push_back(output.sigma_pol_2nd);
				output.sigma_2nd_vect.push_back(output.sigma_2nd);
			}
		}
		if (param.flag[param.brems] == 3) {
			output.ev_brems_2nd.push_back(output.sigma_unpol_inelastic_2nd / 
                                          output.sigma_unpol_2nd);
			output.ev_brems_2nd_l1k1.push_back(sigma_hp_2nd_l1k1 / 
                                               output.sigma_unpol_inelastic_2nd);
			output.ev_brems_2nd_l1k2.push_back(sigma_hp_2nd_l1k2 / 
                                               output.sigma_unpol_inelastic_2nd);
			output.ev_brems_2nd_l2k1.push_back(sigma_hp_2nd_l2k1 / 
                                               output.sigma_unpol_inelastic_2nd);
			output.ev_brems_2nd_l2k2.push_back(sigma_hp_2nd_l2k2 / 
                                               output.sigma_unpol_inelastic_2nd);
			output.sigma_unpol_2nd_vect.push_back(output.sigma_unpol_2nd);
			if (param.flag[param.asymmetry] == 1) {
				output.sigma_pol_2nd_vect.push_back(output.sigma_pol_2nd);
				output.sigma_2nd_vect.push_back(output.sigma_2nd);
			}
		}
	}
// End section (param.en <= 10.)
	}


// Carbon target 
	if (param.flag[param.target] == 1) {

		if (param.flag[param.LO] == 1) {

			std::cout << "Numerical integration for the LO unpolarized cross-section ... ";

			llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_born_carbon,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_LO, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);

			std::cout <<"completed\n";

			output.sigma_unpol_born = integral[CP.comp]*nb;
			errors.sigma_unpol_born = error[CP.comp]*nb;
		}

		if (param.flag[param.order] >= 0) {

			std::cout << "Numerical integration for the soft-photon unpolarized cross-section ... ";

			llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_elastic_carbon,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_LO, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_elastic, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);

				errors.sigma_unpol_elastic_1st = error[CP.comp]*nb;
				output.sigma_unpol_elastic_1st = integral[CP.comp]*nb;

			if (param.flag[param.brems] == 0) {
				output.sigma_unpol_1st = output.sigma_unpol_elastic_1st;
				errors.sigma_unpol_1st = errors.sigma_unpol_elastic_1st;
			}
			std::cout <<"completed\n";
		}

		if (param.flag[param.asymmetry] == 1) {
			if (param.flag[param.LO] == 1) {

				std::cout << "Numerical integration for the LO polarized cross-section ... ";

				llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_interf_born_carbon,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_LO, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);


				output.sigma_pol_born = param.P*integral[CP.comp]*nb;
				errors.sigma_pol_born = param.P*error[CP.comp]*nb;
				output.asymm_born = output.sigma_pol_born / output.sigma_unpol_born;
				errors.asymm_born = abs(output.sigma_pol_born/output.sigma_unpol_born) *
						sqrt(pow(errors.sigma_unpol_born/output.sigma_unpol_born,2.)
								+ pow(errors.sigma_pol_born/output.sigma_pol_born,2.));
				output.sigma_born = output.sigma_pol_born + output.sigma_unpol_born;
				errors.sigma_born = sqrt(pow(output.sigma_pol_born,2.)
													+ pow(output.sigma_unpol_born,2.));

				std::cout <<"completed\n";
			}

			if (param.flag[param.order] >= 0) {

				std::cout << "Numerical integration for the soft-photon polarized cross-section ... ";

				llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_interf_elastic_carbon,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_LO, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);

					errors.sigma_pol_elastic_1st = param.P*error[CP.comp]*nb;
					output.sigma_pol_elastic_1st = param.P*integral[CP.comp]*nb;

				if (param.flag[param.brems] == 0) {
					output.sigma_pol_1st = output.sigma_pol_elastic_1st;
					errors.sigma_pol_1st = errors.sigma_pol_elastic_1st;
				}
			}

			std::cout <<"completed\n";

		}

		if (param.flag[param.brems] == 1) {

			std::cout << "Numerical integration for one hard-photon unpolarized cross-section ... ";

			if (param.flag[param.int_method] == 0) {

					llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_carbon,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
							CP.MINEVAL, 1e7, CP.NSTART, CP.NINCREASE, CP.NBATCH,
							CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);

					llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_carbon,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
							CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
							CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
			}

			if (param.flag[param.int_method] == 1) {

				llSuave(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_carbon,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}

			output.sigma_unpol_inelastic_1st = integral[CP.comp]*nb;
			errors.sigma_unpol_inelastic_1st = error[CP.comp]*nb;
			output.sigma_unpol_1st = output.sigma_unpol_inelastic_1st + output.sigma_unpol_elastic_1st;
			errors.sigma_unpol_1st = errors.sigma_unpol_inelastic_1st + errors.sigma_unpol_elastic_1st;

		}

		std::cout <<"completed\n";

		if (param.flag[param.asymmetry] == 1 && param.flag[param.brems] == 1) {

			std::cout << "Numerical integration for one hard-photon polarized cross-section ... ";

			if (param.flag[param.int_method] == 0) {

				llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_interf_brems_1st_carbon,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e7, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_interf, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);

				llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_interf_brems_1st_carbon,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_interf, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}

			if (param.flag[param.int_method] == 1) {

				llSuave(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_interf_brems_1st_carbon,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}

				output.sigma_pol_inelastic_1st = param.P*integral[CP.comp]*nb;
				errors.sigma_pol_inelastic_1st = param.P*error[CP.comp]*nb;
				output.sigma_pol_1st = output.sigma_pol_inelastic_1st + output.sigma_pol_elastic_1st;
				errors.sigma_pol_1st = errors.sigma_pol_inelastic_1st + errors.sigma_pol_elastic_1st;
				output.asymm_1st = output.sigma_pol_1st / output.sigma_unpol_1st;
				errors.asymm_1st = abs(output.asymm_1st) * sqrt(pow(errors.sigma_pol_1st / output.sigma_pol_1st, 2.)
						+ pow(errors.sigma_unpol_1st / output.sigma_unpol_1st, 2.));
				output.sigma_1st = output.sigma_pol_1st + output.sigma_unpol_1st;
				errors.sigma_1st = sqrt(pow(errors.sigma_pol_1st,2.) + pow(errors.sigma_unpol_1st,2.));

			std::cout << "completed\n";
		}

		if (param.en <= 10.) {
			if (param.flag[param.brems] == 0) {
				output.sigma_unpol_1st_vect.push_back(output.sigma_unpol_1st);
				output.ev_brems_1st.push_back(0.);
				output.sigma_unpol_1st_elastic_vect.push_back(output.sigma_unpol_elastic_1st);
				if (param.flag[param.asymmetry] == 1) {
					output.sigma_pol_1st_vect.push_back(output.sigma_pol_1st);
					output.sigma_1st_vect.push_back(output.sigma_1st);
				}
			}
			if (param.flag[param.brems] == 1) {
				output.ev_brems_1st.push_back(output.sigma_unpol_inelastic_1st / output.sigma_unpol_1st);
				output.sigma_unpol_1st_vect.push_back(output.sigma_unpol_1st);
				output.sigma_unpol_1st_elastic_vect.push_back(output.sigma_unpol_elastic_1st);
				output.sigma_unpol_1st_inelastic_vect.push_back(output.sigma_unpol_inelastic_1st);
				if (param.flag[param.asymmetry] == 1) {
					output.sigma_pol_1st_vect.push_back(output.sigma_pol_1st);
					output.sigma_1st_vect.push_back(output.sigma_1st);
				}
			}
		}
	}


// Electron target
	if (param.flag[param.target] == 2) {

			if (param.flag[param.LO] == 1) {

				std::cout << "Numerical integration for the LO unpolarized cross-section ... ";

				llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_born_e,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_LO, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);

				std::cout <<"completed\n";

				output.sigma_unpol_born = integral[CP.comp]*nb;
				errors.sigma_unpol_born = error[CP.comp]*nb;
			}


			if (param.flag[param.brems] == 1) {

				std::cout << "Numerical integration for one hard-photon unpolarized cross-section ... ";

				if (param.flag[param.int_method] == 0) {

						llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_e,
								USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
								CP.MINEVAL, 1e7, CP.NSTART, CP.NINCREASE, CP.NBATCH,
								CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
								&CP.neval, &CP.fail, integral, error, prob);

						llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_e,
								USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems, CP.SEED,
								CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
								CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
								&CP.neval, &CP.fail, integral, error, prob);
				}

				if (param.flag[param.int_method] == 1) {

					llSuave(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_e,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
							CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
							CP.STATEFILE, CP.SPIN,
							&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
				}

				output.sigma_unpol_inelastic_1st = integral[CP.comp]*nb;
				errors.sigma_unpol_inelastic_1st = error[CP.comp]*nb;
//				output.sigma_unpol_1st = output.sigma_unpol_inelastic_1st + output.sigma_unpol_elastic_1st;
//				errors.sigma_unpol_1st = errors.sigma_unpol_inelastic_1st + errors.sigma_unpol_elastic_1st;

			}

			std::cout <<"completed\n";
	}

	write_output();

	std::cout <<"\n";

	return 0;

}


int PES::events(){

	param.final_param(input);
	CP.set_param(&param);
	integrands.set_param(&param, &interpolation, &FS);

	void *USERDATA = &integrands;

	cubacores(1, 10000);

	param.flag[param.ps] = 1;
	param.aux = FS.E_prime_l;

	double integral[CP.NCOMP], error[CP.NCOMP], prob[CP.NCOMP];

	double r2 = rand.uniform();
//	double r1 = rand.uniform();
	double r0 = rand.uniform();

	if (param.flag[param.brems] == 0) {
		while (FS.E_prime_l == param.aux) {
				llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_elastic,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
						0, 1, 1, 0, 1,
						CP.GRIDNO_elastic, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				// FS.seed += 1e7;
				FS.event_type = 0;
				FS.avg_weight = FS.weight*nb / output.sigma_unpol_1st_elastic_vect[int(param.en/param.Delta_E) -
						                                          int(param.min[param.E]/param.Delta_E)];
			if (FS.E_prime_l == param.aux) FS.failed_ev++;
		}
	}

	if (param.flag[param.brems] == 1 && param.flag[param.order] == 0) {
		while (FS.E_prime_l == param.aux) {
			if (param.flag[param.PS] == 0) {
				llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
						0, 1, 1, 0, 1,
						CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}
			if (param.flag[param.PS] == 1) {
				llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_ps2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
						0, 1, 1, 0, 1,
						CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}
			// FS.seed += 1e7;
			FS.avg_weight = FS.weight*nb / output.sigma_unpol_1st_inelastic_vect[int(param.en/param.Delta_E) -
					                                          int(param.min[param.E]/param.Delta_E)];
			if (FS.E_prime_l == param.aux) FS.failed_ev++;
		}
	}

	if (param.flag[param.brems] == 1 && param.flag[param.order] != 0) {
		while (FS.E_prime_l == param.aux) {
			if (r0 >= output.ev_brems_1st[int(param.en/param.Delta_E) -
			                                          int(param.min[param.E]/param.Delta_E)]) {
				llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_elastic,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
						0, 1, 1, 0, 1,
						CP.GRIDNO_elastic, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				// FS.seed += 1e7;
				FS.event_type = 0;
				FS.avg_weight = FS.weight*nb / output.sigma_unpol_1st_elastic_vect[int(param.en/param.Delta_E) -
						                                          int(param.min[param.E]/param.Delta_E)];
			}
			else {
				if (param.flag[param.PS] == 0) {
					llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
				}
				if (param.flag[param.PS] == 1) {
					llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_ps2,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
				}
				// FS.seed += 1e7;
				FS.event_type = 1;
				FS.avg_weight = FS.weight*nb / output.sigma_unpol_1st_inelastic_vect[int(param.en/param.Delta_E) -
						                                          int(param.min[param.E]/param.Delta_E)];
			}
			if (FS.E_prime_l == param.aux) FS.failed_ev++;
		}
	}


	if (param.flag[param.brems] == 2 && param.flag[param.order] != 0) {
		while (FS.E_prime_l == param.aux) {
			if (r0 >= output.ev_brems_1st[int(param.en/param.Delta_E) -
			                                          int(param.min[param.E]/param.Delta_E)]
			          + output.ev_brems_2nd[int(param.en/param.Delta_E) -
                                            int(param.min[param.E]/param.Delta_E)]) {
				llVegas(CP.NDIM_ELASTIC, CP.NCOMP, integrands.cuba_integrand_elastic,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
						0, 1, 1, 0, 1,
						CP.GRIDNO_elastic, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				// FS.seed += 1e7;
				FS.event_type = 0;
			}
			else if (r0 >= output.ev_brems_2nd[int(param.en/param.Delta_E) -
			                                               int(param.min[param.E]/param.Delta_E)]) {
				if (param.flag[param.PS] == 0) {
					llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
				}
				if (param.flag[param.PS] == 1) {
					llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_ps2,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
				}
				// FS.seed += 1e7;
				FS.event_type = 1;
			}
			else {
				if (r2 <= output.ev_brems_2nd_l1k1[int(param.en/param.Delta_E) -
				                                               int(param.min[param.E]/param.Delta_E)]){
					llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k1,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_2nd_l1k1, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
					// FS.seed += 1e7;
					FS.event_type = 2;
				}
				else if (r2 <= output.ev_brems_2nd_l1k1[int(param.en/param.Delta_E) -
				                                                    int(param.min[param.E]/param.Delta_E)]
				        + output.ev_brems_2nd_l1k2[int(param.en/param.Delta_E) -
				                                                    int(param.min[param.E]/param.Delta_E)]){
					llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k2,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_2nd_l1k2, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
					// FS.seed += 1e7;
					FS.event_type = 2;
				}
				else if (r2 <= output.ev_brems_2nd_l1k1[int(param.en/param.Delta_E) -
				                                               int(param.min[param.E]/param.Delta_E)]
				        + output.ev_brems_2nd_l1k2[int(param.en/param.Delta_E) -
				                                                    int(param.min[param.E]/param.Delta_E)]
				        + output.ev_brems_2nd_l2k1[int(param.en/param.Delta_E) -
				                                                    int(param.min[param.E]/param.Delta_E)]){
					llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k1,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_2nd_l2k1, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
					// FS.seed += 1e7;
					FS.event_type = 2;

				}
				else {
					llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k2,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_2nd_l2k2, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
					// FS.seed += 1e7;
					FS.event_type = 2;
				}
			}
			if (FS.E_prime_l == param.aux) FS.failed_ev++;
		}
	}


	if (param.flag[param.brems] == 2 && param.flag[param.order] == 0) {
		while (FS.E_prime_l == param.aux) {
			if (r0 >= output.ev_brems_2nd[int(param.en/param.Delta_E) -
			                              int(param.min[param.E]/param.Delta_E)]) {
				if (param.flag[param.PS] == 0) {
					llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
				}
				if (param.flag[param.PS] == 1) {
					llVegas(CP.NDIM_brems_1st, CP.NCOMP, integrands.cuba_integrand_brems_1st_ps2,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_1st, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
				}
				// FS.seed += 1e7;
				FS.event_type = 1;
			}
			else if (r2 <= output.ev_brems_2nd_l1k1[int(param.en/param.Delta_E) -
			                                   int(param.min[param.E]/param.Delta_E)]){
				llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k1,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
						0, 1, 1, 0, 1,
						CP.GRIDNO_brems_2nd_l1k1, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				// FS.seed += 1e7;
				FS.event_type = 3;
			}
			else if (r2 <= output.ev_brems_2nd_l1k1[int(param.en/param.Delta_E) -
			                                        int(param.min[param.E]/param.Delta_E)]
			                                        + output.ev_brems_2nd_l1k2[int(param.en/param.Delta_E) -
			                                                                   int(param.min[param.E]/param.Delta_E)]){
				llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
						0, 1, 1, 0, 1,
						CP.GRIDNO_brems_2nd_l1k2, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				// FS.seed += 1e7;
				FS.event_type = 4;
			}
			else if (r2 <= output.ev_brems_2nd_l1k1[int(param.en/param.Delta_E) -
			                                        int(param.min[param.E]/param.Delta_E)]
			                                        + output.ev_brems_2nd_l1k2[int(param.en/param.Delta_E) -
			                                                                   int(param.min[param.E]/param.Delta_E)]
			                                                                   + output.ev_brems_2nd_l2k1[int(param.en/param.Delta_E) -
			                                                                                              int(param.min[param.E]/param.Delta_E)]){
				llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k1,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
						0, 1, 1, 0, 1,
						CP.GRIDNO_brems_2nd_l2k1, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				// FS.seed += 1e7;
				FS.event_type = 5;

			}
			else {
				llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
						0, 1, 1, 0, 1,
						CP.GRIDNO_brems_2nd_l2k2, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				// FS.seed += 1e7;
				FS.event_type = 6;
			}
			if (FS.E_prime_l == param.aux) FS.failed_ev++;
		}
	}

	if (param.flag[param.brems] == 3 && param.flag[param.order] == 0) {
		while (FS.E_prime_l == param.aux) {
				if (r2 <= output.ev_brems_2nd_l1k1[int(param.en/param.Delta_E) -
				                                               int(param.min[param.E]/param.Delta_E)]){
					llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k1,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_2nd_l1k1, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
					// FS.seed += 1e7;
					FS.event_type = 3;
				}
				else if (r2 <= output.ev_brems_2nd_l1k1[int(param.en/param.Delta_E) -
				                                                    int(param.min[param.E]/param.Delta_E)]
				        + output.ev_brems_2nd_l1k2[int(param.en/param.Delta_E) -
				                                                    int(param.min[param.E]/param.Delta_E)]){
					llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l1k2,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_2nd_l1k2, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
					// FS.seed += 1e7;
					FS.event_type = 4;
				}
				else if (r2 <= output.ev_brems_2nd_l1k1[int(param.en/param.Delta_E) -
				                                               int(param.min[param.E]/param.Delta_E)]
				        + output.ev_brems_2nd_l1k2[int(param.en/param.Delta_E) -
				                                                    int(param.min[param.E]/param.Delta_E)]
				        + output.ev_brems_2nd_l2k1[int(param.en/param.Delta_E) -
				                                                    int(param.min[param.E]/param.Delta_E)]){
					llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k1,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_2nd_l2k1, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
					// FS.seed += 1e7;
					FS.event_type = 5;

				}
				else {
					llVegas(CP.NDIM_brems_2nd, CP.NCOMP, integrands.cuba_integrand_brems_2nd_l2k2,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, FS.seed,
							0, 1, 1, 0, 1,
							CP.GRIDNO_brems_2nd_l2k2, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);
					// FS.seed += 1e7;
					FS.event_type = 6;
				}
			if (FS.E_prime_l == param.aux) FS.failed_ev++;
		}
	}

	set_final_state();
	if (FS.seed <= 1e7) FS.seed += 1e7;

	return 0;

}


int PES::shiftQ2(const double thl_deg) {

	double  average_deltaQ2, sigma_thl, inelastic_thl,
	elastic_thl, error_average, error_inelastic_thl;

	param.final_param(input);
	CP.set_param(&param);
	integrands.set_param(&param, &interpolation);
	param.set_thl(thl_deg);

	void *USERDATA = &integrands;

	cubacores(CP.no_cores, 10000);

	double integral[CP.NCOMP], error[CP.NCOMP], prob[CP.NCOMP];

	int grid_thl = CP.GRIDNO_brems + int(thl_deg);

	if (param.flag[param.int_method] == 0) {
		llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_shiftQ2,
				USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
				CP.MINEVAL, 1e5, CP.NSTART, CP.NINCREASE, CP.NBATCH,
				grid_thl, CP.STATEFILE, CP.SPIN,
				&CP.neval, &CP.fail, integral, error, prob);

		llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_shiftQ2,
				USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
				CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
				grid_thl, CP.STATEFILE, CP.SPIN,
				&CP.neval, &CP.fail, integral, error, prob);
	}

	if (param.flag[param.int_method] == 1) {
		llSuave(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_shiftQ2,
				USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
				CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
				CP.STATEFILE, CP.SPIN,
				&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);

			}

		average_deltaQ2 = integral[CP.comp];
		error_average = error[CP.comp];

		if (param.flag[param.int_method] == 0) {
			if (param.flag[param.PS] == 0) {
				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e5, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);

				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}

			if (param.flag[param.PS] == 1) {
				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl_ps2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e5, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);

				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl_ps2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}
		}

		if (param.flag[param.int_method] == 1) {
					if (param.flag[param.PS] == 0) {
						llSuave(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl,
								USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
								CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
								CP.STATEFILE, CP.SPIN,
								&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
					}

					if (param.flag[param.PS] == 1) {
						llSuave(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl_ps2,
								USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
								CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
								CP.STATEFILE, CP.SPIN,
								&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
					}
				}

		inelastic_thl = integral[CP.comp];
		error_inelastic_thl = error[CP.comp];

		elastic_thl = integrands.CS.crsect_elastic_thl(param.Q2);

		sigma_thl = inelastic_thl + elastic_thl;

		output.shiftQ2 = average_deltaQ2 / sigma_thl;
		errors.shiftQ2 = output.shiftQ2 * sqrt(pow(error_average/average_deltaQ2,2.)
				+ pow(error_inelastic_thl/inelastic_thl,2.));

		output.rel_shiftQ2 = -100. * output.shiftQ2 / param.Q2;
		errors.rel_shiftQ2 = -100. * errors.shiftQ2 / param.Q2;

		return 0;
}


int PES::sigma_diff_Omega_l(const double thl_deg) {

	param.final_param(input);
	CP.set_param(&param);
	integrands.set_param(&param, &interpolation);
	param.set_thl(thl_deg);
	output.Q2 = param.Q2;

	void *USERDATA = &integrands;

	cubacores(CP.no_cores, 10000);

	double integral[CP.NCOMP], error[CP.NCOMP], prob[CP.NCOMP];

	if (param.flag[param.target] == 0) {

		output.sigma_unpol_born =  integrands.CS.crsect_born_thl(param.Q2) * nb / (2.*pi);

		if (param.flag[param.order] == 1 || param.flag[param.hadr_corr] != 0) {
			output.sigma_unpol_elastic_1st = integrands.CS.crsect_elastic_thl(param.Q2);
		}

		if (param.flag[param.order] == 2) {
			output.sigma_unpol_elastic_2nd = integrands.CS.crsect_elastic_thl(param.Q2);
		}

	if (param.flag[param.brems] == 1 || param.flag[param.brems] == 2) {

		if (param.flag[param.int_method] == 0) {
			if (param.flag[param.PS] == 0) {
				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e5, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);

				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}

			if (param.flag[param.PS] == 1) {
				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl_ps2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e5, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);

				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl_ps2,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}
		}

		if (param.flag[param.int_method] == 1) {
					if (param.flag[param.PS] == 0) {
						llSuave(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl,
								USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
								CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
								CP.STATEFILE, CP.SPIN,
								&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
					}

					if (param.flag[param.PS] == 1) {
						llSuave(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl_ps2,
								USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
								CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
								CP.STATEFILE, CP.SPIN,
								&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
					}
				}

		if (param.flag[param.order] == 0 || param.flag[param.order] == 1) {
			output.sigma_unpol_inelastic_1st = integral[CP.comp];
			errors.sigma_unpol_inelastic_1st = error[CP.comp];
		}

		if (param.flag[param.order] == 2) {
			output.sigma_unpol_inelastic_loop = integral[CP.comp];
			errors.sigma_unpol_inelastic_loop = error[CP.comp];
		}
	}

	if (param.flag[param.brems_hadr] == 1 || param.flag[param.brems_hadr] == 3) {

		if (param.flag[param.int_method] == 0) {

				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_hadr_interf_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e5, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_hadr_interf, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);

				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_hadr_interf_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_hadr_interf, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
		}

		if (param.flag[param.int_method] == 1) {
						llSuave(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_hadr_interf_thl,
								USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
								CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
								CP.STATEFILE, CP.SPIN,
								&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
		}

		output.sigma_unpol_inelastic_1st_hadr_interf = integral[CP.comp];
		errors.sigma_unpol_inelastic_1st_hadr_interf = error[CP.comp];

		if (param.flag[param.order] == 0 || param.flag[param.order] == 1) {
			output.sigma_unpol_inelastic_1st += output.sigma_unpol_inelastic_1st_hadr_interf;
			errors.sigma_unpol_inelastic_1st = sqrt(pow(errors.sigma_unpol_inelastic_1st,2.)
										+ pow(errors.sigma_unpol_inelastic_1st_hadr_interf,2.));
		}

		if (param.flag[param.order] == 2) {
			output.sigma_unpol_inelastic_loop += output.sigma_unpol_inelastic_1st_hadr_interf;
			errors.sigma_unpol_inelastic_loop = sqrt(pow(errors.sigma_unpol_inelastic_1st,2.)
												+ pow(errors.sigma_unpol_inelastic_1st_hadr_interf,2.));
		}
	}

	if (param.flag[param.brems_hadr] == 2 || param.flag[param.brems_hadr] == 3) {

		if (param.flag[param.int_method] == 0) {

				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_hadr_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e5, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_hadr, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);

				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_hadr_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_hadr, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
		}

		if (param.flag[param.int_method] == 1) {
						llSuave(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_hadr_thl,
								USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
								CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
								CP.STATEFILE, CP.SPIN,
								&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
		}

		output.sigma_unpol_inelastic_1st_hadr = integral[CP.comp];
		errors.sigma_unpol_inelastic_1st_hadr = error[CP.comp];

		if (param.flag[param.order] == 0 || param.flag[param.order] == 1) {
			output.sigma_unpol_inelastic_1st += output.sigma_unpol_inelastic_1st_hadr;
			errors.sigma_unpol_inelastic_1st = sqrt(pow(errors.sigma_unpol_inelastic_1st,2.)
										+ pow(errors.sigma_unpol_inelastic_1st_hadr,2.));
		}

		if (param.flag[param.order] == 2) {
			output.sigma_unpol_inelastic_loop += output.sigma_unpol_inelastic_1st_hadr;
			errors.sigma_unpol_inelastic_loop = sqrt(pow(errors.sigma_unpol_inelastic_1st,2.)
												+ pow(errors.sigma_unpol_inelastic_1st_hadr,2.));
		}
	}

	if (param.flag[param.order] == 0 || param.flag[param.order] == 1) {
		output.sigma_unpol_1st = output.sigma_unpol_inelastic_1st + output.sigma_unpol_elastic_1st;
		errors.sigma_unpol_1st = errors.sigma_unpol_inelastic_1st;
	}


	if (param.flag[param.order] == 2 && param.flag[param.brems_add] == 1) {

		double sigma_add_1, sigma_add_2, error_sigma_add_1, error_sigma_add_2;

//		llSuave(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_thl,
//				USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
//				CP.MINEVAL, CP.MAXEVAL_2nd_add, CP.NNEW, CP.NMIN, CP.FLATNESS,
//				CP.STATEFILE, CP.SPIN,
//				&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
//
//		output.sigma_unpol_2nd_add = integral[CP.comp];
//		errors.sigma_unpol_2nd_add = error[CP.comp];

		llSuave(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_thl_1,
				USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
				CP.MINEVAL, CP.MAXEVAL_2nd_add, CP.NNEW, CP.NMIN, CP.FLATNESS,
				CP.STATEFILE, CP.SPIN,
				&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);

		sigma_add_1 = integral[CP.comp];
		error_sigma_add_1 = error[CP.comp];

		llSuave(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_add_thl_2,
				USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
				CP.MINEVAL, CP.MAXEVAL_2nd_add, CP.NNEW, CP.NMIN, CP.FLATNESS,
				CP.STATEFILE, CP.SPIN,
				&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);

		sigma_add_2 = integral[CP.comp];
		error_sigma_add_2 = error[CP.comp];

		output.sigma_unpol_2nd_add = sigma_add_1 + sigma_add_2;
		errors.sigma_unpol_2nd_add = sqrt(pow(error_sigma_add_1,2.)+pow(error_sigma_add_2,2.));

		output.sigma_unpol_inelastic_loop += output.sigma_unpol_2nd_add;
		errors.sigma_unpol_inelastic_loop = sqrt(pow(errors.sigma_unpol_2nd_add,2.) +
										pow(errors.sigma_unpol_inelastic_loop,2.));
	}

	if (param.flag[param.order] == 2 && param.flag[param.brems] == 1) {
		output.sigma_unpol_2nd = output.sigma_unpol_inelastic_loop + output.sigma_unpol_elastic_2nd;
		errors.sigma_unpol_2nd = errors.sigma_unpol_inelastic_loop;
	}

	if (param.flag[param.brems] == 2 || param.flag[param.brems] == 3) {

		double sigma_2hp_1, error_sigma_2hp_1, sigma_2hp_2, error_sigma_2hp_2;

		if (param.flag[param.int_method] == 0) {
//			llVegas(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_thl,
//					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
//					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
//					CP.GRIDNO_brems_2nd, CP.STATEFILE, CP.SPIN,
//					&CP.neval, &CP.fail, integral, error, prob);
//
//			llVegas(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_thl,
//					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
//					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
//					CP.GRIDNO_brems_2nd, CP.STATEFILE, CP.SPIN,
//					&CP.neval, &CP.fail, integral, error, prob);

			llVegas(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_thl_1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l1k1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);

			llVegas(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_thl_1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l1k1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);

			sigma_2hp_1 = integral[CP.comp];
			error_sigma_2hp_1 = error[CP.comp];

			llVegas(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_thl_2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l1k2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);

			llVegas(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_thl_2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_2nd_l1k2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);

			sigma_2hp_2 = integral[CP.comp];
			error_sigma_2hp_2 = error[CP.comp];

		}

		if (param.flag[param.int_method] == 1) {
//			llSuave(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_thl,
//					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
//					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NNEW, CP.NMIN, CP.FLATNESS,
//					CP.STATEFILE, CP.SPIN,
//					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			llSuave(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_thl_1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);

			sigma_2hp_1 = integral[CP.comp];
			error_sigma_2hp_1 = error[CP.comp];

			llSuave(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_thl_2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);

			sigma_2hp_2 = integral[CP.comp];
			error_sigma_2hp_2 = error[CP.comp];

			}

//		output.sigma_unpol_inelastic_2nd = integral[CP.comp];
//		errors.sigma_unpol_inelastic_2nd = error[CP.comp];
		output.sigma_unpol_inelastic_2nd = sigma_2hp_1 + sigma_2hp_2;
		errors.sigma_unpol_inelastic_2nd = sqrt(pow(error_sigma_2hp_1,2.) + pow(error_sigma_2hp_2,2.));
		output.sigma_unpol_2nd = output.sigma_unpol_inelastic_2nd +
				output.sigma_unpol_inelastic_loop + output.sigma_unpol_elastic_2nd;
		errors.sigma_unpol_2nd = sqrt(pow(errors.sigma_unpol_inelastic_loop,2.) +
										pow(errors.sigma_unpol_inelastic_2nd,2.));
	}

		if (param.flag[param.asymmetry] == 1) {
		output.sigma_pol_born = param.P*integrands.CS.interf_born_thl(param.Q2) * nb / (2.*pi);
		output.asymm_born = output.sigma_pol_born/output.sigma_unpol_born;
		output.sigma_born = output.sigma_pol_born + output.sigma_unpol_born;

		if (param.flag[param.brems] == 1 || param.flag[param.brems] == 2) {

			if (param.flag[param.int_method] == 0) {
			llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_1st_thl,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e5, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_interf, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);

			llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_1st_thl,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_interf, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);
			}
			if (param.flag[param.int_method] == 1) {
				llSuave(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_1st_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}

			if (param.flag[param.order] == 1) {
				output.sigma_pol_inelastic_1st = param.P*integral[CP.comp];
				errors.sigma_pol_inelastic_1st = param.P*error[CP.comp];
		}

			if (param.flag[param.order] == 2) {
				output.sigma_pol_inelastic_loop = param.P*integral[CP.comp];
				errors.sigma_pol_inelastic_loop = param.P*error[CP.comp];
			}
		}

		if (param.flag[param.order] == 1) {
			output.sigma_pol_elastic_1st = param.P*integrands.CS.interf_elastic_thl(param.Q2);
			output.sigma_pol_1st = output.sigma_pol_inelastic_1st + output.sigma_pol_elastic_1st;
			errors.sigma_pol_1st = errors.sigma_pol_inelastic_1st;
			output.asymm_1st = output.sigma_pol_1st / output.sigma_unpol_1st;
			output.rel_asymm_1st = - 100.*(output.asymm_1st - output.asymm_born) / output.asymm_born;
			errors.asymm_1st = abs(output.asymm_1st) * sqrt(pow(errors.sigma_pol_1st / output.sigma_pol_1st, 2.)
					+ pow(errors.sigma_unpol_1st / output.sigma_unpol_1st, 2.));
			output.sigma_1st = output.sigma_pol_1st + output.sigma_unpol_1st;
			errors.sigma_1st = sqrt(pow(errors.sigma_pol_1st,2.) + pow(errors.sigma_unpol_1st,2.));
		}

		if (param.flag[param.order] == 2) {
			output.sigma_pol_elastic_2nd = param.P*integrands.CS.interf_elastic_thl(param.Q2);
		}

		if (param.flag[param.order] == 2 && param.flag[param.brems_add] == 1) {
			llSuave(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_brems_2nd_pol_add_thl,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd_add, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);

			output.sigma_pol_2nd_add = param.P*integral[CP.comp];
			errors.sigma_pol_2nd_add = param.P*error[CP.comp];
			output.sigma_pol_inelastic_loop += output.sigma_pol_2nd_add;
			errors.sigma_pol_inelastic_loop = sqrt(pow(errors.sigma_pol_inelastic_loop,2.) +
											pow(errors.sigma_pol_2nd_add,2.));
		}
	}

	if (param.flag[param.asymmetry] == 1 && (param.flag[param.brems] == 2 || param.flag[param.brems] == 3)) {

		double sigma_pol_2hp_1, error_sigma_pol_2hp_1, sigma_pol_2hp_2, error_sigma_pol_2hp_2;

		if (param.flag[param.int_method] == 0) {
			llVegas(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_thl_1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_interf_2nd_1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);

			llVegas(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_thl_1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_interf_2nd_1, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);

			sigma_pol_2hp_1 = integral[CP.comp];
			error_sigma_pol_2hp_1 = error[CP.comp];

			llVegas(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_thl_2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_interf_2nd_2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);

			llVegas(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_thl_2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_2nd, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems_interf_2nd_2, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);

			sigma_pol_2hp_2 = integral[CP.comp];
			error_sigma_pol_2hp_2 = error[CP.comp];
		}

		if (param.flag[param.int_method] == 1) {
			llSuave(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_thl_1,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);

			sigma_pol_2hp_1 = integral[CP.comp];
			error_sigma_pol_2hp_1 = error[CP.comp];

			llSuave(CP.NDIM_brems_2nd_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_2nd_thl_2,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
					CP.STATEFILE, CP.SPIN,
					&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);

			sigma_pol_2hp_2 = integral[CP.comp];
			error_sigma_pol_2hp_2 = error[CP.comp];
		}

		output.sigma_pol_inelastic_2nd = param.P*(sigma_pol_2hp_1 + sigma_pol_2hp_2);
		errors.sigma_pol_inelastic_2nd = param.P*sqrt(pow(error_sigma_pol_2hp_1,2.)
													+ pow(error_sigma_pol_2hp_2,2.));

		output.sigma_pol_2nd = output.sigma_pol_inelastic_2nd +
				output.sigma_pol_inelastic_loop + output.sigma_pol_elastic_2nd;
		errors.sigma_pol_2nd = sqrt(pow(errors.sigma_pol_inelastic_loop,2.)
				+ pow(errors.sigma_pol_inelastic_2nd,2.));

		output.asymm_2nd = output.sigma_pol_2nd / output.sigma_unpol_2nd;
		output.rel_asymm_2nd = -100.*(output.asymm_2nd - output.asymm_born) / output.asymm_born;
		errors.asymm_2nd = abs(output.asymm_2nd) * sqrt(pow(errors.sigma_pol_2nd / output.sigma_pol_2nd, 2.)
				+ pow(errors.sigma_unpol_2nd / output.sigma_unpol_2nd, 2.));
		output.sigma_2nd = output.sigma_pol_2nd + output.sigma_unpol_2nd;
		errors.sigma_2nd = sqrt(pow(errors.sigma_pol_2nd,2.) + pow(errors.sigma_unpol_2nd,2.));

	}

		if (param.flag[param.brems] == 2 && param.flag[param.order] == 2) {

			param.flag[param.order] = 1;

			llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);

			llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_thl,
					USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
					CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
					CP.GRIDNO_brems, CP.STATEFILE, CP.SPIN,
					&CP.neval, &CP.fail, integral, error, prob);

			output.sigma_unpol_elastic_1st = integrands.CS.crsect_elastic_thl(param.Q2);
			output.sigma_unpol_inelastic_1st = integral[CP.comp];
			errors.sigma_unpol_inelastic_1st = error[CP.comp];
			output.sigma_unpol_1st = output.sigma_unpol_inelastic_1st + output.sigma_unpol_elastic_1st;
			errors.sigma_unpol_1st = errors.sigma_unpol_inelastic_1st;

			if (param.flag[param.asymmetry] == 1) {

					llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_1st_thl,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
							CP.MINEVAL, 1e8, CP.NSTART, CP.NINCREASE, CP.NBATCH,
							CP.GRIDNO_brems_interf, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);

					llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_1st_thl,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
							CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
							CP.GRIDNO_brems_interf, CP.STATEFILE, CP.SPIN,
							&CP.neval, &CP.fail, integral, error, prob);

					output.sigma_pol_elastic_1st = param.P*integrands.CS.interf_elastic_thl(param.Q2);
					output.sigma_pol_inelastic_1st = param.P*integral[CP.comp];
					errors.sigma_pol_inelastic_1st = param.P*error[CP.comp];
					output.sigma_pol_1st = output.sigma_pol_inelastic_1st + output.sigma_pol_elastic_1st;
					errors.sigma_pol_1st = errors.sigma_pol_inelastic_1st;
					output.asymm_1st = output.sigma_pol_1st / output.sigma_unpol_1st;
					output.asymm_1st *= param.P;
					output.rel_asymm_1st = - 100.*(output.asymm_1st - output.asymm_born) / output.asymm_born;
					errors.asymm_1st = abs(output.asymm_1st) * sqrt(pow(errors.sigma_pol_1st / output.sigma_pol_1st, 2.)
							+ pow(errors.sigma_unpol_1st / output.sigma_unpol_1st, 2.));
					output.sigma_1st = output.sigma_pol_1st + output.sigma_unpol_1st;
					errors.sigma_1st = sqrt(pow(errors.sigma_pol_1st,2.) + pow(errors.sigma_unpol_1st,2.));
				}

			param.flag[param.order] = 2;
		}
	}

	if (param.flag[param.target] == 1) {
		output.sigma_unpol_born =  integrands.CS.crsect_born_carbon_thl(param.Q2) * nb / (2.*pi);
		output.sigma_unpol_elastic_1st = integrands.CS.crsect_elastic_carbon_thl(param.Q2);

//		std::cout << integrands.CS.asymm_born_carbon_test(param.Q2) <<"\n";
//		std::cout << integrands.CS.asymm_born_carbon_test(param.Q2)*nb/(2.*pi)/output.sigma_unpol_born <<"\n";

		if (param.flag[param.brems] == 1) {

			if (param.flag[param.int_method] == 0) {
				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_carbon_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e5, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);

				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_carbon_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
			}

			if (param.flag[param.int_method] == 1) {
				llSuave(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_brems_1st_carbon_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
						CP.STATEFILE, CP.SPIN,
						&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
			}

			output.sigma_unpol_inelastic_1st = integral[CP.comp];
			errors.sigma_unpol_inelastic_1st = error[CP.comp];

		}
		output.sigma_unpol_1st = output.sigma_unpol_inelastic_1st + output.sigma_unpol_elastic_1st;
		errors.sigma_unpol_1st = errors.sigma_unpol_inelastic_1st;

		if (param.flag[param.asymmetry] == 1) {

			output.sigma_pol_born = param.P*integrands.CS.interf_born_carbon_thl(param.Q2) * nb / (2.*pi);
			output.sigma_pol_elastic_1st = param.P*integrands.CS.interf_elastic_carbon_thl(param.Q2);
			output.sigma_born = output.sigma_pol_born + output.sigma_unpol_born;

			if (param.flag[param.brems] == 1) {

				if (param.flag[param.int_method] == 0) {
				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_1st_carbon_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, 1e5, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_interf, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);

				llVegas(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_1st_carbon_thl,
						USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags, CP.SEED,
						CP.MINEVAL, CP.MAXEVAL_1st, CP.NSTART, CP.NINCREASE, CP.NBATCH,
						CP.GRIDNO_brems_interf, CP.STATEFILE, CP.SPIN,
						&CP.neval, &CP.fail, integral, error, prob);
				}
				if (param.flag[param.int_method] == 1) {
					llSuave(CP.NDIM_brems_1st_1diff, CP.NCOMP, integrands.cuba_integrand_interf_brems_1st_carbon_thl,
							USERDATA, CP.NVEC, CP.EPSREL, CP.EPSABS, CP.flags_brems | 4, CP.SEED,
							CP.MINEVAL, CP.MAXEVAL_1st, CP.NNEW, CP.NMIN, CP.FLATNESS,
							CP.STATEFILE, CP.SPIN,
							&CP.nregions, &CP.neval, &CP.fail, integral, error, prob);
				}
			}
				output.sigma_pol_inelastic_1st = param.P*integral[CP.comp];
				errors.sigma_pol_inelastic_1st = param.P*error[CP.comp];
				output.sigma_pol_1st = output.sigma_pol_inelastic_1st + output.sigma_pol_elastic_1st;
				errors.sigma_pol_1st = errors.sigma_pol_inelastic_1st;
				output.asymm_born = output.sigma_pol_born / output.sigma_unpol_born;
				output.asymm_1st = output.sigma_pol_1st / output.sigma_unpol_1st;
				errors.asymm_1st = abs(output.asymm_1st) * sqrt(pow(errors.sigma_pol_1st / output.sigma_pol_1st, 2.)
						+ pow(errors.sigma_unpol_1st / output.sigma_unpol_1st, 2.));
				output.sigma_1st = output.sigma_pol_1st + output.sigma_unpol_1st;
				errors.sigma_1st = sqrt(pow(errors.sigma_pol_1st,2.) + pow(errors.sigma_unpol_1st,2.));
				output.rel_asymm_1st = - 100.*(output.asymm_1st - output.asymm_born) / output.asymm_born;
		}

	}
		output.sigma_unpol_elastic_1st *= nb / (2.*pi);
		output.sigma_unpol_inelastic_1st *= nb / (2.*pi);
		errors.sigma_unpol_inelastic_1st *= nb / (2.*pi);
		output.sigma_unpol_inelastic_1st_hadr *= nb / (2.*pi);
		errors.sigma_unpol_inelastic_1st_hadr *= nb / (2.*pi);
		output.sigma_unpol_inelastic_1st_hadr_interf *= nb / (2.*pi);
		errors.sigma_unpol_inelastic_1st_hadr_interf *= nb / (2.*pi);
		output.sigma_unpol_inelastic_loop *= nb / (2.*pi);
		errors.sigma_unpol_inelastic_loop *= nb / (2.*pi);
		output.sigma_unpol_1st *= nb / (2.*pi);
		errors.sigma_unpol_1st *= nb / (2.*pi);
		output.sigma_unpol_elastic_2nd *= nb / (2.*pi);
		output.sigma_unpol_inelastic_2nd *= nb / (2.*pi);
		errors.sigma_unpol_inelastic_2nd *= nb / (2.*pi);
		output.sigma_unpol_2nd_add *= nb / (2.*pi);
		errors.sigma_unpol_2nd_add *= nb / (2.*pi);
		output.sigma_unpol_2nd *= nb / (2.*pi);
		errors.sigma_unpol_2nd *= nb / (2.*pi);
		output.sigma_pol_elastic_1st *= nb / (2.*pi);
		output.sigma_pol_inelastic_1st *= nb / (2.*pi);
		errors.sigma_pol_inelastic_1st *= nb / (2.*pi);
		output.sigma_pol_1st *= nb / (2.*pi);
		errors.sigma_pol_1st *= nb / (2.*pi);
		output.sigma_pol_elastic_2nd *= nb / (2.*pi);
		output.sigma_pol_inelastic_loop *= nb / (2.*pi);
		errors.sigma_pol_inelastic_loop *= nb / (2.*pi);
		output.sigma_pol_inelastic_2nd *= nb / (2.*pi);
		errors.sigma_pol_inelastic_2nd *= nb / (2.*pi);
		output.sigma_pol_2nd_add *= nb / (2.*pi);
		errors.sigma_pol_2nd_add *= nb / (2.*pi);
		output.sigma_pol_2nd *= nb / (2.*pi);
		errors.sigma_pol_2nd *= nb / (2.*pi);
		output.sigma_1st *= nb / (2.*pi);
		errors.sigma_1st *= nb / (2.*pi);
		output.sigma_2nd *= nb / (2.*pi);
		errors.sigma_2nd *= nb / (2.*pi);

	write_output();

	return 0;
}


double PES::delta(const double thl_deg) {

	param.final_param(input);
	CP.set_param(&param);
	integrands.set_param(&param, &interpolation);
	param.set_thl(thl_deg);

	double delta, delta_comp;

	FS.Q2 = param.Q2;

	Form_factors ff;

	double gpe, gpm, gpze, gpzm, gae, f1, f2, f1z, f2z;
	double kappa = integrands.CS.VC.kappa_weak(param.Q2);
	double ga = -1./2.;
	double gv = -1./2. + 2.*kappa*param.sw2;

	ff.ffz(param.Q2, gpze, gpzm, kappa);
	ff.ffgae (param.Q2, gae);
	ff.ffactp (param.Q2, gpe, gpm);

// Pauli and Dirac Form Factors
	double tau = param.Q2/(4.*M2);
	f2=(gpm-gpe)/(1.+tau);
	f1=(gpe+tau*gpm)/(1.+tau);
//	f2z=(gpzm-gpze)/(1.+tau);
//	f1z=(gpze+tau*gpzm)/(1.+tau);

	delta = integrands.CS.VC.d_brems_ee(param.Q2) + integrands.CS.VC.d_vert(param.Q2, f1, f2);

	if (param.flag[param.vac_pol] == 1)
		delta += integrands.CS.VC.d_vac_1st(param.Q2,me);
	if (param.flag[param.vac_pol] == 2)
		delta += integrands.CS.VC.d_vac_1st(param.Q2,me) + integrands.CS.VC.d_vac_1st(param.Q2, m_mu)
		+ integrands.CS.VC.d_vac_1st(param.Q2,m_tau);
	if (param.flag[param.vac_pol] == 3)
		delta += interpolation.d_vac_hadr(param.Q2);
	if (param.flag[param.vac_pol] == 4)
		delta += interpolation.d_vac_hadr(param.Q2)
				+ integrands.CS.VC.d_vac_2nd(param.Q2,me) + integrands.CS.VC.d_vac_2nd(param.Q2,m_mu)
				+ integrands.CS.VC.d_vac_2nd(param.Q2,m_tau);
	if (param.flag[param.tpe] == 1) {
		if (param.flag[param.lepton] == 0)
			delta += interpolation.tpe(param.Q2, param.eps);
		if (param.flag[param.lepton] == 1)
			delta -= interpolation.tpe(param.Q2, param.eps);
	}

	if (param.flag[param.brems] == 2) {
		delta += integrands.CS.VC.d_2nd_total(param.Q2, f1, f2);
	}

	return delta;
}


double PES::running_sw2(const double Q2) {

	param.final_param(input);
	CP.set_param(&param);
	integrands.set_param(&param, &interpolation);

	output.asymm_born = integrands.CS.interf_born_thl(Q2) / integrands.CS.crsect_born_thl(Q2);

	return integrands.CS.VC.kappa_weak(Q2)*param.sw2;
}


double PES::get_total_unpol_cross_section(const double E) {

	param.en = round(E/param.Delta_E) * param.Delta_E;

	return output.sigma_unpol_1st_vect[int(param.en/param.Delta_E) - int(param.min[param.E]/param.Delta_E)];
}


void PES::set_final_state() {

	double l2 = sqrt(pow(FS.E_prime_l,2.) - m2);

	//FS.phi_l = rand.uniform() * 2.*pi;
	FS.phi_l = 0.;
	FS.E_p = param.en - FS.E_prime_l - FS.E_gamma - FS.E_gamma_prime + M;
	FS.Q2 = 2.*M*(param.en - FS.E_prime_l - FS.E_gamma - FS.E_gamma_prime);

	double p2 = sqrt(pow(FS.E_p,2.) - M2);

//	if (FS.E_gamma_prime != 0 && FS.phi_gamma < 0.) FS.phi_gamma += 2.*pi;
//	if (FS.E_gamma_prime != 0 && FS.phi_gamma_prime < 0.) FS.phi_gamma_prime += 2.*pi;

//	if (FS.E_gamma != 0 && FS.E_gamma_prime != 0)
//		if (rand.uniform() > 0.5) {
//			FS.phi_gamma_prime = 2.*pi - FS.phi_gamma_prime;
//			FS.phi_gamma = 2.*pi - FS.phi_gamma;
//		}
	if (FS.E_gamma != 0 && FS.E_gamma_prime == 0)
		if (rand.uniform() > 0.5) FS.phi_gamma = 2.*pi - FS.phi_gamma;

	FS.E = param.en;

	FS.l_1[0] = param.en;
	FS.l_1[1] = 0.;
	FS.l_1[2] = 0.;
	FS.l_1[3] = param.l1;

	FS.p_1[0] = M;
	FS.p_1[1] = 0.;
	FS.p_1[2] = 0.;
	FS.p_1[3] = 0.;

	FS.l_2[0] = FS.E_prime_l;
	FS.l_2[1] = l2 * sin(FS.theta_l) * cos(FS.phi_l);
	FS.l_2[2] = l2 * sin(FS.theta_l) * sin(FS.phi_l);
	FS.l_2[3] = l2 * cos(FS.theta_l);

	FS.k_1[0] = FS.E_gamma;
	FS.k_1[1] = FS.E_gamma * sin(FS.theta_gamma) * cos(FS.phi_gamma);
	FS.k_1[2] = FS.E_gamma * sin(FS.theta_gamma) * sin(FS.phi_gamma);
	FS.k_1[3] = FS.E_gamma * cos(FS.theta_gamma);

	FS.k_2[0] = FS.E_gamma_prime;
	FS.k_2[1] = FS.E_gamma_prime * sin(FS.theta_gamma_prime) * cos(FS.phi_gamma_prime);
	FS.k_2[2] = FS.E_gamma_prime * sin(FS.theta_gamma_prime) * sin(FS.phi_gamma_prime);
	FS.k_2[3] = FS.E_gamma_prime * cos(FS.theta_gamma_prime);

	FS.p_2[0] = param.en - FS.l_2[0] - FS.k_1[0] - FS.k_2[0] + M;
	FS.p_2[1] = - FS.l_2[1] - FS.k_1[1] - FS.k_2[1];
	FS.p_2[2] = - FS.l_2[2] - FS.k_1[2] -  FS.k_2[2];
	FS.p_2[3] = param.l1 - FS.l_2[3] - FS.k_1[3] - FS.k_2[3];

	double test = m2 + M*(param.en-FS.E_prime_l-FS.E_gamma-FS.E_gamma_prime) - param.en*(FS.E_prime_l+FS.E_gamma+FS.E_gamma_prime)
				+ param.l1*l2*cos(FS.theta_l) + param.l1*FS.E_gamma*cos(FS.theta_gamma)
				+ param.l1*FS.E_gamma_prime*cos(FS.theta_gamma_prime) + FS.E_prime_l*(FS.E_gamma+FS.E_gamma_prime)
				- l2*FS.E_gamma*sin(FS.theta_l)*sin(FS.theta_gamma)*cos(FS.phi_gamma)
				- l2*FS.E_gamma_prime*sin(FS.theta_l)*sin(FS.theta_gamma_prime)*cos(FS.phi_gamma_prime)
				- l2*FS.E_gamma*cos(FS.theta_l)*cos(FS.theta_gamma)
				- l2*FS.E_gamma_prime*cos(FS.theta_l)*cos(FS.theta_gamma_prime)
				+ FS.E_gamma*FS.E_gamma_prime*(1.
						-sin(FS.theta_gamma)*cos(FS.phi_gamma)*sin(FS.theta_gamma_prime)*cos(FS.phi_gamma_prime)
						-sin(FS.theta_gamma)*sin(FS.phi_gamma)*sin(FS.theta_gamma_prime)*sin(FS.phi_gamma_prime)
						-cos(FS.theta_gamma)*cos(FS.theta_gamma_prime));

	if (std::abs(test) > 1e-10)	std::cout<<"Warning! The momentum conservation is not fulfilled"<<"\n";

	FS.event_no++;
}


void PES::write_output() {

	time_t present_time;
	time(&present_time);
	char* dt = ctime(&present_time);

	if (param.output_file != "0") {
		std::ofstream out;
		param.output_file += ".out";
		out.precision(7);

// Output to file
		out.open(param.output_file.c_str());

		out <<"##########################################################################\n"
            <<"##                               POLARES "<< VERSION << "\n"
            <<"##\n"
            <<"##       Radiative Corrections for Polarized Electron-Proton Scattering   \n"
            <<"##\n"
            <<"##                              R.-D. Bucoveanu                           \n"
            <<"##\n"
            <<"##                         "<<dt
            <<"##\n"
            <<"## If you use POLARES please cite R.-D. Bucoveanu and H. Spiesberger,     \n"
            <<"## Eur. Phys. J. A (2019) 55: 57, arXiv:1811.04970 [hep-ph].              \n"
            <<"## Copyright (c) Razvan Bucoveanu, 2019. E-mail: rabucove@uni-mainz.de    \n"
            <<"##########################################################################\n";
        out <<"##                                  Input                                 \n"
            <<"##\n"
            <<"## [General Input]                                                        \n"
            <<"## Type of incident lepton = ";
		if (param.flag[param.lepton] == 0) out<<"electron\n";
		if (param.flag[param.lepton] == 1) out<<"positron\n";
		if (param.flag[param.lepton] == 2) out<<"muon\n";
		if (param.flag[param.lepton] == 3) out<<"anti-muon\n";
		out	<<"## Type of target particle = ";
		if (param.flag[param.target] == 0) out<<"proton\n";
		if (param.flag[param.target] == 1) out<<"carbon-12\n";
		out <<"## Incident lepton energy = "<<param.en<<" GeV\n";
		out <<"## Lower cut-off value for the photon energy (Delta) = "
				<<param.min[param.E_gamma] * 1000.<<" MeV\n";
		out <<"## Type of cuts for elastic scattering = ";
		if (param.flag[param.cuts_born] == 0) {
			out <<"Scattering angle (theta_l) cuts\n";
			out <<"## theta_l min = "<<param.min[param.theta_l_deg]<<" degrees\n";
			out <<"## theta_l max = "<<param.max[param.theta_l_deg]<<" degrees\n";
		}
		if (param.flag[param.cuts_born] == 1) {
			out <<"Q^2 cuts\n";
			out <<"## Q^2 min = "<<param.min[param.Q2_elastic]<<" GeV^2\n";
			out <<"## Q^2 max = "<<param.max[param.Q2_elastic]<<" GeV^2\n";
		}
		out <<"## Form factor parametrization = ";
		if (param.flag[param.form_factors] == 0) out <<"Simple Dipole\n";
		if (param.flag[param.form_factors] == 1) out <<"Dipole x Polynomial\n";
		if (param.flag[param.form_factors] == 2) out <<"Friedrich Walcher\n";
		if (param.flag[param.form_factors] == 3) out <<"Static Limit\n";
		if (param.flag[param.form_factors] == 4) out <<"User defined\n";
        out <<"## Degree of Polarization = "<<param.P*100.<<"%\n";
		out <<"## Calculate the asymmetry = ";
		if (param.flag[param.asymmetry] == 0) out<<"no\n";
		if (param.flag[param.asymmetry] == 1) {
			out<<"yes\n";
			out<<"## sin2thetaW = "<<param.sw2<<"\n";
			out<<"## Running weak mixing angle = ";
			if (param.flag[param.kappa_weak] == 0) out<<"0 - no running, fixed sin2thetaW\n";
			if (param.flag[param.kappa_weak] == 1) out<<"1 - effective sin2thetaW from CM\n";
		}
		out <<"## Maximum number of evaluations for 1st order bremsstrahlung = "<<param.MAXEVAL_1st<<"\n";
		if (param.flag[param.brems] == 2 || param.flag[param.brems] == 3)
			out <<"## Maximum number of evaluations for 2nd order bremsstrahlung = "<<param.MAXEVAL_2nd<<"\n";
		//	out <<"## Minimum number of evaluations = "<<param.MINEVAL<<"\n";
		//	out <<"## Relative accuracy = "<<param.EPSREL<<"\n";
		//	out <<"## Number of cores = "<<param.no_cores<<"\n";
		//	out <<"## Seed = "<<param.SEED<<"\n";
		out<<"##\n";
		out <<"## [E_gamma < Delta]\n";
		out <<"## Vacuum Polarization = ";
		if (param.flag[param.vac_pol] == 0) out<<"not included\n";
		if (param.flag[param.vac_pol] == 1) out<<"Only electron-positron loops\n";
		if (param.flag[param.vac_pol] == 2) out<<"Full leptonic contributions\n";
		if (param.flag[param.vac_pol] == 3) 
            out<<"Including hadronic contributions (Ignatov)\n";
		if (param.flag[param.vac_pol] == 4) 
            out<<"Including adronic contributions (Jegerlehner)\n";
		if (param.flag[param.vac_pol] == 5) 
            out<<"Including hadronic contributions (NKT18))\n";
		out <<"## Two-photon exchange correction (TPE) = ";
		if (param.flag[param.tpe] == 0) out<<"not included\n";
		if (param.flag[param.tpe] == 1) out<<"Calculation for a point like particle (Feshbach term)\n";
		if (param.flag[param.tpe] == 2) out<<"Tomalak total calculation\n";
		out<<"##\n";
		out <<"## [E_gamma > Delta]\n";
		out <<"## Type of hard-photon bremsstrahlung = ";
		if (param.flag[param.brems] == 0) out<<"no hard photon contribution\n";
		if (param.flag[param.brems] == 1) out<<"1st order\n";
		if (param.flag[param.brems] == 2) out<<"1st and 2nd order\n";
		if (param.flag[param.brems] == 3) out<<"2nd order\n";
        out <<"## Hadronic sof+virtual corrections = ";
        if (param.flag[param.hadr_corr] == 0) out<<"no hadronic corrections\n";
        if (param.flag[param.hadr_corr] == 1) out<<"only tpe and lep-had interference\n";
        if (param.flag[param.hadr_corr] == 2) out<<"only hadronic corrections\n";
        if (param.flag[param.hadr_corr] == 3) out<<"complete tpe + interference + hadronic\n";
		out <<"## Hadronic radiation = ";
		if (param.flag[param.brems_hadr] == 0) out<<"no contribution from hadronic radiation \n";
		if (param.flag[param.brems_hadr] == 1) out<<"only lep-had interference\n";
        if (param.flag[param.brems_hadr] == 2) out<<"only hadronic radiation\n";
        if (param.flag[param.brems_hadr] == 3) out<<"complete interference + hadronic\n";
        if (param.flag[param.brems_hadr] != param.flag[param.hadr_corr]) {
            out<<"***** WARNING: ***** \n"
               <<"***** Input for Hadronic corrections and Hadronic Radiation \n" 
               <<"***** do not match. Result unphysical and for testing only! \n";
        }
		out <<"## E_gamma max = "<<param.max[param.E_gamma] <<" GeV\n";
		out <<"## E' min = "<<param.min[param.E_prime] <<" GeV\n";
		out <<"## E' max = "<<param.max[param.E_prime] <<" GeV\n";
		out <<"## theta_gamma min = "<<param.min[param.theta_gamma_deg]<<" degrees\n";
		out <<"## theta_gamma max = "<<param.max[param.theta_gamma_deg]<<" degrees\n";
		//	out <<"## Q'^2 min = "<<param.min[param.Q2_prime]<<" GeV^2\n";
		//	out <<"## Q'^2 max = "<<param.max[param.Q2_prime]<<" GeV^2\n";
		out <<"\n"
            <<"##########################################################################\n"
            <<"                     Results of the numerical integration                 \n"
            <<"\n";

        out << "Leading-order results:\n";
        out << "Sigma unpol Born =                        " 
            << output.sigma_unpol_born
            << " +- "<< errors.sigma_unpol_born << " nb\n";
        out << "Sigma pol Born =                          " 
            << output.sigma_pol_born
            << " +- "<< errors.sigma_pol_born << " nb\n";
        out << "\n";
        out << "\n";
        
        out << "First-order corrections:\n";
        out << "Sigma unpol Born + soft-photon + 1-loop = " 
            << output.sigma_unpol_elastic_1st
            << " +- "<< errors.sigma_unpol_elastic_1st << " nb\n";
        out << "Sigma unpol 1 hard photon hadr interf=    "
            << output.sigma_unpol_inelastic_1st_hadr_interf
            << " +- " << errors.sigma_unpol_inelastic_1st_hadr_interf << " nb\n";
        out << "Sigma unpol 1 hard photon hadronic =      "
            << output.sigma_unpol_inelastic_1st_hadr
            << " +- " << errors.sigma_unpol_inelastic_1st_hadr << " nb\n";
        out << "Sigma unpol 1 hard photon (incl.hadr) =   " 
            << output.sigma_unpol_inelastic_1st
            << " +- "<< errors.sigma_unpol_inelastic_1st << " nb\n";
        out << "Sigma unpol 1st order (complete) =        " 
            << output.sigma_unpol_1st
            << " +- "<< errors.sigma_unpol_1st << " nb\n";
        out << "\n";
        out << "Sigma pol Born + soft-photon + 1-loop =   " 
            << output.sigma_pol_elastic_1st
            << " +- "<< errors.sigma_pol_elastic_1st << " nb\n";
        out << "Sigma pol 1 hard photon (leptonic only) = " 
            << output.sigma_pol_inelastic_1st
            << " +- "<< errors.sigma_pol_inelastic_1st << " nb\n";
        out << "Sigma pol 1st order (complete) =          " 
            << output.sigma_pol_1st
            << " +- " << errors.sigma_pol_1st << " nb\n";
        out << "\n";
        out << "\n";

        out << "Second-order corrections:\n";
        out << "Sigma unpol Born + 1st + 2 soft photons + 2-loop = " 
            << output.sigma_unpol_elastic_2nd
            << " +- "<< errors.sigma_unpol_elastic_2nd << " nb\n";
        out << "Sigma unpol 1 hard photon with loop =              " 
            << output.sigma_unpol_inelastic_loop
            << " +- "<< errors.sigma_unpol_inelastic_loop << " nb\n";
        out << "Sigma unpol 2 hard photons =                       " 
            << output.sigma_unpol_inelastic_2nd
            << " +- "<< errors.sigma_unpol_inelastic_2nd << " nb\n";
        out << "Sigma unpol hard photon + soft photon (fin) =      " 
            << output.sigma_unpol_2nd_add
            << " +- " << errors.sigma_unpol_2nd_add << " nb\n";
        out << "Sigma unpol 2nd order (complete) =                 " 
            << output.sigma_unpol_2nd
            << " +- "<< errors.sigma_unpol_2nd << " nb\n";
        out << "\n"; 

        out << "Sigma pol Born + 1st + 2 soft photons + 2-loop =   " 
            << output.sigma_pol_elastic_2nd
            << " +- "<< errors.sigma_pol_elastic_2nd << " nb\n";
        out << "Sigma pol 1 hard photon + with loop =              " 
            << output.sigma_pol_inelastic_loop
            << " +- "<< errors.sigma_pol_inelastic_loop << " nb\n";
        out << "Sigma pol 2 hard photons =                         " 
            << output.sigma_pol_inelastic_2nd
            << " +- "<< errors.sigma_pol_inelastic_2nd << " nb\n";
        out << "Sigma pol hard photon + soft photon (fin) =        " 
            << output.sigma_pol_2nd_add
            << " +- " << errors.sigma_pol_2nd_add << " nb\n";
        out << "Sigma pol 2nd order (complete) =                   " 
            << output.sigma_pol_2nd
            << " +- " << errors.sigma_pol_2nd << " nb\n";
        out << "\n"; 
        out << "\n"; 

        out << "Polarization asymmetry:\n";
        out << "PV asymmetry at leading order =                " 
            << output.asymm_born 
            << " +- "<< errors.asymm_born << "\n";
        out << "PV asymmetry including 1st order corrections = " 
            << output.asymm_1st 
            << " +- "<< errors.asymm_1st << "\n";
        out << "PV asymmetry including 2nd order corrections = " 
            << output.asymm_2nd 
            << " +- "<< errors.asymm_2nd << "\n";
        out << "1st order correction in percent = " 
            << output.rel_asymm_1st << "\n";
        out << "2nd order correction in percent = " 
            << output.rel_asymm_2nd << "\n";
        out << "\n"; 
        out << "\n"; 

		out <<"##########################################################################";


		out.close();
		param.output_file.erase (param.output_file.end()-4, param.output_file.end());
	}

// Echo output at terminal
	std::cout <<"\n"
              <<"##########################################################################\n"
			  <<"##                               POLARES "<< VERSION << "\n"
			  <<"##\n"
			  <<"##       Radiative Corrections for Polarized Electron-Proton Scattering   \n"
			  <<"##\n"
			  <<"##                              R.-D. Bucoveanu                           \n"
			  <<"##\n"
			  <<"##                         "<<dt
	          <<"##\n"
			  <<"## If you use POLARES please cite R.-D. Bucoveanu and H. Spiesberger,     \n"
			  <<"## Eur. Phys. J. A (2019) 55: 57, arXiv:1811.04970 [hep-ph].              \n"
	          <<"## Copyright (c) Razvan Bucoveanu, 2019. E-mail: rabucove@uni-mainz.de    \n";

	if (param.flag[param.echo_input] == 1) {

		std::cout <<"##########################################################################\n"
				<<"##                                  Input                                   \n"
				<<"##\n"
				<<"## [General Input]                                                          \n"
				<<"## Type of incident lepton = ";
		if (param.flag[param.lepton] == 0) std::cout<<"electron\n";
		if (param.flag[param.lepton] == 1) std::cout<<"positron\n";
		std::cout	<<"## Type of target particle = ";
		if (param.flag[param.target] == 0) std::cout<<"proton\n";
		if (param.flag[param.target] == 1) std::cout<<"carbon-12\n";
		std::cout <<"## Incident lepton energy = "<<param.en<<" GeV\n";
		std::cout <<"## Lower cut-off value for the photon energy (Delta) = "<<param.min[param.E_gamma] * 1000.<<" MeV\n";
		std::cout <<"## Type of cuts for elastic scattering = ";
		if (param.flag[param.cuts_born] == 0) {
			std::cout <<"Scattering angle (theta_l) cuts\n";
			std::cout <<"## theta_l min = "<<param.min[param.theta_l_deg]<<" degrees\n";
			std::cout <<"## theta_l max = "<<param.max[param.theta_l_deg]<<" degrees\n";
		}
		if (param.flag[param.cuts_born] == 1) {
			std::cout <<"Q^2 cuts\n";
			std::cout <<"Q^2 min = "<<param.min[param.Q2_elastic]<<" GeV^2\n";
			std::cout <<"Q^2 max = "<<param.max[param.Q2_elastic]<<" GeV^2\n";
		}
		std::cout <<"## Form factor parametrization = ";
		if (param.flag[param.form_factors] == 0) std::cout <<"Simple Dipole\n";
		if (param.flag[param.form_factors] == 1) std::cout <<"Dipole x Polynomial (Bernauer PhD thesis pp. 181)\n";
		if (param.flag[param.form_factors] == 2) std::cout <<"Friedrich Walcher\n";
		if (param.flag[param.form_factors] == 3) std::cout <<"Static Limit\n";
		if (param.flag[param.form_factors] == 4) std::cout <<"User defined\n";
		std::cout <<"## Degree of Polarization = "<<param.P*100.<<"%\n";
		std::cout <<"## Calculate the asymmetry = ";
		if (param.flag[param.asymmetry] == 0) std::cout<<"no\n";
		if (param.flag[param.asymmetry] == 1) {
			std::cout<<"yes\n";
			std::cout<<"## sin2thetaW = "<<param.sw2<<"\n";
			std::cout<<"## Running weak mixing angle = ";
			if (param.flag[param.kappa_weak] == 0) std::cout<<"0 - no running, fixed sin2thetaW\n";
			if (param.flag[param.kappa_weak] == 1) std::cout<<"1 - effective sin2thetaW from CM\n";
		}
        std::cout <<"## Maximum number of evaluations for 1st order bremsstrahlung = "<<param.MAXEVAL_1st<<"\n";
		if (param.flag[param.brems] == 2 || param.flag[param.brems] == 3)
			std::cout <<"## Maximum number of evaluations for 2nd order bremsstrahlung = "<<param.MAXEVAL_2nd<<"\n";
//		std::cout <<"## Minimum number of evaluations = "<<param.MINEVAL<<"\n";
//		std::cout <<"## Relative accuracy = "<<param.EPSREL<<"\n";
//		std::cout <<"## Number of cores = "<<param.no_cores<<"\n";
//		std::cout <<"## Seed = "<<param.SEED<<"\n";
		std::cout<<"##\n";
		std::cout <<"## [E_gamma < Delta]\n";
		std::cout <<"## Vacuum Polarization = ";
		if (param.flag[param.vac_pol] == 0) std::cout<<"not included\n";
		if (param.flag[param.vac_pol] == 1) std::cout<<"Only electron-positron loops\n";
		if (param.flag[param.vac_pol] == 2) std::cout<<"Full leptonic contributions\n";
		if (param.flag[param.vac_pol] == 3) 
            std::cout<<"Including hadronic contributions (Ignatov)\n";
        if (param.flag[param.vac_pol] == 4) 
            std::cout<<"Including hadronic contributions (Jegerlehner)\n";
        if (param.flag[param.vac_pol] == 5) {
            std::cout<<"Including hadronic contributions (KNT18)\n"; 
            std::cout<<"***** Grid file needed ! \n";
        }
		std::cout <<"## Two-photon exchange correction (TPE) = ";
		if (param.flag[param.tpe] == 0) std::cout<<"not included\n";
		if (param.flag[param.tpe] == 1) std::cout<<"Calculation for a point like particle (Feshbach term)\n";
		std::cout<<"##\n";
		std::cout <<"## [E_gamma > Delta]\n";
		std::cout <<"## Type of hard-photon bremsstrahlung = ";
		if (param.flag[param.brems] == 1) std::cout<<"1st order\n";
		if (param.flag[param.brems] == 2) std::cout<<"1st and 2nd order\n";
		if (param.flag[param.brems] == 3) std::cout<<"2nd order\n";
        std::cout <<"## Hadronic soft+virtual corrections = ";
        if (param.flag[param.hadr_corr] == 0) std::cout<<"no hadronic corrections\n";
        if (param.flag[param.hadr_corr] == 1) std::cout<<"only tpe and lep-had interference\n";
        if (param.flag[param.hadr_corr] == 2) std::cout<<"only hadronic corrections\n";
        if (param.flag[param.hadr_corr] == 3) std::cout<<"complete tpe + interference + hadronic\n";
		std::cout <<"## Hadronic radiation = ";
        if (param.flag[param.brems_hadr] == 0) std::cout<<"no contribution from hadronic radiation \n";
        if (param.flag[param.brems_hadr] == 1) std::cout<<"only lep-had interference\n";
        if (param.flag[param.brems_hadr] == 2) std::cout<<"only hadronic radiation\n";
        if (param.flag[param.brems_hadr] == 3) std::cout<<"complete interference + hadronic\n";
        if (param.flag[param.brems_hadr] != param.flag[param.hadr_corr]) {
            std::cout<<"***** WARNING: ***** \n"
                     <<"***** Input for Hadronic corrections and Hadronic Radiation \n" 
                     <<"***** do not match. Result unphysical and for testing only! \n";
        }
		std::cout <<"## E_gamma max = "<<param.max[param.E_gamma] <<" GeV\n";
		std::cout <<"## E' min = "<<param.min[param.E_prime] <<" GeV\n";
		std::cout <<"## E' max = "<<param.max[param.E_prime] <<" GeV\n";
		std::cout <<"## theta_gamma min = "<<param.min[param.theta_gamma_deg]<<" degrees\n";
		std::cout <<"## theta_gamma max = "<<param.max[param.theta_gamma_deg]<<" degrees\n";
	}
}


bool PES::change_energy_initialization(const double E){
	bool E_is_OK(true);

	if (E > 10.) {
		E_is_OK = false;
		std::cout << "Warning! Energy has to be smaller than 10 GeV for generating events\n\n";
	}

	param.en = E;
	return E_is_OK;
}


bool PES::change_energy_events(const double E){
	bool E_is_OK(true);

	param.en = round(E/param.Delta_E) * param.Delta_E;

	if (param.en - param.min[param.E] < 0. ||
			int(param.en/param.Delta_E) - int(param.min[param.E]/param.Delta_E)
			> output.ev_brems_1st.size()) {
		E_is_OK = false;
		std::cout << " Warning! The energy "<<param.en<<" GeV was not initialized\n";
	}

	return E_is_OK;
}


bool PES::set_child_process(const int child_process) {
	bool seed_is_OK(false);

	if (int(child_process) == child_process && child_process < 100) {
		seed_is_OK = true;
		if (!seed_is_set) {
			FS.seed += child_process;
			seed_is_set = true;
		}
	}

	return seed_is_OK;
}


double PES::get_E() {

	return param.en;
}
