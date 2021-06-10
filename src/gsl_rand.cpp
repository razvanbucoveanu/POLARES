/*
 * gsl_rand.cpp
 *
 *  Created on: Dec 11, 2015
 *      Author: razvan
 */

#include "gsl_rand.h"
#ifdef __cplusplus
extern "C" {
#endif
#include <gsl/gsl_rng.h>
#ifdef __cplusplus
}
#endif

using POLARES::Rand;

void Rand::_init(){
	switch (type) {
		case GSL_RAND_TYPE_mt19937:
			gsl_T	= gsl_rng_mt19937;
			break;
		case GSL_RAND_TYPE_ranlxs0:
			gsl_T	= gsl_rng_ranlxs0;
			break;
		case GSL_RAND_TYPE_ranlxs1:
			gsl_T	= gsl_rng_ranlxs1;
			break;
		case GSL_RAND_TYPE_ranlxs2:
			gsl_T	= gsl_rng_ranlxs2;
			break;
		case GSL_RAND_TYPE_ranlxd1:
			gsl_T	= gsl_rng_ranlxd1;
			break;
		case GSL_RAND_TYPE_ranlxd2:
			gsl_T	= gsl_rng_ranlxd2;
			break;
		case GSL_RAND_TYPE_ranlux:
			gsl_T	= gsl_rng_ranlux;
			break;
		case GSL_RAND_TYPE_ranlux389:
			gsl_T	= gsl_rng_ranlux389;
			break;
		case GSL_RAND_TYPE_cmrg:
			gsl_T	= gsl_rng_cmrg;
			break;
		case GSL_RAND_TYPE_mrg:
			gsl_T	= gsl_rng_mrg;
			break;
		case GSL_RAND_TYPE_taus:
			gsl_T	= gsl_rng_taus;
			break;
		case GSL_RAND_TYPE_taus2:
			gsl_T	= gsl_rng_taus2;
			break;
		case GSL_RAND_TYPE_gfsr4:
			gsl_T	= gsl_rng_gfsr4;
			break;
		default:
			gsl_T	= NULL;
			break;
	}
	if(gsl_T==NULL){
		std::cerr	<<"Wrong GSL Random Number Generator Type -> exit(1)"<<std::endl;
		std::cerr	<<"Available generators:"<<std::endl;
		for(const gsl_rng_type** t = gsl_rng_types_setup(); *t != 0; t++){
			std::cerr	<<"\t"<<(*t)->name<<std::endl;
		}
		exit(1);
	}
	gsl_G	= static_cast<void *>(gsl_rng_alloc(static_cast<const gsl_rng_type*>(gsl_T)));
	gsl_rng_set(static_cast<gsl_rng*>(gsl_G),seed);
}

void Rand::_exit(){
	if((*stack)>0){
		--(*stack);
	}else{
		delete stack;
		gsl_rng_free(static_cast<gsl_rng*>(gsl_G));
	}
}

double Rand::uniform()const{
	return gsl_rng_uniform(static_cast<gsl_rng*>(gsl_G));
}
double Rand::uniform_pos()const{
	return gsl_rng_uniform_pos(static_cast<gsl_rng*>(gsl_G));
}

long int Rand::uniform_int(const long int max)const{
	return gsl_rng_uniform_int(static_cast<gsl_rng*>(gsl_G),max+1);
}

std::ostream& Rand::to_stream(std::ostream& os)const{
	return os << "Random generator type:" << gsl_rng_name(static_cast<gsl_rng*>(gsl_G)) << "(" << seed << ")";
}
