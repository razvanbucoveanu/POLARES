/*
 * gsl_rand.h
 *
 *  Created on: Dec 11, 2015
 *      Author: razvan
 */

#ifndef SRC_GSL_RAND_H_
#define SRC_GSL_RAND_H_


#include <iostream>

namespace POLARES {

class Rand{
public:
	enum GSL_RAND_TYPE {
		GSL_RAND_TYPE_mt19937
		,GSL_RAND_TYPE_ranlxs0,GSL_RAND_TYPE_ranlxs1,GSL_RAND_TYPE_ranlxs2
		,GSL_RAND_TYPE_ranlxd1,GSL_RAND_TYPE_ranlxd2
		,GSL_RAND_TYPE_ranlux,GSL_RAND_TYPE_ranlux389
		,GSL_RAND_TYPE_cmrg
		,GSL_RAND_TYPE_mrg
		,GSL_RAND_TYPE_taus,GSL_RAND_TYPE_taus2
		,GSL_RAND_TYPE_gfsr4
	};
protected:
	unsigned long int seed;
	GSL_RAND_TYPE type;

	unsigned int* stack;
	const void * gsl_T;	// gsl_random_number_generator Type
	mutable void * gsl_G;	// gsl_random_number_generator

	void _init();
	void _exit();


public:

	inline unsigned int get_seed()const{return seed;}
	inline GSL_RAND_TYPE get_type()const{return type;}
	void change_seed(unsigned long int in_seed){seed = in_seed;}

	inline Rand(const unsigned long int seed=time(NULL), const GSL_RAND_TYPE type=GSL_RAND_TYPE_ranlux)
	:seed(seed),type(type)
	,stack(new unsigned int(0))
	,gsl_T(NULL),gsl_G(NULL){
		_init();
	}

	inline Rand(const Rand& R)
	:seed(R.seed),type(R.type)
	,stack(R.stack)
	,gsl_T(R.gsl_T),gsl_G(R.gsl_G){
		++(*stack);
	}

	inline Rand& operator=(const Rand& R){
		_exit();

		seed	= (R.seed);
		type	= (R.type);
		stack	= (R.stack);
		gsl_T	= (R.gsl_T);
		gsl_G	= (R.gsl_G);
		++(*stack);

		return *this;
	}

	inline ~Rand(){
		_exit();
	}

	// returns double in [0,1)
	double uniform()const;
	// returns double in [min,max)
	inline double uniform(const double min, const double max)const{
		return (min+(max-min)*uniform());
	}
	// returns double in (0,1)
	double uniform_pos()const;
	// returns double in (min,max)
	inline double uniform_pos(const double min, const double max)const{
		return (min+(max-min)*uniform_pos());
	}

	// returns uniform_pos (0,1)
	inline double operator ()()const{
		return uniform_pos();
	}
	// returns uniform_pos (min,max)
	inline double operator ()(const double min, const double max)const{
		return uniform_pos(min,max);
	}

	// returns int in [0,max]
	long int uniform_int(const long int max)const;

	std::ostream& to_stream(std::ostream& os)const;
};

}  // namespace POLARES

#endif //GSL_RAND_H_
