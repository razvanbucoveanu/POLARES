#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <string>
//#include "gsl/gsl_interp.h"
//#include "gsl/gsl_spline.h"
namespace POLARES {
class Interpolation {
protected:
//	gsl_interp_accel *acc;
//	gsl_spline *spline;
	void *gsl_interp_accel_ptr;
	void *gsl_interp_accel_ptr_tpe_Tomalak;
	void *gsl_spline_ptr;
	void *gsl_spline_ptr_tpe_Tomalak;
	void *gsl_spline2d_ptr;
	void *gsl_xacc_ptr;
	void *gsl_yacc_ptr;
//	Interpolation& operator=(const Interpolation &I);

public:
	void init_d_vac_hadr(const std::string& filename);
	void init_tpe_Tomalak(const std::string& filename);
	void init_tpe(const std::string& filename);
	double d_vac_hadr(const double Q2) const;
	double tpe_Tomalak(const double Q2) const;
	double tpe(const double Q2, const double eps) const;

	Interpolation();
//	:acc(gsl_interp_accel_alloc())
//	,spline(NULL){
//
//	}

	Interpolation(const Interpolation &I);
//	:acc(gsl_interp_accel_alloc())
//	,spline(NULL){
//
//	}


	~Interpolation();
//	{
//		gsl_interp_accel_free(acc);
//		if(spline) gsl_spline_free(spline);
//	}

};
}
#endif
