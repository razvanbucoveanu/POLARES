#include "interpolation.h"
#include <fstream>
#include <iostream>
#include <vector>

#include <gsl/gsl_math.h>
#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

using POLARES::Interpolation;

Interpolation::Interpolation()
:gsl_interp_accel_ptr(static_cast<void*>(gsl_interp_accel_alloc()))
,gsl_interp_accel_ptr_tpe_Tomalak(static_cast<void*>(gsl_interp_accel_alloc()))
,gsl_spline_ptr(NULL)
,gsl_spline_ptr_tpe_Tomalak(NULL)
,gsl_xacc_ptr(static_cast<void*>(gsl_interp_accel_alloc()))
,gsl_yacc_ptr(static_cast<void*>(gsl_interp_accel_alloc()))
,gsl_spline2d_ptr(NULL){
}

Interpolation::~Interpolation(){
	gsl_interp_accel_free(static_cast<gsl_interp_accel*>(gsl_interp_accel_ptr));
	if(gsl_spline_ptr_tpe_Tomalak)	gsl_spline_free(static_cast<gsl_spline*>(gsl_spline_ptr_tpe_Tomalak));
	gsl_interp_accel_free(static_cast<gsl_interp_accel*>(gsl_interp_accel_ptr_tpe_Tomalak));
	if(gsl_spline_ptr)	gsl_spline_free(static_cast<gsl_spline*>(gsl_spline_ptr));
	gsl_interp_accel_free(static_cast<gsl_interp_accel*>(gsl_xacc_ptr));
	gsl_interp_accel_free(static_cast<gsl_interp_accel*>(gsl_yacc_ptr));
	if(gsl_spline2d_ptr) gsl_spline2d_free(static_cast<gsl_spline2d*>(gsl_spline2d_ptr));
}

void Interpolation::init_tpe(const std::string& filename){

	// std::cout <<"size: "<< gsl_interp2d_type_min_size(gsl_interp2d_bicubic) <<"\n";

	std::ifstream infile(filename.c_str());

	const gsl_interp2d_type *T = gsl_interp2d_bicubic;

	double Q2, eps, tpe;
	unsigned int i;
	double aux_Q2, aux_eps;
	std::vector<double> x, y, z;

	infile >> aux_Q2 >> aux_eps >> tpe;
	infile.clear();
	infile.seekg(0, std::ios::beg);

	while (infile >> Q2 >> eps >> tpe) {
		if (Q2 == aux_Q2) y.push_back(eps);
		if (eps == aux_eps) x.push_back(Q2);
		z.push_back(tpe);
	}
	infile.close();

	double xa[x.size()], ya[y.size()];

	// std::cout <<"test: "<<y.size()<<"\n";

	for (i = 0; i < x.size() ; i++)
		xa[i] = x[i];
	for (i = y.size(); i >= 0 ; i--)
		ya[i] = y[i];

	const size_t nx = sizeof(xa) / sizeof(double); /* x grid points */
	const size_t ny = sizeof(ya) / sizeof(double); /* y grid points */
	double *za = (double*) malloc(nx * ny * sizeof(double));

	for (i = 0; i < z.size() ; i++)
		za[i] = z[i];

	if(gsl_spline2d_ptr) gsl_spline2d_free(static_cast<gsl_spline2d*>(gsl_spline2d_ptr));
	gsl_spline2d_ptr = static_cast<void*>(gsl_spline2d_alloc (T, nx, ny));
	gsl_spline2d_init(static_cast<gsl_spline2d*>(gsl_spline2d_ptr), xa, ya, za, nx, ny);
}

void Interpolation::init_d_vac_hadr(const std::string& filename){

	std::ifstream infile(filename.c_str());

	std::vector<double> x, y;
	double Q2, vp;

	while (infile >> Q2 >> vp) {
		x.push_back(Q2);
		y.push_back(vp);
	}
	infile.close();

	double Q2_vp[x.size()], vp_h[y.size()];

	for (unsigned int i = 0; i < x.size(); i++) {
		Q2_vp[i] = x[i];
		vp_h[i] = y[i];
	}

//	if(spline) gsl_spline_free(spline);
//	spline = gsl_spline_alloc (gsl_interp_cspline, N_points);
//	gsl_spline_init (spline, Q2_vp, vp_h, N_points);
	if(gsl_spline_ptr)	gsl_spline_free(static_cast<gsl_spline*>(gsl_spline_ptr));
	gsl_spline_ptr	= static_cast<void*>(gsl_spline_alloc (gsl_interp_cspline, x.size()));
	gsl_spline_init(static_cast<gsl_spline*>(gsl_spline_ptr), Q2_vp, vp_h, x.size());
}

void Interpolation::init_tpe_Tomalak(const std::string& filename){

	std::ifstream infile(filename.c_str());

	std::vector<double> x, y;
	double Q2, eps, tpe;

	while (infile >> Q2 >> eps >> tpe) {
		x.push_back(Q2);
		y.push_back(tpe);
	}
	infile.close();

	double Q2_array[x.size()], tpe_array[y.size()];

	for (unsigned int i = 0; i < x.size(); i++) {
		Q2_array[i] = x[i];
		tpe_array[i] = y[i];
		}

//	if(spline) gsl_spline_free(spline);
//	spline = gsl_spline_alloc (gsl_interp_cspline, N_points);
//	gsl_spline_init (spline, Q2_vp, vp_h, N_points);
	if (gsl_spline_ptr_tpe_Tomalak)	gsl_spline_free(static_cast<gsl_spline*>(gsl_spline_ptr_tpe_Tomalak));
	gsl_spline_ptr_tpe_Tomalak = static_cast<void*>(gsl_spline_alloc (gsl_interp_cspline, x.size()));
	gsl_spline_init(static_cast<gsl_spline*>(gsl_spline_ptr_tpe_Tomalak), Q2_array, tpe_array, x.size());
}

double Interpolation::d_vac_hadr(const double Q2)const {
//	return 2.*gsl_spline_eval (spline, Q2, acc);
	return 2.*gsl_spline_eval (static_cast<gsl_spline*>(gsl_spline_ptr), Q2,static_cast<gsl_interp_accel*>(gsl_interp_accel_ptr));
}

double Interpolation::tpe_Tomalak(const double Q2)const {
	return gsl_spline_eval (static_cast<gsl_spline*>(gsl_spline_ptr_tpe_Tomalak), Q2,static_cast<gsl_interp_accel*>
																					(gsl_interp_accel_ptr_tpe_Tomalak))/100.;
}

double Interpolation::tpe(const double Q2, const double eps)const {
	return gsl_spline2d_eval (static_cast<gsl_spline2d*>(gsl_spline2d_ptr), Q2, eps,
			static_cast<gsl_interp_accel*>(gsl_xacc_ptr), static_cast<gsl_interp_accel*>(gsl_yacc_ptr));
}
