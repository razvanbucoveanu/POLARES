#ifndef CONST_H
#define CONST_H

#include <cmath>
//----------------------------------------------------------------------------------------------
// File with constant values for POLARES
//==============================================================================================
namespace constants {
//----------------------------------------------------------------------------------------------
// Some mathematical and physical constants:
const double alpha = 1/137.035999084;       //fine-structure constant
const double pi = acos(-1.);                // value of pi
const double pi2 = pow(pi,2.);
const double mpr = 0.938272;                //proton mass
const double mpr2 = pow(mpr,2.);
const double mpr3 = pow(mpr,3.);
const double mpr4 = pow(mpr,4.);
const double mpr6 = pow(mpr,6.);
const double me = 0.0005109989;             // electron/positron mass
const double me2 = pow(me,2.);
const double me3 = pow(me,3.);
const double me4 = pow(me,4.);
const double me6 = pow(me,6.);
const double me8 = pow(me,8.);
const double me10 = pow(me,10.);
const double me12 = pow(me,12.);
const double m_mu = 105.658372e-3;         // muon/antimuon mass
const double m_tau = 1.77682;              // tau/antitau mass
const double m_carbon12 = 0.931494*12. - 6.*me;		// mass of carbon-12 in GeV
const double nb = 0.389379304*1E6;          //conversion to nanobarns factor
const double muep = 2.7928473;              //proton magnetic moment
const double muen = -1.9130427;             //neutron magnetic moment
const double gf = 0.00001166364;            //Fermi constant
const double sw2_msbar = 0.23122; //sinus of the weak mixing angle squared (MS-bar scheme)
//const double ga = -1./2.;        		    // Z^0 coupling constants
//const double gv = -1./2. + 2.*sw2;
const double mw = 80.385;
const double mz = 91.1871;
const double mz2 = pow(mz,2.);
const double mw2 = pow(mw,2.);
const double fm = 1/sqrt(nb/1e7); 					//fm to GeV^-1 conversion

//weak corrections
//double const rho_ew = 0.9878, k_ew = 1.0027;
}
#endif
