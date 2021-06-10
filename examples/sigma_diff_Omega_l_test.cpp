/*
 * main.cpp
 *
 *  Created on: Mar 25, 2015
 *  Author: Razvan
 *  Modified: HS, Mar 15-16, 2017
 *
 */
#include "POLARES.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>
#include <sstream>
#include <cstdio>
#include <sys/time.h>

using std::cout;
using std::cin;
using std::endl;

int main()
{
	time_t time_start;
	time_t time_end;

	time(&time_start);
	char* dt = ctime(&time_start); // Current time

	//Output
	// Print to file?
	//	std::ofstream out;
	//	out.open("this-run-output.dat");

	//Number of digits for output
	cout.precision(10);

	//Print header
	cout<<"Test sigma_diff_thl \n\n";

	using namespace POLARES;

	Input input;

	input.thl_min = 25.;
	input.thl_max = 45.;
	input.E_prime_min = 0.045;
	input.E_min = 0.150;
	input.E_max = 0.160;
	input.Delta_E = 0.001;

	PES pes;
	pes.set_input(input);
	pes.initialization();

	cout << "Summary total cross sections: \n";

	if (pes.output.sigma_unpol_elastic_1st != 0) cout << "Sigma soft-photon 1st order = " << pes.output.sigma_unpol_elastic_1st << endl;
	if (pes.output.sigma_unpol_elastic_2nd != 0) cout << "Sigma soft-photon 2nd order = " << pes.output.sigma_unpol_elastic_2nd << endl;
	if (pes.output.sigma_unpol_inelastic_1st != 0) cout << "Sigma hard-photon unpol 1 gamma = " << pes.output.sigma_unpol_inelastic_1st
			<< " +- " << pes.errors.sigma_unpol_inelastic_1st << endl;
	if (pes.output.sigma_unpol_inelastic_loop != 0) cout << "Sigma hard-photon + 1-loop = " << pes.output.sigma_unpol_inelastic_loop
			<< " +- " << pes.errors.sigma_unpol_inelastic_loop << endl;
	if (pes.output.sigma_unpol_inelastic_2nd != 0) cout << "Sigma hard-photon unpol 2 gamma = " << pes.output.sigma_unpol_inelastic_2nd
			<< " +- " << pes.errors.sigma_unpol_inelastic_2nd << endl;
	if (pes.output.sigma_unpol_1st != 0) cout << "Sigma 1st order = " << pes.output.sigma_unpol_1st
			<< " +- " << pes.errors.sigma_unpol_1st << endl;
	if (pes.output.sigma_unpol_2nd != 0) cout << "Sigma 2nd order = " << pes.output.sigma_unpol_2nd
			<< " +- " << pes.errors.sigma_unpol_2nd << endl;
	if (pes.output.sigma_unpol_born != 0) cout << "Sigma Born = " << pes.output.sigma_unpol_born << endl;
	if (pes.output.sigma_pol_1st != 0) cout << "Sigma pol 1st order = " << pes.output.sigma_pol_1st
			<< " +- " << pes.errors.sigma_pol_1st << endl;
	if (pes.output.sigma_pol_2nd != 0) cout << "Sigma pol 2nd order = " << pes.output.sigma_pol_2nd
			<< " +- " << pes.errors.sigma_pol_2nd << endl;
	if (pes.output.asymm_1st != 0) cout << "Asymm 1st order = " << pes.output.asymm_1st <<
			" +- " << pes.errors.asymm_1st << endl;
	if (pes.output.asymm_2nd != 0) cout << "Asymm 2nd order = " << pes.output.asymm_2nd <<
			" +- " << pes.errors.asymm_2nd << endl;
	if (pes.output.asymm_born != 0) cout << "Asymm Born = " << pes.output.asymm_born << endl;


	time(&time_end);
	cout<<"Time for initialization: "<<difftime(time_end,time_start)<<" seconds "<<"\n\n";

	double thetai;
	for (thetai = input.thl_min; thetai <= input.thl_max; thetai += 5.) {
		pes.sigma_diff_Omega_l(thetai);
		cout << "Scattering angle = " << thetai << " degrees" << endl;
		cout << "Sigma Born thl = " << pes.output.sigma_unpol_born << endl;
		cout << "Sigma at 1st order thl = " << pes.output.sigma_unpol_1st << endl;
		cout << "Sigma at 2nd order thl = " << pes.output.sigma_unpol_2nd << endl;
		cout << "Asymmetry Born thl = " << pes.output.asymm_born << endl;
		cout << "Asymmetry at 1st order thl = " << pes.output.asymm_1st<< endl;
		cout << "Asymmetry at 2nd order thl = " << pes.output.asymm_2nd<< endl;
		cout << "Relative correction of asymmetry 1st = " << pes.output.rel_asymm_1st << endl;
		cout << endl;
	}

	return 0;
}
