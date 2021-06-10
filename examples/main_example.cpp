/*
 * main.cpp
 *
 *  Created on: Mar 25, 2015
 *      Author: Razvan
 */
#include "POLARES.h" // Main header file
#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include <sys/time.h>
#include "config.h"
#include "const.h"
#include "string"

using std::cout;
using std::cin;
using std::endl;

//-------------------------------------------------------------------------------------------------
//Main program
int main()
{
	time_t time_start;
	time_t time_end;
	time_t time_events;

	cout<<"\nNumber of events?\n";
	int no_events;
	cin>>no_events;

	cout<<"\nName of the input file?\n";
	std::string input_file;
	cin>>input_file;

	time(&time_start);
	char* dt = ctime(&time_start); // Current time
//	Output
	//-------------------------------------------------------------------------------------------------
//	Precision of output
	using namespace POLARES;

//	input class declaration
	Input input;
	input.input_file = input_file;
//	input.thl_min = 25.;
//	input.thl_max = 45.;
//	input.E = 0.155;
//	input.E_prime_min = 0.045;
//	input.Delta = 1e-13;
//	input.E_prime_max = 1e10;
//	input.E_gamma_max = 0.11;
//	input.no_cores = 1;
//	input.thg_min = 0.;
//	input.thg_max = 180.;
//	input.flag[input.asymmetry] = 0;
//	input.flag[input.order] = 1;
//	input.flag[input.brems] = 4;
//	input.flag[input.GL] = 0.;
//	input.flag[input.LO] = 0;
//	input.flag[input.brems_add] = 0;
//	input.flag[input.int_method] = 0;
//	input.no_eval_1st = 1e7;
//	input.no_eval_gamma_loop = 1e7;
//	input.no_eval_2nd_add = 1e7;
//	input.lambda = constants::me2;
//	input.Delta_eps = 0.;
//	input.thl_deg = 35.;
//	input.thg_deg = 35.;

	PES pes;
	pes.set_input(input);
	pes.initialization();
//	pes.sigma_diff_Omega_l(35.);

	cout.precision(14);
	cout <<"################################################################################\n"
			<<"                     Numerical integration results                          \n\n";
	if (pes.output.sigma_unpol_elastic_1st != 0) cout << "Sigma unpol soft-photon 1st order = " << pes.output.sigma_unpol_elastic_1st
			<< " +- " << pes.errors.sigma_unpol_elastic_1st << endl;
	if (pes.output.sigma_unpol_elastic_2nd != 0) cout << "Sigma unpol soft-photon 2nd order = " << pes.output.sigma_unpol_elastic_2nd
			<< " +- " << pes.errors.sigma_unpol_elastic_2nd<< endl;
	if (pes.output.sigma_unpol_inelastic_1st != 0) cout << "Sigma unpol hard-photon 1 gamma = " << pes.output.sigma_unpol_inelastic_1st
			<< " +- " << pes.errors.sigma_unpol_inelastic_1st << endl;
	if (pes.output.sigma_unpol_inelastic_1st_hadr_interf != 0) cout << "Sigma unpol hard-photon 1 gamma hadr interf= "
			<< pes.output.sigma_unpol_inelastic_1st_hadr_interf
			<< " +- " << pes.errors.sigma_unpol_inelastic_1st_hadr_interf << endl;
	if (pes.output.sigma_unpol_inelastic_1st_hadr != 0) cout << "Sigma unpol hard-photon 1 gamma purely hadr = "
			<< pes.output.sigma_unpol_inelastic_1st_hadr
			<< " +- " << pes.errors.sigma_unpol_inelastic_1st_hadr << endl;
	if (pes.output.sigma_unpol_inelastic_loop != 0) cout << "Sigma unpol hard-photon + 1-loop = " << pes.output.sigma_unpol_inelastic_loop
			<< " +- " << pes.errors.sigma_unpol_inelastic_loop << endl;
	if (pes.output.sigma_unpol_inelastic_2nd != 0) cout << "Sigma unpol hard-photon 2 gamma = " << pes.output.sigma_unpol_inelastic_2nd
			<< " +- " << pes.errors.sigma_unpol_inelastic_2nd << endl;
	if (pes.output.sigma_unpol_1st != 0) cout << "Sigma unpol 1st order = " << pes.output.sigma_unpol_1st
			<< " +- " << pes.errors.sigma_unpol_1st << endl;
	if (pes.output.sigma_unpol_2nd != 0) cout << "Sigma unpol 2nd order = " << pes.output.sigma_unpol_2nd
			<< " +- " << pes.errors.sigma_unpol_2nd << endl;
	if (pes.output.sigma_unpol_2nd_add != 0) cout << "Sigma hard-photon + soft photon additional = " << pes.output.sigma_unpol_2nd_add
			<< " +- " << pes.errors.sigma_unpol_2nd_add << endl;
	if (pes.output.sigma_unpol_born != 0) cout << "Sigma unpol Born = " << pes.output.sigma_unpol_born
			<< " +- " << pes.errors.sigma_unpol_born<< endl;
	if (pes.output.sigma_pol_elastic_1st != 0) cout << "Sigma pol soft-photon 1st order = " << pes.output.sigma_pol_elastic_1st
			<< " +- " << pes.errors.sigma_pol_elastic_1st << endl;
	if (pes.output.sigma_pol_elastic_2nd != 0) cout << "Sigma pol soft-photon 2nd order = " << pes.output.sigma_pol_elastic_2nd
			<< " +- " << pes.errors.sigma_pol_elastic_2nd << endl;
	if (pes.output.sigma_pol_inelastic_1st != 0) cout << "Sigma pol hard-photon 1 gamma = " << pes.output.sigma_pol_inelastic_1st
			<< " +- " << pes.errors.sigma_pol_inelastic_1st << endl;
	if (pes.output.sigma_pol_inelastic_loop != 0) cout << "Sigma pol hard-photon + 1-loop = " << pes.output.sigma_pol_inelastic_loop
			<< " +- " << pes.errors.sigma_pol_inelastic_loop << endl;
	if (pes.output.sigma_pol_inelastic_2nd != 0) cout << "Sigma pol hard-photon 2 gamma = " << pes.output.sigma_pol_inelastic_2nd
			<< " +- " << pes.errors.sigma_pol_inelastic_2nd << endl;
	if (pes.output.sigma_pol_1st != 0) cout << "Sigma pol 1st order = " << pes.output.sigma_pol_1st
			<< " +- " << pes.errors.sigma_pol_1st << endl;
	if (pes.output.sigma_pol_2nd != 0) cout << "Sigma pol 2nd order = " << pes.output.sigma_pol_2nd
			<< " +- " << pes.errors.sigma_pol_2nd << endl;
	if (pes.output.sigma_pol_born != 0) cout << "Sigma pol Born = " << pes.output.sigma_pol_born << endl;
	if (pes.output.sigma_pol_2nd_add != 0) cout << "Sigma hard-photon + soft photon pol additional = " << pes.output.sigma_pol_2nd_add
			<< " +- " << pes.errors.sigma_pol_2nd_add << endl;
	if (pes.output.asymm_1st != 0) cout << "Asymm 1st order = " << pes.output.asymm_1st <<
			" +- " << pes.errors.asymm_1st << endl;
	if (pes.output.asymm_2nd != 0) cout << "Asymm 2nd order = " << pes.output.asymm_2nd <<
			" +- " << pes.errors.asymm_2nd << endl;
	if (pes.output.asymm_born != 0) cout << "Asymm Born = " << pes.output.asymm_born << endl;
	if (pes.output.rel_asymm_1st != 0) cout << "Relative Asymmetry 1st = " << pes.output.rel_asymm_1st << endl;
	if (pes.output.rel_asymm_2nd != 0) cout << "Relative Asymmetry 2nd = " << pes.output.rel_asymm_2nd << endl;

	time(&time_events);
	cout<<"\nInitialization Time: "<<difftime(time_events,time_start)<<" seconds "<<"\n\n";

	cout << "NLO = " << pes.output.sigma_unpol_1st - pes.output.sigma_unpol_born <<"\n";

//	Warning!
	//Current version of event generator works only for energies between 20 and 200 MeV with an increment of 0.1 MeV
	if (no_events != 0) {

		FILE* fout;
		fout = fopen("events.txt", "w");

		double E = pes.get_E();
		double avg_wgt = 0.;

		for (int i = 0; i < no_events; i++) {

			if (pes.change_energy_events(E))
			pes.events();

			fprintf(fout, "%d %8.4lf %8.7lf %6.4lf %7.4lf %8.4lf %6.4lf %7.4lf %8.4lf %6.4lf %7.4lf "
					"%8.4lf %6.4lf %7.4lf %4.10lf \n",
					pes.FS.event_no,
					pes.FS.E, pes.FS.E_prime_l, pes.FS.theta_l, pes.FS.phi_l,
					pes.FS.E_p, pes.FS.theta_p, pes.FS.phi_p,
					pes.FS.E_gamma, pes.FS.theta_gamma, pes.FS.phi_gamma,
					pes.FS.E_gamma_prime, pes.FS.theta_gamma_prime, pes.FS.phi_gamma_prime,
					pes.FS.weight);
		}

		time(&time_end);
		cout << "failed events: " << pes.FS.failed_ev << "\n";
		cout << "efficiency: " <<100 - pes.FS.failed_ev / double (no_events) * 100. << "\n";
		cout << std::setprecision(15)<< "Events Time: "<<difftime(time_end,time_events)<<" seconds "<<"\n\n";
	}
//	}
	return 0;
}
