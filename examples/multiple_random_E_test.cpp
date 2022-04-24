/*
 * multiple_random_E_test.cpp
 *
 *  Created on: Mar 23, 2017
 *      Author: razvan
 */

#include "POLARES.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>
#include <sstream>
#include <cstdio>
#include <stdlib.h>
#include <sys/time.h>

using std::cout;
using std::cin;
using std::endl;

int main()
{
	time_t time_start;
	time_t time_end;
	time_t time_events;

	cout<<"\nNumber of events?\n";
	int no_events;
	cin>>no_events;

//	cout<<"\nName of the input file?\n";
//	std::string input_file;
//	cin>>input_file;

	time(&time_start);
	char* dt = ctime(&time_start); // Current time

	cout.precision(14);

	cout<<"Beware: initialisation for this run will take a while (2 to 3 min)\n\n";

	using namespace POLARES;

	Input input;
//	input.input_file = input_file;
	input.thl_min = 25.;
	input.thl_max = 45.;
	input.Delta = 0.01;
	input.E_prime_min = 0.045;
	input.E_min = 0.155;
	input.E_max = 0.155;
	input.Delta_E = 0.001;
	input.flag[input.order] = 2;
	input.flag[input.brems] = 2;
	input.flag[input.asymmetry] = 0;
	input.flag[input.brems_add] = 0;
	input.flag[input.GL] = 1;

	PES pes;
	pes.set_input(input);

	double E;

	cout <<"Initialisation \n\n";
	for (E = input.E_min; E <= input.E_max; E += input.Delta_E) {
		if (pes.change_energy_initialization(E))
			pes.initialization();

		cout <<"Incident electron energy = "<< E <<" GeV "<< endl;
		if (pes.output.sigma_unpol_elastic_1st != 0) cout << "Sigma soft-photon 1st order = " << pes.output.sigma_unpol_elastic_1st << endl;
		if (pes.output.sigma_unpol_elastic_2nd != 0) cout << "Sigma soft-photon 2nd order = " << pes.output.sigma_unpol_elastic_2nd << endl;
		if (pes.output.sigma_unpol_inelastic_1st != 0) cout << "Sigma hard-photon unpol 1 gamma = " << pes.output.sigma_unpol_inelastic_1st
				<< " +- " << pes.errors.sigma_unpol_inelastic_1st << endl;
		if (pes.output.sigma_unpol_inelastic_loop != 0) cout << "Sigma hard-photon + 1-loop = " << pes.output.sigma_unpol_inelastic_1st
				<< " +- " << pes.errors.sigma_unpol_inelastic_1st << endl;
		if (pes.output.sigma_unpol_inelastic_2nd != 0) cout << "Sigma hard-photon unpol 2 gamma = " << pes.output.sigma_unpol_inelastic_2nd
				<< " +- " << pes.errors.sigma_unpol_inelastic_2nd << endl;
		if (pes.output.sigma_unpol_1st != 0) cout << "Sigma 1st order = " << pes.output.sigma_unpol_1st
				<< " +- " << pes.errors.sigma_unpol_1st << endl;
		if (pes.output.sigma_unpol_2nd != 0) cout << "Sigma 2nd order = " << pes.output.sigma_unpol_2nd
				<< " +- " << pes.errors.sigma_unpol_2nd << endl;
		if (pes.output.sigma_unpol_born != 0) cout << "Sigma Born = " << pes.output.sigma_unpol_born << endl;

		cout << endl;

	}

	time(&time_events);
	cout<<"Time for initialization: "<<difftime(time_events,time_start)<<" seconds "<<"\n\n";

	cout <<"Event generation \n\n";
	if (no_events != 0) {

		FILE* fout;
		fout = fopen("events.txt", "w");

		for (int i = 0; i < no_events; i++) {

			double Ei=input.E_min+(input.E_max-input.E_min)*((double) rand() / (RAND_MAX));
			if (pes.change_energy_events(Ei))
				pes.events();

				if (pes.FS.k_2[0] != 0.)
			fprintf(fout,
					" %7ll   %8.14lf  %8.14lf  %8.14lf %8.14lf %8.14lf %8.14lf  %8.14lf %8.14lf %8.14lf %8.14lf  %8.14lf %8.14lf %8.14lf %8.14lf %8.14lf %8.14lf %8.14lf %8.14lf %4.10lf\n",
					pes.FS.event_no,
					Ei, pes.FS.l_1[0],
					pes.FS.l_2[0],pes.FS.l_2[1],pes.FS.l_2[2],pes.FS.l_2[3],
					pes.FS.p_2[0],pes.FS.p_2[1],pes.FS.p_2[2],pes.FS.p_2[3],
					pes.FS.k_1[0],pes.FS.k_1[1],pes.FS.k_1[2],pes.FS.k_1[3],
					pes.FS.k_2[0],pes.FS.k_2[1],pes.FS.k_2[2],pes.FS.k_2[3],
					pes.FS.weight
			);

		}

		time(&time_end);
		cout<<"Time for event generation: "<<difftime(time_end,time_events)<<" seconds "<<"\n\n";
	}

	return 0;
}
