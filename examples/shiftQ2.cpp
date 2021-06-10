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
#include <sys/time.h>
#include "config.h"
#include "const.h"

using std::cout;
using std::cin;
using std::endl;

int main()
{
	time_t time_start;
	time_t time_end;

	time(&time_start);
	char* dt = ctime(&time_start); // Current time

	cout.precision(14);

	std::ofstream out;
	out.precision(6);
	out.open("Q2shift_155.dat");
	out << "thl\tshift_45\tshift_50\tshift_55\n";

	using namespace POLARES;

	Input input;

	double E_min = 0.045;
	double E_max = 0.055;
	double Delta_E = 0.005;

	for (double thl=1.;thl<=178.;thl++) {

		cout << thl << "\n";
		out << thl << "\t";
		for (double E = E_min; E <= E_max; E += Delta_E) {
			cout << E << "\n";
			PES pes;
			input.E_prime_min = E;
			pes.set_input(input);
			pes.shiftQ2(thl);
			cout << "Shift Q2 = " << pes.output.shiftQ2 <<" +- "<< pes.errors.shiftQ2 << endl;
			cout << "Relative Shift Q2 = " << pes.output.rel_shiftQ2 <<" +- "<<
					pes.errors.rel_shiftQ2 << endl;
			out << pes.output.rel_shiftQ2 <<"\t";
		}
		out << "\n";
		cout << endl;
	}

	time(&time_end);
	cout<<"\nTime: "<<difftime(time_end,time_start)<<" seconds "<<"\n\n";

	return 0;
}
