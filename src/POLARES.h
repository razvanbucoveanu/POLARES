#ifndef POLARES_H
#define POLARES_H

#include "interpolation.h"
#include "IO_classes.h"
#include "gsl_rand.h"
#include "parameters.h"
#include "cuba_param.h"
#include "integrands.h"

/**
\dot
digraph G {
main -> POLARES namespace;
POLARES namespace -> PES (POLARES Class)
}
\enddot
*/

namespace POLARES {

/**
 * main class
 */
class PES {
//Hidden variables
protected:
	/**
	 * Storing for all the required parameters
	 */
	Parameters param;
	/**
	 * Interpolation of vacuum polarization and two-photon exchange corrections
	 */
	Interpolation interpolation;
	/**
	 * Class for storing user defined input
	 */
	Input input;
	/**
	 * GSL random number generator
	 */
	Rand rand;
	/**
	 * Cuba related parameters
	 */
	Cuba_parameters CP;
	/**
	 * Class that contains all the integrands
	 */
	Integrands integrands;
	/**
	 * Lepton and target particle masses
	 */
	mutable double m, M, m2, M2;

//Public variables
public:
	/**
	 * A class in which the results of the numerical integration are placed
	 */
	Output output;
	/**
	 * A class in which the uncertainties of the numerical integration are placed
	 */
	Output errors;

	/**
	 * A class in which all the information about the final state particles
	 * is stored when events are generated
	 */
	Final_State FS;

//Hidden functions
protected:

	//Copy Constructor
	PES(const PES&);
	//Copy operator
	PES& operator=(const PES&);
	/**
	 * A function that sets the final values for the final state particles
	 */
	void set_final_state();
	void write_output();
	bool seed_is_set;

//Public functions
public:
	/**
	 * Constructor that can receive external user-defined input
	 * sets the input
	 * @param input
	 */
	PES(const Input& input);
	/**
	 * Standard Constructor
	 * sets the input
	 */
	PES();
	//Destructor
	~PES();
	/**
	 * Function that can be used to change the energy before the initialization
	 * @param
	 * E (energy of the incoming electron)
	 * @return
	 * true if the value can be used to generate events
	 */
	bool change_energy_initialization(const double E);
	/**
	 * Function that can be used to change the energy for generating events
	 * @param
	 * E (energy of the incoming electron)
	 * @return
	 * true if the value was initialized
	 */
	bool change_energy_events(const double E);
	/**
	 * Function for changing the events seed in case of multi-core usage
	 * @param child_process
	 * @return true if the value is different from other seeds
	 */
	bool set_child_process(const int child_process);
//	const Input& get_input()const{
//		return input;
//	}
    /**
     * Function that can be used for the initialization of the grid in case
	   events are generated or simply to obtain the cross sections or asymmetries
     * @return results and errors
     */
	int initialization();
	/**
	 * Function that generates events. Can be used only after the initialization
	was performed. Warning! The current version works only with energies between 20-200 MeV
	with an increment of 0.1 MeV
	 * @return phase_space
	 */
	int events();
	/**
	 * Function that calculates the shift in Q2 due to first order hard-photon bremsstrahlung for a given scattering angle theta _l
	 * @return results and errors
	 */
	int shiftQ2(const double thl_deg);
	/**
	 * Function that calculates the one-fold differential cross-section and asymmetry in respect
	to the scattering angle of the final electron
	 * @return results and errors
	 */
	int sigma_diff_Omega_l(const double thl_deg);
	/**
	 * Function that calculates the total (virtual + real) first order radiative corrections
	 * @return results and errors
	 */
	double delta(const double thl_deg);
	/**
	 * Function that calculates the running of the weak mixing angle
	 * @return results
	 */
	double running_sw2(const double Q2);
	/**
	 * Function that can be used to replace the input from "POLARES.in"
	 * @return void
	 */
	void set_input(const Input& input);
	/**
	 * Function that gives the total cross section for a given energy
	 * @return sigma_unpol_vect
	 */
	double get_total_unpol_cross_section(const double E);
	/**
	 * Function that echoes the energy provided by the user
	 * @return sigma_unpol_vect
	 */
	double get_E();
};

}  // namespace POLARES

#endif
