//#include <vector>
//#include <string>
/*#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/base/Solution.h"*/
/*#include <cantera/kinetics/Kinetics.h>
#include <cantera/thermo/IdealGasPhase.h>
#include <cantera/zerodim.h>
#include <cantera/transport/MultiTransport.h>*/
#include "cantera/clib/ct.h"
using namespace std;

extern "C"
void reaction(double& pressure, double& temperature,
							//const double* mole_fractions,
							const char* mole_proportions,
											size_t& species_len, //const char**& species_data,
											//double rate_time,
												double*& net_production_rates,
												/*size_t& reactions_len,
													const char**& equations,
													double*& equilibrium_constants,
													double*& forward_rates_of_progress,
													double*& reverse_rates_of_progress,
											double state_time,
											double*& concentrations*/) {//try {
	//using namespace Cantera;
	//auto mechanism = newSolution("gri30.yaml", "gri30", "mixture-averaged"/*Multi*/);
	//auto thermo = mechanism->thermo();
	//auto phase = newPhase("gri30.yaml", "gri30");
	auto phase = thermo_newFromFile("gri30.yaml", "gri30");
	//species_len = phase->nSpecies();
	species_len = thermo_nSpecies(phase);
	auto species = new const char*[species_len];
	for(auto k=0; k<species_len; k++) { species[k] = new char[8]; thermo_getSpeciesName(phase, k, 8, species[k]); }
	//species_data = species->data();
	//phase->setState_TPX(temperature, pressure, mole_proportions);
	thermo_setTemperature(phase, temperature);
	thermo_setPressure(phase, pressure);
	//thermo_setMoleFractions(phase, species_len, mole_fractions);
	thermo_setMoleFractionsByName(phase, mole_proportions);
#if 0
	IdealGasReactor reactor;
	//IdealGasConstPressureReactor reactor;
	reactor.insert(mechanism);
	ReactorNet system;
	//system.setTolerances(relative_tolerance, absolute_tolerance);
	system.setTolerances(/*relative_tolerance:*/ 1e-8, /*absolute_tolerance:*/ 1e-14);
	system.addReactor(reactor);
	system.advance(rate_time);
#endif

	net_production_rates = new double[species_len];
	/*std::vector<ThermoPhase*> phases;
	phases.push_back(phase);
	//sol->setKinetics();
	auto kinetics = newKinetics(phases, "gri30.yaml", "gri30"); //mechanism->kinetics();*/
	auto kinetics = kin_newFromFile("gri30.yaml", "gri30", phase, 0, 0, 0, 0);
	//kinetics->getNetProductionRates(net_production_rates);
	kin_getNetProductionRates(kinetics, species_len, net_production_rates);

	/*reactions_len = kinetics->nReactions();
	auto reactions = new std::vector<const char*>();
	//for(auto k=0; k<kinetics->nReactions(); k++) { reactions->push_back((new std::string(kinetics->reaction(k)->equation()))->data()); }
	equations_data = reactions->data();
	auto equilibrium_constants = new std::vector<double>();
	equilibrium_constants->resize(kinetics->nReactions());
	kinetics->getEquilibriumConstants(equilibrium_constants->data());
	equilibrium_constants_data = equilibrium_constants->data();
	auto forward_rates_of_progress = new std::vector<double>();
	forward_rates_of_progress->resize(kinetics->nReactions());
	kinetics->getFwdRatesOfProgress(forward_rates_of_progress->data());
	forward_rates_of_progress_data = forward_rates_of_progress->data();
	auto reverse_rates_of_progress = new std::vector<double>();
	reverse_rates_of_progress->resize(kinetics->nReactions());
	kinetics->getRevRatesOfProgress(reverse_rates_of_progress->data());
	reverse_rates_of_progress_data  = reverse_rates_of_progress->data();

	system.advance(state_time);

	temperature = phase->temperature();
	pressure = phase->pressure();
	auto concentrations = new std::vector<double>();
	concentrations->resize(kinetics->nTotalSpecies());
	phase->getConcentrations(concentrations->data());
	concentrations_data = concentrations->data();*/
} //catch (std::exception& err) { std::cerr << err.what() << std::endl; }

#if 0
extern "C"
void transport(double pressure, double temperature, const char* mole_proportions, double& viscosity, double& thermal_conductivity, size_t& species_len, const char**& species_data, double*& mixture_averaged_thermal_diffusion_coefficients_data) try {
	using namespace Cantera;
	auto mechanism = newSolution("gri30.yaml", "gri30", "mixture-averaged"/*Multi*/);
	auto kinetics = mechanism->kinetics();
	species_len = kinetics->nTotalSpecies();
	auto species = new std::vector<const char*>();
	for(auto k=0; k<kinetics->nTotalSpecies(); k++) { species->push_back((new std::string(kinetics->kineticsSpeciesName(k)))->data()); }
	species_data = species->data();
	auto phase = mechanism->thermo();
	phase->setState_TPX(temperature, pressure, mole_proportions);
	auto transport = mechanism->transport();
	viscosity = transport->viscosity();
	thermal_conductivity = transport->thermalConductivity();
	auto mixture_averaged_thermal_diffusion_coefficients = new std::vector<double>();
	mixture_averaged_thermal_diffusion_coefficients->resize(kinetics->nTotalSpecies());
	transport->getThermalDiffCoeffs(mixture_averaged_thermal_diffusion_coefficients->data());
	mixture_averaged_thermal_diffusion_coefficients_data = mixture_averaged_thermal_diffusion_coefficients->data();
} catch (std::exception& err) { std::cerr << err.what() << std::endl; }
#endif
