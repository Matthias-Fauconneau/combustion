#include <vector>
#include <string>
#include <cantera/kinetics/Kinetics.h>
#include <cantera/thermo/IdealGasPhase.h>
#include <cantera/zerodim.h>
#include <cantera/transport/MultiTransport.h>
using namespace std;

extern "C"
void reaction(double& pressure, double& temperature, const char* mole_proportions,
											size_t& species_len, const char**& species_data,
											double rate_time,
												double*& net_production_rates_data,
												size_t& reactions_len,
													const char**& equations_data,
													double*& equilibrium_constants_data,
													double*& forward_rates_of_progress_data,
													double*& reverse_rates_of_progress_data,
											double state_time,
											double*& concentrations_data) try {
	using namespace Cantera;
	auto mechanism = newSolution("gri30.yaml", "gri30", "mixture-averaged"/*Multi*/);
	auto kinetics = mechanism->kinetics();
	species_len = kinetics->nTotalSpecies();
	auto species = new std::vector<const char*>();
	for(auto k=0; k<kinetics->nTotalSpecies(); k++) { species->push_back((new std::string(kinetics->kineticsSpeciesName(k)))->data()); }
	species_data = species->data();
	auto phase = mechanism->thermo();
	phase->setState_TPX(temperature, pressure, mole_proportions);
	IdealGasConstPressureReactor reactor;
	reactor.insert(mechanism);
	ReactorNet system;
	//system.setTolerances(relative_tolerance, absolute_tolerance);
	system.setTolerances(/*relative_tolerance:*/ 1e-8, /*absolute_tolerance:*/ 1e-14);
	system.addReactor(reactor);
	system.advance(rate_time);

	auto net_production_rates = new std::vector<double>();
	net_production_rates->resize(kinetics->nTotalSpecies());
	kinetics->getNetProductionRates(net_production_rates->data());
	net_production_rates_data = net_production_rates->data();

	reactions_len = kinetics->nReactions();
	auto reactions = new std::vector<const char*>();
	for(auto k=0; k<kinetics->nReactions(); k++) { reactions->push_back((new std::string(kinetics->reaction(k)->equation()))->data()); }
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
	concentrations_data = concentrations->data();
} catch (std::exception& err) { std::cerr << err.what() << std::endl; }

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
