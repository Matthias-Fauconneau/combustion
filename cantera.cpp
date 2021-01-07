#include <vector>
#include <string>
#include <cantera/kinetics/Kinetics.h>
#include <cantera/thermo/IdealGasPhase.h>
#include <cantera/transport/MultiTransport.h>
using namespace std;

extern "C"
void cantera(double pressure, double temperature, const char* mole_proportions, double& viscosity, double& thermal_conductivity, size_t& species_len, const char**& species_data, double*& mixture_averaged_thermal_diffusion_coefficients_data) try {
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
