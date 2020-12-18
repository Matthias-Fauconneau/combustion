#include <vector>
#include <string>
using namespace std;

struct IndexIterator { size_t index;

	size_t operator*() { return index; }
	void operator++() { index += 1; }
	bool operator !=(const IndexIterator& o) const { return index != o.index; }
};

struct Range {
	public: size_t start;
	public: size_t end_;

	inline IndexIterator begin() const { return {start}; }
	inline IndexIterator end() const { return {end_}; }
};
auto Range_new(size_t size) -> Range { return {.start=0, .end_=size}; }

#include <cantera/zerodim.h>
#include <cantera/thermo/IdealGasPhase.h>
#include <cantera/transport/MultiTransport.h>

extern "C"
void cantera(double pressure, double temperature, const char* mole_proportions, double& viscosity, double& thermal_conductivity, size_t& species_len, const char**& species_data, double*& thermal_diffusion_coefficients_data) try {
	using namespace Cantera;
	auto mechanism = newSolution("gri30.yaml", "gri30", "Multi");
	auto kinetics = mechanism->kinetics();
	species_len = kinetics->nTotalSpecies();
	auto species = new std::vector<const char*>();
	for(auto k: Range_new(kinetics->nTotalSpecies())) { species->push_back((new std::string(kinetics->kineticsSpeciesName(k)))->data()); }
	species_data = species->data();
	auto phase = mechanism->thermo();
	phase->setState_TPX(temperature, pressure, mole_proportions);
	auto transport = mechanism->transport();
	viscosity = transport->viscosity();
	thermal_conductivity = transport->thermalConductivity();
	auto thermal_diffusion_coefficients = new std::vector<double>();
	thermal_diffusion_coefficients->resize(kinetics->nTotalSpecies());
	transport->getThermalDiffCoeffs(thermal_diffusion_coefficients->data());
	thermal_diffusion_coefficients_data = thermal_diffusion_coefficients->data();
} catch (std::exception& err) { std::cerr << err.what() << std::endl; }
