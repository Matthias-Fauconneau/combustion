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

extern "C"
auto cantera(double relative_tolerance, double absolute_tolerance, double& temperature, double& pressure, const char* mole_proportions, double time_step, size_t& species_len, const char**& species_data, double*& concentrations_data) -> void try {
	using namespace Cantera;
	auto mechanism = newSolution("gri30.yaml", "gri30", "None");
	auto phase = mechanism->thermo();
	phase->setState_TPX(temperature, pressure, mole_proportions);
	IdealGasConstPressureReactor reactor;
	reactor.insert(mechanism);
	ReactorNet system;
	system.setTolerances(relative_tolerance, absolute_tolerance);
	system.addReactor(reactor);
	auto kinetics = mechanism->kinetics();
	species_len = kinetics->nTotalSpecies();
	auto species = new std::vector<const char*>();
	for(auto k: Range_new(kinetics->nTotalSpecies())) { species->push_back((new std::string(kinetics->kineticsSpeciesName(k)))->data()); }
	species_data = species->data();
	system.advance(time_step);
	temperature = phase->temperature();
	pressure = phase->pressure();
	auto concentrations = new std::vector<double>();
	concentrations->resize(kinetics->nTotalSpecies());
	phase->getConcentrations(concentrations->data());
	concentrations_data = concentrations->data();
} catch (std::exception& err) { std::cerr << err.what() << std::endl; }
//catch(CanteraError& cterr) { std::cerr << cterr.what()) << std::endl; }
