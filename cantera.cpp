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
auto cantera(double relative_tolerance, double absolute_tolerance, double& temperature, double& pressure, const char* mole_proportions, double time_step, size_t& species_len, const char**& species_data,
											double*& net_production_rates_data, double*& concentrations_data,
											size_t& reactions_len, const char**& reactions_data, double*& equilibrium_constants_data, double*& forward_rates_of_progress_data, double*& reverse_rates_of_progress_data
						) -> void try {
	using namespace Cantera;
	auto mechanism = newSolution("h2o2.yaml", "ohmech", "None");
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
	auto net_production_rates = new std::vector<double>();
	net_production_rates->resize(kinetics->nTotalSpecies());
	kinetics->getNetProductionRates(net_production_rates->data());
	net_production_rates_data = net_production_rates->data();

	reactions_len = kinetics->nReactions();
	auto reactions = new std::vector<const char*>();
	for(auto k: Range_new(kinetics->nReactions())) { reactions->push_back((new std::string(kinetics->reaction(k)->equation()))->data()); }
	reactions_data = reactions->data();
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

	system.advance(time_step);

	temperature = phase->temperature();
	pressure = phase->pressure();
	auto concentrations = new std::vector<double>();
	concentrations->resize(kinetics->nTotalSpecies());
	phase->getConcentrations(concentrations->data());
	concentrations_data = concentrations->data();

} catch (std::exception& err) { std::cerr << err.what() << std::endl; }
