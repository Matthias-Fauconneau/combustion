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

extern "C"  auto cantera(double rtol, double atol, double T, double P, const char* X, size_t& len, const char**& species_data, double*& dtw_data) -> void {
	using namespace Cantera;
	auto mechanism = newSolution("gri30.yaml", "gri30", "None");
	auto phase = mechanism->thermo();
	phase->setState_TPX(T, P, X);
	IdealGasConstPressureReactor reactor;
	reactor.insert(mechanism);
	ReactorNet system;
	system.setTolerances(rtol, atol);
	system.addReactor(reactor);
	//system.advance(dt);
	auto kinetics = mechanism->kinetics();
	//new (&species) std::vector<std::string>;
	//for(auto k: Range_new(kinetics->nTotalSpecies())) { species.push_back(kinetics->kineticsSpeciesName(k)); }
	auto species = new std::vector<const char*>();
	for(auto k: Range_new(kinetics->nTotalSpecies())) { species->push_back((new std::string(kinetics->kineticsSpeciesName(k)))->data()); }
	species_data = species->data();
	len = kinetics->nTotalSpecies();
	//new (&dtw) std::vector<double>;
	auto dtw = new std::vector<double>();
	dtw->resize(kinetics->nTotalSpecies());
	kinetics->getNetProductionRates(dtw->data());
	dtw_data = dtw->data();
}
