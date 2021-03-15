#include <stddef.h>
int thermo_newFromFile(const char* filename, const char* phasename);
size_t thermo_nSpecies(int n);
int thermo_setTemperature(int n, double t);
int thermo_setMoleFractions(int n, size_t lenx, const double* x, int norm);
int thermo_getSpeciesName(int n, size_t m, size_t lennm, /*out*/ char* nm);
int thermo_setPressure(int n, double p);
//int thermo_setDensity(int n, double rho);
int kin_newFromFile(const char* filename, const char* phasename, int reactingPhase, int neighbor1, int neighbor2, int neighbor3, int neighbor4);
int kin_getNetProductionRates(int n, size_t len, /*out*/ double* wdot);
