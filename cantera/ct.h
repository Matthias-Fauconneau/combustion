#include <stddef.h>
int thermo_newFromFile(const char* file_name, const char* phase_name);
size_t thermo_nSpecies(int n);
int thermo_setTemperature(int n, double t);
int thermo_setMoleFractions(int n, size_t len, const double* x, int norm);
int thermo_getSpeciesName(int n, size_t m, size_t len, /*out*/ char* buffer);
int thermo_setPressure(int n, double p);
int kin_newFromFile(const char* file_name, const char* phase_name, int reactingPhase, int neighbor0, int neighbor1, int neighbor2, int neighbor3);
int kin_getEquilibriumConstants(int n, size_t len, double* kc);
int kin_getNetProductionRates(int n, size_t len, double* w_dot);
//int kin_getCreationRates(int n, size_t len, /*out*/ double* cdot);
//int kin_getDestructionRates(int n, size_t len, /*out*/ double* ddot);
int kin_getNetRatesOfProgress(int n, size_t len, double* net_ROP);
int kin_getReactionString(int n, size_t i, size_t len, /*out*/ char* buffer);
