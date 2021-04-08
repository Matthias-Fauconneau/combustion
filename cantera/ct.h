#include <stddef.h>
int thermo_newFromFile(const char* file_name, const char* phase_name);
size_t thermo_nSpecies(int n);
int thermo_getMolecularWeights(int n, size_t lenm, double* mw);
int thermo_setTemperature(int n, double t);
int thermo_setMoleFractions(int n, size_t len, const double* x, int norm);
int thermo_getSpeciesName(int n, size_t m, size_t len, /*out*/ char* buffer);
int thermo_setPressure(int n, double p);
int thermo_getEnthalpies_RT(int n, size_t len, double* H_RT);
int thermo_getCp_R(int n, size_t len, double* Cp_R);
int kin_newFromFile(const char* file_name, const char* phase_name, int reactingPhase, int neighbor0, int neighbor1, int neighbor2, int neighbor3);
int kin_getEquilibriumConstants(int n, size_t len, double* kc);
int kin_getNetProductionRates(int n, size_t len, double* w_dot);
int kin_getFwdRatesOfProgress(int n, size_t len, double* rate);
int kin_getRevRatesOfProgress(int n, size_t len, double* rate);
int kin_getNetRatesOfProgress(int n, size_t len, double* rate);
int kin_getReactionString(int n, size_t i, size_t len, /*out*/ char* buffer);
int trans_newDefault(int th, int loglevel);
double trans_viscosity(int n);
double trans_thermalConductivity(int n);
int trans_getThermalDiffCoeffs(int n, int ldt, double* dt); // Mass-averaged
int trans_getBinDiffCoeffs(int n, int ld, double* d);
