__device__ double sq(double x) { return x*x; }
__device__ double dot5(double a[5], double b[5]) { double sum = 0.; for(uint k=0;k<5;k++) { sum += a[k]*b[k]; } return sum; }
__device__ double eval_poly(double P[5], double x) {
	double generate_k[5]; for(uint k=0;k<5;k++) generate_k[k] = pow(x, k);
	return dot5(P, generate_k);
}

__device__ double dot(double a[SPECIES], double b[SPECIES]) { double sum = 0.; for(uint k=0;k<SPECIES;k++) { sum += a[k]*b[k]; } return sum; }

struct TransportPolynomials {
	double sqrt_viscosity_T14[SPECIES][5];
	double thermal_conductivity_T12[SPECIES][5];
	double binary_thermal_diffusion_coefficients_T32[SPECIES][SPECIES][5];
};

struct Species {
	double molar_mass[SPECIES];
	double thermodynamics[SPECIES][2][7];
	TransportPolynomials transport_polynomials;
};

__device__ double sqrt_viscosity(Species species, uint a, double T) { return sqrt(sqrt(T)) * eval_poly(species.transport_polynomials.sqrt_viscosity_T14[a], log(T)); }
__device__ double thermal_conductivity(Species species, uint a, double T) { return sqrt(T) * eval_poly(species.transport_polynomials.thermal_conductivity_T12[a], log(T)); }
__device__ double binary_thermal_diffusion_coefficient(Species species, uint a, uint b, double T) { return pow(T,3./2.) * eval_poly(species.transport_polynomials.binary_thermal_diffusion_coefficients_T32[a>b?a:b][a>b?b:a], log(T)); }

extern "C" __global__ void rates_transport(const uint len, const double pressure_R, double* temperature, double* _amounts, double* d_temperature, double* d_amounts, double* viscosity, double* _thermal_conductivity, double* mixture_averaged_thermal_diffusion_coefficients) {
const Species species = Species{
#include "species.h"
};
const double volume = 1.;
for (uint i = blockIdx.x * /*workgroup size*/blockDim.x + /*SIMD lane*/threadIdx.x; i < len; i += /*workgroup size*/blockDim.x * /*SIMD width*/gridDim.x) {
	double T = temperature[i];
	double amounts[SPECIES]; for(int k=0;k<SPECIES;k++) amounts[k] = _amounts[k*len+i];
	double generate_k[SPECIES]; for(uint k=0;k<SPECIES;k++) {
		double generate_j[SPECIES];
		for(uint j=0;j<SPECIES;j++) {
			generate_j[j] = sq(1. + (sqrt_viscosity(species, k, T)/sqrt_viscosity(species, j, T)) * sqrt(sqrt(species.molar_mass[j]/species.molar_mass[k]))) /
			(sqrt(8.) * sqrt(1. + species.molar_mass[k]/species.molar_mass[j]));
		}
		generate_k[k] = sq(sqrt_viscosity(species, k, T)) / dot(amounts, generate_j);
	}
	viscosity[i] = dot(amounts, generate_k);
#if 0
	double amount = pressure_R / T * volume;
	double thermal_conductivities[SPECIES]; for(uint k=0;k<SPECIES;k++) { thermal_conductivities[k] = thermal_conductivity(species, k, T); }
	double rcp_thermal_conductivities[SPECIES]; for(uint k=0;k<SPECIES;k++) { rcp_thermal_conductivities[k] = 1. / thermal_conductivity(species, k, T); }
	_thermal_conductivity[i] = 1./2. * (dot(amounts, thermal_conductivities) / amount + amount / dot(amounts, rcp_thermal_conductivities));
	for(uint k=0;k<SPECIES;k++) {
		double generate_j[SPECIES];
		for(uint j=0;j<SPECIES;j++) {
			if (j != k) {
				generate_j[j] = 1. / binary_thermal_diffusion_coefficient(species, k, j, T);
			} else {
				generate_j[j] = 0.;
			}
		}
		mixture_averaged_thermal_diffusion_coefficients[k*len+i] = (1. - amounts[k]/amount) / dot(amounts, generate_j);
	}
#endif
}
}
